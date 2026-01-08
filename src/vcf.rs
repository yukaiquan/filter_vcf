use anyhow::{Context, Result};
use clap::Parser;
use regex::Regex;
use std::collections::HashMap;
use std::io::{BufRead};
use crate::args::Args;

/// 基因型统计结果结构体
#[derive(Debug, Default)]
struct GenotypeStats {
    a_count: u32,    // 0/0 纯合参考基因型数量
    b_count: u32,    // 1/1 纯合替代基因型数量
    h_count: u32,    // 0/1 杂合基因型数量
    n_count: u32,    // ./. 缺失基因型数量
    present: u32,    // 有效基因型数量（a + b + h）
}

/// 从字符串中提取INFO字段的DP数值
fn extract_info_dp(input: &str, pattern: &Regex) -> u32 {
    pattern
        .captures(input)
        .and_then(|caps| caps.get(1))
        .map(|m| m.as_str().parse::<u32>().unwrap_or(0))
        .unwrap_or(0)
}

/// 动态解析FORMAT字段，返回字段名到索引的映射
fn parse_format_fields(format_str: &str) -> HashMap<&str, usize> {
    let mut field_map = HashMap::new();
    for (idx, field) in format_str.split(':').enumerate() {
        field_map.insert(field, idx);
    }
    field_map
}

/// 提取样本的DP、r值，同时保留所有原始字段
fn extract_sample_info(
    gt_str: &str,
    format_map: &HashMap<&str, usize>
) -> (u32, f64, Vec<String>) {
    // 关键修复：将&str迭代器转为String迭代器后再collect
    let all_fields: Vec<String> = gt_str.split(':')
        .map(|s| s.to_string()) // 每个&str转为String
        .collect();

    // 提取DP（用于判断基因型，动态找DP索引）
    let dp = format_map.get("DP")
        .and_then(|&idx| all_fields.get(idx))
        .and_then(|s| if s == "." { None } else { s.parse::<u32>().ok() })
        .unwrap_or(0);

    // 提取AD计算r值（仅用于判断，不修改AD字段）
    let r = if dp > 0 {
        let dv = format_map.get("AD")
            .and_then(|&idx| all_fields.get(idx))
            .map(|ad_str| {
                if ad_str == "." {
                    0
                } else {
                    let ad_parts: Vec<&str> = ad_str.split(',').collect();
                    match ad_parts.len() {
                        1 => ad_parts[0].parse::<u32>().unwrap_or(0),
                        _ => ad_parts[1].parse::<u32>().unwrap_or(0),
                    }
                }
            })
            .unwrap_or(0);
        dv as f64 / dp as f64
    } else {
        0.0
    };

    (dp, r, all_fields)
}

/// 处理单个样本的基因型
fn process_sample_gt(
    sample_str: &str,
    alt_base: &str,
    format_map: &HashMap<&str, usize>,
    dphom: u32,
    dphet: u32,
    tol: f64,
    stats: &mut GenotypeStats
) -> String {
    // 提取DP、r值，以及样本列的所有原始字段
    let (dp, r, mut all_fields) = extract_sample_info(sample_str, format_map);

    // 找到GT字段的索引
    let gt_idx = match format_map.get("GT") {
        Some(&idx) => idx,
        None => {
            // 无GT字段则返回原始样本字符串
            return sample_str.to_string();
        }
    };

    // 确保索引有效
    if gt_idx >= all_fields.len() {
        return sample_str.to_string();
    }

    // 仅修改GT字段的值，其他字段完全保留
    let new_gt = if dp == 0 && alt_base == "." {
        stats.a_count += 1;
        "0/0"
    } else if dp >= dphom && r <= tol {
        stats.a_count += 1;
        "0/0"
    } else if dp >= dphom && r >= 1.0 - tol {
        stats.b_count += 1;
        "1/1"
    } else if dp >= dphet && r >= 0.5 - tol && r <= 0.5 + tol {
        stats.h_count += 1;
        "0/1"
    } else {
        stats.n_count += 1;
        "./."
    };

    // 替换GT字段的值
    all_fields[gt_idx] = new_gt.to_string();

    // 重新拼接所有字段（保留所有原始字段，仅修改GT）
    all_fields.join(":")
}

/// 核心处理逻辑：FORMAT列保持原始，样本列仅修改GT字段，保留所有其他字段
pub fn process_vcf_line(
    line: &str,
    args: &Args,
    dp_re: &Regex,
) -> Result<Option<String>> {
    let parts: Vec<&str> = line.split('\t').collect();
    if parts.len() < 10 {
        return Ok(None); // 列数不足，跳过无效行
    }

    // 1. 基础字段提取
    let ref_base = parts[3];
    let alt_base = parts[4];
    let qual = parts[5].parse::<f64>().unwrap_or(0.0);
    let info = parts[7];          // 保留原始INFO
    let format_str = parts[8];    // 原始FORMAT列
    let info_dp = extract_info_dp(info, dp_re);

    // 前置过滤（多碱基、低QUAL等）这里把indel类型去除了，仅对SNP进行过滤
    if ref_base.len() > 1
        || alt_base.len() > 1
        || qual < args.minqual
        || ref_base == "N"
        || info_dp < args.mindp
    {
        return Ok(None);
    }

    // 2. 解析FORMAT（仅用于找GT/DP/AD索引）
    let format_map = parse_format_fields(format_str);
    // 必须包含GT和DP字段
    if !format_map.contains_key("GT") || !format_map.contains_key("DP") {
        return Ok(None);
    }

    // 3. 处理样本列
    let mut stats = GenotypeStats::default();
    let mut modified_samples = Vec::new();
    for sample_str in parts.iter().skip(9) {
        let new_sample = process_sample_gt(
            sample_str,
            alt_base,
            &format_map,
            args.dphom,
            args.dphet,
            args.tol,
            &mut stats
        );
        modified_samples.push(new_sample);
    }

    // 4. 群体/MAF过滤
    stats.present = stats.a_count + stats.b_count + stats.h_count;
    let present_ratio = if (stats.present + stats.n_count) > 0 {
        stats.present as f64 / (stats.present + stats.n_count) as f64
    } else {
        0.0
    };
    let hom_ratio = if stats.present > 0 {
        (stats.a_count + stats.b_count) as f64 / stats.present as f64
    } else {
        0.0
    };

    if stats.present == 0
        || stats.a_count < args.minhomn
        || stats.b_count < args.minhomn
        || present_ratio < args.minpresent
        || hom_ratio < args.minhomp
    {
        return Ok(None);
    }

    let maf = if stats.b_count > stats.a_count {
        (2 * stats.a_count + stats.h_count) as f64 / (2 * stats.present) as f64
    } else {
        (2 * stats.b_count + stats.h_count) as f64 / (2 * stats.present) as f64
    };
    if maf < args.minmaf {
        return Ok(None);
    }

    // 5. 构建最终行（FORMAT列保持原始，样本列保留所有字段）
    let mut new_parts = Vec::with_capacity(parts.len());
    // 前8列：保留原始值（转为String）
    for part in &parts[0..8] {
        new_parts.push(part.to_string());
    }
    // FORMAT列：100%保留原始值
    new_parts.push(format_str.to_string());
    // 样本列：仅修改GT，保留所有其他字段
    new_parts.extend(modified_samples);

    Ok(Some(new_parts.join("\t")))
}

/// 生成过滤规则注释行
pub fn generate_filter_comment(args: &Args) -> String {
    format!(
        "##FilterRule=<ID=CustomFilter,Description=\"dphom={}, dphet={}, tol={}, minqual={}, mindp={}, minhomn={}, minpresent={}, minhomp={}, minmaf={}\">",
        args.dphom, args.dphet, args.tol, args.minqual, args.mindp, args.minhomn, args.minpresent, args.minhomp, args.minmaf
    )
}