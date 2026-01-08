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

/// 动态解析FORMAT字段，返回字段名到索引的映射（仅用于找DP/AD）
fn parse_format_fields(format_str: &str) -> HashMap<&str, usize> {
    let mut field_map = HashMap::new();
    for (idx, field) in format_str.split(':').enumerate() {
        field_map.insert(field, idx);
    }
    field_map
}

/// 提取样本中的DP和用于计算r值的DV（仅用于逻辑判断，不输出）
fn extract_sample_dp_and_r_value(gt_str: &str, format_map: &HashMap<&str, usize>) -> (u32, f64) {
    let gt_parts: Vec<&str> = gt_str.split(':').collect();

    // 提取DP（从FORMAT映射中动态找，兼容任意位置）
    let dp = format_map.get("DP")
        .and_then(|&idx| gt_parts.get(idx))
        .and_then(|s| if *s == "." { None } else { s.parse::<u32>().ok() })
        .unwrap_or(0);

    // 提取AD计算r值（仅用于判断基因型，不输出）
    let r = if dp > 0 {
        let dv = format_map.get("AD")
            .and_then(|&idx| gt_parts.get(idx))
            .map(|ad_str| {
                if *ad_str == "." {
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

    (dp, r)
}

/// 处理单个样本的基因型（仅修改GT部分，保留DP，不输出DV）
fn process_sample_gt(
    sample_str: &str,
    alt_base: &str,
    format_map: &HashMap<&str, usize>,
    dphom: u32,
    dphet: u32,
    tol: f64,
    stats: &mut GenotypeStats
) -> String {
    let (dp, r) = extract_sample_dp_and_r_value(sample_str, format_map);

    // 仅修改GT部分，DP保留原始值，不输出DV
    let new_gt = if dp == 0 && alt_base == "." {
        stats.a_count += 1;
        format!("0/0:{}", dp)
    } else if dp >= dphom && r <= tol {
        stats.a_count += 1;
        format!("0/0:{}", dp)
    } else if dp >= dphom && r >= 1.0 - tol {
        stats.b_count += 1;
        format!("1/1:{}", dp)
    } else if dp >= dphet && r >= 0.5 - tol && r <= 0.5 + tol {
        stats.h_count += 1;
        format!("0/1:{}", dp)
    } else {
        stats.n_count += 1;
        format!("./.:{}", dp)
    };

    new_gt
}

/// 核心处理逻辑：FORMAT列完全保留原始值，修复类型不匹配问题
pub fn process_vcf_line(
    line: &str,
    args: &Args,
    dp_re: &Regex,
) -> Result<Option<String>> {
    let parts: Vec<&str> = line.split('\t').collect();
    if parts.len() < 10 {
        return Ok(None); // 列数不足，跳过无效行
    }

    // 1. 基础字段提取（FORMAT列直接用原始值）
    let ref_base = parts[3];
    let alt_base = parts[4];
    let qual = parts[5].parse::<f64>().unwrap_or(0.0);
    let info = parts[7];          // 保留原始INFO
    let format_str = parts[8];    // 【关键】原始FORMAT列，全程不修改
    let info_dp = extract_info_dp(info, dp_re);

    // 前置过滤（多碱基、低QUAL等）
    if ref_base.len() > 1
        || alt_base.len() > 1
        || qual < args.minqual
        || ref_base == "N"
        || info_dp < args.mindp
    {
        return Ok(None);
    }

    // 2. 解析FORMAT（仅用于找DP/AD，不修改FORMAT本身）
    let format_map = parse_format_fields(format_str);
    if !format_map.contains_key("DP") {
        return Ok(None); // 无DP字段则跳过
    }

    // 3. 处理样本列（仅修改GT，FORMAT列不变）
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

    // 4. 群体/MAF过滤（逻辑不变）
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

    // ========== 将new_parts改为Vec<String> ==========
    let mut new_parts = Vec::with_capacity(parts.len());

    // 前8列：将&str转为String存入
    for part in &parts[0..8] {
        new_parts.push(part.to_string());
    }

    // FORMAT列：原始值转为String存入（保持原列）
    new_parts.push(format_str.to_string());

    // 样本列：直接扩展String类型的modified_samples（无类型冲突）
    new_parts.extend(modified_samples);

    // 拼接最终行
    Ok(Some(new_parts.join("\t")))
}

/// 生成过滤规则注释行（无MQ/DV相关）
pub fn generate_filter_comment(args: &Args) -> String {
    let dphom = args.dphom;
    let dphet = args.dphet;
    let tol = args.tol;

    format!(
        "##FilterRule=<ID=CustomFilter,Description=\"dphom={}, dphet={}, tol={}, minqual={}, mindp={}, minhomn={}, minpresent={}, minhomp={}, minmaf={}\">",
        dphom, dphet, tol, args.minqual, args.mindp, args.minhomn, args.minpresent, args.minhomp, args.minmaf
    )
}