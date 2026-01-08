use anyhow::{Context, Result};
use clap::Parser;
use regex::Regex;
use std::collections::HashMap;
use std::io::{ BufRead};
use crate::args::Args;

/// 基因型统计结果结构体（仅保留计数，不处理字段）
#[derive(Debug, Default)]
struct GenotypeStats {
    a_count: u32,    // 0/0 纯合参考基因型数量
    b_count: u32,    // 1/1 纯合替代基因型数量
    h_count: u32,    // 0/1 杂合基因型数量
    n_count: u32,    // ./. 缺失基因型数量
    present: u32,    // 有效基因型数量（a + b + h）
}

/// 从INFO字段中提取指定数值
fn extract_info_value(info: &str, pattern: &Regex) -> Result<u32> {
    pattern
        .captures(info)
        .and_then(|caps| caps.get(1))
        .map(|m| m.as_str().parse::<u32>())
        .unwrap_or(Ok(0))
        .context("Failed to extract value from INFO field")
}

/// 动态解析FORMAT字段，返回字段名到索引的映射
fn parse_format_fields(format_str: &str) -> HashMap<&str, usize> {
    let mut field_map = HashMap::new();
    for (idx, field) in format_str.split(':').enumerate() {
        field_map.insert(field, idx);
    }
    field_map
}

/// 从样本基因型字符串中提取DP和DV（仅用于过滤计算）
fn extract_dp_dv(gt_str: &str, format_map: &HashMap<&str, usize>) -> (u32, u32) {
    let gt_parts: Vec<&str> = gt_str.split(':').collect();

    // 提取DP（处理.的情况）
    let dp = format_map.get("DP")
        .and_then(|&idx| gt_parts.get(idx))
        .and_then(|s| if *s == "." { None } else { s.parse::<u32>().ok() })
        .unwrap_or(0);

    // 提取DV（从AD字段，仅用于计算r值）
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

    (dp, dv)
}

/// 核心过滤逻辑：仅判断是否符合条件，符合则返回原始行，否则返回None
pub fn process_vcf_line(
    line: &str,
    args: &Args,
    dp_re: &Regex,
) -> Result<Option<String>> {
    let parts: Vec<&str> = line.split('\t').collect();
    if parts.len() < 10 {
        return Ok(None); // 列数不足，跳过无效行
    }

    // 1. 基础过滤条件
    let ref_base = parts[3];
    let alt_base = parts[4];
    let qual = parts[5].parse::<f64>().unwrap_or(0.0);
    let info = parts[7];
    let format_str = parts[8];

    // 提取INFO中的总DP（用于过滤）
    let info_dp = extract_info_value(info, dp_re)?;

    // 基础过滤：多碱基变异、低质量、参考碱基为N、低总DP
    if ref_base.len() > 1
        || alt_base.len() > 1
        || qual < args.minqual
        || ref_base == "N"
        || info_dp < args.mindp
    {
        return Ok(None);
    }

    // 2. 解析FORMAT（必须包含DP字段）
    let format_map = parse_format_fields(format_str);
    if !format_map.contains_key("DP") {
        return Ok(None);
    }

    // 3. 统计样本基因型（仅计数，不修改字段）
    let mut stats = GenotypeStats::default();
    for i in 9..parts.len() {
        let gt_str = parts[i];
        let (dp, dv) = extract_dp_dv(gt_str, &format_map);
        let r = if dp > 0 { dv as f64 / dp as f64 } else { 0.0 };

        // 基因型分类计数（仅用于过滤，不修改原始GT）
        if dp == 0 && alt_base == "." {
            stats.a_count += 1;
        } else if dp >= args.dphom && r <= args.tol {
            stats.a_count += 1;
        } else if dp >= args.dphom && r >= 1.0 - args.tol {
            stats.b_count += 1;
        } else if dp >= args.dphet && r >= 0.5 - args.tol && r <= 0.5 + args.tol {
            stats.h_count += 1;
        } else {
            stats.n_count += 1;
        }
    }

    // 计算有效样本数
    stats.present = stats.a_count + stats.b_count + stats.h_count;

    // 4. 群体水平过滤
    let present_ratio = if stats.present + stats.n_count > 0 {
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

    // 5. MAF过滤
    let maf = if stats.b_count > stats.a_count {
        (2 * stats.a_count + stats.h_count) as f64 / (2 * stats.present) as f64
    } else {
        (2 * stats.b_count + stats.h_count) as f64 / (2 * stats.present) as f64
    };

    if maf < args.minmaf {
        return Ok(None);
    }

    // 所有过滤条件都满足，返回原始行
    Ok(Some(line.to_string()))
}

/// 生成过滤规则注释行（用于写入Header最后面）
pub fn generate_filter_comment(args: &Args) -> String {
    format!(
        "##FilterRule=<ID=CustomFilter,Description=\"dphom={}, dphet={}, tol={}, minqual={}, mindp={}, minhomn={}, minpresent={}, minhomp={}, minmaf={}\">",
        args.dphom, args.dphet, args.tol, args.minqual, args.mindp, args.minhomn, args.minpresent, args.minhomp, args.minmaf
    )
}
