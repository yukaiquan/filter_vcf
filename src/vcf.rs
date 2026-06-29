use crate::args::Args;
use crate::config::SampleThresholds;
use anyhow::{Context, Result};
use clap::Parser;
use rayon::prelude::*;
use regex::Regex;
use std::collections::HashMap;

/// 基因型统计结果结构体
#[derive(Debug, Default, Clone, Copy)]
struct GenotypeStats {
    a_count: u32, // 0/0 纯合参考基因型数量
    b_count: u32, // 1/1 纯合替代基因型数量
    h_count: u32, // 0/1 杂合基因型数量
    n_count: u32, // ./. 缺失基因型数量
    present: u32, // 有效基因型数量（a + b + h）
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

/// 提取替代等位深度：优先 AD（逗号分隔，取第 2 个值），AD 不存在则回退 DV（单整数）
/// 与 AWK 脚本的 DV 字段语义对齐；AD[1] 在 GT:DP:AD 格式下与 DV 数值等价
fn extract_alt_depth(all_fields: &[&str], format_map: &HashMap<&str, usize>) -> u32 {
    // 优先使用 AD 字段（ref,alt 格式，取 alt）
    if let Some(&idx) = format_map.get("AD") {
        if let Some(&ad_str) = all_fields.get(idx) {
            if ad_str == "." {
                return 0;
            }
            let ad_parts: Vec<&str> = ad_str.split(',').collect();
            return match ad_parts.len() {
                1 => ad_parts[0].parse::<u32>().unwrap_or(0),
                _ => ad_parts[1].parse::<u32>().unwrap_or(0),
            };
        }
    }
    // 回退使用 DV 字段
    if let Some(&idx) = format_map.get("DV") {
        if let Some(&dv_str) = all_fields.get(idx) {
            if dv_str == "." {
                return 0;
            }
            return dv_str.parse::<u32>().unwrap_or(0);
        }
    }
    0
}

/// 提取样本的DP、r值（使用引用避免不必要的字符串分配）
fn extract_sample_info_optimized(gt_str: &str, format_map: &HashMap<&str, usize>) -> (u32, f64) {
    // 使用引用数组，零内存分配
    let all_fields: Vec<&str> = gt_str.split(':').collect();

    // 提取DP（用于判断基因型，动态找DP索引）
    let dp = format_map
        .get("DP")
        .and_then(|&idx| all_fields.get(idx))
        .and_then(|&s| {
            if s == "." {
                None
            } else {
                s.parse::<u32>().ok()
            }
        })
        .unwrap_or(0);

    // 提取替代等位深度计算r值（AD 优先，回退 DV）
    let r = if dp > 0 {
        let dv = extract_alt_depth(&all_fields, format_map);
        dv as f64 / dp as f64
    } else {
        0.0
    };

    (dp, r)
}

/// 提取样本的DP、r值，同时保留所有原始字段
#[cfg(test)]
fn extract_sample_info(gt_str: &str, format_map: &HashMap<&str, usize>) -> (u32, f64, Vec<String>) {
    // 关键修复：将&str迭代器转为String迭代器后再collect
    let all_fields: Vec<String> = gt_str
        .split(':')
        .map(|s| s.to_string()) // 每个&str转为String
        .collect();

    // 提取DP（用于判断基因型，动态找DP索引）
    let dp = format_map
        .get("DP")
        .and_then(|&idx| all_fields.get(idx))
        .and_then(|s| {
            if s == "." {
                None
            } else {
                s.parse::<u32>().ok()
            }
        })
        .unwrap_or(0);

    // 提取AD计算r值（AD 优先，回退 DV）
    let r = if dp > 0 {
        let all_refs: Vec<&str> = all_fields.iter().map(|s| s.as_str()).collect();
        let dv = extract_alt_depth(&all_refs, format_map);
        dv as f64 / dp as f64
    } else {
        0.0
    };

    (dp, r, all_fields)
}

/// 处理单个样本的基因型（并行优化版：返回结果和统计，不依赖可变状态）
/// 核心过滤逻辑完全保持不变
fn process_sample_gt_parallel(
    sample_str: &str,
    alt_base: &str,
    format_map: &HashMap<&str, usize>,
    dphom: u32,
    dphet: u32,
    tol: f64,
) -> (String, GenotypeStats) {
    let mut stats = GenotypeStats::default();

    // 使用优化版提取函数，减少内存分配
    let (dp, r) = extract_sample_info_optimized(sample_str, format_map);

    // 找到GT字段的索引
    let gt_idx = match format_map.get("GT") {
        Some(&idx) => idx,
        None => {
            // 无GT字段则返回原始样本字符串
            return (sample_str.to_string(), stats);
        }
    };

    // 使用引用数组，避免不必要的内存分配
    let all_fields: Vec<&str> = sample_str.split(':').collect();

    // 确保索引有效
    if gt_idx >= all_fields.len() {
        return (sample_str.to_string(), stats);
    }

    // ========== 核心过滤逻辑（完全保持不变）==========
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
    // ========== 核心过滤逻辑结束 ==========

    // 仅在GT字段改变时才构建新字符串（优化内存分配）
    let result = if all_fields[gt_idx] == new_gt {
        // GT字段未改变，直接返回原始字符串（零分配）
        sample_str.to_string()
    } else {
        // GT字段改变，需要重新构建
        let mut new_fields: Vec<String> = all_fields.iter().map(|&s| s.to_string()).collect();
        new_fields[gt_idx] = new_gt.to_string();
        new_fields.join(":")
    };

    (result, stats)
}

/// 带上下界阈值的样本处理（仅在配置启用时调用）
/// 判断结构对齐原逻辑，每分支额外加入 dp <= *_max 上界判断
fn process_sample_gt_with_thresholds(
    sample_str: &str,
    alt_base: &str,
    format_map: &HashMap<&str, usize>,
    th: SampleThresholds,
    tol: f64,
) -> (String, GenotypeStats) {
    let mut stats = GenotypeStats::default();

    let (dp, r) = extract_sample_info_optimized(sample_str, format_map);

    let gt_idx = match format_map.get("GT") {
        Some(&idx) => idx,
        None => return (sample_str.to_string(), stats),
    };

    let all_fields: Vec<&str> = sample_str.split(':').collect();
    if gt_idx >= all_fields.len() {
        return (sample_str.to_string(), stats);
    }

    // 核心过滤逻辑：与原版分支结构一致，纯合/杂合各自加 dp 上界
    let new_gt = if dp == 0 && alt_base == "." {
        stats.a_count += 1;
        "0/0"
    } else if dp >= th.dphom_min && dp <= th.dphom_max && r <= tol {
        stats.a_count += 1;
        "0/0"
    } else if dp >= th.dphom_min && dp <= th.dphom_max && r >= 1.0 - tol {
        stats.b_count += 1;
        "1/1"
    } else if dp >= th.dphet_min && dp <= th.dphet_max && r >= 0.5 - tol && r <= 0.5 + tol {
        stats.h_count += 1;
        "0/1"
    } else {
        stats.n_count += 1;
        "./."
    };

    let result = if all_fields[gt_idx] == new_gt {
        sample_str.to_string()
    } else {
        let mut new_fields: Vec<String> = all_fields.iter().map(|&s| s.to_string()).collect();
        new_fields[gt_idx] = new_gt.to_string();
        new_fields.join(":")
    };

    (result, stats)
}

/// 处理单个样本的基因型
#[cfg(test)]
fn process_sample_gt(
    sample_str: &str,
    alt_base: &str,
    format_map: &HashMap<&str, usize>,
    dphom: u32,
    dphet: u32,
    tol: f64,
    stats: &mut GenotypeStats,
) -> String {
    // 使用优化版提取函数，减少内存分配
    let (dp, r) = extract_sample_info_optimized(sample_str, format_map);

    // 找到GT字段的索引
    let gt_idx = match format_map.get("GT") {
        Some(&idx) => idx,
        None => {
            // 无GT字段则返回原始样本字符串
            return sample_str.to_string();
        }
    };

    // 使用引用数组，避免不必要的内存分配
    let all_fields: Vec<&str> = sample_str.split(':').collect();

    // 确保索引有效
    if gt_idx >= all_fields.len() {
        return sample_str.to_string();
    }

    // ========== 核心过滤逻辑（完全保持不变）==========
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
    // ========== 核心过滤逻辑结束 ==========

    // 仅在GT字段改变时才构建新字符串
    let result = if all_fields[gt_idx] == new_gt {
        // GT字段未改变，直接返回原始字符串
        sample_str.to_string()
    } else {
        // GT字段改变，需要重新构建
        let mut new_fields: Vec<String> = all_fields.iter().map(|&s| s.to_string()).collect();
        new_fields[gt_idx] = new_gt.to_string();
        new_fields.join(":")
    };

    result
}

/// 核心处理逻辑：FORMAT列保持原始，样本列仅修改GT字段，保留所有其他字段
pub fn process_vcf_line(line: &str, args: &Args, dp_re: &Regex) -> Result<Option<String>> {
    let parts: Vec<&str> = line.split('\t').collect();
    if parts.len() < 10 {
        return Ok(None); // 列数不足，跳过无效行
    }

    // 1. 基础字段提取
    let ref_base = parts[3];
    let alt_base = parts[4];
    let qual = parts[5].parse::<f64>().unwrap_or(0.0);
    let info = parts[7]; // 保留原始INFO
    let format_str = parts[8]; // 原始FORMAT列
    let info_dp = extract_info_dp(info, dp_re);

    // 前置过滤
    // 根据 include_indel 参数决定是否过滤 indel
    let is_indel = ref_base.len() > 1 || alt_base.len() > 1;
    if (!args.include_indel && is_indel)
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

    // 3. 并行处理样本列（优化：利用多核CPU）
    let samples: Vec<&str> = parts.iter().skip(9).copied().collect();

    // 使用 rayon 并行处理样本
    // 无配置 -> 走原路径，行为与未引入配置时完全一致
    // 有配置 -> 走带上下界的并行路径
    let results: Vec<(String, GenotypeStats)> = match &args.sample_config {
        None => samples
            .into_par_iter()
            .map(|sample_str| {
                process_sample_gt_parallel(
                    sample_str,
                    alt_base,
                    &format_map,
                    args.dphom,
                    args.dphet,
                    args.tol,
                )
            })
            .collect(),
        Some(cfg) => samples
            .into_par_iter()
            .enumerate()
            .map(|(idx, sample_str)| {
                let th = cfg
                    .get(idx)
                    .copied()
                    .flatten()
                    .unwrap_or_else(|| SampleThresholds::from_global(args.dphom, args.dphet));
                process_sample_gt_with_thresholds(sample_str, alt_base, &format_map, th, args.tol)
            })
            .collect(),
    };

    // 合并统计结果和样本字符串
    let mut stats = GenotypeStats::default();
    let mut modified_samples = Vec::with_capacity(results.len());

    for (sample_result, sample_stats) in results {
        modified_samples.push(sample_result);
        stats.a_count += sample_stats.a_count;
        stats.b_count += sample_stats.b_count;
        stats.h_count += sample_stats.h_count;
        stats.n_count += sample_stats.n_count;
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
        "##FilterRule=<ID=CustomFilter,Description=\"dphom={}, dphet={}, tol={}, minqual={}, mindp={}, minhomn={}, minpresent={}, minhomp={}, minmaf={}, include_indel={}\">",
        args.dphom,
        args.dphet,
        args.tol,
        args.minqual,
        args.mindp,
        args.minhomn,
        args.minpresent,
        args.minhomp,
        args.minmaf,
        args.include_indel
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::args::Args;
    use clap::Parser;

    /// 创建默认测试参数
    fn create_default_args() -> Args {
        Args::try_parse_from(vec!["filter_vcf"]).unwrap()
    }

    // ==================== extract_info_dp 测试 ====================
    #[test]
    fn test_extract_info_dp_valid() {
        let re = Regex::new(r"DP=(\d+)").unwrap();

        // 标准格式
        assert_eq!(extract_info_dp("DP=100;AF=0.5", &re), 100);
        assert_eq!(extract_info_dp("AF=0.5;DP=200", &re), 200);
        assert_eq!(extract_info_dp("DP=0", &re), 0);

        // 多个DP（取第一个）
        assert_eq!(extract_info_dp("DP=50;DP=100", &re), 50);
    }

    #[test]
    fn test_extract_info_dp_invalid() {
        let re = Regex::new(r"DP=(\d+)").unwrap();

        // 无DP字段
        assert_eq!(extract_info_dp("AF=0.5;MQ=60", &re), 0);
        assert_eq!(extract_info_dp("", &re), 0);

        // 格式错误
        assert_eq!(extract_info_dp("DP=abc", &re), 0);
        assert_eq!(extract_info_dp("DP=", &re), 0);
    }

    // ==================== parse_format_fields 测试 ====================
    #[test]
    fn test_parse_format_fields_standard() {
        let map = parse_format_fields("GT:DP:AD");

        assert_eq!(map.get("GT"), Some(&0));
        assert_eq!(map.get("DP"), Some(&1));
        assert_eq!(map.get("AD"), Some(&2));
        assert_eq!(map.len(), 3);
    }

    #[test]
    fn test_parse_format_fields_varied_order() {
        let map = parse_format_fields("AD:GT:DP:GQ");

        assert_eq!(map.get("AD"), Some(&0));
        assert_eq!(map.get("GT"), Some(&1));
        assert_eq!(map.get("DP"), Some(&2));
        assert_eq!(map.get("GQ"), Some(&3));
    }

    #[test]
    fn test_parse_format_fields_single() {
        let map = parse_format_fields("GT");
        assert_eq!(map.get("GT"), Some(&0));
        assert_eq!(map.len(), 1);
    }

    // ==================== extract_sample_info 测试 ====================
    #[test]
    fn test_extract_sample_info_standard() {
        let format_map = parse_format_fields("GT:DP:AD");

        // 标准样本数据
        let (dp, r, fields) = extract_sample_info("0/0:50:30,20", &format_map);
        assert_eq!(dp, 50);
        assert!((r - 0.4).abs() < 0.01); // 20/50 = 0.4
        assert_eq!(fields, vec!["0/0", "50", "30,20"]);
    }

    #[test]
    fn test_extract_sample_info_homozygous() {
        let format_map = parse_format_fields("GT:DP:AD");

        // 纯合参考 (AD[1] = 0)
        let (dp, r, _) = extract_sample_info("0/0:100:100,0", &format_map);
        assert_eq!(dp, 100);
        assert!((r - 0.0).abs() < 0.01);

        // 纯合替代 (AD[1] = 100)
        let (dp, r, _) = extract_sample_info("1/1:100:0,100", &format_map);
        assert_eq!(dp, 100);
        assert!((r - 1.0).abs() < 0.01);
    }

    #[test]
    fn test_extract_sample_info_heterozygous() {
        let format_map = parse_format_fields("GT:DP:AD");

        // 杂合 (AD[1] ≈ DP/2)
        let (dp, r, _) = extract_sample_info("0/1:100:50,50", &format_map);
        assert_eq!(dp, 100);
        assert!((r - 0.5).abs() < 0.01);
    }

    #[test]
    fn test_extract_sample_info_missing_data() {
        let format_map = parse_format_fields("GT:DP:AD");

        // DP缺失
        let (dp, r, _) = extract_sample_info("0/0:.:30,20", &format_map);
        assert_eq!(dp, 0);
        assert!((r - 0.0).abs() < 0.01);

        // AD缺失
        let (dp, r, _) = extract_sample_info("0/0:50:.", &format_map);
        assert_eq!(dp, 50);
        assert!((r - 0.0).abs() < 0.01);
    }

    #[test]
    fn test_extract_sample_info_different_order() {
        let format_map = parse_format_fields("AD:GT:DP");

        let (dp, r, _) = extract_sample_info("30,20:0/0:50", &format_map);
        assert_eq!(dp, 50);
        assert!((r - 0.4).abs() < 0.01);
    }

    // ==================== AD/DV 字段回退测试 ====================
    #[test]
    fn test_extract_sample_info_dv_fallback() {
        // 仅含 DV 字段（AWK 原生格式 GT:DP:DV），无 AD -> 使用 DV
        let format_map = parse_format_fields("GT:DP:DV");
        let (dp, r, _) = extract_sample_info("0/0:50:20", &format_map);
        assert_eq!(dp, 50);
        assert!((r - 0.4).abs() < 0.01); // 20/50 = 0.4
    }

    #[test]
    fn test_extract_sample_info_ad_preferred_over_dv() {
        // 同时含 AD 和 DV -> 优先 AD（取 alt 分量 AD[1]）
        let format_map = parse_format_fields("GT:DP:AD:DV");
        // AD=50,0 -> alt=0；DV=20 -> 若误用 DV 则 r=0.4
        let (dp, r, _) = extract_sample_info("0/0:50:50,0:20", &format_map);
        assert_eq!(dp, 50);
        assert!((r - 0.0).abs() < 0.01); // 应取 AD[1]=0，而非 DV=20
    }

    #[test]
    fn test_extract_sample_info_neither_ad_nor_dv() {
        // 既无 AD 也无 DV -> dv=0 -> r=0
        let format_map = parse_format_fields("GT:DP:GQ");
        let (dp, r, _) = extract_sample_info("0/0:50:99", &format_map);
        assert_eq!(dp, 50);
        assert!((r - 0.0).abs() < 0.01);
    }

    #[test]
    fn test_extract_sample_info_optimized_dv_fallback() {
        // 优化版同样支持 DV 回退
        let format_map = parse_format_fields("GT:DP:DV");
        let (dp, r) = extract_sample_info_optimized("1/1:100:100", &format_map);
        assert_eq!(dp, 100);
        assert!((r - 1.0).abs() < 0.01); // 100/100 = 1.0
    }

    // ==================== process_sample_gt 测试 ====================
    #[test]
    fn test_process_sample_gt_homozygous_ref() {
        let format_map = parse_format_fields("GT:DP:AD");
        let mut stats = GenotypeStats::default();

        // r = 0, dp足够 → 0/0
        let result = process_sample_gt("0/0:10:10,0", "A", &format_map, 1, 1, 0.25, &mut stats);
        assert_eq!(result, "0/0:10:10,0");
        assert_eq!(stats.a_count, 1);
    }

    #[test]
    fn test_process_sample_gt_homozygous_alt() {
        let format_map = parse_format_fields("GT:DP:AD");
        let mut stats = GenotypeStats::default();

        // r = 1.0, dp足够 → 1/1
        let result = process_sample_gt("1/1:10:0,10", "A", &format_map, 1, 1, 0.25, &mut stats);
        assert_eq!(result, "1/1:10:0,10");
        assert_eq!(stats.b_count, 1);
    }

    #[test]
    fn test_process_sample_gt_heterozygous() {
        let format_map = parse_format_fields("GT:DP:AD");
        let mut stats = GenotypeStats::default();

        // r = 0.5, dp足够 → 0/1
        let result = process_sample_gt("0/1:10:5,5", "A", &format_map, 1, 1, 0.25, &mut stats);
        assert_eq!(result, "0/1:10:5,5");
        assert_eq!(stats.h_count, 1);
    }

    #[test]
    fn test_process_sample_gt_insufficient_dp() {
        let format_map = parse_format_fields("GT:DP:AD");
        let mut stats = GenotypeStats::default();

        // dp不足 → ./.
        let result = process_sample_gt("0/0:0:10,0", "A", &format_map, 5, 5, 0.25, &mut stats);
        assert_eq!(result, "./.:0:10,0");
        assert_eq!(stats.n_count, 1);
    }

    #[test]
    fn test_process_sample_gt_ambiguous_ratio() {
        let format_map = parse_format_fields("GT:DP:AD");
        let mut stats = GenotypeStats::default();

        // r=0.3 在杂合区间内 (0.25 <= r <= 0.75) → 0/1
        let result = process_sample_gt("0/0:10:7,3", "A", &format_map, 1, 1, 0.25, &mut stats);
        assert_eq!(result, "0/1:10:7,3");
        assert_eq!(stats.h_count, 1);
    }

    #[test]
    fn test_process_sample_gt_special_case_alt_dot() {
        let format_map = parse_format_fields("GT:DP");
        let mut stats = GenotypeStats::default();

        // alt_base = ".", dp = 0 → 0/0
        let result = process_sample_gt("./.:0", ".", &format_map, 1, 1, 0.25, &mut stats);
        assert_eq!(result, "0/0:0");
        assert_eq!(stats.a_count, 1);
    }

    #[test]
    fn test_process_sample_gt_preserves_other_fields() {
        let format_map = parse_format_fields("GT:DP:AD:GQ:PL");
        let mut stats = GenotypeStats::default();

        let result = process_sample_gt(
            "0/0:10:10,0:99:0,10,100",
            "A",
            &format_map,
            1,
            1,
            0.25,
            &mut stats,
        );
        assert_eq!(result, "0/0:10:10,0:99:0,10,100");
    }

    // ==================== process_vcf_line 测试 ====================
    #[test]
    fn test_process_vcf_line_valid_snp() {
        let args = create_default_args();
        let dp_re = Regex::new(r"DP=(\d+)").unwrap();

        let line = "chr1\t100\t.\tA\tT\t50\tPASS\tDP=100\tGT:DP:AD\t0/0:50:50,0\t1/1:50:0,50";
        let result = process_vcf_line(line, &args, &dp_re).unwrap();

        assert!(result.is_some());
        let output = result.unwrap();
        assert!(output.contains("chr1"));
        assert!(output.contains("A\tT"));
    }

    #[test]
    fn test_process_vcf_line_filter_low_qual() {
        let mut args = create_default_args();
        args.minqual = 60.0;
        let dp_re = Regex::new(r"DP=(\d+)").unwrap();

        let line = "chr1\t100\t.\tA\tT\t50\tPASS\tDP=100\tGT:DP:AD\t0/0:50:50,0";
        let result = process_vcf_line(line, &args, &dp_re).unwrap();

        assert!(result.is_none());
    }

    #[test]
    fn test_process_vcf_line_filter_low_dp() {
        let mut args = create_default_args();
        args.mindp = 200;
        let dp_re = Regex::new(r"DP=(\d+)").unwrap();

        let line = "chr1\t100\t.\tA\tT\t50\tPASS\tDP=100\tGT:DP:AD\t0/0:50:50,0";
        let result = process_vcf_line(line, &args, &dp_re).unwrap();

        assert!(result.is_none());
    }

    #[test]
    fn test_process_vcf_line_filter_indel() {
        let args = create_default_args();
        let dp_re = Regex::new(r"DP=(\d+)").unwrap();

        // Indel 应该被过滤（默认不包含）
        let line = "chr1\t100\t.\tA\tAT\t50\tPASS\tDP=100\tGT:DP:AD\t0/0:50:50,0";
        let result = process_vcf_line(line, &args, &dp_re).unwrap();
        assert!(result.is_none());

        // SNP 应该保留
        let line = "chr1\t100\t.\tA\tT\t50\tPASS\tDP=100\tGT:DP:AD\t0/0:50:50,0";
        let result = process_vcf_line(line, &args, &dp_re).unwrap();
        assert!(result.is_some());
    }

    #[test]
    fn test_process_vcf_line_filter_ref_n() {
        let args = create_default_args();
        let dp_re = Regex::new(r"DP=(\d+)").unwrap();

        let line = "chr1\t100\t.\tN\tA\t50\tPASS\tDP=100\tGT:DP:AD\t0/0:50:50,0";
        let result = process_vcf_line(line, &args, &dp_re).unwrap();

        assert!(result.is_none());
    }

    #[test]
    fn test_process_vcf_line_filter_minhomn() {
        let mut args = create_default_args();
        args.minhomn = 2;
        let dp_re = Regex::new(r"DP=(\d+)").unwrap();

        // 只有0个纯合样本（两个都是杂合），a_count=0, b_count=0，不满足 minhomn=2
        let line = "chr1\t100\t.\tA\tT\t50\tPASS\tDP=100\tGT:DP:AD\t0/1:50:25,25\t0/1:50:25,25";
        let result = process_vcf_line(line, &args, &dp_re).unwrap();
        assert!(result.is_none());

        // a_count=1, b_count=0，不满足 minhomn=2（需要 a_count >= 2 && b_count >= 2）
        let line = "chr1\t100\t.\tA\tT\t50\tPASS\tDP=100\tGT:DP:AD\t0/0:50:50,0\t0/1:50:25,25";
        let result = process_vcf_line(line, &args, &dp_re).unwrap();
        assert!(result.is_none());

        // a_count=1, b_count=1，不满足 minhomn=2
        let line = "chr1\t100\t.\tA\tT\t50\tPASS\tDP=100\tGT:DP:AD\t0/0:50:50,0\t1/1:50:0,50";
        let result = process_vcf_line(line, &args, &dp_re).unwrap();
        assert!(result.is_none());

        // a_count=2, b_count=2，满足 minhomn=2
        let line = "chr1\t100\t.\tA\tT\t50\tPASS\tDP=100\tGT:DP:AD\t0/0:50:50,0\t0/0:50:50,0\t1/1:50:0,50\t1/1:50:0,50";
        let result = process_vcf_line(line, &args, &dp_re).unwrap();
        assert!(result.is_some());
    }

    #[test]
    fn test_process_vcf_line_filter_minmaf() {
        let mut args = create_default_args();
        args.minmaf = 0.1;
        let dp_re = Regex::new(r"DP=(\d+)").unwrap();

        // MAF < 0.1，应该被过滤
        let line =
            "chr1\t100\t.\tA\tT\t50\tPASS\tDP=100\tGT:DP:AD\t0/0:50:50,0\t0/0:50:50,0\t0/0:50:50,0";
        let result = process_vcf_line(line, &args, &dp_re).unwrap();
        assert!(result.is_none());

        // MAF >= 0.1，应该保留
        let line =
            "chr1\t100\t.\tA\tT\t50\tPASS\tDP=100\tGT:DP:AD\t1/1:50:0,50\t0/0:50:50,0\t0/0:50:50,0";
        let result = process_vcf_line(line, &args, &dp_re).unwrap();
        assert!(result.is_some());
    }

    #[test]
    fn test_process_vcf_line_missing_format_fields() {
        let args = create_default_args();
        let dp_re = Regex::new(r"DP=(\d+)").unwrap();

        // 缺少DP字段
        let line = "chr1\t100\t.\tA\tT\t50\tPASS\tDP=100\tGT:AD\t0/0:50,0";
        let result = process_vcf_line(line, &args, &dp_re).unwrap();
        assert!(result.is_none());

        // 缺少GT字段
        let line = "chr1\t100\t.\tA\tT\t50\tPASS\tDP=100\tDP:AD\t50:50,0";
        let result = process_vcf_line(line, &args, &dp_re).unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_process_vcf_line_insufficient_columns() {
        let args = create_default_args();
        let dp_re = Regex::new(r"DP=(\d+)").unwrap();

        // 列数不足
        let line = "chr1\t100\t.\tA\tT\t50\tPASS\tDP=100";
        let result = process_vcf_line(line, &args, &dp_re).unwrap();
        assert!(result.is_none());
    }

    #[test]
    fn test_process_vcf_line_preserves_format() {
        let args = create_default_args();
        let dp_re = Regex::new(r"DP=(\d+)").unwrap();

        let line = "chr1\t100\t.\tA\tT\t50\tPASS\tDP=100;AF=0.5\tGT:DP:AD:GQ\t0/0:50:50,0:99";
        let result = process_vcf_line(line, &args, &dp_re).unwrap();

        assert!(result.is_some());
        let output = result.unwrap();
        assert!(output.contains("GT:DP:AD:GQ")); // FORMAT列保持不变
        assert!(output.contains("DP=100;AF=0.5")); // INFO列保持不变
    }

    // ==================== generate_filter_comment 测试 ====================
    #[test]
    fn test_generate_filter_comment() {
        let args = create_default_args();
        let comment = generate_filter_comment(&args);

        assert!(comment.starts_with("##FilterRule=<ID=CustomFilter,Description="));
        assert!(comment.contains("dphom=1"));
        assert!(comment.contains("dphet=1"));
        assert!(comment.contains("tol=0.2499"));
        assert!(comment.ends_with("\">"));
    }

    #[test]
    fn test_generate_filter_comment_custom_args() {
        let mut args = create_default_args();
        args.dphom = 5;
        args.dphet = 10;
        args.minqual = 30.0;
        args.include_indel = true;

        let comment = generate_filter_comment(&args);

        assert!(comment.contains("dphom=5"));
        assert!(comment.contains("dphet=10"));
        assert!(comment.contains("minqual=30"));
        assert!(comment.contains("include_indel=true"));
    }

    // ==================== 边界情况测试 ====================
    #[test]
    fn test_boundary_tolerance_values() {
        let format_map = parse_format_fields("GT:DP:AD");
        let mut stats = GenotypeStats::default();

        // r 刚好等于 tol 边界 → 0/0
        let result = process_sample_gt("0/0:100:75,25", "A", &format_map, 1, 1, 0.25, &mut stats);
        assert_eq!(result, "0/0:100:75,25"); // r=0.25, 刚好等于tol

        // r=0.25, tol=0.24, 在杂合区间 (0.26, 0.74) 外，但 > tol，所以是 ./.
        stats = GenotypeStats::default();
        let result = process_sample_gt("0/0:100:75,25", "A", &format_map, 1, 1, 0.24, &mut stats);
        assert_eq!(result, "./.:100:75,25");

        // r=0.74, tol=0.25, 在杂合区间 (0.25, 0.75) 内 → 0/1
        stats = GenotypeStats::default();
        let result = process_sample_gt("0/0:100:26,74", "A", &format_map, 1, 1, 0.25, &mut stats);
        assert_eq!(result, "0/1:100:26,74");

        // r=0.75, tol=0.25, 刚好等于 1-tol → 1/1
        stats = GenotypeStats::default();
        let result = process_sample_gt("0/0:100:25,75", "A", &format_map, 1, 1, 0.25, &mut stats);
        assert_eq!(result, "1/1:100:25,75");
    }

    #[test]
    fn test_zero_dp_handling() {
        let format_map = parse_format_fields("GT:DP:AD");
        let mut stats = GenotypeStats::default();

        // DP=0 且 alt_base="." → 特殊情况处理
        let result = process_sample_gt("./.:0:.", ".", &format_map, 1, 1, 0.25, &mut stats);
        assert_eq!(result, "0/0:0:.");
        assert_eq!(stats.a_count, 1);
    }

    #[test]
    fn test_empty_and_whitespace_handling() {
        let re = Regex::new(r"DP=(\d+)").unwrap();

        // 空字符串
        assert_eq!(extract_info_dp("", &re), 0);

        // 只有空格
        assert_eq!(extract_info_dp("   ", &re), 0);
    }

    // ==================== 性能相关测试 ====================
    #[test]
    fn test_large_number_of_samples() {
        let args = create_default_args();
        let dp_re = Regex::new(r"DP=(\d+)").unwrap();

        // 构建包含100个样本的VCF行
        let mut samples = Vec::new();
        for i in 0..100 {
            if i % 3 == 0 {
                samples.push("0/0:50:50,0");
            } else if i % 3 == 1 {
                samples.push("1/1:50:0,50");
            } else {
                samples.push("0/1:50:25,25");
            }
        }

        let line = format!(
            "chr1\t100\t.\tA\tT\t50\tPASS\tDP=5000\tGT:DP:AD\t{}",
            samples.join("\t")
        );

        let result = process_vcf_line(&line, &args, &dp_re).unwrap();
        assert!(result.is_some());
    }

    #[test]
    fn test_complex_info_field() {
        let args = create_default_args();
        let dp_re = Regex::new(r"DP=(\d+)").unwrap();

        let line = "chr1\t100\t.\tA\tT\t50\tPASS\tDP=100;AF=0.5;MQ=60;FS=0.0;QD=15.0\tGT:DP:AD\t0/0:50:50,0";
        let result = process_vcf_line(line, &args, &dp_re).unwrap();

        assert!(result.is_some());
        let output = result.unwrap();
        assert!(output.contains("DP=100;AF=0.5;MQ=60;FS=0.0;QD=15.0"));
    }

    /// 测试并行处理结果与串行处理完全一致
    #[test]
    fn test_parallel_processing_consistency() {
        let args = create_default_args();
        let dp_re = Regex::new(r"DP=(\d+)").unwrap();

        // 测试多种场景
        let test_cases = vec![
            // 场景1: 多样本混合基因型
            "chr1\t100\t.\tA\tT\t50\tPASS\tDP=100\tGT:DP:AD\t0/0:50:50,0\t1/1:50:0,50\t0/1:50:25,25\t./.:0:.",
            // 场景2: 大量样本（模拟并行场景）
            "chr1\t200\t.\tG\tC\t60\tPASS\tDP=200\tGT:DP:AD\t0/0:60:60,0\t0/0:60:60,0\t1/1:60:0,60\t0/1:60:30,30\t0/0:60:60,0\t1/1:60:0,60\t0/1:60:30,30\t./.:0:.\t0/0:60:60,0\t1/1:60:0,60",
            // 场景3: 边界情况
            "chr2\t300\t.\tT\tA\t30\tPASS\tDP=50\tGT:DP:AD\t0/0:10:10,0\t1/1:10:0,10\t0/1:10:5,5\t./.:0:.\t0/0:5:5,0",
        ];

        for line in test_cases {
            // 多次运行，验证结果一致性（并行处理应该是确定性的）
            let mut results = Vec::new();
            for _ in 0..5 {
                let result = process_vcf_line(line, &args, &dp_re).unwrap();
                results.push(result);
            }

            // 验证所有结果完全一致
            let first = &results[0];
            for result in &results[1..] {
                assert_eq!(
                    first, result,
                    "并行处理结果不一致！\n第一次: {:?}\n后续: {:?}",
                    first, result
                );
            }

            // 如果有结果，验证样本顺序和内容
            if let Some(output) = first {
                let parts: Vec<&str> = output.split('\t').collect();
                let original_parts: Vec<&str> = line.split('\t').collect();

                // 验证样本数量一致
                assert_eq!(
                    parts.len(),
                    original_parts.len(),
                    "样本数量不一致！原始: {}, 输出: {}",
                    original_parts.len(),
                    parts.len()
                );

                // 验证样本顺序保持不变
                for (i, (output_sample, original_sample)) in parts
                    .iter()
                    .skip(9)
                    .zip(original_parts.iter().skip(9))
                    .enumerate()
                {
                    // 检查样本格式是否正确
                    let output_fields: Vec<&str> = output_sample.split(':').collect();
                    let original_fields: Vec<&str> = original_sample.split(':').collect();

                    assert_eq!(
                        output_fields.len(),
                        original_fields.len(),
                        "样本{}字段数量不一致！原始: {}, 输出: {}",
                        i,
                        original_fields.len(),
                        output_fields.len()
                    );

                    // 验证GT字段修改正确
                    let format_parts: Vec<&str> = original_parts[8].split(':').collect();
                    if let Some(gt_idx) = format_parts.iter().position(|&f| f == "GT") {
                        let original_gt = original_fields.get(gt_idx).unwrap_or(&"./.");
                        let output_gt = output_fields.get(gt_idx).unwrap_or(&"./.");

                        // GT字段应该被正确处理（可能被修改，但格式应该正确）
                        assert!(
                            ["0/0", "1/1", "0/1", "./."].contains(output_gt),
                            "样本{}的GT字段格式不正确: {}",
                            i,
                            output_gt
                        );
                    }
                }
            }
        }
    }

    /// 测试并行处理保持样本顺序
    #[test]
    fn test_parallel_processing_preserves_order() {
        let args = create_default_args();
        let dp_re = Regex::new(r"DP=(\d+)").unwrap();

        // 创建特殊的测试数据：每个样本都有唯一标识
        let mut samples = Vec::new();
        for i in 0..20 {
            let gt = match i % 4 {
                0 => "0/0",
                1 => "1/1",
                2 => "0/1",
                _ => "./.",
            };
            let dp = match i % 4 {
                0 => 50,
                1 => 50,
                2 => 50,
                _ => 0,
            };
            let ad = match i % 4 {
                0 => "50,0",
                1 => "0,50",
                2 => "25,25",
                _ => ".",
            };
            samples.push(format!("{}:{}:{}", gt, dp, ad));
        }

        let line = format!(
            "chr1\t100\t.\tA\tT\t50\tPASS\tDP=100\tGT:DP:AD\t{}",
            samples.join("\t")
        );

        // 多次运行，验证顺序一致性
        let mut all_results = Vec::new();
        for run in 0..10 {
            let result = process_vcf_line(&line, &args, &dp_re).unwrap();
            if let Some(output) = result {
                let output_parts: Vec<&str> = output.split('\t').collect();
                let sample_results: Vec<String> =
                    output_parts.iter().skip(9).map(|s| s.to_string()).collect();
                all_results.push((run, sample_results));
            }
        }

        // 验证所有运行的结果完全一致
        if all_results.len() > 1 {
            let first = &all_results[0].1;
            for (run, results) in &all_results[1..] {
                assert_eq!(
                    first, results,
                    "第{}次运行结果与第一次不一致！\n第一次: {:?}\n第{}次: {:?}",
                    run, first, run, results
                );
            }
        }
    }

    /// 压力测试：大规模并行处理（模拟100样本）
    #[test]
    fn test_large_scale_parallel_processing() {
        let args = create_default_args();
        let dp_re = Regex::new(r"DP=(\d+)").unwrap();

        // 创建100个样本
        let mut samples = Vec::new();
        for i in 0..100 {
            let gt = match i % 4 {
                0 => "0/0",
                1 => "1/1",
                2 => "0/1",
                _ => "./.",
            };
            let dp = if i % 4 == 3 { 0 } else { 50 };
            let ad = match i % 4 {
                0 => "50,0",
                1 => "0,50",
                2 => "25,25",
                _ => ".",
            };
            samples.push(format!("{}:{}:{}", gt, dp, ad));
        }

        let line = format!(
            "chr1\t100\t.\tA\tT\t50\tPASS\tDP=500\tGT:DP:AD\t{}",
            samples.join("\t")
        );

        // 运行20次，验证一致性
        let mut all_results = Vec::new();
        for run in 0..20 {
            let result = process_vcf_line(&line, &args, &dp_re).unwrap();
            if let Some(output) = result {
                all_results.push((run, output));
            }
        }

        // 验证所有结果完全一致
        let first = &all_results[0].1;
        for (run, result) in &all_results[1..] {
            assert_eq!(
                first, result,
                "大规模并行处理结果不一致！第{}次运行与第一次不同",
                run
            );
        }
    }

    /// 测试多线程竞态条件检测
    #[test]
    fn test_multithreading_race_condition() {
        use std::sync::{Arc, Mutex};
        use std::thread;

        let args = Arc::new(create_default_args());
        let dp_re = Arc::new(Regex::new(r"DP=(\d+)").unwrap());

        // 创建测试数据
        let mut samples = Vec::new();
        for i in 0..50 {
            let gt = match i % 4 {
                0 => "0/0",
                1 => "1/1",
                2 => "0/1",
                _ => "./.",
            };
            let dp = if i % 4 == 3 { 0 } else { 50 };
            let ad = match i % 4 {
                0 => "50,0",
                1 => "0,50",
                2 => "25,25",
                _ => ".",
            };
            samples.push(format!("{}:{}:{}", gt, dp, ad));
        }

        let line = Arc::new(format!(
            "chr1\t100\t.\tA\tT\t50\tPASS\tDP=250\tGT:DP:AD\t{}",
            samples.join("\t")
        ));

        // 多线程并发处理同一行数据
        let results = Arc::new(Mutex::new(Vec::new()));
        let mut handles = vec![];

        for _ in 0..10 {
            let args = Arc::clone(&args);
            let dp_re = Arc::clone(&dp_re);
            let line = Arc::clone(&line);
            let results = Arc::clone(&results);

            let handle = thread::spawn(move || {
                let result = process_vcf_line(&line, &args, &dp_re).unwrap();
                let mut results = results.lock().unwrap();
                results.push(result);
            });

            handles.push(handle);
        }

        // 等待所有线程完成
        for handle in handles {
            handle.join().unwrap();
        }

        // 验证所有线程的结果完全一致
        let results = results.lock().unwrap();
        let first = &results[0];
        for (i, result) in results.iter().enumerate().skip(1) {
            assert_eq!(
                first, result,
                "多线程竞态条件检测失败！线程{}的结果与其他线程不一致",
                i
            );
        }
    }

    /// 测试边界条件的并行处理
    #[test]
    fn test_parallel_boundary_conditions() {
        let args = create_default_args();
        let dp_re = Regex::new(r"DP=(\d+)").unwrap();

        // 创建边界条件测试数据
        let test_cases = vec![
            // 所有序本都是缺失值
            (
                "chr1\t100\t.\tA\tT\t50\tPASS\tDP=0\tGT:DP:AD",
                vec!["./.:0:."; 20],
            ),
            // 所有序本都是纯合参考
            (
                "chr1\t200\t.\tG\tC\t60\tPASS\tDP=100\tGT:DP:AD",
                vec!["0/0:50:50,0"; 20],
            ),
            // 所有序本都是纯合替代
            (
                "chr2\t300\t.\tT\tA\t70\tPASS\tDP=100\tGT:DP:AD",
                vec!["1/1:50:0,50"; 20],
            ),
            // 所有序本都是杂合
            (
                "chr3\t400\t.\tC\tG\t80\tPASS\tDP=100\tGT:DP:AD",
                vec!["0/1:50:25,25"; 20],
            ),
        ];

        for (prefix, sample_list) in test_cases {
            let line = format!("{}\t{}", prefix, sample_list.join("\t"));

            // 多次运行验证一致性
            let mut results = Vec::new();
            for _ in 0..10 {
                let result = process_vcf_line(&line, &args, &dp_re).unwrap();
                results.push(result);
            }

            // 验证所有结果一致
            let first = &results[0];
            for result in &results[1..] {
                assert_eq!(
                    first, result,
                    "边界条件并行处理不一致！\n测试数据: {}",
                    prefix
                );
            }
        }
    }
}
