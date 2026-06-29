use anyhow::{Context, Result, bail};
use std::collections::{HashMap, HashSet};
use std::fs;

/// 单个样本的深度阈值（纯合/杂合的下限与上限）
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SampleThresholds {
    pub dphom_min: u32,
    pub dphom_max: u32,
    pub dphet_min: u32,
    pub dphet_max: u32,
}

impl SampleThresholds {
    /// 用全局阈值构造默认阈值（上界为 u32::MAX，等价于无上界）
    /// 当配置未覆盖某样本时回退到此值
    pub fn from_global(dphom: u32, dphet: u32) -> Self {
        Self {
            dphom_min: dphom,
            dphom_max: u32::MAX,
            dphet_min: dphet,
            dphet_max: u32::MAX,
        }
    }
}

/// 按样本列顺序对齐的阈值表；下标对应该样本在 #CHROM 行中的列序号
/// 元素为 None 表示该样本未在配置中指定（回退到全局阈值）
pub type SampleConfig = Vec<Option<SampleThresholds>>;

/// 解析单个阈值单元格：`.` 或空表示不设置（None），否则解析为 u32
fn parse_cell(s: &str, field: &str, line_no: usize) -> Result<Option<u32>> {
    let s = s.trim();
    if s.is_empty() || s == "." {
        return Ok(None);
    }
    s.parse::<u32>().map(Some).with_context(|| {
        format!(
            "Invalid value for {} on line {} of config file: {}",
            field, line_no, s
        )
    })
}

/// 加载样本阈值配置文件，并按 VCF 的样本名顺序对齐为 SampleConfig
///
/// 配置文件格式（TSV，含表头）：
/// ```text
/// sample	dphom_min	dphom_max	dphet_min	dphet_max
/// SampleA	5	100	3	50
/// SampleB	.	.	.	.
/// ```
/// - 第 1 列为样本名（须与 #CHROM 行样本列一致）
/// - 其余 4 列为纯合下限/上限、杂合下限/上限；`.` 或留空表示不设置（回退全局）
/// - 未在配置中出现的样本回退到全局 --dphom/--dphet 阈值
pub fn load_config(path: &str, sample_names: &[String]) -> Result<SampleConfig> {
    let content = fs::read_to_string(path)
        .with_context(|| format!("Failed to read config file: {}", path))?;

    // 临时按样本名建索引
    let mut map: HashMap<String, SampleThresholds> = HashMap::new();
    let mut header_seen = false;

    for (line_idx, raw) in content.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        // 跳过表头：首个非空非注释行，且第一列为 "sample"
        if !header_seen
            && line
                .split('\t')
                .next()
                .map(|s| s.trim() == "sample")
                .unwrap_or(false)
        {
            header_seen = true;
            continue;
        }
        header_seen = true;

        let cols: Vec<&str> = line.split('\t').collect();
        if cols.len() < 5 {
            bail!(
                "Insufficient columns on line {} of config file (expected 5: sample,dphom_min,dphom_max,dphet_min,dphet_max)",
                line_idx + 1
            );
        }
        let sample = cols[0].trim().to_string();
        let dphom_min = parse_cell(cols[1], "dphom_min", line_idx + 1)?;
        let dphom_max = parse_cell(cols[2], "dphom_max", line_idx + 1)?;
        let dphet_min = parse_cell(cols[3], "dphet_min", line_idx + 1)?;
        let dphet_max = parse_cell(cols[4], "dphet_max", line_idx + 1)?;

        // 校验：下限 <= 上限
        if let (Some(lo), Some(hi)) = (dphom_min, dphom_max) {
            if lo > hi {
                bail!(
                    "dphom_min({}) > dphom_max({}) on line {} of config file",
                    lo,
                    hi,
                    line_idx + 1
                );
            }
        }
        if let (Some(lo), Some(hi)) = (dphet_min, dphet_max) {
            if lo > hi {
                bail!(
                    "dphet_min({}) > dphet_max({}) on line {} of config file",
                    lo,
                    hi,
                    line_idx + 1
                );
            }
        }

        map.insert(
            sample,
            SampleThresholds {
                dphom_min: dphom_min.unwrap_or(0),
                dphom_max: dphom_max.unwrap_or(u32::MAX),
                dphet_min: dphet_min.unwrap_or(0),
                dphet_max: dphet_max.unwrap_or(u32::MAX),
            },
        );
    }

    if map.is_empty() {
        bail!(
            "Config file {} contains no valid sample configuration rows",
            path
        );
    }

    // 按样本名顺序对齐
    let cfg: SampleConfig = sample_names
        .iter()
        .map(|name| map.get(name).copied())
        .collect();

    // 报告未匹配的配置项（仅警告，不影响运行）
    let matched: HashSet<&str> = sample_names.iter().map(|s| s.as_str()).collect();
    for name in map.keys() {
        if !matched.contains(name.as_str()) {
            eprintln!(
                "[warn] sample '{}' in config file not found in VCF, ignored",
                name
            );
        }
    }

    Ok(cfg)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_config(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        f
    }

    #[test]
    fn test_load_config_basic() {
        let f = write_config(
            "sample\tdphom_min\tdphom_max\tdphet_min\tdphet_max\n\
             SampleA\t5\t100\t3\t50\n\
             SampleB\t10\t200\t5\t80\n",
        );
        let names = vec!["SampleA".to_string(), "SampleB".to_string()];
        let cfg = load_config(f.path().to_str().unwrap(), &names).unwrap();

        assert_eq!(cfg.len(), 2);
        let a = cfg[0].unwrap();
        assert_eq!(
            a,
            SampleThresholds {
                dphom_min: 5,
                dphom_max: 100,
                dphet_min: 3,
                dphet_max: 50
            }
        );
        let b = cfg[1].unwrap();
        assert_eq!(b.dphom_min, 10);
        assert_eq!(b.dphom_max, 200);
    }

    #[test]
    fn test_load_config_partial_samples_fallback() {
        // 配置只覆盖部分样本，未覆盖的回退 None
        let f = write_config(
            "sample\tdphom_min\tdphom_max\tdphet_min\tdphet_max\n\
                             OnlyA\t5\t100\t3\t50\n",
        );
        let names = vec!["OnlyA".to_string(), "Missing".to_string()];
        let cfg = load_config(f.path().to_str().unwrap(), &names).unwrap();

        assert!(cfg[0].is_some());
        assert!(cfg[1].is_none()); // 回退全局
    }

    #[test]
    fn test_load_config_dot_means_unset() {
        // `.` 表示不设置：下限回退 0，上限回退 u32::MAX
        let f = write_config(
            "sample\tdphom_min\tdphom_max\tdphet_min\tdphet_max\n\
                             A\t.\t.\t.\t.\n",
        );
        let names = vec!["A".to_string()];
        let cfg = load_config(f.path().to_str().unwrap(), &names).unwrap();
        let a = cfg[0].unwrap();
        assert_eq!(a.dphom_min, 0);
        assert_eq!(a.dphom_max, u32::MAX);
        assert_eq!(a.dphet_min, 0);
        assert_eq!(a.dphet_max, u32::MAX);
    }

    #[test]
    fn test_load_config_partial_bounds() {
        // 只设下限，上限留空（用 . 显式占位）-> 回退 u32::MAX
        let f = write_config(
            "sample\tdphom_min\tdphom_max\tdphet_min\tdphet_max\n\
                             A\t5\t.\t3\t.\n",
        );
        let names = vec!["A".to_string()];
        let cfg = load_config(f.path().to_str().unwrap(), &names).unwrap();
        let a = cfg[0].unwrap();
        assert_eq!(a.dphom_min, 5);
        assert_eq!(a.dphom_max, u32::MAX);
        assert_eq!(a.dphet_min, 3);
        assert_eq!(a.dphet_max, u32::MAX);
    }

    #[test]
    fn test_load_config_invalid_lower_gt_upper() {
        let f = write_config(
            "sample\tdphom_min\tdphom_max\tdphet_min\tdphet_max\n\
                             A\t100\t5\t3\t50\n",
        );
        let names = vec!["A".to_string()];
        let result = load_config(f.path().to_str().unwrap(), &names);
        assert!(result.is_err());
    }

    #[test]
    fn test_load_config_invalid_value() {
        let f = write_config(
            "sample\tdphom_min\tdphom_max\tdphet_min\tdphet_max\n\
                             A\tabc\t5\t3\t50\n",
        );
        let names = vec!["A".to_string()];
        let result = load_config(f.path().to_str().unwrap(), &names);
        assert!(result.is_err());
    }

    #[test]
    fn test_load_config_too_few_columns() {
        let f = write_config("sample\tdphom_min\tdphom_max\nA\t5\t100\n");
        let names = vec!["A".to_string()];
        let result = load_config(f.path().to_str().unwrap(), &names);
        assert!(result.is_err());
    }

    #[test]
    fn test_load_config_empty_file() {
        let f = write_config("");
        let names = vec!["A".to_string()];
        let result = load_config(f.path().to_str().unwrap(), &names);
        assert!(result.is_err());
    }

    #[test]
    fn test_load_config_nonexistent_file() {
        let result = load_config("/nonexistent/path/config.tsv", &["A".to_string()]);
        assert!(result.is_err());
    }

    #[test]
    fn test_load_config_skips_comments() {
        // 注释行（# 开头）被跳过，表头被跳过，数据行正常解析
        let f = write_config(
            "# 这是一个注释\n\
             # 另一行注释\n\
             sample\tdphom_min\tdphom_max\tdphet_min\tdphet_max\n\
             A\t5\t100\t3\t50\n",
        );
        let names = vec!["A".to_string()];
        let cfg = load_config(f.path().to_str().unwrap(), &names).unwrap();
        assert!(cfg[0].is_some());
        assert_eq!(cfg[0].unwrap().dphom_min, 5);
    }

    #[test]
    fn test_sample_thresholds_from_global() {
        let th = SampleThresholds::from_global(7, 4);
        assert_eq!(th.dphom_min, 7);
        assert_eq!(th.dphom_max, u32::MAX);
        assert_eq!(th.dphet_min, 4);
        assert_eq!(th.dphet_max, u32::MAX);
    }
}
