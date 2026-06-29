use clap::Parser;

/// VCF Filtering Tool (https://bitbucket.org/ipk_dg_public/vcf_filtering/)
/// kaiquanyu@icloud.com
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
pub struct Args {
    /// 输入VCF文件路径（支持.gz/bgzip压缩，不指定则从标准输入/管道读取）
    /// Input VCF file path (supports .gz/bgzip compression; reads from stdin/pipeline if not specified)
    #[arg(short, long)]
    pub input: Option<String>,

    /// 输出VCF文件路径（支持.gz/bgzip压缩，不指定则输出到标准输出）
    /// Output VCF file path (supports .gz/bgzip compression; writes to stdout if not specified)
    #[arg(short, long)]
    pub output: Option<String>,

    /// 纯合基因型最小DP阈值 [默认: 1]
    /// Minimum DP threshold for homozygous genotypes [default: 1]
    #[arg(long, default_value_t = 1)]
    pub dphom: u32,

    /// 杂合基因型最小DP阈值 [默认: 1]
    /// Minimum DP threshold for heterozygous genotypes [default: 1]
    #[arg(long, default_value_t = 1)]
    pub dphet: u32,

    /// 频率容差阈值 [默认: 0.2499]
    /// Frequency tolerance threshold [default: 0.2499]
    #[arg(long, default_value_t = 0.2499)]
    pub tol: f64,

    /// 最小质量值阈值 [默认: 0]
    /// Minimum quality score threshold [default: 0]
    #[arg(long, default_value_t = 0.0)]
    pub minqual: f64,

    /// INFO字段最小DP阈值 [默认: 0]
    /// Minimum DP threshold in INFO field [default: 0]
    #[arg(long, default_value_t = 0)]
    pub mindp: u32,

    /// 最小纯合样本数阈值 [默认: 0]
    /// Minimum number of homozygous samples threshold [default: 0]
    #[arg(long, default_value_t = 0)]
    pub minhomn: u32,

    /// 有效样本占比阈值 (present/(present+n)) [默认: 0.0]
    /// Valid sample ratio threshold (present/(present+n)) [default: 0.0]
    #[arg(long, default_value_t = 0.0)]
    pub minpresent: f64,

    /// 纯合样本占有效样本比阈值 (A+B/present) [默认: 0.0]
    /// Homozygous sample ratio in valid samples threshold (A+B/present) [default: 0.0]
    #[arg(long, default_value_t = 0.0)]
    pub minhomp: f64,

    /// 最小MAF阈值 [默认: 0.0]
    /// Minimum MAF (Minor Allele Frequency) threshold [default: 0.0]
    #[arg(long, default_value_t = 0.0)]
    pub minmaf: f64,

    /// 是否包含indel [默认: false，仅处理SNP]
    /// Whether to include indels [default: false, SNP only]
    #[arg(long, default_value_t = false)]
    pub include_indel: bool,

    /// 压缩级别 (1-9, 6=平衡) [默认: 6，仅对文件输出生效]
    /// Compression level (1-9, 6=balanced) [default: 6, only effective for file output]
    #[arg(long, default_value_t = 6)]
    pub compress_level: u32,

    /// 输入是否为压缩格式（仅对标准输入有效，自动检测则不指定）
    /// Whether the input is compressed (only valid for stdin; auto-detect if not specified)
    #[arg(long)]
    pub input_compressed: bool,

    /// 并行处理线程数 [默认: 1=自动使用1个CPU核心]
    /// Number of threads for parallel processing [default: 1=auto-use 1 CPU cores]
    #[arg(long, default_value_t = 1)]
    pub threads: usize,

    /// 样本深度配置文件路径（TSV：sample,dphom_min,dphom_max,dphet_min,dphet_max）
    /// Sample depth config file (TSV). 不指定则所有样本使用全局 --dphom/--dphet
    #[arg(long)]
    pub config: Option<String>,

    /// [内部字段] 按样本列对齐的已加载阈值表，非命令行参数
    /// [internal] per-sample loaded thresholds aligned by column index, not a CLI arg
    #[arg(skip)]
    pub sample_config: Option<crate::config::SampleConfig>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_args() {
        let args = Args::try_parse_from(vec!["filter_vcf"]).unwrap();

        assert_eq!(args.dphom, 1);
        assert_eq!(args.dphet, 1);
        assert!((args.tol - 0.2499).abs() < 0.0001);
        assert_eq!(args.minqual, 0.0);
        assert_eq!(args.mindp, 0);
        assert_eq!(args.minhomn, 0);
        assert_eq!(args.minpresent, 0.0);
        assert_eq!(args.minhomp, 0.0);
        assert_eq!(args.minmaf, 0.0);
        assert!(!args.include_indel);
        assert_eq!(args.compress_level, 6);
        assert!(!args.input_compressed);
        assert_eq!(args.threads, 1);
    }

    #[test]
    fn test_parse_threads() {
        // 测试默认值（1=自动）
        let args = Args::try_parse_from(vec!["filter_vcf"]).unwrap();
        assert_eq!(args.threads, 1);

        // 测试指定线程数
        let args = Args::try_parse_from(vec!["filter_vcf", "--threads", "4"]).unwrap();
        assert_eq!(args.threads, 4);

        let args = Args::try_parse_from(vec!["filter_vcf", "--threads", "8"]).unwrap();
        assert_eq!(args.threads, 8);
    }

    #[test]
    fn test_parse_input_output() {
        let args = Args::try_parse_from(vec!["filter_vcf", "-i", "input.vcf", "-o", "output.vcf"])
            .unwrap();

        assert_eq!(args.input, Some("input.vcf".to_string()));
        assert_eq!(args.output, Some("output.vcf".to_string()));
    }

    #[test]
    fn test_parse_long_options() {
        let args = Args::try_parse_from(vec![
            "filter_vcf",
            "--input",
            "input.vcf.gz",
            "--output",
            "output.vcf.gz",
        ])
        .unwrap();

        assert_eq!(args.input, Some("input.vcf.gz".to_string()));
        assert_eq!(args.output, Some("output.vcf.gz".to_string()));
    }

    #[test]
    fn test_parse_dphom() {
        let args = Args::try_parse_from(vec!["filter_vcf", "--dphom", "5"]).unwrap();

        assert_eq!(args.dphom, 5);
    }

    #[test]
    fn test_parse_dphet() {
        let args = Args::try_parse_from(vec!["filter_vcf", "--dphet", "10"]).unwrap();

        assert_eq!(args.dphet, 10);
    }

    #[test]
    fn test_parse_tol() {
        let args = Args::try_parse_from(vec!["filter_vcf", "--tol", "0.3"]).unwrap();

        assert!((args.tol - 0.3).abs() < 0.0001);
    }

    #[test]
    fn test_parse_minqual() {
        let args = Args::try_parse_from(vec!["filter_vcf", "--minqual", "30.5"]).unwrap();

        assert!((args.minqual - 30.5).abs() < 0.01);
    }

    #[test]
    fn test_parse_mindp() {
        let args = Args::try_parse_from(vec!["filter_vcf", "--mindp", "100"]).unwrap();

        assert_eq!(args.mindp, 100);
    }

    #[test]
    fn test_parse_minhomn() {
        let args = Args::try_parse_from(vec!["filter_vcf", "--minhomn", "2"]).unwrap();

        assert_eq!(args.minhomn, 2);
    }

    #[test]
    fn test_parse_minpresent() {
        let args = Args::try_parse_from(vec!["filter_vcf", "--minpresent", "0.9"]).unwrap();

        assert!((args.minpresent - 0.9).abs() < 0.01);
    }

    #[test]
    fn test_parse_minhomp() {
        let args = Args::try_parse_from(vec!["filter_vcf", "--minhomp", "0.8"]).unwrap();

        assert!((args.minhomp - 0.8).abs() < 0.01);
    }

    #[test]
    fn test_parse_minmaf() {
        let args = Args::try_parse_from(vec!["filter_vcf", "--minmaf", "0.05"]).unwrap();

        assert!((args.minmaf - 0.05).abs() < 0.001);
    }

    #[test]
    fn test_parse_include_indel() {
        // 测试不指定参数时的默认值
        let args = Args::try_parse_from(vec!["filter_vcf"]).unwrap();
        assert!(!args.include_indel);

        // 测试指定参数
        let args = Args::try_parse_from(vec!["filter_vcf", "--include-indel"]).unwrap();
        assert!(args.include_indel);
    }

    #[test]
    fn test_parse_compress_level() {
        let args = Args::try_parse_from(vec!["filter_vcf", "--compress-level", "9"]).unwrap();

        assert_eq!(args.compress_level, 9);
    }

    #[test]
    fn test_parse_input_compressed() {
        // 测试不指定参数时的默认值
        let args = Args::try_parse_from(vec!["filter_vcf"]).unwrap();
        assert!(!args.input_compressed);

        // 测试指定参数
        let args = Args::try_parse_from(vec!["filter_vcf", "--input-compressed"]).unwrap();
        assert!(args.input_compressed);
    }

    #[test]
    fn test_parse_multiple_parameters() {
        let args = Args::try_parse_from(vec![
            "filter_vcf",
            "-i",
            "input.vcf.gz",
            "-o",
            "output.vcf.gz",
            "--dphom",
            "5",
            "--dphet",
            "10",
            "--minqual",
            "30.0",
            "--mindp",
            "100",
            "--minmaf",
            "0.01",
            "--include-indel",
        ])
        .unwrap();

        assert_eq!(args.input, Some("input.vcf.gz".to_string()));
        assert_eq!(args.output, Some("output.vcf.gz".to_string()));
        assert_eq!(args.dphom, 5);
        assert_eq!(args.dphet, 10);
        assert!((args.minqual - 30.0).abs() < 0.01);
        assert_eq!(args.mindp, 100);
        assert!((args.minmaf - 0.01).abs() < 0.001);
        assert!(args.include_indel);
    }

    #[test]
    fn test_help_flag() {
        let result = Args::try_parse_from(vec!["filter_vcf", "--help"]);
        assert!(result.is_err());
    }

    #[test]
    fn test_version_flag() {
        let result = Args::try_parse_from(vec!["filter_vcf", "--version"]);
        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_number() {
        let result = Args::try_parse_from(vec!["filter_vcf", "--dphom", "abc"]);
        assert!(result.is_err());
    }

    #[test]
    fn test_invalid_float() {
        let result = Args::try_parse_from(vec!["filter_vcf", "--minqual", "not_a_number"]);
        assert!(result.is_err());
    }
}
