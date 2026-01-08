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

    /// 压缩级别 (1-9, 6=平衡) [默认: 6，仅对文件输出生效]
    /// Compression level (1-9, 6=balanced) [default: 6, only effective for file output]
    #[arg(long, default_value_t = 6)]
    pub compress_level: u32,

    /// 输入是否为压缩格式（仅对标准输入有效，自动检测则不指定）
    /// Whether the input is compressed (only valid for stdin; auto-detect if not specified)
    #[arg(long)]
    pub input_compressed: bool,
}