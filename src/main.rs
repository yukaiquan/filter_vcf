use anyhow::{Context, Result};
use clap::Parser;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use flate2::Compression;
use regex::Regex;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

/// VCF Filtering Tool (https://bitbucket.org/ipk_dg_public/vcf_filtering/)
/// kaiquanyu@icloud.com
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// 输入VCF文件路径（支持.gz/bgzip压缩，不指定则从标准输入/管道读取）
    /// Input VCF file path (supports .gz/bgzip compression; reads from stdin/pipeline if not specified)
    #[arg(short, long)]
    input: Option<String>,

    /// 输出VCF文件路径（支持.gz/bgzip压缩，不指定则输出到标准输出）
    /// Output VCF file path (supports .gz/bgzip compression; writes to stdout if not specified)
    #[arg(short, long)]
    output: Option<String>,

    /// 纯合基因型最小DP阈值 [默认: 1]
    /// Minimum DP threshold for homozygous genotypes [default: 1]
    #[arg(long, default_value_t = 1)]
    dphom: u32,

    /// 杂合基因型最小DP阈值 [默认: 1]
    /// Minimum DP threshold for heterozygous genotypes [default: 1]
    #[arg(long, default_value_t = 1)]
    dphet: u32,

    /// 频率容差阈值 [默认: 0.2499]
    /// Frequency tolerance threshold [default: 0.2499]
    #[arg(long, default_value_t = 0.2499)]
    tol: f64,

    /// 最小质量值阈值 [默认: 0]
    /// Minimum quality score threshold [default: 0]
    #[arg(long, default_value_t = 0.0)]
    minqual: f64,

    /// INFO字段最小DP阈值 [默认: 0]
    /// Minimum DP threshold in INFO field [default: 0]
    #[arg(long, default_value_t = 0)]
    mindp: u32,

    /// 最小纯合样本数阈值 [默认: 0]
    /// Minimum number of homozygous samples threshold [default: 0]
    #[arg(long, default_value_t = 0)]
    minhomn: u32,

    /// 有效样本占比阈值 (present/(present+n)) [默认: 0.0]
    /// Valid sample ratio threshold (present/(present+n)) [default: 0.0]
    #[arg(long, default_value_t = 0.0)]
    minpresent: f64,

    /// 纯合样本占有效样本比阈值 (A+B/present) [默认: 0.0]
    /// Homozygous sample ratio in valid samples threshold (A+B/present) [default: 0.0]
    #[arg(long, default_value_t = 0.0)]
    minhomp: f64,

    /// 最小MAF阈值 [默认: 0.0]
    /// Minimum MAF (Minor Allele Frequency) threshold [default: 0.0]
    #[arg(long, default_value_t = 0.0)]
    minmaf: f64,

    /// 压缩级别 (1-9, 6=平衡) [默认: 6，仅对文件输出生效]
    /// Compression level (1-9, 6=balanced) [default: 6, only effective for file output]
    #[arg(long, default_value_t = 6)]
    compress_level: u32,

    /// 输入是否为压缩格式（仅对标准输入有效，自动检测则不指定）
    /// Whether the input is compressed (only valid for stdin; auto-detect if not specified)
    #[arg(long)]
    input_compressed: bool,
}

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
fn process_vcf_line(
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
fn generate_filter_comment(args: &Args) -> String {
    format!(
        "##FilterRule=<ID=CustomFilter,Description=\"dphom={}, dphet={}, tol={}, minqual={}, mindp={}, minhomn={}, minpresent={}, minhomp={}, minmaf={}\">",
        args.dphom, args.dphet, args.tol, args.minqual, args.mindp, args.minhomn, args.minpresent, args.minhomp, args.minmaf
    )
}

/// 打开输入（支持文件/标准输入，压缩自动检测）
fn open_input(args: &Args) -> Result<Box<dyn BufRead>> {
    match &args.input {
        Some(path) => {
            let file = File::open(path).context("Failed to open input file")?;
            if path.ends_with(".gz") {
                let decoder = MultiGzDecoder::new(file);
                Ok(Box::new(BufReader::new(decoder)))
            } else {
                Ok(Box::new(BufReader::new(file)))
            }
        }
        None => {
            let stdin = io::stdin();
            let stdin_lock = stdin.lock();
            if args.input_compressed {
                let decoder = MultiGzDecoder::new(stdin_lock);
                Ok(Box::new(BufReader::new(decoder)))
            } else {
                Ok(Box::new(BufReader::new(stdin_lock)))
            }
        }
    }
}

/// 打开输出（支持文件/标准输出，压缩可选）
fn open_output(output_path: &Option<String>, compress_level: u32) -> Result<Box<dyn Write>> {
    match output_path {
        Some(path) => {
            let file = File::create(path).context("Failed to create output file")?;
            // 修正：Compression::new需要u8类型
            let compress_level = Compression::new(compress_level as u32);
            if path.ends_with(".gz") {
                let encoder = GzEncoder::new(file, compress_level);
                Ok(Box::new(BufWriter::new(encoder)))
            } else {
                Ok(Box::new(BufWriter::new(file)))
            }
        }
        None => {
            let stdout = io::stdout();
            Ok(Box::new(BufWriter::new(stdout.lock())))
        }
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    println!("Using parameters: {:?}", args);

    // 预编译DP正则（仅提取INFO的DP）
    let dp_re = Regex::new(r"DP=(\d+)").context("Failed to compile DP regex")?;
    // 生成过滤规则注释行
    let filter_comment = generate_filter_comment(&args);

    // 打开输入输出
    let mut reader = open_input(&args)?;
    let mut writer = open_output(&args.output, args.compress_level)?;

    // 打印输入输出来源
    match (&args.input, &args.output) {
        (Some(in_path), Some(out_path)) => {
            println!("Reading from file: {}, Writing to file: {}", in_path, out_path);
        }
        (Some(in_path), None) => {
            println!("Reading from file: {}, Writing to stdout", in_path);
        }
        (None, Some(out_path)) => {
            println!("Reading from stdin, Writing to file: {}", out_path);
        }
        (None, None) => {
            println!("Reading from stdin, Writing to stdout (full pipeline mode)");
        }
    }

    // ========== 关键修改：收集所有Header行，在#CHROM前插入过滤规则 ==========
    let mut buf = String::new();
    let mut header_lines = Vec::new(); // 收集所有原始Header行
    let mut chrom_header = String::new(); // 存储#CHROM行
    let mut header_finished = false;

    while reader.read_line(&mut buf)? > 0 {
        let line = buf.trim_end().to_string();
        if line.is_empty() {
            buf.clear();
            continue;
        }

        if line.starts_with('#') {
            if line.starts_with("#CHROM") {
                // 分离#CHROM行，后续最后写入
                chrom_header = line;
                header_finished = true;
            } else {
                // 收集所有非#CHROM的Header行
                header_lines.push(line);
            }
            buf.clear();
            continue;
        }

        // Header处理完成，跳出循环
        if header_finished {
            break;
        }

        buf.clear();
    }

    // 1. 写入所有原始非#CHROM Header行
    for header_line in header_lines {
        writeln!(writer, "{}", header_line)?;
    }

    // 2. 写入过滤规则注释行（Header最后面，#CHROM之前）
    writeln!(writer, "{}", filter_comment)?;

    // 3. 写入#CHROM行（Header的最后一行）
    writeln!(writer, "{}", chrom_header)?;

    // ========== 处理数据行（逻辑不变） ==========
    let mut line_num = 0;
    loop {
        buf.clear();
        let bytes_read = reader.read_line(&mut buf)?;
        if bytes_read == 0 {
            break;
        }
        line_num += 1;

        let line = buf.trim_end();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }

        // 过滤逻辑：符合条件则输出原始行
        match process_vcf_line(line, &args, &dp_re) {
            Ok(Some(output_line)) => {
                writeln!(writer, "{}", output_line)?;
                writer.flush()?;
            }
            Ok(None) => {} // 不符合条件，跳过
            Err(e) => {
                eprintln!("Failed to process line {}: {:?}, skipping", line_num, e);
            }
        }
    }

    // 刷新缓冲区
    writer.flush()?;

    // 完成提示
    match (&args.input, &args.output) {
        (_, Some(path)) => println!("Processing completed! Output file: {}", path),
        (_, None) => println!("Processing completed! Output to stdout"),
    }

    Ok(())
}