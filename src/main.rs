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

/// 高性能VCF过滤工具（High-performance VCF filtering tool）
/// 核心功能：保留原始Header、INFO总DP，从AD提取DV，支持管道
/// https://bitbucket.org/ipk_dg_public/vcf_filtering
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
struct Args {
    /// 输入VCF文件路径（支持.gz/bgzip压缩，不指定则从标准输入/管道读取）
    /// Input VCF file path (supports .gz/bgzip compression, read from stdin/pipeline if not specified)
    #[arg(short, long)]
    input: Option<String>,

    /// 输出VCF文件路径（支持.gz/bgzip压缩，不指定则输出到标准输出）
    /// Output VCF file path (supports .gz/bgzip compression, output to stdout if not specified)
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
    /// Minimum quality threshold [default: 0.0]
    #[arg(long, default_value_t = 0.0)]
    minqual: f64,

    /// INFO字段最小DP阈值 [默认: 0]
    /// Minimum DP threshold in INFO field [default: 0]
    #[arg(long, default_value_t = 0)]
    mindp: u32,

    /// 最小纯合样本数阈值 [默认: 0]
    /// Minimum number of homozygous samples [default: 0]
    #[arg(long, default_value_t = 0)]
    minhomn: u32,

    /// 有效样本占比阈值 (present/(present+n)) [默认: 0.0]
    /// Valid sample ratio threshold (present/(present+n)) [default: 0.0]
    #[arg(long, default_value_t = 0.0)]
    minpresent: f64,

    /// 纯合样本占有效样本比阈值 (A+B/present) [默认: 0.0]
    /// Homozygous ratio in valid samples (A+B/present) [default: 0.0]
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
    /// Whether input is compressed (only for stdin, auto-detect if not specified)
    #[arg(long)]
    input_compressed: bool,
}

/// 基因型统计结果结构体
/// 存储每个位点的基因型分类计数和处理后的样本基因型字符串
#[derive(Debug, Default)]
struct GenotypeStats {
    a_count: u32,    // 0/0 纯合参考基因型数量
    b_count: u32,    // 1/1 纯合替代基因型数量
    h_count: u32,    // 0/1 杂合基因型数量
    n_count: u32,    // ./. 缺失基因型数量
    present: u32,    // 有效基因型数量（a + b + h）
    sample_gts: Vec<String>, // 处理后的样本基因型字符串（GT:DP:DV）
}

/// 从INFO字段中提取指定数值（如DP、MQ）
/// 输入：INFO字段字符串、匹配正则表达式
/// 输出：提取的数值（失败时返回0）
fn extract_info_value(info: &str, pattern: &Regex) -> Result<u32> {
    pattern
        .captures(info)
        .and_then(|caps| caps.get(1))
        .map(|m| m.as_str().parse::<u32>())
        .unwrap_or(Ok(0))
        .context("Failed to extract value from INFO field")
}

/// 动态解析FORMAT字段，返回字段名到索引的映射
/// 适配不同的FORMAT字段顺序（如GT:PL:DP:AD或GT:DP:AD等）
fn parse_format_fields(format_str: &str) -> HashMap<&str, usize> {
    let mut field_map = HashMap::new();
    for (idx, field) in format_str.split(':').enumerate() {
        field_map.insert(field, idx);
    }
    field_map
}

/// 从样本基因型字符串中提取DP和DV
/// DP：直接从FORMAT的DP字段提取
/// DV：从AD字段提取（AD格式为"REF深度,ALT深度"，取第二个值）
fn extract_dp_dv(gt_str: &str, format_map: &HashMap<&str, usize>) -> (u32, u32) {
    let gt_parts: Vec<&str> = gt_str.split(':').collect();

    // 提取DP（动态索引，适配不同FORMAT顺序）
    let dp = format_map.get("DP")
        .and_then(|&idx| gt_parts.get(idx))
        .and_then(|s| s.parse::<u32>().ok())
        .unwrap_or(0);

    // 提取DV（从AD字段的第二个值获取，AD格式：ref_depth,alt_depth）
    let dv = format_map.get("AD")
        .and_then(|&idx| gt_parts.get(idx))
        .and_then(|ad_str| {
            ad_str.split(',')
                .nth(1)  // 取ALT等位基因的深度
                .and_then(|s| s.parse::<u32>().ok())
        })
        .unwrap_or(0);

    (dp, dv)
}

/// 处理单行VCF数据（核心过滤逻辑）
/// 步骤：1. 基础过滤 → 2. 动态解析FORMAT → 3. 基因型分类 → 4. 群体水平过滤 → 5. MAF过滤 → 6. 组装输出行
fn process_vcf_line(
    line: &str,
    args: &Args,
    dp_re: &Regex,
    mq_re: &Regex,
) -> Result<Option<String>> {
    let parts: Vec<&str> = line.split('\t').collect();
    if parts.len() < 10 {
        return Ok(None); // 列数不足，跳过无效行
    }

    // 1. 基础过滤：多碱基变异、低质量、参考碱基为N、低总DP
    let ref_base = parts[3];
    let alt_base = parts[4];
    let qual = parts[5].parse::<f64>().unwrap_or(0.0);
    let info = parts[7];
    let format_str = parts[8];

    // 提取INFO中的总DP（用于过滤和保留）
    let info_dp = extract_info_value(info, dp_re)?;

    // 过滤条件：
    // - 参考/替代碱基长度>1（多碱基变异）
    // - 质量值<阈值
    // - 参考碱基为N
    // - INFO总DP<阈值
    if ref_base.len() > 1
        || alt_base.len() > 1
        || qual < args.minqual
        || ref_base == "N"
        || info_dp < args.mindp
    {
        return Ok(None);
    }

    // 2. 动态解析FORMAT字段（必须包含DP字段，否则跳过）
    let format_map = parse_format_fields(format_str);
    if !format_map.contains_key("DP") {
        return Ok(None); // 无DP字段，跳过该位点
    }

    // 3. 处理每个样本的基因型，进行分类计数
    let mut stats = GenotypeStats::default();
    for i in 9..parts.len() {
        let gt_str = parts[i];
        let (dp, dv) = extract_dp_dv(gt_str, &format_map);

        // 计算替代等位基因频率（r = DV/DP）
        let r = if dp > 0 { dv as f64 / dp as f64 } else { 0.0 };
        let mut gt_result = String::new();

        // 基因型分类逻辑（完全复刻原AWK脚本）
        if dp == 0 && alt_base == "." {
            // DP=0且无替代等位基因 → 纯合参考（0/0）
            gt_result = format!("0/0:{dp}:{dv}");
            stats.a_count += 1;
        } else if dp >= args.dphom && r <= args.tol {
            // DP≥纯合阈值且r≤容差 → 纯合参考（0/0）
            gt_result = format!("0/0:{dp}:{dv}");
            stats.a_count += 1;
        } else if dp >= args.dphom && r >= 1.0 - args.tol {
            // DP≥纯合阈值且r≥1-容差 → 纯合替代（1/1）
            gt_result = format!("1/1:{dp}:{dv}");
            stats.b_count += 1;
        } else if dp >= args.dphet && r >= 0.5 - args.tol && r <= 0.5 + args.tol {
            // DP≥杂合阈值且r在0.5±容差范围内 → 杂合（0/1）
            gt_result = format!("0/1:{dp}:{dv}");
            stats.h_count += 1;
        } else {
            // 不满足任何条件 → 缺失基因型（./.）
            gt_result = format!("./.:{dp}:{dv}");
            stats.n_count += 1;
        }

        stats.sample_gts.push(gt_result);
    }

    // 计算有效基因型数量（present = a_count + b_count + h_count）
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

    // 过滤条件：
    // - 无有效基因型（present=0）
    // - 纯合参考/替代样本数<阈值
    // - 有效样本占比<阈值
    // - 纯合样本占有效样本比<阈值
    if stats.present == 0
        || stats.a_count < args.minhomn
        || stats.b_count < args.minhomn
        || present_ratio < args.minpresent
        || hom_ratio < args.minhomp
    {
        return Ok(None);
    }

    // 5. MAF（最小等位基因频率）过滤
    let maf = if stats.b_count > stats.a_count {
        // 替代等位基因为次要等位基因
        (2 * stats.a_count + stats.h_count) as f64 / (2 * stats.present) as f64
    } else {
        // 参考等位基因为次要等位基因
        (2 * stats.b_count + stats.h_count) as f64 / (2 * stats.present) as f64
    };

    if maf < args.minmaf {
        return Ok(None); // MAF<阈值，跳过
    }

    // 6. 提取MQ并组装输出行
    let mq = extract_info_value(info, mq_re)?;
    let new_info = format!("DP={info_dp};MQ={mq}"); // 保留INFO总DP + 添加MQ

    // 组装输出行（CHROM → 样本基因型）
    let mut output_parts = vec![
        parts[0].to_string(),
        parts[1].to_string(),
        parts[2].to_string(),
        parts[3].to_string(),
        parts[4].to_string(),
        parts[5].to_string(),
        parts[6].to_string(),
        new_info,
        "GT:DP:DV".to_string(),
    ];
    output_parts.extend(stats.sample_gts);

    Ok(Some(output_parts.join("\t")))
}

/// 生成过滤规则注释行（写入VCF头文件，便于实验追溯）
/// 格式：##FilterRule=<ID=CustomFilter,Description="参数说明">
fn generate_filter_comment(args: &Args) -> String {
    format!(
        "##FilterRule=<ID=CustomFilter,Description=\"dphom={}, dphet={}, tol={}, minqual={}, mindp={}, minhomn={}, minpresent={}, minhomp={}, minmaf={}\">",
        args.dphom, args.dphet, args.tol, args.minqual, args.mindp, args.minhomn, args.minpresent, args.minhomp, args.minmaf
    )
}

/// 打开输入（支持文件/标准输入，文件支持.gz/bgzip压缩，stdin支持手动指定压缩）
fn open_input(args: &Args) -> Result<Box<dyn BufRead>> {
    match &args.input {
        Some(path) => {
            // 输入来自文件
            let file = File::open(path).context("Failed to open input file")?;
            if path.ends_with(".gz") {
                let decoder = MultiGzDecoder::new(file);
                Ok(Box::new(BufReader::new(decoder)))
            } else {
                Ok(Box::new(BufReader::new(file)))
            }
        }
        None => {
            // 输入来自标准输入（管道）
            let stdin = io::stdin();
            let stdin_lock = stdin.lock();
            if args.input_compressed {
                // 手动指定输入为压缩格式（如zcat input.vcf.gz | ./vcf_filter --input-compressed）
                let decoder = MultiGzDecoder::new(stdin_lock);
                Ok(Box::new(BufReader::new(decoder)))
            } else {
                // 非压缩格式（如cat input.vcf | ./vcf_filter）
                Ok(Box::new(BufReader::new(stdin_lock)))
            }
        }
    }
}

/// 打开输出（支持文件/标准输出，文件支持.gz/bgzip压缩）
fn open_output(output_path: &Option<String>, compress_level: u32) -> Result<Box<dyn Write>> {
    match output_path {
        Some(path) => {
            // 输出到文件
            let file = File::create(path).context("Failed to create output file")?;
            let compress_level = Compression::new(compress_level as u32);
            if path.ends_with(".gz") {
                let encoder = GzEncoder::new(file, compress_level);
                Ok(Box::new(BufWriter::new(encoder)))
            } else {
                Ok(Box::new(BufWriter::new(file)))
            }
        }
        None => {
            // 输出到标准输出（禁用压缩）
            let stdout = io::stdout();
            Ok(Box::new(BufWriter::new(stdout.lock())))
        }
    }
}

fn main() -> Result<()> {
    let args = Args::parse();
    println!("Using parameters: {:?}", args);

    // 预编译正则表达式（避免重复编译，提升处理性能）
    let dp_re = Regex::new(r"DP=(\d+)").context("Failed to compile DP regex")?;
    let mq_re = Regex::new(r"MQ=(\d+)").context("Failed to compile MQ regex")?;

    // 打开输入（支持文件/标准输入）和输出（支持文件/标准输出）
    let mut reader = open_input(&args)?;
    let mut writer = open_output(&args.output, args.compress_level)?;

    // 打印输入输出来源提示
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

    // 第一步：读取并写入所有原始Header（完整保留，仅添加过滤规则注释）
    let mut buf = String::new();
    let mut filter_comment_written = false;
    let mut header_finished = false;

    while reader.read_line(&mut buf)? > 0 {
        let line = buf.trim_end(); // 去除换行符
        if line.is_empty() {
            buf.clear();
            continue;
        }

        // 处理Header行（以#开头）
        if line.starts_with('#') {
            // 在##fileformat行后插入过滤规则注释（保证VCF格式合法）
            if !filter_comment_written && line.starts_with("##fileformat=") {
                writeln!(writer, "{}", line)?;
                let filter_comment = generate_filter_comment(&args);
                writeln!(writer, "{}", filter_comment)?;
                filter_comment_written = true;
            } else {
                // 原样写入其他所有原始Header行
                writeln!(writer, "{}", line)?;
            }

            // 检测Header结束（#CHROM行是最后一个Header行）
            if line.starts_with("#CHROM") {
                header_finished = true;
            }
            buf.clear();
            continue;
        }

        // Header处理完成，跳出循环开始处理数据行
        if header_finished {
            break;
        }

        buf.clear();
    }

    // 第二步：逐行处理数据并实时输出（边处理边输出，适配管道流）
    let mut line_num = 0;
    loop {
        buf.clear();
        let bytes_read = reader.read_line(&mut buf)?;
        if bytes_read == 0 {
            break; // 读取完毕，退出循环
        }
        line_num += 1;

        let line = buf.trim_end();
        if line.is_empty() || line.starts_with('#') {
            continue; // 跳过空行或遗漏的Header行
        }

        // 处理数据行，符合条件则立即输出（管道流实时处理）
        match process_vcf_line(line, &args, &dp_re, &mq_re) {
            Ok(Some(output_line)) => {
                writeln!(writer, "{}", output_line)?;
                // 强制刷新缓冲区（管道模式下确保实时输出，避免数据滞留）
                writer.flush()?;
            }
            Ok(None) => {} // 不满足过滤条件，跳过
            Err(e) => {
                eprintln!("Failed to process line {}: {:?}, skipping", line_num, e);
            }
        }
    }

    // 刷新输出缓冲区，确保所有数据写入
    writer.flush()?;

    // 打印完成提示
    match (&args.input, &args.output) {
        (_, Some(path)) => println!("Processing completed! Output file: {}", path),
        (_, None) => println!("Processing completed! Output to stdout"),
    }

    Ok(())
}