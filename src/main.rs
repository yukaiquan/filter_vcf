mod args;
mod config;
mod input;
mod vcf;
use crate::args::Args;
use crate::config::load_config;
use crate::input::{open_input, open_output};
use crate::vcf::{generate_filter_comment, process_vcf_line};
use anyhow::{Context, Result};
use clap::Parser;
use rayon::ThreadPoolBuilder;
use regex::Regex;
use std::io::{BufRead, Read, Write};

fn main() -> Result<()> {
    let mut args = Args::parse();
    // println!("Using parameters: {:?}", args);

    // 配置并行处理线程数
    if args.threads > 0 {
        ThreadPoolBuilder::new()
            .num_threads(args.threads)
            .build_global()
            .context("Failed to configure thread pool")?;
    }

    // 预编译DP正则
    let dp_re = Regex::new(r"DP=(\d+)").context("Failed to compile DP regex")?;
    // 生成过滤规则注释行
    let filter_comment = generate_filter_comment(&args);

    // 打开输入输出
    let mut reader = open_input(&args)?;
    let mut writer = open_output(&args.output, args.compress_level)?;

    // 打印输入输出来源
    match (&args.input, &args.output) {
        (Some(in_path), Some(out_path)) => {
            // println!("Reading from file: {}, Writing to file: {}", in_path, out_path);
        }
        (Some(in_path), None) => {
            // println!("Reading from file: {}, Writing to stdout", in_path);
        }
        (None, Some(out_path)) => {
            // println!("Reading from stdin, Writing to file: {}", out_path);
        }
        (None, None) => {
            // println!("Reading from stdin, Writing to stdout (full pipeline mode)");
        }
    }

    // ========== 关键修改：收集所有Header行，在#CHROM前插入过滤规则 ==========
    let mut buf = String::new();
    let mut header_lines = Vec::new(); // 收集所有原始Header行
    let mut chrom_header = String::new(); // 存储#CHROM行
    let mut header_finished = false;
    // 保存 #CHROM 后的第一条数据行（while 循环 break 时 buf 中残留的行）
    let mut first_data_line: Option<String> = None;

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

        // Header处理完成，跳出循环（此时 buf 中是第一条数据行）
        if header_finished {
            first_data_line = Some(line);
            break;
        }

        buf.clear();
    }

    // 解析样本名（#CHROM 行第 10 列起），并按需加载按样本深度配置
    {
        let cols: Vec<&str> = chrom_header.split('\t').collect();
        let sample_names: Vec<String> = cols.iter().skip(9).map(|s| s.to_string()).collect();
        if let Some(cfg_path) = &args.config {
            args.sample_config = Some(load_config(cfg_path, &sample_names)?);
        }
    }

    // 1. 写入所有原始非#CHROM Header行
    for header_line in header_lines {
        writeln!(writer, "{}", header_line)?;
    }

    // 2. 写入过滤规则注释行
    writeln!(writer, "{}", filter_comment)?;

    // 3. 写入#CHROM行
    writeln!(writer, "{}", chrom_header)?;

    // ========== 处理数据行 ==========
    let mut line_num = 0;

    // 先处理 while 循环中残留的第一条数据行
    if let Some(ref first_line) = first_data_line {
        line_num += 1;
        if !first_line.is_empty() && !first_line.starts_with('#') {
            match process_vcf_line(first_line, &args, &dp_re) {
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
    }

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
    // match (&args.input, &args.output) {
    //     (_, Some(path)) => println!("Processing completed! Output file: {}", path),
    //     (_, None) => println!("Processing completed! Output to stdout"),
    // }

    Ok(())
}
