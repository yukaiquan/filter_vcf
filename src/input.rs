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
use crate::args::Args;



/// 打开输入（支持文件/标准输入，压缩自动检测）
pub fn open_input(args: &Args) -> anyhow::Result<Box<dyn BufRead>> {
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
pub fn open_output(output_path: &Option<String>, compress_level: u32) -> anyhow::Result<Box<dyn Write>> {
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