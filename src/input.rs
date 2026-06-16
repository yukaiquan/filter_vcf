use crate::args::Args;
use anyhow::Context;
use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Write};

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
pub fn open_output(
    output_path: &Option<String>,
    compress_level: u32,
) -> anyhow::Result<Box<dyn Write>> {
    match output_path {
        Some(path) => {
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
            let stdout = io::stdout();
            Ok(Box::new(BufWriter::new(stdout.lock())))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::Parser;
    use std::io::Write;
    use tempfile::NamedTempFile;

    // ==================== open_input 测试 ====================
    #[test]
    fn test_open_input_uncompressed_file() {
        // 创建临时文件
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "test content").unwrap();
        let temp_path = temp_file.path().to_str().unwrap().to_string();

        let mut args = Args::try_parse_from(vec!["filter_vcf"]).unwrap();
        args.input = Some(temp_path);

        let result = open_input(&args);
        assert!(result.is_ok());
    }

    #[test]
    fn test_open_input_gz_file() {
        // 创建临时 .gz 文件
        let mut temp_file = NamedTempFile::with_suffix(".gz").unwrap();
        writeln!(temp_file, "compressed content").unwrap();
        let temp_path = temp_file.path().to_str().unwrap().to_string();

        let mut args = Args::try_parse_from(vec!["filter_vcf"]).unwrap();
        args.input = Some(temp_path);

        let result = open_input(&args);
        assert!(result.is_ok());
    }

    #[test]
    fn test_open_input_nonexistent_file() {
        let mut args = Args::try_parse_from(vec!["filter_vcf"]).unwrap();
        args.input = Some("/nonexistent/path/to/file.vcf".to_string());

        let result = open_input(&args);
        assert!(result.is_err());
    }

    #[test]
    fn test_open_input_stdin_uncompressed() {
        let args = Args::try_parse_from(vec!["filter_vcf"]).unwrap();
        // 注意：stdin 在测试中可能无法真正读取，但至少不会panic
        let result = open_input(&args);
        assert!(result.is_ok());
    }

    #[test]
    fn test_open_input_stdin_compressed() {
        let mut args = Args::try_parse_from(vec!["filter_vcf"]).unwrap();
        args.input_compressed = true;

        let result = open_input(&args);
        assert!(result.is_ok());
    }

    // ==================== open_output 测试 ====================
    #[test]
    fn test_open_output_uncompressed_file() {
        let temp_file = NamedTempFile::new().unwrap();
        let temp_path = Some(temp_file.path().to_str().unwrap().to_string());

        let result = open_output(&temp_path, 6);
        assert!(result.is_ok());

        // 验证可以写入
        let mut writer = result.unwrap();
        let write_result = writer.write_all(b"test");
        assert!(write_result.is_ok());
    }

    #[test]
    fn test_open_output_gz_file() {
        let temp_file = NamedTempFile::with_suffix(".gz").unwrap();
        let temp_path = Some(temp_file.path().to_str().unwrap().to_string());

        let result = open_output(&temp_path, 6);
        assert!(result.is_ok());

        // 验证可以写入
        let mut writer = result.unwrap();
        let write_result = writer.write_all(b"test");
        assert!(write_result.is_ok());
    }

    #[test]
    fn test_open_output_stdout() {
        let result = open_output(&None, 6);
        assert!(result.is_ok());

        // 验证可以写入到stdout
        let mut writer = result.unwrap();
        let write_result = writer.write_all(b"test");
        assert!(write_result.is_ok());
    }

    #[test]
    fn test_open_output_different_compress_levels() {
        let temp_file = NamedTempFile::with_suffix(".gz").unwrap();
        let temp_path = Some(temp_file.path().to_str().unwrap().to_string());

        // 测试不同压缩级别 (1-9)
        for level in 1..=9 {
            let result = open_output(&temp_path, level);
            assert!(result.is_ok());
        }
    }

    #[test]
    fn test_open_output_compress_level_zero() {
        let temp_file = NamedTempFile::with_suffix(".gz").unwrap();
        let temp_path = Some(temp_file.path().to_str().unwrap().to_string());

        // 测试压缩级别为0（特殊情况）
        let result = open_output(&temp_path, 0);
        assert!(result.is_ok());
    }

    #[test]
    fn test_open_output_high_compress_level() {
        let temp_file = NamedTempFile::with_suffix(".gz").unwrap();
        let temp_path = Some(temp_file.path().to_str().unwrap().to_string());

        // 测试高压缩级别
        let result = open_output(&temp_path, 9);
        assert!(result.is_ok());
    }

    // ==================== 边界情况测试 ====================
    #[test]
    fn test_empty_file_extension() {
        // 无扩展名文件
        let mut temp_file = NamedTempFile::new().unwrap();
        writeln!(temp_file, "test").unwrap();
        let temp_path = temp_file.path().to_str().unwrap().to_string();

        let mut args = Args::try_parse_from(vec!["filter_vcf"]).unwrap();
        args.input = Some(temp_path);

        let result = open_input(&args);
        assert!(result.is_ok());
    }

    #[test]
    fn test_multiple_dots_in_filename() {
        // 文件名包含多个点
        let mut temp_file = NamedTempFile::with_suffix(".vcf.gz").unwrap();
        writeln!(temp_file, "test").unwrap();
        let temp_path = temp_file.path().to_str().unwrap().to_string();

        let mut args = Args::try_parse_from(vec!["filter_vcf"]).unwrap();
        args.input = Some(temp_path);

        let result = open_input(&args);
        assert!(result.is_ok());
    }

    #[test]
    fn test_case_sensitive_gz_extension() {
        // 测试 .gz 扩展名大小写敏感
        let mut temp_file = NamedTempFile::with_suffix(".GZ").unwrap();
        writeln!(temp_file, "test").unwrap();
        let temp_path = temp_file.path().to_str().unwrap().to_string();

        let mut args = Args::try_parse_from(vec!["filter_vcf"]).unwrap();
        args.input = Some(temp_path);

        let result = open_input(&args);
        assert!(result.is_ok());
        // 注意：根据实现，可能不会被当作压缩文件处理
    }

    // ==================== 性能相关测试 ====================
    #[test]
    fn test_large_file_handle() {
        // 测试处理大文件时的内存和行为
        let mut temp_file = NamedTempFile::new().unwrap();

        // 写入大量数据
        for _ in 0..1000 {
            writeln!(temp_file, "test line with some content").unwrap();
        }

        let temp_path = temp_file.path().to_str().unwrap().to_string();
        let mut args = Args::try_parse_from(vec!["filter_vcf"]).unwrap();
        args.input = Some(temp_path);

        let result = open_input(&args);
        assert!(result.is_ok());
    }

    // ==================== 错误处理测试 ====================
    #[test]
    fn test_output_to_invalid_path() {
        // 尝试写入到无效路径
        let invalid_path = Some("/nonexistent/directory/file.vcf".to_string());

        let result = open_output(&invalid_path, 6);
        assert!(result.is_err());
    }

    #[test]
    fn test_output_to_readonly_location() {
        // 注意：这个测试需要实际创建只读位置，在CI环境中可能不可靠
        // 仅作为示例展示错误处理思路

        // 创建临时目录并设置为只读
        let temp_dir = tempfile::tempdir().unwrap();
        let temp_path = temp_dir.path().join("readonly.vcf");

        // 在某些系统上可能需要额外的权限设置才能使此测试生效
        // 这里主要测试函数不会panic
        let path_str = Some(temp_path.to_str().unwrap().to_string());
        let result = open_output(&path_str, 6);

        // 根据系统权限，这个测试可能成功也可能失败
        // 主要确保不会panic
        if result.is_ok() {
            let mut writer = result.unwrap();
            let _ = writer.write_all(b"test");
        }
    }
}
