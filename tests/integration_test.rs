use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::process::{Command, Stdio};
use tempfile::NamedTempFile;

/// 创建测试用的VCF文件
fn create_test_vcf() -> NamedTempFile {
    let mut temp_file = NamedTempFile::with_suffix(".vcf").unwrap();
    
    let vcf_content = r#"##fileformat=VCFv4.2
##fileDate=20260602
##source=example
##reference=GRCh38
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2	Sample3
chr1	100	.	A	T	50	PASS	DP=100	GT:DP:AD	0/0:50:50,0	1/1:50:0,50	0/1:50:25,25
chr1	200	.	C	G	30	PASS	DP=80	GT:DP:AD	0/0:40:40,0	0/0:40:40,0	0/0:40:40,0
chr1	300	.	G	A	60	PASS	DP=120	GT:DP:AD	1/1:60:0,60	0/1:60:30,30	0/0:60:60,0
chr2	100	.	N	A	50	PASS	DP=100	GT:DP:AD	0/0:50:50,0	1/1:50:0,50	0/1:50:25,25
chr2	200	.	TG	T	50	PASS	DP=100	GT:DP:AD	0/0:50:50,0	1/1:50:0,50	0/1:50:25,25
chr2	300	.	A	T	20	PASS	DP=50	GT:DP:AD	0/0:25:25,0	1/1:25:0,25	0/1:25:12,13
"#;
    
    write!(temp_file, "{}", vcf_content).unwrap();
    temp_file
}

/// 测试基本的VCF过滤功能
#[test]
fn test_basic_filtering() {
    let input_file = create_test_vcf();
    let output_file = NamedTempFile::with_suffix(".vcf").unwrap();
    
    let status = Command::new("./target/debug/filter_vcf")
        .args([
            "-i", input_file.path().to_str().unwrap(),
            "-o", output_file.path().to_str().unwrap(),
        ])
        .status()
        .expect("Failed to execute command");
    
    assert!(status.success());
    
    // 验证输出文件存在且有内容
    let output_metadata = std::fs::metadata(output_file.path()).unwrap();
    assert!(output_metadata.len() > 0);
    
    // 读取输出文件，验证格式正确
    let file = File::open(output_file.path()).unwrap();
    let reader = BufReader::new(file);
    let lines: Vec<String> = reader.lines().map(|l| l.unwrap()).collect();
    
    // 应该包含header行
    assert!(lines.iter().any(|l| l.starts_with("##fileformat=VCFv4.2")));
    assert!(lines.iter().any(|l| l.starts_with("#CHROM")));
    
    // 应该包含过滤规则注释
    assert!(lines.iter().any(|l| l.contains("FilterRule")));
}

/// 测试最小质量值过滤
#[test]
fn test_minqual_filtering() {
    let input_file = create_test_vcf();
    let output_file = NamedTempFile::with_suffix(".vcf").unwrap();
    
    let status = Command::new("./target/debug/filter_vcf")
        .args([
            "-i", input_file.path().to_str().unwrap(),
            "-o", output_file.path().to_str().unwrap(),
            "--minqual", "40",
        ])
        .status()
        .expect("Failed to execute command");
    
    assert!(status.success());
    
    // 验证输出文件
    let file = File::open(output_file.path()).unwrap();
    let reader = BufReader::new(file);
    let data_lines: Vec<String> = reader
        .lines()
        .map(|l| l.unwrap())
        .filter(|l| !l.starts_with('#'))
        .collect();
    
    // 只有QUAL >= 40的SNP应该保留
    for line in data_lines {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() > 5 {
            let qual: f64 = parts[5].parse().unwrap();
            assert!(qual >= 40.0);
        }
    }
}

/// 测试最小DP过滤
#[test]
fn test_mindp_filtering() {
    let input_file = create_test_vcf();
    let output_file = NamedTempFile::with_suffix(".vcf").unwrap();
    
    let status = Command::new("./target/debug/filter_vcf")
        .args([
            "-i", input_file.path().to_str().unwrap(),
            "-o", output_file.path().to_str().unwrap(),
            "--mindp", "100",
        ])
        .status()
        .expect("Failed to execute command");
    
    assert!(status.success());
    
    // 验证输出文件
    let file = File::open(output_file.path()).unwrap();
    let reader = BufReader::new(file);
    let data_lines: Vec<String> = reader
        .lines()
        .map(|l| l.unwrap())
        .filter(|l| !l.starts_with('#'))
        .collect();
    
    // 只有INFO DP >= 100的SNP应该保留
    for line in data_lines {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() > 7 {
            let info = parts[7];
            // 提取DP值
            if let Some(dp_str) = info.split("DP=").nth(1) {
                if let Some(dp_value) = dp_str.split(';').next() {
                    let dp: u32 = dp_value.parse().unwrap_or(0);
                    assert!(dp >= 100);
                }
            }
        }
    }
}

/// 测试indel过滤
#[test]
fn test_indel_filtering() {
    let input_file = create_test_vcf();
    let output_file = NamedTempFile::with_suffix(".vcf").unwrap();
    
    // 测试不包含indel（默认）
    let status = Command::new("./target/debug/filter_vcf")
        .args([
            "-i", input_file.path().to_str().unwrap(),
            "-o", output_file.path().to_str().unwrap(),
        ])
        .status()
        .expect("Failed to execute command");
    
    assert!(status.success());
    
    // 验证输出文件不包含indel
    let file = File::open(output_file.path()).unwrap();
    let reader = BufReader::new(file);
    let data_lines: Vec<String> = reader
        .lines()
        .map(|l| l.unwrap())
        .filter(|l| !l.starts_with('#'))
        .collect();
    
    for line in &data_lines {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() > 4 {
            let ref_base = parts[3];
            let alt_base = parts[4];
            // REF和ALT应该都是单个碱基
            assert_eq!(ref_base.len(), 1);
            assert_eq!(alt_base.len(), 1);
        }
    }
}

/// 测试包含indel
#[test]
fn test_include_indel() {
    let input_file = create_test_vcf();
    let output_file = NamedTempFile::with_suffix(".vcf").unwrap();
    
    let status = Command::new("./target/debug/filter_vcf")
        .args([
            "-i", input_file.path().to_str().unwrap(),
            "-o", output_file.path().to_str().unwrap(),
            "--include-indel",
        ])
        .status()
        .expect("Failed to execute command");
    
    assert!(status.success());
    
    // 验证输出文件包含indel
    let file = File::open(output_file.path()).unwrap();
    let reader = BufReader::new(file);
    let data_lines: Vec<String> = reader
        .lines()
        .map(|l| l.unwrap())
        .filter(|l| !l.starts_with('#'))
        .collect();
    
    // 应该有包含indel的行（TG -> T）
    let has_indel = data_lines.iter().any(|line| {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() > 4 {
            parts[3].len() > 1 || parts[4].len() > 1
        } else {
            false
        }
    });
    
    assert!(has_indel);
}

/// 测试MAF过滤
#[test]
fn test_minmaf_filtering() {
    let input_file = create_test_vcf();
    let output_file = NamedTempFile::with_suffix(".vcf").unwrap();
    
    let status = Command::new("./target/debug/filter_vcf")
        .args([
            "-i", input_file.path().to_str().unwrap(),
            "-o", output_file.path().to_str().unwrap(),
            "--minmaf", "0.1",
        ])
        .status()
        .expect("Failed to execute command");
    
    assert!(status.success());
    
    // 验证输出文件存在
    let output_metadata = std::fs::metadata(output_file.path()).unwrap();
    assert!(output_metadata.len() > 0);
}

/// 测试管道模式（stdin到stdout）
#[test]
fn test_pipeline_mode() {
    let input_file = create_test_vcf();
    
    // 使用管道模式：从文件读取，输出到stdout
    let mut child = Command::new("./target/debug/filter_vcf")
        .args(["-i", input_file.path().to_str().unwrap()])
        .stdout(Stdio::piped())
        .spawn()
        .expect("Failed to execute command");
    
    let output = child.wait_with_output().expect("Failed to read stdout");
    
    assert!(output.status.success());
    assert!(!output.stdout.is_empty());
    
    // 验证输出是有效的VCF格式
    let output_str = String::from_utf8(output.stdout).unwrap();
    assert!(output_str.contains("##fileformat=VCFv4.2"));
    assert!(output_str.contains("#CHROM"));
}

/// 测试压缩输出
#[test]
fn test_gz_output() {
    let input_file = create_test_vcf();
    let output_file = NamedTempFile::with_suffix(".vcf.gz").unwrap();
    
    let status = Command::new("./target/debug/filter_vcf")
        .args([
            "-i", input_file.path().to_str().unwrap(),
            "-o", output_file.path().to_str().unwrap(),
        ])
        .status()
        .expect("Failed to execute command");
    
    assert!(status.success());
    
    // 验证输出文件存在
    let output_metadata = std::fs::metadata(output_file.path()).unwrap();
    assert!(output_metadata.len() > 0);
    
    // 验证文件是gzip格式（通过魔数）
    let mut file = File::open(output_file.path()).unwrap();
    let mut buffer = [0u8; 2];
    std::io::Read::read_exact(&mut file, &mut buffer).unwrap();
    
    // Gzip魔数: 0x1f 0x8b
    assert_eq!(buffer[0], 0x1f);
    assert_eq!(buffer[1], 0x8b);
}

/// 测试多个参数组合
#[test]
fn test_combined_filters() {
    let input_file = create_test_vcf();
    let output_file = NamedTempFile::with_suffix(".vcf").unwrap();
    
    let status = Command::new("./target/debug/filter_vcf")
        .args([
            "-i", input_file.path().to_str().unwrap(),
            "-o", output_file.path().to_str().unwrap(),
            "--minqual", "30",
            "--mindp", "80",
            "--minmaf", "0.05",
            "--dphom", "2",
            "--dphet", "4",
        ])
        .status()
        .expect("Failed to execute command");
    
    assert!(status.success());
    
    // 验证输出文件
    let output_metadata = std::fs::metadata(output_file.path()).unwrap();
    assert!(output_metadata.len() > 0);
}

/// 测试REF=N过滤
#[test]
fn test_ref_n_filtering() {
    let input_file = create_test_vcf();
    let output_file = NamedTempFile::with_suffix(".vcf").unwrap();
    
    let status = Command::new("./target/debug/filter_vcf")
        .args([
            "-i", input_file.path().to_str().unwrap(),
            "-o", output_file.path().to_str().unwrap(),
        ])
        .status()
        .expect("Failed to execute command");
    
    assert!(status.success());
    
    // 验证输出文件不包含REF=N的SNP
    let file = File::open(output_file.path()).unwrap();
    let reader = BufReader::new(file);
    let data_lines: Vec<String> = reader
        .lines()
        .map(|l| l.unwrap())
        .filter(|l| !l.starts_with('#'))
        .collect();
    
    for line in data_lines {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() > 3 {
            let ref_base = parts[3];
            assert_ne!(ref_base, "N");
        }
    }
}

/// 测试文件格式保留
#[test]
fn test_format_preservation() {
    let input_file = create_test_vcf();
    let output_file = NamedTempFile::with_suffix(".vcf").unwrap();
    
    let status = Command::new("./target/debug/filter_vcf")
        .args([
            "-i", input_file.path().to_str().unwrap(),
            "-o", output_file.path().to_str().unwrap(),
        ])
        .status()
        .expect("Failed to execute command");
    
    assert!(status.success());
    
    let file = File::open(output_file.path()).unwrap();
    let reader = BufReader::new(file);
    let lines: Vec<String> = reader.lines().map(|l| l.unwrap()).collect();
    
    // 验证header格式
    assert!(lines.iter().any(|l| l.starts_with("##fileformat=VCFv4.2")));
    assert!(lines.iter().any(|l| l.starts_with("##fileDate=")));
    assert!(lines.iter().any(|l| l.starts_with("##source=")));
    assert!(lines.iter().any(|l| l.starts_with("##reference=")));
    
    // 验证#CHROM行格式
    let chrom_header = lines.iter().find(|l| l.starts_with("#CHROM")).unwrap();
    let parts: Vec<&str> = chrom_header.split('\t').collect();
    assert!(parts.len() >= 10);
    assert_eq!(parts[0], "#CHROM");
    assert_eq!(parts[1], "POS");
    assert_eq!(parts[2], "ID");
    assert_eq!(parts[3], "REF");
    assert_eq!(parts[4], "ALT");
}
