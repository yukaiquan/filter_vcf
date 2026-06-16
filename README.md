# VCF Filter (Rust Version)

A high-performance Rust reimplementation of the AWK VCF filtering script, with enhanced compatibility, performance and maintainability.

## Reference 参考

M. Mascher (2022-08-17): Filtration script for genetic variant matrices in Variant Call Format (VCF). DOI:10.5447/ipk/2022/15
Original AWK script source: [https://bitbucket.org/ipk_dg_public/vcf_filtering/](https://bitbucket.org/ipk_dg_public/vcf_filtering/)

## Input File Format 输入文件格式

The program reads VCF/VCF.gz files (supports BGZIP compression) and has the following assumptions:  
该程序读取VCF/VCF.gz文件（支持BGZIP压缩），并做如下适配：

1. The genotype field (FORMAT column) contains DP subfield (dynamic position, no longer fixed to 3rd position);  
   基因型字段（FORMAT列）包含DP子字段（动态位置，不再固定为第3位）；
2. Fully compatible with standard VCFv4.0+ format, retains all original header information.  
   完全兼容标准VCFv4.0+格式，保留所有原始头文件信息。
3. Suport stdin/stdout
   支持管道

## Installation & Compilation 安装与编译

### Prerequisites 前置条件

- Rust programming environment (1.60+): [https://www.rust-lang.org/tools/install](https://www.rust-lang.org/tools/install)

### Install 安装

```bash

wget https://github.com/yukaiquan/filter_vcf/releases/download/v0.01/filter_vcf

chmod 775 ./filter_vcf

./filter_vcf --help

```

### Compilation 编译

```bash
# Clone or download the source code
# 克隆或下载源代码
git clone https://github.com/yukaiquan/filter_vcf.git


# Compile in release mode (optimized for performance)
# 以发布模式编译
cargo build --release

# The compiled binary is located at
# 编译后的可执行文件位于
./target/release/vcf_filter
```

## Usage Example 使用示例

```bash
# Basic usage (same parameter logic as AWK script)
# 基础用法（与AWK脚本参数逻辑一致）
./vcf_filter \
  -i input.vcf.gz \
  -o output.vcf.gz \
  --dphom 2 \
  --dphet 4 \
  --minqual 40.0 \
  --mindp 100 \
  --minhomn 1 \
  --minhomp 0.9 \
  --tol 0.2 \
  --minmaf 0.01 \
  --minpresent 0.9 \
  --compress-level 6

# Simplified usage with default parameters
# 使用默认参数的简化用法
./vcf_filter -i input.vcf -o output.vcf --min-dp 50 --min-maf 0.05

# ～v～
bcftools view -V indels -i 'MAC > 3' -m2 -M2 GBS_chr1A.vcf.gz | filter_vcf --dphom 2 --dphet 4 --minpresent 0.5 --minmaf 0.05 | bgzip -c > GBS_chr1A_M2m2_dp2het4miss50maf005.vcf.gz
```

## Parameters 参数说明

| Parameter          | English Description                                                                                                                                   | 中文说明                                                                                                                                                         |
| ------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `-i/--input`       | Input VCF/VCF.gz file path (required)                                                                                                                 | 输入VCF/VCF.gz文件路径（必填）                                                                                                                                   |
| `-o/--output`      | Output VCF/VCF.gz file path (required)                                                                                                                | 输出VCF/VCF.gz文件路径（必填）                                                                                                                                   |
| `--dphom`          | Minimum read depth to accept a homozygous genotype call (default: 1)                                                                                  | 接受纯合基因型调用的最小读取深度（默认：1）                                                                                                                      |
| `--dphet`          | Minimum read depth to accept a heterozygous genotype call (default: 1)                                                                                | 接受杂合基因型调用的最小读取深度（默认：1）                                                                                                                      |
| `--minqual`        | Minimum SNP quality (6th column of the VCF) (default: 0.0)                                                                                            | 最低SNP质量值（VCF第6列）（默认：0.0）                                                                                                                           |
| `--mindp`          | Minimum total depth (parsed from DP subfield of INFO field) - SNPs below this value are discarded (default: 0)                                        | 最小总测序深度（从INFO字段的DP子字段解析）- 低于该值的SNP会被丢弃（默认：0）                                                                                     |
| `--minhomn`        | Minimum number of homozygous calls (REF/ALT) - SNPs with fewer calls are discarded (default: 0)                                                       | 纯合基因型（REF/ALT）的最小调用数 - 低于该值的SNP会被丢弃（默认：0）                                                                                             |
| `--minhomp`        | Maximum fraction of heterozygous calls (1 - minhomp) - SNPs exceeding this ratio are discarded (default: 0.0)                                         | 杂合基因型调用的最大比例（1 - minhomp）- 超过该比例的SNP会被丢弃（默认：0.0）                                                                                    |
| `--tol`            | Tolerance threshold for allele frequency ratio (DV/DP): - DV/DP ≤ tol → 0/0 - 0.5-tol ≤ DV/DP ≤ 0.5+tol → 0/1 - DV/DP ≥ 1-tol → 1/1 (default: 0.2499) | 等位基因频率比（DV/DP）的容差阈值： - DV/DP ≤ tol → 0/0（纯合参考） - 0.5-tol ≤ DV/DP ≤ 0.5+tol → 0/1（杂合） - DV/DP ≥ 1-tol → 1/1（纯合替代） （默认：0.2499） |
| `--minmaf`         | Minimum minor allele frequency (MAF) - SNPs below this value are discarded (default: 0.0)                                                             | 最小次要等位基因频率（MAF）- 低于该值的SNP会被丢弃（默认：0.0）                                                                                                  |
| `--minpresent`     | Minimum fraction of present (non-missing) genotypes - SNPs below this ratio are discarded (default: 0.0)                                              | 有效（非缺失）基因型的最小比例 - 低于该比例的SNP会被丢弃（默认：0.0）                                                                                            |
| `--threads`        | Number of threads for parallel processing (0=auto-use all CPU cores) (default: 0)                                                                     | 并行处理线程数（0=自动使用所有CPU核心）（默认：0）                                                                                                               |
| `--compress-level` | BGZIP compression level (1-9, 6=balanced) (default: 6)                                                                                                | BGZIP压缩级别（1-9，6为平衡值）（默认：6）                                                                                                                       |
