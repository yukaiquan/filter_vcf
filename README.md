# VCF Filter & Genotype Calling (Rust Version)
A high-performance Rust reimplementation of the AWK VCF filtering script, with enhanced compatibility, performance and maintainability.  
原AWK VCF过滤脚本的高性能Rust重实现版本，提升了兼容性、性能和可维护性。

## Reference 参考
Original AWK script source: [https://bitbucket.org/ipk_dg_public/vcf_filtering/](https://bitbucket.org/ipk_dg_public/vcf_filtering/)  
原始AWK脚本来源：[https://bitbucket.org/ipk_dg_public/vcf_filtering/](https://bitbucket.org/ipk_dg_public/vcf_filtering/)

## Input File Format 输入文件格式
The program reads VCF/VCF.gz files (supports BGZIP compression) and has the following assumptions:  
该程序读取VCF/VCF.gz文件（支持BGZIP压缩），并做如下适配：

1. The genotype field (FORMAT column) contains DP subfield (dynamic position, no longer fixed to 3rd position);  
基因型字段（FORMAT列）包含DP子字段（动态位置，不再固定为第3位）；
2. The alternative allele depth (DV) is extracted from the AD subfield (format: REF_depth,ALT_depth, takes the 2nd value);  
替代等位基因深度（DV）从AD子字段提取（格式：REF_depth,ALT_depth，取第二个值）；
3. The INFO field (8th column) contains DP (total depth) and MQ (mapping quality) subfields;  
INFO字段（第8列）包含DP（总测序深度）和MQ（映射质量）子字段；
4. Fully compatible with standard VCFv4.0+ format, retains all original header information.  
完全兼容标准VCFv4.0+格式，保留所有原始头文件信息。

## Installation & Compilation 安装与编译
### Prerequisites 前置条件
+ Rust programming environment (1.60+): [https://www.rust-lang.org/tools/install](https://www.rust-lang.org/tools/install)  
Rust编程环境（1.60及以上版本）：[https://www.rust-lang.org/tools/install](https://www.rust-lang.org/tools/install)

### Compilation 编译
```bash
wget https://github.com/yukaiquan/filter_vcf/releases/download/v0.01/filter_vcf
chmod 775 ./filter_vcf

./filter_vcf --help
```



### Compilation 编译
```bash
# Clone or download the source code
# 克隆或下载源代码
git clone https://github.com/yukaiquan/filter_vcf.git # Or copy the main.rs and Cargo.toml files
# 或直接复制main.rs和Cargo.toml文件

# Compile in release mode (optimized for performance)
# 以发布模式编译（优化性能）
cargo build --release

# The compiled binary is located at
# 编译后的可执行文件位于
./target/release/vcf_filter
```

## Usage Example 使用示例
```bash
# Basic usage (same parameter logic as AWK script)
# 基础用法（与AWK脚本参数逻辑一致）
./target/release/vcf_filter \
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
./target/release/vcf_filter -i input.vcf -o output.vcf --min-dp 50 --min-maf 0.05
```

## Parameters 参数说明
| Parameter | English Description | 中文说明 |
| --- | --- | --- |
| `-i/--input` | Input VCF/VCF.gz file path (required) | 输入VCF/VCF.gz文件路径（必填） |
| `-o/--output` | Output VCF/VCF.gz file path (required) | 输出VCF/VCF.gz文件路径（必填） |
| `--dphom` | Minimum read depth to accept a homozygous genotype call (default: 1) | 接受纯合基因型调用的最小读取深度（默认：1） |
| `--dphet` | Minimum read depth to accept a heterozygous genotype call (default: 1) | 接受杂合基因型调用的最小读取深度（默认：1） |
| `--minqual` | Minimum SNP quality (6th column of the VCF) (default: 0.0) | 最低SNP质量值（VCF第6列）（默认：0.0） |
| `--mindp` | Minimum total depth (parsed from DP subfield of INFO field) - SNPs below this value are discarded (default: 0) | 最小总测序深度（从INFO字段的DP子字段解析）- 低于该值的SNP会被丢弃（默认：0） |
| `--minhomn` | Minimum number of homozygous calls (REF/ALT) - SNPs with fewer calls are discarded (default: 0) | 纯合基因型（REF/ALT）的最小调用数 - 低于该值的SNP会被丢弃（默认：0） |
| `--minhomp` | Maximum fraction of heterozygous calls (1 - minhomp) - SNPs exceeding this ratio are discarded (default: 0.0) | 杂合基因型调用的最大比例（1 - minhomp）- 超过该比例的SNP会被丢弃（默认：0.0） |
| `--tol` | Tolerance threshold for allele frequency ratio (DV/DP):   - DV/DP ≤ tol → 0/0   - 0.5-tol ≤ DV/DP ≤ 0.5+tol → 0/1   - DV/DP ≥ 1-tol → 1/1   (default: 0.2499) | 等位基因频率比（DV/DP）的容差阈值：   - DV/DP ≤ tol → 0/0（纯合参考）   - 0.5-tol ≤ DV/DP ≤ 0.5+tol → 0/1（杂合）   - DV/DP ≥ 1-tol → 1/1（纯合替代）   （默认：0.2499） |
| `--minmaf` | Minimum minor allele frequency (MAF) - SNPs below this value are discarded (default: 0.0) | 最小次要等位基因频率（MAF）- 低于该值的SNP会被丢弃（默认：0.0） |
| `--minpresent` | Minimum fraction of present (non-missing) genotypes - SNPs below this ratio are discarded (default: 0.0) | 有效（非缺失）基因型的最小比例 - 低于该比例的SNP会被丢弃（默认：0.0） |
| `--compress-level` | BGZIP compression level (1-9, 6=balanced) (default: 6) | BGZIP压缩级别（1-9，6为平衡值）（默认：6） |


## Core Features 核心特性
### English
1. **High Performance**: Rust-native implementation with buffered I/O and precompiled regex, significantly faster than AWK for large VCF files.
2. **Enhanced Compatibility**: 
    - Supports BGZIP/GZIP compressed input/output (compatible with tabix/bcftools)
    - Dynamically parses FORMAT fields (no fixed position dependency for DP/DV)
    - Automatically extracts DV from AD field (compatible with standard VCF without explicit DV)
3. **Header Preservation**: Retains all original VCF header information, adds filter rule annotation for traceability.
4. **Robust Error Handling**: Skips invalid lines with clear error messages, no program crash.
5. **Full AWK Compatibility**: Exact replication of the original AWK script's filtering and genotype calling logic.

### 中文
1. **高性能**：Rust原生实现，结合缓冲IO和正则预编译，处理大型VCF文件速度远超AWK脚本。
2. **增强兼容性**：
    - 支持BGZIP/GZIP压缩的输入/输出（兼容tabix/bcftools等工具）
    - 动态解析FORMAT字段（不再依赖DP/DV的固定位置）
    - 自动从AD字段提取DV（兼容无显式DV字段的标准VCF）
3. **头文件保留**：完整保留所有原始VCF头文件信息，添加过滤规则注释便于实验追溯。
4. **健壮的错误处理**：跳过无效行并输出清晰的错误提示，程序不会崩溃。
5. **完全兼容AWK**：精准复刻原AWK脚本的过滤和基因型调用逻辑。

## Output Format 输出格式
The output file is a standard VCF/VCF.gz file with:  
输出文件为标准VCF/VCF.gz格式，包含：

1. All original header lines + additional `##FilterRule` annotation (records all filter parameters)  
所有原始头文件行 + 新增`##FilterRule`注释（记录所有过滤参数）
2. INFO field: `DP=total_depth;MQ=mapping_quality` (retains total DP and MQ)  
INFO字段：`DP=总深度;MQ=映射质量`（保留总DP和MQ）
3. FORMAT field: `GT:DP:DV` (standardized genotype format)  
FORMAT字段：`GT:DP:DV`（标准化的基因型格式）
4. Genotype calls: Strictly follows the DV/DP ratio rules (0/0, 0/1, 1/1, ./. )  
基因型调用：严格遵循DV/DP比率规则（0/0、0/1、1/1、./.）

