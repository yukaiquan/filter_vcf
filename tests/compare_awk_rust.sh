#!/bin/bash
# 对比 AWK 脚本与 Rust 程序在不同过滤参数下的输出差异
# 用法: bash tests/compare_awk_rust.sh
# 前置: cargo build --release 已完成

set -u
cd "$(dirname "$0")/.."

VCF=tests/test.vcf.gz
RUST=./target/release/filter_vcf
AWK=awk
AWK_SCRIPT=tests/filter_vcf_bsd.awk
TMP=/tmp/cmp_vcf
mkdir -p "$TMP"

# 检查文件
[ -f "$VCF" ] || { echo "缺少 $VCF"; exit 1; }
[ -x "$RUST" ] || { echo "缺少 $RUST，请先 cargo build --release"; exit 1; }
[ -f "$AWK_SCRIPT" ] || { echo "缺少 $AWK_SCRIPT"; exit 1; }

# 对比函数：参数1=参数描述，后续参数=传给两者的参数（AWK 用 -v，Rust 用 --）
# 这里我们手动指定每种参数组合
compare() {
    local desc="$1"; shift
    local awk_vars=()
    local rust_args=()
    local k v
    # 解析 key=val 对
    while [ $# -gt 0 ]; do
        local kv="$1"; shift
        local k="${kv%%=*}"
        local v="${kv#*=}"
        awk_vars+=("-v" "$k=$v")
        rust_args+=("--$k" "$v")
    done

    echo "=============================="
    echo "参数组合: $desc"
    echo "  AWK vars: ${awk_vars[*]:-（默认）}"
    echo "  Rust args: ${rust_args[*]:-（默认）}"

    # 运行 AWK（输出 GT:DP:DV 格式）
    if [ ${#awk_vars[@]} -gt 0 ]; then
        zcat "$VCF" | $AWK -f "$AWK_SCRIPT" "${awk_vars[@]}" > "$TMP/awk.out" 2>"$TMP/awk.err"
    else
        zcat "$VCF" | $AWK -f "$AWK_SCRIPT" > "$TMP/awk.out" 2>"$TMP/awk.err"
    fi
    local awk_rc=$?
    local awk_lines=$(grep -vc "^#" "$TMP/awk.out" | tr -d ' ')

    # 运行 Rust（输出到 stdout，保留原 FORMAT）
    if [ ${#rust_args[@]} -gt 0 ]; then
        zcat "$VCF" | $RUST "${rust_args[@]}" > "$TMP/rust.out" 2>"$TMP/rust.err"
    else
        zcat "$VCF" | $RUST > "$TMP/rust.out" 2>"$TMP/rust.err"
    fi
    local rust_rc=$?
    local rust_lines=$(grep -vc "^#" "$TMP/rust.out" | tr -d ' ')

    echo "  AWK : rc=$awk_rc 数据行=$awk_lines"
    echo "  Rust: rc=$rust_rc 数据行=$rust_lines"

    if [ "$awk_rc" -ne 0 ]; then
        echo "  [AWK 报错]"; head -3 "$TMP/awk.err"
        return
    fi
    if [ "$rust_rc" -ne 0 ]; then
        echo "  [Rust 报错]"; head -3 "$TMP/rust.err"
        return
    fi

    if [ "$awk_lines" -eq "$rust_lines" ]; then
        echo "  ✅ 行数一致 ($awk_lines)"
    else
        echo "  ❌ 行数不一致: AWK=$awk_lines Rust=$rust_lines (差 $((rust_lines-awk_lines)))"
    fi

    # 提取两者的样本 GT 矩阵做逐位对比（只比较 GT 字段，忽略 DP/DV/INFO/FORMAT 差异）
    # AWK 输出: CHROM POS ID REF ALT QUAL FILTER INFO FORMAT sample... 样本格式 GT:DP:DV
    # Rust 输出: 保留原列，样本 GT 在 FORMAT 指定位置
    # 统一提取: 染色体+位置 作为 key，每个样本的 GT 作为 value
    extract_gt_awk() {
        awk -F'\t' '
        /^#/ {next}
        {
            key=$1":"$2
            gt_str=""
            for(i=10;i<=NF;i++){
                split($i,a,":")
                gt_str=gt_str (i>10?",":"") a[1]
            }
            print key"\t"gt_str
        }' "$1"
    }
    extract_gt_rust() {
        # Rust 保留原 FORMAT，GT 位置需动态找；这里假设 GT 是 FORMAT 第一个字段（常见）
        # 更稳妥：用 Rust 输出本身的 FORMAT 列解析
        awk -F'\t' '
        NR==1 || /^#/ {next}
        {
            key=$1":"$2
            # 找 FORMAT 中 GT 的位置
            fmt=$9
            n=split(fmt,fa,":")
            gtpos=1
            for(j=1;j<=n;j++) if(fa[j]=="GT") gtpos=j
            gt_str=""
            for(i=10;i<=NF;i++){
                split($i,a,":")
                gt_str=gt_str (i>10?",":"") a[gtpos]
            }
            print key"\t"gt_str
        }' "$1"
    }

    extract_gt_awk "$TMP/awk.out" | sort > "$TMP/awk.gt"
    extract_gt_rust "$TMP/rust.out" | sort > "$TMP/rust.gt"

    local diff_out
    diff_out=$(diff "$TMP/awk.gt" "$TMP/rust.gt" | head -20)
    if [ -z "$diff_out" ]; then
        echo "  ✅ 所有位点 GT 完全一致"
    else
        local diff_cnt
        diff_cnt=$(diff "$TMP/awk.gt" "$TMP/rust.gt" | grep -c "^[<>]")
        echo "  ⚠️  GT 差异行数: $diff_cnt (双向)"
        echo "$diff_out" | head -10 | sed 's/^/      /'
    fi
    echo ""
}

echo "========================================"
echo "AWK vs Rust 过滤参数对比测试"
echo "VCF: $VCF ($(zcat $VCF | grep -vc '^#') 条记录, $(zcat $VCF | grep '^#CHROM' | awk -F'\t' '{print NF-9}') 个样本)"
echo "========================================"
echo ""

# 测试组1: 默认参数
compare "默认参数 (dphom=1,dphet=1,tol=0.2499)"

# 测试组2: 提高纯合深度阈值
compare "dphom=3" dphom=3

# 测试组3: 提高杂合深度阈值
compare "dphet=5" dphet=5

# 测试组4: 收紧 tol
compare "tol=0.15" tol=0.15

# 测试组5: minqual 过滤
compare "minqual=30" minqual=30

# 测试组6: mindp 过滤
compare "mindp=50" mindp=50

# 测试组7: minhomn 过滤
compare "minhomn=2" minhomn=2

# 测试组8: minpresent 过滤
compare "minpresent=0.8" minpresent=0.8

# 测试组9: minhomp 过滤
compare "minhomp=0.9" minhomp=0.9

# 测试组10: minmaf 过滤
compare "minmaf=0.05" minmaf=0.05

# 测试组11: 组合参数
compare "组合: dphom=2 dphet=3 minqual=20 minmaf=0.01" dphom=2 dphet=3 minqual=20 minmaf=0.01

# 测试组12: 严格组合
compare "严格: dphom=5 dphet=5 tol=0.1 minhomn=3 minpresent=0.9 minhomp=0.95 minmaf=0.1" dphom=5 dphet=5 tol=0.1 minhomn=3 minpresent=0.9 minhomp=0.95 minmaf=0.1

echo "========================================"
echo "对比测试完成"
echo "========================================"
