#!/usr/bin/awk -f
# BSD awk 兼容版：用 match+substr 替代 gensub，逻辑与 filter_vcf.awk 完全一致

BEGIN{
 OFS=FS="\t"

 print "##fileformat=VCFv4.0"
 print "##FILTER=<ID=PASS,Description=\"All filters passed\">"
 print "##INFO=<ID=MQ,Number=1,Type=Integer,Description=\"Average mapping quality\">"
 print "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
 print "##FORMAT=<ID=DP,Number=.,Type=Integer,Description=\"Read depth\">"
 print "##FORMAT=<ID=DV,Number=.,Type=Integer,Description=\"Read depth of the alternative allele\">"

 if(!dphom)
  dphom=1
 if(!dphet)
  dphet=1
 if(!tol)
  tol=0.2499
}

/^#(#reference=|CHROM)/ {
 print; next
}

/^#/ || length($4) > 1 || length($5) > 1 || $6 < minqual || $4 == "N" || 0+extract_dp($8) < mindp {
 next
}

{
 o=""
 n=A=B=H=0

 for(i = 10; i <= NF; i++){
  split($i, a, ":")
  dp = 0+a[2]
  dv = 0+a[3]
  if(dp > 0)
   r=dv/dp
  else if($5 == "."){
   o = o"\t0/0:"dp":"dv
   A++
  }
  if(dp >= dphom && r <= tol ){
   o = o"\t0/0:"dp":"dv
   A++
  }
  else if(dp >= dphom && r >= 1-tol){
   o = o"\t1/1:"dp":"dv
   B++
  }
  else if(dp >= dphet && r >= 0.5-tol && r <= 0.5+tol){
   o = o"\t0/1:"dp":"dv
   H++
  }
  else{
   o = o"\t./.:"dp":"dv
   n++
  }
 }
 present = A + B + H
}

!present || A < minhomn || B < minhomn || present < minpresent * (present + n) || A + B < minhomp * present {
 next
}

{
 if(B > A)
  m = (2*A + H) / 2 / present
 else
  m = (2*B + H) / 2 / present
}

m >= minmaf {
 mq = extract_mq($8)
 print $1, $2, $3, $4, $5, $6, $7, "MQ="mq, "GT:DP:DV"o
}

# 从 INFO 字段提取 DP=xxx 的数值（替代 gensub）
function extract_dp(s,    m, p) {
 p = index(s, "DP=")
 if(p == 0) return 0
 m = substr(s, p+3)
 # 取连续数字
 return 0+extract_num(m)
}

# 从 INFO 字段提取 MQ=xxx 的数值（替代 gensub）
function extract_mq(s,    p, m) {
 p = index(s, "MQ=")
 if(p == 0) return 0
 m = substr(s, p+3)
 return 0+extract_num(m)
}

# 从字符串起始提取连续数字
function extract_num(s,    i, c, n) {
 n = ""
 for(i = 1; i <= length(s); i++){
  c = substr(s, i, 1)
  if(c >= "0" && c <= "9")
   n = n c
  else
   break
 }
 return (n == "" ? 0 : n)
}
