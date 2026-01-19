#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# APOE e2/e3/e4 inference from per-sample gVCF using bcftools.
# No new VCF is written.
#
# Correct handling for gVCF SNP lines where ALT is like "C,<NON_REF>".
# Correct severity labeling regardless of genotype order (e4/e3 == e3/e4).
#
# OUTPUT TSV columns:
#   sample  rs429358_alleles rs429358_source rs7412_alleles rs7412_source apoe_genotype severity
#
# USAGE:
#   bash apoe_from_gvcf_bcftools_fixed.sh \
#     --gvcf_dir /path/to/gvcfs \
#     --ref_fa /path/to/GRCh38.fasta \
#     --out apoe_report.tsv
#
# Optional:
#   --pattern "*.g.vcf.gz"
# ============================================================

GVCF_DIR=""
REF_FA=""
OUT="apoe_report.tsv"
PATTERN="*.g.vcf.gz"

# GRCh38 APOE-defining SNPs
RS1_CHR="chr19"; RS1_POS="44908684"   # rs429358 (T>C)
RS2_CHR="chr19"; RS2_POS="44908822"   # rs7412   (C>T)

die() { echo "[ERROR] $*" >&2; exit 1; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    --gvcf_dir) GVCF_DIR="$2"; shift 2 ;;
    --ref_fa)   REF_FA="$2"; shift 2 ;;
    --out)      OUT="$2"; shift 2 ;;
    --pattern)  PATTERN="$2"; shift 2 ;;
    *) die "Unknown argument: $1" ;;
  esac
done

[[ -d "${GVCF_DIR}" ]] || die "--gvcf_dir not found: ${GVCF_DIR}"
[[ -f "${REF_FA}" ]]   || die "--ref_fa not found: ${REF_FA}"
[[ -f "${REF_FA}.fai" ]] || die "Reference FASTA index missing: ${REF_FA}.fai (run: samtools faidx ${REF_FA})"

command -v bcftools >/dev/null 2>&1 || die "bcftools not in PATH"
command -v samtools >/dev/null 2>&1 || die "samtools not in PATH"

detect_vcf_chr_style() {
  local vcf="$1"
  local contig
  contig="$(bcftools view -h "$vcf" | awk -F'[=,>]' '/^##contig=<ID=/{print $3; exit}')"
  if [[ -n "$contig" && "$contig" == chr* ]]; then
    echo "chr"
  else
    echo "nochr"
  fi
}

norm_chr() {
  local chr="$1" style="$2"
  if [[ "$style" == "chr" ]]; then
    [[ "$chr" == chr* ]] && echo "$chr" || echo "chr${chr}"
  else
    [[ "$chr" == chr* ]] && echo "${chr#chr}" || echo "$chr"
  fi
}

fa_style_from_fai() {
  awk 'NR==1{print ($1 ~ /^chr/ ? "chr" : "nochr"); exit}' "${REF_FA}.fai"
}

fa_base() {
  local chr="$1" pos="$2"
  samtools faidx "$REF_FA" "${chr}:${pos}-${pos}" | awk 'NR==2{print toupper($0)}'
}

# ------------------------------------------------------------
# site_call:
#  - Prefer explicit record at POS where ALT is NOT exactly "<NON_REF>"
#    (ALT may be "C,<NON_REF>" which is a true SNP line)
#  - Else if a <NON_REF> block covers POS (END>=POS), treat as 0/0 and use FASTA base
# Returns: "<A1>/<A2>\t<SOURCE>"
# SOURCE: SITE | BLOCK | MISSING | SITE_NO_GT | BLOCK_NO_REF
# ------------------------------------------------------------
site_call() {
  local vcf="$1" chr="$2" pos="$3" vcf_style="$4"

  local chr_vcf
  chr_vcf="$(norm_chr "$chr" "$vcf_style")"

  local recs
  recs="$(bcftools view -r "${chr_vcf}:${pos}-${pos}" -H "$vcf" 2>/dev/null || true)"
  if [[ -z "$recs" ]]; then
    echo -e "./.\tMISSING"
    return 0
  fi

  # Explicit record: POS==pos and ALT != "<NON_REF>"
  # This includes ALT like "C,<NON_REF>".
  local explicit
  explicit="$(echo "$recs" | awk -v p="$pos" '
    ($2==p && $5 != "<NON_REF>") {print; exit}
  ')"

  if [[ -n "$explicit" ]]; then
    echo "$explicit" | awk -v OFS="\t" '
      function idx(fmt, key,   n,i,a){
        n=split(fmt,a,":");
        for(i=1;i<=n;i++) if(a[i]==key) return i;
        return 0;
      }
      {
        ref=toupper($4);
        alt=toupper($5);
        fmt=$9;
        samp=$10;

        gi=idx(fmt,"GT");
        if(gi==0){ print "./.","SITE_NO_GT"; exit }

        n=split(samp,sv,":");
        gt=sv[gi];
        gsub(/[|]/,"/",gt);

        if(gt=="." || gt=="./."){ print "./.","SITE_NO_GT"; exit }

        # ALT list can be "C,<NON_REF>" etc.
        na=split(alt,alts,",");

        split(gt,gta,"/");
        if(gta[1]=="" || gta[2]==""){ print "./.","SITE_NO_GT"; exit }

        # allele index: 0=REF, 1..=ALT(s) in order
        a1 = (gta[1]==0)?ref:alts[gta[1]];
        a2 = (gta[2]==0)?ref:alts[gta[2]];

        if(a1=="" || a2==""){ print "./.","SITE_NO_GT"; exit }

        print a1"/"a2,"SITE";
      }'
    return 0
  fi

  # Otherwise, look for a <NON_REF> block that covers POS using END=
  local block
  block="$(echo "$recs" | awk -v p="$pos" '
    ($5 == "<NON_REF>") {
      end=$2;
      if (match($8, /END=[0-9]+/)) end=substr($8, RSTART+4, RLENGTH-4)+0;
      if ($2 <= p && end >= p) {print; exit}
    }
  ')"

  if [[ -n "$block" ]]; then
    local fa_style fa_chr base
    fa_style="$(fa_style_from_fai)"
    fa_chr="$(norm_chr "$chr" "$fa_style")"
    base="$(fa_base "$fa_chr" "$pos")"
    if [[ -z "$base" ]]; then
      echo -e "./.\tBLOCK_NO_REF"
    else
      echo -e "${base}/${base}\tBLOCK"
    fi
    return 0
  fi

  echo -e "./.\tMISSING"
}

infer_apoe() {
  local a429="$1" a741="$2"
  [[ "$a429" != "./." && "$a741" != "./." ]] || { echo "NA"; return 0; }

  local a1 a2 b1 b2
  a1="${a429%/*}"; a2="${a429#*/}"
  b1="${a741%/*}"; b2="${a741#*/}"

  hap() {
    local x="$1" y="$2"
    if [[ "$x" == "T" && "$y" == "T" ]]; then echo "e2"; return; fi
    if [[ "$x" == "T" && "$y" == "C" ]]; then echo "e3"; return; fi
    if [[ "$x" == "C" && "$y" == "C" ]]; then echo "e4"; return; fi
    echo ""
  }

  local h1 h2 g1 g2
  h1="$(hap "$a1" "$b1")"; h2="$(hap "$a2" "$b2")"
  if [[ -n "$h1" && -n "$h2" ]]; then
    # normalize ordering with e4 first if present, else lexical
    if [[ "$h1" == "e4" && "$h2" != "e4" ]]; then echo "e4/$h2"; return; fi
    if [[ "$h2" == "e4" && "$h1" != "e4" ]]; then echo "e4/$h1"; return; fi
    [[ "$h1" < "$h2" ]] && echo "$h1/$h2" || echo "$h2/$h1"
    return
  fi

  g1="$(hap "$a1" "$b2")"; g2="$(hap "$a2" "$b1")"
  if [[ -n "$g1" && -n "$g2" ]]; then
    if [[ "$g1" == "e4" && "$g2" != "e4" ]]; then echo "e4/$g2"; return; fi
    if [[ "$g2" == "e4" && "$g1" != "e4" ]]; then echo "e4/$g1"; return; fi
    [[ "$g1" < "$g2" ]] && echo "$g1/$g2" || echo "$g2/$g1"
    return
  fi

  echo "AMBIG"
}

# Order-safe severity labeling (e4/e3 == e3/e4, e4/e2 == e2/e4, etc.)
severity_label() {
  local g="$1"
  case "$g" in
    e4/e4)       echo "Very_high_risk" ;;
    e4/e3|e3/e4) echo "High_risk" ;;
    e2/e4|e4/e2) echo "Moderate_mixed" ;;
    e2/e3|e3/e2) echo "Mild_protective" ;;
    e3/e3)       echo "Baseline" ;;
    e2/e2)       echo "Protective" ;;
    AMBIG)       echo "Ambiguous" ;;
    NA)          echo "Missing" ;;
    *)           echo "Unknown" ;;
  esac
}

# ---------------------- MAIN ----------------------
shopt -s nullglob
GVCFS=( "${GVCF_DIR}"/${PATTERN} )
[[ ${#GVCFS[@]} -gt 0 ]] || die "No gVCFs matched: ${GVCF_DIR}/${PATTERN}"

echo -e "sample\trs429358_alleles\trs429358_source\trs7412_alleles\trs7412_source\tapoe_genotype\tseverity" > "${OUT}"

for vcf in "${GVCFS[@]}"; do
  base="$(basename "$vcf")"
  sample="${base%.g.vcf.gz}"   # keep version in name (e.g., GARDWGSN00376_v1.1.0)

  style="$(detect_vcf_chr_style "$vcf")"

  c1="$(site_call "$vcf" "$RS1_CHR" "$RS1_POS" "$style")"
  a429="$(echo "$c1" | awk '{print $1}')"
  s429="$(echo "$c1" | awk '{print $2}')"

  c2="$(site_call "$vcf" "$RS2_CHR" "$RS2_POS" "$style")"
  a741="$(echo "$c2" | awk '{print $1}')"
  s741="$(echo "$c2" | awk '{print $2}')"

  apoe="$(infer_apoe "$a429" "$a741")"
  sev="$(severity_label "$apoe")"

  echo -e "${sample}\t${a429}\t${s429}\t${a741}\t${s741}\t${apoe}\t${sev}" >> "${OUT}"
done

echo "[OK] Wrote: ${OUT}"
