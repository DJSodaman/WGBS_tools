#!/bin/bash
# NCBI WT 专用重跑脚本 — 逐个样本稳健处理
set -euo pipefail
export PATH="/data/home/zhength/miniconda3/envs/wgbs_final/bin:$PATH"

REF="/data/home/zhength/work/2026.6.11_CaMethy_suppleFig/ref"
FNA="${REF}/GCF_000182925.2_NC12_genomic.fna"
OUT="/data/home/zhength/work/2026.6.11_CaMethy_suppleFig/WGBS_tools_new_result_v2"
WT="/data/home/zhength/work/2026.6.11_CaMethy_suppleFig/rawdata_NCBI/wildtype"
PY="/data/home/zhength/miniconda3/envs/wgbs_final/bin/python"
log() { echo "[$(date '+%H:%M:%S')] $*"; }

SAMPLES=(
  "SRR1566116:${WT}/SRR1566116_1.fastq.gz"
  "SRR3106959:${WT}/SRR3106959_1.fastq.gz"
  "SRR3106960:${WT}/SRR3106960_1.fastq.gz"
  "SRR3476867:${WT}/SRR3476867_1.fastq.gz"
  "DRR001162:${WT}/DRR001162_1.fastq.gz"
  "DRR001163:${WT}/DRR001163_1.fastq.gz"
  "DRR001164:${WT}/DRR001164_1.fastq.gz"
  "DRR001165:${WT}/DRR001165_1.fastq.gz"
  "DRR001166:${WT}/DRR001166_1.fastq.gz"
)

for item in "${SAMPLES[@]}"; do
  IFS=':' read -r name fastq <<< "$item"
  SD="$OUT/$name"; BD="$SD/bam"; TD="$SD/trimmed"; ED="$SD/meth_extract"
  mkdir -p "$BD" "$TD" "$ED"
  log "=== $name ==="

  # Step 1: Trim
  TRIMMED="${TD}/${name}_trimmed_1.fq.gz"
  if [ -f "$TRIMMED" ] && [ "$(stat -c%s "$TRIMMED" 2>/dev/null)" -gt 1000000 ]; then
    log "  Trim: done"
  else
    log "  Trim: running..."
    fastp --in1 "$fastq" --out1 "$TRIMMED" \
      --cut_front_window_size=1 --cut_front_mean_quality=3 \
      --cut_tail_window_size=1 --cut_tail_mean_quality=3 \
      --cut_right_window_size=4 --cut_right_mean_quality=15 \
      --length_required=36 --thread 4 \
      --html "${TD}/${name}_fastp.html" --json "${TD}/${name}_fastp.json" 2>&1 | tail -2
  fi

  # Step 2: Bismark (skip if valid BAM exists)
  RAW_BAM=$(find "$BD" -maxdepth 1 -name "*_bismark_bt2.bam" ! -name "*.temp.*" ! -name "*.sorted.*" ! -name "*.dedup*" 2>/dev/null | head -1)
  if [ -n "$RAW_BAM" ] && [ "$(stat -c%s "$RAW_BAM" 2>/dev/null)" -gt 1000000 ]; then
    log "  Bismark: done ($(du -sh "$RAW_BAM" | cut -f1))"
  else
    rm -f "$BD"/*.temp.* "$BD"/*_bismark_bt2* 2>/dev/null
    log "  Bismark: running..."
    bismark --bowtie2 --genome_folder "$REF" --score_min L,0,-0.2 \
      --parallel 5 -p 7 --bam -o "$BD" "$TRIMMED" 2>&1 | tail -2
    RAW_BAM=$(find "$BD" -maxdepth 1 -name "*_bismark_bt2.bam" ! -name "*.temp.*" 2>/dev/null | head -1)
    [ -z "$RAW_BAM" ] && { log "  FAIL: Bismark"; continue; }
  fi

  # Step 3: Extract (dedup not needed for SE — skip and go straight to extract)
  if [ "$(ls "$ED"/*.txt.gz 2>/dev/null | wc -l)" -ge 3 ]; then
    log "  Extract: done"
  else
    log "  Extract: running..."
    bismark_methylation_extractor --multicore 8 --comprehensive --gzip \
      --single-end --ignore 5 --output_dir "$ED" "$RAW_BAM" 2>&1 | tail -2
  fi

  # Step 4: Caller
  MC="$SD/${name}.mC_level_Identification_stat.txt"
  if [ -f "$MC" ]; then
    log "  Caller: done ($(wc -l < "$MC") sites)"
  else
    log "  Caller: running..."
    $PY -c "
from wgbs_tools.modules.caller import run_methylation_calling
run_methylation_calling('$ED', '$SD', '$name', '$FNA', cores=20, min_mc=2, min_depth=5)
" 2>&1 | tail -2
  fi

  # Step 5: Summary
  SF="$SD/${name}.methylation_summary.stat.txt"
  if [ -f "$SF" ]; then
    log "  Summary: done"
  else
    $PY -c "
from wgbs_tools.modules.summary import generate_summary
from wgbs_tools.modules.genome import compute_and_cache_c_stats
c=compute_and_cache_c_stats('$REF','GCF_000182925.2_NC12_genomic.fna')
generate_summary('$MC',c,'$SD','$name')
" 2>&1 | tail -1
  fi
  log "  DONE: $(wc -l < "$MC" 2>/dev/null) sites"
done
log "ALL NCBI WT DONE!"
