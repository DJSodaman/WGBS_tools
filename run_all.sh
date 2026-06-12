#!/bin/bash
# WGBS_tools 批量运行脚本 (72核)
# 处理: NC2M_1/2/3 + NCBI WT (9) = 12 样本
# 用法: nohup bash run_all.sh > run_all.log 2>&1 &
set -euo pipefail

source /data/home/zhength/miniconda3/etc/profile.d/conda.sh
conda activate /data/home/zhength/miniconda3/envs/wgbs_final
export PATH="/data/home/zhength/miniconda3/envs/wgbs_final/bin:$PATH"

# ========== 配置 ==========
REF_DIR="/data/home/zhength/work/2026.6.11_CaMethy_suppleFig/ref"
FASTA="${REF_DIR}/GCF_000182925.2_NC12_genomic.fna"
OUT_DIR="/data/home/zhength/work/2026.6.11_CaMethy_suppleFig/WGBS_tools_new_result"
RAWDATA_USER="/data/home/zhength/work/2026.6.11_CaMethy_suppleFig/rawdata"
RAWDATA_NCBI="/data/home/zhength/work/2026.6.11_CaMethy_suppleFig/rawdata_NCBI"

BISMARK_PARALLEL=3
BISMARK_BT2_THREADS=7
METH_EXTRACT_CORES=8
FASTP_THREADS=4
CALLER_CORES=40

# ========== 工具函数 ==========
log() { echo "[$(date '+%m-%d %H:%M:%S')] $*"; }

# 安全文件查找 (兼容 pipefail)
safe_find() { find "$1" -maxdepth 1 -name "$2" 2>/dev/null | head -1 || true; }
safe_find_count() { find "$1" -maxdepth 1 -name "*.txt.gz" 2>/dev/null | wc -l || echo 0; }

run_fastp_se() {
    local r1=$1 outdir=$2 name=$3
    mkdir -p "$outdir"
    local out1="${outdir}/${name}_trimmed_1.fq.gz"
    if [ -f "$out1" ]; then log "   修剪已完成, 跳过"; echo "$out1"; return; fi
    log "   Fastp SE: $name"
    fastp --in1 "$r1" --out1 "$out1" \
        --cut_front_window_size=1 --cut_front_mean_quality=3 \
        --cut_tail_window_size=1 --cut_tail_mean_quality=3 \
        --cut_right_window_size=4 --cut_right_mean_quality=15 \
        --length_required=36 --thread $FASTP_THREADS \
        --html "${outdir}/${name}_fastp.html" --json "${outdir}/${name}_fastp.json" 2>&1 | tail -3 || true
    echo "$out1"
}

# ========== 主流程 ==========
log "=========================================="
log "WGBS_tools 批量分析 (72核, 12样本)"
log "=========================================="

# 基因组C统计
log "基因组C统计..."
python -c "
from wgbs_tools.modules.genome import compute_and_cache_c_stats
s = compute_and_cache_c_stats('$REF_DIR', 'GCF_000182925.2_NC12_genomic.fna')
print(f'C: total={s[\"totals\"][\"total_C\"]:,} CG={s[\"totals\"][\"C_CG\"]:,}')
" || { log "基因组统计失败!"; exit 1; }

# 样本列表
SAMPLES=(
    "NC2M_1:pe:${RAWDATA_USER}/NC2M_1_clean_1.fq.gz:${RAWDATA_USER}/NC2M_1_clean_2.fq.gz:skip_trim"
    "NC2M_2:pe:${RAWDATA_USER}/NC2M_2_clean_1.fq.gz:${RAWDATA_USER}/NC2M_2_clean_2.fq.gz:skip_trim"
    "NC2M_3:pe:${RAWDATA_USER}/NC2M_3_clean_1.fq.gz:${RAWDATA_USER}/NC2M_3_clean_2.fq.gz:skip_trim"
    "SRR1566116:se:${RAWDATA_NCBI}/wildtype/SRR1566116_1.fastq.gz::trim"
    "SRR3106959:se:${RAWDATA_NCBI}/wildtype/SRR3106959_1.fastq.gz::trim"
    "SRR3106960:se:${RAWDATA_NCBI}/wildtype/SRR3106960_1.fastq.gz::trim"
    "SRR3476867:se:${RAWDATA_NCBI}/wildtype/SRR3476867_1.fastq.gz::trim"
    "DRR001162:se:${RAWDATA_NCBI}/wildtype/DRR001162_1.fastq.gz::trim"
    "DRR001163:se:${RAWDATA_NCBI}/wildtype/DRR001163_1.fastq.gz::trim"
    "DRR001164:se:${RAWDATA_NCBI}/wildtype/DRR001164_1.fastq.gz::trim"
    "DRR001165:se:${RAWDATA_NCBI}/wildtype/DRR001165_1.fastq.gz::trim"
    "DRR001166:se:${RAWDATA_NCBI}/wildtype/DRR001166_1.fastq.gz::trim"
)

TOTAL=${#SAMPLES[@]}
PASSED=0; FAILED=0
START_TIME=$(date +%s)

for ((i=0; i<TOTAL; i++)); do
    IFS=':' read -r name mode r1 r2 trim_flag <<< "${SAMPLES[$i]}"
    idx=$((i+1))
    log ""
    log "=========================================="
    log "[$idx/$TOTAL] $name ($mode)"
    log "=========================================="

    SAMPLE_DIR="${OUT_DIR}/${name}"
    BAM_DIR="${SAMPLE_DIR}/bam"
    TRIM_DIR="${SAMPLE_DIR}/trimmed"
    EXTRACT_DIR="${SAMPLE_DIR}/meth_extract"
    mkdir -p "$SAMPLE_DIR" "$BAM_DIR" "$TRIM_DIR" "$EXTRACT_DIR"

    t0=$(date +%s)

    # --- Trimming ---
    if [ "$trim_flag" = "trim" ] && [ "$mode" = "se" ]; then
        r1=$(run_fastp_se "$r1" "$TRIM_DIR" "$name")
    fi

    # --- Bismark ---
    EXISTING_BAM=$(safe_find "$BAM_DIR" "*_bismark_bt2*.bam")
    if [ -n "$EXISTING_BAM" ]; then
        BAM_FILE="$EXISTING_BAM"
        log "   Bismark已完成: $(basename $BAM_FILE)"
    else
        log "   Bismark $mode: $name"
        if [ "$mode" = "pe" ]; then
            bismark --bowtie2 --genome_folder "$REF_DIR" \
                --score_min L,0,-0.2 --parallel $BISMARK_PARALLEL -p $BISMARK_BT2_THREADS \
                --bam --quiet -o "$BAM_DIR" \
                -1 "$r1" -2 "$r2" -X 700 --dovetail 2>&1 | tail -3 || true
        else
            bismark --bowtie2 --genome_folder "$REF_DIR" \
                --score_min L,0,-0.2 --parallel $BISMARK_PARALLEL -p $BISMARK_BT2_THREADS \
                --bam --quiet -o "$BAM_DIR" "$r1" 2>&1 | tail -3 || true
        fi
        BAM_FILE=$(safe_find "$BAM_DIR" "*_bismark_bt2*.bam")
        if [ -z "$BAM_FILE" ]; then log "   ✗ Bismark失败!"; FAILED=$((FAILED+1)); continue; fi
    fi

    # --- Sort BAM (Bismark --parallel 不自动排序) ---
    SORTED_BAM="${BAM_DIR}/$(basename "${BAM_FILE%.bam}").sorted.bam"
    if [ -f "$SORTED_BAM" ]; then
        log "   Sort已完成: $(basename $SORTED_BAM)"
        BAM_FILE="$SORTED_BAM"
    elif [ "$(samtools view -H "$BAM_FILE" 2>/dev/null | grep SO: | head -1)" = "@HD	VN:1.0	SO:coordinate" ]; then
        log "   Sort已排序 (SO:coordinate)"
    else
        log "   Sort: $name"
        samtools sort -@ 4 -o "$SORTED_BAM" "$BAM_FILE" 2>&1 | tail -3 || true
        BAM_FILE="$SORTED_BAM"
    fi

    # --- Dedup ---
    DEDUP_BAM=$(safe_find "$BAM_DIR" "*.deduplicated.bam")
    if [ -n "$DEDUP_BAM" ]; then
        log "   Dedup已完成: $(basename $DEDUP_BAM)"
    else
        log "   Dedup $mode: $name"
        if [ "$mode" = "pe" ]; then
            deduplicate_bismark --paired --output_dir "$BAM_DIR" --bam "$BAM_FILE" 2>&1 | tail -3 || true
        else
            deduplicate_bismark --single --output_dir "$BAM_DIR" --bam "$BAM_FILE" 2>&1 | tail -3 || true
        fi
        DEDUP_BAM=$(safe_find "$BAM_DIR" "*.deduplicated.bam")
        if [ -z "$DEDUP_BAM" ]; then log "   ✗ Dedup失败!"; FAILED=$((FAILED+1)); continue; fi
    fi

    # --- Methylation Extraction ---
    N_EXT=$(safe_find_count "$EXTRACT_DIR")
    if [ "$N_EXT" -ge 3 ]; then
        log "   Extract已完成 ($N_EXT files)"
    else
        log "   Extract $mode: $name"
        if [ "$mode" = "pe" ]; then
            bismark_methylation_extractor --multicore $METH_EXTRACT_CORES \
                --comprehensive --gzip --paired-end --no_overlap \
                --ignore 5 --ignore_r2 5 --output_dir "$EXTRACT_DIR" "$DEDUP_BAM" 2>&1 | tail -3 || true
        else
            bismark_methylation_extractor --multicore $METH_EXTRACT_CORES \
                --comprehensive --gzip --single-end \
                --ignore 5 --output_dir "$EXTRACT_DIR" "$DEDUP_BAM" 2>&1 | tail -3 || true
        fi
    fi

    # --- Methylation Calling ---
    MC_FILE="${SAMPLE_DIR}/${name}.mC_level_Identification_stat.txt"
    if [ -f "$MC_FILE" ]; then
        N_LINES=$(wc -l < "$MC_FILE" 2>/dev/null || echo 0)
        log "   Caller已完成 ($N_LINES lines)"
    else
        log "   Caller: $name"
        python -c "
from wgbs_tools.modules.caller import run_methylation_calling
run_methylation_calling('$EXTRACT_DIR', '$SAMPLE_DIR', '$name', '$FASTA', cores=$CALLER_CORES)
" 2>&1 | tail -5 || { log "   ✗ Caller失败!"; FAILED=$((FAILED+1)); continue; }
    fi

    # --- Summary ---
    SUM_FILE="${SAMPLE_DIR}/${name}.methylation_summary.stat.txt"
    if [ -f "$SUM_FILE" ]; then
        log "   Summary已完成"
    else
        log "   Summary: $name"
        python -c "
from wgbs_tools.modules.summary import generate_summary
from wgbs_tools.modules.genome import compute_and_cache_c_stats
c_stats = compute_and_cache_c_stats('$REF_DIR', 'GCF_000182925.2_NC12_genomic.fna')
generate_summary('$MC_FILE', c_stats, '$SAMPLE_DIR', '$name')
" 2>&1 | tail -3 || { log "   ✗ Summary失败!"; FAILED=$((FAILED+1)); continue; }
    fi

    elapsed=$(($(date +%s) - t0))
    log "   ✓ 完成! (${elapsed}s)"
    PASSED=$((PASSED+1))
done

TOTAL_TIME=$(($(date +%s) - START_TIME))
log ""
log "=========================================="
log "全部完成! 成功: $PASSED/$TOTAL"
log "总耗时: $((TOTAL_TIME/60))m $((TOTAL_TIME%60))s"
log "=========================================="
