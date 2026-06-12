#!/bin/bash
# WGBS_tools v2: sorted BAM + parallel caller (up to 4 concurrent)
# 用法: nohup bash run_all_v2.sh > run_all_v2.log 2>&1 &
set -euo pipefail

# 直接使用完整路径，避免 conda activate hang
export PATH="/data/home/zhength/miniconda3/envs/wgbs_final/bin:$PATH"
export PYTHON="/data/home/zhength/miniconda3/envs/wgbs_final/bin/python"
BISMARK_CMD="/data/home/zhength/miniconda3/envs/wgbs_final/bin/bismark"
SAMTOOLS_CMD="/data/home/zhength/miniconda3/envs/wgbs_final/bin/samtools"
DEDUP_CMD="/data/home/zhength/miniconda3/envs/wgbs_final/bin/deduplicate_bismark"
EXTRACT_CMD="/data/home/zhength/miniconda3/envs/wgbs_final/bin/bismark_methylation_extractor"
FASTP_CMD="/data/home/zhength/miniconda3/envs/wgbs_final/bin/fastp"

# ========== CONFIG ==========
REF_DIR="/data/home/zhength/work/2026.6.11_CaMethy_suppleFig/ref"
FASTA="${REF_DIR}/GCF_000182925.2_NC12_genomic.fna"
OUT_DIR="/data/home/zhength/work/2026.6.11_CaMethy_suppleFig/WGBS_tools_new_result_v2"
RAWDATA_USER="/data/home/zhength/work/2026.6.11_CaMethy_suppleFig/rawdata"
RAWDATA_NCBI="/data/home/zhength/work/2026.6.11_CaMethy_suppleFig/rawdata_NCBI"
OLD_OUT_DIR="/data/home/zhength/work/2026.6.11_CaMethy_suppleFig/WGBS_tools_new_result"

# 资源: ~72核总量
BISMARK_PARALLEL=5; BISMARK_BT2_THREADS=7   # Bismark: 5*8=40核
SORT_THREADS=4
METH_EXTRACT_CORES=8
CALLER_CORES=20                              # caller: 最多4样本并行 → 4*20=80, 但实际上受3文件限制
FASTP_THREADS=4
MAX_PARALLEL_CALLERS=4

# ========== HELPERS ==========
log() { echo "[$(date '+%m-%d %H:%M:%S')] $*"; }
safe_find() { find "$1" -maxdepth 1 -name "$2" 2>/dev/null | head -1 || true; }
safe_count() { find "$1" -maxdepth 1 -name "*.txt.gz" 2>/dev/null | wc -l || echo 0; }

# ========== FUNCTIONS ==========
run_fastp_se() {
    local r1=$1 d=$2 n=$3; mkdir -p "$d"
    local o="${d}/${n}_trimmed_1.fq.gz"
    [ -f "$o" ] && { echo "$o"; return; }
    log "    Fastp: $n" >&2
    $FASTP_CMD --in1 "$r1" --out1 "$o" \
        --cut_front_window_size=1 --cut_front_mean_quality=3 \
        --cut_tail_window_size=1 --cut_tail_mean_quality=3 \
        --cut_right_window_size=4 --cut_right_mean_quality=15 \
        --length_required=36 --thread $FASTP_THREADS \
        --html "${d}/${n}_fastp.html" --json "${d}/${n}_fastp.json" >/dev/null 2>&1
    echo "$o"
}
run_sort_bam() {
    local bam=$1 d=$2 n=$3
    local sorted="${d}/$(basename "${bam%.bam}").sorted.bam"
    [ -f "$sorted" ] && { echo "$sorted"; return; }
    local so_tag
    so_tag=$(samtools view -H "$bam" 2>/dev/null | grep "^@HD" | grep -o "SO:coordinate" || true)
    [ -n "$so_tag" ] && { echo "$bam"; return; }
    log "    Sort: $n"
    $SAMTOOLS_CMD sort -@ $SORT_THREADS -o "$sorted" "$bam" >/dev/null 2>&1
    echo "$sorted"
}
run_dedup_v2() {
    local bam=$1 d=$2 n=$3 pe=$4
    local mode; mode=$([ "$pe" = "pe" ] && echo "PE" || echo "SE")
    local flag; flag=$([ "$pe" = "pe" ] && echo "--paired" || echo "--single")
    log "    Dedup $mode: $n"
    $DEDUP_CMD $flag --output_dir "$d" --bam "$bam" >/dev/null 2>&1
    local out; out=$(safe_find "$d" "*.deduplicated.bam")
    [ -z "$out" ] && { log "    ✗ Dedup 失败!"; return 1; }
    echo "$out"
}

run_bismark_pe() {
    local r1=$1 r2=$2 d=$3 n=$4
    log "    Bismark PE: $n"
    $BISMARK_CMD --bowtie2 --genome_folder "$REF_DIR" \
        --score_min L,0,-0.2 --parallel $BISMARK_PARALLEL -p $BISMARK_BT2_THREADS \
        --bam --quiet -o "$d" -1 "$r1" -2 "$r2" -X 700 --dovetail 2>&1 | tail -2 || true
}

run_bismark_se() {
    local r1=$1 d=$2 n=$3
    log "    Bismark SE: $n"
    $BISMARK_CMD --bowtie2 --genome_folder "$REF_DIR" \
        --score_min L,0,-0.2 --parallel $BISMARK_PARALLEL -p $BISMARK_BT2_THREADS \
        --bam --quiet -o "$d" "$r1" 2>&1 | tail -2 || true
}

run_extract() {
    local bam=$1 d=$2 n=$3 pe=$4
    mkdir -p "$d"
    local nf; nf=$(safe_count "$d")
    [ "$nf" -ge 3 ] && { log "    Extract 已完成 ($nf files)"; return 0; }
    local args mode; mode=$([ "$pe" = "pe" ] && echo "PE" || echo "SE")
    if [ "$pe" = "pe" ]; then
        args="--paired-end --no_overlap --ignore 5 --ignore_r2 5"
    else
        args="--single-end --ignore 5"
    fi
    log "    Extract $mode: $n"
    $EXTRACT_CMD --multicore $METH_EXTRACT_CORES \
        --comprehensive --gzip $args --output_dir "$d" "$bam" 2>&1 | tail -2 || true
}

run_caller() {
    local ext_dir=$1 out_dir=$2 n=$3
    local mc="${out_dir}/${n}.mC_level_Identification_stat.txt"
    [ -f "$mc" ] && { local l; l=$(wc -l < "$mc"); log "    Caller 已完成 ($l lines)"; return 0; }
    log "    Caller: $n"
    $PYTHON -c "
from wgbs_tools.modules.caller import run_methylation_calling
run_methylation_calling('$ext_dir', '$out_dir', '$n', '$FASTA', cores=$CALLER_CORES, min_mc=2, min_depth=5)
" 2>&1 | tail -3
    [ -f "$mc" ] || { log "    ✗ Caller 失败!"; return 1; }
}

run_summary() {
    local mc=$1 out_dir=$2 n=$3
    local sf="${out_dir}/${n}.methylation_summary.stat.txt"
    [ -f "$sf" ] && { log "    Summary 已完成"; return 0; }
    log "    Summary: $n"
    $PYTHON -c "
from wgbs_tools.modules.summary import generate_summary
from wgbs_tools.modules.genome import compute_and_cache_c_stats
c = compute_and_cache_c_stats('$REF_DIR', 'GCF_000182925.2_NC12_genomic.fna')
generate_summary('$mc', c, '$out_dir', '$n')
" 2>&1 | tail -2
}

# ========== SAMPLE LIST ==========
SAMPLES=(
    "NC2M_1:pe:${RAWDATA_USER}/NC2M_1_clean_1.fq.gz:${RAWDATA_USER}/NC2M_1_clean_2.fq.gz:skip"
    "NC2M_2:pe:${RAWDATA_USER}/NC2M_2_clean_1.fq.gz:${RAWDATA_USER}/NC2M_2_clean_2.fq.gz:skip"
    "NC2M_3:pe:${RAWDATA_USER}/NC2M_3_clean_1.fq.gz:${RAWDATA_USER}/NC2M_3_clean_2.fq.gz:skip"
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
START_TIME=$(date +%s)
PASSED=0; FAILED=0

# ========== PHASE 1: 基因组准备 ==========
log "=========================================="
log "WGBS_tools v2 (sorted BAM + parallel caller)"
log "=========================================="
python -c "
from wgbs_tools.modules.genome import compute_and_cache_c_stats
s = compute_and_cache_c_stats('$REF_DIR', 'GCF_000182925.2_NC12_genomic.fna')
print(f'C: total={s[\"totals\"][\"total_C\"]:,}')
"

# ========== PHASE 2: Bismark→Sort→Dedup→Extract (顺序, 每个样本) ==========
log ""
log "=========================================="
log "PHASE 2: 比对+排序+去重+提取 (顺序执行)"
log "=========================================="

declare -A SAMPLE_INFO  # store extract_dir for later caller

for ((i=0; i<TOTAL; i++)); do
    IFS=':' read -r name mode r1 r2 trim_flag <<< "${SAMPLES[$i]}"
    idx=$((i+1))
    log ""
    log "--- [$idx/$TOTAL] $name ($mode) ---"

    SD="${OUT_DIR}/${name}"; mkdir -p "$SD"
    BD="${SD}/bam"; mkdir -p "$BD"
    TD="${SD}/trimmed"; mkdir -p "$TD"
    ED="${SD}/meth_extract"; mkdir -p "$ED"
    SAMPLE_INFO["${name}_extract"]="$ED"

    t0=$(date +%s)

    # Trimming
    [ "$trim_flag" = "trim" ] && [ "$mode" = "se" ] && r1=$(run_fastp_se "$r1" "$TD" "$name")

    # Bismark alignment (reuse from old output dir if available)
    # 排除 temp.* 文件 (残留的不完整分片)
    EXISTING_RAW=$(find "$BD" -maxdepth 1 -name "*_bismark_bt2*.bam" ! -name "*.temp.*" 2>/dev/null | head -1 || true)
    OLD_RAW=$(find "${OLD_OUT_DIR}/${name}/bam" -maxdepth 1 -name "*_bismark_bt2_pe.bam" ! -name "*.temp.*" 2>/dev/null | head -1 || true)
    [ -z "$OLD_RAW" ] && OLD_RAW=$(find "${OLD_OUT_DIR}/${name}/bam" -maxdepth 1 -name "*_bismark_bt2.bam" ! -name "*.temp.*" 2>/dev/null | head -1 || true)
    if [ -n "$EXISTING_RAW" ]; then
        BAM_FILE="$EXISTING_RAW"
        log "    Bismark 已有: $(basename $BAM_FILE)"
    elif [ -n "$OLD_RAW" ]; then
        cp "$OLD_RAW" "$BD/" 2>/dev/null || ln "$OLD_RAW" "$BD/" 2>/dev/null
        BAM_FILE=$(safe_find "$BD" "*_bismark_bt2*.bam")
        log "    Bismark 复用旧 BAM: $(basename $BAM_FILE)"
    else
        [ "$mode" = "pe" ] && run_bismark_pe "$r1" "$r2" "$BD" "$name" || run_bismark_se "$r1" "$BD" "$name"
        BAM_FILE=$(find "$BD" -maxdepth 1 -name "*_bismark_bt2*.bam" ! -name "*.temp.*" 2>/dev/null | head -1 || true)
        [ -z "$BAM_FILE" ] && { log "    ✗ Bismark失败!"; FAILED=$((FAILED+1)); continue; }
    fi

    # Sort
    BAM_FILE=$(run_sort_bam "$BAM_FILE" "$BD" "$name")

    # Dedup
    DEDUP_BAM=$(safe_find "$BD" "*.sorted.deduplicated.bam")
    [ -z "$DEDUP_BAM" ] && DEDUP_BAM=$(safe_find "$BD" "*.deduplicated.bam")
    if [ -n "$DEDUP_BAM" ] && [ "$(samtools view -H "$DEDUP_BAM" 2>/dev/null | grep -c SO:coordinate)" -gt 0 ]; then
        log "    Dedup 已完成 (sorted): $(basename $DEDUP_BAM)"
    else
        rm -f "$BD"/*.deduplicated* 2>/dev/null  # remove old unsorted dedup
        DEDUP_BAM=$(run_dedup_v2 "$BAM_FILE" "$BD" "$name" "$mode") || { FAILED=$((FAILED+1)); continue; }
    fi

    # Methylation extraction (reuse from old output if available)
    OLD_EXT_DIR="${OLD_OUT_DIR}/${name}/meth_extract"
    if [ "$(safe_count "$ED")" -ge 3 ]; then
        log "    Extract 已完成 ($(safe_count "$ED") files)"
    elif [ -d "$OLD_EXT_DIR" ] && [ "$(safe_count "$OLD_EXT_DIR")" -ge 3 ]; then
        cp -l "$OLD_EXT_DIR"/*.txt.gz "$ED/" 2>/dev/null || cp "$OLD_EXT_DIR"/*.txt.gz "$ED/" 2>/dev/null
        log "    Extract 复用旧数据 ($(safe_count "$ED") files)"
    else
        rm -rf "$ED"/* 2>/dev/null
        mkdir -p "$ED"
        run_extract "$DEDUP_BAM" "$ED" "$name" "$mode"
    fi

    elapsed=$(($(date +%s) - t0))
    log "  ✓ $name phase2 完成 (${elapsed}s)"
    PASSED=$((PASSED+1))
done

BISMARK_TIME=$(($(date +%s) - START_TIME))
log ""
log "PHASE 2 完成! 耗时: $((BISMARK_TIME/60))m"

# ========== PHASE 3: Caller + Summary (并行, 最多4个) ==========
log ""
log "=========================================="
log "PHASE 3: 甲基化调用+摘要 (并行, max $MAX_PARALLEL_CALLERS)"
log "=========================================="

# 用临时文件追踪每个样本的状态
TMP_DIR="/tmp/wgbs_caller_$$"
mkdir -p "$TMP_DIR"

run_caller_and_summary() {
    local name=$1 extract_dir=$2 sample_dir=$3
    (
        if run_caller "$extract_dir" "$sample_dir" "$name"; then
            local mc="${sample_dir}/${name}.mC_level_Identification_stat.txt"
            run_summary "$mc" "$sample_dir" "$name"
            echo "DONE" > "${TMP_DIR}/${name}"
            log "  ✓ $name 完成!"
        else
            echo "FAILED" > "${TMP_DIR}/${name}"
            log "  ✗ $name 失败!"
        fi
    )
}

# 收集需要 run caller 的样本
CALLER_JOBS=()
for ((i=0; i<TOTAL; i++)); do
    IFS=':' read -r name mode r1 r2 trim_flag <<< "${SAMPLES[$i]}"
    ED="${OUT_DIR}/${name}/meth_extract"
    SD="${OUT_DIR}/${name}"
    MC="${SD}/${name}.mC_level_Identification_stat.txt"
    if [ -f "$MC" ]; then
        log "  $name: Caller 已跳过 (文件存在)"
    else
        CALLER_JOBS+=("$name:$ED:$SD")
    fi
done

# 并行执行 caller, 每次最多 MAX_PARALLEL_CALLERS 个
ACTIVE=0; JOB_COUNT=${#CALLER_JOBS[@]}
log "待处理 caller jobs: $JOB_COUNT"

for job in "${CALLER_JOBS[@]}"; do
    IFS=':' read -r jname jed jsd <<< "$job"

    # Wait if we've reached max parallel
    while [ "$(jobs -r | wc -l)" -ge $MAX_PARALLEL_CALLERS ]; do
        sleep 10
    done

    log "  启动 caller: $jname"
    run_caller_and_summary "$jname" "$jed" "$jsd" &
    ACTIVE=$((ACTIVE+1))
done

# Wait for all
log "等待所有 caller 完成..."
wait

# ========== SUMMARY ==========
TOTAL_TIME=$(($(date +%s) - START_TIME))
FINAL_PASSED=0; FINAL_FAILED=0
for ((i=0; i<TOTAL; i++)); do
    IFS=':' read -r name mode r1 r2 trim_flag <<< "${SAMPLES[$i]}"
    MC="${OUT_DIR}/${name}/${name}.mC_level_Identification_stat.txt"
    SF="${OUT_DIR}/${name}/${name}.methylation_summary.stat.txt"
    if [ -f "$MC" ] && [ -f "$SF" ]; then
        N=$(wc -l < "$MC")
        log "  ✓ $name: $N sites"
        FINAL_PASSED=$((FINAL_PASSED+1))
    else
        log "  ✗ $name: 缺少输出文件"
        FINAL_FAILED=$((FINAL_FAILED+1))
    fi
done

rm -rf "$TMP_DIR"
log ""
log "=========================================="
log "全部完成! 成功: $FINAL_PASSED/$TOTAL, 失败: $FINAL_FAILED"
log "总耗时: $((TOTAL_TIME/3600))h $(((TOTAL_TIME%3600)/60))m"
log "输出: $OUT_DIR"
log "=========================================="
