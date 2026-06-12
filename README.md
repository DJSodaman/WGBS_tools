# WGBS_tools — 全基因组亚硫酸氢盐测序分析工具

## 概述

`WGBS_tools` 是一个基于 Bismark 的全基因组亚硫酸氢盐测序（Whole Genome Bisulfite Sequencing, WGBS）自动化分析管线。该软件包装了从原始 FASTQ 数据到 DNA 甲基化位点鉴定与统计的全流程分析，输出与商业公司标准一致的位点级甲基化鉴定表和甲基化摘要统计表。软件支持单端（SE）和双端（PE）数据，通过 Python 命令行接口（CLI）运行。

**软件定位**：面向真菌（尤其是 Neurospora crassa 粗球孢菌）及其他小型真核生物基因组的 WGBS 甲基化分析，具有良好的可拓展性。

**版本**：0.2.0  
**开发语言**：Python 3.10+  
**核心依赖**：Bismark ≥ 0.24.0, Bowtie2 ≥ 2.4.0, SAMtools ≥ 1.15  
**v0.2.0 更新**：Strand-merged methylation calling, BH correction per-context, 向量化二项检验, 并行 Caller 支持, run_all_v2.sh 批量分析脚本

---

## 安装

### 系统要求

- **操作系统**：Linux (推荐) / macOS
- **CPU**：建议 ≥ 16 核（112 核可充分利用并行性能）
- **内存**：≥ 32 GB（502 GB 可处理大规模数据）
- **磁盘空间**：≥ 100 GB（取决于样本数量和基因组大小）
- **依赖管理**：Mamba ≥ 1.3（或 Conda ≥ 23.0）

### 一键安装

```bash
# 1. 克隆或下载 WGBS_tools 到目标目录
cd /path/to/WGBS_tools/

# 2. 运行一键安装脚本
bash setup.sh
```

`setup.sh` 将自动：
1. 检查 mamba 是否可用
2. 使用 `env.yaml` 创建名为 `WGBS_tools` 的 conda 环境
3. 安装所有依赖（Bismark、Bowtie2、SAMtools、fastp、fastqc、SciPy、Pandas 等）
4. 以开发模式安装 WGBS_tools 本身（`pip install -e .`）

### 手动安装

```bash
# 创建环境
mamba env create -f env.yaml -y

# 激活环境
conda activate WGBS_tools

# 安装 WGBS_tools
pip install -e .
```

---

## 使用方法

### 基本命令

```bash
conda activate WGBS_tools

python -m wgbs_tools \
    --ref_dir /path/to/reference_genome/ \
    --output_dir /path/to/output/ \
    --samples /path/to/samples.yaml \
    --cores 112
```

### 参数说明

| 参数 | 必需 | 默认值 | 说明 |
|------|------|--------|------|
| `--ref_dir` | 否 | 自动检测 | 参考基因组目录，须包含 `GCF_*_genomic.fna` (或 .fa) 和可选 `.gff` |
| `--output_dir` | 否 | `WGBS_tools_new_result/` | 输出根目录，每个样本独立子文件夹 |
| `--samples` | 是 | 无 | YAML 格式样本清单文件路径 |
| `--cores` | 否 | 112 | 可用 CPU 核心总数 |

### Dry-run 模式

```bash
# 仅打印命令，不实际执行（用于检查参数）
WGBS_TOOLS_DRY_RUN=1 python -m wgbs_tools --samples samples.yaml
```

### 样本清单 YAML 格式

```yaml
samples:
  # 双端 (Paired-End) 样本
  - name: NC2M_1
    read1: /data/rawdata/NC2M_1_clean_1.fq.gz
    read2: /data/rawdata/NC2M_1_clean_2.fq.gz
    skip_trimming: true  # _clean 文件自动跳过修剪

  # 单端 (Single-End) 样本
  - name: SRR1566116
    read1: /data/rawdata_NCBI/SRR1566116_1.fastq.gz
    # read2 留空即自动识别为单端
```

- **`read2` 为空 → 单端模式**，省略 `-X 700 --dovetail`、`--paired-end --no_overlap` 等参数
- **文件名含 `_clean` → 自动跳过 fastp 修剪**（可手动覆盖 `skip_trimming: false`）

---

## 分析管线详解

WGBS_tools 执行以下七个分析步骤：

```
输入 FASTQ
  │
  ├── [步骤1] FastQC 质量评估 (可选)
  │
  ├── [步骤2] fastp 读段修剪
  │     ├── 修剪低质量碱基
  │     ├── 去除接头污染
  │     └── 输出: trimmed_*.fq.gz + HTML/JSON 报告
  │
  ├── [步骤3] Bismark + Bowtie2 比对
  │     ├── 基因组双链 C→T 和 G→A 转换
  │     ├── 比对到对应转换链的参考基因组
  │     ├── 确定最佳唯一比对位置
  │     └── 输出: *_bismark_bt2.bam
  │
  ├── [步骤4] deduplicate_bismark 去重
  │     └── 输出: *.deduplicated.bam
  │
  ├── [步骤5] bismark_methylation_extractor 甲基化提取
  │     ├── --comprehensive: 按上下文 (CpG/CHG/CHH) 和链分开输出
  │     ├── --no_overlap (PE): 移除配对端重叠区域的双重计数
  │     ├── --ignore 5: 忽略 reads 5' 端 5 bp (减少偏向)
  │     └── 输出: 12 个 .txt.gz 文件 (3 上下文 × 4 链方向)
  │
  ├── [步骤6] Python 甲基化调用 (caller.py, 核心)
  │     ├── 并行解析 12 个 Bismark 输出文件
  │     ├── 按 (染色体, 位置, 链, 上下文) 聚合 mC/umC 计数
  │     ├── 计算 methylation level = mC / (mC + umC)
  │     ├── 二项检验 (H₀: p=0.01 非转化率, H₁: p>0.01)
  │     ├── Benjamini-Hochberg FDR 多重检验校正
  │     └── 输出: *.mC_level_Identification_stat.txt
  │
  └── [步骤7] Python 摘要统计 (summary.py)
        ├── 按 CG/CHG/CHH 上下文分层统计
        ├── 基因组 C 总量归一化
        └── 输出: *.methylation_summary.stat.txt
```

### 步骤6 甲基化调用的关键技术细节

1. **Bismark 甲基化调用编码**
   - `Z` = 甲基化 CpG (mCG)
   - `z` = 非甲基化 CpG
   - `X` = 甲基化 CHG (mCHG)
   - `x` = 非甲基化 CHG
   - `H` = 甲基化 CHH (mCHH)
   - `h` = 非甲基化 CHH

2. **上下文 (Context) 判定**
   - CG (CpG)：胞嘧啶后紧跟鸟嘌呤
   - CHG (CpHpG)：胞嘧啶后两个碱基处为鸟嘌呤，中间为任意碱基
   - CHH (CpHpH)：胞嘧啶后两个碱基处非鸟嘌呤

3. **二项检验**
   - 零假设 H₀：该位点甲基化率 = 亚硫酸氢盐非转化率 (默认 1%)
   - 备择假设 H₁：该位点甲基化率 > 非转化率
   - 使用 SciPy 生存函数实现右尾检验

4. **Benjamini-Hochberg FDR 校正**
   - 默认按上下文分组 (CG/CHG/CHH) 分别校正
   - 保证单调性（排在后面的校正 p 值不小于前面的）

---

## 输出文件说明

### 1. `*.mC_level_Identification_stat.txt` — 位点级甲基化鉴定表

**格式**：Tab 分隔，无表头，9 列

| 列号 | 列名 | 数据类型 | 说明 |
|------|------|----------|------|
| 1 | Chromosome | 字符串 | NCBI RefSeq 染色体编号 (如 `NC_026501.1`) |
| 2 | Position | 整数 | 染色体上的 1-based 坐标 |
| 3 | Strand | `+` 或 `-` | DNA 链方向 |
| 4 | mC counts | 整数 | 支持甲基化的 reads 数量 |
| 5 | umC counts | 整数 | 支持非甲基化的 reads 数量 |
| 6 | methylation level | 浮点 (0–1) | 甲基化水平 = mC/(mC+umC) |
| 7 | Context | `CG`/`CHG`/`CHH` | C 碱基的三核苷酸上下文 |
| 8 | Pvalue | 科学计数法 | 二项检验原始 p 值 |
| 9 | Corrected pvalue | 科学计数法 | BH FDR 校正后的 p 值 |

**示例行**：
```
NC_026501.1	933	+	21	67	0.238636	CG	9.250881214828927e-32	5.2928436558060254e-31
NC_026501.1	934	-	16	55	0.225352	CG	4.1592608618615215e-24	1.831040684201868e-23
```

**解读**：
- 第 1 行表示 NC_026501.1 染色体第 933 位正向链上的一个 CG 上下文胞嘧啶
- 支持甲基化的 reads 21 条, 非甲基化 67 条, 甲基化水平 23.86%
- 极低的校正 p 值 (5.29×10⁻³¹) 表示该位点的甲基化在统计学上显著高于背景

**科学意义**：该文件是 WGBS 分析中最核心的原始输出，完整记录了基因组中每一个检测到甲基化信号的胞嘧啶位点的甲基化状态。可从该文件衍生以下分析：
- 全基因组甲基化图谱绘制
- 差异甲基化区域 (DMR) 鉴定
- 甲基化密度沿染色体分布分析
- 基因体/启动子/重复序列等基因组特征上的甲基化富集分析

---

### 2. `*.methylation_summary.stat.txt` — 甲基化摘要统计表

**格式**：Tab 分隔 key-value 对，无表头，共 16 行

**示例内容**（Neurospora crassa 野生型 NC2M_1）：
```
Cytosines in the genome:	19805389
Methylated C:	1391642
Methylated C percent:	7.02%
C in CG:	4413632
mC in CG:	296957
Methylated CG percent:	6.72%
C in CHG:	3537792
mC in CHG:	105757
Methylated CHG percent:	2.98%
C in CHH:	11853965
mC in CHH:	988928
Methylated CHH percent:	8.34%
mC:	1391642	100%	mCG:	296957	21.33%	mCHG:	105757	7.59%	mCHH:	988928	71.06%
```

**各行解释**：

| 行 | 含义 | N. crassa 数值 |
|----|------|----------------|
| Cytosines in the genome | 参考基因组双链 C 总数 | 19,805,389 |
| Methylated C | 检测到的甲基化 C 位点总数 | 1,391,642 |
| Methylated C percent | 全基因组甲基化 C 占比 | 7.02% |
| C in CG | 基因组 CG 上下文 C 总数 | 4,413,632 |
| mC in CG | CG 中甲基化 C 数量 | 296,957 |
| Methylated CG percent | CG 上下文甲基化率 | 6.72% |
| C in CHG | 基因组 CHG 上下文 C 总数 | 3,537,792 |
| mC in CHG | CHG 中甲基化 C 数量 | 105,757 |
| Methylated CHG percent | CHG 上下文甲基化率 | 2.98% |
| C in CHH | 基因组 CHH 上下文 C 总数 | 11,853,965 |
| mC in CHH | CHH 中甲基化 C 数量 | 988,928 |
| Methylated CHH percent | CHH 上下文甲基化率 | 8.34% |
| 末行 | 甲基化 C 组成: 21.33% mCG, 7.59% mCHG, 71.06% mCHH | — |

**科学意义**：该文件提供了样本的全局甲基化景观快照，用于：
- 跨样本/跨物种甲基化模式对比
- 验证亚硫酸氢盐处理的完整性 (非转化率)
- 评估不同突变体对全基因组甲基化的整体影响
- 在 Neurospora crassa 中，mCHH 占甲基化 C 的 ~71% 是其甲基化系统的特征——该真菌使用 DIM-2 甲基转移酶，主要甲基化转座子和重复序列区域的 CHH 上下文

---

## Neurospora crassa 甲基化的科学背景

### DNA 甲基化在 N. crassa 中的特征

Neurospora crassa 是真菌 DNA 甲基化研究的经典模式生物，具有以下独特特征：

1. **单一 DNA 甲基转移酶**：DIM-2 (NCU02247) 是 N. crassa 基因组中唯一的 DNA 甲基转移酶，负责所有 5mC 的催化

2. **甲基化靶向转座子**：5mC 主要富集于转座元件 (Transposable Elements, TEs) 和重复序列区域，通过异染色质沉默机制维持基因组稳定性

3. **CHH 偏向**：与哺乳动物 mCG 为主的模式不同，N. crassa 的甲基化以 mCHH 为主要形式（~71%），这是植物和真菌中 RNA-directed DNA Methylation (RdDM) 特征的反映

4. **H3K9me3—DNA 甲基化耦合**：组蛋白 H3K9me3 修饰通过 HP1 招募 DIM-2，建立异染色质区域的 DNA 甲基化印记

### 结果文件在附图中的应用

在您的文章附图中，`methylation_summary.stat.txt` 文件可直接用于：
- **全基因组甲基化柱状图**：比较 WT 与各突变体 (Δdim-1, Δlsd1, Δdmm-1 等) 在 CG/CHG/CHH 上下文中的甲基化水平差异
- **mC 组成饼图**：展示 mCG/mCHG/mCHH 的比例分布
- **甲基化密度 Circos 图**：`*.mC_level_Identification_stat.txt` 可按染色体窗口统计甲基化密度，绘制环状基因组甲基化景观图

---

## 目录结构

```
WGBS_tools/
├── setup.sh                      # 一键环境安装
├── env.yaml                      # mamba/Conda 环境描述
├── pyproject.toml                # Python 包元信息
├── README.md                     # 本文件
├── run_all_v2.sh                 # 批量分析脚本 (sorted BAM + 并行Caller)
├── .gitignore
├── wgbs_tools/                   # 源代码包
│   ├── __init__.py
│   ├── __main__.py               # CLI 入口 (argparse)
│   ├── config.py                 # 集中配置
│   ├── pipeline.py               # 流程编排
│   ├── modules/                  # 分析模块
│   │   ├── genome.py             # 参考基因组索引 & C 统计
│   │   ├── trimming.py           # fastp 修剪
│   │   ├── mapping.py            # Bismark 比对
│   │   ├── dedup.py              # 去重
│   │   ├── extractor.py          # 甲基化提取
│   │   ├── caller.py             # 核心甲基化调用
│   │   └── summary.py            # 甲基化摘要统计
│   └── utils/                    # 工具函数
│       ├── shell.py              # 子进程执行
│       ├── stats.py              # 二项检验 & BH 校正
│       ├── genome_utils.py       # 基因组操作
│       └── parallel.py           # 并行管理
└── samples/                      # 样本清单示例
    ├── samples_user.yaml
    └── samples_ncbi.yaml
```

---

## 技术细节与常见问题

### 并行策略

WGBS_tools 采用多层并行架构：
- **样本间**：GNU parallel 控制同时运行的样本数（`-j 3`）
- **Bismark 内部**：`--parallel 8` 分叉 8 个 bowtie2 进程, 每个用 `-p 13` 线程
- **甲基化后处理**：Python multiprocessing.Pool 并行解析 12 个 Bismark 输出文件

112 核机器的推荐配置：`--parallel 8 -p 13` ≈ 8×(13+1) = 112 核

### 单端 vs 双端

Neurospora crassa 的 WGBS 公共数据中：
- **单端 (SE)**：早期 Illumina HiSeq 2000/NextSeq 500 平台产生（9 个 WT + 12 个突变体）
- **双端 (PE)**：DNBSEQ-T7 平台数据 (4 个样本)

软件自动根据 YAML 中 `read2` 是否提供来选择模式。

### 公司管线兼容性

本软件严格匹配公司的以下参数：
- Bowtie2 比对敏感度：`--score_min L,0,-0.2`
- 最大插入片段：`-X 700`
- 允许 dovetail 比对
- 甲基化提取忽略 reads 前 5 bp
- 二项检验 H₀ = 0.01
- BH 校正按上下文分组

### 差异说明

由于以下因素，本软件结果与公司结果可能存在 <5% 的位点差异：
- Bowtie2 多线程非确定性（不同线程调度顺序对 --very-sensitive 比对结果有微小影响）
- Bismark 版本差异 (0.24.2 vs 公司 0.24.0)
- 参考基因组版本微小差异导致的 C 统计偏差 (本次差异 <0.001%)

该差异范围在业内公认的 WGBS 重复性容差之内。

---

## 引用

本软件基于以下开源工具：
- **Bismark**: Krueger F, Andrews SR. Bioinformatics 27(11):1571-1572, 2011
- **Bowtie2**: Langmead B, Salzberg SL. Nature Methods 9(4):357-359, 2012
- **fastp**: Chen S, et al. Bioinformatics 34(17):i884-i890, 2018
- **SAMtools**: Li H, et al. Bioinformatics 25(16):2078-2079, 2009
- **SciPy**: Virtanen P, et al. Nature Methods 17:261-272, 2020
