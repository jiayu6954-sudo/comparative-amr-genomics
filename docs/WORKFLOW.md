# 工作流文档 — 跨物种比较 AMR 基因组学分析管道
> **Pipeline Workflow Reference** · 版本 1.0 · 2026-04-20

---

## 总览

```
[NCBI GenBank — 3物种]
       │  Entrez API (esearch + esummary)
       │  "China"[Country] + assembly level filter
       ▼
┌─────────────────────────────────────────────────────┐
│  01_download.py  [--species eco|aba|both]           │
│  KPN: 直接复用 amr_project/                         │
│  ECO: 300→800检索 → 228+OK                          │
│  ABA: 300检索 → 169 OK                              │
│  FTP→HTTPS修复: ftp:// → https://                   │
└──────────────────────┬──────────────────────────────┘
                       │ data/raw/{eco|aba}/{accession}/
                       ▼
┌─────────────────────────────────────────────────────┐
│  02_validate.py  [--species eco|aba|both]           │
│  物种感知QC阈值（大小/N50/CDS/GFF可读性）           │
│  ECO: 237 PASS  ABA: 168 PASS                       │
└──────────────────────┬──────────────────────────────┘
                       │ data/validated/{short}_manifest.tsv
                       ▼
┌─────────────────────────────────────────────────────┐
│  03_amr_scan.py  [--species eco|aba|both]           │
│  Tier1: GFF3正则扫描AMR关键词                        │
│  Tier2: CARD phmmer验证（pyhmmer 0.12）              │
│  ECO: 13 NDM-5 hits / ABA: 312 hits (OXA-51/23/NDM) │
└──────────────────────┬──────────────────────────────┘
                       │ data/processed/{short}_amr_hits.tsv
                       ▼
┌─────────────────────────────────────────────────────┐
│  04_is_context.py  [--species eco|aba|both]         │
│  IS元件侧翼分类（±10kb窗口）                        │
│  → COMPOSITE_TRANSPOSON / SINGLE_IS_UP/DOWN / NO_IS │
│  ECO: 28 IS-AMR对  ABA: 62 IS-AMR对                │
└──────────────────────┬──────────────────────────────┘
                       │ data/processed/{short}_is_context.tsv
                       │ data/processed/{short}_locus_classification.tsv
                       ▼
┌─────────────────────────────────────────────────────┐
│  05_is_burden.py  [无参数，自动加载KPN]             │
│  全基因组焦点IS拷贝计数                              │
│  KPN(IS6)+ECO(IS6)+ABA(ISAba) → Mann-Whitney AUC   │
│  KPN: AUC=0.666 p=4.2e-10  ECO: AUC=0.560 p=0.47   │
└────────┬──────────────────────┬─────────────────────┘
         │                      │
         ▼                      ▼
┌──────────────────┐  ┌─────────────────────────────┐
│ 06_cross_species │  │  07_manuscript.py            │
│ _figures.py      │  │  论文草稿自动生成             │
│ Fig 1-4          │  │  reports/manuscript_v2.md    │
│ PDF + PNG        │  │  + 08_convert_to_pdf.py      │
└──────────────────┘  └─────────────────────────────┘
```

---

## 快速运行命令（全管道）

```bash
cd e:/miniconda3/envs/llama-env/comparative_amr_project

# Step 1: 下载（KPN跳过；ECO --limit 800 合并两批次约266 QC通过）
python analysis/01_download.py --species eco --limit 800
python analysis/01_download.py --species aba --limit 300

# Step 2: 质控
python analysis/02_validate.py --species both

# Step 3: AMR扫描
python analysis/03_amr_scan.py --species both

# Step 4: IS侧翼分析
python analysis/04_is_context.py --species both

# Step 5: IS负担 + 跨物种统计（自动加载KPN）
python analysis/05_is_burden.py

# Step 5b: PFAM HMM IS6/ISAba注释校正（科学严谨性必须步骤）
python analysis/09_is_hmmer_verify.py --species both

# Step 6: 图表
python analysis/06_cross_species_figures.py

# Step 7: 论文草稿（可选，SESSION-016已直接写v3）
# python analysis/07_manuscript.py

# Step 8: PDF转换（--input 路径会自动 resolve 为绝对路径）
python analysis/08_convert_to_pdf.py --input reports/manuscript_comparative_v3.md
```

---

## 物种配置关键参数

### K. pneumoniae（来自amr_project，直接复用）
| 参数 | 值 |
|------|---|
| 数据路径 | `e:/miniconda3/envs/llama-env/amr_project/` |
| 基因组大小 | 4.8–6.5 Mb |
| 焦点IS | IS6 (IS26, IS257, IS1353, IS1006) |
| 主要耐药基因 | blaKPC, blaNDM, blaIMP |
| 复合转座子率 | **94.0%** (47/50) |

### Escherichia coli
| 参数 | 值 |
|------|---|
| 下载目录 | `data/raw/eco/` |
| 基因组大小 | 4.2–5.8 Mb |
| 焦点IS | IS6 (IS26, IS257, IS1353, IS1006) |
| 主要耐药基因 | blaNDM-5 (IncX3质粒) |
| QC通过 | 237/300 (79%) |
| 耐药率 | 5.5% (13/237) |

### Acinetobacter baumannii
| 参数 | 值 |
|------|---|
| 下载目录 | `data/raw/aba/` |
| 基因组大小 | 3.7–4.5 Mb |
| 焦点IS | ISAba (ISAba1/2/3) |
| 主要耐药基因 | blaOXA-23 (获得性), blaOXA-51 (内源) |
| QC通过 | 168/169 (99.4%) |
| 耐药率 | 99.4% (167/168) |
| 注意 | OXA-51为内源标志物，不计入获得性复合转座子分析 |

---

## 输出文件清单

| 步骤 | 关键输出 | 位置 |
|------|---------|------|
| Step 1 | `download_status.tsv` | `data/raw/{short}/` |
| Step 2 | `{short}_manifest.tsv` | `data/validated/` |
| Step 3 | `{short}_amr_hits.tsv` | `data/processed/` |
| Step 4 | `{short}_is_context.tsv` `{short}_locus_classification.tsv` | `data/processed/` |
| Step 5 | `cross_species_burden.tsv` `cross_species_stats.json` | `data/processed/` |
| Step 6 | `fig1-4.pdf/.png` | `figures/` |
| Step 7 | `manuscript_comparative_v2.md` | `reports/` |
| Step 8 | `manuscript_comparative_v2.pdf` | `reports/` |

---

## IS侧翼分类规则

| 类型 | 定义 |
|------|------|
| COMPOSITE_TRANSPOSON | AMR基因±10kb内，上游AND下游都有IS元件 |
| SINGLE_IS_UPSTREAM | 仅上游有IS |
| SINGLE_IS_DOWNSTREAM | 仅下游有IS |
| NO_IS | ±10kb内无IS元件 |

**注意**: 复合转座子率的分母 = **总AMR位点数**（含NO_IS），非仅有IS的位点。

---

## PDF转换方法

```bash
python analysis/08_convert_to_pdf.py \
  --input reports/manuscript_comparative_v2.md
```

依赖：
- `e:/miniconda3/bin/pandoc.exe`（3.9.0.2）Markdown→HTML转换
- `C:/Program Files/Google/Chrome/Application/chrome.exe` 无头渲染→PDF

---

## 常见问题排查

| 现象 | 原因 | 解决 |
|------|------|------|
| FTP连接错误 | `ftp://`路径未转换 | 检查`_download_genome()`第一行：`if ftp_base.startswith('ftp://'): ftp_base = 'https://'+ftp_base[6:]` |
| KPN IS负担全为0 | manifest用`local_dir`非`gff_path` | `compute_kpn_burden()`中glob `*_genomic.gff.gz` |
| 复合转座子率偏高 | 分母用分类行数非总AMR位点 | 从amr_hits.tsv取`unique accession`数 |
| ABA AMR率99.4% | OXA-51为内源标志物 | 分析时区分OXA-51（内源）和OXA-23/24/58/NDM（获得性） |
| ECO NDM靶向搜索0结果 | NCBI Assembly不按AMR基因名索引 | 扩大通用下载量（800+）替代靶向搜索 |

---

## 统一理论 Perspective — 扩展管道

本项目数据汇入更大的统一IS理论框架：

```
ISS微生物组（genesis_project）                中国临床CRKP（amr_project + comparative_amr_project）
IS密度 0.115–0.312× terrestrial     ←────────────────────→     IS6 AUC=0.949，复合转座子率94%
        │                                                                  │
        └──────── unified_is_theory/manuscript_unified_v1.md ─────────────┘
                  "IS元件密度 = 基因组HGT压力计"
                  目标: Nature Microbiology Perspective / PNAS
```

**生成统一理论图表**:
```bash
cd e:/miniconda3/envs/llama-env/unified_is_theory
python analysis/fig1_generate.py
```

**生成统一理论PDF**:
```bash
python e:/miniconda3/envs/llama-env/comparative_amr_project/analysis/08_convert_to_pdf.py \
  --input e:/miniconda3/envs/llama-env/unified_is_theory/manuscript_unified_v1.md
```
