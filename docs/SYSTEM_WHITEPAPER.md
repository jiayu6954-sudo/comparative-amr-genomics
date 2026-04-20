# 系统白皮书 — 跨物种比较 AMR 基因组自动化科学发现系统
> **System Whitepaper: Cross-Species Comparative AMR Genomics Pipeline**  
> 版本 1.0 · 2026-04-20 · 作者: ZJY  
> 项目路径: `e:/miniconda3/envs/llama-env/comparative_amr_project/`

---

## 一、系统定位与价值主张

### 1.1 解决什么问题

IS元件介导的碳青霉烯耐药基因传播研究长期局限于单物种分析，缺乏跨物种的系统比较。三个关键问题从未被同时回答：

1. IS介导的复合转座子是所有革兰氏阴性菌的普适机制，还是物种特异性现象？
2. IS6(IS26)的双侧翼括弧功能与ISAba1的启动子激活功能在基因组规模有何量化差异？
3. IS拷贝数作为耐药预测生物标志物，是否可以跨物种推广？

**本系统用8个脚本、自动化分析897个基因组，首次同时回答了这三个问题。**

核心发现颠覆了"复合转座子是普适机制"的预设假说，揭示了深刻的物种特异性IS生态学分歧。

### 1.2 系统能做什么

```
输入: NCBI GenBank（三物种中国临床基因组）
输出: 跨物种比较论文 + 4张发表级图表 + PDF

自动化流程:
  ① 三物种Entrez检索下载（含缓存/断点续传）
  ② 物种感知质控（种特异性基因组大小/N50/CDS阈值）
  ③ 多层次AMR基因检测（GFF3扫描 + CARD验证）
  ④ IS元件侧翼架构分类（COMPOSITE/SINGLE/NO_IS）
  ⑤ 跨物种IS负担定量（Mann-Whitney AUC + Cliff's δ）
  ⑥ 4张发表级跨物种比较图表
  ⑦ 论文草稿自动生成（数字实时替换）
  ⑧ pandoc+Chrome PDF渲染
```

### 1.3 系统边界

- 不做基因组组装（依赖NCBI已组装基因组）
- 不做实验室验证（MIC、接合实验等）
- IS检测依赖GFF3注释质量（有偏差，需PFAM HMM纠正）
- OXA-51为ABA内源标志，与获得性AMR分析需区分

---

## 二、系统架构

### 2.1 目录结构

```
comparative_amr_project/
├── analysis/
│   ├── config.py                    # 三物种中心配置
│   ├── 01_download.py               # Entrez下载（ECO/ABA）
│   ├── 01b_download_eco_cr.py       # ECO耐药株靶向补充（Entrez All Fields）
│   ├── 02_validate.py               # 物种感知质控
│   ├── 03_amr_scan.py               # AMR扫描+CARD验证
│   ├── 04_is_context.py             # IS侧翼分析
│   ├── 05_is_burden.py              # IS负担+AUC（三物种合并）
│   ├── 06_cross_species_figures.py  # 4张跨物种图表
│   ├── 07_manuscript.py             # 论文草稿生成器
│   └── 08_convert_to_pdf.py        # pandoc+Chrome PDF转换
├── data/
│   ├── raw/eco/          # ECO下载基因组
│   ├── raw/aba/          # ABA下载基因组
│   ├── validated/        # QC通过manifest
│   └── processed/        # 分析结果
├── figures/              # fig1-4.pdf/.png
├── reports/              # manuscript_v2.md + .pdf
├── logs/                 # 各步骤运行日志
└── docs/                 # 本文档目录
    ├── OPERATION_MEMORY.md
    ├── WORKFLOW.md
    └── SYSTEM_WHITEPAPER.md
```

### 2.2 数据流

```
NCBI GenBank
    │
    ▼ Entrez esearch + esummary
01_download.py ──→ data/raw/{eco,aba}/{accession}/{gff,fna,stats}.gz
    │
    ▼ 物种特异性QC
02_validate.py ──→ data/validated/{short}_manifest.tsv
    │
    ▼ 两层检测
03_amr_scan.py ──→ data/processed/{short}_amr_hits.tsv
    │
    ▼ IS窗口分析
04_is_context.py ──→ {short}_is_context.tsv + {short}_locus_classification.tsv
    │
    ▼ 跨物种合并（+KPN来自amr_project/）
05_is_burden.py ──→ cross_species_burden.tsv + cross_species_stats.json
    │
    ├──→ 06_figures ──→ figures/fig1-4.pdf/.png
    └──→ 07_manuscript → reports/manuscript_v2.md
                              └──→ 08_pdf → manuscript_v2.pdf
```

---

## 三、关键科学发现（SESSION-016 — 最终版）

### 3.1 核心比较结果

| 维度 | K. pneumoniae | E. coli | A. baumannii |
|------|--------------|---------|--------------|
| 基因组数 | 492 | **266** | 168 |
| 耐药率 | 40.7% | **10.5%** | 99.4% |
| 复合转座子率（获得性AMR位点） | **94.0%** | **62.1%** | **6.8%** |
| IS家族（最多） | IS26 | IS26(23%), Tn3(20%), IS5(13%), IS1(9.5%) | ISAba1 |
| IS功能类型 | 结构动员 | 多元（IS26结构+IncX3平台） | 转录激活 |
| IS burden AUC (PFAM校正) | 0.667 ✓ | **0.718 ✓** | N/A (n=1 susc) |
| IS burden p | 4.2×10⁻¹⁰ | **2.3×10⁻⁶** | 0.097 |

### 3.2 论文叙事（SESSION-016最终版 — 叙事更新）

**标题**: *IS Element Ecology Shapes Divergent Carbapenem Resistance Gene Mobilisation Strategies Across Clinical Gram-Negative Pathogens: A Three-Species Comparative Genomic Analysis*

**梯度模型（更新）**:
- **KPN + ECO（肠杆菌科）**: 共享IS26复合转座子（94% vs 62%），IS6负担均显著
- **ABA**: ISAba1上游启动子（92.3% upstream）— 与肠杆菌科完全不同
- **核心梯度**: KPN(94%) > ECO(62.1%) > ABA(6.8%)，三对比较均显著

**科学价值（更新）**: 
1. IS26复合转座子架构是**肠杆菌科共享策略**，非KPN专有
2. ECO IS生态更多样（IS26+Tn3+IS5+IS1），ABA独特使用启动子捕获
3. IS6负担是**肠杆菌科（KPN+ECO）的通用耐药预测标志物**
4. PFAM校正是ECO IS6分析的必要步骤（9.9×注释偏差）

---

## 四、技术实现要点

### 4.1 Entrez查询构建
```python
query = (f'"{species_name}"[Organism] AND "China"[Country] AND '
         f'latest[filter] AND '
         f'("Complete Genome"[Assembly Level] OR '
         f'"Chromosome"[Assembly Level] OR '
         f'"Scaffold"[Assembly Level])')
```

### 4.2 FTP→HTTPS转换（关键修复）
```python
# 在_download_genome()最开头
if ftp_base.startswith('ftp://'):
    ftp_base = 'https://' + ftp_base[6:]
```

### 4.3 KPN数据复用
```python
# config.py
KPN_PROJECT  = Path('e:/miniconda3/envs/llama-env/amr_project')
KPN_MANIFEST = KPN_PROJECT / 'data/validated/genome_manifest.tsv'
# manifest使用'local_dir'列，需glob GFF文件
local_dir = Path(row.get('local_dir', ''))
if local_dir.is_dir():
    gff_path = sorted(local_dir.glob('*_genomic.gff.gz'))[0]
```

### 4.4 物种感知QC阈值
| 物种 | 基因组大小 | N50 | CDS |
|------|-----------|-----|-----|
| KPN | 4.8–6.5 Mb | ≥50 kb | ≥3,000 |
| ECO | 4.2–5.8 Mb | ≥50 kb | ≥3,000 |
| ABA | 3.7–4.5 Mb | ≥50 kb | ≥2,500 |

### 4.5 ABA内源vs获得性AMR区分
```python
# ABA OXA-51为内源（所有ABA均有），不计入复合转座子分析
ACQUIRED_ABA = ['OXA-23', 'OXA-24', 'OXA-40', 'OXA-58', 'NDM-1', 'NDM-5']
# 复合转座子率分母 = 获得性AMR位点数（n=146），非总hits（n=312）
```

### 4.6 PDF转换
```python
# pandoc MD→HTML，Chrome headless HTML→PDF
pandoc = 'e:/miniconda3/bin/pandoc.exe'
chrome = 'C:/Program Files/Google/Chrome/Application/chrome.exe'
# Chrome命令
['--headless', '--print-to-pdf=output.pdf', '--no-pdf-header-footer',
 '--virtual-time-budget=8000', f'file:///tmp_file.html']
```

---

## 五、与其他项目的关系

| 项目 | 路径 | 关系 |
|------|------|------|
| Genesis ISS项目 | `genesis_project/` | 独立；共享pandoc/Chrome PDF方法 |
| AMR K. pneumoniae项目 | `amr_project/` | 本项目**复用**KPN基因组数据 |
| 本跨物种项目 | `comparative_amr_project/` | 扩展AMR项目 + 三物种比较 |
| **统一IS理论** | `unified_is_theory/` | 本项目 + AMR项目 + ISS项目数据**汇总**为Perspective论文 |

---

## 五-A、统一IS理论框架（2026-04-20 新增）

**核心命题**: IS元件密度是"基因组水平基因转移压力计"（Genomic HGT Pressure Gauge）

| 数据来源 | HGT压力端 | IS密度观测 | 来源论文 |
|---------|----------|----------|---------|
| ISS表面微生物 | **极低**（封闭/无抗生素/少噬菌体） | 0.115–0.312× terrestrial | Genesis/ISS项目 |
| 中国临床敏感株 | 中等 | 参考基线（1.0×） | amr_project |
| 中国临床CRKP/CREC | **极高**（高抗生素压力，HGT密集） | IS6 AUC=0.949；复合转座子94% | comparative + amr项目 |

**统一理论论文**:
- 路径: `e:/miniconda3/envs/llama-env/unified_is_theory/manuscript_unified_v1.md`
- 格式: Nature Microbiology Perspective (~2,931词)
- 图表: `unified_is_theory/figures/fig1_is_pressure_gauge.pdf` + `fig2_is_ecology_panels.pdf`
- 投稿策略: 两篇主文投稿/预印本后，Perspective并行投稿

---

## 六、Zenodo存档计划

| 项目 | DOI | 状态 |
|------|-----|------|
| Genesis ISS项目 | `https://doi.org/10.5281/zenodo.19638104` | ✅ 已发布 |
| AMR K. pneumoniae项目 | 待申请 | 论文投稿前 |
| **本跨物种项目** | **https://doi.org/10.5281/zenodo.19665193** | ✅ DOI已预留，发布前最后确认 |
| 统一IS理论 | 与主文同步 | 主文发表后 |

**comparative_amr_project Zenodo元数据（待填写）**:
```
Title: IS Element Ecology Shapes Divergent Carbapenem Resistance Gene Mobilisation Strategies
       Across Clinical Gram-Negative Pathogens
Authors: ZJY
Description: Analysis code, processed data, and manuscript for three-species comparative
             genomic study of IS element architecture (926 genomes: KPN/ECO/ABA)
Keywords: insertion sequences, IS26, carbapenem resistance, K. pneumoniae, E. coli, A. baumannii
License: MIT
Related identifiers: https://github.com/jiayu6954-sudo/comparative-amr-genomics
Upload: figures/*.pdf, reports/manuscript_comparative_v3.pdf, requirements.txt
```

---

## 七、参考文献（系统设计）

1. Harmer CJ, Hall RM. IS26-mediated formation of transposons. *mSphere*. 2016.
2. Mugnier PD, et al. Phylogeny of ISAba1. *PLoS One*. 2009.
3. Corvec S, et al. ISAba1 insertion upstream of OXA-51. *AAC*. 2007.
4. NCBI Entrez Programming Utilities (eutils.ncbi.nlm.nih.gov)
5. CARD v3.3: McArthur AG, et al. *Antimicrob Agents Chemother*. 2013.
