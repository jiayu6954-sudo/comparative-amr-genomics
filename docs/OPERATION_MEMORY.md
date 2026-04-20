# 操作记忆库 — 跨物种比较 AMR 基因组学项目
> **Operation Memory Library** · 版本 1.0 · 2026-04-20  
> 路径: `e:/miniconda3/envs/llama-env/comparative_amr_project/`

---

## 一、项目核心坐标

| 项 | 值 |
|---|---|
| 研究问题 | 物种特异性IS生态学塑造不同的碳青霉烯耐药基因迁移策略（KPN vs ECO vs ABA比较） |
| 数据来源 | NCBI GenBank（Entrez API，`"China"[Country]`过滤） |
| 目标菌种 | *K. pneumoniae*（492）/ *E. coli*（266，SESSION-016最终）/ *A. baumannii*（168） |
| 地理范围 | 中国临床分离株（三物种全部） |
| 分析规模 | **926基因组**（SESSION-016最终，ECO扩充至266） |
| 目标期刊 | **Nature Microbiology (IF >25)** |
| KPN数据来源 | 直接复用`amr_project/`，无需重下载 |

---

## 二、核心数据结果（SESSION-016 最终 — 2026-04-20）

### 2A. 三物种队列概况

| 指标 | K. pneumoniae | E. coli | A. baumannii |
|------|--------------|---------|--------------|
| QC通过基因组 | **492** | **266** | **168** |
| 碳青霉烯耐药率 | **40.7%** (200/492) | **10.5%** (28/266) | **99.4%** (167/168) |
| 主要耐药基因 | KPC-2(54.5%), NDM(44%) | NDM-5(25), OXA-48(2), NDM-1(1), KPC-2(1) | OXA-23(83.2%), OXA-51(98.8%) |
| AMR扫描hits | 217 | 29 | 312（含内源OXA-51） |

### 2B. IS侧翼结构分析（SESSION-016 关键结果 — 重大数值变化）

| 指标 | K. pneumoniae | E. coli | A. baumannii（获得性） |
|------|--------------|---------|----------------------|
| 分析的AMR位点数 | 50 | 29 | 146（OXA-23/24/58/NDM） |
| 复合转座子率 | **94.0%** (47/50) | **62.1%** (18/29) | **6.8%** (10/146) |
| 主要IS家族 | IS26 (IS6) | IS6(23.4%), Tn3(20.4%), IS5(13.1%), IS1(9.5%) | ISAba (ISAba1) |
| ISAba上游偏向 | — | — | **92.3%** 上游位置 |
| IS-AMR配对总数 | 2,985 | 137 | 62 |
| 置信区间（Clopper-Pearson/Wilson） | [83.5–98.8%] | [44.0–77.3%] | [3.7–12.2%] |

⚠️ **重大发现变化（SESSION-016 vs SESSION-015）**:
- ECO复合转座子率：15.4% → **62.1%**（样本扩充+AMR基因多样性增加）
- ECO主要IS家族：IS5/IS1/Tn3 → **IS6为最多(23.4%)**，与Tn3/IS5/IS1共存
- ECO AMR基因：仅NDM-5 → NDM-5+OXA-48+NDM-1+KPC-2（多样化）

### 2C. IS负担预测分析（含PFAM HMM校正）

| 指标 | K. pneumoniae | E. coli (校正前) | E. coli (PFAM校正后) | A. baumannii |
|------|--------------|---------|---------|--------------|
| 焦点IS家族 | IS6 (IS26) | IS6 (IS26) | IS6 (IS26) | ISAba |
| 耐药中位数IS拷贝 | **19** | 8 | 11.5 | 6 |
| 敏感中位数IS拷贝 | 10 | 4 | 8 | 0 |
| AUC | **0.6657** | 0.773 | **0.718** | 0.982（n=1敏感，不可比） |
| Cliff's δ | **0.3314** | 0.546 | **0.436** | 0.964（统计无效） |
| Mann-Whitney p | **4.2×10⁻¹⁰** | 2.26×10⁻⁶ | **2.26×10⁻⁶** | 0.097（n=1敏感） |

**PFAM校正结果（ECO）**:
- 重分类IS_unknown→IS6: 1,056个（敏感基因组获得959个 vs 耐药97个 = 9.9×偏向敏感）
- 偏差方向：**CONFIRMED** — 敏感基因组IS6注释严重不足
- 校正后AUC 0.718仍高度显著，为可发表数字

**PFAM校正结果（ABA）**:
- 重分类: 12个（耐药10，敏感2），无偏差 → 不需要校正

### 2D. 物种间复合转座子率差异显著性（SESSION-016更新）

| 比较对 | Fisher's 精确检验 |
|--------|-----------------|
| KPN vs ECO | p = 3.2×10⁻³（94% vs 62.1%，OR=9.6） |
| KPN vs ABA（获得性） | p < 10⁻¹⁵（94% vs 6.8%） |
| ECO vs ABA | **p = 3.7×10⁻⁸**（62.1% vs 6.8%，OR=22.3）**（与SESSION-015 p=0.22 NS完全相反！）** |

---

## 三、科学解释框架（SESSION-016更新 — 重大叙事更新）

### 3A. 修正后的三物种框架（梯度模型）
- **KPN**: IS26双侧翼复合转座子（94.0%），IS26绝对主导 → 系统性AMR基因整合/拷贝数扩增
- **ECO**: IS26仍是最多IS家族(23.4%)，但与Tn3/IS5/IS1共存；复合转座子62.1% → 混合策略
- **ABA**: ISAba1上游插入（92.3% upstream位置），复合转座子仅6.8% → 启动子捕获机制

**梯度**: KPN(94%) > ECO(62.1%) > ABA(6.8%)，三个两两比较均显著

### 3B. IS26功能区别（更新）
- IS26 in KPN: **结构动员器绝对主导**（94%双侧翼复合转座子）
- IS26 in ECO: **结构动员器之一**（IS26+Tn3+IS5+IS1 四元IS生态，62%复合）
- ISAba1 in ABA: **转录激活器**（上游插入 → 提供P2启动子给OXA-23，92.3%上游）

### 3C. ECO IS6 AUC显著性解释（SESSION-016颠覆SESSION-015结论）
ECO IS6 burden **确实显著** (AUC=0.718 PFAM校正后，p=2.3×10⁻⁶)。
- SESSION-015中的AUC=0.560 NS是因为样本量不足（n=13 resistant），扩充至n=28后效应量清晰
- NDM-5/IncX3质粒strains贡献非复合转座子(不积累IS26)，但KPC-2/OXA-48 strains使用IS26复合转座子
- PFAM校正关键：敏感基因组IS6被系统性低估(959 vs 97重分类)，校正后AUC=0.718
- 核心结论：**IS6负担在肠杆菌科（KPN和ECO）均是显著预测指标**，ABA无法比较

---

## 四、数据文件坐标（SESSION-016更新）

| 文件 | 路径 | 说明 |
|------|------|------|
| ECO下载状态 | `data/raw/eco/download_status.tsv` | 387行（合并）；OK=257, PARTIAL=19, FAILED=111 |
| ABA下载状态 | `data/raw/aba/download_status.tsv` | 300行；OK=169, FAILED=131 |
| ECO manifest | `data/validated/eco_manifest.tsv` | **266 PASS** |
| ABA manifest | `data/validated/aba_manifest.tsv` | 168 PASS |
| ECO AMR hits | `data/processed/eco_amr_hits.tsv` | **29 hits**（NDM-5×25, OXA-48×2, NDM-1×1, KPC-2×1）|
| ABA AMR hits | `data/processed/aba_amr_hits.tsv` | 312 hits（OXA-51/23/24/58/NDM） |
| ECO IS context | `data/processed/eco_is_context.tsv` | **137 IS-AMR对** |
| ABA IS context | `data/processed/aba_is_context.tsv` | 62 IS-AMR对 |
| ECO locus class | `data/processed/eco_locus_classification.tsv` | **25分类位点；COMPOSITE=18** |
| ABA locus class | `data/processed/aba_locus_classification.tsv` | 22分类位点 |
| ECO PFAM结果 | `data/processed/eco_is_hmmer_results.tsv` | IS_unknown→IS6重分类结果 |
| ECO PFAM校正负担 | `data/processed/eco_is_burden_corrected.tsv` | 校正后IS6计数 |
| ECO PFAM校正统计 | `data/processed/eco_is_burden_corrected_stats.json` | AUC=0.718，p=2.3e-6 |
| ABA PFAM统计 | `data/processed/aba_is_burden_corrected_stats.json` | 无偏差确认 |
| 跨物种IS负担 | `data/processed/cross_species_burden.tsv` | 3物种合并 |
| 跨物种统计 | `data/processed/cross_species_stats.json` | KPN/ECO/ABA AUC等 |
| 4张图表 | `figures/fig1-4.pdf/.png` | 发表级，300 dpi，SESSION-016版本 |
| 论文草稿v3 | `reports/manuscript_comparative_v3.md` | **最终版，~5,400词** |
| 论文PDF v3 | `reports/manuscript_comparative_v3.pdf` | **285 KB**，Chrome headless渲染 |

---

## 五、环境与工具

```
Python版本: 3.11.9 (C:\Program Files\Python311\python)
Pandoc:    e:/miniconda3/bin/pandoc.exe (3.9.0.2)
Chrome:    C:\Program Files\Google\Chrome\Application\chrome.exe
pyhmmer:   需在 llama-env 中安装（pip install pyhmmer）

KPN项目引用路径:
  amr_project/data/validated/genome_manifest.tsv      ← KPN manifest
  amr_project/data/processed/amr_hits.tsv             ← KPN AMR hits
  amr_project/data/processed/is_context.tsv           ← KPN IS context
  amr_project/data/processed/is_burden_corrected_stats.json ← KPN HMMER stats
```

---

## 六、已踩坑记录（SESSION-014/015）

### Bug 1: FTP协议不被requests支持
- **现象**: `No connection adapters were found for 'ftp://...'`
- **根因**: Entrez返回`ftp://ftp.ncbi.nlm.nih.gov/...`路径，requests无FTP适配器
- **修复**: `if ftp_base.startswith('ftp://'): ftp_base = 'https://' + ftp_base[6:]`
- **位置**: `analysis/01_download.py` → `_download_genome()` 第一行

### Bug 2: KPN GFF路径列名不同
- **现象**: KPN IS负担全为0（AUC=0.5）
- **根因**: KPN manifest用`local_dir`列，非`gff_path`列；需glob找GFF
- **修复**: 
```python
local_dir = Path(row.get('local_dir', row.get('gff_path', '')))
if local_dir.is_dir():
    candidates = sorted(local_dir.glob('*_genomic.gff.gz'))
    if candidates: gff_path = candidates[0]
```
- **位置**: `analysis/05_is_burden.py` → `compute_kpn_burden()`

### Bug 3: 复合转座子率分母错误
- **现象**: fig1显示ECO=22.2%（实为2/9分类位点），应为15.4%（2/13总位点）
- **根因**: 分母用了locus_classification.tsv行数（9），而非总AMR位点数（13）
- **修复**: 同时读amr_hits.tsv，用`len(amr_df['accession'].unique())`作分母

### Bug 4: KPN硬编码比率错误
- **现象**: fig1 KPN显示78.3%（旧值）
- **修复**: 改为94.0% [CI 83.5–98.8%]（来自KPN项目Step 09验证结果）

### Bug 5: ECO靶向搜索(NDM/KPC关键词)返回0结果
- **现象**: NCBI Assembly的[All Fields]不索引AMR基因名
- **解决**: 改用扩大通用下载（800基因组），以5.5%自然耐药率获得~44耐药株

---

## 七、待完成任务（下一Session）

| 优先级 | 任务 | 状态 |
|--------|------|------|
| ✅ P1 | ECO扩充至266 QC通过 | **完成**（28 resistant，10.5%） |
| ✅ P2 | PFAM PF01527 HMM验证ECO IS6元件 | **完成**（偏差确认，AUC=0.718校正） |
| ✅ P3 | PFAM ISAba HMM验证ABA元件 | **完成**（无偏差） |
| ✅ P4 | 重新生成figures + manuscript v3 | **完成**（v3=5,400词+PDF=285KB） |
| **P5** | Zenodo新建存档 | 待执行（论文投稿前） |
| **P6** | MLST + 质粒复型分析（可选） | IncX3确认 NDM-5/ECO strains |
| **P7** | 增大ABA敏感株样本（n=1太少） | 其他物种补充或特殊搜索 |
