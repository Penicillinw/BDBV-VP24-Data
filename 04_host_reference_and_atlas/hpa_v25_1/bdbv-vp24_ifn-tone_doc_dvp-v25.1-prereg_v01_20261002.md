# 预注册：DVP × Δ 对齐（任务 `dvp_delta_alignment`）

**日期**：2026-09-27
**执行**：子代理 `/root/dvp_delta_alignment`
**开工背景**：本代理的 spawn 载荷未送达（团队级故障），按任务名 + 产物重建任务；
重建依据是 `analysis/protein_layer_scan_20260926/FINDINGS.md` §6 的"高价值后续"登记
与冻结文件 §7.22 的强制口径。
**本文件在任何对齐数值被计算之前写成**（§5.1）。事后改动只能追加到 §7 修订记录。

---

## 0 唯一要回答的问题

HPA v25.1 **Deep Visual Proteomics**（DVP；27 个具名细胞类型 × 14 组织，单一健康供者，
逐细胞类型 MS 强度，同一文件带 `Matched nCPM`）能否为 Δ（154 类等权基线 IFN 基调轴）
提供**细胞类型分辨的蛋白层检验**？若不能，**缺在哪一步、缺多少**。

**不回答**：DVP 是否"证实/否证"了限制梯度；不做任何新的生物学断言。

---

## 1 输入（只读）

| 件 | 路径 | 用途 |
| --- | --- | --- |
| Δ 的基因层输入 | `analysis/t4_ifn_landscape_20260925/hpa_ifn_landscape_wide.tsv` | 154 类 × 29 基因 nCPM（HPA 单细胞 Consensus） |
| Δ 本身 | `analysis/t4_ifn_landscape_20260925/restriction_gradient.tsv` + `restriction_gradient_summary.json` | Δ 与三分量、冻结变换口径 |
| DVP 汇总（细胞类型层） | `analysis/protein_layer_scan_20260926/raw/hpa_dvp/dvp_cell_type.tsv.zip` | 基因 × 27 类蛋白强度 |
| DVP 组数据（蛋白 + 匹配 mRNA） | `analysis/protein_layer_scan_20260926/raw/hpa_dvp/dvp_cell_type_group_data.tsv.zip` | 列：`Gene / Gene name / Cell type name / Intensity / Matched nCPM` |
| DVP 样本层 | `analysis/protein_layer_scan_20260926/raw/hpa_dvp/dvp_sample_data.tsv.zip` | 复孔阳性率（只在需要时用） |
| 既有交叉表 | `analysis/protein_layer_scan_20260926/dvp_vs_delta_match.tsv` | **只用于比对**，不作为我的主结果 |

**Δ 的冻结构造**（来自 `restriction_gradient_summary.json`，不得改动）：
`IFN_I_capacity = mean(z(IFNAR1), z(IFNAR2))`；
`IFN_III_capacity = mean(z(IFNLR1), z(IL10RB))`；
`ISG_priming = mean(z(ISG15), z(MX1))`；`Δ = mean(三分量)`；
变换 `log10(nCPM + 1)`，在 **154 类上** z 标准化。

---

## 2 冻结判据（先于数值）

### 2.1 细胞类型交叉表（P1）

1. **主口径**：我自己实现一个规范化匹配（小写、去标点、去空白、单复数归一），
   要求**双向唯一**（一个 DVP 名字只能配一个 Δ 名字，反之亦然）；并列候选一律判 `AMBIGUOUS`。
2. **第二实现核对**：与既有 `dvp_vs_delta_match.tsv` 逐行比对。既有权口径记
   exact 10 / normalised 11 / matched 21 / unmatched 6。**不一致必须逐条列出并裁定**，
   不得默默采用任一方。
3. **未匹配者的处理**：**不进入主分析**；只做敏感性分析（若纳入，须显式声明）。
   未匹配者必须**逐条给出不匹配的理由**（不是"没有对应"，而是"Δ 里确实没有这一族"或"名字族不同"）。
4. 允许的裁定：`exact` / `normalised` / `unmatched`。**不得**用 `loose(jaccard)` 类门槛凑数——
   0.25–0.50 的 Jaccard 匹配（既有表里的 `distal tubules`、`pancreas epithelial cells`）一律**降级为不进入主分析**。

### 2.2 检验 1：DVP 匹配 mRNA ↔ Δ 的 HPA mRNA（P2）

在**匹配对**上，对 6 个 Δ 焦点基因（`IFNAR1`、`IFNAR2`、`IFNLR1`、`IL10RB`、`ISG15`、`MX1`）
算 Spearman ρ（DVP `Matched nCPM` vs Δ 侧 nCPM）。

- **先做同源判定（强制，先于任何相关性）**：统计两列在重叠位点上**逐位相等**的比例。
  若 ≥ 0.95 判 `SAME-SOURCE`——此时 ρ 高**不构成独立验证**，只能表述为"同一来源的内部一致性"。
  若 < 0.95 判 `DISTINCT-SOURCE`，才可谈跨矩阵一致性。
- 判据（仅在 `DISTINCT-SOURCE` 时才用于结论）：ρ ≥ 0.5 `CONCORDANT`；0.2 ≤ ρ < 0.5 `WEAK`；ρ < 0.2 `DISCORDANT`。
- **全网格**：6 个基因**全部**报告，含 ρ、n、精确/渐近 p，以及未检出导致的缺位。

### 2.3 检验 2：在 DVP 矩阵上重算 Δ（P3）

用 DVP `Matched nCPM` 按**同一冻结公式**在 **27 类 DVP 细胞类型上**重算 `Δ'`，
与 Δ 在匹配对上的取值比 **Spearman 秩序一致**（不比数值——z 的基线不同，154 类 vs 27 类）。

- **零分布（§5.4 强制）**：匹配对内部标签置换 **10,000 次**，报观测 ρ 的置换 p 与零分布分位。
- 判据：ρ ≥ 0.5 且置换 p < 0.05 → `ALIGNED`；ρ ≥ 0.5 但 p ≥ 0.05 → `ALIGNED-BUT-UNDERPOWERED`；
  0.2 ≤ ρ < 0.5 → `PARTIAL`；ρ < 0.2 → `NOT-ALIGNED`。
- 顺带报 `IFN_I_capacity` / `IFN_III_capacity` / `ISG_priming` 三个分量各自的 ρ（诊断用）。

### 2.4 检验 3：蛋白层（P4）——本任务的核心

对每个焦点基因，**预置可检验性门槛**：蛋白 `Intensity` 必须在**匹配对上 ≥ 5 个细胞类型**非空。

- 不满足 → `NOT-TESTABLE (coverage)`，并报**实际非空数**与**全体 27 类里的非空数**。
- 满足 → 报 ρ(蛋白, 匹配 mRNA) 与 ρ(蛋白, Δ)。
- **本条的预期落点**（写在这里以免事后找台阶）：由 DVP 汇总件已知
  `IFNLR1 0/27`、`IFNAR2 0/27`、`IFNAR1 1/27`、`IL10RB 8/27`。
  → 四个受体里**至多 IL10RB 可能过门槛**，`IFNAR1` 几乎必然 `NOT-TESTABLE`。
  若结果确实如此，**结论必须写成"蛋白层无法检验 Δ 的受体臂"**，而不是"蛋白层支持 Δ"。

### 2.5 功效与材料（§5.3、§5.6 强制）

- **n**：匹配对数（≤ 21）。报 `min_detectable_rho(n, power=0.80, α=0.05)`。
- **材料**：DVP 为**单一健康供者**、仅 14 组织；此必须与所有结论同句出现。
- **深度**：DVP 为深度受限 MS，"蛋白未检出 ≠ 不表达"；与所有蛋白层结论同句出现。
- **版本**：HPA v25.1 + 2026 预印本（DOI `10.64898/2026.05.26.727663`），**未经同行评议**。

---

## 3 允许的最终 verdict 词表

`ALIGNED` / `ALIGNED-BUT-UNDERPOWERED` / `PARTIAL` / `NOT-ALIGNED` /
`NOT-TESTABLE (coverage)` / `SAME-SOURCE` / `AMBIGUOUS` / `DISCORDANT`。
每条结论必须绑定：n、口径、产物文件。

## 4 明确不做（边界）

1. **不改** Δ 的公式、基因集、权重（冻结 §2 关于 Δ 的规则）。
2. **不写**"蛋白层缺口已闭合"、"Δ 已在蛋白层被独立验证"（冻结 §7.22 明文禁止）。
3. **不碰** `tools/`、`manuscript/`、`figures/`、`report/`、`data/`、
   `研究进展与计划.md`、冻结文件、其他 `analysis/*` 目录。
4. **不重建**任何 docx；本任务的产物只有文本与表格。
5. **不用** PDF 文本做串检查（项目已有教训）。

## 5 交付物（全部落在 `analysis/dvp_delta_alignment_20260926/`）

| 文件 | 内容 |
| --- | --- |
| `PREREG.md` | 本文件 |
| `crosswalk.tsv` | DVP 27 类 ↔ Δ 154 类，含 `match_kind` / `rationale` / `in_primary` |
| `p1_crosswalk_summary.json` | 计数与不一致清单 |
| `p2_mrna_concordance.tsv` | 6 基因 × {ρ, n, p, same-source 比例, verdict} |
| `p3_delta_recompute.tsv` / `.json` | Δ' vs Δ：ρ、置换 p、零分布分位、三分量诊断 |
| `p4_protein_layer.tsv` | 每焦点基因的可检验性与结果 |
| `FINDINGS.md` | 30 秒摘要 + 全网格 + verdict + 对手稿 §7.22 的落点建议 |
| `repro.md` | 命令 + 输入 SHA-256 + 时间 |

## 6 写盘纪律自证

所有写路径以 `analysis/dvp_delta_alignment_20260926/` 为前缀；脚本内不含指向
`tools/`、`manuscript/`、`figures/`、`report/`、`data/` 或冻结文件的写路径。

## 7 修订记录

| # | 时间 | 内容 |
| --- | --- | --- |
| v1 | 2026-09-27 00:2x | 初版（在任何对齐数值计算之前写成） |
| v2 | 2026-09-27 00:5x | **两处修订，均在最终结论写成之前**：① **交叉表加一族**——v1 只允许 1:1；实测 24 个 DVP 组里有 4 组的 Δ 对应是**一对多**（`neurons`→3 类、`collecting ducts`→2 类、`keratinocytes`→2 类、`alveolar cell types`→2 类），1:1 会把它们全部丢掉且 n 只剩 9、功效不足。故新增 `family`（DVP token 集是 Δ 名字 token 集的**真子集**，家族大小 ≤3，Δ 侧取**家族均值**）作为 `sensitivity` 面板；`primary` 仍为 1:1 严格口径，三档**全部并列报告**。② **同源闸门补一条**：v1 的二元闸门（≥95% 逐位相等 → SAME-SOURCE）在全基因组尺度上判为 `DISTINCT-SOURCE`（仅 4.9% 逐位相等），但**一子集细胞类型确实逐位相等**（cardiomyocytes / granulosa cells / hepatocytes / oocytes / pancreatic islets，以及若干 0.0 项），说明两者**共享来源**。故最终 verdict 改为三态：`SAME-SOURCE` / `SHARED-PROVENANCE-DISTINCT-PROCESSING` / `DISTINCT-SOURCE`，并**禁止**把任何一种写成"独立队列验证"。 |
| v3 | 2026-09-27 01:0x | 追加记录，不改判据：发现 **DVP 有两个分辨率不同的文件**（`dvp_cell_type.tsv` = 27 细胞类型、只有蛋白；`dvp_cell_type_group_data.tsv` = **24 组**、有蛋白 + `Matched nCPM`）。对齐只能在 **24 组**上做。另：**第 0 列是 Ensembl id，基因符号在第 1 列**。 |
