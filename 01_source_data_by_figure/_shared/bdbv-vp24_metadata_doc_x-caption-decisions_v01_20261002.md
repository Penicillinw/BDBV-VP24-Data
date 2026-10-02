# X_caption_decisions.md — 五条不一致的裁定与可直接粘贴的替换句

> 轮次 X（2026-10-02）｜角色 agent_analysis
> 对应产物：`X1_p83_distance_sweep.tsv`、`X2_contacts_authoritative.tsv`、
> `X3_figs2d_provenance.tsv`、`X4_md_system_inventory.tsv`、`X4_md_recomputed_metrics.tsv`、
> `X5_fig1a_denominator_recompute.tsv`（同一目录）。
> 复现：`python _coord/_round_X_20261002/analysis/x_run_all.py`（日志 `X_console.txt`）。

---

## X1｜Fig. 1b「55.0 Å from the second VP24 copy」

**裁定：删去该分句（该数字不可复算）。**

穷尽扫描 66 条度量（4U2X 的 A–B / A–C / B–C 三种拷贝组合；出厂同源模型
`vp24_kpna_complex.pdb` 的 A–B；`data/template_transfer_20260919/` 的两个 BDBV VP24 模型；
度量含任意重原子最小距离、CA–CA、CA–质心、残基质心–残基质心、整链质心–整链质心、
对应主链原子 N/C/O/CB、以及每个组合中 P83 到字面链 B 最近原子的距离）中，
**没有任何一条落在 55.0 ± 0.1 Å**。

「P83 到第二个 VP24 拷贝」字面定义的全部候选值：

| 定义（模型 = 出厂同源模型 A/B；4U2X = 实验坐标） | 值 (Å) |
| --- | ---: |
| 模型：P83 到链 B 最近重原子 | 22.902 |
| 模型：P83 CA–CA（第二拷贝同号残基） | 80.169 |
| 模型：P83 残基质心–B83 残基质心 | 80.505 |
| 模型：整链 A 质心–整链 B 质心 | 36.196 |
| 4U2X：链 A 质心–链 B 质心 | 48.548 |
| 4U2X：链 A 质心–链 C 质心 | 56.616 |
| 4U2X：链 B 质心–链 C 质心 | 45.776 |
| 4U2X：P83(A) CA–CA(B) / CA–CA(C) | 24.257 / 15.441 |

唯一落在 55.0 ± 0.1 Å 的工作区量是**第二个拷贝自身的 P83 到最近 KPNA5 原子的距离**：
模型链 B 的 P83 到链 Z 最近重原子 = **54.983 Å**（`X1_p83_distance_sweep.tsv` 行
`model(vp24_kpna_complex) | B-Z-per-copy-p83`）。也就是说，这个数字几乎可以肯定是
「copy 2 的 P83 距 KPNA5 约 55 Å」被误写成了「P83 距第二个拷贝 55 Å」。

**建议替换句（英文，可直接粘贴；与 v13 语气一致）：**

> dashed yellow lines mark six inter-chain polar contacts below 3.6 Å, and P83 lies 11.7 Å from the nearest KPNA5 atom in the model and 11.4 Å in the experimental coordinates, with no partner atom within 8 Å.

若作者希望保留一个第二拷贝的数字（备选，第二句可整段替换上式末句）：

> dashed yellow lines mark six inter-chain polar contacts below 3.6 Å, and P83 lies 11.7 Å from the nearest KPNA5 atom in the model and 11.4 Å in the experimental coordinates, with no partner atom within 8 Å; the second VP24 copy of the model is 36.2 Å away by chain centroid, and its own P83 lies 55.0 Å from the nearest KPNA5 atom.

---

## X2｜Fig. 1b 接触数 6 vs 7

**裁定：图注口径正确（6 条），PyMOL 渲染脚本口径过宽；改渲染脚本，不改图注计数。**

从零重算（不读任何既有接触表）：

* **面板自身的结构**（`figures/pymol_publication/complexes/vp24_kpna_complex.pdb`，
  chain A = BDBV 同源模型 / chain Z = KPNA5），四个变异位点 135/140/141/184 与配偶链之间的
  **重原子 N/O···N/O、< 3.6 Å** 接触恰为 **6 条**，与 `Fig1b_geometry.json` 逐条一致：

  | 病毒 | 病毒原子 | 配偶 | 配偶原子 | 距离 (Å) |
  | --- | --- | --- | --- | ---: |
  | Q135 | O | D480 | OD1 | 2.733 |
  | Q135 | OE1 | E483 | OE1 | 3.227 |
  | R184 | NH1 | K481 | NZ | 2.969 |
  | R184 | NE | K481 | NZ | 3.234 |
  | R184 | NH2 | K481 | NZ | 3.513 |
  | R184 | O | R398 | NH1 | 3.538 |

  归属为 Q135 两条、R184 四条；H140 与 A141 在 3.6 Å 内没有链间极性接触。
* **渲染脚本**（`figures/main_v1_20260926/_fig1b/render_fig1b.py`）的 7 条键是在
  **另一个结构**（4U2X 的实验 EBOV 侧链，位置 140 是 R 而不是 H）上、用**更宽的 4.0 Å**
  截断枚举的。在 4U2X 上按同一判据重算，四个位点在 < 4.0 Å 内是 **10 条**、在 < 3.6 Å 内是
  **8 条**——脚本列出的那 7 条既不是前者也不是后者（漏了 `135 O–480 OD1 3.203` 与
  `184 NE2–398 NH2 3.552`，却收了 `184 NE2–431 OD2 3.846`）。
* **争议对**：4U2X 中 `R140 NH2 ··· D475 O` = **2.950 Å**，两端是 N 与 O，**是**极性原子对，
  且 < 3.6 Å；但面板画的是 BDBV 模型，位置 140 是组氨酸，没有 NH2，因此这条键在
  所画结构里**不存在**。
* 参考量（非四变异位点的全界面）：4U2X 三条同源拷贝的极性对 < 3.6 Å 分别为 22 / 16 / 15 条；
  模型 A–Z 为 17 条（其中含一对非物理的 1.65 Å 接触 `185 ND2 ··· 434 OG1`，不得引用）。

**建议替换句（把判据与归属写进图注）：**

> dashed yellow lines mark the six inter-chain heavy-atom N/O···N/O contacts below 3.6 Å between the four variant positions and KPNA5 (two from Q135 and four from R184), and P83 lies 11.7 Å from the nearest KPNA5 atom in the model and 11.4 Å in the experimental coordinates, with no partner atom within 8 Å.

**代码侧动作（不经我改）**：把 `figures/main_v1_20260926/_fig1b/render_fig1b.py` 的 `BONDS`
换成上表 6 条（或把该脚本标为 superseded 后删除），使渲染输入与出厂面板同结构、同截断。

---

## X3｜Fig. S2d 需求矩阵的凭证化

**裁定：35 格中 34 格有凭证证实，1 格被凭证推翻——`342661 | Publicly released` 应改为 1。**

* 35 个 1/0 单元格全部给出 `voucher_file / voucher_field / voucher_value`，
  见 `X3_figs2d_provenance.tsv`；34 格与冻结网格一致，1 格反驳。
* 被推翻的格子：冻结网格给 `342661 = 0`，但凭证显示 GSE342661 的**处理后 TPM 矩阵已于
  2026-09-25 在 GEO 放行**（`GSE342661_counts_tpm.matrix.gz`，10,325,209 B，且已落在
  `analysis/gse342661_matrix_20260926/raw/`）；`analysis/gse342661_matrix_20260926/PUBLICATION_STATUS_20260926.md`
  亦明确记录「矩阵放行后，‘不可取得’的说法已成事实错误」。据此该格应为 1。
* 五个分数的逐项对应（冻结值 → 凭证校正后）：

  | 数据集 | 冻结分数 | 逐项（7 需求） | 校正后 |
  | --- | ---: | --- | ---: |
  | 114905 | 4/7 | BDBV✔ EBOV✔ same-cell✔ host-tx✘ IFN✘ ≥2rep✘ public✔ | 4/7 |
  | 309699 | 5/7 | BDBV✘ EBOV✔ same-cell✘ host-tx✔ IFN✔ ≥2rep✔ public✔ | 5/7 |
  | 226106 | 5/7 | BDBV✘ EBOV✔ same-cell✘ host-tx✔ IFN✔ ≥2rep✔ public✔ | 5/7 |
  | 46599 | 4/7 | BDBV✘ EBOV✘ same-cell✘ host-tx✔ IFN✔ ≥2rep✔ public✔ | 4/7 |
  | 342661 | 3/7 | BDBV✘ EBOV✘ same-cell✘ host-tx✔ IFN✔ ≥2rep✔ **public✘→✔** | **4/7** |

**建议替换句（英文）：**

> Panel d: 114905 (4/7), 309699 (5/7), 226106 (5/7), 46599 (4/7) and 342661 (4/7) requirements met; pale = absent or not applicable; 309699 = human keratinocytes (IFN-α/β/γ/λ, EBOV read-outs); no dataset satisfies all seven.

（同步需要把 `FigS2_check_and_gaps.py` 的 `panel_d` 网格第 5 列 `publicly released` 由 0 改为 1，
使图与图注一致。若作者坚持保留冻结时的 0，则必须在图注写明该格是 2026-09-25 之前的状态。）

---

## X4｜试点 MD 数据的保留裁定

**裁定：静态坐标判「入包」；由轨迹导出的量判「无法复算」，不得冒充。**

* `openmm_results.json` 报告的界面接触数可从同目录坐标**逐位复算**：以「重原子、跨链、
  残基对 < 4.5 Å」为口径，WT 55→53、N135A 55→52、R140A 53→52、4x_BDBVres 52→53、
  BDBV_interface5 49→55、BDBV_iface5_plus83 49→55、RESTV 46→49、SUDV 46→49，
  与 JSON **完全一致**。删掉坐标，这 8 组整数就无法复算 → 坐标必须入包。
* 静态可复算量：77–90 CA-RMSD 相对 WT（模型内自拟合 0.083–0.402 Å；全 CA 拟合
  0.137–0.207 Å）、0.8 nm 天然接触分数（天然集 240 残基对，1.2 nm 容许下全部体系 = 1.000）、
  界面接触数（同上）。
* **无法复算**（工作区只有最小化坐标，无任何轨迹文件）：`local_rmsd_nm`、`Q_fraction`、
  `minimised_energy_kcal_mol`、`void_volume_nm3`，`X4_md_recomputed_metrics.tsv` 中标
  `not_recomputable`（86 行中 41 行）。`TAFV_interface` 连最小化坐标都没有（力场报错），
  其 77–90 RMSD 亦为 `not_recomputable`。
* 附带登记（非本条要求）：Methods 的 500 步最小化属 pilot 脚本
  （`sim.minimizeEnergy(maxIterations=500)`），`openmm_results.json` 的 300 步属另一批单点扫描，
  二者不是同一次运行，不必改。

**图注/正文：no change。** 正文只用「两个 1 ns 副本未能区分 WT 与 S83」这类定性表述，
未落具体数；若要落具体数，必须先存轨迹。

---

## X5｜Fig. 1a 分母链

**裁定：no change——图注分母链可逐位复现。**

从 `analysis/vp24_pan_species_20260925/coordinator_crosscheck/vp24_six_species_positions.tsv`
从零重算（`x5_denominator.py`）：

| species | retrieved | excluded_fragment | assessed | uncalled_at_some_position | called_all_four |
| --- | ---: | ---: | ---: | ---: | ---: |
| BDBV | 32 | 1 | 31 | 0 | 31 |
| EBOV | 3302 | 31 | 3271 | 18 | 3253 |
| SUDV | 147 | 0 | 147 | 1 | 146 |
| TAFV | 4 | 0 | 4 | 0 | 4 |
| RESTV | 27 | 0 | 27 | 0 | 27 |
| BOMV | 10 | 0 | 10 | 0 | 10 |
| **TOTAL** | **3522** | **32** | **3490** | **19** | **3471** |

非 Bundibugyo 的 `called_all_four` = 3471 − 31 = **3440**。与图注
`3,522 / 32 / 3,490 / 19 / 3,471 / 3,440` **完全一致，无差值**。

**图注：no change。**
