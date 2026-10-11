# single-cell

常规单细胞 RNA-seq skill。先用 `workflow/workflow-sc-analysis` 判断阶段，再按 Python、R 或跨语言桥接选择实现 skill。

| 分组 | Skill | 用途 |
|------|-------|------|
| workflow | `workflow-sc-analysis` | scRNA-seq 预处理、参考构建、h5ad 交接和下游边界。 |
| python | `sc-preprocessing-qc` | 数据接收、按捕获文库 cell calling/QC、双细胞与按需环境 RNA。 |
| python | `sc-integration-clustering` | counts/log/HVG/PCA、整合评估、邻居图、分辨率与 UMAP。 |
| python | `sc-marker-annotation` | 全可用基因 marker、竞争谱系、注释版本与身份/状态。 |
| python | `sc-sample-level-comparison` | 生物重复设计、edgeR pseudobulk、比例、富集与评分。 |
| python | `sc-cell-communication` | 按问题触发的配体受体预测与样本支持。 |
| python | `sc-subcluster-analysis` | 子集重聚类、marker 驱动注释、污染判断。 |
| python | `sc-velocity-paga-cellrank-analysis` | scVelo、PAGA、CellRank 轨迹和命运分析。 |
| python | `sc-pyscenic-regulon-analysis` | pySCENIC/AUCell regulon activity 分析。 |
| python | `sc-multiomics-tcr-vdj-analysis` | scRNA 与 TCR/BCR/VDJ 整合。 |
| r | `sc-r-monocle2-analysis` | R/Seurat 中 Monocle2 拟时序。 |
| shared | `h5ad-to-r-bridge` | AnnData h5ad 到 R/Seurat/spacexr/RCTD 的桥接。 |

主流程采用 Python，R 负责统计及专项方法。共享矩阵、唯一键与验收契约见 [data-contract.md](workflow/workflow-sc-analysis/references/data-contract.md)。靶向空间平台保留自身规则；专项无需每次全部运行。结构审计不能证明方法正确；无 annotation registry 的普通项目报告不适用，有 registry 的项目继续验证其契约。
