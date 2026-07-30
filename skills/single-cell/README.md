# single-cell

常规单细胞 RNA-seq skill。先用 `workflow/workflow-sc-analysis` 判断阶段，再按 Python、R 或跨语言桥接选择实现 skill。

| 分组 | Skill | 用途 |
|------|-------|------|
| workflow | `workflow-sc-analysis` | scRNA-seq 预处理、参考构建、h5ad 交接和下游边界。 |
| python | `sc-subcluster-analysis` | 子集重聚类、marker 驱动注释、污染判断。 |
| python | `sc-velocity-paga-cellrank-analysis` | scVelo、PAGA、CellRank 轨迹和命运分析。 |
| python | `sc-pyscenic-regulon-analysis` | pySCENIC/AUCell regulon activity 分析。 |
| python | `sc-multiomics-tcr-vdj-analysis` | scRNA 与 TCR/BCR/VDJ 整合。 |
| r | `sc-r-monocle2-analysis` | R/Seurat 中 Monocle2 拟时序。 |
| shared | `h5ad-to-r-bridge` | AnnData h5ad 到 R/Seurat/spacexr/RCTD 的桥接。 |
