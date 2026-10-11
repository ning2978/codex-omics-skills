# _meta

全局路由和环境策略，不承载具体分析流程。

| Skill | 用途 |
|-------|------|
| `omics-analysis-runtime-preflight` | 精确 Python/counts/GPU 运行时与合成训练预检。 |
| `skills-taxonomy` | 按 single-cell / spatial、空间分辨率、分析阶段和算法选择 skill。 |
| `environment-installation-policy` | Conda/pip 环境与可重复性策略（勿用 mamba/micromamba）。 |

主分析 skill 已归入：

- `single-cell/`：常规单细胞 RNA-seq。
- `spatial/`：空间组学；先分 low-resolution / high-resolution，再按平台和算法路由。
