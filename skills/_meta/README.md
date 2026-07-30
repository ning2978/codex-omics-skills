# _meta

全局路由和环境策略，不承载具体分析流程。

| Skill | 用途 |
|-------|------|
| `skills-taxonomy` | 按 single-cell / spatial、空间分辨率、分析阶段和算法选择 skill。 |
| `omics-analysis-runtime-preflight` | 审计 Python、AnnData、CUDA/PyTorch 与真实 GPU 训练可用性。 |
| `environment-installation-policy` | Conda/pip 环境与可重复性策略（勿用 mamba/micromamba）。 |

主分析 skill 已归入：

- `single-cell/`：常规单细胞 RNA-seq。
- `spatial/`：空间组学；先分 low-resolution / high-resolution，再按平台和算法路由。
