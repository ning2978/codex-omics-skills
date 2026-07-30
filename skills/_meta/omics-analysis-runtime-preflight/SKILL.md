---
name: omics-analysis-runtime-preflight
description: Audit single-cell and spatial-omics analysis readiness across exact Python environments, AnnData data contracts, NVIDIA/CUDA/PyTorch layers, and real synthetic GPU training for scVI or cell2location. Use before heavy scRNA-seq, Xenium, Visium, deconvolution, latent-model, or other accelerator-backed runs; when users ask whether a server GPU actually works; or when diagnosing imports, CUDA visibility, driver failures, hangs, OOM risk, and environment drift.
---

# Omics Analysis Runtime Preflight

## Role

Prove whether an analysis is runnable before launching a costly job. Keep four claims separate:

1. the physical accelerator exists;
2. the driver exposes it to user space;
3. the exact analysis environment can execute CUDA tensors;
4. the named model actually trains on that CUDA device.

Never infer a later claim from an earlier one. A CUDA-enabled PyTorch wheel, `lspci` entry, or successful import alone does not prove GPU training.

## Workflow

### 1. Fix the scope

- Identify the exact Python executable for every analysis environment. Test scVI and cell2location separately when they use different environments.
- Identify representative scRNA-seq, Xenium, or Visium h5ad inputs and their intended raw-count, batch, sample, label, and spatial-coordinate keys.
- Keep project paths in a repository-local adapter or report; do not encode them in this global skill.
- Use `skills-taxonomy` plus the relevant workflow or algorithm skill for method semantics. This skill owns runtime evidence only.

### 2. Inventory the host and each environment

Run the probe with the exact interpreter that will launch the analysis:

```bash
<analysis-python> <skill-dir>/scripts/runtime_probe.py \
  --require-gpu \
  --import-module scanpy \
  --import-module anndata \
  --import-module scvi \
  --output <report-dir>/runtime_scvi.json
```

For a cell2location environment, replace `scvi` or add `cell2location`. Run separate reports rather than merging package versions from different interpreters.

Treat timeouts as failures with evidence. Do not retry an indefinitely hanging CUDA probe inside the same process.

### 3. Validate representative AnnData contracts

Use backed reads and a bounded matrix sample:

```bash
<analysis-python> <skill-dir>/scripts/h5ad_contract.py \
  --input <representative.h5ad> \
  --kind xenium \
  --count-source layer:counts \
  --obs-key sample_roi \
  --output <report-dir>/xenium_h5ad_contract.json
```

Use `--kind scrna`, `xenium`, `visium`, or `generic`. For training inputs, require sampled values to be finite, nonnegative, and integer-like. For Xenium/Visium, require a finite two-column-or-wider `obsm["spatial"]`. Report that the integer check is sampled, not a proof over the full matrix.

### 4. Execute named-model smoke tests

Run scVI with a short timeout wrapper supplied by the host shell:

```bash
timeout 300 <scvi-python> <skill-dir>/scripts/scvi_smoke.py \
  --accelerator gpu --epochs 2 \
  --output <report-dir>/scvi_gpu_smoke.json
```

Run both cell2location training stages independently of scVI:

```bash
timeout 420 <cell2location-python> <skill-dir>/scripts/cell2location_smoke.py \
  --accelerator gpu --epochs 2 \
  --output <report-dir>/cell2location_gpu_smoke.json
```

GPU PASS requires the trainer root device to be CUDA and a nonzero CUDA allocation during training. If CUDA is unavailable, fail at the precondition and do not silently fall back to CPU. A CPU smoke test may isolate package/API health, but label it CPU and never use it as GPU evidence.

### 5. Report and hand off

Read `references/evidence-contract.md` before interpreting results or declaring readiness. Summarize each environment and model as PASS, WARN, FAIL, or SKIP, with:

- exact interpreter and package versions;
- physical GPU, driver, device-node, CUDA runtime, and tensor-probe evidence;
- representative data-contract results;
- requested and observed trainer device;
- elapsed time, finite-loss/latent checks, and peak allocated GPU memory;
- limitations and the next corrective action.

After PASS, hand off analysis parameters and full-run completion checks to the relevant modality or algorithm skill. A synthetic smoke test proves plumbing, not dataset-scale memory capacity, convergence, or biological validity.

## Safety and Environment Rules

- Keep the default workflow read-only except for explicit JSON reports.
- Do not start full-data training as a preflight.
- Do not install or change drivers or packages merely to make a probe pass. Diagnose first and request authority for a material system repair.
- If an environment must change, apply `environment-installation-policy`: use conda and pip only, never mamba or micromamba.
- Preserve raw counts and stable sample boundaries before any model handoff.
- Run one main cell2location training process per GPU unless allocation is deliberate.

## Routed Method Contracts

- Use `workflow-sc-analysis` for scRNA-seq stages and count/annotation handoffs.
- Use `workflow-st-analysis` for spatial stage and platform contracts.
- Use `algorithm-cell2location` plus `workflow-deconvolution-run` for cell2location parameters, launch, monitoring, posterior export, and final validation.
- Use the repository's Xenium integration adapter for project-specific scVI inputs, batch keys, checkpoints, and output validation.
