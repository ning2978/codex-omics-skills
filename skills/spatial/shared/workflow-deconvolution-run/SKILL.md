---
name: workflow-deconvolution-run
description: Use this skill for algorithm-agnostic deconvolution execution flow, including run preparation, launch, monitoring, and completion checks.
---

# Skill: workflow-deconvolution-run

## Purpose

Use this skill to manage the execution lifecycle of deconvolution jobs, independent of algorithm.

Stages:

1. preflight contract check
2. launch strategy
3. runtime monitoring
4. completion validation

## Stage 1: Preflight

- confirm inputs exist and are versioned
- confirm output directory and log directory
- confirm selected algorithm skill (`spatial/low-resolution/visium-spot/algorithm-cell2location`, `spatial/low-resolution/visium-spot/algorithm-rctd`, etc.)
- confirm environment activation command

## Stage 2: Launch

- prefer non-interactive CLI command
- long runs on server should use `nohup` or scheduler
- print PID and log path at launch

## Stage 3: Monitor

- check process alive
- verify log growth
- detect hard failures (`Traceback`, `Killed`, `OOM`, `KeyboardInterrupt`)
- avoid duplicate concurrent runs unless explicitly requested

## Stage 4: Validate Completion

- main expected artifacts exist
- summary table exists
- logs contain completion markers or final save events

## Generic Deconvolution Contract

For Visium-like deconvolution, keep these paths explicit:

- QCed spatial `.h5ad` input
- single-cell reference `.h5ad` input
- label key and batch key in the reference
- output root
- per-sample result object, commonly `spatial.h5ad`
- run summary table

Prefer a CLI-style runner with explicit arguments over notebook-only execution for server runs.

## Related Skills

- `spatial/low-resolution/visium-spot/algorithm-cell2location`
- `spatial/low-resolution/visium-spot/algorithm-rctd`
- `spatial/shared/workflow-st-analysis`
- `single-cell/workflow/workflow-sc-analysis`
