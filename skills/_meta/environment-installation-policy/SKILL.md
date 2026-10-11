---
name: environment-installation-policy
description: "Enforce the user's environment policy for installs or environment changes: use conda and pip; never use mamba or micromamba unless explicitly requested in the current turn."
---

# Skill: environment-installation-policy

## Purpose

Apply the user's environment installation preference consistently.

## Hard Rule

When creating, installing, updating, or documenting a software environment:

- Use `conda` for Conda environment creation, environment updates, and Conda package installation.
- Use `pip` for Python package installation when appropriate.
- Do not use `mamba` or `micromamba` for installation commands, environment creation, environment updates, examples, scripts, setup docs, or troubleshooting steps.

## Command Preferences

Prefer commands like:

```bash
conda create -n ENV_NAME python=3.11
conda install -n ENV_NAME PACKAGE
conda env create -f environment.yml
conda env update -n ENV_NAME -f environment.yml
pip install PACKAGE
```

Avoid commands like:

```bash
mamba create ...
mamba install ...
mamba env create ...
mamba env update ...
micromamba create ...
micromamba install ...
```

## When Existing Files Mention Mamba

If editing existing project files that already mention `mamba` or `micromamba`:

- Prefer replacing those instructions with equivalent `conda` commands when the task touches environment setup.
- If a file must retain historical or third-party text, do not add new `mamba` usage.
- If a tool output suggests `mamba`, translate the actionable command to `conda` before presenting it to the user.

## Conflict Handling

Only use `mamba` or `micromamba` if the user explicitly asks for it in the current turn. Otherwise, treat `conda` plus `pip` as the required installation method.
