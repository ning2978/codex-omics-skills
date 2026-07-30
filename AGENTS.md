# Repository instructions

These rules apply to the entire repository.

- Keep every published skill portable. Do not commit patient data, clinical
  metadata, analysis outputs, project-specific annotation registries, private
  cohort identifiers, local absolute paths, credentials, logs, or caches.
- Keep project-specific paths, labels, and frozen defaults in repository-local
  adapters outside this collection.
- Each skill must have exactly `name` and `description` in `SKILL.md`
  frontmatter. The folder name must match `name`.
- Keep `SKILL.md` concise and below 500 lines. Put detailed contracts in
  `references/` and deterministic implementations in `scripts/`.
- Keep `agents/openai.yaml` aligned with `SKILL.md`; its default prompt must
  explicitly mention `$skill-name`.
- Add, rename, or retire biomedical skills through
  `skills/_meta/skills-taxonomy/` and update its migration map when relevant.
- Regenerate `catalog.json` with `python scripts/build_catalog.py`.
- Run `python scripts/validate_repo.py` before committing.
- Never promote candidate, review, exploratory, unresolved, or
  pending-approval labels into formal annotations through a global skill.
- Treat segmented-cell Xenium expression reclustering as counts-first
  expression analysis. Do not conflate it with spatial-composition
  clustering such as RCN/CN, target-cell NMF, or CellCharter niches.
