#!/usr/bin/env python3
"""Audit active biomedical skills and optional repository-local adapters."""

from __future__ import annotations

import argparse
import itertools
import re
import sys
from dataclasses import dataclass
from pathlib import Path


NAME_RE = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")
SCRIPT_RE = re.compile(r"scripts/[A-Za-z0-9_./-]+\.(?:py|sh|R)")
REFERENCE_RE = re.compile(r"references/[A-Za-z0-9_./-]+")
WORD_RE = re.compile(r"[a-z][a-z0-9_+-]{2,}")
STOP_WORDS = {
    "after",
    "analysis",
    "and",
    "before",
    "from",
    "input",
    "into",
    "only",
    "output",
    "skill",
    "that",
    "the",
    "this",
    "use",
    "when",
    "with",
    "workflow",
}
RETIRED_ACTIVE_NAMES = {
    "spatial-project-analysis",
    "cell2location-params-official",
    "xenium-scimap-rcn",
    "xenium-scimap-nmf-neighborhood",
    "xenium-harmony-integration",
    "xenium-scvi-integration",
    "xenium-cellcharter-celltype",
    "xenium-cellcharter-latent",
    "xenium-radius-sensitivity",
}


@dataclass(frozen=True)
class Skill:
    scope: str
    path: Path
    name: str
    description: str
    text: str
    line_count: int


def parse_frontmatter(path: Path) -> tuple[str, str, set[str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0].strip() != "---":
        raise ValueError("missing opening YAML delimiter")
    try:
        end = next(i for i, line in enumerate(lines[1:], 1) if line.strip() == "---")
    except StopIteration as exc:
        raise ValueError("missing closing YAML delimiter") from exc

    block = lines[1:end]
    keys = {
        match.group(1)
        for line in block
        if not line.startswith((" ", "\t"))
        if (match := re.match(r"^([A-Za-z0-9_-]+):", line))
    }

    name = ""
    description = ""
    for index, line in enumerate(block):
        if line.startswith("name:"):
            name = line.split(":", 1)[1].strip().strip('"').strip("'")
        if line.startswith("description:"):
            value = line.split(":", 1)[1].strip()
            if value in {">", ">-", "|", "|-"}:
                folded: list[str] = []
                for continuation in block[index + 1 :]:
                    if continuation and not continuation.startswith((" ", "\t")):
                        break
                    if continuation.strip():
                        folded.append(continuation.strip())
                description = " ".join(folded)
            else:
                description = value.strip('"').strip("'")
    return name, description, keys


def discover(global_root: Path, project_root: Path | None) -> list[Skill]:
    candidates: list[tuple[str, Path]] = []
    for category in ("_meta", "single-cell", "spatial"):
        root = global_root / category
        if root.exists():
            candidates.extend(("global", path) for path in root.rglob("SKILL.md"))
    if project_root is not None and (project_root / "skills").exists():
        candidates.extend(("project", path) for path in (project_root / "skills").glob("*/SKILL.md"))

    skills: list[Skill] = []
    for scope, path in sorted(candidates, key=lambda item: str(item[1])):
        text = path.read_text(encoding="utf-8")
        name, description, _ = parse_frontmatter(path)
        skills.append(
            Skill(
                scope=scope,
                path=path,
                name=name,
                description=description,
                text=text,
                line_count=len(text.splitlines()),
            )
        )
    return skills


def content_words(skill: Skill) -> set[str]:
    return {
        word
        for word in WORD_RE.findall(skill.text.lower())
        if word not in STOP_WORDS
    }


def audit(
    skills: list[Skill],
    project_root: Path | None,
) -> tuple[list[str], list[str]]:
    errors: list[str] = []
    warnings: list[str] = []

    for skill in skills:
        try:
            name, description, keys = parse_frontmatter(skill.path)
        except ValueError as exc:
            errors.append(f"{skill.path}: {exc}")
            continue

        if keys != {"name", "description"}:
            errors.append(
                f"{skill.path}: frontmatter keys are {sorted(keys)}; expected name and description only"
            )
        if not name or not NAME_RE.fullmatch(name):
            errors.append(f"{skill.path}: invalid skill name {name!r}")
        if name != skill.path.parent.name:
            errors.append(
                f"{skill.path}: name {name!r} does not match folder {skill.path.parent.name!r}"
            )
        if len(description) < 40:
            warnings.append(f"{skill.path}: description may be too narrow ({len(description)} chars)")
        if skill.line_count > 500:
            errors.append(f"{skill.path}: SKILL.md has {skill.line_count} lines (>500)")

        agents = skill.path.parent / "agents" / "openai.yaml"
        if agents.exists():
            agent_text = agents.read_text(encoding="utf-8")
            first = next((line.strip() for line in agent_text.splitlines() if line.strip()), "")
            if first != "interface:":
                errors.append(f"{agents}: first key must be interface:")
            if ("$" + name) not in agent_text:
                errors.append(f"{agents}: default prompt does not mention $" + name)

        for reference in sorted(set(REFERENCE_RE.findall(skill.text))):
            target = skill.path.parent / reference.rstrip(".,:;)")
            if not target.exists():
                errors.append(f"{skill.path}: missing bundled reference {reference}")

        if skill.scope == "project" and project_root is not None:
            for script in sorted(set(SCRIPT_RE.findall(skill.text))):
                target = project_root / script.rstrip(".,:;)")
                if not target.exists():
                    errors.append(f"{skill.path}: missing project script {script}")

        for retired in RETIRED_ACTIVE_NAMES:
            if retired in skill.text:
                errors.append(f"{skill.path}: contains retired active skill name {retired}")

    by_scope: dict[str, list[Skill]] = {}
    for skill in skills:
        by_scope.setdefault(skill.scope, []).append(skill)

    for scope, scoped_skills in by_scope.items():
        names: dict[str, list[Path]] = {}
        for skill in scoped_skills:
            names.setdefault(skill.name, []).append(skill.path)
        for name, paths in names.items():
            if len(paths) > 1:
                errors.append(f"{scope}: duplicate skill name {name}: {paths}")

        word_sets = {skill.path: content_words(skill) for skill in scoped_skills}
        for left, right in itertools.combinations(scoped_skills, 2):
            union = word_sets[left.path] | word_sets[right.path]
            score = len(word_sets[left.path] & word_sets[right.path]) / max(1, len(union))
            if score >= 0.35:
                warnings.append(
                    f"{scope}: high content similarity {score:.3f}: {left.name} vs {right.name}"
                )

    taxonomy = next(
        (skill for skill in skills if skill.scope == "global" and skill.name == "skills-taxonomy"),
        None,
    )
    if taxonomy is None:
        errors.append("global: skills-taxonomy is missing")
    else:
        global_root = taxonomy.path.parents[2]
        actual_catalog = {
            str(skill.path.parent.relative_to(global_root))
            for skill in skills
            if skill.scope == "global"
            and skill.path.is_relative_to(global_root / "single-cell")
            or skill.scope == "global"
            and skill.path.is_relative_to(global_root / "spatial")
        }
        listed_catalog = set(
            re.findall(
                r"(?m)^\s*-\s+`((?:single-cell|spatial)/[a-z0-9][a-z0-9/-]+)`\s*$",
                taxonomy.text,
            )
        )
        for entry in sorted(listed_catalog):
            if not (global_root / entry / "SKILL.md").exists():
                errors.append(f"{taxonomy.path}: catalog entry does not exist: {entry}")
        missing_entries = actual_catalog - listed_catalog
        extra_entries = listed_catalog - actual_catalog
        if missing_entries:
            errors.append(
                f"{taxonomy.path}: active skills missing from catalog: {sorted(missing_entries)}"
            )
        if extra_entries:
            errors.append(
                f"{taxonomy.path}: non-active catalog entries: {sorted(extra_entries)}"
            )

    return errors, warnings


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--global-root",
        type=Path,
        default=Path.home() / ".codex" / "skills",
        help="Root containing _meta, single-cell, and spatial.",
    )
    parser.add_argument(
        "--project-root",
        type=Path,
        help="Optional repository root containing project-local skills/.",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Return nonzero for warnings as well as errors.",
    )
    args = parser.parse_args()

    try:
        skills = discover(args.global_root.expanduser(), args.project_root)
    except (OSError, ValueError) as exc:
        print(f"ERROR\t{exc}", file=sys.stderr)
        return 2

    errors, warnings = audit(skills, args.project_root)
    for message in errors:
        print(f"ERROR\t{message}")
    for message in warnings:
        print(f"WARN\t{message}")

    global_count = sum(skill.scope == "global" for skill in skills)
    project_count = sum(skill.scope == "project" for skill in skills)
    print(
        f"SUMMARY\tglobal={global_count}\tproject={project_count}"
        f"\terrors={len(errors)}\twarnings={len(warnings)}"
    )
    return 1 if errors or (args.strict and warnings) else 0


if __name__ == "__main__":
    raise SystemExit(main())
