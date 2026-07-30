#!/usr/bin/env python3
"""Validate structure, portability, metadata, and catalog integrity."""

from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True

import build_catalog


MARKDOWN_LINK_RE = re.compile(r"\[[^\]]+\]\(([^)]+)\)")
RESOURCE_RE = re.compile(
    r"\b((?:references|scripts)/[A-Za-z0-9_./-]+\.(?:md|py|sh|R|json|ya?ml))"
)
SECURITY_PATTERNS = {
    "Linux home absolute path": re.compile(r"/home/[A-Za-z0-9._-]+/"),
    "macOS home absolute path": re.compile(r"/Users/[A-Za-z0-9._-]+/"),
    "local account identifier": re.compile(
        r"\b" + "t100" + "737" + r"\b", re.IGNORECASE
    ),
    "private project name": re.compile(
        "prostate_tumor" + "_sc_ST_CODEX_downstream", re.IGNORECASE
    ),
    "private sample identifier": re.compile(
        r"\b(?:" + "0066" + r"301|" + "0085" + r"005)\b"
    ),
    "private key": re.compile(r"BEGIN [A-Z ]*PRIVATE KEY"),
    "GitHub token": re.compile(r"\b(?:ghp|github_pat)_[A-Za-z0-9_]{20,}\b"),
    "OpenAI-style token": re.compile(r"\bsk-[A-Za-z0-9_-]{20,}\b"),
    "AWS access key": re.compile(r"\bAKIA[0-9A-Z]{16}\b"),
}
FORBIDDEN_SUFFIXES = {
    ".h5ad",
    ".h5",
    ".loom",
    ".csv",
    ".tsv",
    ".parquet",
    ".feather",
    ".rds",
    ".rdata",
    ".xlsx",
    ".xls",
    ".log",
    ".pid",
    ".pyc",
}
TEXT_SUFFIXES = {".md", ".py", ".R", ".sh", ".json", ".yaml", ".yml"}


def validate(repo_root: Path) -> list[str]:
    errors: list[str] = []
    skill_files = sorted((repo_root / "skills").rglob("SKILL.md"))
    names: dict[str, Path] = {}

    for skill_file in skill_files:
        try:
            name, description, keys = build_catalog.parse_frontmatter(skill_file)
        except ValueError as exc:
            errors.append(str(exc))
            continue
        if keys != {"name", "description"}:
            errors.append(
                f"{skill_file}: frontmatter keys are {sorted(keys)}; "
                "expected name and description only"
            )
        if skill_file.parent.name != name:
            errors.append(
                f"{skill_file}: name {name!r} does not match folder "
                f"{skill_file.parent.name!r}"
            )
        if len(name) > 63:
            errors.append(f"{skill_file}: name has {len(name)} characters (>63)")
        if len(description) < 40:
            errors.append(f"{skill_file}: description is too narrow")
        line_count = len(skill_file.read_text(encoding="utf-8").splitlines())
        if line_count > 500:
            errors.append(f"{skill_file}: {line_count} lines (>500)")
        if name in names:
            errors.append(f"duplicate skill name {name}: {names[name]} and {skill_file}")
        names[name] = skill_file

        agents_file = skill_file.parent / "agents" / "openai.yaml"
        if not agents_file.exists():
            errors.append(f"{agents_file}: missing")
        else:
            agents_text = agents_file.read_text(encoding="utf-8")
            first_key = next(
                (line.strip() for line in agents_text.splitlines() if line.strip()),
                "",
            )
            if first_key != "interface:":
                errors.append(f"{agents_file}: first key must be interface:")
            if f"${name}" not in agents_text:
                errors.append(
                    f"{agents_file}: default prompt does not mention ${name}"
                )
            short_match = re.search(
                r'(?m)^\s*short_description:\s*["\'](.*?)["\']\s*$',
                agents_text,
            )
            if not short_match:
                errors.append(f"{agents_file}: quoted short_description is missing")
            elif not 25 <= len(short_match.group(1)) <= 64:
                errors.append(
                    f"{agents_file}: short_description length is "
                    f"{len(short_match.group(1))}, expected 25-64"
                )

        skill_text = skill_file.read_text(encoding="utf-8")
        for resource in sorted(set(RESOURCE_RE.findall(skill_text))):
            target = skill_file.parent / resource.rstrip(".,:;)")
            if not target.exists():
                errors.append(f"{skill_file}: missing resource {resource}")

    for path in sorted(repo_root.rglob("*")):
        if ".git" in path.parts:
            continue
        if path.is_symlink():
            errors.append(f"{path}: symlinks are not allowed in the portable kit")
            continue
        if not path.is_file():
            continue
        relative = path.relative_to(repo_root)
        if path.suffix.lower() in {suffix.lower() for suffix in FORBIDDEN_SUFFIXES}:
            errors.append(f"{relative}: forbidden data/cache suffix")
        if path.stat().st_size > 1_000_000:
            errors.append(f"{relative}: file is larger than 1 MB")
        if path.suffix not in TEXT_SUFFIXES and path.name not in {
            ".gitignore",
            "AGENTS.md",
        }:
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            errors.append(f"{relative}: expected UTF-8 text")
            continue
        for label, pattern in SECURITY_PATTERNS.items():
            if pattern.search(text):
                errors.append(f"{relative}: contains {label}")

        if path.suffix == ".md":
            for target_text in MARKDOWN_LINK_RE.findall(text):
                target_text = target_text.strip().strip("<>")
                if (
                    not target_text
                    or target_text.startswith(("#", "http://", "https://", "mailto:"))
                ):
                    continue
                target_path = target_text.split("#", 1)[0]
                if target_path and not (path.parent / target_path).resolve().exists():
                    errors.append(f"{relative}: broken Markdown link {target_text}")

        if path.suffix == ".py":
            try:
                compile(text, str(relative), "exec")
            except SyntaxError as exc:
                errors.append(f"{relative}: Python syntax error: {exc}")

    expected_catalog = build_catalog.render(repo_root)
    catalog_path = repo_root / "catalog.json"
    if not catalog_path.exists():
        errors.append("catalog.json: missing")
    elif catalog_path.read_text(encoding="utf-8") != expected_catalog:
        errors.append("catalog.json: stale; run scripts/build_catalog.py")

    taxonomy_audit = (
        repo_root
        / "skills"
        / "_meta"
        / "skills-taxonomy"
        / "scripts"
        / "audit_bio_skills.py"
    )
    if taxonomy_audit.exists():
        result = subprocess.run(
            [
                sys.executable,
                str(taxonomy_audit),
                "--global-root",
                str(repo_root / "skills"),
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode:
            detail = (result.stdout + result.stderr).strip()
            errors.append(f"taxonomy audit failed:\n{detail}")

    if len(skill_files) != 34:
        errors.append(f"expected 34 skills, discovered {len(skill_files)}")
    return errors


def main() -> int:
    repo_root = Path(__file__).resolve().parents[1]
    errors = validate(repo_root)
    for error in errors:
        print(f"ERROR: {error}")
    if errors:
        print(f"FAILED: {len(errors)} error(s)")
        return 1
    print("OK: 34 skills passed structure, metadata, link, syntax, and safety checks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
