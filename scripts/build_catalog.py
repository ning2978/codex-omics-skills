#!/usr/bin/env python3
"""Build the deterministic catalog for the published Codex omics skills."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


NAME_RE = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")


def parse_frontmatter(path: Path) -> tuple[str, str, set[str]]:
    lines = path.read_text(encoding="utf-8").splitlines()
    if not lines or lines[0].strip() != "---":
        raise ValueError(f"{path}: missing opening YAML delimiter")
    try:
        end = next(
            index
            for index, line in enumerate(lines[1:], start=1)
            if line.strip() == "---"
        )
    except StopIteration as exc:
        raise ValueError(f"{path}: missing closing YAML delimiter") from exc

    block = lines[1:end]
    keys: set[str] = set()
    values: dict[str, str] = {}
    index = 0
    while index < len(block):
        line = block[index]
        match = re.match(r"^([A-Za-z0-9_-]+):(?:\s*(.*))?$", line)
        if not match:
            index += 1
            continue
        key, raw_value = match.group(1), (match.group(2) or "").strip()
        keys.add(key)
        if raw_value in {">", ">-", "|", "|-"}:
            folded: list[str] = []
            index += 1
            while index < len(block) and (
                not block[index].strip()
                or block[index].startswith((" ", "\t"))
            ):
                if block[index].strip():
                    folded.append(block[index].strip())
                index += 1
            values[key] = " ".join(folded)
            continue
        values[key] = raw_value.strip("\"'")
        index += 1

    name = values.get("name", "")
    description = values.get("description", "")
    if not NAME_RE.fullmatch(name):
        raise ValueError(f"{path}: invalid skill name {name!r}")
    if not description:
        raise ValueError(f"{path}: missing description")
    return name, description, keys


def discover(repo_root: Path) -> list[dict[str, str]]:
    skills_root = repo_root / "skills"
    records: list[dict[str, str]] = []
    for skill_file in sorted(skills_root.rglob("SKILL.md")):
        name, description, _ = parse_frontmatter(skill_file)
        relative_dir = skill_file.parent.relative_to(repo_root)
        category_parts = relative_dir.relative_to("skills").parts[:-1]
        category = "/".join(category_parts) if category_parts else "cross-cutting"
        records.append(
            {
                "name": name,
                "category": category,
                "path": relative_dir.as_posix(),
                "description": description,
            }
        )
    return sorted(records, key=lambda item: item["name"])


def render(repo_root: Path) -> str:
    records = discover(repo_root)
    dependencies = json.loads((repo_root / "dependencies.json").read_text())
    names = {record["name"] for record in records}
    if set(dependencies) != names:
        raise ValueError("dependencies.json must declare exactly the catalogued skills")
    for record in records:
        entry = dependencies[record["name"]]
        if set(entry) != {"skills", "resources"} or set(entry["skills"]) - names:
            raise ValueError(f"Invalid dependencies for {record['name']}")
        for resource in entry["resources"]:
            path = (repo_root / resource).resolve()
            if not path.is_relative_to(repo_root / "skills" / "_shared") or not path.is_dir():
                raise ValueError(f"Invalid shared resource: {resource}")
        record["dependencies"] = sorted(set(entry["skills"]))
        record["shared_resources"] = sorted(set(entry["resources"]))
    payload = {
        "schema_version": 2,
        "skill_count": len(records),
        "skills": records,
    }
    return json.dumps(payload, ensure_ascii=False, indent=2) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Fail if catalog.json does not match the discovered skills.",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    catalog_path = repo_root / "catalog.json"
    expected = render(repo_root)
    if args.check:
        if not catalog_path.exists() or catalog_path.read_text(
            encoding="utf-8"
        ) != expected:
            print("ERROR: catalog.json is stale; run scripts/build_catalog.py")
            return 1
        print(f"OK: catalog.json contains {len(discover(repo_root))} skills")
        return 0

    catalog_path.write_text(expected, encoding="utf-8")
    print(f"Wrote {catalog_path} with {len(discover(repo_root))} skills")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
