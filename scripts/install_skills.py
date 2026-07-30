#!/usr/bin/env python3
"""Install selected catalogued skills without deleting unrelated skills."""

from __future__ import annotations

import argparse
import json
import os
import shutil
from pathlib import Path


def default_target() -> Path:
    codex_root = os.environ.get("CODEX_HOME")
    if codex_root:
        return Path(codex_root).expanduser() / "skills"
    return Path.home() / ".codex" / "skills"


def load_catalog(repo_root: Path) -> list[dict[str, str]]:
    payload = json.loads((repo_root / "catalog.json").read_text(encoding="utf-8"))
    records = payload.get("skills")
    if not isinstance(records, list):
        raise ValueError("catalog.json does not contain a skills list")
    return records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--target",
        type=Path,
        default=default_target(),
        help="Destination skills root.",
    )
    parser.add_argument(
        "--name",
        action="append",
        default=[],
        help="Install one named skill; repeat to select multiple skills.",
    )
    parser.add_argument(
        "--category",
        help="Install skills whose category starts with this value.",
    )
    parser.add_argument("--list", action="store_true", help="List available skills.")
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned copies without changing the destination.",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Overlay existing skill directories; never deletes unrelated files.",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    records = load_catalog(repo_root)
    if args.list:
        for record in records:
            print(f"{record['name']}\t{record['category']}\t{record['path']}")
        return 0

    selected = records
    if args.name:
        requested = set(args.name)
        known = {record["name"] for record in records}
        unknown = sorted(requested - known)
        if unknown:
            print(f"ERROR: unknown skill name(s): {', '.join(unknown)}")
            return 2
        selected = [record for record in selected if record["name"] in requested]
    if args.category:
        selected = [
            record
            for record in selected
            if record["category"] == args.category
            or record["category"].startswith(args.category.rstrip("/") + "/")
        ]
    if not selected:
        print("ERROR: selection matched no skills")
        return 2

    target_root = args.target.expanduser().resolve()
    failures = 0
    copied = 0
    for record in selected:
        source = (repo_root / record["path"]).resolve()
        relative = source.relative_to((repo_root / "skills").resolve())
        destination = target_root / relative
        if destination.exists() and not args.force:
            print(f"SKIP: {record['name']} already exists at {destination}")
            failures += 1
            continue
        print(f"{'PLAN' if args.dry_run else 'COPY'}: {source} -> {destination}")
        if not args.dry_run:
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copytree(source, destination, dirs_exist_ok=args.force)
        copied += 1

    print(
        f"SUMMARY: selected={len(selected)} copied_or_planned={copied} "
        f"skipped={failures} target={target_root}"
    )
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
