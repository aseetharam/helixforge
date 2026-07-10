#!/usr/bin/env python3
"""Generate a per-command Markdown CLI reference from the live ``helixforge`` app.

This is the zero-third-party-dependency generator + drift gate described in
``AUDIT_REMEDIATION_AND_DOCS_BLUEPRINT.md`` §3.3. It introspects the click
command tree (``helixforge.cli:main``) so the reference is produced *from the
parser itself* and cannot drift from the code. ``click`` is the only import
(already a ``[cli]`` extra dependency).

The reference is rendered entirely from each command's own definition — its
**docstring** (`cmd.help`: description + the real output columns) and its
**epilog** (`cmd.epilog`: worked examples). The docstring/epilog are therefore
the *single source of truth*: ``--help`` and these Markdown pages render from the
same text, so updating a command updates both.

Layout (one file per command/subcommand, plus an index):

    docs/cli/README.md          # index: one line per command, linked
    docs/cli/reconcile.md
    docs/cli/parallel/plan.md   # group subcommands nest under the group dir
    ...

Usage::

    python scripts/gen_cli_reference.py            # write the docs/cli/ tree
    python scripts/gen_cli_reference.py --check     # CI/pre-commit: fail on drift
    python scripts/gen_cli_reference.py --stdout     # print the whole tree

In ``--check`` mode the script regenerates the tree to memory and compares it to
the committed ``docs/cli/`` tree (content **and** the set of files), exiting
non-zero with a hint on any difference — so adding a flag, renaming an option, or
editing a docstring without regenerating the docs fails the gate.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import click

REPO = Path(__file__).resolve().parents[1]
DEFAULT_OUT = REPO / "docs" / "cli"
INDEX_NAME = "README.md"

# click stores the ``\b`` no-rewrap marker as a literal backspace (\x08) inside
# ``help``/``epilog``. We strip it when rendering Markdown (Markdown controls its
# own wrapping via fenced blocks / hard newlines).
_BACKSPACE = "\b"


# ---------------------------------------------------------------------------
# Text helpers
# ---------------------------------------------------------------------------


def _squeeze_blanks(lines: list[str]) -> list[str]:
    """Collapse runs of blank lines (left behind by removed ``\\b`` markers)."""
    out: list[str] = []
    for line in lines:
        if not line.strip() and out and not out[-1].strip():
            continue
        out.append(line)
    while out and not out[0].strip():
        out.pop(0)
    while out and not out[-1].strip():
        out.pop()
    return out


def _description(cmd: click.Command) -> str:
    """The command description (full docstring), ``\\b`` markers removed."""
    text = (cmd.help or "").replace(_BACKSPACE, "")
    lines = [line.rstrip() for line in text.splitlines()]
    return "\n".join(_squeeze_blanks(lines))


def _short_help(cmd: click.Command) -> str:
    """One-line summary for the index table."""
    return cmd.get_short_help_str(limit=120).strip()


# ---------------------------------------------------------------------------
# Option table
# ---------------------------------------------------------------------------


def _opt_line(param: click.Parameter) -> str:
    """One Markdown table row for a click Option/Argument."""
    decls = ", ".join(f"`{o}`" for o in param.opts + param.secondary_opts)
    if isinstance(param, click.Argument):
        decls = f"`{param.name.upper()}`"
    meta = []
    if getattr(param, "required", False):
        meta.append("required")
    if isinstance(param, click.Option) and param.multiple:
        meta.append("repeatable")
    if isinstance(param.type, click.Choice):
        meta.append("choices: " + "/".join(param.type.choices))
    default = getattr(param, "default", None)
    # Suppress empty/sentinel defaults (e.g. the UNSET sentinel used for
    # repeatable options) so the table shows only meaningful defaults.
    default_repr = repr(default)
    is_sentinel = "Sentinel" in default_repr or "UNSET" in default_repr
    if (
        default is not None
        and default != ()
        and not is_sentinel
        and not getattr(param, "required", False)
    ):
        meta.append(f"default: `{default}`")
    help_text = (getattr(param, "help", "") or "").replace(_BACKSPACE, "").strip()
    detail = "; ".join(meta)
    if detail:
        help_text = f"{help_text} ({detail})" if help_text else f"({detail})"
    # Escape table-breaking pipes inside help text.
    help_text = help_text.replace("|", "\\|")
    return f"| {decls} | {help_text.strip()} |"


def _options_block(cmd: click.Command) -> list[str]:
    params = [p for p in cmd.params if not isinstance(p, click.Option) or not p.hidden]
    if not params:
        return []
    lines = ["## Options", "", "| Option | Description |", "|---|---|"]
    lines.extend(_opt_line(p) for p in params)
    lines.append("")
    return lines


# ---------------------------------------------------------------------------
# Examples block (from the epilog)
# ---------------------------------------------------------------------------


def _examples_block(cmd: click.Command) -> list[str]:
    epilog = getattr(cmd, "epilog", None)
    if not epilog:
        return []
    lines = epilog.replace(_BACKSPACE, "").splitlines()
    lines = _squeeze_blanks([line.rstrip() for line in lines])
    # Drop a leading "Examples:" / "Example:" header — the section heading
    # below already names the block.
    if lines and lines[0].strip().rstrip(":").lower() in ("example", "examples"):
        lines = _squeeze_blanks(lines[1:])
    if not lines:
        return []
    return ["## Examples", "", "```bash", *lines, "```", ""]


# ---------------------------------------------------------------------------
# Per-command page + index
# ---------------------------------------------------------------------------


def _command_page(parts: list[str], cmd: click.Command) -> str:
    full = "helixforge " + " ".join(parts)
    lines = [f"# `{full}`", ""]
    desc = _description(cmd)
    if desc:
        lines.append(desc)
        lines.append("")
    lines.extend(_options_block(cmd))
    lines.extend(_examples_block(cmd))
    return "\n".join(lines).rstrip() + "\n"


def _iter_leaf_commands(name: str, cmd: click.Command, parts: list[str]):
    """Yield ``(parts, leaf_command)`` for every non-group command in the tree.

    Hidden commands and groups (e.g. the dev-only ``benchmark``) are skipped, so
    the generated reference reflects exactly the user-facing command set.
    """
    if cmd.hidden:
        return
    if isinstance(cmd, click.Group):
        for sub_name in sorted(cmd.commands):
            yield from _iter_leaf_commands(
                sub_name, cmd.commands[sub_name], parts + [sub_name]
            )
    else:
        yield parts, cmd


def _rel_path(parts: list[str]) -> str:
    return "/".join(parts) + ".md"


def _index_page(root: click.Group, leaves: list[tuple[list[str], click.Command]]) -> str:
    lines = [
        "# HelixForge — Command-line reference",
        "",
        "> Auto-generated from the `helixforge` click app by "
        "`scripts/gen_cli_reference.py`. Do not edit by hand — run the script "
        "(or let the pre-commit/CI drift gate regenerate it). Each command's "
        "description, options, and examples come straight from its docstring and "
        "epilog, so this reference cannot drift from the code.",
        "",
    ]
    if root.help:
        lines.append(root.help.replace(_BACKSPACE, "").strip().splitlines()[0])
        lines.append("")
    lines.append("| Command | Description |")
    lines.append("|---|---|")
    for parts, cmd in sorted(leaves, key=lambda pc: pc[0]):
        full = "helixforge " + " ".join(parts)
        link = _rel_path(parts)
        desc = _short_help(cmd).replace("|", "\\|")
        lines.append(f"| [`{full}`]({link}) | {desc} |")
    lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def build_tree() -> dict[str, str]:
    """Return ``{relative_path: markdown}`` for the whole reference tree."""
    from helixforge.cli import main  # lazy so --help works without click extras

    leaves = list(_iter_leaf_commands("helixforge", main, []))
    tree: dict[str, str] = {INDEX_NAME: _index_page(main, leaves)}
    for parts, cmd in leaves:
        tree[_rel_path(parts)] = _command_page(parts, cmd)
    return tree


# ---------------------------------------------------------------------------
# Disk I/O + drift gate
# ---------------------------------------------------------------------------


def _read_existing(out_dir: Path) -> dict[str, str]:
    if not out_dir.exists():
        return {}
    return {
        str(p.relative_to(out_dir)): p.read_text()
        for p in out_dir.rglob("*.md")
    }


def _write_tree(out_dir: Path, tree: dict[str, str]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    wanted = set(tree)
    # Remove stale .md files (e.g. a renamed/removed command) before writing.
    for existing in out_dir.rglob("*.md"):
        if str(existing.relative_to(out_dir)) not in wanted:
            existing.unlink()
    for rel, content in tree.items():
        path = out_dir / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
    # Prune now-empty subdirectories.
    for sub in sorted(out_dir.rglob("*"), reverse=True):
        if sub.is_dir() and not any(sub.iterdir()):
            sub.rmdir()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(DEFAULT_OUT),
                    help="output directory for the per-command tree")
    ap.add_argument("--check", action="store_true",
                    help="fail if the committed tree is stale")
    ap.add_argument("--stdout", action="store_true",
                    help="print the whole tree to stdout instead of writing")
    args = ap.parse_args()

    tree = build_tree()
    out_dir = Path(args.out)

    if args.stdout:
        for rel in sorted(tree):
            sys.stdout.write(f"===== {rel} =====\n")
            sys.stdout.write(tree[rel])
            sys.stdout.write("\n")
        return 0

    if args.check:
        existing = _read_existing(out_dir)
        if existing != tree:
            only_disk = sorted(set(existing) - set(tree))
            only_new = sorted(set(tree) - set(existing))
            changed = sorted(
                rel for rel in set(tree) & set(existing)
                if tree[rel] != existing[rel]
            )
            sys.stderr.write(
                f"CLI reference is stale: {out_dir} differs from the live CLI.\n"
            )
            if only_new:
                sys.stderr.write(f"  missing files: {', '.join(only_new)}\n")
            if only_disk:
                sys.stderr.write(f"  extra files:   {', '.join(only_disk)}\n")
            if changed:
                sys.stderr.write(f"  changed files: {', '.join(changed)}\n")
            sys.stderr.write(
                f"Regenerate it with:  python {Path(__file__).name}\n"
            )
            return 1
        return 0

    _write_tree(out_dir, tree)
    sys.stderr.write(f"wrote {len(tree)} files under {out_dir}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
