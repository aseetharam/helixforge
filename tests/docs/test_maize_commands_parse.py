"""Every ``helixforge ...`` command in ``docs/maize/*.md`` parses under the live
click group.

The maize handoff documents the *only* path that scales to a large genome, so a
documented flag that the CLI does not actually have would silently mislead an
external user mid-run. This test extracts every fenced ``helixforge ...`` command
from the maize docs and replays it against the live ``helixforge.cli:main`` click
group with ``--help`` appended: ``--help`` is eager, so it short-circuits before
required-option / file-existence checks (which need real inputs we do not have),
but an **unknown option or subcommand still fails at tokenisation** (exit 2).
A documented flag can therefore never be one the CLI lacks.

No real tools or files are needed — parsing only. (Commands that would need real
inputs to *execute* are out of scope here; this asserts they *parse*.)
"""

from __future__ import annotations

import re
import shlex
from pathlib import Path

import pytest
from click.testing import CliRunner

from helixforge.cli import main

MAIZE_DOCS = Path(__file__).resolve().parents[2] / "docs" / "maize"

# Fenced code blocks: ```...\n<body>\n```
_FENCE = re.compile(r"```[^\n]*\n(.*?)```", re.DOTALL)


def _commands_from(text: str) -> list[str]:
    """Return the ``helixforge ...`` commands inside fenced blocks of ``text``.

    Joins backslash line-continuations and ignores non-helixforge lines (shell
    setup like ``samtools``/``awk``/``sbatch`` is not our CLI to validate).
    """
    commands: list[str] = []
    for block in _FENCE.findall(text):
        # stitch backslash-continued lines into one logical command
        logical = block.replace("\\\n", " ")
        for line in logical.splitlines():
            line = line.strip()
            if line.startswith("helixforge "):
                commands.append(line)
    return commands


def _all_maize_commands() -> list[str]:
    commands: list[str] = []
    for md in sorted(MAIZE_DOCS.glob("*.md")):
        commands.extend(_commands_from(md.read_text()))
    return commands


ALL_COMMANDS = _all_maize_commands()


def test_maize_docs_dir_present():
    """The maize handoff directory and its README exist."""
    assert MAIZE_DOCS.is_dir(), f"missing {MAIZE_DOCS}"
    assert (MAIZE_DOCS / "README.md").exists()


def test_at_least_one_command_extracted():
    """The extractor finds helixforge commands to validate (guards a silent no-op)."""
    assert len(ALL_COMMANDS) >= 5, (
        f"expected to extract >=5 helixforge commands, got {len(ALL_COMMANDS)}"
    )


def test_every_doc_file_present():
    """README + 00..04 are all present."""
    expected = {
        "README.md",
        "00_inputs_checklist.md",
        "01_quickstart.md",
        "02_end_to_end_chunked.md",
        "03_troubleshooting.md",
        "04_caveats_known_limits.md",
    }
    present = {p.name for p in MAIZE_DOCS.glob("*.md")}
    assert expected <= present, f"missing maize docs: {expected - present}"


@pytest.mark.parametrize("command", ALL_COMMANDS, ids=lambda c: c[:60])
def test_documented_command_parses(command: str):
    """Each documented command resolves to a real subcommand with only real flags.

    ``--help`` is appended so the parse stops before required/file checks; an
    unknown option or subcommand still produces a nonzero exit.
    """
    tokens = shlex.split(command)
    assert tokens[0] == "helixforge"
    args = tokens[1:] + ["--help"]
    result = CliRunner().invoke(main, args)
    assert result.exit_code == 0, (
        f"documented command failed to parse (exit {result.exit_code}):\n"
        f"  {command}\n"
        f"--- output ---\n{result.output}"
    )


def test_flags_are_known_to_their_subcommand():
    """Sanity backstop: a fabricated flag on a real subcommand is rejected.

    Confirms the --help-append strategy genuinely catches unknown options (so a
    green parametrized run above is meaningful, not vacuous).
    """
    result = CliRunner().invoke(
        main, ["parallel", "plan", "--genome", "x.fa", "--not-a-real-flag", "--help"]
    )
    assert result.exit_code != 0
