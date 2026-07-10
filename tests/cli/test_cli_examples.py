"""Keep the CLI examples honest as the CLI evolves.

Every command's worked examples live in its click ``epilog`` (the single source
of truth that both ``--help`` and the generated ``docs/cli/`` pages render from).
This module extracts each example command and asserts it parses under the *live*
click CLI: the subcommand path resolves and every option token is a real option
of the resolved command. A stale example that references a removed or renamed
flag therefore fails the test.

Parsing is purely structural — no command is executed and no input file is
required — so nothing here needs the ``integration`` marker.
"""

from __future__ import annotations

import shlex

import click
import pytest

from helixforge import cli


# ---------------------------------------------------------------------------
# Walk the live click tree
# ---------------------------------------------------------------------------


def _leaf_commands(cmd, parts):
    """Yield ``(parts, leaf_command)`` for every non-group command."""
    if isinstance(cmd, click.Group):
        for sub_name in sorted(cmd.commands):
            yield from _leaf_commands(cmd.commands[sub_name], parts + [sub_name])
    else:
        yield parts, cmd


LEAVES = list(_leaf_commands(cli.main, []))


def _example_lines(epilog):
    """Yield each ``helixforge …`` example command line from an epilog.

    Joins shell line-continuations (``\\`` at end of line), then keeps only the
    lines that invoke ``helixforge`` (skipping comments and the occasional
    ``parallel``/``hs`` executor line). Nested ``helixforge`` inside a quoted
    ``--command '…'`` value stays a single shell token and is *not* treated as a
    separate command.
    """
    if not epilog:
        return
    text = epilog.replace("\b", "").replace("\\\n", " ")
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("$ "):  # tolerate a v1-style shell prompt prefix
            line = line[2:].strip()
        if line.startswith("helixforge "):
            yield line


def _resolve(tokens):
    """Walk subcommand tokens to the leaf command; return (command, remaining)."""
    cmd = cli.main
    i = 0
    while i < len(tokens) and isinstance(cmd, click.Group) and not tokens[i].startswith("-"):
        if tokens[i] in cmd.commands:
            cmd = cmd.commands[tokens[i]]
            i += 1
        else:
            break
    return cmd, tokens[i:]


def _known_opts(cmd):
    known = set()
    for param in cmd.params:
        known.update(param.opts)
        known.update(param.secondary_opts)
    return known


# Build the (command-path, example) cases at import time so a bad example shows
# up as its own failing test id.
_CASES = []
for _parts, _cmd in LEAVES:
    for _ex in _example_lines(getattr(_cmd, "epilog", None)):
        _CASES.append(
            pytest.param(" ".join(_parts), _ex, id=f"{'_'.join(_parts)}-{len(_CASES)}")
        )


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("cmd_path,example", _CASES)
def test_example_parses_under_live_cli(cmd_path, example):
    tokens = shlex.split(example)
    assert tokens and tokens[0] == "helixforge", example

    cmd, rest = _resolve(tokens[1:])
    assert not isinstance(cmd, click.Group), (
        f"example did not resolve to a command (stopped at a group): {example}"
    )

    known = _known_opts(cmd)
    for tok in rest:
        if tok.startswith("-") and tok != "-":
            name = tok.split("=", 1)[0]
            assert name in known, (
                f"unknown option {name!r} in example for `helixforge {cmd_path}`:\n"
                f"  {example}"
            )


@pytest.mark.parametrize(
    "cmd_path,cmd",
    [pytest.param(" ".join(p), c, id="_".join(p)) for p, c in LEAVES],
)
def test_every_command_has_at_least_one_example(cmd_path, cmd):
    examples = list(_example_lines(getattr(cmd, "epilog", None)))
    assert examples, f"command `helixforge {cmd_path}` has no worked examples in its epilog"


def test_cases_cover_all_commands():
    """Sanity: we actually collected examples (guards against a parsing regression)."""
    assert len(_CASES) >= len(LEAVES)
