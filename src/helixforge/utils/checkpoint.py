"""Stage-level checkpoint manifest."""

from __future__ import annotations

import hashlib
import json
from collections.abc import Iterable
from pathlib import Path
from typing import Any

from helixforge.utils.atomic import atomic_write
from helixforge.utils.logging import get_logger

_log = get_logger(__name__)

# Bytes sampled from each end of a file for the cheap content hash. Hashing the
# full file would be wasteful for a multi-GB maize loci.gff3; size + head + tail
# detects truncation and content changes cheaply (the common corruption modes).
_HASH_SAMPLE = 65536


def cheap_hash(path: str | Path) -> str:
    """A cheap, content-sensitive signature of ``path``.

    sha256 over the file size plus its first and last ``_HASH_SAMPLE`` bytes.
    Catches truncation (size changes) and edits at either end without reading a
    huge file end-to-end.
    """
    p = Path(path)
    size = p.stat().st_size
    h = hashlib.sha256()
    h.update(str(size).encode())
    with open(p, "rb") as fh:
        h.update(fh.read(_HASH_SAMPLE))
        if size > _HASH_SAMPLE:
            fh.seek(max(0, size - _HASH_SAMPLE))
            h.update(fh.read(_HASH_SAMPLE))
    return h.hexdigest()


class Checkpoint:
    """Per-``work_dir`` stage manifest. Inert unless ``enabled`` (resume) is set."""

    def __init__(self, path: str | Path, *, enabled: bool = False) -> None:
        self.path = Path(path)
        self.enabled = enabled
        self._data: dict[str, Any] = self._load() if enabled else {"stages": {}}

    def _load(self) -> dict[str, Any]:
        if self.path.exists():
            try:
                data: Any = json.loads(self.path.read_text())
            except (json.JSONDecodeError, OSError, ValueError):
                _log.warning(
                    "ignoring unreadable checkpoint %s: re-running", self.path
                )
                return {"stages": {}}
            if isinstance(data, dict) and isinstance(data.get("stages"), dict):
                return data
            _log.warning("ignoring malformed checkpoint %s, re-running", self.path)
        return {"stages": {}}

    def is_complete(self, stage: str, outputs: Iterable[str | Path]) -> bool:
        """True iff ``stage`` is recorded AND every recorded output validates.

        Validation: the recorded output set equals ``outputs``, and each file
        still exists, is non-empty, and matches its recorded cheap hash. Returns
        False when disabled (a cold run never skips anything).
        """
        if not self.enabled:
            return False
        rec = self._data["stages"].get(stage)
        if rec is None:
            return False
        recorded: dict[str, Any] = rec.get("outputs", {})
        if sorted(recorded) != sorted(str(Path(o)) for o in outputs):
            return False
        for path_str, info in recorded.items():
            p = Path(path_str)
            if not p.exists() or p.stat().st_size == 0:
                _log.info(
                    "checkpoint stage %r output missing/empty (%s), re-running",
                    stage,
                    p,
                )
                return False
            if info.get("hash") != cheap_hash(p):
                _log.info(
                    "checkpoint stage %r output changed (%s), re-running", stage, p
                )
                return False
        return True

    def mark(self, stage: str, outputs: Iterable[str | Path]) -> None:
        """Record ``stage`` as complete with the current state of ``outputs``.

        No-op when disabled (a cold run leaves no manifest behind, the default
        behavior is byte-identical to pre-Phase-22).
        """
        if not self.enabled:
            return
        # Materialise the iterable so we can traverse it twice safely.
        output_list = list(outputs)
        # Never record a false-complete: if the stage claims success but an
        # output is missing/empty, leave the manifest unmarked so a re-run
        # regenerates it (Phase 22 "what NOT to do" #2).
        missing = [
            str(Path(o))
            for o in output_list
            if not Path(o).exists() or Path(o).stat().st_size == 0
        ]
        if missing:
            _log.warning(
                "not checkpointing stage %r: outputs missing/empty: %s", stage, missing
            )
            return
        recorded: dict[str, Any] = {}
        for o in output_list:
            p = Path(o)
            recorded[str(p)] = {"size": p.stat().st_size, "hash": cheap_hash(p)}
        self._data["stages"][stage] = {"outputs": recorded}
        with atomic_write(self.path) as fh:
            fh.write(json.dumps(self._data, indent=2, sort_keys=True))
        _log.info("checkpoint: stage %r complete (%d outputs)", stage, len(recorded))
