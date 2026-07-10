"""Interval index over NCLS (Nested Containment List)."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Any

import numpy as np
from ncls import NCLS


def bounds(item: Any) -> tuple[int, int]:
    """Return ``(start, end)`` from an ``.start``/``.end`` object or a 2-tuple.

    Canonical interval-accessor (de-duplicated from
    ``reconcile/{mikado_integrate,as_events}.py``). Objects
    carrying ``.start``/``.end`` (Exon, CDSSegment, Interval) are read by
    attribute; anything else is treated as an indexable ``(start, end, ...)``.
    """
    if hasattr(item, "start"):
        return item.start, item.end
    return item[0], item[1]


class IntervalIndex:
    """A queryable index of intervals carrying arbitrary tuple payloads.

    Each added interval is a tuple ``(start, end, *data)``; the integer id
    returned by :meth:`query` is its insertion index. :meth:`query_with_data`
    returns the full stored tuples for overlapping intervals.
    """

    def __init__(self) -> None:
        self._intervals: list[tuple[Any, ...]] = []
        self._ncls: NCLS | None = None

    def add_intervals(self, intervals: Iterable[Sequence[Any]]) -> None:
        """Add a list of ``(start, end, *data)`` tuples and (re)build the index."""
        for iv in intervals:
            iv = tuple(iv)
            if len(iv) < 2:
                raise ValueError(
                    f"interval must have at least (start, end), got {iv!r}"
                )
            start, end = int(iv[0]), int(iv[1])
            if end <= start:
                raise ValueError(f"require start < end, got {start}-{end}")
            self._intervals.append((start, end) + iv[2:])
        self._rebuild()

    def _rebuild(self) -> None:
        if not self._intervals:
            self._ncls = None
            return
        starts = np.array([iv[0] for iv in self._intervals], dtype=np.int64)
        ends = np.array([iv[1] for iv in self._intervals], dtype=np.int64)
        ids = np.arange(len(self._intervals), dtype=np.int64)
        self._ncls = NCLS(starts, ends, ids)

    def query(self, start: int, end: int) -> list[int]:
        """Return sorted ids of intervals overlapping ``[start, end)``."""
        if self._ncls is None:
            return []
        return sorted(i for _, _, i in self._ncls.find_overlap(int(start), int(end)))

    def query_with_data(self, start: int, end: int) -> list[tuple[Any, ...]]:
        """Return the stored ``(start, end, *data)`` tuples overlapping ``[start, end)``."""
        return [self._intervals[i] for i in self.query(start, end)]

    def __len__(self) -> int:
        return len(self._intervals)

    @property
    def is_empty(self) -> bool:
        return len(self._intervals) == 0
