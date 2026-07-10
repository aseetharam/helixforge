"""Synthetic fixtures for Phase 3 locus + classification tests.

Concrete literal coordinates only. No real biological data.
"""

import h5py
import numpy as np
import pytest

# Helixer GFF3 for load_helixer_loci:
#   gene1 chr1 + 101-200   } overlap same strand -> merge into (100,300)
#   gene2 chr1 + 151-300   }
#   gene3 chr1 + 401-500   separate
#   gene4 chr1 - 121-250   opposite strand (antisense) -> NOT merged
#   gene5 chr2 + 101-200   second scaffold
_LOAD_GFF = """\
##gff-version 3
chr1\tHelixer\tgene\t101\t200\t.\t+\t.\tID=gene1
chr1\tHelixer\tmRNA\t101\t200\t.\t+\t.\tID=gene1.m;Parent=gene1
chr1\tHelixer\texon\t101\t200\t.\t+\t.\tID=gene1.e1;Parent=gene1.m
chr1\tHelixer\tgene\t151\t300\t.\t+\t.\tID=gene2
chr1\tHelixer\tmRNA\t151\t300\t.\t+\t.\tID=gene2.m;Parent=gene2
chr1\tHelixer\texon\t151\t300\t.\t+\t.\tID=gene2.e1;Parent=gene2.m
chr1\tHelixer\tgene\t401\t500\t.\t+\t.\tID=gene3
chr1\tHelixer\tmRNA\t401\t500\t.\t+\t.\tID=gene3.m;Parent=gene3
chr1\tHelixer\texon\t401\t500\t.\t+\t.\tID=gene3.e1;Parent=gene3.m
chr1\tHelixer\tgene\t121\t250\t.\t-\t.\tID=gene4
chr1\tHelixer\tmRNA\t121\t250\t.\t-\t.\tID=gene4.m;Parent=gene4
chr1\tHelixer\texon\t121\t250\t.\t-\t.\tID=gene4.e1;Parent=gene4.m
chr2\tHelixer\tgene\t101\t200\t.\t+\t.\tID=gene5
chr2\tHelixer\tmRNA\t101\t200\t.\t+\t.\tID=gene5.m;Parent=gene5
chr2\tHelixer\texon\t101\t200\t.\t+\t.\tID=gene5.e1;Parent=gene5.m
"""


@pytest.fixture
def load_gff_path(tmp_path):
    p = tmp_path / "helixer_load.gff3"
    p.write_text(_LOAD_GFF)
    return str(p)


# HDF5 with uniform high-CDS predictions -> region_confidence ~= 0.90 everywhere.
_CDS_VEC = [0.05, 0.05, 0.90, 0.00]
_LOAD_SCAFFOLDS = {"chr1": 512, "chr2": 256}
_CHUNK = 128


@pytest.fixture
def load_h5_path(tmp_path):
    p = tmp_path / "helixer_load.h5"
    seqids, start_ends, rows = [], [], []
    for sid, length in _LOAD_SCAFFOLDS.items():
        for low in range(0, length, _CHUNK):
            high = min(low + _CHUNK, length)
            row = np.zeros((_CHUNK, 4), dtype=np.float32)
            for j in range(high - low):
                row[j] = _CDS_VEC
            rows.append(row)
            seqids.append(sid.encode())
            start_ends.append([low, high])
    with h5py.File(p, "w") as f:
        f.create_dataset("predictions", data=np.stack(rows))
        f.create_dataset("seqids", data=np.array(seqids))
        f.create_dataset("start_ends", data=np.array(start_ends, dtype=np.int64))
    return str(p)
