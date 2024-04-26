from __future__ import annotations

import itertools as it
from collections import defaultdict
from typing import Sequence
from snakebids import BidsDataset
from snakebids.utils.containers import MultiSelectDict


def consensus_table(dataset: BidsDataset):
    all_entities: dict[str, set[str]] = defaultdict(set)
    for comp in dataset.values():
        for entity, vals in comp.entities.items():
            all_entities[entity] |= set(vals)
    table: dict[str, Sequence[str]] = MultiSelectDict(
        zip(all_entities.keys(), zip(*it.product(*all_entities.values())))
    )
    keep_indices: set[int] = set()
    for comp in dataset.values():
        keys = tuple(comp.zip_lists.keys())
        rows = set(zip(*comp.zip_lists.values()))
        keep_indices |= {
            i for i, row in enumerate(zip(*table[keys].values())) if row in rows
        }
    return MultiSelectDict(
        {key: [val[i] for i in keep_indices] for key, val in table.items()}
    )
