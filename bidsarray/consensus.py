from __future__ import annotations

import itertools as it
from collections import defaultdict
from typing import Sequence
from snakebids import BidsDataset
from snakebids.utils.containers import MultiSelectDict


def _filter_table(*, dataset: BidsDataset, inverse: bool):
    all_entities: dict[str, set[str]] = defaultdict(set)
    for comp in dataset.values():
        for entity, vals in comp.entities.items():
            all_entities[entity] |= set(vals)
    table = MultiSelectDict(
        zip(all_entities.keys(), zip(*it.product(*all_entities.values())))
    )
    keep_indices: set[int] | None = set() if inverse else None
    missing_from: dict[str, list[BidsDataset]] = defaultdict(list)
    for comp in dataset.values():
        keys = tuple(comp.zip_lists.keys())
        rows = set(zip(*comp.zip_lists.values()))
        new_indices = {
            i
            for i, row in enumerate(zip(*table[keys].values()))
            if (row in rows) ^ inverse
        }
        if inverse:
            keep_indices |= new_indices
            for i in new_indices:
                missing_from[i].append(comp)
        elif keep_indices is None:
            keep_indices = new_indices
        else:
            keep_indices &= new_indices
    if inverse:
        new_indices = {
            i for i, comps in missing_from.items() if len(comps) < len(dataset)
        }
        keep_indices &= new_indices
        missing = [missing_from[i] for i in keep_indices]

    filt_table = MultiSelectDict(
        {key: [val[i] for i in keep_indices] for key, val in table.items()}
    )
    if inverse:
        return filt_table, missing
    return filt_table


def consensus_table(dataset: BidsDataset):
    return _filter_table(dataset=dataset, inverse=False)


def excluded_table(dataset: BidsDataset):
    return _filter_table(dataset=dataset, inverse=True)
