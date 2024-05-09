from __future__ import annotations

import itertools as it
from collections import defaultdict
from typing import Literal, overload
from snakebids import BidsComponent, BidsDataset
from snakebids.utils.containers import MultiSelectDict


@overload
def _filter_table(
    *, dataset: BidsDataset, inverse: Literal[True]
) -> tuple[MultiSelectDict[str, list[str]], list[list[BidsComponent]]]: ...


@overload
def _filter_table(
    *, dataset: BidsDataset, inverse: Literal[False]
) -> MultiSelectDict[str, list[str]]: ...


def _filter_table(*, dataset: BidsDataset, inverse: bool):
    all_entities: dict[str, set[str]] = defaultdict(set)
    for comp in dataset.values():
        for entity, vals in comp.entities.items():
            all_entities[entity] |= set(vals)
    table = MultiSelectDict(
        zip(all_entities.keys(), zip(*it.product(*all_entities.values())))
    )
    keep_indices: set[int] | None = None
    missing_from: dict[int, list[BidsComponent]] = defaultdict(list)
    for comp in dataset.values():
        keys = tuple(comp.zip_lists.keys())
        rows = set(zip(*comp.zip_lists.values()))
        new_indices = {
            i
            for i, row in enumerate(zip(*table[keys].values()))
            if (row in rows) ^ inverse
        }
        if keep_indices is None:
            keep_indices = new_indices
        elif inverse:
            keep_indices |= new_indices
        else:
            keep_indices &= new_indices
        if inverse:
            for i in new_indices:
                missing_from[i].append(comp)
    if keep_indices is None:
        keep_indices = set()

    if inverse:
        new_indices = {
            i for i, comps in missing_from.items() if len(comps) < len(dataset)
        }
        keep_indices &= new_indices

    filt_table = MultiSelectDict(
        {key: [val[i] for i in keep_indices] for key, val in table.items()}
    )
    if inverse:
        missing = [missing_from[i] for i in keep_indices]
        return filt_table, missing
    return filt_table


def consensus_table(dataset: BidsDataset):
    return _filter_table(dataset=dataset, inverse=False)


def excluded_table(dataset: BidsDataset):
    return _filter_table(dataset=dataset, inverse=True)
