#!/usr/bin/env python3
import argparse
from collections import defaultdict
from typing import Any, cast
import itertools as it
import more_itertools as itx
import sys
from snakebids import (
    bidsapp,
    plugins,
    generate_inputs,
    bids,
    set_bids_spec,
)
from snakebids.utils.containers import MultiSelectDict
import json

from bidsarray.component import parse_component
from bidsarray.consensus import consensus_table, excluded_table

set_bids_spec("v0_11_0")


def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir."""
    return app.build_parser().parser


@bidsapp.hookimpl
def add_cli_arguments(parser: argparse.ArgumentParser):
    parser.add_argument(
        "--skipped-wildcards",
        default=False,
        action="store_true",
        help="Print wildcard card combinations not used due to missing components",
    )


@bidsapp.hookimpl
def get_argv(argv: list[str], config: dict[str, Any]):
    main_args, *components = itx.split_at(argv, lambda x: x == ":::")
    config["component_args"] = components
    return main_args


@bidsapp.hookimpl
def finalize_config(config: dict[str, Any]):
    components = config["component_args"]
    app.config["pybids_inputs"] = defaultdict[str | int, Any](dict)
    app.config["outputs"] = defaultdict[str | int, Any](dict)
    for ix, comp in enumerate(components):
        parse_component(ix, comp, app.config)
    pybids_inputs = app.config["pybids_inputs"]

    are_labelled = [
        isinstance(label, str)
        for label in it.chain(pybids_inputs, app.config["outputs"])
    ]
    labelled = all(are_labelled)
    if not labelled and any(are_labelled):
        raise ValueError(
            "Either all components must have a label, or none may be labelled"
        )
    if labelled:
        labels = [""] * len(pybids_inputs)
        for label, comp in pybids_inputs.items():
            labels[comp["order"]] = str(label)
        for label, comp in app.config["outputs"].items():
            labels[comp["order"]] = str(label)
        config["labels"] = labels


@bidsapp.hookimpl
def run(config: dict[str, Any]):
    pybids_inputs = config["pybids_inputs"]
    derivatives = app.config["derivatives"]
    derivatives = (
        False if derivatives is None else True if derivatives == [] else derivatives
    )
    inputs = generate_inputs(
        bids_dir=app.config["bids_dir"],
        pybids_inputs=cast("dict[str, Any]", app.config["pybids_inputs"]),
        pybidsdb_dir=app.config.get("pybidsdb_dir"),
        pybidsdb_reset=app.config.get("pybidsdb_reset", False),
        derivatives=derivatives,
        participant_label=app.config.get("participant_label"),
        exclude_participant_label=app.config.get("exclude_participant_label", None),
    )
    if len(inputs) != len(app.config["pybids_inputs"]):
        exit(1)
    if config["skipped_wildcards"]:
        table, comps = excluded_table(inputs)
        keys = list(table.keys())
        print(
            json.dumps(
                [
                    {
                        "wildcards": (wcards := dict(zip(keys, row))),
                        "paths": [comp.path.format(**wcards) for comp in comps[i]],
                    }
                    for i, row in enumerate(zip(*table.values()))
                ]
            )
        )
        return
    table = consensus_table(inputs)
    grouped_entities = tuple(
        set(
            it.chain.from_iterable(
                pybids_inputs[label]["groupby"] for label in pybids_inputs
            )
        )
    )

    ncomponents = len(pybids_inputs) + len(app.config["outputs"])
    if "labels" in config:
        print("\t".join(config["labels"]))
    for row_ in zip(*table.values()):
        row = MultiSelectDict(zip(table.keys(), row_))
        out = [""] * ncomponents
        for label, comp in inputs.items():
            out[app.config["pybids_inputs"][label]["order"]] = " ".join(
                comp.filter(**row[tuple(pybids_inputs[label]["groupby"])]).expand()
            )
        for comp in app.config["outputs"].values():
            out[comp["order"]] = bids(
                app.config["output_dir"],
                **(row[grouped_entities] | comp["entities"]),
            )
        try:
            print("\t".join(out))
            pass
        except BrokenPipeError:
            sys.exit(0)


app = bidsapp.app(
    [
        plugins.BidsArgs(analysis_level=False),
        plugins.Pybidsdb(),
        plugins.Version(distribution="bidsarray"),
        sys.modules[__name__],
    ],
)
app.parser.add_argument(
    "--derivatives",
    nargs="*",
)


if __name__ == "__main__":
    app.run()
