#!/usr/bin/env python3
import argparse
from collections import defaultdict
from types import EllipsisType
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
from snakebids.plugins.component_edit import FilterParse
from snakebids.utils.utils import text_fold

from bidsarray.consensus import consensus_table

set_bids_spec("v0_11_0")


def make_component_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        nargs="?",
        metavar="LABEL",
        const=...,
        help="Classify the component as an input. Optionally provide a name",
    )
    parser.add_argument(
        "--output",
        nargs="?",
        metavar="LABEL",
        const=...,
        help="Classify the component as an output. Optionally provide a name",
    )
    parser.add_argument(
        "--filter",
        action=FilterParse,
        nargs="+",
        metavar="ENTITY[:METHOD]=VALUE",
        help=text_fold(
            """
            Specify filters to select the component. Only usable for --input components.
            """
        ),
    )
    parser.add_argument(
        "--groupby",
        nargs="+",
        metavar="ENTITY",
        help=text_fold(
            """
            Specify variable entities that should kept distinct in the output. This
            allows parallel running over multiple subjects, sessions, contrasts, etc.
            Only usable for --input components.
            """
        ),
    )
    parser.add_argument(
        "--aggregate",
        nargs="+",
        metavar="ENTITY",
        help=text_fold(
            """
            Specify variable entities that should merged in the output. This can be 
            used e.g. to create an average image. Only usable for --input components.
            """
        ),
    )
    parser.add_argument(
        "--entities",
        nargs="+",
        action=FilterParse,
        metavar="ENTITY=VALUE",
        help=text_fold(
            """
            Specify the entities of the output file. VALUE may be set to '{WILDCARD}',
            where WILDCARD is a wildcard entity specified in the --groupby argument of
            at least one component. Only usable for --output components.
            """
        ),
    )
    return parser


app = bidsapp.app(
    [
        plugins.BidsArgs(analysis_level=False),
        plugins.Pybidsdb(),
        plugins.Version(distribution="bidsarray"),
    ],
)

component_parser = make_component_parser()


def get_parser():
    """Exposes parser for sphinx doc generation, cwd is the docs dir."""
    return app.build_parser().parser


def parse_component(ix: int, component_args: list[str], config: dict[str, Any]):
    parsed = component_parser.parse_args(component_args)
    if parsed.input is None and parsed.output is None:
        msg = "Each component must specify either --input or --ouput"
        raise ValueError(msg)
    if parsed.input is not None and parsed.output is not None:
        msg = "Only one of --input and --output may be specified for each component"
        raise ValueError(msg)
    label_: str | EllipsisType = (
        parsed.input if parsed.input is not None else parsed.output
    )
    label = ix if label_ is ... else label_
    if parsed.input is not None and parsed.entities is not None:
        msg = f"--entities may not be specified in input component '{label}'"
        raise ValueError(msg)
    if parsed.output is not None:
        invalid_args = {
            "--groupby": parsed.groupby is not None,
            "--aggregate": parsed.aggregate is not None,
            "--filter": parsed.filter is not None,
        }
        msg = ", ".join(key for key, val in invalid_args.items() if val)
        if msg:
            msg += f" may not be specified in output component '{label}'"
            raise ValueError(msg)
    if (
        parsed.aggregate is not None
        and parsed.groupby is not None
        and set(parsed.aggregate) & set(parsed.groupby)
    ):
        msg = f"--aggregate and --groupby must specify different entities in component '{label}'"
        raise ValueError(msg)
    if parsed.input is not None:
        comp_config = config["pybids_inputs"][label]
        if parsed.filter is not None:
            comp_config["filters"] = {}
            for entity, filter_ in parsed.filter.items():
                if isinstance(filter_, (bool, str)):
                    comp_config["filters"][entity] = filter_
        comp_config["wildcards"] = list(
            it.chain(parsed.aggregate or [], parsed.groupby or [])  # type: ignore
        )
        comp_config["groupby"] = parsed.groupby or []
    else:
        comp_config = config["outputs"][label]
        comp_config["entities"] = parsed.entities or {}

    comp_config["order"] = ix


def main():
    main_args, *components = itx.split_at(sys.argv, lambda x: x == ":::")
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
    app.parse_args(main_args[1:])
    grouped_entities = tuple(
        set(
            it.chain.from_iterable(
                pybids_inputs[label]["groupby"] for label in pybids_inputs
            )
        )
    )
    for label, comp in app.config["outputs"].items():
        if set(comp["entities"]) & set(grouped_entities):
            msg = (
                "may not directly specify a groupby variable entity as an output "
                f"entity in component {label}"
            )
            raise ValueError(msg)
    inputs = generate_inputs(
        bids_dir=app.config["bids_dir"],
        pybids_inputs=cast("dict[str, Any]", app.config["pybids_inputs"]),
        pybidsdb_dir=app.config.get("pybidsdb_dir"),
        pybidsdb_reset=app.config.get("pybidsdb_reset", False),
        derivatives=app.config.get("derivatives", False),
        participant_label=app.config.get("participant_label"),
        exclude_participant_label=app.config.get("exclude_participant_label", None),
    )
    table = consensus_table(inputs)

    ncomponents = len(pybids_inputs) + len(app.config["outputs"])
    if labelled:
        out = [""] * ncomponents
        for label, comp in pybids_inputs.items():
            out[comp["order"]] = str(label)
        for label, comp in app.config["outputs"].items():
            out[comp["order"]] = str(label)
        print("\t".join(out))
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
        except BrokenPipeError:
            sys.exit(0)


if __name__ == "__main__":
    main()
