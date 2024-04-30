from __future__ import annotations
from types import EllipsisType
from typing import Any
import argparse
from snakebids.utils.utils import text_fold
from snakebids.plugins.component_edit import FilterParse
import functools as ft
import itertools as it


@ft.lru_cache
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


def parse_component(ix: int, component_args: list[str], config: dict[str, Any]):
    component_parser = make_component_parser()
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
