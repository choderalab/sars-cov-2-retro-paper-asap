"""
This script creates the calculations that will be run on a cross-docking dataset
"""

import click
from pathlib import Path

import pandas as pd

from harbor.analysis.cross_docking import (
    DockingDataModel,
    EvaluatorFactory,
    ScaffoldSplitOptions,
)
from harbor.analysis.utils import FileLogger


def save_and_create_evs(
    evf: EvaluatorFactory, data: DockingDataModel, name: str, output: Path, logger
):
    evf.to_yaml_file(output)
    evs = evf.create_evaluators(data)
    logger.info(f"created {len(evs)} for {name}")
    for i, evaluator in enumerate(evs):
        evaluator.to_json_file(output / f"evaluator_{name}_{i}.json")

    df = pd.DataFrame.from_records([ev.get_records() for ev in evs])
    df.to_csv(output / f"{name}_evaluators.csv")


@click.command()
@click.option(
    "-i",
    "--input-parquet",
    required=True,
    help="Path to input parquet file made by DockinDataModel",
    type=click.Path(exists=True, path_type=Path),
)
@click.option(
    "--settings",
    required=True,
    help="Path to serialized EvaluatorFactory settings object",
    type=click.Path(exists=True, path_type=Path),
)
@click.option(
    "-o",
    "--output",
    type=Path,
    required=True,
    help="Path to the output directory where the results will be stored",
)
def main(input_parquet, settings, output):
    # load evaluator factory
    evf = EvaluatorFactory.from_yaml_file(settings)
    output = output / evf.name
    output.mkdir(exist_ok=True, parents=True)

    logger = FileLogger(
        logname="create_evaluators", path=output, logfile="create_evaluators.log"
    ).getLogger()
    logger.info(f"Reading data model from {input_parquet}")
    # load docking model
    data = DockingDataModel.deserialize(input_parquet)
    save_and_create_evs(evf, data, evf.name, output, logger)


if __name__ == "__main__":
    main()
