"""
This script creates the calculations that will be run on a cross-docking dataset
"""

import click
from pathlib import Path
import pandas as pd
import harbor.analysis.cross_docking as cd


@click.command()
@click.option(
    "-i",
    "--input-parquet",
    required=True,
    help="Path to input parquet file made by DockingDataModel",
    type=click.Path(exists=True, path_type=Path),
)
@click.option(
    "-o",
    "--output",
    type=Path,
    required=False,
    default=Path("./"),
    help="Path to the output directory where the results will be stored",
)
def main(input_parquet, output):
    name = "scaffold_datesplit_evaluators"

    data = cd.DockingDataModel.deserialize(input_parquet)

    output = output / name
    output.mkdir(exist_ok=True, parents=True)

    ref_structure_column = "Reference_Structure"

    n_refs_per_scaffold_list = [1, 2, 5, 10, 20, -1]
    n_refs = [1, 2, 5, 10, 20, 30, 40, 50, 100, 137, 200, 300, 403]
    n_poses = [1]
    pose_selectors = [
        cd.PoseSelector(
            name="PoseSelector", variable="Pose_ID", ascending=True, number_to_return=n
        )
        for n in n_poses
    ]

    scorers = [
        cd.POSITScorer(variable="PoseData_docking-confidence-POSIT"),
        cd.RMSDScorer(variable="PoseData_RMSD", cutoff=2),
    ]

    dataset_splits = []
    for n in n_refs:
        dataset_splits.append(
            cd.DateSplit(
                date_column="RefData_Date",
                randomize_by_n_days=1,
                n_reference_structures=n,
                reference_structure_column=ref_structure_column,
            )
        )
        dataset_splits.append(
            cd.RandomSplit(
                reference_structure_column=ref_structure_column,
                n_reference_structures=n,
            )
        )
        for n_refs_per_scaffold in n_refs_per_scaffold_list:
            unique_refs = (
                data.dataframe.sort_values("RefData_Date")
                .groupby("RefData_Scaffold_ID")
                .head(n_refs_per_scaffold)[ref_structure_column]
                .unique()
            )
            if len(unique_refs) >= n:

                # Only add the split if we have enough scaffolds with enough references
                dataset_splits.append(
                    cd.ScaffoldDateSplit(
                        date_column="RefData_Date",
                        scaffold_id_column="RefData_Scaffold_ID",
                        randomize_by_n_days=1,
                        n_reference_structures=n,
                        reference_structure_column=ref_structure_column,
                        n_refs_per_scaffold=n_refs_per_scaffold,
                    )
                )
                if len(unique_refs) > n:
                    # If not exactly one of the n_refs, add another split
                    # with the actual number of unique refs
                    dataset_splits.append(
                        cd.ScaffoldDateSplit(
                            date_column="RefData_Date",
                            scaffold_id_column="RefData_Scaffold_ID",
                            randomize_by_n_days=1,
                            n_reference_structures=len(unique_refs),
                            reference_structure_column=ref_structure_column,
                            n_refs_per_scaffold=n_refs_per_scaffold,
                        )
                    )

    evs = []
    for pose_selector in pose_selectors:
        for scorer in scorers:
            for dataset_split in dataset_splits:
                ev = cd.Evaluator(
                    pose_selector=pose_selector,
                    dataset_split=dataset_split,
                    scorer=scorer,
                    evaluator=cd.BinaryEvaluation(variable="PoseData_RMSD", cutoff=2),
                    n_bootstraps=1000,
                )
                evs.append(ev)

    for i, evaluator in enumerate(evs):
        evaluator.to_json_file(output / f"evaluator_{name}_{i}.json")

    df = pd.DataFrame.from_records([ev.get_records() for ev in evs])
    df.to_csv(output / f"{name}.csv")


if __name__ == "__main__":
    main()
