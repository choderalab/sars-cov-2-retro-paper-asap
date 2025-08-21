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
    name = "reverse_similarity_split"
    output = output / name
    output.mkdir(exist_ok=True, parents=True)
    n_refs = [1, 2, 5, 10, 20, 30, 40, 50, 100, 137, 200, 300, 403]
    n_most_similars = [1, 2, 5, 10, 50]
    ref_structure_column = "Reference_Structure"
    dataset_splits = []
    for n_ref in n_refs:
        dataset_splits.append(
            cd.DateSplit(
                date_column="RefData_Date",
                randomize_by_n_days=1,
                n_reference_structures=n_ref,
                reference_structure_column=ref_structure_column,
            )
        )
        dataset_splits.append(
            cd.RandomSplit(
                reference_structure_column=ref_structure_column,
                n_reference_structures=n_ref,
            )
        )
    similarity_splits = []
    for n_most_similar in n_most_similars:
        similarity_splits.append(
            cd.SimilaritySplit(
                groupby={"MCSData_Type": "MCS"},
                similarity_column="MCSData_Tanimoto",
                query_ligand_column="Query_Ligand",
                reference_ligand_column="RefData_Ligand",
                include_similar=True,
                higher_is_more_similar=True,
                sort_instead_of_threshold=True,
                split_level=1,
                n_similar=n_most_similar,
            )
        )
        similarity_splits.append(
            cd.SimilaritySplit(
                similarity_column="TanimotoComboData_Tanimoto",
                groupby={
                    "TanimotoComboData_Type": "TanimotoCombo",
                    "TanimotoComboData_Aligned": True,
                },
                query_ligand_column="Query_Ligand",
                reference_ligand_column="RefData_Ligand",
                include_similar=True,
                higher_is_more_similar=True,
                sort_instead_of_threshold=True,
                split_level=1,
                n_similar=n_most_similar,
            )
        )
        similarity_splits.append(
            cd.SimilaritySplit(
                similarity_column="ECFPData_Tanimoto",
                groupby={"ECFPData_fingerprint": "ECFP4_2048"},
                query_ligand_column="Query_Ligand",
                reference_ligand_column="RefData_Ligand",
                include_similar=True,
                higher_is_more_similar=True,
                sort_instead_of_threshold=True,
                split_level=1,
                n_similar=n_most_similar,
            )
        )
    scorers = [
        cd.POSITScorer(variable="PoseData_docking-confidence-POSIT"),
        cd.RMSDScorer(variable="PoseData_RMSD", cutoff=2),
    ]
    evs = []
    for scorer in scorers:
        for dataset_split in dataset_splits:
            for similarity_split in similarity_splits:
                # need to do this because N_Reference_Structures is getting overwritten and I don't want to go back through the code and change it yet
                # Just setting it to the n_reference_structures of the dataset_split without copying the similarity split does annoying python pointer things (it ends up setting all the ref splits to the same value
                similarity_split = similarity_split.model_copy()
                similarity_split.n_reference_structures = (
                    dataset_split.n_reference_structures
                )
                ev = cd.Evaluator(
                    pose_selector=cd.PoseSelector(
                        name="PoseSelector",
                        variable="Pose_ID",
                        ascending=True,
                        number_to_return=1,
                    ),
                    dataset_split=dataset_split,
                    similarity_split=similarity_split,
                    dataset_before_similarity=True,
                    scorer=scorer,
                    evaluator=cd.BinaryEvaluation(variable="PoseData_RMSD", cutoff=2),
                    n_bootstraps=1000,
                )
                evs.append(ev)

    for i, evaluator in enumerate(evs):
        evaluator.to_json_file(output / f"evaluator_{name}_{i}.json")

    df = pd.DataFrame.from_records([ev.get_records() for ev in evs])
    df.to_csv(output / f"{name}_evaluators.csv")


if __name__ == "__main__":
    main()
