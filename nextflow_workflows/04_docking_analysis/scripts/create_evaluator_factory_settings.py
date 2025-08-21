"""
This script creates the settings files used to create the analysis
"""

import click
from pathlib import Path
from harbor.analysis.cross_docking import (
    EvaluatorFactory,
    ScaffoldSplitOptions,
)


@click.command()
@click.option(
    "-o",
    "--output",
    type=Path,
    required=False,
    default="./",
    help="Path to the output directory where the results will be stored",
)
def main(output):
    output.mkdir(exist_ok=True, parents=True)

    # update default settings
    default = EvaluatorFactory(name="default")

    default.success_rate_evaluator_settings.use = True
    default.success_rate_evaluator_settings.success_rate_column = "PoseData_RMSD"

    default.scorer_settings.rmsd_scorer_settings.use = True
    default.scorer_settings.rmsd_scorer_settings.rmsd_column_name = "PoseData_RMSD"

    default.scorer_settings.posit_scorer_settings.use = True
    default.scorer_settings.posit_scorer_settings.posit_score_column_name = (
        "PoseData_docking-confidence-POSIT"
    )

    # basic date split cross docking
    evf = default.__deepcopy__()
    evf.name = "reference_split_comparison"
    evf.reference_split_settings.use = True
    evf.reference_split_settings.date_split_settings.use = True
    evf.reference_split_settings.date_split_settings.reference_structure_date_column = (
        "RefData_Date"
    )
    evf.reference_split_settings.random_split_settings.use = True
    evf.reference_split_settings.update_reference_settings.use = True
    evf.reference_split_settings.update_reference_settings.use_logarithmic_scaling = (
        True
    )
    evf.to_yaml_file(output)

    # Scaffold split options
    default_scaffold = default.__deepcopy__()
    default_scaffold.name = "default_scaffold_settings"
    default_scaffold.pairwise_split_settings.use = True
    default_scaffold.pairwise_split_settings.scaffold_split_settings.use = True
    default_scaffold.pairwise_split_settings.scaffold_split_settings.reference_scaffold_id_column = (
        "RefData_Scaffold_ID"
    )
    default_scaffold.pairwise_split_settings.scaffold_split_settings.query_scaffold_id_column = (
        "QueryData_Scaffold_ID"
    )
    default_scaffold.pairwise_split_settings.scaffold_split_settings.reference_scaffold_min_count = (
        1
    )
    default_scaffold.pairwise_split_settings.scaffold_split_settings.query_scaffold_min_count = (
        1
    )

    # x_to_x default
    x_to_x_default = default_scaffold.__deepcopy__()
    x_to_x_default.name = "x_to_x_scaffold_split"
    x_to_x_default.pairwise_split_settings.scaffold_split_settings.scaffold_split_option = (
        ScaffoldSplitOptions.X_TO_X
    )
    x_to_x_default.to_yaml_file(output)

    # x_to_x 5 refs
    evf = x_to_x_default.__deepcopy__()
    evf.name = "x_to_x_scaffold_split_5_refs"
    evf.pairwise_split_settings.scaffold_split_settings.scaffold_split_option = (
        ScaffoldSplitOptions.X_TO_X
    )
    evf.dataset_before_similarity = False
    evf.combine_reference_and_similarity_splits = True
    evf.reference_split_settings.use = True
    evf.reference_split_settings.random_split_settings.use = True
    evf.reference_split_settings.n_reference_structures = [5]
    evf.pairwise_split_settings.scaffold_split_settings.reference_scaffold_min_count = 5
    evf.pairwise_split_settings.scaffold_split_settings.query_scaffold_min_count = 5
    evf.to_yaml_file(output)

    # x_to_y default
    x_to_y_default = default_scaffold.__deepcopy__()
    x_to_y_default.name = "x_to_y_scaffold_split"
    x_to_y_default.pairwise_split_settings.scaffold_split_settings.scaffold_split_option = (
        ScaffoldSplitOptions.X_TO_Y
    )
    x_to_y_default.to_yaml_file(output)

    # x to y scaffold split with 5 refs
    evf = x_to_y_default.__deepcopy__()
    evf.name = "x_to_y_scaffold_split_5_refs"
    evf.dataset_before_similarity = False
    evf.combine_reference_and_similarity_splits = True
    evf.reference_split_settings.use = True
    evf.reference_split_settings.random_split_settings.use = True
    evf.reference_split_settings.n_reference_structures = [5]
    evf.pairwise_split_settings.scaffold_split_settings.reference_scaffold_min_count = 5
    evf.pairwise_split_settings.scaffold_split_settings.query_scaffold_min_count = 5

    evf.to_yaml_file(output)

    # x to not x scaffold split
    x_to_not_x = default_scaffold.__deepcopy__()
    x_to_not_x.name = "x_to_not_x_scaffold_split"
    x_to_not_x.pairwise_split_settings.scaffold_split_settings.scaffold_split_option = (
        ScaffoldSplitOptions.X_TO_NOT_X
    )
    x_to_not_x.dataset_before_similarity = True
    x_to_not_x.combine_reference_and_similarity_splits = True
    x_to_not_x.reference_split_settings.use = True
    x_to_not_x.reference_split_settings.random_split_settings.use = True
    x_to_not_x.reference_split_settings.update_reference_settings.use = True
    x_to_not_x.reference_split_settings.update_reference_settings.use_logarithmic_scaling = (
        True
    )
    x_to_not_x.to_yaml_file(output)

    # not x to x scaffold split
    not_x_to_x = default_scaffold.__deepcopy__()
    not_x_to_x.name = "not_x_to_x_scaffold_split"
    not_x_to_x.pairwise_split_settings.scaffold_split_settings.scaffold_split_option = (
        ScaffoldSplitOptions.NOT_X_TO_X
    )
    not_x_to_x.to_yaml_file(output)

    # not x to x scaffold split w 5 refs
    not_x_to_x = default_scaffold.__deepcopy__()
    not_x_to_x.name = "not_x_to_x_scaffold_split_5_refs"
    not_x_to_x.pairwise_split_settings.scaffold_split_settings.scaffold_split_option = (
        ScaffoldSplitOptions.NOT_X_TO_X
    )
    not_x_to_x.combine_reference_and_similarity_splits = True
    not_x_to_x.dataset_before_similarity = False
    not_x_to_x.reference_split_settings.use = True
    not_x_to_x.reference_split_settings.random_split_settings.use = True
    not_x_to_x.reference_split_settings.n_reference_structures = [5]
    not_x_to_x.pairwise_split_settings.scaffold_split_settings.reference_scaffold_min_count = (
        5
    )
    not_x_to_x.pairwise_split_settings.scaffold_split_settings.query_scaffold_min_count = (
        5
    )
    not_x_to_x.to_yaml_file(output)

    # similarity split
    sim_split = default.__deepcopy__()
    sim_split.name = "increasing_similarity_tanimoto_combo_aligned"
    sim_split.pairwise_split_settings.use = True
    sim_split.pairwise_split_settings.similarity_split_settings.use = True
    sim_split.pairwise_split_settings.similarity_split_settings.similarity_column_name = (
        "TanimotoComboData_Tanimoto"
    )
    sim_split.pairwise_split_settings.similarity_split_settings.include_similar = False
    sim_split.pairwise_split_settings.similarity_split_settings.similarity_groupby_dict = {
        "TanimotoComboData_Type": "TanimotoCombo",
        "TanimotoComboData_Aligned": True,
    }
    sim_split.pairwise_split_settings.similarity_split_settings.update_reference_settings.use = (
        True
    )
    sim_split.pairwise_split_settings.similarity_split_settings.update_reference_settings.use_logarithmic_scaling = (
        True
    )
    sim_split.to_yaml_file(output)

    # MCS
    sim_split = default.__deepcopy__()
    sim_split.name = "increasing_similarity_mcs"
    sim_split.pairwise_split_settings.use = True
    sim_split.pairwise_split_settings.similarity_split_settings.use = True
    sim_split.pairwise_split_settings.similarity_split_settings.include_similar = False
    sim_split.pairwise_split_settings.similarity_split_settings.similarity_column_name = (
        "MCSData_Tanimoto"
    )
    sim_split.pairwise_split_settings.similarity_split_settings.similarity_groupby_dict = {
        "MCSData_Type": "MCS"
    }
    sim_split.pairwise_split_settings.similarity_split_settings.update_reference_settings.use = (
        True
    )
    sim_split.pairwise_split_settings.similarity_split_settings.update_reference_settings.use_logarithmic_scaling = (
        True
    )
    sim_split.to_yaml_file(output)

    # ECFP
    sim_split = default.__deepcopy__()
    sim_split.name = "increasing_similarity_ecfp4"
    sim_split.pairwise_split_settings.use = True
    sim_split.pairwise_split_settings.similarity_split_settings.use = True
    sim_split.pairwise_split_settings.similarity_split_settings.include_similar = False
    sim_split.pairwise_split_settings.similarity_split_settings.similarity_column_name = (
        "ECFPData_Tanimoto"
    )
    sim_split.pairwise_split_settings.similarity_split_settings.similarity_groupby_dict = {
        "ECFPData_fingerprint": "ECFP4_2048"
    }
    sim_split.pairwise_split_settings.similarity_split_settings.update_reference_settings.use = (
        True
    )
    sim_split.pairwise_split_settings.similarity_split_settings.update_reference_settings.use_logarithmic_scaling = (
        True
    )
    sim_split.to_yaml_file(output)


if __name__ == "__main__":
    main()
