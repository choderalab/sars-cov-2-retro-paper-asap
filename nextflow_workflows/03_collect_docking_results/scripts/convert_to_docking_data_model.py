"""
This script combines the results csvs and generates an input CSV that will be used to run the full cross docking evaluation.
"""
import click
from pathlib import Path
from harbor.analysis.cross_docking import DataFrameModel, DataFrameType, DockingDataModel
import pandas as pd

@click.command()
@click.option(
    "-i",
    "--input-csv",
    required=True,
    help="One or more CSV files containing docking results",
    type=click.Path(exists=True, path_type=Path)
)
@click.option("--output-file-prefix",
              required=True,
              help="Output file suffix")
def main(
    input_csv,
        output_file_prefix
):
    """Convert INPUT_CSV into a DockingDataModel"""
    df = pd.read_csv(input_csv)

    # construct Data
    pose_data = DataFrameModel(dataframe=df, name="PoseData", type=DataFrameType.POSE, key_columns=["Query_Ligand", "Reference_Structure", "Pose_ID"])
    ref_data = DataFrameModel(dataframe=df, name="RefData", type=DataFrameType.REFERENCE, key_columns=["Reference_Structure"])
    lig_data = DataFrameModel(dataframe=df, name="LigData", type=DataFrameType.QUERY, key_columns=["Query_Structure"])
    similarity_data = DataFrameModel(dataframe=df, name="SimilarityData", type=DataFrameType.CHEMICAL_SIMILARITY, key_columns=["Query_Structure", "Reference_Structure", "Aligned", "radius", "bitsize", "fingerprint"])
    scaffold_data = DataFrameModel(dataframe=df, name="ScaffoldData", type=DataFrameType.CHEMICAL_SIMILARITY, key_columns=["Query_Ligand", "Reference_Ligand"])
    data = DockingDataModel.from_models([pose_data, ref_data, lig_data, similarity_data, scaffold_data])

    data.serialize(output_file_prefix)


if __name__ == "__main__":
    main()