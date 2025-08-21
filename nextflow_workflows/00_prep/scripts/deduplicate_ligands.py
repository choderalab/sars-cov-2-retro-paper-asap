#!/usr/bin/env python
import click

"""
Example usage:
python deduplicate_ligands.py \
    --fragalysis-dir /data1/choderaj/paynea/asap-datasets/full_cross_dock_v2/mpro_fragalysis-04-01-24_curated \
    --prepped-path /data1/choderaj/paynea/asap-datasets/full_cross_dock_v2/mpro_fragalysis-04-01-24_curated_cache \
    --output-dir /data1/choderaj/paynea/asap-datasets/full_cross_dock_v2/mpro_fragalysis-04-01-24_curated_cache_fixed \
    --remove-covalent
"""
from pathlib import Path
import pandas as pd
from datetime import datetime
import shutil
from asapdiscovery.modeling.protein_prep import PreppedComplex
from harbor.analysis.utils import FileLogger


def get_duplicates(df):
    from itertools import combinations

    return_dict = {}
    for col1, col2 in combinations(df.columns, 2):
        counts = df.groupby(col1).nunique()
        return_dict[f"{col1}_to_{col2}"] = counts[(counts[col2] > 1)].index.unique()
    return return_dict


def date_processor(date_string):
    if isinstance(date_string, str) and date_string != "None":
        try:
            return datetime.strptime(date_string, "%Y-%m-%d %H:%M:%S")
        except ValueError:
            return datetime.strptime(date_string, "%d/%m/%Y %H:%M")
    return None


def process_crystal_data(soaks):
    ddf = soaks.loc[:, ["Sample Name", "Data Collection Date"]]
    ddf["Sanitized_Date"] = ddf["Data Collection Date"].apply(date_processor)
    ddf.columns = ["Structure_Name", "Data_Collection_Date", "Structure_Date"]
    date_dict = ddf.set_index("Structure_Name").to_dict()["Structure_Date"]
    date_dict = {k: str(v) for k, v in date_dict.items() if str(v) != "NaT"}
    return date_dict


def get_records_from_complexes(complexes):
    records = []
    for c in complexes:
        records.append(
            {
                "SMILES": c.ligand.smiles,
                "Compound_Name": c.ligand.compound_name,
                "Target_Name": c.target.target_name,
            }
        )
    return records


@click.command()
@click.option(
    "--fragalysis-dir",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    required=True,
    help="Directory containing Fragalysis data",
)
@click.option(
    "--prepped-path",
    type=click.Path(exists=True, file_okay=False, dir_okay=True),
    required=True,
    help="Directory containing prepped complexes",
)
@click.option(
    "--output-dir",
    type=click.Path(file_okay=False, dir_okay=True),
    required=True,
    help="Output directory for filtered structures",
)
@click.option(
    "--remove-covalent", is_flag=True, default=False, help="Remove covalent ligands"
)
def main(fragalysis_dir, prepped_path, output_dir, remove_covalent):
    """Filter and copy protein structures based on deduplication criteria."""
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    logger = FileLogger(
        "deduplicate_ligands", path=output_path, logfile="deduplicate_ligands.log"
    ).getLogger()

    pcs_to_load = list(Path(prepped_path).glob("./*/*.json"))
    if not pcs_to_load:
        logger.info("No prepped complexes found to load.")
        return

    logger.info(f"Found {len(pcs_to_load)} prepped complexes to load.")

    # Load Fragalysis data
    pcs = [PreppedComplex.from_json_file(f) for f in pcs_to_load]

    # Create initial dataframe
    df = pd.DataFrame.from_records(get_records_from_complexes(pcs))

    # Remove covalent ligands if specified
    if remove_covalent:
        data = pd.read_csv(
            Path(fragalysis_dir) / "extra_files" / "Mpro_compound_tracker_csv.csv"
        )
        relevant_data = data[data["Compound ID"].isin(df.Compound_Name.unique())]
        suspected_covalent = relevant_data[
            relevant_data.why_suspected_SMILES == "Covalent"
        ]["Compound ID"].unique()
        noncovalent = df[~df.Compound_Name.isin(suspected_covalent)]

        covalent_target_names = set(df.Target_Name.unique()) - set(
            noncovalent.Target_Name.unique()
        )
        logger.info(
            f"Removing {len(suspected_covalent)} covalent compounds with {len(covalent_target_names)} associated crystal structures."
        )
        df = noncovalent

    # Load and process dates
    soaks_path = Path(fragalysis_dir) / "extra_files" / "Mpro_soaks.csv"
    soaks = pd.read_csv(soaks_path)
    date_dict = process_crystal_data(soaks)

    # Add dates and find structures to keep
    df["Date"] = df.Target_Name.apply(lambda x: date_dict.get(x[:-3], None))

    duplicates = get_duplicates(df)
    for key, value in duplicates.items():
        logger.info(f"Duplicate {key}: {len(value)} entries")

    deduped = df.sort_values("Date").groupby("SMILES").head(1)
    deduped = deduped.sort_values("Date").groupby("Compound_Name").head(1)

    duplicates = get_duplicates(deduped)

    failed = False
    for key, value in duplicates.items():
        if len(value) > 1:
            # Only report duplicates with more than one entry
            logger.info(f"Duplicate {key}: {len(value)} entries")
            failed = True
    if failed:
        raise ValueError("Deduplication failed due to remaining duplicates.")

    all_targets_to_keep = set(deduped.Target_Name.unique())

    logger.info(
        f"Keeping {len(all_targets_to_keep)} unique targets after deduplication."
    )

    # Copy files, skipping removed targets
    prepped_path = Path(prepped_path)
    copied = 0
    skipped = 0

    for src_path in prepped_path.glob("*/*.json"):
        src_dir = src_path.parent
        target_name = src_dir.name
        # Check if the target should be removed (starts with any name in targets_to_remove)
        if any(target_name.startswith(x) for x in all_targets_to_keep):
            dest_dir = output_path / target_name
            shutil.copytree(src_dir, dest_dir, dirs_exist_ok=True)
            copied += 1
        else:
            skipped += 1

    logger.info(f"Copied {copied} files")
    logger.info(f"Skipped {skipped} files")
    logger.info(f"Total files processed: {copied + skipped}")


if __name__ == "__main__":
    main()
