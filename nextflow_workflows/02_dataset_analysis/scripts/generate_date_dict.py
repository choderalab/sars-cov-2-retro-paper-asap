"""
This script generates a dictionary of dates for the Fragalysis project.
"""

from pathlib import Path
import argparse
from asapdiscovery.data.util.logging import FileLogger
from datetime import datetime
import pandas as pd
import json


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i",
        "--fragalysis-dir",
        type=Path,
        help="Path to fragalysis_download.",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        default="cmpd_date_dict",
        type=Path,
        help="Path to output directory.",
    )
    return parser.parse_args()


def date_processor(date_string):
    if type(date_string) == str and not date_string == "None":
        try:
            return datetime.strptime(date_string, "%Y-%m-%d %H:%M:%S")
        except ValueError:
            return datetime.strptime(date_string, "%d/%m/%Y %H:%M")
    else:
        return None


def process_crystal_data(soaks):
    ddf = soaks.loc[:, ["Sample Name", "Data Collection Date"]]
    ddf["Sanitized_Date"] = ddf["Data Collection Date"].apply(date_processor)
    ddf.columns = ["Structure_Name", "Data_Collection_Date", "Structure_Date"]
    date_dict = ddf.set_index("Structure_Name").to_dict()["Structure_Date"]
    date_dict = {k: str(v) for k, v in date_dict.items() if str(v) != "NaT"}
    structure_to_cmpd_dict = {
        row["Sample Name"]: row["Compound ID"]
        for idx, row in soaks.iterrows()
        if row["Sample Name"] in date_dict
    }

    return date_dict, structure_to_cmpd_dict


def main():
    args = get_args()
    fragalysis_dir = args.fragalysis_dir
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    logger = FileLogger(
        "generate_date_dict", args.output_dir, logfile="generate_date_dict.log"
    ).getLogger()
    if not fragalysis_dir.exists():
        raise FileNotFoundError(f"Fragalysis directory {fragalysis_dir} not found")

    # load soaks first
    soaks_path = fragalysis_dir / "extra_files" / "Mpro_soaks.csv"
    soaks = pd.read_csv(soaks_path)
    date_dict, structure_to_cmpd_dict = process_crystal_data(soaks)

    # now load the rest of the data
    cocrystals_path = fragalysis_dir / "extra_files" / "Mpro_cocrystallisation.csv"
    cocrystals = pd.read_csv(cocrystals_path)
    co_date_dict, co_structure_to_cmpd_dict = process_crystal_data(cocrystals)

    # confirm no overlap between the two structure datasets
    overlap = set(date_dict.keys()).intersection(set(co_date_dict.keys()))
    if len(overlap) > 0:
        raise ValueError(f"Overlap between soaks and cocrystals: {overlap}")
    date_dict.update(co_date_dict)
    structure_to_cmpd_dict.update(co_structure_to_cmpd_dict)
    logger.info(f"Found {len(date_dict)} structures in the date dict")
    logger.info(
        f"Found {len(structure_to_cmpd_dict)} structures in the structure to cmpd dict"
    )

    # now save both as json files
    with open(args.output_dir / "date_dict.json", "w") as f:
        json.dump(date_dict, f)
    with open(args.output_dir / "structure_to_cmpd_dict.json", "w") as f:
        json.dump(structure_to_cmpd_dict, f)


if __name__ == "__main__":
    main()
