"""
Takes the previously calculated chemical similarity data and outputs a CSV file with the rest of the docking results
"""

import pandas as pd
import argparse
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine all relevant data for analysis"
    )

    # add any number of csv files
    parser.add_argument("csvs", nargs="+")
    parser.add_argument("--output-dir", type=Path, required=False, default="./")
    return parser.parse_args()


def main():
    args = parse_args()
    output_dir = args.output_dir
    output_dir.mkdir(exist_ok=True, parents=True)
    dfs = []
    for csv in args.csvs:
        dfs.append(pd.read_csv(csv))
    df = pd.concat(dfs)
    df.to_csv(output_dir / "combined_chemical_similarity_data.csv", index=False)


if __name__ == "__main__":
    main()
