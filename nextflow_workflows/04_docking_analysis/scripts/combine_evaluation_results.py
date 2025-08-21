import pandas as pd
import glob
import os
import click


@click.command()
@click.argument(
    "input-csvs",
    nargs=-1,
    type=click.Path(exists=True),
)
@click.argument("output-file", type=click.Path())
def combine_csv_files(input_csvs, output_file):
    """Combine multiple INPUT_CSVS into a single OUTPUT_FILE."""
    csv_files = list(input_csvs)

    # Create empty list to store dataframes
    dfs = []

    # Read each CSV file and append to list
    for file in csv_files:
        df = pd.read_csv(file)
        dfs.append(df)

    # Combine all dataframes
    combined_df = pd.concat(dfs, ignore_index=True)

    # Save combined dataframe
    combined_df.to_csv(output_file, index=False)
    click.echo(f"Combined {len(csv_files)} files into {output_file}")


if __name__ == "__main__":
    combine_csv_files()
