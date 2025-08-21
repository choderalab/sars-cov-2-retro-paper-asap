"""
Script to generate ECFP fingerprint similarity analysis between reference and query ligands.
"""

from openeye import oegraphsim
from asapdiscovery.data.readers.molfile import MolFileFactory
import pandas as pd
import itertools
import argparse
from pathlib import Path
from asapdiscovery.data.util.logging import FileLogger
from chemical_similarity_schema import ECFPSimilarity


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate ECFP fingerprint similarities between ligands"
    )
    parser.add_argument(
        "--ref-ligand-sdf",
        type=Path,
        required=True,
        help="Path to directory containing prepped reference ligand sdf.",
    )
    parser.add_argument(
        "--output-dir", required=True, type=Path, help="Path to output directory"
    )
    parser.add_argument(
        "--query-ligand-sdf",
        type=Path,
        required=False,
        help="Path to directory containing prepped query ligand sdf. If false, ref-ligand-sdf will be used.",
    )
    parser.add_argument(
        "--settings", type=Path, required=False, help="Path to settings json file"
    )
    return parser.parse_args()


def get_fp(mol, bit_size=2048, radius=2):
    fp = oegraphsim.OEFingerPrint()
    oegraphsim.OEMakeCircularFP(
        fp,
        mol,
        bit_size,
        0,
        radius,
        oegraphsim.OEFPAtomType_DefaultCircularAtom,
        oegraphsim.OEFPBondType_DefaultCircularBond,
    )
    return fp


def calculate_tanimoto(fp1, fp2):
    return oegraphsim.OETanimoto(fp1, fp2)


def main():
    args = parse_args()
    output_dir = args.output_dir
    output_dir.mkdir(exist_ok=True, parents=True)

    logger = FileLogger(
        "calculate_ecfp_tanimoto",
        output_dir,
        logfile="calculate_ecfp_tanimoto.log",
    ).getLogger()

    logger.info("Loading molecules...")
    references = MolFileFactory(filename=args.ref_ligand_sdf).load()
    queries = (
        MolFileFactory(filename=args.query_ligand_sdf).load()
        if args.query_ligand_sdf
        else references.copy()
    )
    logger.info(f"Loaded {len(references)} reference molecules.")
    logger.info(f"Loaded {len(queries)} query molecules.")

    # Create or load settings
    # radii = [2, 3, 4, 5]
    radii = [2, 5]
    # bit_sizes = [1024, 2048]
    bit_sizes = [2048]
    dfs = []
    for radius, bit_size in itertools.product(radii, bit_sizes):
        logger.info(
            f"Calculating similarities for radius {radius} and bit size {bit_size}"
        )

        # Create all pairs and calculate similarities
        logger.info("Calculating fingerprints...")
        ref_fps = {
            mol.compound_name: get_fp(mol.to_oemol(), bit_size, radius)
            for mol in references
        }
        query_fps = {
            mol.compound_name: get_fp(mol.to_oemol(), bit_size, radius)
            for mol in queries
        }

        logger.info("Calculating similarities...")
        similarities = [
            ECFPSimilarity(
                Reference_Ligand=ref,
                Query_Ligand=query,
                Tanimoto=calculate_tanimoto(ref_fps[ref], query_fps[query]),
                radius=radius,
                bitsize=bit_size,
            )
            for ref, query in itertools.product(ref_fps.keys(), query_fps.keys())
        ]
        df = ECFPSimilarity.construct_dataframe(similarities)
        dfs.append(df)

    # Save results
    logger.info("Saving results...")
    df = pd.concat(dfs, ignore_index=True)
    output_path = output_dir / "fingerprint_similarities.csv"
    df.to_csv(output_path, index=False)
    logger.info(f"Results saved to {output_path}")


if __name__ == "__main__":
    main()
