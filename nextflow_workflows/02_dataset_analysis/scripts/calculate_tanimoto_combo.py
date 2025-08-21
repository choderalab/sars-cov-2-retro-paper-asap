"""
Script to calculate the TanimotoCombo between reference and query ligands.
"""

from openeye import oeshape, oechem
from asapdiscovery.data.readers.molfile import MolFileFactory
from asapdiscovery.data.schema.ligand import Ligand
import pandas as pd
import argparse
from pathlib import Path
from asapdiscovery.data.util.logging import FileLogger
import multiprocessing as mp
from chemical_similarity_schema import TanimotoComboSimilarity


def parse_args():
    parser = argparse.ArgumentParser(
        description="Calculate TanimotoCombo between ligands"
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
    return parser.parse_args()


def calculate_one_to_many_tanimoto_oe(
    refmol: oechem.OEMol,
    ref_name: str,
    fitmols: [oechem.OEMol],
    fit_names: [str],
    align: bool = False,
):
    """
    Calculate the Tanimoto coefficient between two molecules using OpenEye's shape toolkit.

    Parameters
    ----------
    refmol : Ligand
        The reference molecule to which the docked molecule is compared.
    fitmol : Ligand
        The docked molecule to be compared to the reference molecule.

    Returns
    -------
    float
        The Tanimoto coefficient between the two molecules.
    """
    similarity_array = []

    if align:
        # OEROCS aligns the molecules before calculating the Tanimoto coefficient
        res = oeshape.OEROCSResult()
        for i, fitmol in enumerate(fitmols):
            oeshape.OEROCSOverlay(res, refmol, fitmol)
            similarity_array.append(
                TanimotoComboSimilarity.from_tanimoto_results(
                    ref=ref_name, query=fit_names[i], results=res, aligned=True
                )
            )
        return similarity_array

    if not align:
        # Prepare reference molecule for calculation
        # With default options this will remove any explicit hydrogens present
        prep = oeshape.OEOverlapPrep()
        prep.Prep(refmol)

        # Get appropriate function to calculate exact shape
        shapeFunc = oeshape.OEOverlapFunc()
        shapeFunc.SetupRef(refmol)

        res = oeshape.OEOverlapResults()
        for i, fitmol in enumerate(fitmols):
            prep.Prep(fitmol)
            shapeFunc.Overlap(fitmol, res)
            similarity_array.append(
                TanimotoComboSimilarity.from_tanimoto_results(
                    ref=ref_name, query=fit_names[i], results=res, aligned=False
                )
            )
        return similarity_array


def parallelize(ref: Ligand, queries: list[Ligand], logger):
    """
    Calculate the MCS between a reference ligand and a list of query ligands.
    :param ref: Reference ligand
    :param queries: List of query ligands
    :return: Dataframe with MCS results
    """
    logger.info(f"Calculating MCS for {ref.compound_name}...")
    refmol = ref.to_oemol()
    query_mols = [query.to_oemol() for query in queries]
    similarity_array_aligned = calculate_one_to_many_tanimoto_oe(
        refmol,
        ref.compound_name,
        query_mols,
        [query.compound_name for query in queries],
        align=True,
    )
    similarity_array_not_aligned = calculate_one_to_many_tanimoto_oe(
        refmol,
        ref.compound_name,
        query_mols,
        [query.compound_name for query in queries],
        align=False,
    )
    df = pd.concat(
        [
            TanimotoComboSimilarity.construct_dataframe(similarity_array_aligned),
            TanimotoComboSimilarity.construct_dataframe(similarity_array_not_aligned),
        ],
        ignore_index=True,
    )
    return df


def main():
    args = parse_args()
    output_dir = args.output_dir
    output_dir.mkdir(exist_ok=True, parents=True)

    logger = FileLogger(
        "calculate_tanimoto_combo",
        output_dir,
        logfile="calculate_tanimoto_combo.log",
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

    logger.info("Calculating similarities...")
    # Parallelize the MCS calculation
    with mp.Pool(mp.cpu_count()) as pool:
        results = pool.starmap(
            parallelize, [(ref, queries, logger) for ref in references]
        )
        results = [result for result in results if result is not None]

    # Save results
    logger.info("Saving results...")
    df = pd.concat(results, ignore_index=True)
    output_path = output_dir / "tanimoto_combo.csv"
    df.to_csv(output_path, index=False)
    logger.info(f"Results saved to {output_path}")


if __name__ == "__main__":
    main()
