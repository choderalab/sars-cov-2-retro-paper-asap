"""
Scriptify the Bemis-Murcko clustering of ligands using a few different definitions of the cluster
Influenced by:
https://github.com/rdkit/rdkit/discussions/6844

Requires asapdiscovery environment
"""

import pandas as pd
from asapdiscovery.data.schema.ligand import Ligand
from asapdiscovery.data.readers.molfile import MolFileFactory
from rdkit import Chem

PATT = Chem.MolFromSmarts("[$([D1]=[*])]")
REPL = Chem.MolFromSmarts("[*]")

from rdkit.Chem.Scaffolds import MurckoScaffold
from collections import defaultdict
from argparse import ArgumentParser
from pathlib import Path

from pydantic import BaseModel
from abc import abstractmethod


def parse_args():
    parser = ArgumentParser(description="Bemis-Murcko clustering of ligands")
    parser.add_argument(
        "--sdf-2d", type=Path, help="Path to the input 2D SDF file", nargs="+"
    )
    parser.add_argument(
        "--output-dir", default="./", type=Path, help="Path to the output directory"
    )
    return parser.parse_args()


class BaseBemisMurckoScaffold(BaseModel):
    name: str

    @abstractmethod
    def run(self, ligand: Ligand) -> str:
        """
        Run the Bemis-Murcko clustering on the ligands
        """
        pass


class DefaultRDKitBemisMurckoScaffold(BaseBemisMurckoScaffold):
    name = "default"

    def run(self, ligand: Ligand) -> str:
        """
        Run the Bemis-Murcko clustering on the ligands
        """
        mol = ligand.to_rdkit()
        scaff = MurckoScaffold.GetScaffoldForMol(mol)
        return Chem.MolToSmiles(scaff)


class BajorathBemisMurckoScaffold(BaseBemisMurckoScaffold):
    name = "bajorath"

    def run(self, ligand: Ligand) -> str:
        """
        Run the Bemis-Murcko clustering on the ligands
        """
        mol = ligand.to_rdkit()
        scaff = MurckoScaffold.GetScaffoldForMol(mol)
        scaff = Chem.rdmolops.DeleteSubstructs(scaff, PATT)
        return Chem.MolToSmiles(scaff)


class GenericBemisMurckoScaffold(BaseBemisMurckoScaffold):
    name = "generic"

    def run(self, ligand: Ligand) -> str:
        """
        Run the Bemis-Murcko clustering on the ligands
        :param ligands:
        :return:
        """
        mol = ligand.to_rdkit()
        scaff = MurckoScaffold.GetScaffoldForMol(mol)
        scaff = MurckoScaffold.MakeScaffoldGeneric(scaff)
        return Chem.MolToSmiles(scaff)


class CSKBemisMurckoScaffold(BaseBemisMurckoScaffold):
    name = "csk"

    def run(self, ligand: Ligand) -> str:
        """
        Run the Bemis-Murcko clustering on the ligands
        :param ligands:
        :return:
        """
        mol = ligand.to_rdkit()
        scaff = MurckoScaffold.GetScaffoldForMol(mol)
        scaff = Chem.rdmolops.ReplaceSubstructs(scaff, PATT, REPL, replaceAll=True)[0]
        scaff = MurckoScaffold.MakeScaffoldGeneric(scaff)
        scaff = MurckoScaffold.GetScaffoldForMol(scaff)

        return Chem.MolToSmiles(scaff)


def split_by_scaffold(ligands, scaffold_type: BaseBemisMurckoScaffold):
    """
    Split ligands by scaffold.
    """

    scaffolds = defaultdict(list)
    for ligand in ligands:
        scaffold = scaffold_type.run(ligand)
        scaffolds[scaffold].append(ligand)
    scaffold_list = [
        {"scaffold": scaffold, "ligands": ligands}
        for scaffold, ligands in scaffolds.items()
    ]
    return sorted(scaffold_list, key=lambda x: len(x["ligands"]), reverse=True)


def main():
    args = parse_args()

    # create output directory
    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    ligands = []
    for sdf in args.sdf_2d:
        mff = MolFileFactory(filename=sdf)
        ligands.extend(mff.load())

    scaffold_types = [
        DefaultRDKitBemisMurckoScaffold(),
        BajorathBemisMurckoScaffold(),
        GenericBemisMurckoScaffold(),
        CSKBemisMurckoScaffold(),
    ]

    for scaffold_type in scaffold_types:
        print(f"Running {scaffold_type.name}")
        scaffold_list = split_by_scaffold(ligands, scaffold_type)
        cluster_labels = []
        for i, scaffold_dict in enumerate(scaffold_list):
            for ligand in scaffold_dict["ligands"]:
                cluster_labels.append(
                    dict(
                        compound_name=ligand.compound_name,
                        cluster_id=i,
                        scaffold_smarts=scaffold_dict["scaffold"],
                        cluster_type=f"{scaffold_type.name}_bemis_murko",
                    )
                )

        cluster_df = pd.DataFrame.from_records(cluster_labels)
        cluster_df.to_csv(
            output_dir / f"{scaffold_type.name}_cluster_labels.csv", index=False
        )


if __name__ == "__main__":
    main()
