from pydantic import BaseModel, Field, confloat, root_validator, ValidationError
from enum import Enum
import json
import pandas as pd


class MoleculeSimilarity(BaseModel):
    """
    MoleculeSimilarity
    """

    Reference_Ligand: str = Field(..., description="Molecule 1")
    Query_Ligand: str = Field(..., description="Molecule 2")
    Type: str = Field(..., description="Type of similarity")
    Tanimoto: confloat(ge=0, le=1) = Field(
        ..., description="Tanimoto similarity from 0 to 1"
    )

    def save(self, path):
        with open(path, "w") as f:
            json.dump(self.dict(), f, indent=4)

    @classmethod
    def load(cls, path):
        with open(path, "r") as f:
            return cls(**json.load(f))

    def __str__(self):
        return f"{self.Type} similarity between {self.Reference_Ligand} and {self.Query_Ligand}: {self.Tanimoto}"

    @classmethod
    def construct_dataframe(cls, similarity_list):
        return pd.DataFrame.from_records(
            [similarity.dict() for similarity in similarity_list]
        )


class ECFPSimilarity(MoleculeSimilarity):
    """
    ECFPSimilarity
    """

    Type: str = Field("ECFP", description="Type of similarity")
    radius: int = 2
    bitsize: int = 2048

    @property
    def Fingerprint(self):
        return f"ECFP{self.radius * 2}_{self.bitsize}"

    # add fingerprint property to output dictionary
    def dict(self, *args, **kwargs):
        return_dict = super().dict()
        return_dict.update({"fingerprint": self.Fingerprint})
        return return_dict


class MCSSimilarity(MoleculeSimilarity):
    """
    MCSSimilarity
    """

    Type: str = Field("MCS", description="Type of similarity")
    N_Atoms_in_MCS: int = Field(..., description="Number of atoms in the MCS")
    N_Atoms_in_Union: int = Field(
        ..., description="Total number of atoms between the two molecules"
    )

    @root_validator(pre=True)
    @classmethod
    def validate_tanimoto(cls, values):
        if "num_atoms_in_mcs" in values and "num_atoms_in_union" in values:
            expected_tanimoto = (
                values["num_atoms_in_mcs"] / values["num_atoms_in_union"]
            )
            assert (
                abs(values["tanimoto"] - expected_tanimoto) < 0.0001
            ), f"Tanimoto {values['tanimoto']} does not match expected value {expected_tanimoto}"
        return values


class TanimotoComboType(str, Enum):
    """
    Enum for the different types of Tanimoto coefficients that can be calculated.
    """

    SHAPE = "TanimotoShape"
    COLOR = "TanimotoColor"
    COMBO = "TanimotoCombo"


class TanimotoComboSimilarity(MoleculeSimilarity):
    """
    TanimotoComboSimilarity
    """

    Type: str = Field("TanimotoCombo", description="Type of similarity")
    Tanimoto_Shape: confloat(ge=0, le=1) = Field(..., description="Shape Tanimoto")
    Tanimoto_Color: confloat(ge=0, le=1) = Field(..., description="Color Tanimoto")
    Aligned: bool = False

    @classmethod
    def from_tanimoto_results(
        cls, ref, query, results, aligned: bool
    ) -> "TanimotoComboSimilarity":
        """
        This uses the OpenEye results class to get the TanimotoCombo similarity.
        Here we divide the TanimotoCombo by 2 to get it in the same range as all the others.
        :param ref:
        :param query:
        :param results:
        :param aligned:
        :return:
        """
        return cls(
            Reference_Ligand=ref,
            Query_Ligand=query,
            Tanimoto=results.GetTanimotoCombo() / 2,
            Tanimoto_Shape=results.GetShapeTanimoto(),
            Tanimoto_Color=results.GetColorTanimoto(),
            Aligned=aligned,
        )


# just test within this file
if __name__ == "__main__":
    ecfp = ECFPSimilarity(
        Reference_Ligand="mol1",
        Query_Ligand="mol2",
        Tanimoto=0.5,
    )
    ecfp.save("ecfp.json")
    ecfp2 = ECFPSimilarity.load("ecfp.json")
    assert ecfp == ecfp2

    mcs = MCSSimilarity(
        Reference_Ligand="mol1",
        Query_Ligand="mol2",
        N_Atoms_in_MCS=5,
        N_Atoms_in_Union=10,
        Tanimoto=0.5,
    )
    mcs.save("mcs.json")
    mcs2 = MCSSimilarity.load("mcs.json")

    try:
        mcs3 = MCSSimilarity(
            Reference_Ligand="mol1",
            Query_Ligand="mol2",
            N_Atoms_in_MCS=5,
            N_Atoms_in_Union=10,
            Tanimoto=0.6,
        )
        raise AssertionError("Tanimoto should not match")
    except ValidationError as e:
        print("Success")

    # test dataframe construction
    df = MoleculeSimilarity.construct_dataframe([ecfp, ecfp2, mcs, mcs2])
    df.to_csv("similarity.csv", index=False)
