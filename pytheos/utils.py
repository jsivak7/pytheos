# General utilities

from ase import Atoms
from pymatgen.core import Structure


def read_structure_to_ase_atoms(file_path: str) -> Atoms:
    """
    Read in structure file to ASE Atoms Object

    Args:
        file_path (str): relative path to structure file

    Returns:
        Atoms: ASE object to perform other operations
    """
    from ase import io

    s = io.read(f"{file_path}")
    print(f"{file_path} read in as ASE Atoms object")
    return s


def write_structure_from_ase_atoms(
    struc: Atoms,
    file_path: str,
    overwrite=False,
    sort=True,
):
    """
    Write ASE Atoms object to structure file
    NOTE that always writes "direct" coordinates

    Args:
        struc (Atoms): structure to be written
        file_path (str): relative path to write structure file, include suffix for desired file type (*.vasp, *.cif, etc.)
        overwrite (bool): override to write over file already made. Defaults to False.
        sort (bool): sort structure elements by electronegativity using Pymatgen. Defaults to True.

    Raises:
        FileExistsError: if given path for output file already exists, ensures that previously generated structures are not overwritten
    """
    import os
    from ase import io

    if os.path.exists(file_path) and overwrite == False:
        raise FileExistsError(file_path)

    if sort == True:
        print("Sorting structure by electronegativity.")
        from pymatgen.core import Structure

        struc_pmg = Structure.from_ase_atoms(struc)
        struc_pmg.sort()
        struc = struc_pmg.to_ase_atoms()

    io.write(f"{file_path}", struc, direct=True)
    print(f"ASE Atoms written to {file_path}")


def rattle_atoms(struc: Atoms, stddev=0.02) -> Atoms:
    """
    Rattles atoms of a given input structure - often is helpful prior to relaxation to break initial symmetry

    Args:
        struc (Atoms): structure to be rattled
        stddev (float, optional): standard deviation for amount of rattling to perform in Angstroms. Defaults to 0.02.

    Returns:
        Atoms: ASE Atoms object for rattled structure
    """
    import random

    struc.rattle(stddev, seed=int(random.uniform(0, 2000)))  # random seed

    return struc
