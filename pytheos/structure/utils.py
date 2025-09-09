# utilities for manipulating structures

import os
import random
from ase import Atoms
from pymatgen.core import Structure


def read_structure(
    filename: str,
) -> Atoms:
    """
    Reads in structure file as ASE Atoms object.

    Args:
        filename (str): Relative path for structure file.

    Returns:
        Atoms: ASE Atoms object for structure.
    """

    from ase.io import read

    structure = read(filename=filename)

    return structure


def write_structure(
    structure: Atoms,
    output_filename: str,
    overwrite: bool = False,
    sort: bool = True,
) -> None:
    """
    Writes ASE Atoms object to file.

    Args:
        structure (Atoms): ASE Atoms object of structure.
        output_filename (str, optional): Output file name - include path and file type. e.g. "../material.vasp".
        overwrite (bool, optional): Overwrite over already made file. Defaults to False.
        sort (bool, optional): Option to sort elements by electronegativity. Defaults to True.

    Raises:
        FileExistsError: If supplied output_filename already exists to ensure that previously
            generated structures are not overwritten.

    Returns:
        None: Writes structure to file.
    """

    from ase.io import write

    if os.path.exists(output_filename) and overwrite == False:
        raise FileExistsError(output_filename)

    if sort == True:
        structure = sort_elements(structure)

    if "vasp" in output_filename or "poscar" in output_filename:
        write(output_filename, structure, direct=True)

    else:
        write(output_filename, structure)

    return None


def sort_elements(
    structure: Atoms,
) -> Atoms:
    """
    Sorts elements by electronegativity using Pymatgen - follows default for VASP POSCAR.

    Args:
        structure (Atoms): ASE Atoms object of structure.

    Returns:
        Atoms: ASE Atoms object with sorted structure.
    """

    structure_pmg = Structure.from_ase_atoms(structure)
    structure_pmg.sort()
    structure = structure_pmg.to_ase_atoms()

    return structure


def rattle_atoms(
    structure: Atoms,
    std_dev=0.01,
) -> Atoms:
    """
    Rattles atoms in structure using ASE - helpful prior to relaxation to break initial symmetry,
    which can enable symmetry-broken atomic arrangements during structure relaxation that are often
    lower energy structures.
    - see Zunger group's publications on exploring a 'polymorphous representation' for more
    details of this approach

    Args:
        structure (Atoms): ASE Atoms object of structure.
        std_dev (float, optional): standard deviation of rattling amount to perform in Angstroms.
            Defaults to 0.01.

    Returns:
        Atoms: Rattled structure.
    """

    structure.rattle(
        std_dev,
        seed=int(random.uniform(0, 2000)),  # random seed to ensure each time is unique
    )

    return structure
