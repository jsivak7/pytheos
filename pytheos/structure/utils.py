# utilities for manipulating structures

import random
from pymatgen.core import Structure


def rattle_atoms(
    struc: Structure,
    std_dev=0.01,
) -> Structure:
    """
    Rattles atoms in structure using ASE - helpful prior to relaxation to break initial symmetry,
    which can enable symmetry-broken atomic arrangements during structure relaxation that are often
    lower energy structures.
    - see Zunger group's publications on exploring a 'polymorphous representation' for more
    details of this approach

    Args:
        structure (Structure): Pymatgen Structure object.
        std_dev (float, optional): standard deviation of rattling amount to perform in Angstroms.
            Defaults to 0.01.

    Returns:
        Structure: Rattled structure as Pymatgen Structure object.
    """

    struc = struc.to_ase_atoms()

    struc.rattle(
        std_dev,
        seed=int(random.uniform(0, 2000)),  # random seed to ensure each time is unique
    )

    struc = Structure.from_ase_atoms(struc)

    return structure
