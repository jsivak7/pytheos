# tools for generating structures

from pymatgen.core import Structure
from icet import ClusterSpace
from icet.tools.structure_generation import (
    generate_sqs_from_supercells,
    _get_sqs_cluster_vector,
    occupy_structure_randomly,
)
from icet.input_output.logging_tools import set_log_config


def make_supercell(
    struc: Structure,
    dimensions: list,
) -> Structure:
    """
    Makes a supercell from structure using ASE.

    Automatically sorts structure again using Pymatgen to match VASP POSCAR format.

    Args:
        struc (Structure): Pymatgen Structure object.
        dimensions (list): [x, y, z] cell multipliers for supercell generation.

    Returns:
        Structure: Supercell structure.
    """

    struc = struc.to_ase_atoms()
    supercell = struc.repeat(dimensions)
    supercell = Structure.from_ase_atoms(supercell)
    supercell.sort()

    return supercell


def make_sqs(
    struc: Structure,
    dimensions: list,
    chemical_symbols,
    concentrations: dict,
    cutoffs: list,
    num_mc_steps: int = 10000,
) -> Structure:
    """
    Generates a special quasirandom structure (SQS) using ICET.
    - see https://gitlab.com/materials-modeling/icet

    Args:
        struc (Structure): Pymatgen Structure object.
        dimensions (list): [x, y, z] cell multipliers for supercell generation.
        chemical_symbols (list): list of lists for allowed elements following same order as structure.
            e.g. [["Sr"], ["Ti", "Cr"], ["O"], ["O"], ["O"]] for a perovskite.
        concentrations (dict): Fractions of different elements for each lattice site.
            Only need to specify those that are not 1.
        cutoffs (list): cutoffs in order of multiplicity (pair, triplet, quadruplet, ...).
        num_mc_steps (int, optional): Number of Monte Carlo steps to run the SQS generation. Defaults to 10000.

    Returns:
        Structure: SQS structure.
    """

    struc = struc.to_ase_atoms()

    set_log_config(level="INFO")

    cluster_space = ClusterSpace(
        structure=struc,
        cutoffs=cutoffs,
        chemical_symbols=chemical_symbols,
    )
    print(cluster_space)

    sqs = generate_sqs_from_supercells(
        cluster_space=cluster_space,
        supercells=[struc.repeat(dimensions)],
        target_concentrations=concentrations,
        n_steps=num_mc_steps,
    )

    trial_cluster_vector = cluster_space.get_cluster_vector(sqs)
    perfectly_random_cluster_vector = _get_sqs_cluster_vector(
        cluster_space=cluster_space,
        target_concentrations=concentrations,
    )

    sqs = Structure.from_ase_atoms(sqs)

    print(f"\nTrial Cluster Vector ->\n{trial_cluster_vector}")
    print(f"\nPerfectly Random Cluster Vector ->\n{perfectly_random_cluster_vector}")

    return sqs


def decorate_randomly(
    struc: Structure,
    dimensions: list,
    chemical_symbols: list,
    concentrations: dict,
) -> Structure:
    """
    Randomly decorates structure using ICET.
    - see https://gitlab.com/materials-modeling/icet

    Args:
        struc (Structure): Pymatgen Structure object.
        dimensions (list): [x, y, z] cell multipliers for supercell generation.
        chemical_symbols (_type_): list of lists for allowed elements following same order as structure.
            e.g. [["Sr"], ["Ti", "Cr"], ["O"], ["O"], ["O"]] for a perovskite.
        concentrations (dict): Fractions of different elements for each lattice site.
            Only need to specify those that are not 1.

    Returns:
        Structure: Randomly decorated structure.
    """

    struc = struc.to_ase_atoms()

    cluster_space = ClusterSpace(
        structure=struc,
        cutoffs=[0],  # just need something for cluster space construction
        chemical_symbols=chemical_symbols,
    )

    randomly_decorated_structure = struc.repeat(dimensions)

    occupy_structure_randomly(
        structure=randomly_decorated_structure,
        cluster_space=cluster_space,
        target_concentrations=concentrations,
    )

    randomly_decorated_structure = Structure.from_ase_atoms(
        randomly_decorated_structure
    )

    return randomly_decorated_structure
