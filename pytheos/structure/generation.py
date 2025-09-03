# tools for generating structures

from ase import Atoms
from icet import ClusterSpace
from icet.tools.structure_generation import (
    generate_sqs_from_supercells,
    _get_sqs_cluster_vector,
    occupy_structure_randomly,
)
from icet.input_output.logging_tools import set_log_config


def make_supercell(
    structure: Atoms,
    dimensions: list,
) -> Atoms:
    """
    Makes a supercell from structure using ASE.

    Args:
        structure (Atoms): ASE Atoms object of structure.
        dimensions (list): [x, y, z] cell multipliers for supercell generation.

    Returns:
        Atoms: Supercell structure.
    """

    supercell = structure.repeat(dimensions)

    return supercell


def make_sqs(
    structure: Atoms,
    dimensions: list,
    chemical_symbols,
    concentrations: dict,
    cutoffs: list,
    num_mc_steps: int = 10000,
) -> Atoms:
    """
    Generates a special quasirandom structure (SQS) using ICET.
    - see https://gitlab.com/materials-modeling/icet

    Args:
        structure (Atoms): ASE Atoms object of structure.
        dimensions (list): [x, y, z] cell multipliers for supercell generation.
        chemical_symbols (list): list of lists for allowed elements following same order as structure.
            e.g. [["Sr"], ["Ti", "Cr"], ["O"], ["O"], ["O"]] for a perovskite.
        concentrations (dict): Fractions of different elements for each lattice site.
            Only need to specify those that are not 1.
        cutoffs (list): cutoffs in order of multiplicity (pair, triplet, quadruplet, ...).
        num_mc_steps (int, optional): Number of Monte Carlo steps to run the SQS generation. Defaults to 10000.

    Returns:
        Atoms: SQS structure.
    """

    set_log_config(level="INFO")

    cluster_space = ClusterSpace(
        structure=structure,
        cutoffs=cutoffs,
        chemical_symbols=chemical_symbols,
    )
    print(cluster_space)

    sqs = generate_sqs_from_supercells(
        cluster_space=cluster_space,
        supercells=[structure.repeat(dimensions)],
        target_concentrations=concentrations,
        n_steps=num_mc_steps,
    )

    trial_cluster_vector = cluster_space.get_cluster_vector(sqs)
    perfectly_random_cluster_vector = _get_sqs_cluster_vector(
        cluster_space=cluster_space,
        target_concentrations=concentrations,
    )

    print(f"\nTrial Cluster Vector ->\n{trial_cluster_vector}")
    print(f"\nPerfectly Random Cluster Vector ->\n{perfectly_random_cluster_vector}")

    return sqs


def decorate_randomly(
    structure: Atoms,
    dimensions: list,
    chemical_symbols: list,
    concentrations: dict,
) -> Atoms:
    """
    Randomly decorates structure using ICET.
    - see https://gitlab.com/materials-modeling/icet

    Args:
        structure (Atoms): ASE Atoms object of structure.
        dimensions (list): [x, y, z] cell multipliers for supercell generation.
        chemical_symbols (_type_): list of lists for allowed elements following same order as structure.
            e.g. [["Sr"], ["Ti", "Cr"], ["O"], ["O"], ["O"]] for a perovskite.
        concentrations (dict): Fractions of different elements for each lattice site.
            Only need to specify those that are not 1.

    Returns:
        Atoms: Randomly decorated structure.
    """

    cluster_space = ClusterSpace(
        structure=structure,
        cutoffs=[0],  # just need something for cluster space construction
        chemical_symbols=chemical_symbols,
    )

    randomly_decorated_structure = structure.repeat(dimensions)

    occupy_structure_randomly(
        structure=randomly_decorated_structure,
        cluster_space=cluster_space,
        target_concentrations=concentrations,
    )

    return randomly_decorated_structure
