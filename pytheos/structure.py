# structure tools

import os
import random
from ase import Atoms
from pymatgen.core import Structure
from icet import ClusterSpace
from icet.tools.structure_generation import (
    generate_sqs_from_supercells,
    _get_sqs_cluster_vector,
    occupy_structure_randomly,
)
from icet.input_output.logging_tools import set_log_config


def read_structure(
    filename: str,
):
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
        sort_elements (bool, optional): Option to sort elements by electronegativity. Defaults to True.

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
    Rattles atoms in structure using ASE - helpful prior to relaxation to break initial symmetry.

    Can enable symmetry-broken atomic arrangements during structure relaxation.
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
        seed=int(random.uniform(0, 2000)),  # random seed
    )

    return structure


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

    structure = supercell

    return structure


def make_sqs(
    structure: Atoms,
    dimensions: list,
    chemical_symbols,
    concentrations: dict,
    cutoffs: list,
    num_mc_steps: int = 10000,
) -> Atoms:
    """
    Generates a special quasi-random structure (SQS) using ICET.
    - https://gitlab.com/materials-modeling/icet

    Args:
        structure (Atoms): ASE Atoms object of structure.
        dimensions (list): [x, y, z] cell multipliers for supercell generation.
        chemical_symbols (_type_): list of lists for allowed elements following same order as structure.
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
    - https://gitlab.com/materials-modeling/icet

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


### analysis tools - subject to change since might make a StructureAnalyzer class...


def get_space_group(
    struc: Atoms,
    symprec: float = 0.01,
    angle_tolerance: float = 5.0,
) -> tuple:
    """
    Get space group symbol and international space group symbol for inputted structure.

    The Pymatgen defaults are used as defaults here as well.

    Args:
        struc (Atoms): structure to get space group
        symprec (float, optional): Tolerance for symmetry search. Defaults to 0.01.
        angle_tolerance (float, optional): Angle tolerance for symmetry search. Defaults to 5.0.

    Returns:
        tuple: (space group symbol, international space group number)
    """
    from pymatgen.core import Structure
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    struc = Structure.from_ase_atoms(struc)
    symbol = SpacegroupAnalyzer(
        structure=struc,
        symprec=symprec,
        angle_tolerance=angle_tolerance,
    ).get_space_group_symbol()

    number = SpacegroupAnalyzer(
        structure=struc,
        symprec=symprec,
        angle_tolerance=angle_tolerance,
    ).get_space_group_number()

    return (symbol, number)


def get_diffraction_pattern(struc: Atoms, scaled=True) -> dict:
    """
    Get simulated diffraction pattern for a structure file using Pymatgen

    Args:
        struc (Atoms): structure to for generating diffraction pattern
        scaled (bool, optional): if intensities should be scaled so that max=100. Defaults to True.

    Returns:
        dict: 2theta values with corresponding intensities
    """
    from pymatgen.core.structure import Structure
    from pymatgen.analysis.diffraction.xrd import XRDCalculator

    struc = Structure.from_ase_atoms(struc)
    calculator = XRDCalculator()
    pattern = calculator.get_pattern(struc, scaled=scaled)
    diffraction_data = {"2theta": pattern.x, "intensity": pattern.y}
    return diffraction_data


def get_lattice_parameters(struc: Atoms) -> tuple:
    """
    Gets lattice parmeters from an ASE Atoms object.

    Args:
        struc (Atoms): Structure input.

    Returns:
        tuple: (a, b, c)
    """

    from pymatgen.core.structure import Structure
    import numpy as np

    struc = Structure.from_ase_atoms(struc)
    lattice_parameters = struc.lattice.abc
    print(f"a = {np.round(lattice_parameters[0], 4)} \u212b")
    print(f"b = {np.round(lattice_parameters[1], 4)} \u212b")
    print(f"c = {np.round(lattice_parameters[2], 4)} \u212b")
    return lattice_parameters


def get_firstNN_bonds(
    struc: Atoms,
    atom_num: int,
    num_NNs: int,
    radius=3.00,
    anion="O",
) -> tuple:
    """
    Gets all first nearest neighbor (NN) bond lengths for a specified cation with surrounding anions within a structure.

    A flexible, self-consistent scheme has been implemented that allows for consistent NN extraction even for highly distorted lattices commonly present in HEOs.
    This is done by applying a 'shift' to the search radius around atom of interest by comparing the expected and actual number of NNs found.

    Args:
        struc (Atoms): ASE Atoms object as structure input.
        atom_num (int): Atom of interest.
        num_NNs (int): Number of NN to extract.
        radius (float, optional): Starting radius for NN search in Angstroms. Defaults to 3.00.
        anion (str, optional): Anion element. Defaults to "O".

    Raises:
        ValueError: If implemented self-consistent scheme still cannot find correct number of NNs.

    Returns:
        tuple: (1D list bondlengths, 1D list anion indices corresponding to bond lengths)
    """

    from pymatgen.core.structure import Structure
    from pymatgen.core import Composition
    import numpy as np

    struc = Structure.from_ase_atoms(struc)

    print(f"\natom #{atom_num} ({struc[atom_num].species})")

    # variables related to the self-consistent scheme I have set up to extract the desired number of NN bonds for highly distorted structures
    shift = 0.01  # Angstroms - initial guess of how much to 'shift' radius and try to find desired number of NNs if no success initially
    counter = 0  # some reasonable default
    max_counter = 100  # max number of 'shifts' before changing shift value

    if struc[atom_num].species != Composition(anion):  # ensure cation
        print(f"radius = {np.round(radius, 6)} Å")
        current_atom_neighbors = struc.get_all_neighbors_py(r=radius)[atom_num]
        indices = []
        distances = []

        for nn in current_atom_neighbors:
            if nn.species == Composition(anion):  # only want anions as first NNs
                indices.append(nn.index)
                distances.append(nn.nn_distance)

        while len(indices) != num_NNs:
            if counter >= max_counter:
                counter = 0
            counter += 1

            # attempt finds less than desired number of NNs -> increase search radius by shift value
            if len(indices) < num_NNs:
                radius += shift

            # attempt finds more than desired number of NNs -> decrease search radius by shift value
            elif len(indices) > num_NNs:
                radius -= shift

            indices = []
            distances = []

            print(f"radius = {np.round(radius, 6)} Å")
            current_atom_neighbors = struc.get_all_neighbors_py(r=radius)[atom_num]

            for nn in current_atom_neighbors:
                if nn.species == Composition(anion):  # only want anion since 1st NN
                    indices.append(nn.index)
                    distances.append(nn.nn_distance)
            print("\t--> {} NNs".format(len(indices)))

            # only do the allowed amount of trials at each radius before making smaller
            if counter >= max_counter:
                # progressively make |shift| smaller for finer resolution
                shift = shift * 0.05

        # fail safe for if the implemented scheme still does not work to due to large deviations expected coordination
        if len(indices) != num_NNs:
            raise ValueError(
                "Number of 1st NN does not equal {}!\n\t--> {}".format(
                    num_NNs, len(indices)
                )
            )

    else:  # should be anion - printing for clarity
        print("Current atom is an anion ({})!".format(anion))

    return (distances, indices)


# TODO NEEDS TO BE FIXED I THINK - CAREFUL USING YET
# TODO I BELIEVE I HAVE AN OPTIMIZED VERSION OF THIS FUNCTION IN PEROVSKITE CONDUCTORS PROJECT FILES
def get_octahedral_bondangles(
    struc_path: str,
    bsite_cations: tuple,
    bondlength_max=2.5,
    bondangle_min=120,
    write_csv=False,
) -> list:
    """
    Given a structure file, extracts the B-O-B bond angles.
    - NOTE that function currently does not take into account PBCs, thus a 2x2x2 supercell is used

    Args:
        struc_path (str): relative path to structure file
        bsite_cations (tuple): b-site cations that will be searched over
        bondlength_max (float, optional): maximum bond length allowed between atom1-atom2 & atom2-atom3. Defaults to 2.5 for a reasonable default.
        bondangle_min (int, optional): minimum bond angle allowed. Defaults to 120 for a reasonable value.
        write_csv (bool, optional): if user wants to write a .csv file with data on bond angles. Defaults to False.
            - "bondangles.csv"

    Returns:
        float: average bond angle
    """
    from pymatgen.core.structure import Structure
    from pymatgen.core.composition import Composition
    import numpy as np
    import pandas as pd
    import time

    start_time = time.time()

    print(f"\n========== Getting average octahedral bond angle ==========")
    print(f"structure path:\t\t{struc_path}")
    print(f"bsite cations:\t\t{bsite_cations}")

    struc = Structure.from_file(struc_path)
    struc.make_supercell((2, 2, 2))  # since PBC are not considered

    len_struc = np.arange(0, len(struc))
    counter = 0
    tracker_list = []
    atoms1 = []
    atoms2 = []
    atoms3 = []
    bond_angles = []
    bsite_compositions = []

    for bsite_cation in bsite_cations:
        bsite_compositions.append(Composition(bsite_cation))

    # loop through all combinations of 3 atoms within structure
    for atom1 in len_struc:
        if struc[atom1].species in bsite_compositions:
            for atom2 in len_struc:
                if atom1 != atom2:
                    if (
                        struc[atom2].species == Composition("O")
                        and struc.get_distance(atom1, atom2) < bondlength_max
                    ):
                        for atom3 in len_struc:
                            if (
                                struc[atom3].species in bsite_compositions
                                and struc.get_distance(atom2, atom3) < bondlength_max
                            ):
                                if (
                                    atom1 != atom3
                                    and atom2 != atom3
                                    and struc[atom3].species in bsite_compositions
                                ):
                                    if [atom1, atom2, atom3] not in tracker_list and [
                                        atom3,
                                        atom2,
                                        atom1,
                                    ] not in tracker_list:
                                        if (
                                            struc.get_angle(atom1, atom2, atom3)
                                            > bondangle_min
                                        ):
                                            print(
                                                f"\t{struc[atom1].species}(#{atom1}) - {struc[atom2].species}(#{atom2}) - {struc[atom3].species}(#{atom3}) -> {struc.get_angle(atom1, atom2, atom3):.1f}\u00b0"
                                            )

                                            # all checks passed so moving forward with getting information
                                            tracker_list.append([atom1, atom2, atom3])
                                            counter += 1

                                            atoms1.append(atom1)
                                            atoms2.append(atom2)
                                            atoms3.append(atom3)
                                            bond_angles.append(
                                                struc.get_angle(atom1, atom2, atom3)
                                            )
    end_time = time.time()

    print(f"num bond angles found:\t{counter}")
    avg_bondangle = np.average(bond_angles)

    if write_csv == True:
        data = pd.DataFrame(
            {"atom1": atoms1, "atom2": atom2, "atom3": atoms3, "bondangle": bond_angles}
        )

        data.to_csv("bondangles.csv", index=False)

    print(f"avg bond angle:\t\t{avg_bondangle:.1f}\u00b0")
    print(f"elapsed time:\t\t{end_time-start_time:.2f} seconds")
    return avg_bondangle
