# tools for analyzing structures

from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.analysis.diffraction.xrd import XRDCalculator
import numpy as np
from pymatgen.core import Composition
import pandas as pd
import time


def get_space_group(
    struc: Structure,
    symprec: float = 0.01,
    angle_tolerance: float = 5.0,
) -> tuple:
    """
    Get space group symbol and international space group number for inputted structure.

    The Pymatgen defaults are used as defaults here as well.

    Args:
        struc (Structure): structure to get space group
        symprec (float, optional): Tolerance for symmetry search. Defaults to 0.01.
        angle_tolerance (float, optional): Angle tolerance for symmetry search. Defaults to 5.0.

    Returns:
        dict: "symbol" key for space group symbol, and "number" key for international space group symbol
    """

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

    space_groups = {"symbol": symbol, "number": number}

    return space_groups


def simulate_diffraction_pattern(struc: Structure, scaled=True) -> dict:
    """
    Get simulated diffraction pattern for a structure file using Pymatgen

    Args:
        struc (Structure): structure to for generating diffraction pattern
        scaled (bool, optional): if intensities should be scaled so that max=100. Defaults to True.

    Returns:
        dict: 2theta values with corresponding intensities
    """

    calculator = XRDCalculator()

    pattern = calculator.get_pattern(struc, scaled=scaled)
    diffraction_data = {"2theta": pattern.x, "intensity": pattern.y}

    return diffraction_data


def extract_lattice_parameters(struc: Structure) -> tuple:
    """
    Gets lattice parameters.

    Args:
        struc (Structure): Structure input.

    Returns:
        tuple: (a, b, c)
    """

    lattice_parameters = struc.lattice.abc

    print(f"a = {np.round(lattice_parameters[0], 4)} \u212b")
    print(f"b = {np.round(lattice_parameters[1], 4)} \u212b")
    print(f"c = {np.round(lattice_parameters[2], 4)} \u212b")

    return lattice_parameters


def extract_firstNN_bonds(
    struc: Structure,
    num_NNs: int,
    radius=3.00,
    anion="O",
) -> dict:
    """
    Gets all first nearest neighbor (NN) cation-anion bond lengths for a specified structure.

    *Note that this has primarily been applied only for rocksalt systems (where all cation-oxygen bonds are equivalent).*
    I plan to implement another methodology for structure containing multiple inequivalent cation sites in the future...

    A flexible, self-consistent scheme has been implemented that allows for consistent NN extraction even for highly distorted lattices commonly present in HEOs.
    This is done by applying a 'shift' to the search radius around atom of interest by comparing the expected and actual number of NNs found.

    Args:
        struc (Structure): Pymatgen Structure input.
        num_NNs (int): Number of NN to extract.
        radius (float, optional): Starting radius for NN search in Angstroms. Defaults to 3.00.
        anion (str, optional): Anion element. Defaults to "O".

    Raises:
        ValueError: If implemented self-consistent scheme still cannot find correct number of NNs.

    Returns:
        dict: {"distances": 1D list bondlengths, "indices": 1D list anion indices corresponding to bond lengths}
    """

    all_data = {
        "species": [],
        "distance": [],
        "anion index": [],
    }

    # loop through all atoms within structure
    for atom_num in range(int(len(struc))):

        cation_species = struc[atom_num].species.chemical_system
        print(f"\natom #{atom_num} ({cation_species})")

        # variables related to the self-consistent scheme I have set up to extract the desired number of NN bonds for highly distorted structures
        shift = 0.01  # Angstroms - initial guess of how much to 'shift' radius and try to find desired number of NNs if no success initially
        counter = 0  # some reasonable default
        max_counter = 100  # max number of 'shifts' before changing shift value

        if struc[atom_num].species != Composition(anion):  # ensure cation
            print(f" radius = {np.round(radius, 6)} Å")
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

                # attempt if finds less than desired number of NNs -> increase search radius by shift value
                if len(indices) < num_NNs:
                    radius += shift

                # attempt if finds more than desired number of NNs -> decrease search radius by shift value
                elif len(indices) > num_NNs:
                    radius -= shift

                indices = []
                distances = []

                print(f" radius = {np.round(radius, 6)} Å")
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

            print(f" bondlengths = {np.round(distances, 3)} Å")

            all_data["species"].extend([cation_species] * 6)
            all_data["distance"].extend(distances)
            all_data["anion index"].extend(indices)

    print(
        f"A total of {len(all_data['distance'])} cation-anion first NN bonds were found."
    )

    return pd.DataFrame.from_dict(all_data)


def extract_octahedral_bondangles(
    struc: Structure,
    bsite_cations: list,
    max_bondlength=2.5,
    min_bondangle=120,
) -> pd.DataFrame:
    """
    Given a structure file, extracts the B-O-B bond angles.
    - This is only tested and usually used for a perovskite structure, however this should (at least in theory) work for other structures...
    - NOTE currently does NOT take into account PBCs, thus you may want to utilize a 2x2x2 supercell of your relaxed calculation...

    Args:
        struc (Structure): Pymatgen Structure object.
        bsite_cations (tuple): b-site cations that will be searched over. Example: ["Ti"], or ["Ti", "V"]
        max_bondlength (float, optional): maximum bond length allowed between atom1-atom2 & atom2-atom3. Defaults to 2.5 for a reasonable default.
        min_bondangle (int, optional): minimum bond angle allowed. Defaults to 120 for a reasonable value.

    Returns:
        Pandas DataFrame: bond lengths for B-O-B combinations including atom numbers and species
    """

    start_time = time.time()

    print(f"Extracting average octahedral bond angles...")
    print(f"bsite cations: {bsite_cations}")

    len_struc = np.arange(0, len(struc))
    counter = 0
    tracker_list = []
    atoms1 = []
    atoms1_species = []
    atoms2 = []
    atoms2_species = []
    atoms3 = []
    atoms3_species = []
    bond_angles = []
    bsite_compositions = []

    for bsite_cation in bsite_cations:
        print(bsite_cation)
        bsite_compositions.append(Composition(bsite_cation))

    # loop through all combinations of 3 atoms within structure
    for atom1 in len_struc:
        if struc[atom1].species in bsite_compositions:
            for atom2 in len_struc:
                if atom1 != atom2:
                    if (
                        struc[atom2].species == Composition("O")
                        and struc.get_distance(atom1, atom2) < max_bondlength
                    ):
                        for atom3 in len_struc:
                            if (
                                struc[atom3].species in bsite_compositions
                                and struc.get_distance(atom2, atom3) < max_bondlength
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
                                            > min_bondangle
                                        ):
                                            print(
                                                f"\t{struc[atom1].species}(#{atom1}) - {struc[atom2].species}(#{atom2}) - {struc[atom3].species}(#{atom3}) -> {struc.get_angle(atom1, atom2, atom3):.1f}\u00b0"
                                            )

                                            # all checks passed so moving forward with getting information
                                            tracker_list.append([atom1, atom2, atom3])
                                            counter += 1

                                            atoms1.append(atom1)
                                            atoms1_species.append(
                                                struc[atom1].species.reduced_formula
                                            )
                                            atoms2.append(atom2)
                                            atoms2_species.append(
                                                struc[atom2].species.reduced_formula
                                            )
                                            atoms3.append(atom3)
                                            atoms3_species.append(
                                                struc[atom3].species.reduced_formula
                                            )
                                            bond_angles.append(
                                                np.round(
                                                    struc.get_angle(
                                                        atom1, atom2, atom3
                                                    ),
                                                    3,
                                                )
                                            )
    end_time = time.time()

    print(f"num bond angles found: {counter}")
    avg_bondangle = np.average(bond_angles)

    data = pd.DataFrame(
        {
            "atom1": atoms1,
            "atom1_species": atoms1_species,
            "atom2": atoms2,
            "atom2_species": atoms2_species,
            "atom3": atoms3,
            "atom3_species": atoms3_species,
            "bondangle": bond_angles,
        }
    )

    print(f"avg bond angle: {avg_bondangle:.1f}\u00b0")
    print(f"elapsed time: {end_time-start_time:.2f} seconds")
    return data
