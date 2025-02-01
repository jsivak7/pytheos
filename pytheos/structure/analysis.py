# for analyzing structures

from ase import Atoms


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
) -> list:
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
        list: 1D list of first NN bond lengths
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
                indices = []
                distances = []
                radius += shift
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

            # attempt finds more than desired number of NNs -> decrease search radius by shift value
            elif len(indices) > num_NNs:
                indices = []
                distances = []
                radius -= shift
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

    return distances


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
