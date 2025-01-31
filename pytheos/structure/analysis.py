# for analyzing structures
# I think I have somewhat optimized versions of these somewhere...???
# TODO look around for these across different projects...


def get_firstNN_bonds(
    atom_num: int, struc_path="output.vasp", radius=3.00, anion="O", num_NNs=6
):
    """
    Gets all first nearest neighbor (NN) bond lengths for a specified cation with surrounding anions within a structure
    - a flexible scheme has been implemented that allows for consistent NN extraction even for highly distorted lattices

    Args:
    - atom_num (int): atom number of interest
    - struc_path (string): relative path to structure file
        - default = 'output.vasp'
    - radius (float): radius in Angstroms to search for NNs
        - default = 3.00
    - anion (string): anion element
        - default = 'O'
    - num_NNs (int): number of NNs to extract
        - default = 6 (octahedral)

    Returns:
        1d list of first near-neighbor bond lengths
    """

    from pymatgen.core.structure import Structure
    from pymatgen.core import Composition
    import numpy as np

    struc = Structure.from_file(struc_path)

    print(f"\natom number {atom_num} ({struc[atom_num].species})")

    # variables related to the scheme I have set up for always getting the desired number of NN bonds
    shift = 0.01  # Angstroms - initial guess of how much to shift radius and try to find desired number of NNs if no success initially
    counter = 0
    max_counter = 100  # max number of 'shifts' before decreasing shift value

    if struc[atom_num].species != Composition(anion):  # check to make sure that
        print("radius = {:.6f} Angstroms".format(radius))
        current_atom_neighbors = struc.get_all_neighbors_py(r=radius)[atom_num]
        indices = []
        distances = []

        for nn in current_atom_neighbors:
            if nn.species == Composition(anion):  # only want anion since 1st NN
                indices.append(nn.index)
                distances.append(nn.nn_distance)

        while len(indices) != num_NNs:
            if counter >= max_counter:
                counter = 0
            counter += 1

            if len(indices) < num_NNs:
                indices = []
                distances = []
                radius += shift
                print("radius = {:.6f} Angstroms".format(radius))
                current_atom_neighbors = struc.get_all_neighbors_py(r=radius)[atom_num]

                for nn in current_atom_neighbors:
                    if nn.species == Composition(anion):  # only want anion since 1st NN
                        indices.append(nn.index)
                        distances.append(nn.nn_distance)
                print("\t--> {} NNs".format(len(indices)))

                if (
                    counter >= max_counter
                ):  # only do the allowed amount of trials at each radius before making smaller
                    shift = shift * 0.05

            elif len(indices) > num_NNs:
                indices = []
                distances = []
                radius -= shift
                print("radius = {:.6f} Angstroms".format(radius))
                current_atom_neighbors = struc.get_all_neighbors_py(r=radius)[atom_num]

                for nn in current_atom_neighbors:
                    if nn.species == Composition(anion):  # only want anion since 1st NN
                        indices.append(nn.index)
                        distances.append(nn.nn_distance)
                print("\t--> {} NNs".format(len(indices)))

                if (
                    counter >= max_counter
                ):  # only do the allowed amount of trials at each radius before making smaller
                    shift = shift * 0.05

        if (
            len(indices) != num_NNs
        ):  # fail safe for if this still does not work to due to large deviations from octahedral coordination
            raise ValueError(
                "Number of 1st NN does not equal {}!\n\t--> {}".format(
                    num_NNs, len(indices)
                )
            )

    else:  # extra fail safe if the automatic scheme above does not work properly
        print("Current atom is an anion ({})!".format(anion))

    return distances


def get_diffraction_pattern(struc_path: str, scaled=True):
    """
    Gets simulated diffraction pattern for a structure file using Pymatgen

    Args:
        struc_path (str): relative path to structure file
        scaled (bool, optional): if intensities should be scaled so that max=100. Defaults to True.

    Returns:
        dict: 2theta values with corresponding intensities
    """
    from pymatgen.core.structure import Structure
    from pymatgen.analysis.diffraction.xrd import XRDCalculator

    struc = Structure.from_file(struc_path)
    calculator = XRDCalculator()
    pattern = calculator.get_pattern(struc, scaled=scaled)
    data = {"2theta": pattern.x, "intensity": pattern.y}
    return data


# TODO
def get_octahedral_bondangles(
    struc_path: str,
    bsite_cations: tuple,
    bondlength_max=2.5,
    bondangle_min=120,
    write_csv=False,
):
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


# TODO
def get_lattice_parameters(struc="output.vasp"):
    """
    Gets lattice parmeters for a CONTCAR structure

    Args:
        structure (string): structure file
            - default = 'CONTCAR'

    Returns:
        1d list of lattice parameters: [a, b, c]
    """

    from pymatgen.core.structure import Structure

    lattice_parameters = struc.lattice.abc
    print("a = {:.4f} \u212b".format(lattice_parameters[0]))
    print("b = {:.4f} \u212b".format(lattice_parameters[1]))
    print("c = {:.4f} \u212b".format(lattice_parameters[2]))
    return lattice_parameters
