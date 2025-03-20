# facilitate Vienna Ab initio Simulation Package (VASP) calculation inputs
# https://www.vasp.at/wiki/index.php/The_VASP_Manual

# NOTE that it is assumed here that the KSPACING INCAR tag is used rather than a KPOINTS file
# https://www.vasp.at/wiki/index.php/KSPACING

from ase import Atoms
from pymatgen.io.vasp.inputs import Incar, Poscar, Kpoints, Potcar


def write_inputs(
    struc: Atoms,
    output_dir: str = "relaxation",
    user_incar_changes=None,
    functional="r2scan",
) -> None:
    """
    Makes VASP input files from a structure as an ASE Atoms object.
    Usually is a relaxation, however can handle other calculation types if desired.

    Args:
        struc (Atoms): structure to be relaxed
        output_dir (str): relative path to output generated files. Defaults to "relaxation".
        user_incar_changes (dict, optional): changes that deviate from default
            INCAR parameters given in pbe/r2scan.yaml files. Defaults to None.
        functional (str, optional): Exchange-Correlation (XC) functional to be used.
            Defaults to "r2scan". Options can be found in /pytheos/vasp/inputs_sets/.

    Raises:
        FileExistsError: if output_dir already exists
        ValueError: if an invalid functional is given

    Returns:
        None: New directory is made for a VASP calculation -> ./output_dir
    """

    import os
    import yaml
    from pymatgen.core import Structure
    from pymatgen.io.vasp.sets import MPRelaxSet, MPScanRelaxSet
    from pymatgen.io.vasp.inputs import Kpoints

    print(f"Writing general VASP calculation.")

    os.mkdir(output_dir)

    # read in structure file as ASE Atoms object
    s = Structure.from_ase_atoms(struc)

    # get specific functional (xc.yaml) file for incar settings
    module_dir = os.path.dirname(__file__)  # get path to module
    with open(f"{module_dir}/input_sets/{functional}.yaml", "r") as f:
        incar_settings = yaml.load(f, Loader=yaml.SafeLoader)

    # only make further changes if user gives them
    if bool(user_incar_changes) == True:
        incar_settings.update(user_incar_changes)

    if functional in ["pbe", "pbesol"]:
        calc = MPRelaxSet(
            structure=s,
            user_incar_settings=incar_settings,
            user_potcar_functional="PBE",
            sort_structure=True,
        )
    elif functional in ["r2scan", "scan"]:
        calc = MPScanRelaxSet(
            structure=s,
            user_incar_settings=incar_settings,
            user_potcar_functional="PBE_54",
            sort_structure=True,
        )
    else:
        raise ValueError(
            f"{functional} is not available. Choose either a GGA ('pbe', 'pbesol') or a MetaGGA ('r2scan', 'scan') functional."
        )

    calc.write_input(f"{output_dir}")

    return None


def get_magorder(
    struc_path: str,
    magorder_name: str,
    magmom_values: dict,
    rattle_amount: float = 0,
    round_coords: int = 0,
) -> list:
    """
    Gets the magnetic ordering for an arbitrary structure.

    NOTE that a *.yaml file must exist for the given magorder_name arguement in pytheos/vasp/mag_orders/

    Magnetic moments can be "rattled" some amount (arg: rattle_amount) around the specified magnetic moment (arg: magmom_values).
    This can be beneficial to break the initial magnetic symmetry and assists in electronic convergence.
    - random value between magmom_value +/- rattle_amount

    Args:
        struc_path (str): relative path to structure file
        magorder_name (str): magnetic ordering name to apply to structure
        magmom_values (dict): magnetic moment values (in Bohr-Magneton) for different elements. Example -> {"Ni": 3, "O": 0}
        rattle_amount (float, optional): amount to distort initial magnetic moment. Defaults to 0.
        round_coords (int, optional): number of decimals to round atom positions - sometimes helpful for distorted structures. Defaults to 0.

    Raises:
        ValueError: If atom position cannot be identified as spin-up or spin-down.

    Returns:
        list: Initial magnetic moment values to be passed to VASP input set
    """

    from pytheos.structure import utils
    import numpy as np
    import os
    import yaml
    import random

    print(f"----- Getting AFM magnetic ordering -----")

    # load afm_order.yaml per "magorder_name" argument
    module_path = os.path.dirname(__file__)
    with open(f"{module_path}/mag_orders/{magorder_name}.yaml", "r") as f:
        spins = yaml.load(f, Loader=yaml.SafeLoader)

    # split up spin up and spin down positions
    spin_up_positions = np.round(np.array(spins["spin_up"]), round_coords)
    spin_down_positions = np.round(np.array(spins["spin_down"]), round_coords)

    s = utils.read_to_ase_atoms(file_path=struc_path)
    s = utils.convert_ase_atoms_to_pmg_structure(s)

    magmoms = []
    counts = {"up": 0, "down": 0}

    atom_index = 0
    for atom in s:
        spin_label = None
        rounded_coords = np.round(atom.coords, round_coords)

        try:
            magmom_values[atom.label]
        except KeyError:
            print(f"{atom.label} not found in magmom_values dictionary.")

        magmom_value = magmom_values[atom.label]
        magmom = np.round(
            random.uniform(
                magmom_value - rattle_amount,
                magmom_value + rattle_amount,
            ),
            2,
        )

        if atom.label == "O":
            spin_label = None  # only care about cations

        else:  # cations
            position = np.round(atom.coords, round_coords)
            identified_position = False

            for spin_up_position in spin_up_positions:
                if np.array_equal(position, spin_up_position):
                    identified_position = True
                    spin_label = "spin up"
                    counts["up"] += 1
                    magmom = +magmom

            for spin_down_position in spin_down_positions:
                if np.array_equal(position, spin_down_position):
                    identified_position = True
                    spin_label = "spin down"
                    counts["down"] += 1
                    magmom = -magmom

            if identified_position == False:
                raise ValueError(
                    f"Spin not properly identified for atom #{atom_index} ({atom.label}) @ {rounded_coords}"
                )

        magmoms.append(magmom)
        print(
            f"#{atom_index:<4} {atom.label:<4} {str(rounded_coords):<20}  -> {str(spin_label):<12} {str(magmom):<4}"
        )
        atom_index += 1

    if counts["up"] != counts["down"]:
        question = input(
            "WARNING! The number of atoms identified as spin-up and spin-down are not the same - do you want to proceed? (y/n)\n"
        )
        if question.lower() in ["yes", "y"]:
            return magmoms
        else:
            print("ABORTING!")
    else:
        return magmoms


def update_incar_with_magorder(path: str, magmom_list: list) -> None:
    """
    Updates INCAR file with new MAGMOM magnetic ordering.

    Args:
        path (str): relative path to INCAR file.
        magmom_list (list): magnetic moments that will update MAGMOM in INCAR file

    Returns:
        None
    """
    from pymatgen.io.vasp.inputs import Incar

    print(f"----- Updating INCAR with magnetic ordering -----")

    i = Incar.from_file(filename=path)
    i["MAGMOM"] = magmom_list
    i.write_file(filename=path)

    return None


def load_incar(path: str = "INCAR") -> Incar:
    """
    Loads INCAR file.

    Args:
        path (str, optional): Relative path for INCAR file. Defaults to "INCAR".

    Returns:
        Incar: Pymatgen Incar object
    """

    i = Incar.from_file(filename=path)

    return i


def load_poscar(path: str = "POSCAR") -> Poscar:
    """
    Loads POSCAR file.

    Args:
        path (str, optional): Relative path for POSCAR file. Defaults to "POSCAR".

    Returns:
        Poscar: Pymatgen Poscar object
    """

    p = Poscar.from_file(filename=path)

    return p


def load_potcar(path: str = "POTCAR") -> Potcar:
    """
    Loads POTCAR file.

    Args:
        path (str, optional): Relative path for POTCAR file. Defaults to "POTCAR".

    Returns:
        Potcar: Pymatgen Potcar object
    """

    pot = Potcar.from_file(filename=path)

    return pot


def load_kpoints(path: str = "KPOINTS") -> Kpoints:
    """
    Loads KPOINTS file.

    Args:
        path (str, optional): Relative path for KPOINTS file. Defaults to "KPOINTS".

    Returns:
        Kpoints: Pymatgen Kpoints object
    """

    kpts = Kpoints.from_file(filename=path)

    return kpts
