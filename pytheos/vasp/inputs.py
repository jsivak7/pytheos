# For creating VASP input files and directories

from ase import Atoms


def write_relax_inputs(
    struc: Atoms,
    output_dir: str,
    user_incar_changes=None,
    functional="r2scan",
):
    """
    Makes VASP input files (INCAR, KPOINTS, POSCAR, POTCAR) using Pymatgen for a relaxation calculation

    Args:
        struc (Atoms): structure to be relaxed
        output_dir (str): relative path to output generated files
        user_incar_changes (dict, optional): changes that deviate from default INCAR parameters given in *.yaml files. Defaults to None.
        functional (str, optional): XC functional to be used. Defaults to "r2scan".

    Raises:
        FileExistsError: if output_dir already exists
        ValueError: if an invalid functional is given
    """

    import os
    import yaml
    from pymatgen.core import Structure
    from pymatgen.io.vasp.sets import MPRelaxSet, MPScanRelaxSet
    from pymatgen.io.vasp.inputs import Kpoints

    if os.path.exists(output_dir):
        raise FileExistsError(output_dir)
    else:
        os.mkdir(output_dir)

    print(f"Writing VASP input files to {output_dir}")

    module_dir = os.path.dirname(__file__)
    s = Structure.from_ase_atoms(struc)

    with open(f"{module_dir}/{functional}.yaml", "r") as f:
        incar_settings = yaml.load(f, Loader=yaml.SafeLoader)

    if bool(user_incar_changes) != False:
        incar_settings.update(user_incar_changes)

    if functional == "r2scan":
        calc = MPScanRelaxSet(
            structure=s,
            user_incar_settings=incar_settings,
            user_potcar_functional="PBE_54",
            sort_structure=True,
        )
    elif functional == "pbe":
        calc = MPRelaxSet(
            structure=s,
            user_incar_settings=incar_settings,
            user_potcar_functional="PBE",
            sort_structure=True,
        )
    else:
        raise ValueError(f"{functional} is not available.")

    calc.write_input(f"{output_dir}")
