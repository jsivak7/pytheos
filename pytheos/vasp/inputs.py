# module to facilitate Vienna Ab initio Simulation Package (VASP) calculation inputs
# see https://www.vasp.at/wiki/index.php/The_VASP_Manual

# NOTE that it is assumed here that the KSPACING INCAR tag is used rather than a KPOINTS file
# --> see https://www.vasp.at/wiki/index.php/KSPACING

from ase import Atoms
from pymatgen.io.vasp.inputs import Incar


def set_up_relaxation(
    struc: Atoms,
    output_dir: str = "relaxation",
    user_incar_changes=None,
    functional="r2scan",
) -> None:
    """
    Sets up VASP input files for a ground-state (0K) relaxation calculation from a structure as an ASE Atoms object.

    Args:
        struc (Atoms): structure to be relaxed
        output_dir (str): relative path to output generated files. Defaults to "relaxation".
        user_incar_changes (dict, optional): changes that deviate from default INCAR parameters given in pbe/r2scan.yaml files. Defaults to None.
        functional (str, optional): Exchange-Correlation (XC) functional to be used. Defaults to "r2scan". Options can be found in /pytheos/vasp/inputs_sets/.

    Raises:
        FileExistsError: if output_dir already exists
        ValueError: if an invalid functional is given

    Returns:
        None: New directory is made for a VASP relaxation calculation -> output_dir
    """

    import os
    import yaml
    from pymatgen.core import Structure
    from pymatgen.io.vasp.sets import MPRelaxSet, MPScanRelaxSet
    from pymatgen.io.vasp.inputs import Kpoints

    print(f"----- Setting up VASP relaxation -----")

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
            f"{functional} is not available. Choose either a GGA ('pbe', 'pbesol') or a MetaGGA ('r2scan', 'scan')."
        )

    calc.write_input(f"{output_dir}")

    return None


def set_up_dos(
    source_dir: str = "relaxation",
    output_dir: str = "dos",
    user_incar_changes: dict = None,
) -> None:
    """
    Sets up VASP input files for a DOS calculation from a previous calculation (usually a relaxation calculation).

    The source directory should ideally contain CHGCAR and WAVECAR files.

    Some DOS-specific changes are made to the INCAR by default from the source calculation ->
    - NELECT = NELECT - 1 + 0.999999 (https://www.vasp.at/forum/viewtopic.php?t=17981)
    - EMIN = EFERMI - 8 (more reasonable energy window)
    - EMAX = EFERMI + 6 (more reasonable energy window)
    - ISMEAR = -5 (tetrahedron method with Blöchl corrections)
    - ALGO = Normal (since ISMEAR = -5 can often fail with ALGO = All)
    - NSW = 0 (single ionic step)
    - LCHARG = False (don't need usually)
    - LWAVE = False (don't need usually)
    - NEDOS = 501 (higher resolution DOS sampling)

    Args:
        source_dir (str, optional): Directory to source previous VASP files. Defaults to "relaxation".
        output_dir (str, optional): Directory to output new VASP files for DOS calculation. Defaults to "dos".
        user_incar_changes (dict, optional): Additional changes to INCAR that can supply (ex. {"NCORE": 24}). Defaults to None.

    Returns:
        None: New directory is made for a VASP DOS calculation -> output_dir
    """
    import os
    from pytheos.vasp import load_vasprun
    from pymatgen.io.vasp.outputs import Eigenval
    from pymatgen.io.vasp.inputs import Incar

    print(f"----- Setting up VASP density of states -----")

    # to get back to where this was originally called
    og_path = os.path.abspath(".")

    # load vasprun.xml (also performs a check for sufficient convergence)
    v = load_vasprun(
        path=f"{source_dir}/vasprun.xml",
        parse_dos_eigen=True,
    )

    # make directory for DOS calculation
    os.mkdir(f"{output_dir}")

    # get necessary files from source calculation
    os.system(
        f"cp {source_dir}/INCAR {source_dir}/CONTCAR {source_dir}/POTCAR {source_dir}/CHGCAR {source_dir}/WAVECAR {output_dir}"
    )

    # due to issue with fermi level placement
    e = Eigenval(f"{source_dir}/EIGENVAL")
    num_elec = e.nelect
    new_num_elec = num_elec - 1 + 0.999999

    # for a more reasonable energy window
    efermi = v.efermi
    emin = efermi - 8
    emax = efermi + 6

    os.chdir(f"{output_dir}")
    os.system("mv CONTCAR POSCAR")  # to ensure final structure is used

    # read in previous incar + update with general DOS flags
    incar = Incar.from_file("INCAR")
    incar.update(
        {
            "NELECT": new_num_elec,  # per https://www.vasp.at/forum/viewtopic.php?t=17981
            "EMIN": emin,  # for more reasonable energy window
            "EMAX": emax,  # for more reasonable energy window
            "ISMEAR": -5,  # proper ISMEAR for DOS -> https://www.vasp.at/wiki/index.php/ISMEAR
            "ALGO": "Normal",  # since ALGO = All can often fail with ISMEAR = -5
            "NSW": 0,  # static calculation
            "LCHARG": False,  # usually not needed
            "LWAVE": False,  # usually not needed
            "NEDOS": 501,  # higher resolution DOS
        }
    )

    # only make further changes if user gives them
    if bool(user_incar_changes) == True:
        incar.update(user_incar_changes)

    incar.write_file("INCAR")  # overwrite INCAR with new flags

    # get back to original directory where this was called
    os.chdir(og_path)

    return None


def set_up_bader(
    source_dir: str = "relaxation",
    output_dir: str = "bader",
    user_incar_changes: dict = None,
) -> None:
    """
    Sets up VASP input files for a Bader calculation from a previous calculation (usually a relaxation calculation).
    - details on Bader charge analysis -> https://theory.cm.utexas.edu/henkelman/code/bader/

    The source directory should ideally contain CHGCAR and WAVECAR files.

    Some Bader-specific changes are made to the INCAR by default from the source calculation ->

    - ISMEAR = 0 (Gaussian smearing)
    - NSW = 0 (single ionic step)
    - LAECHG = True (needed for Bader charge analysis - https://www.vasp.at/wiki/index.php/LAECHG)
    - LCHARG = True (needed for Bader charge analysis)

    Args:
        source_dir (str, optional): Directory to source previous VASP files. Defaults to "relaxation".
        output_dir (str, optional): Directory to output new VASP files for DBaderOS calculation. Defaults to "dos".
        user_incar_changes (dict, optional): Additional changes to INCAR that can supply (ex. {"NCORE": 24}). Defaults to None.

    Returns:
        None: New directory is made for a VASP Bader calculation -> output_dir
    """

    import os
    from pytheos.vasp.outputs import load_vasprun
    from pymatgen.io.vasp.inputs import Incar

    print(f"----- Setting up VASP Bader charge -----")

    # to get back to where this was originally called
    og_path = os.path.abspath(".")

    # load vasprun.xml (also performs a check for sufficient convergence)
    v = load_vasprun(
        path=f"{source_dir}/vasprun.xml",
        parse_dos_eigen=False,
    )

    # make directory for bader calculation
    os.mkdir(f"{output_dir}")

    # get necessary files from source calculation
    os.system(
        f"cp {source_dir}/INCAR {source_dir}/CONTCAR {source_dir}/POTCAR {source_dir}/CHGCAR {source_dir}/WAVECAR {output_dir}"
    )

    os.chdir(f"{output_dir}")
    os.system("mv CONTCAR POSCAR")  # to ensure final structure is used

    # read in previous incar + update with general DOS flags
    incar = Incar.from_file("INCAR")
    incar.update(
        {
            "ISMEAR": 0,  # Gaussian smearing
            "NSW": 0,  # static calculation (important for VASP Bader -> https://theory.cm.utexas.edu/henkelman/code/bader/)
            "LAECHG": True,  # needed for Bader charge analysis (https://www.vasp.at/wiki/index.php/LAECHG)
            "LCHARG": True,  # needed for Bader charge analysis
        }
    )

    # only make further changes if user gives them
    if bool(user_incar_changes) == True:
        incar.update(user_incar_changes)

    incar.write_file("INCAR")  # overwrite INCAR with new flags

    # get back to original directory where this was called
    os.chdir(og_path)

    return None


def write_custodian_static(
    output_dir: str,
    vasp_cmd: tuple = ("srun", "vasp_std"),
) -> None:
    """
    Write a generic static Custodian script for calculation workflow and error handling.
    - https://github.com/materialsproject/custodian

    Args:
        output_dir (str): relative path to write submission file.

    Returns:
        None: {output_dir}/cstdn.py script written.
    """

    print(f"----- Writing static Custodian script -----")

    cstdn_script = f"# Custodian static script.\n\nvasp_cmd = {vasp_cmd}\n"

    cstdn_script += """import os
from custodian.custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler
from custodian.vasp.jobs import VaspJob

subset = list(VaspErrorHandler.error_msgs.keys())
handlers = [VaspErrorHandler(errors_subset_to_catch=subset)]

step1 = VaspJob(
    vasp_cmd=vasp_cmd,
    final=True,
)

jobs = [step1]
c = Custodian(handlers, jobs, max_errors=3)
c.run()"""
    with open(f"{output_dir}/cstdn.py", "w+") as f:
        f.writelines(cstdn_script)


def write_custodian_doublerelax(
    output_dir: str,
    vasp_cmd: tuple = ("srun", "vasp_std"),
    kspacing=0.25,
    half_kmesh_first_relax=True,
) -> None:
    """
    Write a generic double relaxation Custodian script for calculation workflow and error handling.
    - https://github.com/materialsproject/custodian

    Restricted to using the "KSPACING" INCAR flag as commonly used for high-throughput calculations
    - KSPACING = 0.25 used as a reasonably safe default -> do you own convergence testing!!

    Using half_kmesh_first_relax = True can be very beneficial to get closer to a structural minima at a drastically reduced cost.
    - see discussion -> https://github.com/hackingmaterials/atomate/issues/19

    Args:
        output_dir (str): relative path to write submission file
        kspacing (float, optional): K-point mesh spacing with VASP KSPACING tag (https://www.vasp.at/wiki/index.php/KSPACING). Defaults to 0.25.
        half_kmesh_first_relax (bool, optional): Use more sparse k-mesh for initial relax. Defaults to True.

    Returns:
        None: {output_dir}/cstdn.py script written.
    """

    print(f"----- Writing double-relax Custodian script -----")

    cstdn_script = f"# Custodian double-relaxation script.\n\nvasp_cmd = {vasp_cmd}\n"

    if half_kmesh_first_relax == True:
        cstdn_script += (
            f"""kspacing_initial = {kspacing*2}\nkspacing = {kspacing}\n\n"""
        )
    else:
        cstdn_script += f"""kspacing_initial = {kspacing}\nkspacing = {kspacing}\n\n"""

    cstdn_script += """import os
from custodian.custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler
from custodian.vasp.jobs import VaspJob

subset = list(VaspErrorHandler.error_msgs.keys())
handlers = [VaspErrorHandler(errors_subset_to_catch=subset)]

step1 = VaspJob(
    vasp_cmd=vasp_cmd,
    final=False,
    suffix=".1",
    settings_override=[
        {
            "dict": "INCAR",
            "action": {
                "_set": {
                    "KSPACING": kspacing_initial
                }
            },
        },
    ],
)

step2 = VaspJob(
    vasp_cmd=vasp_cmd,
    final=True,
    settings_override=[
        {
            "dict": "INCAR",
            "action": {
                "_set": {
                    "KSPACING": kspacing
                }
            },
        },
        {"file": "CONTCAR", "action": {"_file_copy": {"dest": "POSCAR"}}},
    ],
)


jobs = [step1, step2]
c = Custodian(handlers, jobs, max_errors=3)
c.run()"""
    with open(f"{output_dir}/cstdn.py", "w+") as f:
        f.writelines(cstdn_script)


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


# TODO just copied in for now from the perovskite HEO project - still needs to be fixed
def set_up_optical_calc(sigma=0.1) -> None:
    """
    For making a new directory for optical calculations (should be called in directory where "02_static" exists).
    - NBANDS is doubled from the default value of 02_static calculation

    Args:
        sigma (float, optional): Smearing parameter than is usually increased for optical calculations compared to relaxation. Defaults to 0.2.
    """
    import os
    from pymatgen.io.vasp.outputs import Eigenval, Vasprun

    print(f"-> Making optical directory...")

    os.mkdir("04_optical")
    os.system(
        f"cp 02_static/INCAR 02_static/POSCAR 02_static/KPOINTS 02_static/POTCAR 02_static/CHGCAR 02_static/WAVECAR 02_static/submitvasp ./04_optical"
    )
    os.chdir("04_optical")

    e = Eigenval("../02_static/EIGENVAL")

    nbands = e.nbands
    new_nbands = nbands * 2  # doubling NBANDS
    os.system(f"echo NBANDS = {new_nbands} >> INCAR")
    print(f"number of bands: {nbands} -> {new_nbands}")

    os.system("perl -pi -e 's/LWAVE = True/LWAVE = False/g' INCAR")
    os.system("perl -pi -e 's/LCHARG = True/LCHARG = False/g' INCAR")
    os.system(f"echo NEDOS = 20000 >> INCAR")  # for adequate sampling
    os.system(f"perl -pi -e 's/SIGMA = 0.05/SIGMA = {sigma}/g' INCAR")
    os.system("perl -pi -e 's/ALGO = All/ALGO = Normal/g' INCAR")
    os.system("perl -pi -e 's/EDIFF = 1e-06/EDIFF = 1e-08/g' INCAR")
    os.system("perl -pi -e 's/LREAL = Auto/LREAL = False/g' INCAR")
    os.system(f"echo LOPTICS = True >> INCAR")
    os.system(f"echo CSHIFT = 0.01 >> INCAR")


####################################


# TODO just copied in for now from the perovskite HEO project - still needs to be fixed
def set_up_bandstructure_calc():
    """For making a new directory for density of states calculations (should be called in directory where "02_static" exists)
    - NBANDS is 1.5x from the default value of 02_static calculation

    Additional, manual steps are required for this function still due to band structure process:
        NOTE that this is for the 'regular' band structure calculations (i.e., not unfolding)
        1. run 01_pbe calculation
        2. get high-symmetry k-path with sumo (> sumo-kgen --hybrid --symprec 0.1 --pymatgen)
        3. copy KPOINTS_band to KPOINTS
        4. copy WAVECAR from finished 01_pbe calculation to the 'bandstructure' directory
        5. now can run actual metaGGA band structure calculation
    """

    import os
    from pymatgen.io.vasp.outputs import Eigenval

    print(f"-> Making bandstructure directory...")

    os.mkdir("05_bandstructure")
    os.system(
        f"cp 02_static/INCAR 02_static/IBZKPT 02_static/POSCAR 02_static/KPOINTS 02_static/POTCAR 02_static/CHGCAR 02_static/WAVECAR 02_static/submitvasp ./05_bandstructure"
    )
    os.chdir("05_bandstructure")

    e = Eigenval("../02_static/EIGENVAL")

    nbands = e.nbands
    new_nbands = int(nbands * 1.5)  # increasing NBANDS by 1.5x
    os.system(f"echo NBANDS = {new_nbands} >> INCAR")
    print(f"number of bands: {nbands} -> {new_nbands}")

    os.system("perl -pi -e 's/LCHARG = True/LCHARG = False/g' INCAR")
    os.system(f"echo NEDOS = 5001 >> INCAR")  # for adequate sampling
    os.system("perl -pi -e 's/LREAL = Auto/LREAL = False/g' INCAR")

    os.mkdir("01_pbe")
    os.system("cp * 01_pbe")
    os.system("rm CHGCAR WAVECAR")
    os.chdir("01_pbe")

    os.system("perl -pi -e 's/METAGGA = R2scan/GGA = PE/g' INCAR")
    os.system(f"echo ICHARG = 11 >> INCAR")
