# For creating VASP input files and directories
# Assumes KSPACING INCAR tag is used rather than a KPOINTS file (https://www.vasp.at/wiki/index.php/KSPACING)

from ase import Atoms


def set_up_relaxation(
    struc: Atoms,
    output_dir: str = "relaxation",
    user_incar_changes=None,
    functional="r2scan",
) -> None:
    """
    Sets up VASP input files for a relaxation calculation from a structure as an ASE Atoms object.

    Args:
        struc (Atoms): structure to be relaxed
        output_dir (str): relative path to output generated files. Defaults to "relaxation".
        user_incar_changes (dict, optional): changes that deviate from default INCAR parameters given in pbe/r2scan.yaml files. Defaults to None.
        functional (str, optional): Exchange-Correlation (XC) functional to be used. Defaults to "r2scan".

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

    print(f"\nSetting up VASP relaxation inputs...")

    os.mkdir(output_dir)

    # read in structure file as ASE Atoms object
    s = Structure.from_ase_atoms(struc)

    print(f"-> structure = {s.composition.formula}")
    print(f"-> output directory = {output_dir}")
    print(f"-> user incar changes = {user_incar_changes}\n")

    # get specific functional (xc.yaml) file for incar settings
    module_dir = os.path.dirname(__file__)  # get path to module
    with open(f"{module_dir}/{functional}.yaml", "r") as f:
        incar_settings = yaml.load(f, Loader=yaml.SafeLoader)

    # only make further changes if user gives them
    if bool(user_incar_changes) == True:
        incar_settings.update(user_incar_changes)

    if functional == "pbe":
        calc = MPRelaxSet(
            structure=s,
            user_incar_settings=incar_settings,
            user_potcar_functional="PBE",
            sort_structure=True,
        )
    elif functional == "r2scan":
        calc = MPScanRelaxSet(
            structure=s,
            user_incar_settings=incar_settings,
            user_potcar_functional="PBE_54",
            sort_structure=True,
        )
    else:
        raise ValueError(
            f"{functional} is not available. Choose either 'pbe' or 'r2scan'."
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
    from pytheos.vasp.outputs import load_vasprun
    from pymatgen.io.vasp.outputs import Eigenval
    from pymatgen.io.vasp.inputs import Incar

    print(f"\nSetting up VASP density of states inputs...")
    print(f"-> source directory = {source_dir}")
    print(f"-> output directory = {output_dir}")
    print(f"-> user incar changes = {user_incar_changes}\n")

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

    print(f"\nSetting up VASP Bader charge inputs...")
    print(f"-> source directory = {source_dir}")
    print(f"-> output directory = {output_dir}")
    print(f"-> user incar changes = {user_incar_changes}\n")

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


###################### TODO #######################


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


# TODO need to ensure this is properly implemented
def get_magorder(calc_path: str, magorder_name: str, coordinate_rounding=0):

    from ase.io import read
    import numpy as np
    import random
    import os
    import yaml

    s = read(f"{calc_path}/POSCAR")

    coords = s.get_positions()
    coords = np.around(
        coords, coordinate_rounding
    )  # optional to change coordinate rounding of inputted structure

    elems = s.get_chemical_symbols()

    module_dir = os.path.dirname(__file__)

    with open(
        f"{module_dir}/incar_defaults/r2scan.yaml", "r"
    ) as f:  # need to make this more general...
        incar_settings = yaml.load(f, Loader=yaml.SafeLoader)
    magmoms = incar_settings["MAGMOM"]

    with open(f"{module_dir}/magnetic_orders/{magorder_name}.yaml", "r") as f:
        spin_ordering = yaml.load(f, Loader=yaml.SafeLoader)

    spin_up = np.around(np.array(spin_ordering["spin_up"]), 1)
    spin_down = np.around(np.array(spin_ordering["spin_down"]), 1)

    up_index = []
    down_index = []

    for i in range(len(coords)):
        for j in range(
            len(spin_up)
        ):  # NOTE that this needs to be an equal amount of spin up and down with the current implementation
            if (
                np.array_equal(spin_up[j], coords[i]) == True
            ):  # this tests for if the atomic coordinates match what the user defines for a spin-up atom
                up_index.append(
                    i
                )  # if it matches, then append the atom number to the up_index list
            elif (
                np.array_equal(spin_down[j], coords[i]) == True
            ):  # this tests for if the atomic coordinates match what the user defines for a spin-up atom
                down_index.append(
                    i
                )  # if it matches, then append the atom number to the down_index list

    magmom = []  # empty list that will be in the same order as poscar file
    # loop for creating a list that will correspond to the final magmom tag that will be added to incar file
    for i in range(
        len(coords)
    ):  # loops through all of the atomic numbers within the inputted structure
        if i in up_index:
            magmom.append(
                +np.round(
                    random.uniform(magmoms[elems[i]] - 0.5, magmoms[elems[i]] + 0.5), 2
                )
            )
        elif i in down_index:
            magmom.append(
                -np.round(
                    random.uniform(magmoms[elems[i]] - 0.5, magmoms[elems[i]] + 0.5), 2
                )
            )
        elif elems[i] == "O":
            magmom.append(magmoms["O"])
        elif elems[i] == "Sr":
            magmom.append(magmoms["Sr"])
        elif elems[i] == "La":
            magmom.append(magmoms["La"])
        else:
            print(f"WARNING!!!! An atom was not identified. Atom #{i}")

    print(f"\nnumber spin up atoms:\t{len(up_index)}")
    print(f"number spin down atoms:\t{len(down_index)}")

    print(f"magnetic ordering:\t{magmom}")
    if len(up_index) != len(down_index):
        print("WARNING!!!! The up index does not equal the down index.")

    # write the magmom list (the actual magnetic ordering) to a temporary file "tmp_ORDER.txt"
    with open("tmp_ORDER.txt", "w") as output:
        output.write(str(magmom))

    # the following are so that the magnetic ordering is added to the MAGMOM tag in the proper format
    os.system("perl -pi -e 's/\\[//g' tmp_ORDER.txt")  # get rid of bracket
    os.system("perl -pi -e 's/\\]//g' tmp_ORDER.txt")  # get rid of bracket
    os.system("perl -pi -e 's/,//g' tmp_ORDER.txt")  # get rid of commas
    os.system("echo 'MAGMOM = CHANGE' >> MAGMOM.txt")
    os.system(f"perl -pi -e 's/MAGMOM.*/MAGMOM = CHANGE/g' {calc_path}/INCAR")
    os.system(
        f"perl -pe 's/CHANGE/`cat tmp_ORDER.txt`/e' -i {calc_path}/INCAR"
    )  # this replaces the "CHANGE" that we just added with the magnetic ordering for the corresponding POSCAR file
    os.system("rm MAGMOM.txt ")
    os.system(f"mv tmp_ORDER.txt {calc_path}/magorder")
