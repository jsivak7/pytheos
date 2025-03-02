# module to facilitate Vienna Ab initio Simulation Package (VASP) calculations
# see https://www.vasp.at/wiki/index.php/The_VASP_Manual
# NOTE that it is assumed that the KSPACING incar tag is used rather than a KPOINTS file
# --> see https://www.vasp.at/wiki/index.php/KSPACING

from ase import Atoms
from pymatgen.io.vasp import Vasprun
from pandas import DataFrame


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


def load_vasprun(
    path: str = "vasprun.xml",
    parse_dos_eigen: bool = False,
) -> Vasprun:
    """
    Loads vasprun.xml file for further analyses. This is a very powerful class that can be used to get most common data from VASP calculations.
    - https://github.com/materialsproject/pymatgen/blob/master/src/pymatgen/io/vasp/outputs.py

    Args:
        path (str, optional): Relative path for *.xml file. Defaults to "vasprun.xml".
        parse_dos_eigen (bool, optional): Option to parse DOS and Eigenvalues. Can save significant time if not interested in those data.. Defaults to False.

    Raises:
        ValueError: If sufficient convergence has not been achieved.

    Returns:
        Vasprun: Pymatgen.io.vasp.Vasprun object
    """

    from pymatgen.io.vasp import Vasprun

    v = Vasprun(filename=path)

    # to ensure that calculation has reached convergence
    if v.converged == False:

        # to give the user some more specifics
        if v.converged_ionic == False:
            print("Ionic convergence has not been achieved!")
        if v.converged_electronic == False:
            print("Electronic convergence has not been achieved!")

        # don't allow object to be returned since not sufficiently converged.
        raise ValueError(
            f"VASP calculation has not converged ({path}). See specifics above."
        )

    return v


def bader_convert_ACF_dat_to_csv() -> None:
    """
    Converts the ACF.dat output file from VASP Bader charge analysis to a .csv file with the name ACF.csv
    - https://theory.cm.utexas.edu/henkelman/code/bader/
    - NOTE: call in the same location as VASP output files.
    """

    import os

    os.system("cp ACF.dat ACF.csv")

    os.system("perl -pi -e 's/    /,/g' ACF.csv")  # four spaces
    os.system("perl -pi -e 's/   /,/g' ACF.csv")  # three spaces
    os.system("perl -pi -e 's/  /,/g' ACF.csv")  # two spaces
    os.system("perl -pi -e 's/ /,/g' ACF.csv")  # one space

    os.system("perl -pi -e 's/-//g' ACF.csv")  # dashes removed

    os.system("perl -pi -e 's/,,,/,/g' ACF.csv")
    os.system("perl -pi -e 's/,,/,/g' ACF.csv")

    os.system("sed -i 's/^.\{1\}//g' ACF.csv")

    os.system("perl -pi -e 's/#/ATOM_NUM/g' ACF.csv")

    os.system("perl -pi -e 's/MIN,DIST/MIN_DIST/g' ACF.csv")
    os.system("perl -pi -e 's/ATOMIC,VOL/ATOMIC_VOL/g' ACF.csv")

    os.system("perl -pi -e 's/VACUUM.*//g' ACF.csv")  # ending removed
    os.system("perl -pi -e 's/NUMBER.*//g' ACF.csv")  # ending removed

    os.system("sed -i '/^$/d' ACF.csv")  # remove all blank lines


def extract_optical_data(run="vasprun.xml", anisotropic=False) -> DataFrame:
    """
    Uses Sumo (https://github.com/SMTG-Bham/sumo) to extract dielectric function from VASP calculation
    and then converts it to the real part of the refractive index.

    This is currently for calculations that have been run with LOPTICS = True
    - default sumo output is a .dat file, but this function converts it to a .csv file
    - also gets energy to wavelength conversion
    - NOTE needs to be run in directory where VASP loptics calculation was run

    Extracts the following properties:
        - absorption
        - loss
        - eps_real
        - eps_imag
        - n_real
        - n_imag

    Args:
        - run (str, optional): relative location for vasprun.xml file. Defaults to "vasprun.xml".
            - default = './vasprun.xml'
        - anisotropic (bool)

    Returns:
        Pandas DataFrame with optical data
    """

    import os
    import pandas as pd

    properties = [
        "absorption",  # optical absorption
        "loss",  # energy-loss function -Im(1/eps)
        "eps_real",  # real part of dielectric function
        "eps_imag",  # imaginary part of dielectric function
        "n_real",  # real part of refractive index
        "n_imag",  # imaginary part of refractive index
    ]

    data_all = {}  # will be populated

    print("\nStarting Optical Property Extraction...")

    if anisotropic == True:
        print("**anisotropic**")

    for property in properties:
        print("--> {}".format(property))

        if anisotropic == True:
            os.system(
                "sumo-optplot {} --anisotropic --filenames {}".format(property, run)
            )
        else:
            os.system("sumo-optplot {} --filenames {}".format(property, run))

        # cleaning up .dat file to .csv
        os.system("perl -pi -e 's/# //g' {}.dat".format(property))
        os.system("perl -pi -e 's/ /,/g' {}.dat".format(property))

        # change 'alpha' to actual property (unsure of why always defaults to alpha in sumo)
        os.system("perl -pi -e 's/alpha/{}/g' {}.dat".format(property, property))
        os.system("mv {}.dat {}.csv".format(property, property))

        # dont want all of these since plotting myself (for some reason always defaults to absorption.pdf in sumo)
        os.system("rm absorption.pdf")

        # getting energy in wavelength (nm) and making new column in .csv file
        data = pd.read_csv("./{}.csv".format(property))

        # convert energy to joules
        energy_joules = data["energy(eV)"] * 1.602176634e-19

        # calculate wavelength from energy
        wavelength = (6.626e-34 * 2.998e8) / energy_joules
        wavelength = wavelength * (1e9)  # convert to nm

        data["wavelength(nm)"] = wavelength
        # data = data[[col for col in data.columns if col != property] + [property]] # move property column to the end of the dataframe
        data.to_csv("{}.csv".format(property), index=False)

        # collect all data
        df = pd.read_csv("{}.csv".format(property))

        if property == "absorption":  # so that we only get energy and wavelength once
            data_all["energy_eV"] = df["energy(eV)"]
            data_all["wavelength_nm"] = df["wavelength(nm)"]

        if anisotropic == True:
            data_all["{}_xx".format(property)] = df["{}_xx".format(property)]
            data_all["{}_yy".format(property)] = df["{}_yy".format(property)]
            data_all["{}_zz".format(property)] = df["{}_zz".format(property)]
        else:
            data_all[property] = df[property]

        # get rid of all of the individual data files
        os.system("rm {}.csv".format(property))

    df_final = pd.DataFrame().from_dict(data_all)

    return df_final
