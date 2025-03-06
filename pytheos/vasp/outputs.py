# module to facilitate Vienna Ab initio Simulation Package (VASP) calculation outputs
# see https://www.vasp.at/wiki/index.php/The_VASP_Manual

from pymatgen.io.vasp import Vasprun
from pandas import DataFrame


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


def convert_baderACF_from_dat_to_csv() -> None:
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
