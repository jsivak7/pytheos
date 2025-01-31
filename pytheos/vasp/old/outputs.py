"""
For extracting data from VASP calculation outputs
"""


def get_formula(run="vasprun.xml"):
    """
    Extracts the chemical formula for a VASP calculation

    Args:
        run (str, optional): relative location for vasprun.xml file. Defaults to "vasprun.xml".

    Returns:
        tuple: (formula, reduced_formula, number_atoms)
    """
    from pymatgen.io.vasp import Vasprun

    v = Vasprun(run)
    formula = v.final_structure.composition.formula
    reduced_formula = v.final_structure.composition.reduced_formula
    num_atoms = v.final_structure.composition.num_atoms

    return (formula, reduced_formula, num_atoms)


def get_total_energy(run="vasprun.xml"):
    """
    Extracts total energy (in eV) for a VASP calculation

    Args:
        run (str, optional): relative location for vasprun.xml file. Defaults to "vasprun.xml".

    Returns:
        float: total energy in eV
    """
    from pymatgen.io.vasp import Vasprun

    v = Vasprun(run)
    en_tot = v.final_energy
    return en_tot


# TODO
def get_mp2020compat_energy(run="vasprun.xml"):
    """
    Returns the energy per atom for a given PBE calculation using the MP2020Compatibility correction scheme
    """
    from pymatgen.io.vasp import Vasprun
    from pymatgen.entries.compatibility import MaterialsProject2020Compatibility

    v = Vasprun(run).get_computed_entry()
    en_og = v.energy / len(v.structure)
    print("\noriginal energy = {:.6f}/atom".format(en_og))
    v.energy_adjustments = MaterialsProject2020Compatibility().get_adjustments(v)
    en_corr = v.energy / len(v.structure)
    print("corrected energy = {:.6f}/atom\n".format(en_corr))
    return en_corr


def get_lattice(run="vasprun.xml"):
    """
    Extract volume and lattice parameters from a VASP calculation.

    Args:
        run (str, optional): relative location for vasprun.xml file. Defaults to "vasprun.xml".

    Returns:
        tuple: (volume, a, b, c, alpha, beta, gamma)
    """
    from pymatgen.io.vasp import Vasprun
    import numpy as np

    v = Vasprun(run)
    volume = v.final_structure.volume
    a = v.final_structure.lattice.a
    b = v.final_structure.lattice.b
    c = v.final_structure.lattice.c
    alpha = v.final_structure.lattice.alpha
    beta = v.final_structure.lattice.beta
    gamma = v.final_structure.lattice.gamma

    return (
        np.round(volume, 2),
        np.round(a, 4),
        np.round(b, 4),
        np.round(c, 4),
        np.round(alpha, 2),
        np.round(beta, 2),
        np.round(gamma, 2),
    )


def get_bandgap(run="vasprun.xml"):
    """
    Extracts electronic band gap from a VASP calculation.
    - NOTE to be careful given the recently issues with the ISMEAR = -5 issues with Fermi energy placement
        (https://www.vasp.at/forum/viewtopic.php?t=18407)
        - recommendation from JTSivak for this issue -> (https://www.vasp.at/forum/viewtopic.php?t=17981)
    - one should also be careful to have an adequately sampled DOS (NEDOS = large, odd value)

    Args:
        run (str, optional): relative location for vasprun.xml file. Defaults to "vasprun.xml".

    Returns:
        float: electronic band gap in eV
    """
    from pymatgen.io.vasp import Vasprun
    import numpy as np

    v = Vasprun(run)
    bandgap = v.eigenvalue_band_properties[0]

    return np.round(bandgap, 2)


def get_optical(run="vasprun.xml", anisotropic=False):
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
        os.system(
            "perl -pi -e 's/# //g' {}.dat".format(property)
        )  # cleaning up .dat file to .csv
        os.system(
            "perl -pi -e 's/ /,/g' {}.dat".format(property)
        )  # cleaning up .dat file to .csv
        os.system(
            "perl -pi -e 's/alpha/{}/g' {}.dat".format(property, property)
        )  # change 'alpha' to actual property (unsure of why always defaults to alpha)
        os.system("mv {}.dat {}.csv".format(property, property))
        os.system(
            "rm absorption.pdf"
        )  # dont want all of these since plotting myself (for some reason always defaults to absorption.pdf in sumo)

        # getting energy in wavelength (nm) and making new column in .csv file
        data = pd.read_csv("./{}.csv".format(property))
        energy_joules = data["energy(eV)"] * 1.602176634e-19  # convert energy to joules
        wavelength = (
            6.626e-34 * 2.998e8
        ) / energy_joules  # calculate wavelength from energy
        wavelength = wavelength * (1e9)  # get in nm
        data["wavelength(nm)"] = wavelength
        # data = data[[col for col in data.columns if col != property] + [property]] # move property column to the end of the dataframe
        data.to_csv("{}.csv".format(property), index=False)

        # collect
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

        os.system(
            "rm {}.csv".format(property)
        )  # get rid of all of the individual data files

    df_final = pd.DataFrame().from_dict(data_all)

    return df_final
