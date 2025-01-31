# for VASP calculation outputs

from pymatgen.io.vasp import Vasprun
from ase import Atoms


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


def calc_mp2020compat_energy(run: Vasprun):
    """_summary_

    Args:
        run (Vasprun): _description_

    Returns:
        _type_: _description_
    """
    from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
    import numpy as np

    v = run.get_computed_entry()

    # get original energy in eV/atom
    energy_og = v.energy / len(v.structure)
    print(f"\noriginal energy = {np.round(energy_og, 4)}/atom")

    # calculate corrected energy with MP2020Compatibility corrections
    v.energy_adjustments = MaterialsProject2020Compatibility().get_adjustments(v)
    energy_mp2020 = v.energy / len(v.structure)
    print(f"corrected energy = {np.round(energy_mp2020, 4)}/atom\n")

    return energy_mp2020


def calc_form_decomp_energy(
    struc: Atoms,
    energy_per_atom: float,
    MPApiKey: str,
    xc: str = "GGA_GGA+U",
):
    """
    Calculates the formation and decomposition energies using the Materials Project database.
    - see C.J. Bartel et al. Npj Comput Mater 5 (2019) 4. https://doi.org/10.1038/s41524-018-0143-2 for more details on computation

    Calculation parameters/potcars should be consistent with MPRelaxSet/MPScanRelaxSet for valid computations.
    - https://github.com/materialsproject/pymatgen/tree/master/src/pymatgen/io/vasp
    - Any energy correction/mixing schemes (such as MP2020Compatibility should already have been applied to the total energy)
        - https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/thermodynamic-stability

    Args:
        struc (Atoms): ASE Atoms object for structure.
        energy_per_atom (float): Total energy in eV/atom from VASP calculation.
        MPApiKey (str): Materials Project API Key (https://next-gen.materialsproject.org/api) - user specific.
        xc (str, optional): Exchange-Correlation functional used in calculation. Options -> ["GGA_GGA+U", "R2SCAN", "GGA_GGA+U_R2SCAN"] per https://github.com/materialsproject/emmet/blob/main/emmet-core/emmet/core/thermo.py. Defaults to "GGA_GGA+U".

    Returns:
        dict: A dictionary with the following keys and values:
        - "form_energy": formation energy in eV/atom
        - "decomp_energy": decomposition energy in eV/atom
        - "decomp_rxn": decomposition reaction
    """

    from pymatgen.core import Structure
    from mp_api.client import MPRester
    from pymatgen.analysis import phase_diagram
    from emmet.core.thermo import ThermoType

    struc = Structure.from_ase_atoms(struc)
    chemical_formula = (
        struc.composition.formula
    )  # NOTE might need to correct this, TODO

    # for our system of interest
    target_entry = phase_diagram.PDEntry(
        composition=chemical_formula,
        energy=energy_per_atom,
        name="target",  # need some unique name for identification purposes
    )

    # get all elements from the entry of interest and places them into a list with the required format
    elements = target_entry.composition.elements
    elements_newlist = []
    for element in range(len(elements)):
        elements_newlist.append(elements[element].symbol)

    # get all entries from Materials Project that contain the elements of interest
    with MPRester(MPApiKey) as mpr:
        entries = mpr.get_entries_in_chemsys(
            elements_newlist,
            additional_criteria={"thermo_types": [ThermoType[xc]]},
        )

    # add in our target entry to the entries
    entries.append(target_entry)

    # compute phase diagram
    phasediagram = phase_diagram.PhaseDiagram(entries)
    all_entries = phasediagram.all_entries

    for e in all_entries:
        if e.name == "target":  # only for our entry of interest
            e_form = phasediagram.get_form_energy_per_atom(e)
            e_decomp = phasediagram.get_decomp_and_phase_separation_energy(e)[1]
            decomp = phasediagram.get_decomp_and_phase_separation_energy(e)[0]
            decomp_formula = ""

            # for getting the decomposition reaction in a nicer format
            for decomp_entry in range(len(list(decomp))):
                decomp_formula += "{:.2f}".format(list(decomp.values())[decomp_entry])
                decomp_formula += "({})".format(
                    list(decomp.keys())[decomp_entry].composition.reduced_formula
                )
                if decomp_entry != range(len(list(decomp.keys())))[-1]:
                    decomp_formula += " + "

    return {
        "form_energy": e_form,
        "decomp_energy": e_decomp,
        "decomp_rxn": decomp_formula,
    }


def extract_optical_data(run="vasprun.xml", anisotropic=False):
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
