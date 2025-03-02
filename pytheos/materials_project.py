# for extracting data and performing computations using the Materials Project database
# https://next-gen.materialsproject.org

# TODO add in chemical potential overlap extraction and computation
# TODO add in general database extraction for materials disvoery exploration

from pymatgen.io.vasp.outputs import Vasprun
from ase import Atoms


def calc_mp2020compat_energy(run: Vasprun) -> float:
    """
    Calculates the corrected energy using the MP2020Compatbility correction scheme for GGA/GGA+U and anion mixing calculations.
    - https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/thermodynamic-stability/anion-and-gga-gga+u-mixing

    Calculation parameters/potcars should be consistent with MPRelaxSet for valid computations.
    - https://github.com/materialsproject/pymatgen/blob/master/src/pymatgen/io/vasp/MPRelaxSet.yaml

    Args:
        run (Vasprun): Pymatgen vasprun object. Used preferentially over raw energies to ensure scheme is correctly implemented.

    Returns:
        float: Corrected energy in eV/atom
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
    energy: float,
    MPApiKey: str,
    xc: str = "R2SCAN",
) -> tuple:
    """
    Calculates the formation and decomposition energies using the Materials Project database.
    - see C.J. Bartel et al. Npj Comput Mater 5 (2019) 4. https://doi.org/10.1038/s41524-018-0143-2 for more details on computation

    Calculation parameters/potcars should be consistent with MPRelaxSet/MPScanRelaxSet for valid computations.
    - https://github.com/materialsproject/pymatgen/tree/master/src/pymatgen/io/vasp
    - Any energy correction/mixing schemes (such as MP2020Compatibility should already have been applied to the total energy)
        - https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/thermodynamic-stability

    Args:
        struc (Atoms): ASE Atoms object for structure.
        energy (float): Total energy in eV from VASP calculation.
        MPApiKey (str): Materials Project API Key (https://next-gen.materialsproject.org/api) - user specific.
        xc (str, optional): Exchange-Correlation functional used in calculation. Options -> ["GGA_GGA_U", "R2SCAN", "GGA_GGA_U_R2SCAN"] per https://github.com/materialsproject/emmet/blob/main/emmet-core/emmet/core/thermo.py. Defaults to "R2SCAN".

    Returns:
        tuple: (formation energy in eV/atom, decomposition energy in eV/atom, decomposition reaction entries)
    """

    from pymatgen.core import Structure
    from mp_api.client import MPRester
    from pymatgen.analysis import phase_diagram
    from emmet.core.thermo import ThermoType

    struc = Structure.from_ase_atoms(struc)
    chemical_formula = struc.composition.formula

    # for our system of interest
    target_entry = phase_diagram.PDEntry(
        composition=chemical_formula,
        energy=energy,
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
            decomp_entries = phasediagram.get_decomp_and_phase_separation_energy(e)[0]

    return (e_form, e_decomp, decomp_entries)


def cleanup_decomp_rxn(decomp_entries: dict) -> str:
    """
    Cleans up decomposition reaction from Pymatgen ComputedStructureEntries to a simpler, human-readable string.

    Args:
        decomp_entries (dict): Dictionaries of ComputedStructureEntries outputted the calc_form_decom_energy function "decomp_entries"

    Returns:
        string: Decomposition reaction - ex. 0.25(MgO) + 0.75(NiO)
    """
    from pymatgen.analysis import phase_diagram

    decomp_formula = ""

    for decomp_entry in range(len(list(decomp_entries))):
        decomp_formula += "{:.2f}".format(list(decomp_entries.values())[decomp_entry])
        decomp_formula += "({})".format(
            list(decomp_entries.keys())[decomp_entry].composition.reduced_formula
        )
        if decomp_entry != range(len(list(decomp_entries.keys())))[-1]:
            decomp_formula += " + "

    return decomp_formula
