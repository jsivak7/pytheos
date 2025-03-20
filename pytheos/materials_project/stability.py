# evaluate stability using MP database
# see https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability

from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram


def make_phasediagram_entry(
    composition: str,
    energy: float,
    name="my_PDEntry",
) -> PDEntry:
    """
    Makes new phase diagram entry (PDEntry) object..

    Usually this is added to queried entries and used to generate PhaseDiagram and determine
    relative stability for your PDEntry. PDEntries contain all of the data needed to construct
    a PhaseDiagram

    NOTE that you have to be very careful that calculation parameters match those you will be
    comparing this entry to!! Remember that MP applies different mixing schemes across their
    database -> https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/thermodynamic-stability.

    More info on PDEntries can be found at -> https://github.com/materialsproject/pymatgen/blob/master/src/pymatgen/analysis/phase_diagram.py.

    Args:
        composition (str): Material composition - include all atoms. Example -> "MgCoNiO3".
            Preferentially provide Pymatgen Composition to minimize possible errors.
        energy (float): Energy (in eV) corresponding to composition.
        name (str, optional): Name for entry - need some unique name for identification purposes. Defaults to "my_PDEntry".

    Returns:
        PDEntry: Pymatgen PDEntry.
    """

    my_PDEntry = PDEntry(
        composition=composition,
        energy=energy,
        name=name,
    )
    return my_PDEntry


def generate_phase_diagram(
    entries: list,
) -> PhaseDiagram:
    """
    Generates a phase diagram for a given set of entries using Pymatgen.

    One can find the algorithm used in the implemented PhaseDiagram class.
    - https://github.com/materialsproject/pymatgen/blob/master/src/pymatgen/analysis/phase_diagram.py

    Args:
        entries (list): list of ComputedEntries/ComputedStructureEntries/PDEntries.

    Returns:
        PhaseDiagram: Pymatgen PhaseDiagram class
    """
    phase_diagram = PhaseDiagram(entries=entries)

    return phase_diagram


def calc_form_energy(
    phasediagram: PhaseDiagram,
    name_target: str = "my_PDEntry",
) -> float:
    """
    Calculates the formation energy of targetted entry using supplied PhaseDiagram.

    The formation energy is the energy relative to the elemental references.

    Any energy correction schemes should already have been applied to the total energy
    - see /pytheos/materials_project/corrections for more detail

    Args:
        phasediagram (PhaseDiagram): Pymatgen PhaseDiagram with target entry added.
        name_target (str, optional): PDEntry name for target. Defaults to "my_PDEntry".

    Returns:
        float: Calculated formation energy in eV/atom
    """

    for entry in phasediagram.all_entries:
        if entry.name == name_target:
            form_energy = phasediagram.get_form_energy_per_atom(entry=entry)

    return form_energy


def calc_decomp_energy(
    phasediagram: PhaseDiagram,
    name_target: str = "my_PDEntry",
) -> float:
    """
    Calculates the decomposition energy of targetted entry using supplied PhaseDiagram.

    The decomposition energy is a measure of instability with respect to phase separation.
    - see C.J. Bartel et al. Npj Comput Mater 5 (2019) 4. https://doi.org/10.1038/s41524-018-0143-2
    - appears to be very useful for determining enthalpy cost of materials, especially for HEOs

    Any energy correction schemes should already have been applied to the total energy
    - see /pytheos/materials_project/corrections for more detail

    Args:
        phasediagram (PhaseDiagram): Pymatgen PhaseDiagram with target entry added.
        name_target (str, optional): PDEntry name for target. Defaults to "my_PDEntry".

    Returns:
        float: Calculated decomposition energy in eV/atom
    """

    for entry in phasediagram.all_entries:
        if entry.name == name_target:
            decomp_energy = phasediagram.get_decomp_and_phase_separation_energy(
                entry=entry
            )[1]

    return decomp_energy


def calc_decomp_rxn(
    phasediagram: PhaseDiagram,
    name_target: str = "my_PDEntry",
) -> str:
    """
    Calculates the decomposition reaction corresponding to the decomposition energy for
    the targetted entry using supplied PhaseDiagram.

    These reactions are quite useful in comparing to experimental synthesis attempts and to
    determine what phases might inhibit single-phase formation of a material.

    Any energy correction schemes should already have been applied to the total energy
    - see /pytheos/materials_project/corrections for more detail

    Args:
        phasediagram (PhaseDiagram): Pymatgen PhaseDiagram with target entry added.
        name_target (str, optional): PDEntry name for target. Defaults to "my_PDEntry".

    Returns:
        float: Calculated decomposition energy in eV/atom
    """
    for entry in phasediagram.all_entries:
        if entry.name == name_target:
            decomp_rxn = phasediagram.get_decomp_and_phase_separation_energy(
                entry=entry
            )[0]

    cleaned_decomp_rxn = clean_up_decomp_rxn(decomp_entries=decomp_rxn)

    return cleaned_decomp_rxn


def clean_up_decomp_rxn(
    decomp_entries: dict,
) -> str:
    """
    Cleans up decomposition reaction from dictionary of ComputedStructureEntriesto a simpler,
    more human-readable string with relative entry amounts.

    Args:
        decomp_entries (dict): Dictionaries of ComputedStructureEntries outputted.

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
