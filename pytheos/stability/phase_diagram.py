# for constructing phase diagrams

from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram


def make_new_pd_entry(composition: str, energy: float, name="my_PDEntry") -> PDEntry:
    """
    Makes new phase diagram entry (PDEntry) object - useful for adding your own calculation
    to a queried entry list before generating a phase diagram.

    NOTE that you have to be very careful that calculation parameters match those you will be
    comparing this entry to!! Remember that MP applies different mixing schemes across their database
    - https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/thermodynamic-stability.

    Args:
        composition (str): Material composition - include all atoms. e.g. "MgCoNiO3".
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


def generate_phase_diagram(entries: list) -> PhaseDiagram:
    """
    Generates a phase diagram for a given set of entries using Pymatgen.

    https://github.com/materialsproject/pymatgen/blob/master/src/pymatgen/analysis/phase_diagram.py

    Args:
        entries (list): list of ComputedEntries/ComputedStructureEntries/PDEntries.

    Returns:
        PhaseDiagram: Pymatgen PhaseDiagram class
    """

    phase_diagram = PhaseDiagram(entries=entries)

    return phase_diagram
