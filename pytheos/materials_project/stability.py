# evaluate stability using MP database
# see https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability

from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram


def make_phasediagram_entry(
    composition: str,
    energy: float,
    name="my_PDEntry",
) -> PDEntry:
    """
    Makes new phase diagram entry (PDEntry) object.

    Usually this is added to queried entries and used to generate PhaseDiagram and determine
    relative stability for your PDEntry.

    NOTE that you have to be very careful that calculation parameters match those you will be
    comparing this entry to!! Remember that MP applies different mixing schemes across their
    database -> https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/thermodynamic-stability.

    More info on PDEntries can be found at -> https://github.com/materialsproject/pymatgen/blob/master/src/pymatgen/analysis/phase_diagram.py.

    Args:
        composition (str): Material composition - can also provide Pymatgen Composition. Include all atoms. Example -> "MgCoNiO3".
        energy (float): Energy (in eV) corresponding to composition.
        name (str, optional): Name for entry - need some unique name for identification purposes. Defaults to "my_PDEntry".

    Returns:
        PDEntry: Pymatgen PDEntry object with all data needed to make a phase diagram.
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

    One can find the algorithm used in the parent PhaseDiagram class in Pymatgen.

    Args:
        entries (list): list of ComputedEntries/ComputedStructureEntries/PDEntries.

    Returns:
        PhaseDiagram: Pymatgen PhaseDiagram class
    """
    phase_diagram = PhaseDiagram(entries=entries)

    return phase_diagram
