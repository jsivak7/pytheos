# for constructing phase diagrams

from pymatgen.analysis.phase_diagram import PDEntry, PhaseDiagram


def make_new_entry(composition: str, energy: float, name="my_PDEntry") -> PDEntry:
    """
    Makes new phase diagram entry (PDEntry) object. Useful for adding your own calculation
    to a queried entry list before generating a phase diagram.

    PDEntries contain all of the data needed to construct a PhaseDiagram.

    NOTE that you have to be very careful that calculation parameters match those you will be
    comparing this entry to!! Remember that MP applies different mixing schemes across their database
    - https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/thermodynamic-stability.

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


def generate(entries: list) -> PhaseDiagram:
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
