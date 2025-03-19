# evaluate stability using MP database
# see https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability

class PhaseDiagram:
    


def generate_phase_diagram(
    MPApiKey: str,
    elements: tuple,
    xc: str = "GGA_GGA_U",
):

    from mp_api.client import MPRester
    from emmet.core.thermo import ThermoType
    from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter

    with MPRester(MPApiKey) as mpr:

        # Obtain only corrected GGA and GGA+U ComputedStructureEntry objects
        entries = mpr.get_entries_in_chemsys(
            elements=elements,
            additional_criteria={"thermo_types": [ThermoType[xc]]},
        )
        # Construct phase diagram
        pd = PhaseDiagram(entries)

        # Plot phase diagram
        plot = PDPlotter(pd).get_plot()
