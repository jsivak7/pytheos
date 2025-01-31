# for analysis of VASP results


def calc_form_decomp_energy(
    composition: str, total_energy: float, MPApiKey: str, xc="R2SCAN"
):
    """
    Calculates the formation and decomposition energy for a single VASP calculation using Materials Project
    - calculation parameters should be ~equivalent with MP(Scan)RelaxSet to be valid

    This could also be done using the vasprun.xml file (probably a bit simpler), but opted for this as it is more flexible

    Args:
        composition (str): Chemical composition ("Mg4Ni4O8", or "Sr2Ti1Nb1O6")
        total_energy (float): total energy corresponding to the composition
        MPApiKey (str): Materials Project API Key (https://next-gen.materialsproject.org/api)
        xc (str): exchange correlation functional used. Defaults to "R2SCAN".
            - other current option is "GGA/GGA+U"

    Returns:
        tuple: (e_form, e_decomp, decomp_formula)
    """
    from pymatgen.analysis.phase_diagram import PDEntry
    from mp_api.client import MPRester
    from pymatgen.analysis.phase_diagram import PhaseDiagram
    from emmet.core.thermo import ThermoType

    print(
        f"\n========== Calculating {xc} ΔHf and ΔHd from Materials Project =========="
    )

    print(f"composition:\t\t{composition}")
    print(f"total energy:\t\t{total_energy} eV")
    system = PDEntry(
        composition=composition,
        energy=total_energy,
        name="my_calc",  # need some identifying name
    )

    # get all elements from the run and places them into a list with the required format
    elements = system.composition.elements
    elements_newlist = []
    for element in range(len(elements)):
        elements_newlist.append(elements[element].symbol)

    # get all entries from Materials Project that contain the elements of interest
    if xc == "R2SCAN":
        with MPRester(MPApiKey) as mpr:
            entries = mpr.get_entries_in_chemsys(
                elements_newlist,
                additional_criteria={"thermo_types": [ThermoType["R2SCAN"]]},
            )
    elif xc == "GGA/GGA+U":
        with MPRester(MPApiKey) as mpr:
            entries = mpr.get_entries_in_chemsys(elements_newlist)

    entries.append(system)
    phasediagram = PhaseDiagram(entries)

    all_entries = phasediagram.all_entries

    for e in all_entries:
        if e.name == "my_calc":  # to only get run
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

    print(f"formation energy:\t{e_form:.3f} eV/atom")
    print(f"decomposition energy:\t{e_decomp:.3f} eV/atom")
    print(f"decomposition rxn:\t{decomp_formula}")
    return (e_form, e_decomp, decomp_formula)
