# for energy correction schemes from the Materials Project

from pymatgen.io.vasp.outputs import Vasprun


def apply_mp2020compat(run: Vasprun) -> float:
    """
    Applies MP2020Compatbility correction scheme for GGA/GGA+U and anion mixing calculations.
    - https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/thermodynamic-stability/anion-and-gga-gga+u-mixing

    Calculation parameters/potcars should be consistent with MPRelaxSet for valid computations.
    - https://github.com/materialsproject/pymatgen/blob/master/src/pymatgen/io/vasp/MPRelaxSet.yaml

    Args:
        run (Vasprun): Pymatgen vasprun object. Used preferentially over raw energies to ensure
            scheme is implemented correctly for a given material system and calculation specs.

    Returns:
        float: Corrected energy in eV/atom
    """
    from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
    import numpy as np

    v = run.get_computed_entry()

    # get original energy in eV/atom
    energy_og = v.energy / len(v.structxure)
    print(f"original energy = {np.round(energy_og, 4)}/atom")

    # calculate corrected energy with MP2020Compatibility corrections
    v.energy_adjustments = MaterialsProject2020Compatibility().get_adjustments(v)
    energy_mp2020 = v.energy / len(v.structure)
    print(f"corrected energy = {np.round(energy_mp2020, 4)}/atom")

    return energy_mp2020
