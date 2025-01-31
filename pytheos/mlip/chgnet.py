# for running CHGNet
# TODO needs to be updated - just was a function from a previous project I copied over


def run_chgnet_relax(struc_path, max_force=0.050):
    """
    General function to run CHGNet relaxation at 0K for a single structure file
    - NOTE that energy includes the MP2020Correction for Materials Project

    Args:
        struc_path (str): Structure file to relax
        max_force (float, optional): force convergence criteria for relaxation in eV/Angstrom
            - slightly lower than chgnet default of 100 meV/Angstrom

    Returns:
        tuple: (relaxed_structure, total energy in eV/atom, and list of magnetic moments)
    """

    from chgnet.model import StructOptimizer
    from chgnet.model.model import CHGNet
    from pymatgen.core import Structure
    from pymatgen.io.vasp import Poscar
    import pandas as pd
    import numpy as np

    structure = Structure.from_file(struc_path)

    chgnet = CHGNet.load()
    relaxer = StructOptimizer()

    result = relaxer.relax(
        structure, steps=10000, fmax=max_force
    )  # using a large number of steps that should never be reached
    final_structure = result[
        "final_structure"
    ]  # pmg structure object, can be written to .vasp file
    print("CHGNet relaxed structure", result["final_structure"])

    prediction = chgnet.predict_structure(final_structure)
    energy = float(prediction["e"])  # eV/atom
    magmoms = final_structure.site_properties["magmom"]  # list of magnetic moments

    return (final_structure, energy, magmoms)
