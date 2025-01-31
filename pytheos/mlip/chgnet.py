# for running CHGNet simulations

from ase import Atoms


def run_chgnet_relax(
    structure: Atoms,
    max_force=0.050,
    # TODO can you run constrained relaxations (i.e. like ISIF = 8 in VASP?)
):
    """
    For running CHGNet relaxation at 0K for a single structure file
    - NOTE that resulting energy includes the MP2020Compatibility corrections from Materials Project
        - (https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/thermodynamic-stability/anion-and-gga-gga+u-mixing)

    Args:
        structure (Atoms): ASE Atoms object for structure.
        max_force (float, optional): Force convergence criteria for relaxation in eV/Angstrom.
            - Slightly lower than CHGNet default of 100 meV/Angstrom, good middle ground

    Returns:
        tuple: (Atoms object of relaxed structure, total energy in eV/atom w/ MP2020Compatibility corrections, and list of magnetic moments)
    """

    from chgnet.model import StructOptimizer
    from chgnet.model.model import CHGNet
    from pymatgen.core import Structure

    print(f"Running CHGNet relaxation.")

    s = Structure.from_ase_atoms(structure)

    chgnet = CHGNet.load()
    relaxer = StructOptimizer()

    # using a large number of steps to ensure force convergence criteria is reached
    # seems that for HEOs sometimes the default of 500 can be reached
    result = relaxer.relax(s, steps=10000, fmax=max_force)

    # Pymatgen structure object, can be written to .vasp file
    s_final_pmg = result["final_structure"]
    print(f"CHGNet relaxed structure -->\n{s_final_pmg}")
    s_final = s_final_pmg.to_ase_atoms()

    # to get predictions from CHGNet on relaxed structure
    prediction = chgnet.predict_structure(final_structure)

    # eV/atom with MP2020Compatibility corrections
    energy = float(prediction["e"])

    # same order as inputted structure
    magmoms = final_structure.site_properties["magmom"]

    return (s_final, energy, magmoms)
