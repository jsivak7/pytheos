# Machine Learning Interatomic Potentials

from ase import Atoms
from chgnet.model.dynamics import StructOptimizer
from chgnet.model.model import CHGNet
from pymatgen.core import Structure


def relax_with_OG_chgnet_MPTrj(
    structure: Atoms,
    max_force=0.050,
    optimizer="FIRE",
):
    """
    For running Crystal Hamiltonian Graph Neural Network (CHGNet) relaxation at 0K for an inputted ASE Atoms object
    - NOTE that resulting energy includes the MP2020Compatibility corrections from Materials Project
        - (https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/thermodynamic-stability/anion-and-gga-gga+u-mixing)

    Args:
        structure (Atoms): ASE Atoms object for structure.
        max_force (float, optional): Force convergence criteria for relaxation in eV/Angstrom.
            - Slightly lower than CHGNet default of 100 meV/Angstrom, good middle ground
        optimizer (string): ASE optimizer to use. Defaults to "FIRE"
            - Sometimes faster relaxation can be achieved with other optimizers like "BFGS" - do your own testing!
            - https://wiki.fysik.dtu.dk/ase/ase/optimize.html
            - https://gpaw.readthedocs.io/devel/ase_optimize/ase_optimize.html

    Returns:
        dict: results = {
            "final_structure": s_final, # ASE Atoms object
            "energy_eVperAtom": energy, # includes MP2020Compatibility corrections
            "magmoms": magmoms, # list of magnetic moments following structure ordering
        }
    """

    print(
        f"Running CHGNet relaxation with force criteria = {max_force} meV/Angstrom..."
    )

    s = Structure.from_ase_atoms(structure)

    chgnet = CHGNet.load()
    relaxer = StructOptimizer(optimizer_class=optimizer)

    # using a large number of steps to ensure force convergence criteria is reached
    # seems that for HEOs sometimes the default of 500 can be reached
    result = relaxer.relax(s, steps=10000, fmax=max_force)

    # Pymatgen structure object, can be written to .vasp file
    s_final_pmg = result["final_structure"]
    print(f"\nCHGNet relaxed structure -->\n\n{s_final_pmg}")
    s_final = s_final_pmg.to_ase_atoms()

    # to get predictions from CHGNet on relaxed structure
    prediction = chgnet.predict_structure(s_final_pmg)

    # eV/atom with MP2020Compatibility corrections
    energy = float(prediction["e"])

    # same order as inputted structure
    magmoms = s_final_pmg.site_properties["magmom"]

    results = {
        "final_structure": s_final,
        "energy_eVperAtom": energy,
        "magmoms": magmoms,
    }

    return results
