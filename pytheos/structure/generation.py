# structure generation tools

from ase import Atoms
import os
from pymatgen.core import Structure
from icet import ClusterSpace
from icet.tools.structure_generation import (
    generate_sqs_from_supercells,
    _get_sqs_cluster_vector,
    occupy_structure_randomly,
)
from icet.input_output.logging_tools import set_log_config


class StructureGenerator:
    """
    Class for structure generation - usually used to make a supercell, generate special
    quasi-random structure (SQS), generate randomly decorated structure, or rattle atoms.

    NOTE that as methods are called, the `structure` attribute is updated.

    Attributes:
        structure (Atoms): ASE Atoms structure object.
    """

    def __init__(
        self,
        structure: Atoms,
    ):
        """
        Args:
            structure (Atoms): ASE Atoms structure object
                Usually want to start with a unit cell.
        """

        self.structure = structure

    def rattle_atoms(self, std_dev=0.02) -> Atoms:
        """
        Rattles atoms in structure using ASE - helpful prior to relaxation to break initial symmetry.

        Can enable symmetry-broken atomic arrangements during structure relaxation.
        - see Zunger group's publications on exploring a 'polymorphous representation' for more
        details of this approach

        Args:
            std_dev (float, optional): standard deviation of rattling amount to perform in Angstroms.
                Defaults to 0.02.

        Returns:
            Atoms: Rattled structure.
        """

        import random

        self.structure.rattle(
            std_dev,
            seed=int(random.uniform(0, 2000)),  # random seed
        )

        return self.structure

    def make_supercell(
        self,
        dimensions: list,
    ) -> Atoms:
        """
        Makes a supercell from structure using ASE.

        Args:
            dimensions (list): [x, y, z] cell multipliers for supercell generation.

        Returns:
            Atoms: Supercell structure.
        """

        supercell = self.structure.repeat(dimensions)

        self.structure = supercell

        return self.structure

    def make_sqs(
        self,
        dimensions: list,
        chemical_symbols,
        concentrations: dict,
        cutoffs: list,
        num_mc_steps: int = 10000,
    ):
        """
        Generates a special quasi-random structure (SQS) using ICET.
        - https://gitlab.com/materials-modeling/icet

        Args:
            dimensions (list): [x, y, z] cell multipliers for supercell generation.
            chemical_symbols (_type_): list of lists for allowed elements following same order as structure.
                e.g. [["Sr"], ["Ti", "Cr"], ["O"], ["O"], ["O"]] for a perovskite.
            concentrations (dict): Fractions of different elements for each lattice site.
                Only need to specify those that are not 1.
            cutoffs (list): cutoffs in order of multiplicity (pair, triplet, quadruplet, ...).
            num_mc_steps (int, optional): Number of Monte Carlo steps to run the SQS generation. Defaults to 10000.

        Returns:
            Atoms: SQS structure.
        """

        set_log_config(level="INFO")

        cluster_space = ClusterSpace(
            structure=self.structure,
            cutoffs=cutoffs,
            chemical_symbols=chemical_symbols,
        )
        print(cluster_space)

        sqs = generate_sqs_from_supercells(
            cluster_space=cluster_space,
            supercells=[self.structure.repeat(dimensions)],
            target_concentrations=concentrations,
            n_steps=num_mc_steps,
        )

        trial_cluster_vector = cluster_space.get_cluster_vector(sqs)
        perfectly_random_cluster_vector = _get_sqs_cluster_vector(
            cluster_space=cluster_space,
            target_concentrations=concentrations,
        )

        print(f"\nTrial Cluster Vector ->\n{trial_cluster_vector}")
        print(
            f"\nPerfectly Random Cluster Vector ->\n{perfectly_random_cluster_vector}"
        )

        self.structure = sqs

        return self.structure

    def decorate_randomly(
        self,
        dimensions: list,
        chemical_symbols: list,
        concentrations: dict,
    ):
        """
        Randomly decorates structure using ICET.
        - https://gitlab.com/materials-modeling/icet

        Args:
            dimensions (list): [x, y, z] cell multipliers for supercell generation.
            chemical_symbols (_type_): list of lists for allowed elements following same order as structure.
                e.g. [["Sr"], ["Ti", "Cr"], ["O"], ["O"], ["O"]] for a perovskite.
            concentrations (dict): Fractions of different elements for each lattice site.
                Only need to specify those that are not 1.

        Returns:
            Atoms: Randomly decorated structure.
        """

        cluster_space = ClusterSpace(
            structure=self.structure,
            cutoffs=[0],  # just need something for cluster space construction
            chemical_symbols=chemical_symbols,
        )

        randomly_decorated_structure = self.structure.repeat(dimensions)

        occupy_structure_randomly(
            structure=randomly_decorated_structure,
            cluster_space=cluster_space,
            target_concentrations=concentrations,
        )

        self.structure = randomly_decorated_structure

        return self.structure

    def write_structure(
        self,
        output_file: str,
        sort: bool = True,
        overwrite: bool = False,
    ) -> None:
        """
        Writes structure to file.

        Args:
            output_file (str, optional): Output file name - include path and file type. e.g. "../material.vasp".
            sort (bool, optional): Sort elements in by electronegativity using Pymatgen. Defaults to True.
            overwrite (bool, optional): Overwrite over already made file. Defaults to False.

        Raises:
            FileExistsError: If supplied file_name already exists to ensure that previously
                generated structures are not overwritten.

        Returns:
            None: Writes structure to file.
        """

        from ase.io import write

        if os.path.exists(output_file) and overwrite == False:
            raise FileExistsError(output_file)

        if sort == True:
            struc_pmg = Structure.from_ase_atoms(self.structure)
            struc_pmg.sort()
            self.structure = struc_pmg.to_ase_atoms()

        if "vasp" in output_file or "poscar" in output_file:
            write(output_file, self.structure, direct=True)

        else:
            write(output_file, self.structure)

        return None
