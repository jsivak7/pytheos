# for generating VASP input files

"""
Probably will get a `BadInputSetWarning` during POTCAR generation,
just be sure it matches the expected POTCAR type in warning and move on

Not 100% sure why this has been occuring, but has been noted by others on github
"""

from ase import Atoms
import os
import yaml
from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPRelaxSet, MPScanRelaxSet
from pymatgen.io.vasp.inputs import Kpoints
import numpy as np
import random
from pytheos import utils


class CalcInputs:
    """
    Class for VASP calculation input file generation INCAR, POSCAR, POTCAR, KPOINTS (if specified).

    The `KSPACING` INCAR flag is used for k-point generation instead of a KPOINTS file
    to faciliate high-throughput calculations. If a `kpoint_mesh` is supplied, a KPOINTS file can still
    be written, however, and the `KSPACING` tag in the INCAR will be removed.

    Attributes:
        structure (Atoms): Initial structure as ASE Atoms.
        mp_input_set (str): Materials Project input set to describe VASP calculation inputs. Currently
            only implemented (& customized) for "MPRelaxSet" and "MPScanRelaxSet". Defaults to "MPScanRelaxSet".
        incar_changes (dict, optional): Specific changes to the INCAR file. Defaults to None.
        kpoint_mesh (tuple, optional): K-point mesh (if desired over `KSPACING` in INCAR. Always Gamma-centered to be safe.
                `KSPACING` is removed from INCAR if this is modified. Defaults to None.
        incar (Incar): Pymatgen INCAR object.
        poscar (Poscar): Pymatgen POSCAR object.
        potcar (Potcar): Pymatgen POTCAR object.
        kpoints (Kpoints): Pymatgen KPOINTS object (if `kpoint_mesh` is specified)
    """

    def __init__(
        self,
        structure: Atoms,
        mp_input_set: str = "MPScanRelaxSet",
        incar_changes: dict = None,
        kpoint_mesh: tuple = None,
    ):
        """
        Args:
            structure (Atoms): Initial structure as ASE Atoms.
            mp_input_set (str): Materials Project VASP input set to describe calculation inputs. Currently
                only implemented (& customized) for "MPRelaxSet" and "MPScanRelaxSet". Defaults to "MPScanRelaxSet".
            incar_changes (dict, optional): Specific changes to the INCAR file. Defaults to None.
            kpoint_mesh (tuple, optional): K-point mesh (if desired over `KSPACING` in INCAR. Gamma-centered.
                `KSPACING` is removed from INCAR if this is modified. Defaults to None.

        Raises:
            Exception: If supplied `mp_input_set` has not been implemented.
        """

        self.structure = structure
        self.mp_input_set = mp_input_set
        self.incar_changes = incar_changes
        self.kpoint_mesh = kpoint_mesh

        # convert ASE Atoms -> PMG Structure (sorted)
        self.structure = Structure.from_ase_atoms(self.structure).sort()

        # get path to this module + get customized pytheos MP set
        module_dir = os.path.dirname(__file__)

        try:
            with open(f"{module_dir}/{self.mp_input_set.lower()}.yaml", "r") as f:
                incar_settings = yaml.load(f, Loader=yaml.SafeLoader)
        except:
            raise Exception(
                f"The input set you provided ({self.mp_input_set}) is not implemented in pytheos.\nCurrent options are: MPRelaxSet (pbe) or MPScanRelaxSet (r2scan)."
            )

        if self.incar_changes:
            incar_settings.update(self.incar_changes)

        kpoint_settings = None
        if self.kpoint_mesh:
            kpoint_settings = Kpoints(kpts=self.kpoint_mesh)
            incar_settings.update({"KSPACING": None})

        # make input files based on specified input set
        if mp_input_set.lower() == "MPRelaxSet".lower():
            input_set = MPRelaxSet(
                structure=self.structure,
                user_incar_settings=incar_settings,
                user_potcar_functional="PBE",
                user_kpoints_settings=kpoint_settings,
                sort_structure=True,
            ).get_input_set()

        elif mp_input_set.lower() == "MPScanRelaxSet".lower():
            input_set = MPScanRelaxSet(
                structure=self.structure,
                user_incar_settings=incar_settings,
                user_potcar_functional="PBE_54",
                user_kpoints_settings=kpoint_settings,
                sort_structure=True,
            ).get_input_set()

        self.incar = input_set.incar
        self.poscar = input_set.poscar
        self.potcar = input_set.potcar

        if kpoint_mesh:
            self.kpoints = input_set.kpoints

    def apply_mag_order(
        self,
        magmom_values: dict,
        mag_order_file: str = "magorder.yaml",
        rattle_amount: float = 0.5,
        coord_decimals: int = 1,
        anion: str = "O",
    ) -> list:
        """
        Applies magnetic ordering to structure using coordinates supplied as spin-up and spin-down in `mag_order_path`.

        Magnetic moments can be "rattled" some amount around the specified magnetic moment value.
        Can be beneficial to break the initial magnetic symmetry and assists in electronic convergence.

        NOTE that rattled values are between magmom_value +/- rattle_amount

        Args:
            magmom_values (dict): Magnetic moment values for different elements in Bohr-Magneton.
            mag_order_file (str, optional): Relative path to .yaml file with desired magnetic ordering.
                Defaults to "magorder.yaml".
            rattle_amount (float, optional): Amount to rattle (randomly distort) magnetic moments in
                Bohr-Magneton. Defaults to 0.5.
            coord_decimals (int, optional): Number of decimals to round atom positions in Angstroms.
                Helpful for distorted structures and rounding errors. Defaults to 1.
            anion (str, optional): Anion in structure. Defaults to "O".

        Raises:
            Exception: If atomic coordinates cannot be identified as spin-up or spin-down.

        Returns:
            list: Magnetic moments for updated `MAGMOM` in `CalcInputs.incar`.
        """

        # load magnetic ordering positions from yaml
        with open(f"{mag_order_file}", "r") as f:
            mag_order_coords = yaml.load(f, Loader=yaml.SafeLoader)

        # split supplied magnetic order positions into spin up and spin down
        spin_up_coords = np.round(np.array(mag_order_coords["spin up"]), coord_decimals)
        spin_down_coords = np.round(
            np.array(mag_order_coords["spin down"]), coord_decimals
        )

        magmoms = []  # as will be added to INCAR
        spin_counts = {"up": 0, "down": 0}
        atom_index = 0

        for atom in self.structure:

            # round atom coordinates
            atom_coords = np.round(atom.coords, coord_decimals)

            # get magmom value from supplied dictionary and rattle it
            magmom_value = magmom_values[atom.label]
            magmom_min_value = magmom_value - rattle_amount
            magmom_max_value = magmom_value + rattle_amount
            magmom_value_rattled = np.round(
                random.uniform(magmom_min_value, magmom_max_value), 2
            )

            # only apply ordering to cations
            if atom.label != anion:
                coord_identified = False

                # try to find in spin-up coordinates
                for spin_up_coord in spin_up_coords:
                    if np.array_equal(atom_coords, spin_up_coord):
                        coord_identified = True
                        spin_counts["up"] += 1
                        magmom_value_rattled = +magmom_value_rattled

                # try to find in spin-down coordinates if not found in spin-up
                if coord_identified == False:
                    for spin_down_coord in spin_down_coords:
                        if np.array_equal(atom_coords, spin_down_coord):
                            coord_identified = True
                            spin_counts["down"] += 1
                            magmom_value_rattled = -magmom_value_rattled

                # if we still cannot ID coordinates
                if coord_identified == False:
                    raise Exception(
                        f"Atom Index {atom_index} ({atom.label}) spin could not be identified."
                    )

            magmoms.append(magmom_value_rattled)

            atom_index += 1

        # check that number of spin-up atoms = spin-down atoms
        # can pass with user input if this is the expected behavior
        if spin_counts["up"] != spin_counts["down"]:
            print("\nWARNING!!!")
            print(
                f"\nThe number of spin-up ({spin_counts["up"]}) and spin-down ({spin_counts["down"]}) atoms are not equal."
            )
            print("Sometimes this is expected (e.g. in cation defect calculations).")
            utils.check_with_user()

        # update incar attribute with magnetic ordering
        self.incar.update({"MAGMOM": magmoms})

        return self.incar.as_dict()["MAGMOM"]

    def write_files(
        self,
        output_dir: str,
    ) -> None:
        """
        Writes generated VASP input files: INCAR, POSCAR, POTCAR, KPOINTS (if specified).

        Args:
            output_dir (str): Relative directory path to write input files.
                A new directory is created.

        Returns:
            None: VASP input files written to `output_dir`.
        """

        os.mkdir(output_dir)

        self.incar.write_file(f"{output_dir}/INCAR")
        self.poscar.write_file(f"{output_dir}/POSCAR")
        self.potcar.write_file(f"{output_dir}/POTCAR")

        if self.kpoints:
            self.kpoints.write_file(f"{output_dir}/KPOINTS")

        return None
