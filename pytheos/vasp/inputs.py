# for VASP input file generation

from ase import Atoms
import os
import yaml
from pymatgen.core import Structure
from pymatgen.io.vasp.sets import MPRelaxSet, MPScanRelaxSet
from pymatgen.io.vasp.inputs import Kpoints, VaspInput
import numpy as np
import random
from pytheos import utils


class CalcInputs:
    """
    Class for generating VASP calculation input files.

    The KSPACING tag is used instead of a KPOINTS file by default. A KPOINTS file can still be used by calling the `use_kpoints_file()` method.

    You will likely get a "BadInputSetWarning", but just ignore once you check you are using the correct POTCARs (see https://matsci.org/t/potcar-warning-msg/34562).

    Attributes:
        structure (Atoms): Initial structure as sorted Pymatgen Structure object.
        incar (INCAR): Pymatgen INCAR object.
        poscar (POSCAR): Pymatgen POSCAR object.
        potcar (POTCAR): Pymatgen POTCAR object.
    """

    def __init__(self, structure: Atoms, mp_input_set: str = "MPScanRelaxSet"):
        """
        Args:
            structure (Atoms): Initial structure as ASE Atoms.
            mp_input_set (str): Materials Project VASP input set for VASP calculations. Options are: "MPRelaxSet" or "MPScanRelaxSet". Defaults to "MPScanRelaxSet".

        Raises:
            Exception: If supplied `mp_input_set` has not been implemented.
        """

        # convert ASE Atoms to Pymatgen Structure and sort
        self.structure = Structure.from_ase_atoms(structure).sort()

        # get path to pytheos module + get customized pytheos MP set
        module_dir = os.path.dirname(__file__)
        try:
            with open(f"{module_dir}/input_sets/{mp_input_set.lower()}.yaml", "r") as f:
                incar_settings = yaml.load(f, Loader=yaml.SafeLoader)
        except:
            raise Exception(
                f"The input set you provided ({mp_input_set}) is not implemented in pytheos."
            )

        # make inputs based on specified MP input set
        if mp_input_set.lower() == "MPRelaxSet".lower():
            input_set = self._get_mprelaxset_inputs(incar_settings=incar_settings)
        elif mp_input_set.lower() == "MPScanRelaxSet".lower():
            input_set = self._get_mpscanrelaxset_inputs(incar_settings=incar_settings)

        self.incar = input_set.incar
        self.poscar = input_set.poscar
        self.potcar = input_set.potcar

    def update_incar(self, changes: dict) -> None:
        """
        Updates the generated INCAR for desired user changes.

        Args:
            changes (dict): changes for INCAR (ex. {"ALGO": "Normal"})

        Returns:
            None
        """

        self.incar.update(changes)

        return None

    def use_kpoints_file(self, kpoint_mesh: list) -> None:
        """
        Uses a KPOINTS file instead of KSPACING tag.

        Args:
            kpoint_mesh (list): k-points mesh - always Gamma centered to be safe (ex. [4, 4, 4])

        Returns:
            None
        """

        kpoints = Kpoints(kpts=kpoint_mesh)
        self.kpoints = kpoints

        self.incar.pop("KSPACING")  # to ensure both are not used in tandem

        return None

    def apply_mag_order(
        self,
        magmom_values: dict,
        mag_order_file: str = "magorder.yaml",
        rattle_amount: float = 0.5,
        coord_decimals: int = 1,
    ) -> list:
        """
        Applies magnetic ordering to structure using coordinates supplied as spin-up and spin-down in `mag_order_file`.

        Magnetic moments can be "rattled" some amount around the specified magnetic moment value.
        Can be beneficial to break the initial magnetic symmetry and assists in electronic convergence.

        NOTE that rattled values are between magmom_value +/- rattle_amount.

        Args:
            magmom_values (dict): Magnetic moment values for different elements in Bohr-Magneton.
            mag_order_file (str, optional): Relative path to .yaml file with desired magnetic ordering.
                Defaults to "magorder.yaml".
            rattle_amount (float, optional): Amount to rattle (randomly distort) magnetic moments in
                Bohr-Magneton. Defaults to 0.5.
            coord_decimals (int, optional): Number of decimals to round atom positions in Angstroms.
                Helpful for distorted structures and rounding errors. Defaults to 1.

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
            if atom.label != "O":  # only want cations
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
            print(f"\nThe number of spin-up and spin-down atoms are not equal.")
            print("Sometimes this is expected (e.g. in cation defect calculations).")
            utils.check_with_user()

        # update incar attribute with magnetic ordering
        self.incar.update({"MAGMOM": magmoms})

        return self.incar.as_dict()["MAGMOM"]

    def write_files(self, output_dir: str) -> None:
        """
        Writes generated VASP input files: INCAR, POSCAR, POTCAR, KPOINTS (if specified).

        Args:
            output_dir (str): Relative directory path to write input files. A new directory is created.

        Returns:
            None: VASP input files written to `output_dir`.
        """

        os.mkdir(output_dir)

        self.incar.write_file(f"{output_dir}/INCAR")
        self.poscar.write_file(f"{output_dir}/POSCAR")
        self.potcar.write_file(f"{output_dir}/POTCAR")

        try:
            self.kpoints.write_file(f"{output_dir}/KPOINTS")
        except AttributeError:
            pass

        return None

    def _get_mprelaxset_inputs(self, incar_settings: dict) -> VaspInput:
        """
        Gets VASP input sets for MPRelaxSet using Pymatgen.

        https://github.com/materialsproject/pymatgen/blob/master/src/pymatgen/io/vasp/MPRelaxSet.yaml

        Args:
            incar_settings (dict): INCAR settings for VASP run.

        Returns:
            VaspInput: Pymatgen objects for VASP inputs (INCAR, POSCAR, POTCAR)
        """

        input_set = MPRelaxSet(
            structure=self.structure,
            user_incar_settings=incar_settings,
            user_potcar_functional="PBE",
        ).get_input_set()

        return input_set

    def _get_mpscanrelaxset_inputs(self, incar_settings) -> VaspInput:
        """
        Gets VASP input sets for MPScanRelaxSet using Pymatgen.

        https://github.com/materialsproject/pymatgen/blob/master/src/pymatgen/io/vasp/MPSCANRelaxSet.yaml

        Args:
            incar_settings (dict): INCAR settings for VASP run.

        Returns:
            VaspInput: Pymatgen objects for VASP inputs (INCAR, POSCAR, POTCAR)
        """
        input_set = MPScanRelaxSet(
            structure=self.structure,
            user_incar_settings=incar_settings,
            user_potcar_functional="PBE_54",
        ).get_input_set()

        return input_set


def write_submission_script(
    job_name: str,
    output_dir: str,
    num_nodes: int = 1,
    num_cpu: int = 64,
    mem_per_cpu: str = "3500MB",
    runtime: int = 72,
    allocation: str = "sbs5563_bc",
) -> None:
    """
    Writes a SLURM submission script called "submitvasp".

    Assumed calculations are being run via a `cstdn.py` script.

    Specific for PennState RoarCollab HPC.

    Args:
        job_name (str): Name for job.
        output_dir (str): Relative directory to write submission script.
        num_nodes (int, optional): Number of nodes. Defaults to 1.
        num_cpu (int, optional): Number of CPUs per node. Defaults to 64.
        mem_per_cpu (str, optional): Amount of memory per CPU. Defaults to "3500MB".
        runtime (int, optional): Run time for job. Defaults to 72.
        allocation (str, optional): Allocation to run on. Defaults to "sbs5563_bc".

    Returns:
        None
    """

    slurm_script = f"""#!/bin/bash
#SBATCH --nodes={num_nodes}
#SBATCH --ntasks-per-node={num_cpu}
#SBATCH --mem-per-cpu={mem_per_cpu}
#SBATCH --time={runtime}:00:00
#SBATCH --partition=sla-prio
#SBATCH --account={allocation}
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --job-name={job_name}

export UCX_TLS=all

cd $SLURM_SUBMIT_DIR
module purge
module use /storage/icds/RISE/sw8/modules/
module load vasp/vasp-6.4.1v

eval "$(conda shell.bash hook)"
conda activate pytheos

python cstdn.py"""

    with open(f"{output_dir}/submitvasp", "w+") as f:
        f.writelines(slurm_script)

    return None
