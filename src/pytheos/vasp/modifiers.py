# for modifying already-run VASP calculations into other calculation types

from pytheos.vasp.outputs import load_vasprun
from pymatgen.io.vasp.inputs import Incar, Poscar, Potcar, Kpoints
from pymatgen.io.vasp.outputs import Eigenval
import os
from pytheos import utils


class CalcModifier:
    """
    Class to load previous VASP calculation files to modify them for subsequent, more specific calculations.

    Simpler modifications can also be performed by just changing `CalcModifier` attributes.
    - e.g. `CalcModifier.incar.update({"NCORE": 24})`

    To perform multiple modifications from a single source calc, be sure to reinitialize this class to
    ensure unintended changes were not applied. A sanity check is performed if the calculation type is
    not from the "source" to help avoid these mistakes

    Attributes:
        source_dir (str): relative path to source calc directory.
        calc_type (str): current calculation type - for monitoring progress. Starts with "source".
        vasprun (Vasprun): Vasprun object from source calc.
        incar (Incar): Incar object from source calc.
        poscar (Poscar): Relaxed structure (CONTCAR) from source calc as a Poscar object.
        potcar (Potcar): Potcar object from source calc.
        eigenval (Eigenval): Eigenval loaded from source calc.
        kpoints (Kpoints): Kpoints loaded from source calc. Only assigned if `KSPACING` does not exist in source INCAR.
        chgcar (bool): Only assigned by certain functions if needed for modifications.
        wavecar (bool): Only assigned by certain methods if needed for modifications.
    """

    def __init__(self, source_dir: str) -> None:
        """
        Args:
            source_dir (str): Relative path to source directory from which VASP files will be loaded.
        """

        self.source_dir: str = source_dir
        self.calc_type = "source"  # for monitoring

        # contains built-in convergence check (ionic & electronic)
        self.vasprun = load_vasprun(f"{source_dir}/vasprun.xml")

        self.incar = Incar.from_file(f"{source_dir}/INCAR")
        self.poscar = Poscar.from_file(f"{source_dir}/CONTCAR")  # relaxed structure
        self.potcar = Potcar.from_file(f"{source_dir}/POTCAR")
        self.eigenval = Eigenval(f"{source_dir}/EIGENVAL")

        # if `KSPACING` not found in INCAR, try to find KPOINTS
        if "KSPACING" not in self.incar:
            self.kpoints = Kpoints.from_file(f"{source_dir}/KPOINTS")

    def to_dos(
        self,
        incar_changes: dict = None,
        increase_nbands: float = 1.2,
        energy_window: float = 6.0,
    ) -> None:
        """
        Modifies VASP calculation objects for a density of states (DOS) calculation.

        Args:
            incar_changes (dict, optional): Additional changes to INCAR that the user can supply as a dictionary. Defaults to None.
            increase_nbands (float, optional): Factor to increase number of bands from source calculation.
                Helpful for DOS convergence. Defaults to 1.2.
            energy_window (float, optional): Energy window to calculation DOS from the Fermi level in +/- directions. Defaults to 6.0.
        """

        self._run_sanity_check()

        self.calc_type = "dos"

        # due to issue with fermi level placement
        num_elec = self.eigenval.nelect
        new_num_elec = num_elec - 1 + 0.999999

        # increase number of bands
        nbands = self.eigenval.nbands
        new_nbands = int(nbands * increase_nbands)

        # for a more reasonable energy window
        efermi = self.vasprun.efermi
        emin = efermi - energy_window
        emax = efermi + energy_window

        self.chgcar = True
        self.wavecar = True

        self.incar.update(
            {
                "NELECT": new_num_elec,  # per https://www.vasp.at/forum/viewtopic.php?t=17981
                "EMIN": emin,
                "EMAX": emax,
                "ISMEAR": -5,  # proper ISMEAR for DOS -> https://www.vasp.at/wiki/index.php/ISMEAR
                "ALGO": "Normal",  # since ALGO = All can often time fail with ISMEAR = -5
                "NSW": 0,  # static calculation
                "LCHARG": False,  # usually not needed
                "LWAVE": False,  # usually not needed
                "NEDOS": 501,  # higher resolution DOS than default
                "NBANDS": new_nbands,
            }
        )

        if incar_changes:
            self.incar.update(incar_changes)

    def to_bader(
        self,
        incar_changes: dict = None,
    ) -> None:
        """
        Modifies VASP calculation objects for a Bader charge calculation.
        - details on Bader charge analysis: https://theory.cm.utexas.edu/henkelman/code/bader/
        - for VASP, you need to run the `chgsum.pl` script prior to running `bader` (see in link above)

        Args:
            incar_changes (dict, optional): Additional changes to INCAR that the user can supply as a dictionary. Defaults to None.
        """

        self._run_sanity_check()

        self.calc_type = "bader"

        self.chgcar = True
        self.wavecar = True

        self.incar.update(
            {
                "ISMEAR": 0,  # Gaussian smearing
                "NSW": 0,  # static calculation (important for VASP Bader -> https://theory.cm.utexas.edu/henkelman/code/bader/)
                "LAECHG": True,  # needed for Bader charge analysis (https://www.vasp.at/wiki/index.php/LAECHG)
                "LCHARG": True,  # needed for Bader charge analysis
            }
        )

        if incar_changes:
            self.incar.update(incar_changes)

    def to_bandstructure(
        self,
        incar_changes: dict = None,
        increase_nbands: float = 2,
        sumo_kgen_cmd: str = "sumo-kgen --pymatgen --hybrid",
    ) -> None:
        """
        Modifies VASP calculation objects for a band structure calculation.
        - `sumo` package is used to generate the high-symmetry k-point path (https://github.com/SMTG-Bham/sumo)
        - remove `--hybrid` from the default sumo_kgen_cmd argument for GGA band structures
        - `sumo-bandplot` can be used to plot band structures (https://smtg-bham.github.io/sumo/sumo-bandplot.html)
        - `sumo-bandstats` can be used to evaluate band gaps and effective masses (https://smtg-bham.github.io/sumo/sumo-bandstats.html)

        Perform band structure calculations in the following steps ->
        1. run NSCF calculation (will be GGA for meta-GGA band structure)
        2. copy WAVECAR from finished NSCF calculation
        3. run BANDSTRUCTURE calculation

        Args:
            incar_changes (dict, optional): Additional changes to INCAR that the user can supply as a dictionary. Defaults to None.
            increase_nbands (float, optional): Factor to increase number of bands from source calculation. Defaults to 2.
            sumo_kgen_cmd (str, optional): Sumo command for high-symmetry k-point path generation. Defaults to "sumo-kgen --pymatgen".

        Raises:
            ValueError: if an invalid combination of INCAR flags and sumo_kgen_cmd is given
        """

        self._run_sanity_check()

        self.calc_type = "bandstructure"

        # increase number of bands
        nbands = self.eigenval.nbands
        new_nbands = int(nbands * increase_nbands)

        # check to see if calculation uses meta-GGA XC functional has proper sumo-kgen command
        if "METAGGA" in self.incar:
            is_metaGGA = True
            print(
                "Detected a meta-GGA calculation -> NSCF calculation will use `GGA = PE`"
            )

            self.ibzkpts = Kpoints.from_file(
                f"{self.source_dir}/IBZKPT"
            )  # need for meta-GGA high-symm kpoint generation

            if "hybrid" not in sumo_kgen_cmd:
                raise ValueError(
                    f"Cannot run meta-GGA band structure calculation without hybrid-type k-point path ({sumo_kgen_cmd})"
                )
        else:
            is_metaGGA = False
            if "hybrid" in sumo_kgen_cmd:
                raise ValueError(
                    f"Cannot run GGA band structure calculation with hybrid-type k-point path ({sumo_kgen_cmd})"
                )

        self.chgcar = True

        self.incar.update(
            {
                "NBANDS": new_nbands,
                "ISMEAR": 0,  # necessary for bandstructure
                "SIGMA": 0.001,  # small smearing value
                "NEDOS": 2001,  # adequate sampling
                "NSW": 0,  # static calculation
                "LCHARG": False,
                "LWAVE": True,
                "NEDOS": 501,  # higher resolution
                "LREAL": False,
            }
        )

        # for determining high-symmetry k-path
        self.kgen_cmd = sumo_kgen_cmd

        if incar_changes:
            self.incar.update(incar_changes)

    # NOTE that this has been less extensively tested/used than other methods!!
    def to_dielectric(
        self,
        incar_changes: dict = None,
        increase_nbands: float = 2,
    ) -> None:
        """
        Modifies VASP calculation objects for a dielectric calculation using the independent-particle approximation.
        - see https://www.vasp.at/wiki/index.php/Dielectric_properties_of_SiC for more information

        Args:
            incar_changes (dict, optional): Additional changes to INCAR that the user can supply as a dictionary. Defaults to None.
            increase_nbands (float, optional): Factor to increase number of bands from source calculation. Defaults to 2.
        """

        self._run_sanity_check()

        self.calc_type = "dielectric"

        self.chgcar = True
        self.wavecar = True

        self.incar.update(
            {
                "NBANDS": new_nbands,
                "LWAVE": False,
                "LCHARG": False,
                "NEDOS": 2000,  # for adequate sampling
                "SIGMA": 0.1,
                "ALGO": "Normal",
                "EDIFF": 1e-8,
                "LREAL": False,
                "LOPTICS": True,
                "CSHIFT": 0.01,
            }
        )

        if incar_changes:
            self.incar.update(incar_changes)

    def write_files(
        self,
        output_dir: str,
        copy_submit_script: bool = True,
    ) -> None:
        """
        Write current CalcModifier objects to VASP files.
        - A more complex scheme is used for 'bandstructure' type calculations due to need of NSCF calculation.

        Args:
            output_dir (str): Relative directory path to output files.
            copy_submit_script (bool): If want to copy same "submitvasp" from source calc. Defaults to True.
        """

        os.mkdir(output_dir)

        self.incar.write_file(f"{output_dir}/INCAR")
        self.poscar.write_file(f"{output_dir}/POSCAR")
        self.potcar.write_file(f"{output_dir}/POTCAR")
        if hasattr(self, "kpoints"):
            self.kpoints.write_file(f"{output_dir}/KPOINTS")

        # specific scheme for bandstructure calculations
        if self.calc_type == "bandstructure":

            self.ibzkpts.write_file(f"{output_dir}/IBZKPT")

            og_path = os.path.abspath("")  # to get back where we started
            os.chdir(output_dir)
            os.system(self.kgen_cmd)
            os.system("cp KPOINTS_band KPOINTS")
            os.chdir(og_path)

            # make separate dir for NSCF calculation
            os.mkdir(f"{output_dir}/nscf")

            if "METAGGA" in self.incar:
                self.incar.update(
                    {"GGA": "PE", "ICHARG": 11}
                )  # run nscf with GGA, ICHARG = 11
                self.incar.pop("METAGGA")  # remove METAGGA
            else:
                self.incar.update({"ICHARG": 11})  # run nscf with GGA, ICHARG = 11

            self.incar.write_file(f"{output_dir}/nscf/INCAR")
            self.poscar.write_file(f"{output_dir}/nscf/POSCAR")
            self.potcar.write_file(f"{output_dir}/nscf/POTCAR")
            if hasattr(self, "kpoints"):
                self.kpoints.write_file(f"{output_dir}/nscf/KPOINTS")

            # CHGCAR is necessary for NSCF calculation
            os.system(f"cp {self.source_dir}/CHGCAR {output_dir}/nscf/CHGCAR")

        # only copies over CHGCAR/WAVECAR if added as boolean attributes to CalcModifier class
        else:
            if hasattr(self, "chgcar"):
                os.system(f"cp {self.source_dir}/CHGCAR {output_dir}/CHGCAR")

            if hasattr(self, "wavecar"):
                os.system(f"cp {self.source_dir}/WAVECAR {output_dir}/WAVECAR")

        # only copy over `submitvasp` script if specified
        if copy_submit_script:
            os.system(f"cp {self.source_dir}/submitvasp {output_dir}/submitvasp")

            if self.calc_type == "bandstructure":
                os.system(
                    f"cp {self.source_dir}/submitvasp {output_dir}/nscf/submitvasp"
                )

        return None

    def add_chgcar(self) -> None:
        """Specify chgcar be copied to modified calculation."""
        self.chgcar = True

    def add_wavecar(self) -> None:
        """Specify wavecar be copied for modified calculation."""
        self.wavecar = True

    def _run_sanity_check(self) -> None:
        """Run a sanity check when CalcModifier methods are attempted for any calc_type other than the inputted 'source' calc"""

        if self.calc_type != "source":
            print("\nWARNING!!!")
            print("You are not modifying from the source calculation initially loaded.")
            print("You probably want to reinitialize a new CalcModifier...")

            utils.check_with_user()
