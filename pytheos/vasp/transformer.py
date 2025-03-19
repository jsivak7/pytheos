# for transforming already-run VASP calculations into other calculation types

from pytheos.vasp.outputs import load_vasprun
from pymatgen.io.vasp.inputs import Incar, Poscar, Potcar, Kpoints
from pymatgen.io.vasp.outputs import Eigenval, Chgcar, Wavecar
import os


class Transformer:
    """
    Loads previous VASP calculation files as objects so that they can be transformed into subsequent calculations.
    - Some commonly used transformations are provided as methods with reasonable default changes.
    - Simple transformations can also be performed by just changing this object's attributes.
        - ex. Transformer.incar.update({"NCORE": 24})

    If one wishes to perform multiple conversions from a single source calculation, this should be reinitialized following each conversion to ensure unintended changes were not applied.
    - a sanity check is performed if the calc type is not from the "source" to help avoid mistakes

    Attributes:
        source_dir (str): relative path to source calc directory.
        type (str): current calculation type - for monitoring progress. Starts with "source".
        vasprun (Vasprun): Vasprun object from source calc.
        incar (Incar): Incar object from source calc.
        poscar (Poscar): Relaxed structure (CONTCAR) from source calc as a Poscar object.
        potcar (Potcar): Potcar object from source calc.
        eigenval (Eigenval): Eigenval loaded from source calc.
        kpoints (Kpoints): Kpoints loaded from source calc. Only assigned if `KSPACING` does not exist in source INCAR.
        chgcar (Chgcar): Chgcar loaded from source calc. Only assigned by certain functions if needed for output calc.
        wavecar (Wavecar): Wavecar loaded from source calc. Only assigned by certain functions if needed for output calc.
    """

    def __init__(self, source_dir: str) -> None:
        """
        Args:
            source_dir (str): Relative path to source directory where VASP files will be loaded.
        """
        self.source_dir: str = source_dir
        print(f"Loading 'source' calc from ./{self.source_dir}/")
        self.type = "source"  # for monitoring

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
        user_incar_changes: dict = None,
        energy_window: float = 6.0,
    ) -> None:
        """
        Transforms VASP calculation objects for density of states (DOS) calculation.

        Args:
            user_incar_changes (dict, optional): Additional changes to INCAR that the user can supply as a dictionary - example. Defaults to None.
            energy_window (float, optional): Energy window to calculation DOS from the Fermi level in +/- directions. Defaults to 6.0.
        """

        run_sanity_check(self.type)

        print(f"Transforming '{self.type}' calc --> 'dos' calc")
        self.type = "dos"

        # due to issue with fermi level placement
        num_elec = self.eigenval.nelect
        new_num_elec = num_elec - 1 + 0.999999

        # for a more reasonable energy window
        efermi = self.vasprun.efermi
        emin = efermi - energy_window
        emax = efermi + energy_window

        # chgcar and wavecar should be read in
        self.chgcar = Chgcar.from_file(f"{self.source_dir}/CHGCAR")
        self.wavecar = Wavecar(f"{self.source_dir}/WAVECAR")

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
            }
        )

        # check if user gave any additional changes
        if bool(user_incar_changes) == True:
            incar.update(user_incar_changes)

    def to_bader(
        self,
        user_incar_changes: dict = None,
    ) -> None:
        """
        Transforms VASP calculation objects for Bader charge calculation.

        Details on Bader charge analysis: https://theory.cm.utexas.edu/henkelman/code/bader/

        Args:
            user_incar_changes (dict, optional): Additional changes to INCAR that the user can supply as a dictionary. Defaults to None.
        """

        run_sanity_check(self.type)

        print(f"Transforming '{self.type}' calc --> 'bader' calc")
        self.type = "bader"

        # chgcar and wavecar should be read in
        self.chgcar = Chgcar.from_file(f"{self.source_dir}/CHGCAR")
        self.wavecar = Wavecar(f"{self.source_dir}/WAVECAR")

        self.incar.update(
            {
                "ISMEAR": 0,  # Gaussian smearing
                "NSW": 0,  # static calculation (important for VASP Bader -> https://theory.cm.utexas.edu/henkelman/code/bader/)
                "LAECHG": True,  # needed for Bader charge analysis (https://www.vasp.at/wiki/index.php/LAECHG)
                "LCHARG": True,  # needed for Bader charge analysis
            }
        )

        # check if user gave any additional changes
        if bool(user_incar_changes) == True:
            incar.update(user_incar_changes)

    def to_bandstructure(
        self,
        user_incar_changes: dict = None,
        increase_nbands: float = 2,
        sumo_kgen_cmd: str = "sumo-kgen --pymatgen",
    ) -> None:
        """
        Transforms VASP calculation objects for bandstructure calculation.

        The `sumo` package is used to generate the high-symmetry k-point path -> https://github.com/SMTG-Bham/sumo

        One should add the `--hybrid` from the default sumo_kgen_cmd argument for metaGGA band structures.

        Band structures can be plotted with the `sumo-bandplot` command -> https://smtg-bham.github.io/sumo/sumo-bandplot.html

        Band gaps and effective masses can be extracted with the `sumo-bandstats` command -> https://smtg-bham.github.io/sumo/sumo-bandstats.html

        Args:
            user_incar_changes (dict, optional): Additional changes to INCAR that the user can supply as a dictionary. Defaults to None.
            increase_nbands (float, optional): Factor to increase number of bands from source calculation. Defaults to 2.
            sumo_kgen_cmd (str, optional): Sumo command for high-symmetry k-point path generation. Defaults to "sumo-kgen --pymatgen".

        Raises:
            ValueError: if an invalid combination of INCAR flags and sumo_kgen_cmd is given
        """

        run_sanity_check(self.type)

        print(f"Transforming '{self.type}' calc --> 'bandstructure' calc")
        self.type = "bandstructure"

        # increase number of bands
        nbands = self.eigenval.nbands
        new_nbands = int(nbands * increase_nbands)

        # check to see if calculation uses meta-GGA XC functional has proper sumo-kgen command
        if "METAGGA" in self.incar:
            is_metaGGA = True
            print(
                "Detected a meta-GGA calculation... NSCF calculation will use `GGA = PE`"
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

        # chgcar and wavecar should be read into calc
        self.chgcar = Chgcar.from_file(f"{self.source_dir}/CHGCAR")
        self.wavecar = Wavecar(f"{self.source_dir}/WAVECAR")

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

        # check if user gave any additional changes
        if bool(user_incar_changes) == True:
            incar.update(user_incar_changes)

    def write_files(
        self,
        output_dir: str,
    ) -> None:
        """
        Write current Transformer to VASP files.
        - A more complex scheme is used for 'bandstructure' type calculations due to need of NSCF calculation.

        Args:
            output_dir (str): Relative directory path to output files.
        """

        print(f"Writing '{self.type}' calc --> ./{output_dir}/")

        os.mkdir(output_dir)

        self.incar.write_file(f"{output_dir}/INCAR")
        self.poscar.write_file(f"{output_dir}/POSCAR")
        self.potcar.write_file(f"{output_dir}/POTCAR")
        if hasattr(self, "kpoints"):
            self.kpoints.write_file(f"{output_dir}/KPOINTS")

        # specific scheme for bandstructure calculations
        if self.type == "bandstructure":
            print(
                f"""\nPerform bandstructure calculations in the following steps...
    1. run the ./{output_dir}/nscf calculation
    2. copy the WAVECAR from finished ./{output_dir}/nscf calculation to ./{output_dir}
    3. run the ./{output_dir} calculation\n"""
            )

            self.ibzkpts.write_file(f"{output_dir}/IBZKPT")

            og_path = os.path.abspath(".")  # to get back where we started
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

        # only copies over CHGCAR/WAVECAR if added as attributes to Transformer object
        else:
            if hasattr(self, "chgcar"):  # copied, not written
                os.system(f"cp {self.source_dir}/CHGCAR {output_dir}/CHGCAR")
            if hasattr(self, "wavecar"):  # copied, not written
                os.system(f"cp {self.source_dir}/WAVECAR {output_dir}/WAVECAR")

    # NOTE that this has been less extensively tested/used than other methods!!
    def to_dielectric(
        self,
        user_incar_changes: dict = None,
        increase_nbands: float = 2,
    ) -> None:
        """
        Transforms VASP calculation objects for dielectric calculation.

        Args:
            user_incar_changes (dict, optional): Additional changes to INCAR that the user can supply as a dictionary. Defaults to None.
            increase_nbands (float, optional): Factor to increase number of bands from source calculation. Defaults to 2.
        """

        run_sanity_check(self.type)

        print(f"Transforming '{self.type}' calc --> 'dielectric' calc")
        self.type = "dielectric"

        # chgcar and wavecar should be read in
        self.chgcar = Chgcar.from_file(f"{self.source_dir}/CHGCAR")
        self.wavecar = Wavecar(f"{self.source_dir}/WAVECAR")

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

        # check if user gave any additional changes
        if bool(user_incar_changes) == True:
            incar.update(user_incar_changes)

    def load_chgcar(self) -> None:
        """Method to load chgcar for unique transformations"""
        self.chgcar = Chgcar.from_file(f"{self.source_dir}/WAVECAR")

    def load_wavecar(self):
        """Method to load wavecar for unique transformations"""
        self.wavecar = Wavecar(f"{self.source_dir}/WAVECAR")


def run_sanity_check(type) -> None:
    """Run a sanity check when transformer methods are attempted for any type other than the inputted 'source' calc"""

    if type != "source":

        print(
            "WARNING!!!\nYou are not transforming from the 'source' calc that was initially loaded."
        )

        while True:
            answer = input("Do you want to continue? (yes/no): ").lower()
            if answer in ["y", "yes"]:
                print("Continuing...")
                break
            elif answer in ["n", "no"]:
                print("Aborting...")
                exit()
                break
            else:
                print("Invalid input. Please enter 'yes' or 'no'.")
