# module to facilitate FEFF x-ray absorption spectra (XAS) simulations
# see https://feff.phys.washington.edu/index-feffproject.html

from pandas import DataFrame
from pymatgen.core import Structure
from pymatgen.io.feff import sets
import os
import yaml
from pymatgen.io.feff import outputs
import pandas as pd


def write_xanes_inputs(
    struc: Structure,
    absorbing_atom: int,
    output_dir: str,
    user_xanes_changes: dict = None,
) -> None:
    """
    Creates input directory and populates with associated files for running a XANES calculation in FEFF.
    Each FEFF calculation should be run in its own directory to ensure files are not overwritten.

    Args:
        struc (Structure): Pymatgen Structure object.
        absorbing_atom (int): Number in structure for absorbing atom of interest.
        output_dir (str): Relative path to create output directory with associated files.
        user_xanes_changes (dict, optional): Override for default card settings in 'xanes_cards.yaml'. Defaults to None.
            Example -> {"EDGE":"L", "FMS":"5.0 0"} changes this to an L-edge calculation with FMS radius of 5.0 Angstroms

    Raises:
        FileExistsError: If output directory already exists
    """

    print(f"\nSetting up FEFF XANES inputs...")

    if os.path.exists(output_dir):
        raise FileExistsError(output_dir)

    print(f"-> structure = {struc.composition.formula}")
    print(
        f"-> absorbing atom = #{absorbing_atom}-{struc[absorbing_atom].species.reduced_formula}"
    )
    print(f"-> output directory = {output_dir}")
    print(f"-> user xanes changes = {user_xanes_changes}\n")

    # read in default cards from 'xanes_card.yaml' file
    module_dir = os.path.dirname(__file__)
    with open(f"{module_dir}/input_sets/xanes.yaml", "r") as f:
        xanes_settings = yaml.load(f, Loader=yaml.SafeLoader)

    # only update FEFF.inp cards if user asks for changes to be made
    if bool(user_xanes_changes) != False:
        xanes_settings.update(user_xanes_changes)

    inputs = sets.MPXANESSet(
        structure=struc,
        absorbing_atom=absorbing_atom,
        user_tag_settings=xanes_settings,
    )
    inputs.write_input(output_dir=output_dir)


def write_custodian_script(
    feff_cmd: str,
    output_dir: str = ".",
    max_errors=3,
) -> None:
    """
    Writes a python script that can be used to run FEFF using Custodian (https://github.com/materialsproject/custodian).

    This can be used for any FEFF calculation (XANES, EXAFS, etc.).

    One should preferentially use this to run FEFF over just utilizing the command line
    as it automatically backs up files and handles possible errors in a consistent, automated format.

    NOTE that the written script needs to be run in the same directory as the 'feff.inp'.

    Args:
        feff_cmd (str): Command to run FEFF software.
        output_dir (str, optional): Relative path to create 'cstdn_feff.py' script. Defaults to ".".
        max_errors (int, optional): Max number of errors for Custodian to try to correct. Defaults to 3.
    """

    print(f"Writing FEFF Custodian script to {output_dir}")

    custodian_script = f"""# runs the FEFF10 software using Custodian package.

from custodian.custodian import Custodian
from custodian.feff import jobs, handlers

feff_job = jobs.FeffJob(feff_cmd="{feff_cmd}")
c = Custodian(
    handlers=[handlers.UnconvergedErrorHandler()],
    jobs=[feff_job],
    max_errors={max_errors},
)
c.run()"""

    with open(f"{output_dir}/cstdn_feff.py", "w") as write_custodian:
        write_custodian.write(custodian_script)


def xmu_dat_to_df(
    xmu_path: str = "xmu.dat",
    feff_inp_path: str = "feff.inp",
) -> DataFrame:
    """
    Converts the xmu.dat output from FEFF to a Pandas dataframe object.

    Args:
        xmu_path (str, optional): Relative path to xmu.dat data file. Defaults to "xmu.dat".
        feff_inp_path (str, optional): Relative path to feff.inp file. Defaults to "feff.inp".

    Returns:
        pd.DataFrame: Pandas dataframe object of FEFF output data from xmu.dat.
    """

    print(f"Reading {xmu_path} as Pandas DataFrame")

    xmu = outputs.Xmu.from_file(xmu_dat_file=xmu_path, feff_inp_file=feff_inp_path)

    edge = xmu.edge
    absorbing_atom = xmu.absorbing_atom
    material = xmu.material_formula
    calc_type = xmu.calc

    xmu_df = pd.DataFrame(
        {
            "energy_eV": xmu.energies,
            "wavenumber": xmu.wavenumber,
            "mu": xmu.mu,
            "mu0": xmu.mu0,
            "chi": xmu.chi,
        }
    )
    return xmu_df
