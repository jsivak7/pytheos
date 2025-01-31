# For creating FEFF input files and directories

from ase import Atoms


def write_xanes_inputs(
    structure: Atoms,
    absorbing_atom: int,
    output_dir: str,
    user_xanes_changes: dict = None,
):
    """
    Creates input directory and populates with associated files for running a XANES calculation in FEFF.
    Each FEFF calculation should be run in its own directory to ensure files are not overwritten.

    Args:
        structure (Atoms): ASE Atoms object for structure.
        absorbing_atom (int): Number in structure for absorbing atom of interest.
        output_dir (str): Relative path to create output directory with associated files.
        user_xanes_changes (dict, optional): Override for default card settings in 'xanes_cards.yaml'. Defaults to None.
            Example -> {"EDGE":"L", "FMS":"5.0 0"} changes this to an L-edge calculation with FMS radius of 5.0 Angstroms

    Raises:
        FileExistsError: If output directory already exists
    """
    from pymatgen.core import Structure
    from pymatgen.io.feff import sets
    import os
    import yaml

    print(f"Writing XANES input files to {output_dir}")

    if os.path.exists(output_dir):
        raise FileExistsError(output_dir)

    struc = Structure.from_ase_atoms(structure)

    # read in default cards from 'xanes_card.yaml' file
    module_dir = os.path.dirname(__file__)
    with open(f"{module_dir}/xanes_cards.yaml", "r") as f:
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


def write_feff_custodian_script(
    feff_cmd: str,
    output_dir: str = ".",
    max_errors=3,
):
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
