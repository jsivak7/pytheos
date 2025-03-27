# custodian script writer functions for vasp calculations
# https://github.com/materialsproject/custodian


def write_script(
    output_dir: str = "./",
    vasp_cmd: tuple = ("srun", "vasp_std"),
    max_errors: int = 3,
) -> None:
    """
    Writes a generic custodian script for VASP calculation.

    Args:
        output_dir (str, optional): Relative path to write `cstdn.py` file.
            Defaults to "./".
        vasp_cmd (tuple, optional): VASP software command to run calculation.
            Defaults to ("srun", "vasp_std").
        max_errors (int, optional): Maximum number of allowed errors handled by Custodian.
            Defaults to 3.

    Returns:
        None: `cstdn.py` file written to `output_path`.
    """

    cstdn_script = f"""# Generic Custodian script for VASP calculation.

import os
from custodian.custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler
from custodian.vasp.jobs import VaspJob

subset = list(VaspErrorHandler.error_msgs.keys())
handlers = [VaspErrorHandler(errors_subset_to_catch=subset)]

step1 = VaspJob(vasp_cmd={vasp_cmd}, final=True)

jobs = [step1]

c = Custodian(handlers, jobs, max_errors={max_errors})

c.run()
"""

    with open(f"{output_dir}/cstdn.py", "w+") as f:
        f.writelines(cstdn_script)

    return None


def write_double_relax_script(
    output_dir: str = "./",
    vasp_cmd_initial: tuple = ("srun", "vasp_std"),
    vasp_cmd_final: tuple = ("srun", "vasp_std"),
    incar_changes_initial: dict = {"KSPACING": 0.50},
    incar_changes_final: dict = {"KSPACING": 0.25},
    max_errors: int = 3,
    kpoint_meshes: list[tuple, tuple] = None,
) -> None:
    """
    Writes a double-relaxation custodian script for VASP calculations.

    The `KSPACING` INCAR flag is used to facilitate high-throughput calculations.
    - `KSPACING = 0.25` is a reasonably safe default, but do your own testing.
    - There is the option to use KPOINTS instead of `KSPACING` by using the `kpoint_meshes`
        argument. An internal check is performed to ensure that both are not used in tandem,
        which should not be done.

    Using a less dense k-mesh for the initial relaxation (larger KSPACING) can be very
    beneficial to get closer to a structural minima at a drastically reduced cost.
    - https://github.com/hackingmaterials/atomate/issues/19 shows some discussion in this

    Args:
        output_dir (str, optional): Relative path to write `cstdn.py` file.
            Defaults to "./".
        vasp_cmd_initial (tuple, optional): VASP software command to run first relaxation.
            Made to be flexible to allow for `vasp_gam` executable for initial relaxations.
            Defaults to ("srun", "vasp_std").
        vasp_cmd_final (tuple, optional): VASP software command to run second relaxation.
            Defaults to ("srun", "vasp_std").
        incar_changes_initial (_type_, optional): INCAR changes desired for first relaxation.
            Defaults to {"KSPACING": 0.50}.
        incar_changes_final (_type_, optional): INCAR changes desired for second relaxation.
            Defaults to {"KSPACING": 0.25}.
        max_errors (int, optional): Maximum number of allowed errors handled by Custodian.
            Defaults to 3.
        kpoint_meshes (list[tuple, tuple], optional): list of k-point meshes for initial and final
            relaxations. Should not be used with `KSPACING` (i.e. set incar_changes_* to empty
            dictionaries: "{}"). Defaults to None.

    Raises:
        ValueError: if both `KSPACING` and KPOINTS are attempted.

    Returns:
        None: `cstdn.py` file written to `output_path`.
    """

    if kpoint_meshes:
        if "KSPACING" in incar_changes_initial or "KSPACING" in incar_changes_final:
            raise ValueError(
                "Detected use of both KSPACING and KPOINTS file in Custodian writer."
            )

        cstdn_script = f"""# Custodian double-relaxation script for VASP calculation using KPOINTS.

import os
from custodian.custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler
from custodian.vasp.jobs import VaspJob

subset = list(VaspErrorHandler.error_msgs.keys())
handlers = [VaspErrorHandler(errors_subset_to_catch=subset)]

step1 = VaspJob(
    vasp_cmd={vasp_cmd_initial}, 
    final=False,
    suffix=".1",
    settings_override=[
        {{"dict": "INCAR", "action": {{"_set": {incar_changes_initial}}}}},
        {{"dict": "KPOINTS", "action": {{"_set": {{"kpoints": [{kpoint_meshes[0]}]}}}}}},
    ]
)

step2 = VaspJob(
    vasp_cmd={vasp_cmd_final}, 
    final=True,
    settings_override=[
        {{"dict": "INCAR", "action": {{"_set": {incar_changes_final}}}}},
        {{"dict": "KPOINTS", "action": {{"_set": {{"kpoints": [{kpoint_meshes[1]}]}}}}}},
        {{"file": "CONTCAR", "action": {{"_file_copy": {{"dest": "POSCAR"}}}}}}
    ]
)

jobs = [step1, step2]

c = Custodian(handlers, jobs, max_errors={max_errors})

c.run()
"""

    else:
        cstdn_script = f"""# Custodian double-relaxation script for VASP calculation using KSPACING.

import os
from custodian.custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler
from custodian.vasp.jobs import VaspJob

subset = list(VaspErrorHandler.error_msgs.keys())
handlers = [VaspErrorHandler(errors_subset_to_catch=subset)]

step1 = VaspJob(
    vasp_cmd={vasp_cmd_initial}, 
    final=False,
    suffix=".1",
    settings_override=[
        {{"dict": "INCAR", "action": {{"_set": {incar_changes_initial}}}}},
    ]
)

step2 = VaspJob(
    vasp_cmd={vasp_cmd_final}, 
    final=True,
    settings_override=[
        {{"dict": "INCAR", "action": {{"_set": {incar_changes_final}}}}},
        {{"file": "CONTCAR", "action": {{"_file_copy": {{"dest": "POSCAR"}}}}}}
    ]
)

jobs = [step1, step2]

c = Custodian(handlers, jobs, max_errors={max_errors})

c.run()
"""

    with open(f"{output_dir}/cstdn.py", "w+") as f:
        f.writelines(cstdn_script)

    return None
