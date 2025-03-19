# automated custodian writer functions for vasp


def write_static(
    output_dir: str = ".",
    vasp_cmd: tuple = ("srun", "vasp_std"),
) -> None:
    """
    Write a generic static Custodian script for calculation workflow and error handling.
    - https://github.com/materialsproject/custodian

    Args:
        output_dir (str): relative path to write submission file.

    Returns:
        None: {output_dir}/cstdn.py script written.
    """

    print(f"Writing static custodian script")

    cstdn_script = f"# Custodian static script.\n\nvasp_cmd = {vasp_cmd}\n"

    cstdn_script += """import os
from custodian.custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler
from custodian.vasp.jobs import VaspJob

subset = list(VaspErrorHandler.error_msgs.keys())
handlers = [VaspErrorHandler(errors_subset_to_catch=subset)]

step1 = VaspJob(
    vasp_cmd=vasp_cmd,
    final=True,
)

jobs = [step1]
c = Custodian(handlers, jobs, max_errors=3)
c.run()"""
    with open(f"{output_dir}/cstdn.py", "w+") as f:
        f.writelines(cstdn_script)


def write_doublerelax(
    output_dir: str = ".",
    vasp_cmd: tuple = ("srun", "vasp_std"),
    kspacing=0.25,
    half_kmesh_first_relax=True,
) -> None:
    """
    Write a generic double relaxation Custodian script for calculation workflow and error handling.
    - https://github.com/materialsproject/custodian

    Restricted to using the "KSPACING" INCAR flag as commonly used for high-throughput calculations
    - KSPACING = 0.25 used as a reasonably safe default -> do you own convergence testing!!

    Using half_kmesh_first_relax = True can be very beneficial to get closer to a structural minima at a drastically reduced cost.
    - see discussion -> https://github.com/hackingmaterials/atomate/issues/19

    Args:
        output_dir (str): relative path to write submission file
        kspacing (float, optional): K-point mesh spacing with VASP KSPACING tag (https://www.vasp.at/wiki/index.php/KSPACING). Defaults to 0.25.
        half_kmesh_first_relax (bool, optional): Use more sparse k-mesh for initial relax. Defaults to True.

    Returns:
        None: {output_dir}/cstdn.py script written.
    """

    print(f"Writing double-relax custodian script")

    cstdn_script = f"# Custodian double-relaxation script.\n\nvasp_cmd = {vasp_cmd}\n"

    if half_kmesh_first_relax == True:
        cstdn_script += (
            f"""kspacing_initial = {kspacing*2}\nkspacing = {kspacing}\n\n"""
        )
    else:
        cstdn_script += f"""kspacing_initial = {kspacing}\nkspacing = {kspacing}\n\n"""

    cstdn_script += """import os
from custodian.custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler
from custodian.vasp.jobs import VaspJob

subset = list(VaspErrorHandler.error_msgs.keys())
handlers = [VaspErrorHandler(errors_subset_to_catch=subset)]

step1 = VaspJob(
    vasp_cmd=vasp_cmd,
    final=False,
    suffix=".1",
    settings_override=[
        {
            "dict": "INCAR",
            "action": {
                "_set": {
                    "KSPACING": kspacing_initial
                }
            },
        },
    ],
)

step2 = VaspJob(
    vasp_cmd=vasp_cmd,
    final=True,
    settings_override=[
        {
            "dict": "INCAR",
            "action": {
                "_set": {
                    "KSPACING": kspacing
                }
            },
        },
        {"file": "CONTCAR", "action": {"_file_copy": {"dest": "POSCAR"}}},
    ],
)


jobs = [step1, step2]
c = Custodian(handlers, jobs, max_errors=3)
c.run()"""
    with open(f"{output_dir}/cstdn.py", "w+") as f:
        f.writelines(cstdn_script)
