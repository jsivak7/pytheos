# Custodian script to run VASP calculations in a double-relaxation scheme using KSPACING
# https://github.com/materialsproject/custodian


vasp_cmd_run1 = ("srun", "vasp_gam")  # modified for Gamma point only calculation
vasp_cmd_run2 = ("srun", "vasp_std")

incar_changes_run1 = {"KSPACING": 5}  # for Gamma point only
incar_changes_run2 = {"KSPACING": 0.25}

max_errors = 3


from custodian.custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler
from custodian.vasp.jobs import VaspJob

subset = list(VaspErrorHandler.error_msgs.keys())
handlers = [VaspErrorHandler(errors_subset_to_catch=subset)]

step1 = VaspJob(
    vasp_cmd=vasp_cmd_run1,
    final=False,
    suffix=".1",
    settings_override=[
        {
            "dict": "INCAR",
            "action": {
                "_set": incar_changes_run1,
            },
        },
    ],
)

step2 = VaspJob(
    vasp_cmd=vasp_cmd_run2,
    final=True,
    settings_override=[
        {
            "dict": "INCAR",
            "action": {
                "_set": incar_changes_run2,
            },
        },
        {
            "file": "CONTCAR",
            "action": {
                "_file_copy": {
                    "dest": "POSCAR",
                }
            },
        },
    ],
)

jobs = [step1, step2]

c = Custodian(handlers, jobs, max_errors=max_errors)

c.run()
