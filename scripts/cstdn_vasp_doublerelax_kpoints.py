# Custodian script to run VASP calculations in a double-relaxation scheme using KPOINTS
# https://github.com/materialsproject/custodian


vasp_cmd_run1 = ("srun", "vasp_std")
vasp_cmd_run2 = ("srun", "vasp_std")

incar_changes_run1 = {}
incar_changes_run2 = {}

kpoint_mesh_run1 = [2, 2, 2]
kpoint_mesh_run2 = [4, 4, 4]

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
        {
            "dict": "KPOINTS",
            "action": {
                "_set": {
                    "kpoints": [kpoint_mesh_run1],
                }
            },
        },
    ],
)

step2 = VaspJob(
    vasp_cmd=vasp_cmd_run1,
    final=True,
    settings_override=[
        {
            "dict": "INCAR",
            "action": {
                "_set": incar_changes_run2,
            },
        },
        {
            "dict": "KPOINTS",
            "action": {
                "_set": {
                    "kpoints": [kpoint_mesh_run2],
                }
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
