# Custodian script to run VASP calculations in a double-relaxation scheme using KPOINTS
# https://github.com/materialsproject/custodian


vasp_cmds = [("srun", "vasp_std"), ("srun", "vasp_std")]  # [1st calc, 2nd calc]
incar_changes = [{"KSPACING": 0.44}, {"KSPACING": 0.22}]  # [1st calc, 2nd calc]
kpoint_meshes = [(2, 2, 2), (4, 4, 4)]  # [1st calc, 2nd calc]
max_errors = 3


from custodian.custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler
from custodian.vasp.jobs import VaspJob

subset = list(VaspErrorHandler.error_msgs.keys())
handlers = [VaspErrorHandler(errors_subset_to_catch=subset)]

step1 = VaspJob(
    vasp_cmd={vasp_cmds[0]},
    final=False,
    suffix=".1",
    settings_override=[
        {{"dict": "INCAR", "action": {{"_set": {incar_changes[0]}}}}},
        {
            {
                "dict": "KPOINTS",
                "action": {{"_set": {{"kpoints": [{kpoint_meshes[0]}]}}}},
            }
        },
    ],
)

step2 = VaspJob(
    vasp_cmd={vasp_cmds[1]},
    final=True,
    settings_override=[
        {{"dict": "INCAR", "action": {{"_set": {incar_changes[1]}}}}},
        {
            {
                "dict": "KPOINTS",
                "action": {{"_set": {{"kpoints": [{kpoint_meshes[1]}]}}}},
            }
        },
        {{"file": "CONTCAR", "action": {{"_file_copy": {{"dest": "POSCAR"}}}}}},
    ],
)

jobs = [step1, step2]

c = Custodian(handlers, jobs, max_errors={max_errors})

c.run()
