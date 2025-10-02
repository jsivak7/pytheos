# generic Custodian script to run VASP calculations
# https://github.com/materialsproject/custodian


vasp_cmd = ("srun", "vasp_std")
max_errors = 3


from custodian.custodian import Custodian
from custodian.vasp.handlers import VaspErrorHandler
from custodian.vasp.jobs import VaspJob

subset = list(VaspErrorHandler.error_msgs.keys())
handlers = [VaspErrorHandler(errors_subset_to_catch=subset)]

step1 = VaspJob(vasp_cmd=vasp_cmd, final=True)

jobs = [step1]

c = Custodian(handlers, jobs, max_errors=max_errors)

c.run()
