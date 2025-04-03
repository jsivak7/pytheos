# slurm submission script writer function for vasp calculations
# specific for Penn State ROAR COLLAB HPC
# defaults to using standard cores
# no longer handles open allocation submissions - rarely used now since switched slower hardware


def write_vasp_submission(
    job_name: str = "VASP_CALC",
    output_dir: str = "./",
    num_nodes: int = 1,
    num_cpu: int = 48,
    mem_per_cpu: str = "5000MB",
    runtime: int = 72,
    allocation: str = "ixd4_s_sc",
) -> None:
    """
    Writes VASP submission script --> `submitvasp`.

    Assumes calculations are being run via a `cstdn.py` script.

    Args:
        job_name (str, optional): Name for vasp calculation job. Defaults to "VASP_CALC.
        output_dir (str, optional): Relative directory to write submission script. Defaults to "./".
        num_nodes (int, optional): Number of nodes. Defaults to 1.
        num_cpu (int, optional): Number of CPUs per node. Defaults to 48.
        mem_per_cpu (str, optional): Amount of memory per CPU. Defaults to "5000MB".
        runtime (int, optional): Run time for job submission. Defaults to 72.
        allocation (str, optional): Allocation for job to run on. Defaults to "ixd4_s_sc".
    """

    slurm_script = f"""#!/bin/bash
#SBATCH --nodes={num_nodes}
#SBATCH --ntasks-per-node={num_cpu}
#SBATCH --mem-per-cpu={mem_per_cpu}
#SBATCH --time={runtime}:00:00
#SBATCH --partition=sla-prio
#SBATCH --account={allocation}
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --job-name={job_name}

export UCX_TLS=all

cd $SLURM_SUBMIT_DIR
module purge
module use /storage/icds/RISE/sw8/modules/
module load vasp/vasp-6.4.1v

eval "$(conda shell.bash hook)"
conda activate pytheos

python cstdn.py"""

    with open(f"{output_dir}/submitvasp", "w+") as f:
        f.writelines(slurm_script)

    return None
