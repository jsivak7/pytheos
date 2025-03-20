# write a general class here for slurm submission


def write_slurm_vasp(
    jobname,
    output_path,
    num_nodes=1,
    num_cpu=64,
    mem_per_cpu="3500MB",
    runtime=96,
    partition="sla-prio",
    account="ixd4_n_bc",
):
    """
    Slurm submission script writer for Roar Collab VASP jobs

    Args:
        jobname (str):
        output_path (str): relative path to write to
        num_nodes (int, optional): Defaults to 1.
        num_cpu (int, optional): Defaults to 64.
        mem_per_cpu (str, optional): Defaults to "3500MB".
        runtime (int, optional): Defaults to 96.
        partition (str, optional): Defaults to "sla-prio".
        account (str, optional): Defaults to "sbs5563_bc".
    """

    submitvasp = f"""#!/bin/bash
#SBATCH --nodes={num_nodes}
#SBATCH --ntasks-per-node={num_cpu}
#SBATCH --mem-per-cpu={mem_per_cpu}
#SBATCH --time={runtime}:00:00
#SBATCH --partition={partition}
#SBATCH --qos=burst4x
#SBATCH --account={account}
#SBATCH --output=vasp.out
#SBATCH --error=vasp.err
#SBATCH --job-name={jobname}

export UCX_TLS=all

cd $SLURM_SUBMIT_DIR
module purge
module use /storage/icds/RISE/sw8/modules/
module load vasp/vasp-6.4.1v

eval "$(conda shell.bash hook)"
conda activate chgnet-heo-screening 

python cstdn.py"""
