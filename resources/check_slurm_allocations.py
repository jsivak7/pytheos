#!/storage/work/jts6114/software/miniconda3/bin/python
# generic script that can be customized for evaluating use on SLURM allocations
# replace __WORD__ with specific values for your allocations
# lines below can be copied for multiple allocations if desired

import os

print(
    "\033[31;1m__ALLOC_NAME__ (CPU = __TOTAL_NUM_CPU__, MEM = __TOTAL_MEMORY__) ->\033[0m"
)
os.system(
    "squeue --account=__ALLOC_NAME__ -o '%.8i %.8u %.16j %.8T %.11M %.11l %.6D %.6C %.11m  %R'"
)
print()
