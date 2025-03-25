#!/storage/work/jts6114/software/miniconda3/bin/python

# generic script that can be customized for evaluating use on SLURM allocations
# last edited on 2025-03-25

import os

# ixd4_s_sc
print("\033[31;1mixd4_s_sc (CPU = 400, MEM = 4000GB) ->\033[0m")
os.system(
    "squeue --account=ixd4_s_sc -o '%.8i %.8u %.16j %.8T %.11M %.11l %.6D %.6C %.11m  %R'"
)
print()

# sbs5563_p_sc
print("\033[31;1msbs5563_p_sc (CPU = 400, MEM = 4000GB) ->\033[0m")
os.system(
    "squeue --account=sbs5563_p_sc -o '%.8i %.8u %.16j %.8T %.11M %.11l %.6D %.6C %.11m  %R'"
)
print()

# ixd4_n_bc
print("\033[31;1mixd4_n_bc (CPU = 600, MEM = 2400GB) ->\033[0m")
os.system(
    "squeue --account=ixd4_n_bc -o '%.8i %.8u %.16j %.8T %.11M %.11l %.6D %.6C %.11m  %R'"
)
print()

# sbs5563_bc
print("\033[31;1msbs5563_bc (CPU = 720, MEM = 2880GB) ->\033[0m")
os.system(
    "squeue --account=sbs5563_bc -o '%.8i %.8u %.16j %.8T %.11M %.11l %.6D %.6C %.11m  %R'"
)
print()

# sbs5563
print("\033[31;1msbs5563 (CPU = 400, MEM = 1600GB) ->\033[0m")
os.system(
    "squeue --account=sbs5563 -o '%.8i %.8u %.16j %.8T %.11M %.11l %.6D %.6C %.11m %R'"
)
print()

# sbs5563_sc
print("\033[31;1msbs5563_sc (CPU = 100, MEM = 1000GB) ->\033[0m")
os.system(
    "squeue --account=sbs5563_sc -o '%.8i %.8u %.16j %.8T %.11M %.11l %.6D %.6C %.11m  %R'"
)
print()
