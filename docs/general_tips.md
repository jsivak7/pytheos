# General Tips

Careful with using the GGA_GGA+U_R2SCAN mixing scheme (https://github.com/materialsproject/pymatgen/blob/v2025.1.24/src/pymatgen/entries/mixing_scheme.py#L35-L755).
- I have had some energy inconsistencies, therefore I have not used it too much, however perhaps I am not doing something exactly correctly...
- One should preferentially use either GGA/GGA+U or R2SCAN and these are the options I mention in the package.