# General Tips

*Last edited on March 22, 2025 by JTSivak*

## Structure Generation & Analysis

## VASP Calculations

## Materials Project
#### GGA/GGA+U/R2SCAN Materials Project Mixing Scheme
Careful with using the GGA_GGA+U_R2SCAN mixing scheme (https://github.com/materialsproject/pymatgen/blob/v2025.1.24/src/pymatgen/entries/mixing_scheme.py#L35-L755).
- I have had some energy inconsistencies, therefore I have not used it too much, however perhaps I am not doing something exactly correctly...
- One should preferentially use either GGA/GGA+U or R2SCAN and these are the options I mention in the package.

#### Using API keys and ".env" file

.env files can be used to load specific environments, and are especially useful for storing API keys without sharing them to .git or with anyone. These are private and should never be shared with anyone - so be careful adding them to your scripts and sharing them/adding them to a version control system (VCS) like GitHub.

A personal Materials Project API key can be generated at https://next-gen.materialsproject.org/api

The following steps can be taken to use your MP API key in `pytheos` ->

1. Get your Materials Project API key from the above link.
2. Copy  [../resources/example.env](../resources/example.env) --> [../.env](../.env) (i.e. `pytheos` root directory).
3. Add your own API key to [../.env](../.env) file by replacing `your_mp_api_key`.
4. Be sure that file is not in VCS/shared with others!

Whenever you would like to access your MP API key quickly, you can load it before querying the MP API using:
```
from pytheos.materials_project import query

mp_api_key = query.load_api_key()
```

## Machine-Learning Interatomic Potentials (MLIPs)

## FEFF for Simulated X-ray Absorption Spectra