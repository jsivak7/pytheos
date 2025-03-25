# Tips

*Last edited on March 24, 2025*

## Structure Generation & Analysis


## VASP Calculations
### Meta-GGA band structure calculations
Band structures can be calculated for r2scan meta-gga in the following manner (assuming you have the chgcar from the final r2scan calculation):
1. pbe static calculation from r2scan chgcar with `icharg = 11`, writing the wavecar
2. r2scan static calculation with band structure high-symmetry k-point path from pbe wavecar
    - chgcar should *not* be present
    - cannot use `icharg = 11` with meta-ggas

A pbe wavecar is needed for meta-gga bandstructure calculations otherwise there are very spiky bandstructures (which I also found documented on a google search). This methodology is recommended for meta-gga and hybrid functionals on the VASP Wiki, however I added in the ICHARG = 11 portion. This allows the magnetic state that is found from r2scan to still be realized in the pbe wavecar if the electronic structure is qualitatively different between the two functionals which is what we usually have with our heo calculations.

This all is handled internally using the `pytheos.vasp.modifier.CalcModifier` class to make band structure calculation inputs.

## Materials Project

### GGA/GGA+U/R2SCAN Materials Project Mixing Scheme
Careful with using the GGA_GGA+U_R2SCAN mixing scheme (https://github.com/materialsproject/pymatgen/blob/v2025.1.24/src/pymatgen/entries/mixing_scheme.py#L35-L755).
- I have had some energy inconsistencies, therefore I have not used it too much, however perhaps I am not doing something exactly correctly...
- One should preferentially use either GGA/GGA+U or R2SCAN and these are the options I mention in the package.

### Using API keys and ".env" file

.env files can be used to load specific environments, and are especially useful for storing API keys without sharing them to .git or with anyone. These are private and should never be shared with anyone - so be careful adding them to your scripts and sharing them/adding them to a version control system (VCS) like GitHub.

A personal Materials Project API key can be generated at https://next-gen.materialsproject.org/api

The following steps can be taken to use your MP API key in `pytheos` ->

1. Get your Materials Project API key from the above link.
2. Copy  [../examples/example_api_key.env](../examples/example_api_key.env) --> [../.env](../.env) (i.e. `pytheos` root directory).
3. Add your own API key to [../.env](../.env) file by replacing `your_api_key`.
4. Be sure that file is not in VCS/shared with others!

Whenever you would like to access your MP API key quickly, you can load it before querying the MP API using:
```
from pytheos.materials_project import query

mp_api_key = query.load_api_key()
```

## Stability


## Machine-Learning Interatomic Potentials (MLIPs)


## FEFF for Simulated X-ray Absorption Spectra