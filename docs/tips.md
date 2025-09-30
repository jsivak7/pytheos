# Tips

## Structure Generation & Analysis
#TODO

## VASP Calculations
### Meta-GGA band structure calculations
Band structures can be calculated for r2scan meta-gga in the following manner (assuming you have the chgcar from the final r2scan calculation):
1. pbe static calculation from r2scan chgcar with `icharg = 11`, writing the wavecar
2. r2scan static calculation with band structure high-symmetry k-point path from pbe wavecar
    - chgcar should *not* be present
    - cannot use `icharg = 11` with meta-ggas

A pbe wavecar is needed for meta-gga bandstructure calculations otherwise there are very spiky bandstructures (which I also found documented on a google search). This methodology is recommended for meta-gga and hybrid functionals on the VASP Wiki, however I added in the ICHARG = 11 portion. This allows the magnetic state that is found from r2scan to still be realized in the pbe wavecar if the electronic structure is qualitatively different between the two functionals which is what we usually have with our heo calculations.

This all is handled internally using the `pytheos.vasp.modifiers.CalcModifier` class to make band structure calculation inputs.

## Materials Project

### GGA/GGA+U/R2SCAN Materials Project Mixing Scheme
Careful with using the GGA_GGA+U_R2SCAN mixing scheme (https://github.com/materialsproject/pymatgen/blob/v2025.1.24/src/pymatgen/entries/mixing_scheme.py#L35-L755).
- I have had some energy inconsistencies, therefore I have not used it too much, however perhaps I am not doing something exactly correctly...
- One should preferentially use either GGA/GGA+U or R2SCAN and these are the options I mention in the package.

### Easily using MP API keys
A personal Materials Project API key can be generated at https://next-gen.materialsproject.org/api

One can easily use their own personal API key by add the following to their .pmgrc.yaml file:
```
PMG_MAPI_KEY: your_mp_api_key
```

MPRester is able to automatically read this when used, thus you do not have to enter an MP API key when querying from the Materials Project Database.

## Stability
#TODO

## Machine-Learning Interatomic Potentials (MLIPs)
#TODO