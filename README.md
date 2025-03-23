<h1 align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="resources/logo_light.svg">
    <img alt="Logo" src="resources/logo.svg"
height="150">
  </picture>
</h1>

Python toolkit for high-entropy oxide simulations (pytheos) is an open-source package built for exploring chemically disordered crystalline materials with computational models. 

Density-functional theory calculations are the primary focus of the code currently as they provide very accurate thermodynamic and electronic structure predictions, however other types of simulations and machine-learning approaches will be implemented in the future.

A variety of other open-source packages such as [pymatgen](https://github.com/materialsproject/pymatgen), [custodian](https://github.com/materialsproject/custodian), and [sumo](https://github.com/SMTG-Bham/sumo) as well as the [Materials Project database](https://next-gen.materialsproject.org/) are heavily utilized to leverage existing computational machinery and data.

## Installation
1. Clone the repo
```
git clone https://github.com/jsivak7/pytheos
cd pytheos 
```
2. Install the package
- `-e` is for installing in editable mode (optional, but recommended).
```
pip install -e .
```

*Jupyter notebooks are provided to demonstrate how to use the code in the `examples/` folder.*

## Development & Bugs
In case you run into any issues please contact Jacob via email (jts6114@psu.edu) as this package is under active development (and is still changing very often). Thank you for any suggestions and contributions - I hope that this is useful to your work!