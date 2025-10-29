<h1 align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/pytheos_logo_darkmode.svg">
    <img alt="Logo" src="docs/pytheos_logo.svg"
height="150">
  </picture>
</h1>

<h4 align="center">

[![Tests](https://github.com/jsivak7/pytheos/actions/workflows/tests.yml/badge.svg)](https://github.com/jsivak7/pytheos/actions/workflows/tests.yml)

</h4>

Python toolkit for high-entropy oxides (`pytheos`) is an open-source package built for exploring chemically disordered crystalline materials using computational methods. 

Our goal for this package is to consolidate code that has assisted our team in our understanding high-entropy oxide materials. We primarily use density functional theory (DFT) calculations for their ability to provide accurate predictions at the electronic and atomic scale. However, we have also started to implement machine-learning interatomic potentials (MLIPs) into our workflows as they have now become accurate enough to rapidly search over wide compositional spaces with drastically reduced costs compared to DFT. Other open-source packages such as [pymatgen](https://github.com/materialsproject/pymatgen), [custodian](https://github.com/materialsproject/custodian), and [sumo](https://github.com/SMTG-Bham/sumo) as well as the [Materials Project database](https://next-gen.materialsproject.org/) are utilized to leverage existing computational machinery and databases.

A list of publications that have led to the development of `pytheos` can be found [here](docs/publications.md).

Notebooks showing how to use `pytheos` can be found [here](notebooks/).

## Installation
1. (*optional*) Create a new python virtual environment and activate it.
- Requires `conda` (recommend using [Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main) as it is fairly lightweight).
```
conda create --name pytheos python=3.12
conda activate pytheos
```

2. Clone the repository.
```
git clone https://github.com/jsivak7/pytheos
cd pytheos 
```

3. Install the package.
- `-e` is for installing in editable mode (optional, but recommended).
```
pip install -e .
```

## Development & Bugs
In case you run into any issues please contact Jacob (jts6114@psu.edu) as this package is under active development. Thank you for any suggestions and contributions!

## Funding
This package is supported by the Penn State MRSEC [[DMR-2011839](https://www.mrsec.psu.edu)].

<h1 align="center">
  <picture>
    <img alt="Logo" src="docs/mrsec_logo.png"
height="150">
  </picture>
</h1>


