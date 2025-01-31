from pytheo.utils import read, write
from pytheo.structure.generation import make_supercell, make_random

unitcell = read("MgO.poscar")

random = make_random(
    unitcell,
    (10, 5, 5),
    chemical_symbols=[
        ["Mg", "Ni", "Co", "Zn"],
        ["Mg", "Ni", "Co", "Zn"],
        ["Mg", "Ni", "Co", "Zn"],
        ["Mg", "Ni", "Co", "Zn"],
        ["O"],
        ["O"],
        ["O"],
        ["O"],
    ],
    concentrations={"Mg": 0.25, "Ni": 0.25, "Co": 0.25, "Zn": 0.25},
)

write(random, "random.vasp")
