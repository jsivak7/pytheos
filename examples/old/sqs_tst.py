from pytheo.structure import read, write, make_supercell, make_sqs

unitcell = read("MgO.poscar")
sqs = make_sqs(
    unitcell,
    (5, 2, 2),
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
    cutoffs=[4.21],
    concentrations={"Mg": 0.25, "Ni": 0.25, "Co": 0.25, "Zn": 0.25},
    num_steps=15000,
)

write(sqs, "sqs.vasp")
