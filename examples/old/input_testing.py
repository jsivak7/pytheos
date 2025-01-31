from pytheo import utils
from pytheo.dft import inputs
from pytheo.structure import generation

s = utils.read("MgO.poscar")
supercell = generation.make_random(
    s,
    dimensions=(2, 2, 1),
    chemical_symbols=[
        ["Mg", "Co", "Ni", "Zn"],
        ["Mg", "Co", "Ni", "Zn"],
        ["Mg", "Co", "Ni", "Zn"],
        ["Mg", "Co", "Ni", "Zn"],
        ["O"],
        ["O"],
        ["O"],
        ["O"],
    ],
    concentrations={
        "Mg": 1 / 4,
        "Co": 1 / 4,
        "Ni": 1 / 4,
        "Zn": 1 / 4,
    },
)
inputs.make_relax(supercell, "random")
poscar = utils.read("random/POSCAR")
print(poscar)
utils.write(poscar, "random/POSCAR_pristine")
poscar_rattled = generation.rattle_atoms(poscar)
print(poscar_rattled)
utils.write(poscar_rattled, "random/POSCAR", overwrite=True)
inputs.write_psu_roar_collab_submission("random", "random", type="open")
inputs.write_custodian_relax("random")
