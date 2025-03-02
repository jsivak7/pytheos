# TODO just copied in for now from the perovskite HEO project - still needs to be fixed
def set_up_optical_calc(sigma=0.1) -> None:
    """
    For making a new directory for optical calculations (should be called in directory where "02_static" exists).
    - NBANDS is doubled from the default value of 02_static calculation

    Args:
        sigma (float, optional): Smearing parameter than is usually increased for optical calculations compared to relaxation. Defaults to 0.2.
    """
    import os
    from pymatgen.io.vasp.outputs import Eigenval, Vasprun

    print(f"-> Making optical directory...")

    os.mkdir("04_optical")
    os.system(
        f"cp 02_static/INCAR 02_static/POSCAR 02_static/KPOINTS 02_static/POTCAR 02_static/CHGCAR 02_static/WAVECAR 02_static/submitvasp ./04_optical"
    )
    os.chdir("04_optical")

    e = Eigenval("../02_static/EIGENVAL")

    nbands = e.nbands
    new_nbands = nbands * 2  # doubling NBANDS
    os.system(f"echo NBANDS = {new_nbands} >> INCAR")
    print(f"number of bands: {nbands} -> {new_nbands}")

    os.system("perl -pi -e 's/LWAVE = True/LWAVE = False/g' INCAR")
    os.system("perl -pi -e 's/LCHARG = True/LCHARG = False/g' INCAR")
    os.system(f"echo NEDOS = 20000 >> INCAR")  # for adequate sampling
    os.system(f"perl -pi -e 's/SIGMA = 0.05/SIGMA = {sigma}/g' INCAR")
    os.system("perl -pi -e 's/ALGO = All/ALGO = Normal/g' INCAR")
    os.system("perl -pi -e 's/EDIFF = 1e-06/EDIFF = 1e-08/g' INCAR")
    os.system("perl -pi -e 's/LREAL = Auto/LREAL = False/g' INCAR")
    os.system(f"echo LOPTICS = True >> INCAR")
    os.system(f"echo CSHIFT = 0.01 >> INCAR")


# TODO just copied in for now from the perovskite HEO project - still needs to be fixed
def set_up_bandstructure_calc():
    """For making a new directory for density of states calculations (should be called in directory where "02_static" exists)
    - NBANDS is 1.5x from the default value of 02_static calculation

    Additional, manual steps are required for this function still due to band structure process:
        NOTE that this is for the 'regular' band structure calculations (i.e., not unfolding)
        1. run 01_pbe calculation
        2. get high-symmetry k-path with sumo (> sumo-kgen --hybrid --symprec 0.1 --pymatgen)
        3. copy KPOINTS_band to KPOINTS
        4. copy WAVECAR from finished 01_pbe calculation to the 'bandstructure' directory
        5. now can run actual metaGGA band structure calculation
    """

    import os
    from pymatgen.io.vasp.outputs import Eigenval

    print(f"-> Making bandstructure directory...")

    os.mkdir("05_bandstructure")
    os.system(
        f"cp 02_static/INCAR 02_static/IBZKPT 02_static/POSCAR 02_static/KPOINTS 02_static/POTCAR 02_static/CHGCAR 02_static/WAVECAR 02_static/submitvasp ./05_bandstructure"
    )
    os.chdir("05_bandstructure")

    e = Eigenval("../02_static/EIGENVAL")

    nbands = e.nbands
    new_nbands = int(nbands * 1.5)  # increasing NBANDS by 1.5x
    os.system(f"echo NBANDS = {new_nbands} >> INCAR")
    print(f"number of bands: {nbands} -> {new_nbands}")

    os.system("perl -pi -e 's/LCHARG = True/LCHARG = False/g' INCAR")
    os.system(f"echo NEDOS = 5001 >> INCAR")  # for adequate sampling
    os.system("perl -pi -e 's/LREAL = Auto/LREAL = False/g' INCAR")

    os.mkdir("01_pbe")
    os.system("cp * 01_pbe")
    os.system("rm CHGCAR WAVECAR")
    os.chdir("01_pbe")

    os.system("perl -pi -e 's/METAGGA = R2scan/GGA = PE/g' INCAR")
    os.system(f"echo ICHARG = 11 >> INCAR")


# TODO need to ensure this is properly implemented
def get_magorder(calc_path: str, magorder_name: str, coordinate_rounding=0):

    from ase.io import read
    import numpy as np
    import random
    import os
    import yaml

    s = read(f"{calc_path}/POSCAR")

    coords = s.get_positions()
    coords = np.around(
        coords, coordinate_rounding
    )  # optional to change coordinate rounding of inputted structure

    elems = s.get_chemical_symbols()

    module_dir = os.path.dirname(__file__)

    with open(
        f"{module_dir}/incar_defaults/r2scan.yaml", "r"
    ) as f:  # need to make this more general...
        incar_settings = yaml.load(f, Loader=yaml.SafeLoader)
    magmoms = incar_settings["MAGMOM"]

    with open(f"{module_dir}/magnetic_orders/{magorder_name}.yaml", "r") as f:
        spin_ordering = yaml.load(f, Loader=yaml.SafeLoader)

    spin_up = np.around(np.array(spin_ordering["spin_up"]), 1)
    spin_down = np.around(np.array(spin_ordering["spin_down"]), 1)

    up_index = []
    down_index = []

    for i in range(len(coords)):
        for j in range(
            len(spin_up)
        ):  # NOTE that this needs to be an equal amount of spin up and down with the current implementation
            if (
                np.array_equal(spin_up[j], coords[i]) == True
            ):  # this tests for if the atomic coordinates match what the user defines for a spin-up atom
                up_index.append(
                    i
                )  # if it matches, then append the atom number to the up_index list
            elif (
                np.array_equal(spin_down[j], coords[i]) == True
            ):  # this tests for if the atomic coordinates match what the user defines for a spin-up atom
                down_index.append(
                    i
                )  # if it matches, then append the atom number to the down_index list

    magmom = []  # empty list that will be in the same order as poscar file
    # loop for creating a list that will correspond to the final magmom tag that will be added to incar file
    for i in range(
        len(coords)
    ):  # loops through all of the atomic numbers within the inputted structure
        if i in up_index:
            magmom.append(
                +np.round(
                    random.uniform(magmoms[elems[i]] - 0.5, magmoms[elems[i]] + 0.5), 2
                )
            )
        elif i in down_index:
            magmom.append(
                -np.round(
                    random.uniform(magmoms[elems[i]] - 0.5, magmoms[elems[i]] + 0.5), 2
                )
            )
        elif elems[i] == "O":
            magmom.append(magmoms["O"])
        elif elems[i] == "Sr":
            magmom.append(magmoms["Sr"])
        elif elems[i] == "La":
            magmom.append(magmoms["La"])
        else:
            print(f"WARNING!!!! An atom was not identified. Atom #{i}")

    print(f"\nnumber spin up atoms:\t{len(up_index)}")
    print(f"number spin down atoms:\t{len(down_index)}")

    print(f"magnetic ordering:\t{magmom}")
    if len(up_index) != len(down_index):
        print("WARNING!!!! The up index does not equal the down index.")

    # write the magmom list (the actual magnetic ordering) to a temporary file "tmp_ORDER.txt"
    with open("tmp_ORDER.txt", "w") as output:
        output.write(str(magmom))

    # the following are so that the magnetic ordering is added to the MAGMOM tag in the proper format
    os.system("perl -pi -e 's/\\[//g' tmp_ORDER.txt")  # get rid of bracket
    os.system("perl -pi -e 's/\\]//g' tmp_ORDER.txt")  # get rid of bracket
    os.system("perl -pi -e 's/,//g' tmp_ORDER.txt")  # get rid of commas
    os.system("echo 'MAGMOM = CHANGE' >> MAGMOM.txt")
    os.system(f"perl -pi -e 's/MAGMOM.*/MAGMOM = CHANGE/g' {calc_path}/INCAR")
    os.system(
        f"perl -pe 's/CHANGE/`cat tmp_ORDER.txt`/e' -i {calc_path}/INCAR"
    )  # this replaces the "CHANGE" that we just added with the magnetic ordering for the corresponding POSCAR file
    os.system("rm MAGMOM.txt ")
    os.system(f"mv tmp_ORDER.txt {calc_path}/magorder")
