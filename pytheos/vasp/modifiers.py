# For modifying VASP input files for different calculation types (DOS, BS, BADER, etc.)


def make_dos_calc() -> None:
    """
    Moves all relaxation files to a new directory called '01_relax' and gets everything ready for a density of states (DOS) calculation.

    Call in the same location as relaxation files.

    Changes the NELECT value per https://www.vasp.at/forum/viewtopic.php?t=17981 -> NELECT - 1 + 0.999999
    - ensures proper Fermi level placement
    """

    import os
    from pymatgen.io.vasp.outputs import Eigenval, Vasprun

    os.mkdir("01_relax")
    os.system("mv * 01_relax")

    os.system(
        "cp 01_relax/CONTCAR 01_relax/INCAR 01_relax/WAVECAR 01_relax/CHGCAR 01_relax/POTCAR 01_relax/submitvasp ."
    )
    os.system("mv CONTCAR POSCAR")

    # since not using custodian here usually
    os.system("perl -pi -e 's/python cstdn.py/srun vasp_std/g' submitvasp")

    # due to issue with fermi level placement (see https://www.vasp.at/forum/viewtopic.php?t=17981)
    e = Eigenval("01_relax/EIGENVAL")
    num_elec = e.nelect
    new_num_elec = num_elec - 1 + 0.999999
    os.system(f"echo NELECT = {new_num_elec} >> INCAR")
    print(f"number of electrons: {num_elec} -> {new_num_elec}")

    # for a more reasonable energy window to sample over
    v = Vasprun("01_relax/vasprun.xml")
    efermi = v.efermi
    print(f"fermi energy = {efermi} eV")
    emin = efermi - 8
    emax = efermi + 6
    os.system(f"echo EMIN = {emin} >> INCAR")
    os.system(f"echo EMAX = {emax} >> INCAR")

    # change other VASP flags for DOS calculation
    os.system("perl -pi -e 's/ISMEAR = 0/ISMEAR = -5/g' INCAR")
    os.system("perl -pi -e 's/NSW = 250/NSW = 0/g' INCAR")
    os.system("perl -pi -e 's/LCHARG = True/LCHARG = False/g' INCAR")
    os.system("perl -pi -e 's/LWAVE = True/LWAVE = False/g' INCAR")
    os.system(
        "perl -pi -e 's/ALGO.*/ALGO = Normal/g' INCAR"
    )  # usually better for ISMEAR = -5

    os.system(f"echo NEDOS = 501 >> INCAR")


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
