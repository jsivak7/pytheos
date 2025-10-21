from pytheos.vasp.inputs import CalcInputs, write_submission_script
import pytest
from pymatgen.core import SETTINGS, Structure
from pymatgen.io.vasp.inputs import Incar, Poscar, Potcar, Kpoints
from ase.build import bulk
import os
from ase.io import read

current_dir = module_dir = os.path.dirname(__file__)
unitcell = bulk(name="MgO", crystalstructure="rocksalt", a=4.2)
unitcell = Structure.from_ase_atoms(unitcell)


@pytest.fixture
def fake_potcar_patch(monkeypatch):
    """
    Creates a path to the fake POTCAR files with scrambled data for copyright purposes
    """

    monkeypatch.setitem(
        dic=SETTINGS,
        name="PMG_VASP_PSP_DIR",
        value=f"{current_dir}/../files/fake_potcars",
    )


def test_calc_inputs_init(fake_potcar_patch):
    """
    Tests that CalcInputs class has proper attributes during its default initialization.
    """

    calcinputs = CalcInputs(struc=unitcell)

    assert isinstance(calcinputs.struc, Structure)
    assert isinstance(calcinputs.incar, Incar)
    assert isinstance(calcinputs.poscar, Poscar)
    assert isinstance(calcinputs.potcar, Potcar)
    assert not hasattr(calcinputs, "kpoints")


def test_update_incar(fake_potcar_patch):

    calcinputs = CalcInputs(struc=unitcell)

    # just changing random things for testing purposes
    calcinputs.update_incar(
        {
            "ALGO": 10000,
            "KSPACING": 5,
            "ENCUT": 50000,
            "KPAR": 64,
        }
    )

    assert calcinputs.incar["ALGO"] == 10000
    assert calcinputs.incar["KSPACING"] == 5
    assert calcinputs.incar["ENCUT"] == 50000
    assert calcinputs.incar["KPAR"] == 64


def test_use_kpoints_file(fake_potcar_patch):
    """
    Tests that KPOINTS file can be used instead of KSPACING INCAR tag.

    Ensures that both KPOINTS file and KSPACING tag are not both used,
    as this could lead to undesirable results.
    """

    calcinputs = CalcInputs(struc=unitcell)

    calcinputs.use_kpoints_file(kpoint_mesh=[5, 5, 5])

    assert "KSPACING" not in calcinputs.incar
    assert hasattr(calcinputs, "kpoints")
    assert calcinputs.kpoints == Kpoints(kpts=[5, 5, 5])


def test_apply_mag_order(fake_potcar_patch):
    """
    Ensures that the automatic magnetic ordering generator method works properly
    for a given magorder.yaml file. Probably could be more complete, however just
    wanted to make sure this is working as expected.
    """

    supercell = Structure.from_file(f"{current_dir}/../files/NiO_prim_2x2x2.poscar")

    calcinputs = CalcInputs(struc=supercell)

    calcinputs.apply_mag_order(
        magmom_values={
            "Ni": 3,
            "O": 0,
        },
        mag_order_file=f"{current_dir}/../files/magorder_rocksalt_2x2x2.yaml",
        rattle_amount=0,  # to ensure we get expected MAGMOM in INCAR
    )

    assert calcinputs.incar["MAGMOM"] == [
        3.0,
        -3.0,
        -3.0,
        3.0,
        -3.0,
        3.0,
        3.0,
        -3.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
    ]


def test_write_files(fake_potcar_patch, tmp_path):

    calcinputs = CalcInputs(struc=unitcell)

    test_dir_name = "TEST_VASP_CALC"

    calcinputs.write_files(output_dir=tmp_path / test_dir_name)

    expected_dir_path = tmp_path / test_dir_name
    expected_incar_path = expected_dir_path / "INCAR"
    expected_poscar_path = expected_dir_path / "POSCAR"
    expected_potcar_path = expected_dir_path / "POTCAR"

    test_incar_content = """ALGO = All
EDIFF = 1e-06
EDIFFG = -0.02
ENAUG = 1360
ENCUT = 680.0
IBRION = 2
ISIF = 3
ISMEAR = 0
ISPIN = 2
ISYM = 0
KPAR = 1
KSPACING = 0.25
LAECHG = False
LASPH = True
LCHARG = True
LELF = False
LMIXTAU = True
LORBIT = 11
LREAL = Auto
LVTOT = False
LWAVE = True
MAGMOM = 2*0.0
METAGGA = R2scan
NCORE = 32
NELM = 500
NSW = 250
POTIM = 0.25
PREC = Accurate
SIGMA = 0.05
TIME = 0.1
"""

    test_poscar_content = """Mg1 O1
1.0
   0.0000000000000000    2.1000000000000001    2.1000000000000001
   2.1000000000000001    0.0000000000000000    2.1000000000000001
   2.1000000000000001    2.1000000000000001    0.0000000000000000
Mg O
1 1
direct
   0.0000000000000000    0.0000000000000000    0.0000000000000000 Mg
  -0.5000000000000000    0.5000000000000000    0.5000000000000000 O
"""

    assert expected_dir_path.is_dir()
    assert expected_incar_path.is_file()
    assert expected_poscar_path.is_file()
    assert expected_potcar_path.is_file()
    assert expected_poscar_path.read_text() == test_poscar_content
    assert expected_incar_path.read_text() == test_incar_content


def test_get_mprelaxset_inputs(fake_potcar_patch):

    calcinputs = CalcInputs(struc=unitcell, mp_input_set="MPRelaxSet")

    assert calcinputs.incar["ENCUT"] == 520
    assert calcinputs.incar["GGA"] == "Pe"
    assert calcinputs.incar["EDIFF"] == 1e-05


def test_get_mpscanrelaxset_inputs(fake_potcar_patch):

    calcinputs = CalcInputs(struc=unitcell, mp_input_set="MPScanRelaxSet")

    assert calcinputs.incar["ENCUT"] == 680
    assert calcinputs.incar["METAGGA"] == "R2scan"
    assert calcinputs.incar["EDIFF"] == 1e-06


def test_write_submission_script(tmp_path):

    expected_script_path = tmp_path / "submitvasp"

    test_submission_script_content = """#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --mem-per-cpu=3500MB
#SBATCH --time=72:00:00
#SBATCH --partition=sla-prio
#SBATCH --account=sbs5563_bc
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --job-name=test-vasp-calc

export UCX_TLS=all

cd $SLURM_SUBMIT_DIR
module purge
module use /storage/icds/RISE/sw8/modules/
module load vasp/vasp-6.4.1v

eval "$(conda shell.bash hook)"
conda activate pytheos

python cstdn.py"""

    write_submission_script(job_name="test-vasp-calc", output_dir=tmp_path)

    assert expected_script_path.read_text() == test_submission_script_content
