from pytheos.vasp.inputs import CalcInputs, write_submission_script
import pytest
from pymatgen.core import SETTINGS, Structure
from pymatgen.io.vasp.inputs import Incar, Poscar, Potcar, Kpoints, VaspInput
from ase.build import bulk
import os

current_dir = module_dir = os.path.dirname(__file__)
unitcell = bulk(name="MgO", crystalstructure="rocksalt", a=4.2)


@pytest.fixture
def fake_potcar_patch(monkeypatch):
    monkeypatch.setitem(
        dic=SETTINGS,
        name="PMG_VASP_PSP_DIR",
        value=f"{current_dir}/../files/fake_potcars",
    )


def test_calc_inputs_init(fake_potcar_patch):

    calcinputs = CalcInputs(structure=unitcell)

    assert isinstance(calcinputs.structure, Structure)
    assert isinstance(calcinputs.incar, Incar)
    assert isinstance(calcinputs.poscar, Poscar)
    assert isinstance(calcinputs.potcar, Potcar)
    assert not hasattr(calcinputs, "kpoints")


def test_update_incar(fake_potcar_patch):

    calcinputs = CalcInputs(structure=unitcell)

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

    calcinputs = CalcInputs(structure=unitcell)

    calcinputs.use_kpoints_file(kpoint_mesh=[5, 5, 5])

    assert "KSPACING" not in calcinputs.incar
    assert hasattr(calcinputs, "kpoints")
    assert calcinputs.kpoints == Kpoints(kpts=[5, 5, 5])


def test_get_mprelaxset_inputs(fake_potcar_patch):

    calcinputs = CalcInputs(structure=unitcell, mp_input_set="MPRelaxSet")

    assert calcinputs.incar["ENCUT"] == 520
    assert calcinputs.incar["GGA"] == "Pe"
    assert calcinputs.incar["EDIFF"] == 1e-05


def test_get_mpscanrelaxset_inputs(fake_potcar_patch):

    calcinputs = CalcInputs(structure=unitcell, mp_input_set="MPScanRelaxSet")

    assert calcinputs.incar["ENCUT"] == 680
    assert calcinputs.incar["METAGGA"] == "R2scan"
    assert calcinputs.incar["EDIFF"] == 1e-06
