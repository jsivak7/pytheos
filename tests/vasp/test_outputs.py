from pytheos.vasp.outputs import CalcOutputs, load_vasprun
import pytest
import os
from pymatgen.io.vasp.outputs import Vasprun

current_dir = module_dir = os.path.dirname(__file__)


def test_load_vasprun():
    """
    Tests that vasprun.xml files can be properly loaded as Vasprun pymatgen objects.
    """

    vrun = load_vasprun(path=f"{current_dir}/../files/vasprun.xml")

    assert isinstance(vrun, Vasprun)


def test_CalcOutputs():
    """
    Tests that CalcOutputs class properly holds generic output data for further processing/analysis.
    """

    calcoutputs = CalcOutputs(source_dir=f"{current_dir}/../files")

    num_atoms = calcoutputs.num_atoms
    final_energy = calcoutputs.final_energy

    assert num_atoms == 8
    assert final_energy == -64.50440241
    assert calcoutputs.final_energy_per_atom == final_energy / num_atoms
    assert calcoutputs.lattice_parameters == pytest.approx((4.20, 4.20, 4.20), rel=1e-3)
    assert calcoutputs.lattice_angles == pytest.approx((90, 90, 90), rel=1e-4)
