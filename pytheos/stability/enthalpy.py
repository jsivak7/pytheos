# for evaluting enthalpic stability

from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry
from pymatgen.io.vasp import Vasprun
from pymatgen.entries.compatibility import MaterialsProject2020Compatibility
import numpy as np


class EnthalpicStability:
    """
    Class for evaluting the enthalpic stability of a material.

    Useful for calculating relative stability of your target material with respect to the
    other calculations in phase diagram.

    Attributes:
        phase_diagram (PhaseDiagram): Pymatgen PhaseDiagram class with added target entry.
        target_entry_name (str): Name of target entry. Defaults to "my_PDEntry".
        form_enthalpy (float): Calculated formation enthalpy of target material in eV/atom.
        decomp_enthalpy (float): Calculated decomposition enthalpy of target material in eV/atom.
        decomp_rxn (str): Decomposition reaction corresponding to the calculated decomposition enthalpy.
    """

    def __init__(
        self,
        phase_diagram: PhaseDiagram,
        target_entry_name: str = "my_PDEntry",
    ) -> None:
        """
        Args:
            phase_diagram (PhaseDiagram): phase diagram with target entry added.
            target_entry_name (str, optional): Name of target entry. Defaults to "my_PDEntry".
        """

        self.phase_diagram = phase_diagram
        self.target_entry_name = target_entry_name
        self.form_enthalpy = self._get_formation_enthalpy()
        self.decomp_enthalpy = self._get_decomposition_enthalpy()
        self.decomp_rxn = self._get_decomposition_rxn()

    def _find_target_entry(self) -> PDEntry:
        """
        Finds the target entry in supplied PhaseDiagram for stability evaluation.

        Returns:
            PDEntry: target PDEntry.
        """

        for entry in self.phase_diagram.all_entries:

            if entry.name == self.target_entry_name:

                return entry

    def _get_formation_enthalpy(self):
        """
        Calculates the formation entahlpy of target entry using supplied PhaseDiagram.

        The formation enthalpy is the enthalpy relative to the elemental references.

        Returns:
            float: Calculated formation enthalpy in eV/atom.
        """

        target_entry = self._find_target_entry()
        self.form_enthalpy = self.phase_diagram.get_form_energy_per_atom(
            entry=target_entry
        )

        return self.form_enthalpy

    def _get_decomposition_enthalpy(self):
        """
        Calculates the decomposition enthalpy of target entry using supplied PhaseDiagram.

        The decomposition enthalpy is a measure of instability with respect to phase separation.
        - see C.J. Bartel et al. Npj Comput Mater 5 (2019) 4. https://doi.org/10.1038/s41524-018-0143-2
        - has been quite successful in evaluating HEO stability

            Returns:
                float: Calculated decomposition enthalpy in eV/atom.
        """

        target_entry = self._find_target_entry()
        self.decomp_enthalpy = (
            self.phase_diagram.get_decomp_and_phase_separation_energy(
                entry=target_entry
            )[1]
        )

        return self.decomp_enthalpy

    def _get_decomposition_rxn(self):
        """
        Calculates the decomposition reaction corresponding to the decomposition enthalpy for
        target entry using supplied PhaseDiagram. Automatically cleans the dictionary of
        decomposition entries for a more human-readable string with relative amounts.

        These reactions are quite useful in comparing to experimental synthesis attempts and to
        determine what phases might inhibit single-phase formation of a material.

        Returns:
            str: Decomposition reaction.
        """

        target_entry = self._find_target_entry()
        decomp_entries = self.phase_diagram.get_decomp_and_phase_separation_energy(
            entry=target_entry
        )[0]

        decomp_rxn = ""

        for decomp_entry in range(len(list(decomp_entries))):
            decomp_rxn += "{:.2f}".format(list(decomp_entries.values())[decomp_entry])

            decomp_rxn += "({})".format(
                list(decomp_entries.keys())[decomp_entry].composition.reduced_formula
            )

            if decomp_entry != range(len(list(decomp_entries.keys())))[-1]:
                decomp_rxn += " + "

        self.decomp_rxn = decomp_rxn

        return self.decomp_rxn


def apply_mp2020compat(run: Vasprun) -> float:
    """
    Applies MP2020Compatbility correction scheme for GGA/GGA+U and anion mixing calculations.
    - https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/thermodynamic-stability/anion-and-gga-gga+u-mixing

    Calculation parameters/potcars should be consistent with MPRelaxSet for valid computations.
    - https://github.com/materialsproject/pymatgen/blob/master/src/pymatgen/io/vasp/MPRelaxSet.yaml

    Args:
        run (Vasprun): Pymatgen vasprun object. Used preferentially over raw energies to ensure
            scheme is implemented correctly for a given material system and calculation specs.

    Returns:
        float: Corrected energy in eV/atom
    """

    v = run.get_computed_entry()

    # get original energy in eV/atom
    energy_og = v.energy / len(v.structxure)
    print(f"original energy = {np.round(energy_og, 4)}/atom")

    # calculate corrected energy with MP2020Compatibility corrections
    v.energy_adjustments = MaterialsProject2020Compatibility().get_adjustments(v)
    energy_mp2020 = v.energy / len(v.structure)
    print(f"corrected energy = {np.round(energy_mp2020, 4)}/atom")

    return energy_mp2020
