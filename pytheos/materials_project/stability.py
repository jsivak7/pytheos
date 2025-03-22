# for evaluting the stability of a material
# TODO add in mixing energy here via materials propject data base, but not 100% how to do this yet...

from pymatgen.analysis.phase_diagram import PhaseDiagram, PDEntry


class MaterialStability:
    """
    Class for evaluting the stability of a material.

    Useful primarily for calculating relative stability of your target material with respect
    to the other calculations in your phase diagram (i.e. usually from the Materials Project database).

    Attributes:
        phase_diagram (PhaseDiagram): Pymatgen PhaseDiagram class with added target entry.
        target_entry_name (str): Name of target entry. Defaults to "my_PDEntry".
        form_energy (float): Calculated formation energy of target material in eV/atom.
        decomp_energy (float): Calculated decomposition energy of target material in eV/atom.
        decomp_rxn (str): Decomposition reaction corresponding to the calculated decomposition energy.
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
        self.form_energy = self._get_formation_energy()
        self.decomp_energy = self._get_decomposition_energy()
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

    def _get_formation_energy(self):
        """
        Calculates the formation energy of target entry using supplied PhaseDiagram.

        The formation energy is the energy relative to the elemental references.

        Returns:
            float: Calculated formation energy in eV/atom.
        """

        form_energy = self.phase_diagram.get_form_energy_per_atom(
            entry=self._find_target_entry()
        )

        return form_energy

    def _get_decomposition_energy(self):
        """
        Calculates the decomposition energy of target entry using supplied PhaseDiagram.

        The decomposition energy is a measure of instability with respect to phase separation.
        - see C.J. Bartel et al. Npj Comput Mater 5 (2019) 4. https://doi.org/10.1038/s41524-018-0143-2
        - has been quite successful in evaluating HEO stability

            Returns:
                float: Calculated decomposition energy in eV/atom.
        """

        decomp_energy = self.phase_diagram.get_decomp_and_phase_separation_energy(
            entry=self._find_target_entry()
        )[1]

        return decomp_energy

    def _get_decomposition_rxn(self):
        """
        Calculates the decomposition reaction corresponding to the decomposition energy for
        target entry using supplied PhaseDiagram. Automatically cleans the dictionary of
        decomposition entries for a more human-readable string with relative amounts.

        These reactions are quite useful in comparing to experimental synthesis attempts and to
        determine what phases might inhibit single-phase formation of a high-entropy material.

        Returns:
            str: Decomposition reaction.
        """
        decomp_entries = self.phase_diagram.get_decomp_and_phase_separation_energy(
            entry=self._find_target_entry()
        )[0]

        decomp_rxn = ""

        for decomp_entry in range(len(list(decomp_entries))):
            decomp_rxn += "{:.2f}".format(list(decomp_entries.values())[decomp_entry])

            decomp_rxn += "({})".format(
                list(decomp_entries.keys())[decomp_entry].composition.reduced_formula
            )

            if decomp_entry != range(len(list(decomp_entries.keys())))[-1]:
                decomp_rxn += " + "

        return decomp_rxn
