# for evaluating material stability using chemical potential diagrams

from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core import Element
from pandas import DataFrame
import matplotlib.pyplot as plt
from matplotlib.pyplot import Figure
import numpy as np
from pymatgen.util.string import latexify
import warnings


class ChemPotDiagram:
    """
    Class for constructing and evaluating chemical potential diagrams for a given cation-anion formula space.

    The chemical potential diagram is the mathemetical dual to convex hull.
    https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/chemical-potential-diagrams-cpds

    Current implementation is limited to two-dimensional (2D) chemical potential diagrams,
    though these can be generated in much higher dimensions as well if desired. You do sacrifice
    human-interprettability however as you increase dimensions though...

    The original 1999 Yokokawa paper on chemical potential diagrams provides a wealth of
    information: Yokokawa, H.  JPE 20, 258 (1999). https://doi.org/10.1361/105497199770335794

    The Materials Project team has utilized these diagrams in the following paper as well:
    Todd, P. K. et al. Journal of the American Chemical Society, 143(37), 15185-15194.
        https://doi.org/10.1021/jacs.1c06229

    Attributes:
        phase_diagram (PhaseDiagram): Pymatgen PhaseDiagram class.
        cation (str): Cation of interest.
        anion (str): Anion of interest. Defaults to "O".
        target_compound (str): Formula of target compound. e.g. "NiO" or "ZrO2". Defaults to None.
        all_stable_ranges (DataFrame): Pandas DataFrame for all stable chemical potential ranges in the format
        {"formula": values, "cation (eV)": values, "anion (eV)": values}. Defaults to None.
        target_ranges (DataFrame): Pandas DataFrame for stable chemical potential ranges for target compound.
            Same DataFrame format as `all_stable_ranges`. Defaults to None.
        target_anion_range (tuple): Bounds along anion chemical potential axis where target compound is stable.
            Takes the format (min, max, distance) - all in eV. Defaults to None.
        diagram (Figure): Matplotlib figure object of 2D chemical potential diagram with x-axis being cation
            and y-axis being anion. Defaults to None.
    """

    def __init__(
        self,
        phase_diagram: PhaseDiagram,
        cation: str,
        anion: str = "O",
        target_compound: str = None,
    ):
        """
        Args:
            phase_diagram (PhaseDiagram): Pymatgen PhaseDiagram class.
            cation (str): Cation of interest.
            anion (str, optional): Anion of interest. Defaults to "O".
            target_compound (str, optional): Formula of target compound. e.g. "NiO" or "ZrO2".
                Defaults to None.
        """

        self.phase_diagram = phase_diagram
        self.target_compound = target_compound
        self.cation = cation
        self.anion = anion

        if target_compound:
            self._check_target_is_stable()

        self.all_stable_ranges: DataFrame = None
        self.target_ranges: DataFrame = None
        self.target_anion_range: tuple = None
        self.diagram: Figure = None

    def _check_target_is_stable(self) -> None:
        """
        Checks to ensure target compound is stable in supplied phase diagram.

        Automatically called during ChemPotDiagram initialization if `target_compound` is supplied.
        """

        stable_entries = self.phase_diagram.stable_entries

        for stable_entry in stable_entries:
            if stable_entry.composition.reduced_formula == self.target_compound:
                break

        else:
            warnings.warn(
                f"WARNING!!! \nTarget compound ({self.target_compound}) is not stable in the supplied phase diagram, therefore it will not exist on the chemical potential diagram."
            )

    def _get_range_map(self, element: str) -> dict:
        """
        Gets chemical potential range map along specific element axis for supplied phase diagram.

        Args:
            element (str): Element symbol to be considered independent variable.
                e.g. if you wanted to get the ranges for Co-O along Oxygen chemical potential,
                you would provide "O" here.

        Returns:
            dict: {entry: [simplices]}
        """

        range_map = self.phase_diagram.get_chempot_range_map(
            elements=[Element(element)]
        )

        return range_map

    def get_target_anion_range(self) -> tuple:
        """
        Extracts the bounds (max, min, ditance) in anion chemical potential space for which
            the target compound is stable.

        Raises:
            Exception: If target compound was not supplied to ChemPotDiagram instance.

        Returns:
            tuple: (minimum, maximum, distance) in eV
        """

        if self.target_compound is None:
            raise Exception("A target compound has not been supplied.")

        range_map = self._get_range_map(element=self.anion)

        for entry in range_map:
            entry_formula = entry.composition.reduced_formula

            if entry_formula == self.target_compound:
                target_anion_min = range_map[entry][0].coords.min()
                target_anion_max = range_map[entry][0].coords.max()
                target_anion_distance = abs(target_anion_min - target_anion_max)

                self.target_anion_range = (
                    target_anion_min,
                    target_anion_max,
                    target_anion_distance,
                )

        return self.target_anion_range

    def get_all_stable_ranges(self) -> DataFrame:
        """
        Extract all ranges within constructed chemical potential diagram.

        NOTE that there are two data points for each composition, as the composition
        is stable for the line drawn between these two points in chemical potential space.

        Returns:
            DataFrame: {
                "formula": values,
                "cation (eV)": values,
                "anion (eV)": values,
                }
        """

        ranges_dict = {
            "formula": [],
            f"{self.cation} (eV)": [],
            f"{self.anion} (eV)": [],
        }

        for element in [self.cation, self.anion]:
            range_map = self._get_range_map(element=element)

            for entry in range_map:
                entry_formula = entry.composition.reduced_formula

                if element == self.cation:  # only get formulas for one element
                    ranges_dict["formula"] += 2 * [entry_formula]

                ranges_dict[f"{element} (eV)"].append(range_map[entry][0].coords[0][0])
                ranges_dict[f"{element} (eV)"].append(range_map[entry][0].coords[1][0])

        ranges_dict = self._add_elemental_ranges(ranges_dict=ranges_dict)

        ranges_df = DataFrame(ranges_dict)
        ranges_df.sort_values(by=f"{self.cation} (eV)", inplace=True)
        ranges_df.reset_index(drop=True, inplace=True)

        self.all_stable_ranges = ranges_df

        if self.target_compound is not None:
            self._get_target_ranges()

        return self.all_stable_ranges

    def _add_elemental_ranges(self, ranges_dict: dict) -> dict:
        """
        Add in elemental ranges to compound ranges as elements are not explicitly
        included in the `get_all_stable_ranges` method, and therefore need to be added afterwards.

        NOTE that the current implementation fixes the minimum value for elements at -50 eV.

        Args:
            ranges_dict (dict): Dictionary of the following form, but containing only compounds:
                {
                    "formula": values,
                    "cation (eV)": values,
                    "anion (eV)": values,
                }

        Returns:
            dict: Dictionary of compounds with elements
        """

        # cation
        ranges_dict["formula"] += 2 * [self.cation]
        ranges_dict[f"{self.cation} (eV)"].append(0)
        ranges_dict[f"{self.cation} (eV)"].append(0)
        ranges_dict[f"{self.anion} (eV)"].append(min(ranges_dict[f"{self.anion} (eV)"]))
        ranges_dict[f"{self.anion} (eV)"].append(-50)

        # anion
        ranges_dict["formula"] += 2 * [self.anion]
        ranges_dict[f"{self.anion} (eV)"].append(0)
        ranges_dict[f"{self.anion} (eV)"].append(0)
        ranges_dict[f"{self.cation} (eV)"].append(
            min(ranges_dict[f"{self.cation} (eV)"])
        )
        ranges_dict[f"{self.cation} (eV)"].append(-50)

        return ranges_dict

    def _get_target_ranges(self) -> DataFrame:
        """
        Gets chemical potential ranges for target compound.

        Automatically called during `get_all_stable_ranges()` method if target compound is supplied.

        Raises:
            Exception: If target compound was not supplied to ChemPotDiagram instance.

        Returns:
            DataFrame: {
                "formula": values,
                "cation (eV)": values,
                "anion (eV)": values,
                }
        """

        if self.target_compound is None:
            raise Exception("A target compound has not been supplied.")

        self.target_ranges = self.all_stable_ranges[
            self.all_stable_ranges["formula"] == self.target_compound
        ]

        return self.target_ranges

    def plot_diagram(
        self,
        with_target: bool = True,
    ) -> Figure:
        """
        Plots a generic 2D chemical potential diagram across relevant chemical space using supplied ranges
        with cation chemical potential on the x-axis and anion on the y-axis.

        Args:
            with_target (bool, optional): Overlays stable region for target compound in green on diagram.
                Defaults to True.

        Raises:
            Exception: If chemical potential ranges have not yet been obtained.
            Exception: If with_target = True, but target compound has not been supplied.

        Returns:
            Figure: Matplotlib Figure object of 2D chemical potential diagram.
        """

        plt.rcParams["figure.dpi"] = 400
        plt.rcParams["figure.figsize"] = (3, 4)
        plt.rcParams["font.size"] = 12
        plt.rcParams["font.serif"] = "Helvetica"
        plt.rcParams["axes.linewidth"] = 1.65
        plt.rcParams["xtick.direction"] = "in"
        plt.rcParams["ytick.direction"] = "in"
        plt.rcParams["xtick.major.width"] = 1.65
        plt.rcParams["ytick.major.width"] = 1.65
        plt.rcParams["font.weight"] = "bold"

        plt.clf()
        fig = plt.figure()

        if self.all_stable_ranges is None:
            raise Exception(
                "The chemical potential ranges need to be obtained before making a diagram."
            )

        plt.plot(
            self.all_stable_ranges[f"{self.cation} (eV)"],
            self.all_stable_ranges[f"{self.anion} (eV)"],
            marker="s",
            linewidth=2,
            color="k",
        )

        if with_target:
            if self.target_compound is None:
                raise Exception(
                    "A target compound has not been supplied. \nSet `with_target` to `False` if you do not want it included in the diagram."
                )
            else:
                plt.plot(
                    self.target_ranges[f"{self.cation} (eV)"],
                    self.target_ranges[f"{self.anion} (eV)"],
                    marker="s",
                    linewidth=2,
                    color="tab:green",
                )

                plt.text(
                    x=0.92,
                    y=0.92,
                    s=rf"{latexify(self.target_compound)}",
                    transform=plt.gca().transAxes,
                    color="tab:green",
                    ha="right",
                )

        ticks = np.arange(-50, 2, 2)
        plt.xticks(ticks)
        plt.yticks(ticks)

        plt.xlim(
            self.all_stable_ranges[f"{self.cation} (eV)"][1] - 0.75,
            self.all_stable_ranges[f"{self.cation} (eV)"].to_list()[-1] + 0.75,
        )
        plt.ylim(
            self.all_stable_ranges[f"{self.anion} (eV)"].to_list()[-2] - 0.75,
            self.all_stable_ranges[f"{self.anion} (eV)"][1] + 0.75,
        )

        plt.xlabel(
            rf"μ$_{{\rm \bf {self.cation}}}$ - μ$_{{\rm \bf {self.cation}}}^{{\rm \bf o}}$ (eV)",
            weight="bold",
        )
        plt.ylabel(
            rf"μ$_{{\rm \bf {self.anion}}}$ - μ$_{{\rm \bf {self.anion}}}^{{\rm \bf o}}$ (eV)",
            weight="bold",
        )

        plt.title(f"{self.cation} - {self.anion}", weight="bold")
        plt.tight_layout()

        self.diagram = fig
        plt.close()

        return self.diagram
