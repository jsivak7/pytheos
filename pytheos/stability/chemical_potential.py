# for evaluating stability using the chemical potential

from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core import Element
from pandas import DataFrame
import matplotlib.pyplot as plt
from pymatgen.util.string import latexify


class ChemPotDiagram:
    """Currently limited to two-dimensional chemical potential diagrams."""

    def __init__(
        self,
        phase_diagram: PhaseDiagram,
        target_compound: str,
        cation: str,
        anion: str = "O",
    ):
        self.phase_diagram = phase_diagram
        self.target_compound = target_compound
        self.cation = cation
        self.anion = anion

        self._check_target_is_stable()

        self.all_ranges: DataFrame = None
        self.target_ranges: DataFrame = None
        self.target_anion_range: tuple = None

    def _check_target_is_stable(self) -> None:
        """
        Checks to ensure target compound is stable in supplied phase diagram.

        Raises:
            Exception: If target compound is not stable.
        """

        stable_entries = self.phase_diagram.stable_entries

        for stable_entry in stable_entries:
            if stable_entry.composition.reduced_formula == self.target_compound:
                break

        else:
            raise Exception(
                f"Target compound ({self.target_compound}) is not stable in supplied phase diagram, therefore a chemical potential diagram cannot be constructed."
            )

    def _get_range_map(self, element: str) -> dict:
        """
        Gets chemical potential range map along specific element axis for supplied phase diagram.

        Args:
            element (str): Element symbol to be considered independent variable.
                Example: if you wanted to get the ranges for Co-O along Oxygen chemical potential,
                you would provide "O" here.

        Returns:
            dict: {entry: [simplices]}
        """
        range_map = self.phase_diagram.get_chempot_range_map(
            elements=[Element(element)]
        )

        return range_map

    def get_target_anion_range(self) -> tuple:

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

    def get_all_ranges(self) -> DataFrame:

        ranges_dict = {
            "formula": [],
            f"{self.cation} (eV)": [],
            f"{self.anion} (eV)": [],
        }

        for element in [self.cation, self.anion]:
            range_map = self._get_range_map(element=element)

            for entry in range_map:
                entry_formula = entry.composition.reduced_formula

                if element == self.cation:  # only need to do for one element
                    ranges_dict["formula"] += 2 * [entry_formula]

                ranges_dict[f"{element} (eV)"].append(range_map[entry][0].coords[0][0])
                ranges_dict[f"{element} (eV)"].append(range_map[entry][0].coords[1][0])

        ranges_dict = self._add_elemental_ranges(ranges_dict=ranges_dict)

        ranges_df = DataFrame(ranges_dict)
        ranges_df.sort_values(by=f"{self.cation} (eV)", inplace=True)
        ranges_df.reset_index(drop=True, inplace=True)

        self.all_ranges = ranges_df

        self._get_target_ranges()

        return self.all_ranges

    def _add_elemental_ranges(self, ranges_dict: dict) -> dict:

        ### add in elemental data points since not included in range map explicitly
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

        self.target_ranges = self.all_ranges[
            self.all_ranges["formula"] == self.target_compound
        ]

        return self.target_ranges

    def plot_diagram(
        self,
        output_path,
        with_target: bool = True,
    ):

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

        plt.plot(
            self.all_ranges[f"{self.cation} (eV)"],
            self.all_ranges[f"{self.anion} (eV)"],
            marker="s",
            linewidth=2,
            color="k",
        )

        if with_target:
            plt.plot(
                self.target_ranges[f"{self.cation} (eV)"],
                self.target_ranges[f"{self.anion} (eV)"],
                marker="s",
                linewidth=2,
                color="tab:green",
            )

            plt.text(
                x=self.all_ranges[f"{self.cation} (eV)"].to_list()[-1] + 0.5,
                y=self.all_ranges[f"{self.anion} (eV)"][1] + 0.5,
                s=rf"{latexify(self.target_compound)}",
                color="tab:green",
                ha="right",
            )

        plt.xlim(
            self.all_ranges[f"{self.cation} (eV)"][1] - 1,
            self.all_ranges[f"{self.cation} (eV)"].to_list()[-1] + 1,
        )
        plt.ylim(
            self.all_ranges[f"{self.anion} (eV)"].to_list()[-2] - 1,
            self.all_ranges[f"{self.anion} (eV)"][1] + 1,
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
        plt.savefig(fname=output_path)


def _calc_midpoint(point1, point2) -> tuple:
    """
    Calculates the midpoint of a line segment.

    Args:
        point1: A tuple representing the first endpoint (x1, y1).
        point2: A tuple representing the second endpoint (x2, y2).

    Returns:
        A tuple representing the midpoint (x_m, y_m).
    """
    x1, y1 = point1
    x2, y2 = point2
    x_m = (x1 + x2) / 2
    y_m = (y1 + y2) / 2
    return (x_m, y_m)
