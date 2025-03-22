from pymatgen.analysis.phase_diagram import PhaseDiagram
from pymatgen.core import Element
from pandas import DataFrame


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
            element (str): Element to be considered independent variable.

        Returns:
            dict: {entry: [simplices]}
        """
        range_map = self.phase_diagram.get_chempot_range_map(
            elements=[Element(element)]
        )

        return range_map

    def get_stable_anion_bounds(self) -> tuple:

        anion_name = Element(self.anion).long_name

        range_map = self._get_range_map(element=self.anion)

        for entry in range_map:

            entry_formula = entry.composition.reduced_formula

            if entry_formula == self.target_compound:

                self.stable_anion_min = range_map[entry][0].coords.min()
                self.stable_anion_max = range_map[entry][0].coords.max()
                self.stable_anion_distance = abs(self.minimum - self.maximum)

        return (
            self.stable_anion_min,
            self.stable_anion_max,
            self.stable_anion_distance,
        )

    def get_all_regions(self) -> DataFrame:

        regions = {
            "formula": [],
            f"{self.cation} (eV)": [],
            f"{self.anion} (eV)": [],
        }

        for element in [self.cation, self.anion]:
            range_map = self._get_range_map(element=element)

            for entry in range_map:
                entry_formula = entry.composition.reduced_formula

                # only add formulas for one element
                if element == self.cation:
                    regions["formula"] += 2 * [entry_formula]

                regions[f"{element} (eV)"].append(range_map[entry][0].coords[0][0])
                regions[f"{element} (eV)"].append(range_map[entry][0].coords[1][0])

        regions = self._add_elemental_regions(regions=regions)

        df = DataFrame(regions)
        df.sort_values(by=f"{self.cation} (eV)", inplace=True)
        df.reset_index(drop=True, inplace=True)
        self.all_regions = df

    def _add_elemental_regions(self, regions: dict) -> dict:

        ### add in elemental data points since not included in range map explicitly
        # cation
        regions["formula"] += 2 * [self.cation]
        regions[f"{self.cation} (eV)"].append(0)
        regions[f"{self.cation} (eV)"].append(0)
        regions[f"{self.anion} (eV)"].append(min(regions[f"{self.anion} (eV)"]))
        regions[f"{self.anion} (eV)"].append(-50)

        # anion
        regions["formula"] += 2 * [self.anion]
        regions[f"{self.anion} (eV)"].append(0)
        regions[f"{self.anion} (eV)"].append(0)
        regions[f"{self.cation} (eV)"].append(min(regions[f"{self.cation} (eV)"]))
        regions[f"{self.cation} (eV)"].append(-50)

        return regions

    def _find_target_region(self, data_frame: DataFrame) -> DataFrame:

        target_df = data_frame[data_frame["formula"] == self.target_compound]

        return target_df


"""

import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"] = (3, 4)
plt.rcParams["font.size"] = 12

plt.plot(
    df[f"{self.cation} (eV)"],
    df[f"{self.anion} (eV)"],
    marker="s",
    linewidth=2,
    color="k",
)
plt.plot(
    target_df[f"{self.cation} (eV)"],
    target_df[f"{self.anion} (eV)"],
    marker="s",
    linewidth=2,
    color="tab:green",
)
plt.xlim(df[f"{self.cation} (eV)"][1] - 1, df[f"{self.cation} (eV)"].to_list()[-1] + 1)
plt.ylim(-7.75, 0.75)

degrees_symbol = "o"

plt.xlabel(rf"μ$_{{\rm {self.cation}}}$ - μ$_{{\rm {self.cation}}}^{{\rm o}}$ (eV)")
plt.ylabel(
    rf"μ$_{{\rm {self.anion}}}$ - μ$_{{\rm {self.anion}}}^{{\rm o}}$ (eV)"
)

plt.title(f"{self.target_compound}")

plt.show()"""
