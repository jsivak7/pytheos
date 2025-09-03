# plotting tools

import matplotlib.pyplot as plt
import numpy as np
from pymatgen.core import Element
from matplotlib import cm, colors
from typing import TYPE_CHECKING
from typing import Literal


def plot_heatmap_periodic_table(
    elemental_data=None,
    cbar_label="",
    output_file: str = "periodic_table.png",
    cmap="Blues",
    cmap_range=None,
    blank_color="#F0F0F0",
    value_format="%.2f",
    max_row: int = 9,
):
    # adapted from pymatgen/src/pymatgen/util/plotting/periodic_table_heatmap()
    # wanted some specifics changed for my own plots

    """A static method that generates a heat map overlaid on a periodic table.

    Args:
        elemental_data (dict): A dictionary with the element as a key and a
            value assigned to it, e.g. surface energy and frequency, etc.
            Elements missing in the elemental_data will be grey by default
            in the final table elemental_data={"Fe": 4.2, "O": 5.0}.
        cbar_label (str): Label of the color bar. Default is "".
        cmap_range (tuple): Minimum and maximum value of the color map scale.
            If None, the color map will automatically scale to the range of the
            data.
        value_format (str): Formatting string to show values. If None, no value
            is shown. Example: "%.4f" shows float to four decimals.
        cmap (str): Color scheme of the heatmap. Default is 'YlOrRd'.
            Refer to the matplotlib documentation for other options.
        blank_color (str): Color assigned for the missing elements in
            elemental_data. Default is "grey".
        edge_color (str): Color assigned for the edge of elements in the
            periodic table. Default is "white".
        max_row (int): Maximum number of rows of the periodic table to be
            shown. Default is 9, which means the periodic table heat map covers
            the standard 7 rows of the periodic table + 2 rows for the lanthanides
            and actinides. Use a value of max_row = 7 to exclude the lanthanides and
            actinides.

    Returns:
        plt.Axes: matplotlib Axes object
    """

    # Convert primitive_elemental data in the form of numpy array for plotting.
    if cmap_range is not None:
        max_val = cmap_range[1]
        min_val = cmap_range[0]
    else:
        max_val = max(elemental_data.values())
        min_val = min(elemental_data.values())

    max_row = min(max_row, 9)

    if max_row <= 0:
        raise ValueError("The input argument 'max_row' must be positive!")

    value_table = np.empty((max_row, 18)) * np.nan
    blank_value = min_val - 0.01

    for el in Element:
        value = elemental_data.get(el.symbol, blank_value)
        if 57 <= el.Z <= 71:
            plot_row = 8
            plot_group = (el.Z - 54) % 32
        elif 89 <= el.Z <= 103:
            plot_row = 9
            plot_group = (el.Z - 54) % 32
        else:
            plot_row = el.row
            plot_group = el.group
        if plot_row > max_row:
            continue
        value_table[plot_row - 1, plot_group - 1] = value

    fig, ax = plt.subplots()
    plt.gcf().set_size_inches(8, 3.6)

    # We set nan type values to masked values (ie blank spaces)
    data_mask = np.ma.masked_invalid(value_table.tolist())
    heatmap = ax.pcolor(
        data_mask,
        cmap=cmap,
        edgecolors="white",
        linewidths=2.5,
        vmin=min_val - 0.001,
        vmax=max_val + 0.001,
    )

    cbar = fig.colorbar(heatmap)

    # Grey out missing elements in input data
    cbar.cmap.set_under(blank_color)

    # Set the color bar label and tick marks
    cbar.set_label(cbar_label, rotation=270, labelpad=25, size=10)
    cbar.ax.tick_params(labelsize=10)

    # Refine and make the table look nice
    ax.axis("off")
    ax.invert_yaxis()

    # Set the scalarmap for fontcolor
    norm = colors.Normalize(vmin=min_val, vmax=max_val)
    scalar_cmap = cm.ScalarMappable(norm=norm, cmap=cmap)

    # Label each block with corresponding element and value
    for ii, row in enumerate(value_table):
        for jj, el in enumerate(row):
            if not np.isnan(el):
                symbol = Element.from_row_and_group(ii + 1, jj + 1).symbol
                rgba = scalar_cmap.to_rgba(el)
                fontcolor = _decide_fontcolor(rgba)
                plt.text(
                    jj + 0.5,
                    ii + 0.4,
                    symbol,
                    horizontalalignment="center",
                    verticalalignment="center",
                    fontsize=10,
                    color=fontcolor,
                    weight="bold",
                )
                if el != blank_value and value_format is not None:
                    plt.text(
                        jj + 0.5,
                        ii + 0.8,
                        value_format % el,
                        horizontalalignment="center",
                        verticalalignment="center",
                        fontsize=6,
                        color=fontcolor,
                    )

    plt.tight_layout()
    plt.savefig(output_file, dpi=1200)

    return ax


def _decide_fontcolor(rgba: tuple) -> Literal["black", "white"]:
    red, green, blue, _ = rgba
    if red * 0.299 + green * 0.587 + blue * 0.114 > (186 / 255):
        return "black"

    return "white"
