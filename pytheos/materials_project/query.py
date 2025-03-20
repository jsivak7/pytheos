# functions to query data from the MP database
# API info -> https://next-gen.materialsproject.org/api

from mp_api.client import MPRester
import pprint


def get_db_version() -> str:
    """
    Gets Materials Project database version - good to ensure reproducibility.

    Returns:
        str: Materials Project database version (YYYY_MM_DD)
    """

    with MPRester() as mpr:

        mpdb_version = mpr.get_database_version()

        return mpdb_version


def print_query_fields() -> None:
    """
    Prints all available fields for an MP API query in alphabetical order.
    """

    with MPRester() as mpr:
        available_fields = mpr.materials.summary.available_fields
        available_fields.sort()
        pprint.pprint(available_fields)


def get_entries_across_chemsys(
    elements: list,
    MPApiKey: str,
    thermo_type: str = "R2SCAN",
    additional_criteria: dict = None,
) -> list:
    """
    Gets a list of ComputedStructureEntries across entire chemical system for a given set of elements.

    For example, if given elements = ["Mg", "Co", "O"] this will return a list of all entries in the
    parent Mg-Co-O chemical system, as well as all the subsystems (all MgxOy, CoxOy, MgxCoy, MgxCoyOz,
    Mg, Co and O phases). Very useful for create phase diagrams across the entire chemical space
    (and then calculating relative stability using say decomposition enthalpy).

    Straightforward implementation of mp_api.client.MPRester.get_entries_in_chemsys
    (https://github.com/materialsproject/api/blob/main/mp_api/client/mprester.py).

    Args:
        elements (list): List of elements for chemical system to get entries. Example -> ["Mg", "Co", "O"].
        MPApiKey (str): Materials Project API Key - user specific.
        thermo_type (str, optional): Exchange-Correlation functional used to calculate entries.
            Options are ["GGA_GGA+U", "R2SCAN"]. Defaults to "R2SCAN".
        additional_criteria (dict, optional): Additional criteria for entry search. Defaults to None.

    Returns:
        list: list of ComputedStructureEntries
    """

    criteria = {"thermo_types": [thermo_type]}

    if additional_criteria:
        criteria.update(additional_criteria)

    with MPRester(MPApiKey) as mpr:
        entries = mpr.get_entries_in_chemsys(
            elements=elements,
            additional_criteria=criteria,
        )

    return entries


def get_entries_for_formula(
    formula: str,
    MPApiKey: str,
    thermo_type: str = "R2SCAN",
    additional_criteria: dict = None,
) -> list:
    """
    Gets a list of ComputedStructureEntries for a given chemical formula.

    For example, if given formula = "CoO" this will return a list of all entries with this
    chemical formula. This is especially useful when paired with for example
    additional_criteria = {"is_stable": True}, which will only return a single ComputedStructureEntry
    in the list if this formula is stable (on the hull) for the thermo_type specified,
    and an empty list otherwise.

    Straightforward implementation of mp_api.client.MPRester.get_entries.
    (https://github.com/materialsproject/api/blob/main/mp_api/client/mprester.py).

    Args:
        formula (str): Chemical formula to get entries. Example -> "CoO"
            One can also give an mpid. Example -> "mp-1265-r2SCAN"
        MPApiKey (str): Materials Project API Key - user specific.
        thermo_type (str, optional): Exchange-Correlation functional used to calculate entries.
            Options are ["GGA_GGA+U", "R2SCAN"]. Defaults to "R2SCAN".
        additional_criteria (dict, optional): Additional criteria for entry search. Defaults to None.

    Returns:
        list: list of ComputedStructureEntries
    """
    criteria = {"thermo_types": [thermo_type]}

    if additional_criteria:
        criteria.update(additional_criteria)

    with MPRester(MPApiKey) as mpr:
        entries = mpr.get_entries(
            chemsys_formula_mpids=formula,
            additional_criteria=criteria,
        )

    return entries
