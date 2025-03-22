# for querying data from the MP database
# API info -> https://next-gen.materialsproject.org/api

from mp_api.client import MPRester
import pprint
from dotenv import load_dotenv
import os


def _load_api_key():
    """
    Uses `python-dotenv` package (https://github.com/theskumar/python-dotenv) to automatically
    load the user's Materials Project API key assuming the .env file has been set up and exists
    somewhere in the `pytheos` source files.

    To set this correctly up see docs/general_tips.md.

    Raises:
        ValueError: If MP_API_KEY cannot be loaded.

    Returns:
        str: User's MP API key for accessing the Materials Project API.
    """

    load_dotenv()

    mp_api_key = os.getenv("MP_API_KEY")

    if mp_api_key:
        pass
    else:
        raise FileNotFoundError(
            """Materials Project API key could not be loaded. Have you set up a `.env` file??"""
        )

    return mp_api_key


def get_db_version() -> str:
    """
    Gets Materials Project database version.

    *Good practice to note this whenever a query is made for reproducibility.*

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
    thermo_type: str = "R2SCAN",
    additional_criteria: dict = None,
) -> list:
    """
    Gets a list of ComputedStructureEntries across entire chemical system for a given set of elements.

    For example, if given elements = ["Mg", "Co", "O"] this will return a list of all entries in the
    parent Mg-Co-O chemical system, as well as all the subsystems (all MgxOy, CoxOy, MgxCoy, MgxCoyOz,
    Mg, Co and O phases). Very useful for create phase diagrams across the entire chemical space
    (and then calculating relative stability using say decomposition enthalpy).

    *One should ideally supply elements argument using `my_PDEntry.composition.chemical_system` to ensure that
    the correct chemical system is queried.*

    Straightforward implementation of mp_api.client.MPRester.get_entries_in_chemsys
    (https://github.com/materialsproject/api/blob/main/mp_api/client/mprester.py).

    Internally handles getting user's MPAPIKey via `_load_api_key()` function.

    Args:
        elements (list): List of elements for chemical system to get entries. Example -> ["Mg", "Co", "O"].
            Preferentially use -> my_PDEntry.composition.chemical_system
        thermo_type (str, optional): Exchange-Correlation functional used to calculate entries.
            Options are ["GGA_GGA+U", "R2SCAN"]. Defaults to "R2SCAN".
        additional_criteria (dict, optional): Additional criteria for entry search. Defaults to None.

    Returns:
        list: list of ComputedStructureEntries
    """

    apikey = _load_api_key()

    criteria = {"thermo_types": [thermo_type]}

    if additional_criteria:
        criteria.update(additional_criteria)

    with MPRester(api_key=apikey) as mpr:
        entries = mpr.get_entries_in_chemsys(
            elements=elements,
            additional_criteria=criteria,
        )

    return entries


def get_entries_for_formula(
    formula: str,
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

    Internally handles getting user's MPAPIKey via `_load_api_key()` function.

    Args:
        formula (str): Chemical formula to get entries. Example -> "CoO"
            One can also give an mpid. Example -> "mp-1265-r2SCAN"
        thermo_type (str, optional): Exchange-Correlation functional used to calculate entries.
            Options are ["GGA_GGA+U", "R2SCAN"]. Defaults to "R2SCAN".
        additional_criteria (dict, optional): Additional criteria for entry search. Defaults to None.

    Returns:
        list: list of ComputedStructureEntries
    """
    apikey = _load_api_key()

    criteria = {"thermo_types": [thermo_type]}

    if additional_criteria:
        criteria.update(additional_criteria)

    with MPRester(api_key=apikey) as mpr:
        entries = mpr.get_entries(
            chemsys_formula_mpids=formula,
            additional_criteria=criteria,
        )

    return entries
