# general MP utilities


def print_available_query_fields(
    MPApiKey: str,
) -> None:
    """
    Prints all available fields for an MP API query in alphabetical order.

    Args:
        MPApiKey (str): Materials Project API Key (https://next-gen.materialsproject.org/api) - user specific
    """
    from mp_api.client import MPRester
    import pprint

    with MPRester(MPApiKey) as mpr:
        available_fields = mpr.materials.summary.available_fields
        available_fields.sort()
        pprint.pprint(available_fields)
