# for using ionic radii

import yaml
import os


def load_ionic_radii() -> dict:
    """
    Loads `ionic_radii.yaml` file.

    Returns:
        dict: Dictionary of ionic radii
    """

    # get path to this module
    module_dir = os.path.dirname(__file__)

    with open(f"{module_dir}/ionic_radii.yaml", "r") as file:
        ionic_radii = yaml.safe_load(file)

    return ionic_radii


def search_ionic_radii(
    ionic_radii: dict,
    target_valence: str,
    target_coord: int,
    spin_state: str = "high",
) -> dict:
    """
    Searches dictionary of ionic radii for specific valence, coordination, and spin states.

    Args:
        ionic_radii (dict): Dictionary of ionic radii matching. Usually `ionic_radii.yaml` within this module.
        target_valence (str): Target valence state. E.g. "4+".
        target_coord (int): Target coordination. E.g. 6.
        spin_state (str, optional): Use high- or low-spin ionic radii (if applicable). Defaults to "high".

    Returns:
        dict: Elements and ionic radii (in Angstroms) as key:value pairs that match criteria.
    """

    print(f"Searching ionic radii with the following criteria:")
    print(f"\t{target_valence} valence")
    print(f"\t{target_coord} coordination")
    print(f"\t{spin_state} spin (if applicable)")

    target_ionic_radii = {}

    for element in ionic_radii:
        for valence in ionic_radii[element]:
            if valence == target_valence:

                for coord in ionic_radii[element][valence]["coord"]:
                    if coord == target_coord:

                        coordinations = ionic_radii[element][valence]["coord"][coord]

                        # if there is "low_spin" or "high_spin"
                        if isinstance(coordinations, dict):
                            try:
                                target_ionic_radius = coordinations[
                                    f"{spin_state}_spin"
                                ]

                            except:  # if we don't have the spin state
                                pass

                        # if no spin
                        else:
                            target_ionic_radius = coordinations

                        target_ionic_radii[element] = target_ionic_radius

    print(f"\nIdentified {len(target_ionic_radii)} ions matching criteria:")

    for radius in target_ionic_radii:
        print(f"\t{radius}\t{target_ionic_radii[radius]:.3f} Å")

    return target_ionic_radii
