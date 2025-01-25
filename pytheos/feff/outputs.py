# For parsing FEFF output files


def xmu_dat_to_df(
    xmu_path: str = "xmu.dat",
):
    """
    Converts the xmu.dat output from FEFF to a Pandas dataframe object.

    Args:
        xmu_path (str, optional): Relative path to xmu.dat data file. Defaults to "xmu.dat".

    Returns:
        pd.DataFrame: Pandas dataframe object of FEFF output data from xmu.dat.
    """
    from pymatgen.io.feff import outputs
    import pandas as pd

    print(f"Reading {xmu_path} as Pandas DataFrame")

    xmu = outputs.Xmu.from_file(xmu_dat_file=xmu_path)

    edge = xmu.edge
    absorbing_atom = xmu.absorbing_atom
    material = xmu.material_formula
    calc_type = xmu.calc

    xmu_df = pd.DataFrame(
        {
            "energy_eV": xmu.energies,
            "wavenumber": xmu.wavenumber,
            "mu": xmu.mu,
            "mu0": xmu.mu0,
            "chi": xmu.chi,
        }
    )
    return xmu_df
