# General utilities
from ase import Atoms


def read_structure_to_ASEAtoms(file_path: str):
    """
    Read in structure file to ASE Atoms Object

    Args:
        file_path (str): relative path to structure file

    Returns:
        Atoms: modular object to perform other operations
    """
    from ase import io

    s = io.read(f"{file_path}")
    print(f"{file_path} read in as ASE Atoms object")
    return s


def write_structure_from_aseAtoms(struc: Atoms, file_path: str, overwrite=False):
    """
    Write ASE Atoms object to structure file
    NOTE that always writes "direct" coordinates

    Args:
        struc (Atoms): structure to be written
        file_path (str): relative path to write structure file, include suffix for desired file type (*.vasp, *.cif, etc.)
        overwrite (bool): override to write over file already made. Defaults to False.

    Raises:
        FileExistsError: if given path for output file already exists, ensures that previously generated SQSs are not overwritten
    """
    import os
    from ase import io

    if os.path.exists(file_path) and overwrite == False:
        raise FileExistsError(file_path)

    io.write(f"{file_path}", struc, direct=True)
    print(f"ASE Atoms written to {file_path}")
