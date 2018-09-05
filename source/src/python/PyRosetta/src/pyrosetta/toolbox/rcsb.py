#!/usr/bin/env python
# :noTabs=true:


# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   rcsb.py
## @brief
## @author Evan H. Baugh, Johns Hopkins University

from __future__ import print_function

import os, sys
import urllib

from contextlib import contextmanager

if sys.version_info.major > 2:
    import urllib.request

    urllib_urlretrieve = urllib.request.urlretrieve
else:
    urllib_urlretrieve = urllib.urlretrieve


@contextmanager
def download_from_web(url):
    """Get content from the web and make it available within a
    context manager.
    """
    try:
        temp = urllib_urlretrieve(url)[0]
    except:
        raise IOError(
            "Cannot access network resource, please check your Internet access."
        )
    else:
        if os.path.getsize(temp) > 150:
            # Arbitrarily 150... else url was invalid.
            with open(temp) as f:
                file_data = f.readlines()
        else:
            raise IOError("Invalid URL")

        yield file_data
        os.remove(temp)  # Remove temp file.


def load_from_rcsb(pdb_code, pdb_filename=None):
    """Downlaod a PDB file from RCSB and write it to disk.

    Args:
        pdb_code (str): The four-letter accession code of the desired PDB file
        pdb_filename (str): Optional argument for output filename.
            Defaults to <pdb_code>.pdb.

    Examples:
        >>> load_from_rcsb("1YY8")
    """
    pdb_code = pdb_code.upper()
    pdb_url = "http://www.rcsb.org/pdb/files/" + pdb_code + ".pdb"

    if not pdb_filename:
        pdb_filename = pdb_code + ".pdb"

    with download_from_web(pdb_url) as pdb_data:
        with open(pdb_filename, "w") as f:
            f.writelines(pdb_data)


def pose_from_rcsb(pdb_code, ATOM=True, CRYS=False):
    """Downlaod a PDB file from RCSB, write to disk as <pdb_code>.pdb and return
    it as a pose. Optionally calls cleanATOM and/or cleanCRYS.

    Notes:
        Automatic monomer extraction assumes the PDB contains a homodimer.

    Args:
        pdb_code (str): The four-letter accession code of the desired PDB file
        ATOM (bool): Only write ATOM records to disk. Defaults to True.
        CRYS (bool): Attempt to extract a monomer from the target PDB.
            Defaults to False.

    Examples:
        >>> pose = pose_from_rcsb("1YY8")
    """
    from pyrosetta import pose_from_file
    from pyrosetta.toolbox.cleaning import cleanATOM, cleanCRYS

    pdb_code = pdb_code.upper()
    load_from_rcsb(pdb_code)
    if ATOM:
        cleanATOM(pdb_code + ".pdb")
        pdb_code = pdb_code + ".clean"
    if CRYS:
        cleanCRYS(pdb_code + ".pdb")
        pdb_code = pdb_code + ".mono"
    pose = pose_from_file(pdb_code + ".pdb")
    return pose


def load_fasta_from_rcsb(pdb_code, fasta_outfile=None):
    """Downlaod a FASTA file from RCSB and write it to disk.

    Args:
        pdb_code (str): The four-letter accession code of the desired PDB file
        fasta_outfile (str): Optional argument for output filename. Defaults to
        <pdb_code>.fasta.

    Examples:
        >>> load_fasta_from_rcsb("1YY8")
    """
    pdb_code = pdb_code.upper()
    fasta_url = (
        "https://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList="
        + pdb_code
        + "&compressionType=uncompressed"
    )

    if fasta_outfile is None:
        fasta_outfile = pdb_code + ".fasta"

    with download_from_web(fasta_url) as pdb_data:
        with open(fasta_outfile, "w") as f:
            f.writelines(pdb_data)
