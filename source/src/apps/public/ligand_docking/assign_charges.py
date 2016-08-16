#!/usr/bin/env python
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

import sys
from openeye.oechem import * # OpenEye
try:
    from openeye.oequacpac import * # OpenEye, newer releases
except ImportError:
    from openeye.oeproton import * # OpenEye, older releases

ifs = oemolistream()
ifs.open() # stdin
ofs = oemolostream()
ofs.open() # stdout

ifs.SetFormat(OEFormat_MOL2)

# Write MOL2 without renaming atoms to match their input ordering
ofs.SetFormat(OEFormat_MOL2)
ofs.SetFlavor(OEFormat_MOL2, (OEOFlavor_Generic_Default | OEOFlavor_MOL2_Default) & ~OEOFlavor_MOL2_AtomNames)

first = True
for m in ifs.GetOEGraphMols():
    mol = OEGraphMol(m) # otherwise same mol obj is reused
    if first:
        first = False
        OEClearPartialCharges(mol)
        OEAssignPartialCharges(mol, OECharges_AM1BCC)
    OEWriteMolecule(ofs, mol)

