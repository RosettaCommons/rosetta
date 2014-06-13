#!/usr/bin/env python
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   doxygen.py
## @brief  Extract Doxygen comments from C++ files
## @author Brian Weitzner, Sergey Lyskov

from dox_extract import dox_extract

extractor = dox_extract()

def getDoxygenComment( file_name, line_number ):
    return extractor.getDoxygenComment( file_name, line_number )
