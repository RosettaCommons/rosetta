#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file   scorefile_io_utils.py
## @brief  Pandas utils for conversion and analysis of scorefiles
## @author Morgan Nance (@mlnance)

import pandas as pd

def scorefile_to_dataframe(scorefile):
    """Turn a scorefile into a Pandas DataFrame object
    """
    with open(scorefile, "r") as fh:
        scorefile_lines = fh.readlines()
    headers = scorefile_lines[1].split()[1:]
    score_data = [line.split()[1:] for line in scorefile_lines[2:]]

    return pd.DataFrame( score_data, columns=headers, dtype=float )
