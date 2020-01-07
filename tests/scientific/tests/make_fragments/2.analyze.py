#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  cartesian_relax/2.analyze.py
## @brief this script is part of make_fragments scientific test
## @author Danny Farrell

import os, sys, subprocess, math
import numpy as np
import json

from concurrent.futures import ProcessPoolExecutor as PPE
from concurrent.futures import wait
import benchmark
from benchmark.util import quality_measures as qm
from benchmark.util.fragments import Fragments
from benchmark.util.pdb_io_utils import THREE2ONE, ca_pdb_reader


benchmark.load_variables()  # Python black magic: load all variables saved by previous script into
config = benchmark.config()


def get_rmsds(CAs, fragments, ca_olcs, sequence: str):
    rmsds = []
    for i, residue_fragments in enumerate(fragments.per_residue_fragments):
        residue_numbers = [x+1 for x in range(i, i+9)]
        sub_sequence = sequence[i:i+9]
        if len(sub_sequence) < 9:
            continue
        if not all([x in CAs for x in residue_numbers]):
            # Target pdb missing CA for this to work the residue numbering must MATCH
            rmsds.append([])
            continue
        pdb_seq = ''.join(ca_olcs[x] for x in residue_numbers)
        if len(pdb_seq) == 9:
            assert pdb_seq == sub_sequence  # prove that pdb is aligned

        working_CAs = np.asarray([CAs[x] for x in residue_numbers])
        centered_CAs = working_CAs - working_CAs.mean(axis=0)

        sub_rmsds = []
        for j, fragment in enumerate(residue_fragments):
            centered_fragment = np.asarray(fragment.xyzs)
            centered_fragment -= centered_fragment.mean(axis=0)
            _, rmsd = qm.kabsch_align(centered_CAs, centered_fragment)
            sub_rmsds.append(rmsd)
        rmsds.append(sub_rmsds)
    return rmsds


def get_rmsds_wrapper(db_name, target_name):
    with open("1.checkpoint.json") as fh:
        json_data = json.loads(fh.read())
    target_data = json_data[db_name][target_name]
    workdir = target_data["workdir"]
    pdb_file = target_data["pdb_file"]
    sequence = target_data['sequence']
    fragments = Fragments()
    fragments.parse_fragment_file(target_data['final_files'][-1])
    with open(pdb_file) as fh:
        pdb_lines = fh.read().split('\n')
    ca_pdb, ca_olcs = ca_pdb_reader(pdb_lines)
    working_rmsds = get_rmsds(ca_pdb, fragments, ca_olcs, sequence)
    return working_rmsds


checkpoint_filename = "2.checkpoint.json"

with open("1.checkpoint.json", 'w') as fh:
    fh.write(json.dumps(json_data))

if not os.path.isfile(checkpoint_filename):
    with PPE(10) as ppe:
        futures = []
        for db_name, db_data in json_data.items():
            for target_name, target_data in db_data.items():
                working_rmsds = ppe.submit(get_rmsds_wrapper, db_name, target_name)
                futures.append(working_rmsds)
                target_data['rmsds'] = working_rmsds
        wait(futures)
        for db_name, db_data in json_data.items():
            for target_name, target_data in db_data.items():
                    target_data['rmsds'] = target_data['rmsds'].result()

    with open(checkpoint_filename, 'w') as fh:
        fh.write(json.dumps(json_data))
else:
    with open(checkpoint_filename, 'w') as fh:
        json_data = json.loads(fh.read())


benchmark.save_variables('debug targets nstruct working_dir testname results results_file failures json_data')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
