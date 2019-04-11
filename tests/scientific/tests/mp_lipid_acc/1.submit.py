#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  _template_/1.submit.py
## @brief this script is part of the mp_lipid_acc scientific test
## @author JKLeman

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

#==> EDIT HERE
testname    = "mp_lipid_acc"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

#==> EDIT HERE
command_line = '''
-database {rosetta_dir}/database
-in:file:s {rosetta_dir}/tests/scientific/data/{testname}/db_input/{target}__tr.pdb
-mp:setup:spanfiles {rosetta_dir}/tests/scientific/data/{testname}/db_input/{target}__tr.span
-out:path:pdb {prefix}
-run:multiple_processes_writing_to_one_directory 1
-no_color
'''.replace('\n', ' ').replace('  ', ' ')

#==> EDIT HERE
nstruct = 1

#==> EDIT HERE
# targets are sorted by decreasing accuracy
targets = '2XFN 2N2A 2MOF 2MIC 2MFR 2M67 2LZL 2L91 2KSD 2KOG 2KNC 1P49 1FDM 1AFO 3W9T 1EK9 3KVN 3X2R 3CSL 3V8X 1UYN 2X55 4N75 4RDR 3FHH 4MEE 3RFZ 4B7O 3O44 4D5B 3DZM 1XKW 2GUF 2MMU 2MXB 2L2T 4QL0 3EGW 4GX0 5FOK 5IVA 2QTS 4RL8 5HK1 1KMO 4E1S 2KLU 4FQE 2J58 4AFK 5IXM 1FEP 2VDF 4I0U 2BS2 2WJR 4Y25 2WDQ 1KF6 1UUN 2FGQ 1PPJ 2KSE 3QRA 3CX5 4Q35 3SZV 4GEY 7AHL 2YNK 2MPR 3VY8 3DWO 3RLF 2BG9 4C69 2VPZ 3BS0 1OH2 5EKE 1QJP 2ZXE 4RLC 2X27 2X9K 4P02 3RQW 4RDQ 2LHF 2ZFG 4RJW 4HFI 3NE5 4KYT 2POR 3FID 5AJI 4WD8 2GR8 2QKS 4V1G 4UMW 4PL0 3AR4 3TUI 4DX5 4A82 4COF 4PHZ 3D31 4YMU 4CAD 1RZH 3WAJ 4P79 2L35 3GP6 4XES 4EIY 2NQ2 3VW7 4MSW 3AG3 1U19 4NV5 4ZR1 2LCK 4OGQ 5IRX 3ZOJ 2MPN 2K73 4TWK 5FXB 5DQQ 4A01 2KSF 2KS9 3KCU 3WFD 3TDS 4G7V 3WU2 4WW3 1QD6 4RI2 3S8G 3WXW 4O6M 1C17 2F2B 4XP9 4PJ0 4XNV 2JLN 2YEV 5EZM 1U7G 3B9W 3K3F 4X89 2WSW 4RYO 4ENE 4QUV 5A1S 4PD6 4XTL 4U4V 3GIA 4X5M 3RKO 5AWW 4R0C 4PGR 3UG9 2L0J 5IWS 4AL0 2A65 4RP9 4DVE 3ZE3 3PCV 5BZ3 4QTN 2N4X 4MES 4ZW9 3DH4 4Y28 4J05 4K1C 2BL2 4XK8 2AXT 2NWL 2CFQ 4IKV 2XOV 5HYA 4O6Y 4N7W 3M73 4KNF 3W4T 4UVM 4XU4 1H2S 5CKR 1M0K 4O93 5AYN 4KPP 4QND 4TQ3 1JB0 4UC1 4ZP0 4U9N 5I20 3LDC'.split()
targets = targets[:5] if debug else targets

#print(f'extension: {extension}')
#print(f'command_line: {command_line}')
#print(f'config: {benchmark.config()}')
#print(hpc_driver)

hpc_logs = f'{working_dir}/hpc-logs'
if not os.path.exists(hpc_logs): os.makedirs(hpc_logs)
hpc_job_ids = []
for target in targets:
    prefix = f'{working_dir}/output/'
    if not os.path.exists(prefix): os.makedirs(prefix)

    hpc_job_ids.append( hpc_driver.submit_hpc_job(
        name=f'{testname}-{target}',

        #==> EDIT HERE
        executable = f'{rosetta_dir}/source/bin/mp_lipid_acc.{extension}',
        arguments = command_line.format_map(vars()),
        working_dir = prefix,
        jobs_to_queue = min(nstruct, 50),
        log_dir = hpc_logs,
        time=24,
        block=False)
    )


# if not debug:
#     hpc_driver.wait_until_complete(hpc_job_ids, silent=True)
#     time.sleep(64)  # waiting for NFS caching
hpc_driver.wait_until_complete(hpc_job_ids, silent=True)

#==> EDIT HERE
benchmark.save_variables('targets nstruct working_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
