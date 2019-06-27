#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  mp_f19_sequence_recovery/1.submit.py
## @brief this script is part of the franklin2019 sequence recovery test
## @author Sergey Lyskov, Rebecca F. Alford (ralford3@jhu.edu)

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

testname    = "mp_f19_sequence_recovery"

debug       = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

command_line = '''
-database {rosetta_dir}/database
-in:file:s {rosetta_dir}/tests/scientific/data/{testname}/{target}/{target}_tr_ignorechain.pdb
-mp:setup:spanfiles {rosetta_dir}/tests/scientific/data/{testname}/{target}/{target}_tr_ignorechain.span
-nstruct {nstruct}
-score:weights franklin2019
-in:membrane
-out:path:all {prefix}
-out:file:scorefile {prefix}/{target}.score
-no_color
'''.replace('\n', ' ').replace('  ', ' ')

nstruct = 1

targets = '1AFO 1EK9 1FDM 1FEP 1H2S 1KF6 1KMO 1KPK 1M0K 1P49 1QD6 1QJP 1RZH 1U19 1U7G 1UUN 1UYN 1XKW 2A65 2BL2 2BS2 2CFQ 2F2B 2FGQ 2GR8 2GUF 2JLN 2K73 2KIX 2KLU 2KNC 2KOG 2KS9 2KSD 2KSE 2KSF 2L0J 2L2T 2L35 2LCK 2LHF 2LZL 2M67 2MFR 2MIC 2MMU 2MOF 2MPN 2MPR 2MXB 2N2A 2N4X 2NQ2 2NWL 2POR 2QKS 2QTS 2R9R 2VDF 2VPZ 2WCD 2WDQ 2WJR 2WSW 2X27 2X55 2X9K 2XFN 2XOV 2YEV 2YNK 2ZFG 2ZXE 3AG3 3AR4 3B9W 3BS0 3CSL 3D31 3DH4 3DWO 3DZM 3EGW 3EMN 3FHH 3FID 3GIA 3GP6 3K3F 3KCU 3KVN 3LDC 3M73 3NE2 3NE5 3O44 3QRA 3RFZ 3RLF 3RQW 3SZV 3TUI 3UG9 3V8X 3VW7 3VY8 3W4T 3WAJ 3WFD 3WMF 3WXW 3X2R 3ZE3 3ZOJ 4A01 4A82 4AFK 4AL0 4B7O 4C69 4COF 4D5B 4DVE 4DX5 4E1S 4EIY 4ENE 4FQE 4G7V 4GEY 4GX0 4HFI 4I0U 4IKV 4J05 4K1C 4KNF 4KPP 4KYT 4MEE 4MES 4MSW 4N75 4N7W 4NV5 4O6M 4O6Y 4O93 4OGQ 4P02 4P79 4PD6 4PGR 4PHZ 4PL0 4Q35 4QL0 4QND 4QTN 4QUV 4R0C 4RDQ 4RDR 4RI2 4RJW 4RL8 4RLC 4RP9 4RYO 4TNW 4TQ3 4TWK 4U4V 4U9N 4UC1 4UMW 4UVM 4V1G 4WD8 4WW3 4X5M 4X89 4XES 4XNV 4XP9 4XTL 4XU4 4Y25 4ZP0 4ZR1 4ZW9 5A1S 5AJI 5AWW 5AYN 5BZ3 5CKR 5DQQ 5EKE 5EZM 5FOK 5HK1 5HYA 5I20 5IRX 5IVA 5IWS 5IXM 7AHL'.split()
targets = targets[:2] if debug else targets

hpc_logs = f'{working_dir}/hpc-logs'
if not os.path.exists(hpc_logs): os.makedirs(hpc_logs)
hpc_job_ids = []
for target in targets:
    prefix = f'{working_dir}/output/{target}'
    if not os.path.exists(prefix): os.makedirs(prefix)

    hpc_job_ids.append( hpc_driver.submit_hpc_job(
        name=f'{testname}-{target}',

        executable = f'{rosetta_dir}/source/bin/fixbb.{extension}',
        arguments = command_line.format_map(vars()),
        working_dir = prefix,
        jobs_to_queue = min(nstruct, 50),
        log_dir = hpc_logs,
        time=24,
        block=False)
    )

hpc_driver.wait_until_complete(hpc_job_ids, silent=True)

benchmark.save_variables('debug targets nstruct working_dir rosetta_dir extension testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
