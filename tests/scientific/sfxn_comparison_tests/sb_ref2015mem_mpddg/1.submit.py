#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  mpddg/1.submit.py
## @brief this script is part of mpddg scientific test
## @author Julia Koehler Leman, Hope Woods, Johanna Tiemann

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

#==> EDIT HERE
testname	= "sb_ref2015mem_mpddg"

debug	   = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

#==> EDIT HERE
command_line_wt = '''
-database {rosetta_dir}/database
-in:file:s {rosetta_dir}/tests/scientific/data/mp_ddg/{target}/{target}_rm_wt_{mutant}.pdb
-mp:setup:spanfiles {rosetta_dir}/tests/scientific/data/mp_ddg/{target}/{target}_opm_A.span
-nstruct {nstruct}
-parser:protocol {rosetta_dir}/tests/scientific/data/mp_ddg/ddg_wt.xml
-parser:script_vars scorefxn=ref2015_memb
-mp:restore_lazaridis_imm_behavior 1
-packing:resfile {rosetta_dir}/tests/scientific/data/mp_ddg/{target}/{target}_mut_{mutant}.res

-out:file:silent {prefix}/{target}_wt_{mutant}.out 
-out:file:scorefile {prefix}/{target}_wt_{mutant}.sc 

-multiple_processes_writing_to_one_directory
-no_color
'''.replace('\n', ' ').replace('  ', ' ')

command_line_mut = '''
-database {rosetta_dir}/database
-in:file:s {rosetta_dir}/tests/scientific/data/mp_ddg/{target}/{target}_rm_mut_{mutant}.pdb
-mp:setup:spanfiles {rosetta_dir}/tests/scientific/data/mp_ddg/{target}/{target}_opm_A.span
-nstruct {nstruct}
-parser:protocol {rosetta_dir}/tests/scientific/data/mp_ddg/ddg_mut.xml
-parser:script_vars scorefxn=ref2015_memb
-mp:restore_lazaridis_imm_behavior 1
-packing:resfile {rosetta_dir}/tests/scientific/data/mp_ddg/{target}/{target}_mut_{mutant}.res

-out:file:silent {prefix}/{target}_mut_{mutant}.out 
-out:file:scorefile {prefix}/{target}_mut_{mutant}.sc 

-multiple_processes_writing_to_one_directory
-no_color
'''.replace('\n', ' ').replace('  ', ' ')


#==> EDIT HERE
nstruct = 1 if debug else 50

#==> EDIT HERE
targets = '1py6 2k74 2xov'.split()
targets = targets[1:2] if debug else targets

mutants = {
	"1py6" : "107_A 111_A 144_A 144_V 148_A 166_A 167_A 16_A 170_A 181_A 181_F 182_A 185_F 189_A 200_A 208_A 20_A 20_S 20_V 222_A 23_A 26_M 31_A 32_A 33_A 34_A 35_P 36_A 36_P 37_A 37_P 38_A 39_A 39_F 39_P 40_P 41_A 42_A 42_P 42_S 43_A 43_P 44_A 45_A 46_A 47_P 48_A 50_A 51_A 52_A 53_A 54_A 55_A 56_A 57_A 58_A 5_A 75_F 79_A 79_F 86_A 87_A 90_A 92_A 93_A 96_A 9_A".split(), 
	"2k74" : "148_A 152_G 157_G 19_G 26_L 29_G 40_A 57_G 62_G 82_A 83_A 89_A".split(),
	"2xov" : "100_A 104_V 107_A 108_V 109_A 109_V 10_A 110_A 111_A 112_A 113_A 115_A 116_G 117_A 119_V 120_F 124_A 125_V 128_A 129_A 131_A 136_A 139_A 146_G 14_A 14_V 151_A 153_A 157_A 163_V 167_V 171_V 178_A 179_V 21_A 22_A 30_A 44_A 46_A 47_A 48_F 49_A 50_A 51_A 53_A 54_A 55_A 56_A 57_A 5_A 60_A 63_A 64_A 65_A 68_F 70_F 72_V 76_A 78_A 79_A 7_A 80_A 81_A 83_A 84_A 85_A 87_A 88_A 89_A 90_A 94_A".split()
}

mutants = {"2k74" : "148_A 152_G 157_G".split()} if debug else mutants

#print(f'extension: {extension}')
#print(f'command_line: {command_line}')
#print(f'config: {benchmark.config()}')
#print(hpc_driver)

hpc_logs = f'{working_dir}/hpc-logs'
if not os.path.exists(hpc_logs): os.makedirs(hpc_logs)
hpc_job_ids = []
for target in targets:
	prefix = f'{working_dir}/output/{target}'
	if not os.path.exists(prefix): os.makedirs(prefix)
	
	for mutant in mutants[target]:
		
		hpc_job_ids.append( hpc_driver.submit_serial_hpc_job(
			name=f'{testname}-{target}',

			#==> EDIT HERE
			executable = f'{rosetta_dir}/source/bin/rosetta_scripts.{extension}',
			arguments = command_line_wt.format_map(vars()),
			working_dir = prefix,
			jobs_to_queue = min(nstruct, 50),
			log_dir = hpc_logs,
			time=24,
			block=False)
		)

		hpc_job_ids.append( hpc_driver.submit_serial_hpc_job(
			name=f'{testname}-{target}',

			#==> EDIT HERE
			executable = f'{rosetta_dir}/source/bin/rosetta_scripts.{extension}',
			arguments = command_line_mut.format_map(vars()),
			working_dir = prefix,
			jobs_to_queue = min(nstruct, 50),
			log_dir = hpc_logs,
			time=24,
			block=False)
		)

# if not debug:
#	 hpc_driver.wait_until_complete(hpc_job_ids, silent=True)
#	 time.sleep(64)  # waiting for NFS caching
hpc_driver.wait_until_complete(hpc_job_ids, silent=True)

#==> EDIT HERE
benchmark.save_variables('debug targets mutants nstruct working_dir rosetta_dir testname')  # Python black magic: save all listed variable to json file for next script use (save all variables if called without argument)
