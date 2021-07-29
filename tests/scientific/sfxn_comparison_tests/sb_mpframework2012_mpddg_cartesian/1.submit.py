#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

## @file  mpddg_cartesian/1.submit.py
## @brief this script is part of mpddg cartesian scientific test
## @author Johanna Tiemann, Hope Woods, Julia Koehler Leman

import os, sys, time
import benchmark

benchmark.load_variables()  # Python black magic: load all variables saved by previous script into globals
config = benchmark.config()

#==> EDIT HERE
testname	= "sb_mpframework2012_mpddg_cartesian"

debug	   = config['debug']
rosetta_dir = config['rosetta_dir']
working_dir = config['working_dir']
hpc_driver  = benchmark.hpc_driver()
extension   = benchmark.calculate_extension()

#==> EDIT HERE
nstruct = 1 if debug else 5

#==> EDIT HERE
command_line = '''
-database {rosetta_dir}/database 
-s {rosetta_dir}/tests/scientific/data/mp_ddg_cartesian/{target}/{target}_mpfcart.pdb 
-mp:setup:spanfiles {rosetta_dir}/tests/scientific/data/mp_ddg_cartesian/{target}/{target}_opm_A.span 
-mp:lipids:composition DMPC 
-score:weights {rosetta_dir}/tests/scientific/data/mp_ddg_cartesian/mpframework_2012_cart.wts 
-in:membrane 
-fa_max_dis 9 
-missing_density_to_jump 
-has_pore 0 
-ddg:mut_file {rosetta_dir}/tests/scientific/data/mp_ddg_cartesian/{target}/mutfile_{mutant_res} 
-ddg:legacy 0 
-ddg:optimize_proline 1 
-ddg:frag_nbrs 4 
-ddg:bbnbrs 1 
-ddg::dump_pdbs false 
-ddg:iterations {nstruct}
-ddg::cartesian 
-out:prefix {target}_{mutant_res}
'''.replace('\n', ' ').replace('  ', ' ')

#==> EDIT HERE
targets = '1py6 2k74 2xov'.split()
targets = targets[-1:] if debug else targets

mutants = {
	"1py6" : ["E5","L9","M16","T20","F23","K26","S31","D32","P33","D34","A35","K36","K37","F38","Y39","A40","I41","T42","T43","L44","V45","P46","A47","T48","F50","T51","M52","Y53","L54","S55","M56","L57","L58","M64","Y75","Y79","T86","P87","L90","D92","L93","L96","L107","D111","I144","L148","T166","F167","L170","Y181","P182","W185","S189","E200","D208","S222"],
	"2k74" : ["A19","E26","A29","P40","A57","A62","F82","R83","Y89","G148","A152","A157"],
	"2xov" : ["P5","T7","M10","C14","M21","Q22","M30","E44","W46","R47","Y48","F49","T50","H51","L53","M54","H55","F56","S57","H60","F63","N64","L65","W68","Y70","G72","E76","R78","L79","G80","Y81","K83","L84","I85","I87","T88","L89","I90","L94","Q100","G104","F107","G108","G109","L110","S111","G112","V113","Y115","A116","L117","G119","Y120","R124","G125","D128","P129","S131","Q136","L139","W146","W151","D153","M157","A163","H164","G167","G171","D178","S179"],
}

mutants = {"2xov" : ["N64","E76","G109","H164"]} if debug else mutants 

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
	
	for mutant_res in mutants[target]:
		
		hpc_job_ids.append( hpc_driver.submit_serial_hpc_job(
			name=f'{testname}-{target}',

			#==> EDIT HERE
			executable = f'{rosetta_dir}/source/bin/cartesian_ddg.{extension}',
			arguments = command_line.format_map(vars()),
			working_dir = prefix,
			jobs_to_queue = 1, #min(nstruct, 50),
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
