#!/usr/bin/env python
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

import shutil, os, tarfile, subprocess, sys

# this script is responsible for submitting the docking_local_refine
# scientific benchmark to a cluster



# targets.py just sets up a list of input targets. Keep them in a separate
# file so it's easy to re-run just a subset of the benchmark
from targets import targets


# m_vars contains the configuration information for this scientific
# benchmark, such as paths/compilers versions etc. Here we'll add
# additional configuration information needed to setup and run the
# benchmark
m_vars = eval(file('_arguments.py').read())
config_vars = eval(file('config.py').read())

m_vars = dict(m_vars.items() + config_vars.items())

def prepare_input_targets(targets, m_vars):
	""" The zdock set provides bound and unbound versions of the
	ligand and receptors.
	
	To obtain these targets, usually we download the full
	benchmark set from svn, extracts the zip file."""

	print "Preparing input structures..."
	print "\tChecking if the input partners exist..."
	all_input_partners_exist = True
	for target in targets:
		receptor_fname = m_vars["benchmark_data_path"] + "/"
		receptor_fname += m_vars["benchmark_data_set"] + "/"
		receptor_fname += target + m_vars["receptor_target_extension"]
		if not os.path.isfile(receptor_fname):
			if all_input_partners_exist:
				print "\tThe input receptor partner '%s' does not exist." % receptor_fname
				all_input_partners_exist = False

		ligand_fname = m_vars["benchmark_data_path"] + "/"
		ligand_fname += m_vars["benchmark_data_set"] + "/"
		ligand_fname += target + m_vars["ligand_target_extension"]
		if not os.path.isfile(receptor_fname):
			if all_input_partners_exist:
				print "\tThe input ligand partner '%s' does not exist." % ligand_fname
				all_input_partners_exist = False
	
	if all_input_partners_exist:
		print "\tAll input partners exist.\n"
		return
		
	print "\tChecking if the benchmark data set exists..."
	zipped_benchmark_data_set_fname = \
		m_vars["benchmark_data_path"] + "/" + m_vars["benchmark_data_set"] + ".tar.gz"
	if os.path.exists(zipped_benchmark_data_set_fname):
		extract_benchmark_set(m_vars)
		print "The benchmark dataset exists.\n"
		return

	print "\tThe benchmark dataset '%s' does not exist.\n" % (zipped_benchmark_data_set_fname)
	retrieve_benchmark_set(m_vars)
	extract_benchmark_set(m_vars)
	

def retrieve_benchmark_set(m_vars):
	print "\tRetrieving the benchmark set from svn:"
	if not os.path.exists(m_vars["benchmark_data_path"]):
		print "\tCreating benchmark data path '%s' because it does not exist." % m_vars["benchmark_data_path"]
		os.makedirs(m_vars["benchmark_data_path"])

	print "\t\tsvn checkout %s %s " % (m_vars["svn_benchmark_path"], m_vars["benchmark_data_path"])
	p = subprocess.Popen([
			"svn", "checkout", m_vars["svn_benchmark_path"],
			m_vars["benchmark_data_path"]], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
	p.wait()
	stdout, stderr = p.communicate()
	if stdout != '': print stdout
	if stderr != '': print stderr
	if p.returncode != 0:
		print "Unable to retrieve benchmark data set.  Please try to manually retrieve it with the command:"
		print "\n\t\tcd %s && svn checkout %s %s\n" % (
			m_vars["benchmark_data_path"],
			m_vars["svn_benchmark_path"],
			m_vars["benchmark_data_path"])
		print ""
		print "\tSee the README.txt for more details on obtaining the benchmark data set."
		exit(1)
	
def extract_benchmark_set(m_vars):
	print "\tExtracting the benchmark dataset '%s'..." % (
		m_vars["benchmark_data_path"] + "/" +
		m_vars["benchmark_data_set"] + ".tar.gz")
	p = subprocess.Popen([
			'tar', '-xzf', m_vars["benchmark_data_set"] + ".tar.gz"],
			cwd=m_vars["benchmark_data_path"], stderr=subprocess.PIPE)
	stderr = p.communicate()[1]
	if stderr != '': print stderr
	if p.returncode != 0:
		print "Unable to unzip the benchmark set."
		exit(1)

def prepare_output_paths(m_vars):
	print "preparing output paths..."
	run_log_path = os.path.abspath(os.path.join(m_vars["output_dir"], m_vars["output_run_log_path"]))
	if not os.path.exists(run_log_path):
		try:
			os.makedirs(run_log_path)
			print "\tOutput run log path: '%s'" % run_log_path
		except Exception, e:
			print "Unable to create the output run log path '%s':" % run_log_path
			print "\t", e
			exit(1)

	score_path = os.path.abspath(os.path.join(m_vars["output_dir"], m_vars["output_score_path"]))
	if not os.path.exists(score_path):
		try:
			os.makedirs(score_path)
			print "\tOutput score path: '%s'" % os.path.abspath(score_path)
		except Exception, e:
			print "Unable to create the output score path '%s':" % score_path
			print "\t", e
			exit(1)
	
	decoy_path = os.path.abspath(os.path.join(m_vars["output_dir"], m_vars["output_decoy_path"]))
	if not os.path.exists(decoy_path):	
		try:
			os.makedirs(decoy_path)
		except Exception, e:
			print "Unable to create the output decoy path '%s':" % decoy_path
			print "\t", e
			exit(1)


	test_results_log = os.path.abspath(os.path.join(m_vars["output_dir"], m_vars["test_results_log"]))
	if os.path.exists(test_results_log):
		os.remove(test_results_log)
	print "\tTest results log: '%s'" % test_results_log

	test_results_yaml = os.path.abspath(os.path.join(m_vars["output_dir"], m_vars["test_results_yaml"]))
	if os.path.exists(test_results_yaml):
		os.remove(test_results_yaml)
	print "\tTest results yaml: '%s'" % test_results_yaml


def prepare_flags_file(m_vars):
	flags_file = open("flags", 'w')
	flags_file.write('''#Autogenerated by rosetta_tests/cluster/docking_local_refine/submit.py
-in:path %(workdir)s/%(benchmark_data_path)s/%(benchmark_data_set)s
-ignore_unrecognized_res
-use_input_sc
-docking_local_refine
%(extra_chi_rotamers)s
%(energy_function_flags)s
-nstruct %(nstruct)s
-out:path:score %(output_dir)s/%(output_score_path)s
-out:path:all %(output_dir)s/%(output_decoy_path)s
-out:silent_gz
-out:file:silent_struct_type binary
-out:file:silent_preserve_H
''' % m_vars)
	flags_file.close()

def prepare_condor_header(m_vars):
	return '''#Autogenerated by rosetta_tests/cluster/docking_local_refine/submit.py
universe     = vanilla
Notify_user  =
notification = Error
Log          = .condorscript.log
Executable   = %(bin)s/docking_protocol.%(binext)s

Requirements = ( Memory > 512)
GetEnv       = True

''' % m_vars

def prepare_condor_target(target, m_vars):
	m_vars["target_fname"] = '"%s%s %s%s"' % \
			(target, m_vars["receptor_target_extension"],
			target, m_vars["ligand_target_extension"])
	m_vars["native_fname"] = m_vars["target_fname"]
	m_vars["scorefile_fname"] = target + m_vars["output_scorefile_extension"]
	m_vars["silentfile_fname"] = target + m_vars["output_silentfile_extension"]
	
	return '''
Error   = %(output_dir)s/%(output_run_log_path)s/%(target_fname)s.logerr
Output  = %(output_dir)s/%(output_run_log_path)s/%(target_fname)s.logout

arguments = -database %(database)s \\
	-native %(native_fname)s \\
	-s %(target_fname)s \\
	-multiple_processes_writing_to_one_directory  \\
	-out:file:scorefile %(output_dir)s/%(output_score_path)s/%(scorefile_fname)s \\
	-out:file:silent %(output_dir)s/%(output_decoy_path)s/%(silentfile_fname)s \\
        @flags

priority = %(condor_priority)s
queue %(condor_queue)s

''' % m_vars


def prepare_lsf_script(targets, m_vars):
	input_targets_file = open("%(output_dir)s/input.list" % m_vars, 'w')
	for target in targets:
		input_targets_file.write(
			'%s%s %s%s\n' % \
				(target, m_vars["receptor_target_extension"],
				target, m_vars["ligand_target_extension"]))
	input_targets_file.close()


	return '''#!/bin/sh
bsub \\
	-q %(lsf_queue_name)s \\
	-n %(num_cores)s \\
	-J docking_local_refine \\
	-o %(output_dir)s/%(output_run_log_path)s/docking_local_refine_%%J.log \\
	-e %(output_dir)s/%(output_run_log_path)s/docking_local_refine_%%J.err \\
	-a mvapich mpirun \\
	%(bin)s/docking_protocol.%(binext)s \\
	-database %(database)s \\
	-in:file:list %(output_dir)s/input.list \\
	-out:mpi_tracer_to_file %(output_dir)s/%(output_decoy_path)s/docking_local_refine.log \\
	-out:file:scorefile %(output_dir)s/%(output_score_path)s/docking_local_refine%(output_scorefile_extension)s \\
	-out:file:silent_print_all_score_headers \\
	-out:file:silent %(output_dir)s/%(output_decoy_path)s/docking_local_refine.silent \\
	@flags
''' % m_vars


if "run_type" not in m_vars or m_vars["run_type"] == "condor":
	print 'Running submit.py script for cluster docking_local_refine scientific test on a condor cluster...'
	prepare_input_targets(targets, m_vars)
	prepare_output_paths(m_vars)
	prepare_flags_file(m_vars)
	condor_script = prepare_condor_header(m_vars)
	
	for target in targets:
		condor_script += prepare_condor_target(target, m_vars)

	condor_script_fname = '%(output_dir)s/condor_script' % m_vars
	f = file(condor_script_fname, 'w')
	f.write(condor_script)
	f.close()

	p = subprocess.Popen(
		["condor_submit", condor_script_fname],
		stdout = subprocess.PIPE, stderr = subprocess.PIPE)
	p.wait()
	stdout, stderr = p.communicate()
	if stdout != '': print stdout
	if stderr != '': print stderr
	if p.returncode != 0:
        	print "Unable to submit condor script '%s'." % condor_script_fname
		exit(1)
	print "To resubmit this job:"
	print ""
	print "\tcd %(output_dir)s" % m_var
	print "condor_submit condor_script"

			
elif m_vars["run_type"] == "lsf":
	print 'Running submit.py script for cluster docking_local_refine scientific test on an lsf cluster...'
	prepare_input_targets(targets, m_vars)
	prepare_output_paths(m_vars)
	prepare_flags_file(m_vars)
		
	lsf_script = prepare_lsf_script(targets, m_vars)	
	f = file('lsf_script.sh', 'w');  f.write(lsf_script);  f.close()
	p = subprocess.Popen(["sh", "lsf_script.sh"],
		stdout = subprocess.PIPE, stderr = subprocess.PIPE)
	p.wait()
	stdout, stderr = p.communicate()
	if stdout != '': print stdout
	if stderr != '': print stderr
	if p.returncode != 0:
        	print "Unable to submit load sharing facility script 'lsf_script.sh':"
		exit(1)
	print "To resubmit this job:"
	print ""
	print "\tcd %(output_dir)s" % m_var
	print "\tsh lsf_script.sh" 

