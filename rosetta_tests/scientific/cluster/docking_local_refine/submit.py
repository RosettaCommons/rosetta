#!/usr/bin/env python
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

import os, tarfile, subprocess

# This just sets up a list of input targets. Keep them in a separate
# file so it's easy to re-run just a subset of the benchmark
from targets import targets



# m_vars contains the configuration information for this scientific
# benchmark, such as paths/compilers versions etc. Here we'll add
# additional configuration information needed to setup and run the
# benchmark
m_vars = dict(
	eval(file('_arguments.py').read())
	eval(file('config.py').read()))


def prepare_input_targets(targets, m_vars):
    """ The zdock set provides bound and unbound versions of the
    ligand and receptors. For this scientific benchmark we need the
    partners concatenated together to form the input target.

    To obtain these targets, usually we download the full
    benchmark set from svn, extracts the zip file then cats the
    structures together."""

    print "Preparing input_structures..."
    
    print "\tChecking if all the input targets exist..."
    all_input_targets_exist = True
    for target in targets:
       target_fname = m_vars["input_target_path"] + "/"
       target_fname += target + m_vars["input_target_extension"]
       if not os.path.isfile(target_fname):
           print "\tThe input target '%s' does not exist." % target_fname
           all_input_targets_exist = False
           break
       
    if all_input_targets_exist:
        print "\tAll input targets exist.\n"
        return

    print "\tChecking if the input partners exist..."
    all_input_partners_exist = True
    receptor_fnames = []
    ligand_fnames = []
    for target in targets:
        receptor_fname = m_vars["benchmark_data_path"] + "/"
        receptor_fname += m_vars["benchmark_data_set"] + "/"
        receptor_fname += target + m_vars["receptor_target_extension"]
        receptor_fnames.append(receptor_fname)
        if not os.path.isfile(receptor_fname):
            if all_input_partners_exist:
                print "\tThe input receptor partner '%s' does not exist." % receptor_fname
                all_input_partners_exist = False

        ligand_fname = m_vars["benchmark_data_path"] + "/"
        ligand_fname += m_vars["benchmark_data_set"] + "/"
        ligand_fname += target + m_vars["ligand_target_extension"]
        ligand_fnames.append(ligand_fname)
        if not os.path.isfile(receptor_fname):
            if all_input_partners_exist:
                print "\tThe input ligand partner '%s' does not exist." % ligand_fname
                all_input_partners_exist = False
    
    if all_input_partners_exist:
        print "\tAll input partners exist.\n"
        generate_input_targets_from_partners(targets, receptor_fnames, ligand_fnames, m_vars)
        return
    
    print "\tChecking if the benchmark data set exists..."
    zipped_benchmark_data_set_fname = \
        m_vars["benchmark_data_path"] + "/" + m_vars["benchmark_data_set"] + ".tar.gz"
    if os.path.exists(zipped_benchmark_data_set_fname):
        extract_benchmark_set(m_vars)
        print "The benchmark dataset exists.\n"
        generate_input_targets_from_partners(targets, receptor_fnames, ligand_fnames, m_vars)
        return;

    print "\tThe benchmark dataset '%s' does not exist.\n" % (zipped_benchmark_data_set_fname)
    retrieve_benchmark_set(m_vars)
    extract_benchmark_set(m_vars)
    generate_input_targets_from_partners(targets, receptor_fnames, ligand_fnames, m_vars)        
    

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

def generate_input_targets_from_partners(
    targets, receptor_fnames, ligand_fnames, m_vars):
    print "\tGenerating input targets from input partners..."

    if len(receptor_fnames) != len(targets):
        print "WARNING: Requesting to generate '%s' targets from '%s' receptor tructures." % (len(targets), len(receptor_fnames))
    
    if len(ligand_fnames) != len(targets):
        print "WARNING: Requesting to generate '%s' targets from '%s' ligand tructures." % (len(targets), len(ligand_fnames))

    if not os.path.exists(m_vars["input_target_path"]):
        try:
            os.makedirs(m_vars["input_target_path"])
        except OSError, e:
            "Unable to make the input target path '%s':" % m_vars["input_target_path"]
            print e
            exit(1)
            

    for i in range(len(targets)):
        target = targets[i]
        target_fname = m_vars["input_target_path"] + "/" + target + m_vars["input_target_extension"]
        receptor_fname = receptor_fnames[i]
        ligand_fname = ligand_fnames[i]
        
        target_contents = ""
        try:
            receptor_file = open(receptor_fname, 'r')
            target_contents += receptor_file.read()
            receptor_file.close()
        except OSError, e:
            print "Unable to read the receptor file '%s':" % receptor_fname
            print e
            exit(1)

        try:
            ligand_file = open(ligand_fname, 'r')
            target_contents += ligand_file.read()
            ligand_file.close()
        except OSError, e:
            print "Unable to read the ligand file '%s':" % ligand_fname
            print e
            exit(1)

        try:
            target_file = open(target_fname, 'w')
            target_file.write(target_contents)
            target_file.close()
        except OSError, e:
            print "Unable to write the target file '%s':" % target_fname
            print e        
            exit(1)

def prepare_output_paths(m_vars):
    print "preparing output paths..."
    if os.path.exists(m_vars["output_run_log_path"]):
        try:
            os.rmtree(m_vars["output_run_log_path"])
            os.makedirs(m_vars["output_run_log_path"])
            print "\tOutput run log path: '%s'" % os.path.abspath(m_vars["output_run_log_path"])
        except e:
            print "Unable to create the output run log path '%s':" % m_vars["output_run_log_path"]
            print e
            exit(1)

    if os.path.exists(m_vars["output_score_path"]):
        try:
            os.rmtree(m_vars["output_score_path"])
            os.makedirs(m_vars["output_score_path"])
            print "\tOutput score path: '%s'" % os.path.abspath(m_vars["output_score_path"])
        except e:
            print "Unable to create the output score path '%s':" % m_vars["output_score_path"]
            print e
            exit(1)

    if os.path.exists(m_vars["output_decoy_path"]):    
        try:
            os.rmtree(m_vars["output_decoy_path"])
            os.makedirs(m_vars["output_decoy_path"])
        except e:
            print "Unable to create the output decoy path '%s':" % m_vars["output_decoy_path"]
            print e
            exit(1)

    if os.path.exists(m_vars["test_results_log"]):
        os.remove(m_vars["test_results_log"])
    print "\tTest results log: '%s'" % os.path.abspath(m_vars["test_results_log"])

    if os.path.exists(m_vars["test_results_yaml"]):
        os.remove(m_vars["test_results_yaml"])
    print "\tTest results yaml: '%s'" % os.path.abspath(m_vars["test_results_yaml"])


def prepare_condor_header(m_vars):
    return '''
universe     = vanilla
Notify_user  =
notification = Error
Log          = .condorscript.log
Executable   = %(bin)s/docking_protocol.%(binext)s

Requirements = ( Memory > 512)
GetEnv       = True

''' % m_vars

def prepare_condor_target(target, m_vars):
    m_vars["target_fname"] = target + m_vars["input_target_extension"]
    m_vars["native_fname"] = m_vars["target_fname"]
    m_vars["scorefile_fname"] = target + m_vars["output_scorefile_extension"]
    m_vars["silentfile_fname"] = target + m_vars["output_silentfile_extension"]
    
    return '''
Error   = %(output_run_log_path)s/%(target_fname)s.logerr
Output  = %(output_run_log_path)s/%(target_fname)s.logout

arguments = -database %(database)s \\
	-in:path %(input_target_path)s \\
	-native %(native_fname)s \\
	-s %(target_fname)s \\
	-use_input_sc \\
	-docking_local_refine \\
	%(extra_chi_rotamers)s \\
	-nstruct %(nstruct)s \\
	-multiple_processes_writing_to_one_directory  \\
	-out:path:score %(output_score_path)s \\
	-out:file:scorefile %(scorefile_fname)s \\
	-out:path:all %(output_decoy_path)s \\
	-out:silent %(silentfile_fname)s \\
	-out:silent_gz \\
	-out:file:silent_struct_type binary \\
	-out:file:silent_preserve_H

priority = %(condor_priority)s
queue %(condor_queue)s

''' % m_vars




print 'Running submit.py script for cluster docking_local_refine scientific test...'
prepare_input_targets(targets, m_vars)
prepare_output_paths(m_vars)

condor_script = prepare_condor_header(m_vars)

for target in targets:
    condor_script += prepare_condor_target(target, m_vars)

f = file('condor', 'w');  f.write(condor_script);  f.close()
