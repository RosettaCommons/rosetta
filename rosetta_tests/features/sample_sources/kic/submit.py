#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

from shutil import rmtree
from os import path, makedirs, getcwd
from glob import glob
import sys, re
import subprocess
from datetime import date
from optparse import OptionParser

class SampleSource:

    def setup_output_paths(self, mvars):
        output_path = "%(output_dir)s/%(sample_source_id)s" % mvars
        if path.exists(output_path):
            rmtree(output_path)
        try:
            makedirs(output_path)
        except OSError, e:
            print >>sys.stderr, "Cannot make directory:", e

        try:
            makedirs(output_path+"/decoys/plop_set")
        except OSError, e:
            print >>sys.stderr, "Cannot make directory:", e

        try:
            makedirs(output_path+"/decoys/rosetta_set")
        except OSError, e:
            print >>sys.stderr, "Cannot make directory:", e            

    def setup_input_data(self, mvars):
        print "Checkout KIC input dataset..."
        try:
            p = subprocess.Popen(
                args=['svn', 'checkout', "--non-interactive", "--no-auth-cache",
                      'https://svn.rosettacommons.org/source/trunk/mini.data/tests/scientific/cluster/features/sample_sources/kic/input/',
                      'input/'],
                stderr=subprocess.PIPE)
            err = p.communicate()[1]
            if err != '': print err
        except OSError, e:
            print >>sys.stderr, "Execution Failed:", e

        p = subprocess.Popen(args=['tar', '-xzf', 'plop_set.tar.gz'], cwd='input/')
        p.wait()
        
        p = subprocess.Popen(args=['tar', '-xzf', 'rosetta_set.tar.gz'], cwd='input/')
        p.wait()

        mvars["targets"] = [getcwd() + "/" + fname for fname in glob("input/plop_set/start/*_min.pdb") + \
                glob("input/rosetta_set/start/*_min.pdb") ]

        return mvars

    def apply_mvars(self, mvars, in_fname, out_fname):
        contents = file(in_fname).read()
        contents = contents % mvars
        f = file(out_fname, 'w')
        f.write(contents)
        f.close()
    
    def submit_local_job(self, mvars):
        script_fname = "%(output_dir)s/%(sample_source_id)s/submit_local_job.sh" % mvars
        print "Submit local job %s..." % script_fname

        self.apply_mvars(mvars, "submit_local_job.sh.TEMPLATE", script_fname)
        p = subprocess.Popen(["sh", script_fname])
        p.wait()

    def submit_lsf_job(self, mvars):
        script_fname = "%(output_dir)s/%(sample_source_id)s/submit_lsf_job.sh" % mvars
        print "Submit lsf job", script_fname, "to the %(lsf_queue_name)s queue requesting %(num_cores)s cores..." % mvars

        target_queue_stmts_TEMPLATE = """
bsub \\
	-q %(lsf_queue_name)s \\
	-n %(num_cores)s \\
	-J decoys_%(data_set_tag)s_%(target_tag)s \\
	-o %(output_log_fname)s \\
	-e %(error_log_fname)s \\
	-a mvapich mpirun \\
	%(bin)s/loopmodel.%(binext)s \\
	-database %(database)s \\
        @%(decoys_flags)s \\
	-in:file:native %(start_fname)s \\
	-loops:input_pdb %(start_fname)s \\
	-loops:loop_file %(loop_fname)s \\
	-out:file:silent %(silent_fname)s \\
	-out:mpi_tracer_to_file %(tracer_output)s
"""
        target_queue_stmts = ""
        for start_fname in mvars["targets"]:
            data_set_tag = start_fname.split("/")[-3]
            target_tag = start_fname.split("/")[-1][:4]
            target_vars = dict(
                lsf_queue_name = mvars["lsf_queue_name"],
                num_cores = mvars["num_decoy_cores"], 
                bin = mvars["bin"],
                binext = mvars["binext"],
                data_set_tag = data_set_tag,
                target_tag = target_tag,
                error_log_fname = "%(output_dir)s/%(sample_source_id)s/decoys/" + data_set_tag + "/" + target_tag + "_%%J.error.log",
                output_log_fname = "%(output_dir)s/%(sample_source_id)s/decoys/" + data_set_tag + "/" + target_tag + "_%%J.output.log",
                database = mvars["database"],
                decoys_flags = "%(output_dir)s/%(sample_source_id)s/decoys.flags",
                start_fname = start_fname,
                loop_fname = getcwd() + "/input/"+data_set_tag+"/loops/"+target_tag+".loop",
                silent_fname = "%(output_dir)s/%(sample_source_id)s/decoys/" + data_set_tag + "/" + target_tag + ".silent",
                tracer_output = "%(output_dir)s/%(sample_source_id)s/decoys/" + data_set_tag + "/" + target_tag + ".tracer")
            target_queue_stmts += target_queue_stmts_TEMPLATE % target_vars
        target_queue_stmts = target_queue_stmts % mvars

        features_queue_stmt = """
bsub \\
	-q %(lsf_queue_name)s \\
	-n %(num_cores)s \\
	-J features_%(sample_source_id)s \\
	-o %(output_dir)s/%(sample_source_id)s/features_%(sample_source_id)s_%%J.log \\
	-e %(output_dir)s/%(sample_source_id)s/features_%(sample_source_id)s_%%J.err \\
	-w \"done(decoys_plop_*) && done(decoys_rosetta_*)\" \\
	-a mvapich mpirun \\
	%(bin)s/rosetta_scripts.%(binext)s \\
	-database %(database)s \\
	-out:mpi_tracer_to_file %(output_dir)s/%(sample_source_id)s/features_%(sample_source_id)s.tracer \\
	@%(output_dir)s/%(sample_source_id)s/features.flags
"""  % mvars    
        
        self.apply_mvars(
            dict(mvars.items(),
                 target_queue_stmts=target_queue_stmts,
                 features_queue_stmt=features_queue_stmt),
            "submit_lsf_job.sh.TEMPLATE", script_fname)            

        p = subprocess.Popen(['sh', script_fname])
        p.wait()


    def submit_condor_job(self, mvars):
        script_fname = "%(output_dir)s/%(sample_source_id)s/condor_submit_script" % mvars
        print "Submit condor job %s..." % script_fname

        target_queue_stmts_TEMPLATE = """
Error = %(error_log_fname)s 
Output = %(output_log_fname)s
Executable = %(bin)s/loopmodel.%(binext)s
arguments = -database %(database)s \\
	@%(decoys_flags)s \\
	-in:file:native %(start_fname)s \\
	-loops:input_pdb %(start_fname)s \\
	-loops:loop_file %(loop_fname)s \\
	-out:file:silent %(silent_fname)s \\
	-seed_offset $(Process)
queue %(num_cores)s
"""
        target_queue_stmts = ""
        for start_fname in mvars["targets"]:
            data_set_tag = start_fname.split("/")[-3]
            target_tag = start_fname.split("/")[-1][:4]
            target_vars = dict(
                error_log_fname = "%(output_dir)s/%(sample_source_id)s/decoys/" + data_set_tag + "/" + target_tag + ".error.log",
                output_log_fname = "%(output_dir)s/%(sample_source_id)s/decoys/" + data_set_tag + "/" + target_tag + ".output.log",
                bin = mvars["bin"],
                binext = mvars["binext"],
                database = mvars["database"],
                decoys_flags = "%(output_dir)s/%(sample_source_id)s/decoys.flags",
                start_fname = start_fname,
                loop_fname = getcwd() + "/input/"+data_set_tag+"/loops/"+target_tag+".loop",
                silent_fname = "%(output_dir)s/%(sample_source_id)s/decoys/" + data_set_tag + "/" + target_tag + ".silent",
                num_cores = mvars["num_decoy_cores"])

            target_queue_stmts += target_queue_stmts_TEMPLATE % target_vars

        target_queue_stmts = target_queue_stmts % mvars

        features_queue_stmt = """
Error = %(output_dir)s/%(sample_source_id)s/features_%(sample_source_id)s.error.log
Output = %(output_dir)s/%(sample_source_id)s/features_%(sample_source_id)s.output.log
Executable = %(bin)s/rosetta_scripts.%(binext)s
arguments = -database %(database)s \\
	-out:mpi_tracer_to_file %(output_dir)s/%(sample_source_id)s/features_%(sample_source_id)s.log \\
	@%(output_dir)s/%(sample_source_id)s/features.flags \\

queue %(num_cores)s
""" % mvars
            
        self.apply_mvars(
            dict(target_queue_stmts=target_queue_stmts,
                 features_queue_stmt=features_queue_stmt),
            "condor_submit_script.TEMPLATE", script_fname)
        
        p = subprocess.Popen(['condor_submit', script_fname])
        p.wait()
    
    def submit(self,
        mvars,
        submit_type,
        sample_source_id,
        sample_source_description,
        nstruct):

        mvars['date_code'] = date.today().strftime("%y%m%d")
        mvars['sample_source_path'] = getcwd()
        mvars['sample_source_id'] = sample_source_id % mvars
        mvars['sample_source_description'] = sample_source_description % mvars
        mvars['nstruct'] = nstruct
        mvars['num_decoy_cores'] = min(nstruct + 1, int(mvars["num_cores"]))
        mvars['num_features_cores'] = mvars["num_cores"]

        self.setup_output_paths(mvars)
        mvars = self.setup_input_data(mvars)

        print "mvars:"
        for key, value in mvars.iteritems():
            print "\t%s\t'%s'" % (key.rjust(26), value)

        self.apply_mvars(mvars, "decoys.flags.TEMPLATE", "%(output_dir)s/%(sample_source_id)s/decoys.flags" % mvars)

        decoy_silent_fnames_str = ""
        for start_fname in mvars["targets"]:
            data_set_tag = start_fname.split("/")[-3]
            target_tag = start_fname.split("/")[-1][:4]
            silent_fname_base = "%(output_dir)s/%(sample_source_id)s/decoys/" + data_set_tag + "/" + target_tag % mvars
            for i in range(1, mvars["num_decoy_cores"]):
                decoy_silent_fnames_str += " \\\n\t%s_%s.silent" % (silent_fname_base, i) 
        mvars["decoy_silent_fnames"] = decoy_silent_fnames_str % mvars

        self.apply_mvars(mvars, "features.flags.TEMPLATE", "%(output_dir)s/%(sample_source_id)s/features.flags" % mvars)
        self.apply_mvars(mvars, "features.xml.TEMPLATE", "%(output_dir)s/%(sample_source_id)s/features.xml" % mvars)

        if submit_type == "local":
            self.submit_local_job(mvars)
        elif submit_type == "lsf":
            self.submit_lsf_job(mvars)
        else:
            self.submit_condor_job(mvars)

def main(argv):

    ## Parse and validate command line arguments, then call SampleSource.submit()

    parser = OptionParser(usage="Usage: %prog [OPTIONS]")
    parser.add_option("--sample_source_id",
      default="kic_r%(svn_revision)s_%(date_code)s",
      help="The output file is <output_dir>/<sample_source_id>/<sample_source_id>.db3. [default %default]")

    parser.add_option("--sample_source_description",
      default="Rosetta predictions using the kinimatic loop closure protocol.",
      help="This is included into the features database, 'SELECT description from sample_sources;'")

    parser.add_option("--local",
      default=False, action="store_true",
      help="Run this script locally rather submitting the job to run on a cluster.  This can be useful for debugging or one-off jobs. [default %default]")

    parser.add_option("--lsf",
      default=False, action="store_true",
      help="Generate features on an lsf MPI based cluster. [default %default]")

    parser.add_option("--arguments",
      default="../../_arguments.py",
      help="A file that contains a python dictionary of variables initialized by the cluster.py scientific benchmark driver script usually generated by the test/scientific/cluster/cluster.py driver script and deposited into test/scientific/cluster/features/ . [default %default]")

    parser.add_option("--nstruct",
      default=200,
      help="How many predictions should be made for each input structure.")

    (options, args) = parser.parse_args(args=argv)

    mvars = None
    try:
        f = open(options.arguments)        
        mvars = eval(f.read())
        f.close()
    except:
        print "ERROR: The arguments file, %s, could no be loaded. This is usually the _arguments.py file generated by the cluster.py scientific benchmark script." % options.arguments


    submit_type = None
    if options.local:
        submit_type = "local"
    elif options.lsf:
        submit_type = "lsf"
    else:
        submit_type = "condor"

    SampleSource().submit(
      mvars,
      submit_type,
      options.sample_source_id,
      options.sample_source_description,
      options.nstruct)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
