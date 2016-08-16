#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

import os, sys
import subprocess
import glob, shutil

sys.path.insert(0, "/".join(os.getcwd().split("/")[:-1]))

from submit import BaseSampleSource

class SampleSource(BaseSampleSource):


    def __init__(self, argv):
        BaseSampleSource.__init__(self)

        self.sample_source_description_default_value = \
            "Rosetta predictions using the ab initio relax protocol"

        self.initialize_options_parser()

        self.parser.add_option("--nstruct",
            default=200,
            help="How many predictions should be made for each input structure.")

        options = self.parse_options(argv)

        self.mvars['nstruct'] = options.nstruct
        self.mvars['num_decoy_cores'] = \
            min(options.nstruct + 1, int(self.mvars["num_cores"]))
        self.mvars['num_features_cores'] = self.mvars["num_cores"]

        self.setup()




    def setup_input_data(self):
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

        self.mvars["targets"] = [os.getcwd() + "/" + fname for fname in glob.glob("input/plop_set/start/*_min.pdb") + \
                glob.glob("input/rosetta_set/start/*_min.pdb") ]

        decoy_silent_fnames_str = ""
        for start_fname in self.mvars["targets"]:
            data_set_tag = start_fname.split("/")[-3]
            target_tag = start_fname.split("/")[-1][:4]
            silent_fname_base = "%(output_dir)s/%(sample_source_id)s/decoys/" + data_set_tag + "/" + target_tag % self.mvars
            for i in range(1, self.mvars["num_decoy_cores"]):
                decoy_silent_fnames_str += " \\\n\t%s_%s.silent" % (silent_fname_base, i) 
        self.mvars["decoy_silent_fnames"] = decoy_silent_fnames_str % self.mvars

        self.setup_condor_script()
        self.setup_lsf_script()
    
    def setup_lsf_script(self):

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
        for start_fname in self.mvars["targets"]:
            data_set_tag = start_fname.split("/")[-3]
            target_tag = start_fname.split("/")[-1][:4]
            target_vars = dict(
                lsf_queue_name = self.mvars["lsf_queue_name"],
                num_cores = self.mvars["num_decoy_cores"], 
                bin = self.mvars["bin"],
                binext = self.mvars["binext"],
                data_set_tag = data_set_tag,
                target_tag = target_tag,
                error_log_fname = "%(output_dir)s/%(sample_source_id)s/decoys/" + data_set_tag + "/" + target_tag + "_%%J.error.log",
                output_log_fname = "%(output_dir)s/%(sample_source_id)s/decoys/" + data_set_tag + "/" + target_tag + "_%%J.output.log",
                database = self.mvars["database"],
                decoys_flags = "%(output_dir)s/%(sample_source_id)s/decoys.flags",
                start_fname = start_fname,
                loop_fname = os.getcwd() + "/input/"+data_set_tag+"/loops/"+target_tag+".loop",
                silent_fname = "%(output_dir)s/%(sample_source_id)s/decoys/" + data_set_tag + "/" + target_tag + ".silent",
                tracer_output = "%(output_dir)s/%(sample_source_id)s/decoys/" + data_set_tag + "/" + target_tag + ".tracer")
            target_queue_stmts += target_queue_stmts_TEMPLATE % target_vars

        self.mvars["target_queue_stmts"] = target_queue_stmts % self.mvars

        self.mvars["features_queue_stmt"] = """
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
"""  % self.mvars    


    def setup_condor_script(self):
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
        for start_fname in self.mvars["targets"]:
            data_set_tag = start_fname.split("/")[-3]
            target_tag = start_fname.split("/")[-1][:4]
            target_vars = dict(
                error_log_fname = "%(output_dir)s/%(sample_source_id)s/decoys/" + data_set_tag + "/" + target_tag + ".error.log",
                output_log_fname = "%(output_dir)s/%(sample_source_id)s/decoys/" + data_set_tag + "/" + target_tag + ".output.log",
                bin = self.mvars["bin"],
                binext = self.mvars["binext"],
                database = self.mvars["database"],
                decoys_flags = "%(output_dir)s/%(sample_source_id)s/decoys.flags",
                start_fname = start_fname,
                loop_fname = os.getcwd() + "/input/"+data_set_tag+"/loops/"+target_tag+".loop",
                silent_fname = "%(output_dir)s/%(sample_source_id)s/decoys/" + data_set_tag + "/" + target_tag + ".silent",
                num_cores = self.mvars["num_decoy_cores"])

            target_queue_stmts += target_queue_stmts_TEMPLATE % target_vars

        self.mvars["target_queue_stmts"] = target_queue_stmts % self.mvars

        self.mvars["features_queue_stmt"] = """
Error = %(output_dir)s/%(sample_source_id)s/features_%(sample_source_id)s.error.log
Output = %(output_dir)s/%(sample_source_id)s/features_%(sample_source_id)s.output.log
Executable = %(bin)s/rosetta_scripts.%(binext)s
arguments = -database %(database)s \\
	-out:mpi_tracer_to_file %(output_dir)s/%(sample_source_id)s/features_%(sample_source_id)s.log \\
	@%(output_dir)s/%(sample_source_id)s/features.flags \\

queue %(num_cores)s
""" % self.mvars



    def setup_output_paths(self):
        output_path = "%(output_dir)s/%(sample_source_id)s" % self.mvars
        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        try:
            os.makedirs(output_path)
        except OSError, e:
            print >>sys.stderr, "Cannot make directory:", e

        try:
            os.makedirs(output_path+"/decoys/plop_set")
        except OSError, e:
            print >>sys.stderr, "Cannot make directory:", e

        try:
            os.makedirs(output_path+"/decoys/rosetta_set")
        except OSError, e:
            print >>sys.stderr, "Cannot make directory:", e            


def main(argv):

    ss = SampleSource(argv)
    ss.submit()


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
