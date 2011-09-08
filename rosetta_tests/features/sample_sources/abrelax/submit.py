#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

from shutil import rmtree
from os import path, makedirs, getcwd
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

    def setup_input_data(self, mvars):
        print "Checkout abrelax input dataset..."
        try:
            p = subprocess.Popen(
                args=['svn', 'checkout', "--non-interactive", "--no-auth-cache",
                      'https://svn.rosettacommons.org/source/trunk/mini.data/tests/scientific/cluster/features/sample_sources/abrelax/input/abrelax_lr5_input_structures_60_by_103_sequences_100908',
                      'input/'],
                stderr=subprocess.PIPE)
            err = p.communicate()[1]
            if err != '': print err
        except OSError, e:
            print >>sys.stderr, "Execution Failed:", e

        mvars["decoy_input_silent_structures"] = path.join(
            getcwd(), "input",
            "abrelax_lr5_input_structures_60_by_103_sequences_100908.silent")
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
        self.apply_mvars(mvars, "submit_lsf_job.sh.TEMPLATE", script_fname)
        p = subprocess.Popen(["sh", script_fname])
        p.wait()

    def submit_condor_job(self, mvars):
        script_fname = "%(output_dir)s/%(sample_source_id)s/condor_submit_script" % mvars
        print "Submit condor job %s..." % script_fname
        self.apply_mvars(mvars, "condor_submit_script.TEMPLATE", script_fname)
        p = subprocess.Popen(['condor_submit', script_fname])
        p.wait()
    
    def submit(self,
        mvars,
        submit_type,
        sample_source_id,
        sample_source_description):

        mvars['date_code'] = date.today().strftime("%y%m%d")
        mvars['sample_source_path'] = getcwd()
        mvars['sample_source_id'] = sample_source_id % mvars
        mvars['sample_source_description'] = sample_source_description % mvars

        self.setup_output_paths(mvars)
        mvars = self.setup_input_data(mvars)

        print "mvars:"
        for key, value in mvars.iteritems():
            print "\t%s\t'%s'" % (key.rjust(26), value)

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
      default="abrelax_r%(svn_revision)s_%(date_code)s",
      help="The output file is <output_dir>/<sample_source_id>/features_<sample_source_id>.db3. [default %default]")

    parser.add_option("--sample_source_description",
      default="Rosetta predictions using the ab initio relax protocol.",
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
      options.sample_source_description)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
