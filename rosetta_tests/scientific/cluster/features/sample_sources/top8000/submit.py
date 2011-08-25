#!/usr/bin/env python
# :noTabs=true:
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

from shutil import rmtree, copy
from os import path, makedirs, getcwd
import sys, subprocess, glob, re
from datetime import date
from optparse import OptionParser

class SampleSource:

    def setup_input_data(self):
        print "Checkout top8000 dataset..."
#        try:
#	    p = subprocess.Popen(
#	        args=['rsync', '-az'
#	              'https://garin.med.unc.edu/~momeara/top8000/input',
#	              'input'],
#	        stderr=subprocess.PIPE)
#	    err = p.communicate()[1]
#	    if err != '': print err
#        except OSError, e:
#            print >>sys.stderr, "Execution Failed:", e
           

        print "unzipping top8000.tar.gz..."
        try:
            p = subprocess.Popen(['tar', '-xzf', 'top8000.tar.gz'],
                  cwd='input', stderr=subprocess.PIPE)
            err = p.communicate()[1]
            if err != '': print err
        except OSError, e:
            print >>sys.stderr, "Execution Failed:", e
   
        print "generating all_pdbs.list..."
        try:
            f = open("input/all_pdbs.list", 'w')
            p = subprocess.Popen(["find", getcwd() + "/input/top8000_chains_eds_70", "-name", "*pdb"],
                  stdout=f, stderr=subprocess.PIPE)
            err = p.communicate()[1]
            if err != '': print err
            f.close()
        except OSError, e:
            print >>sys.stderr, "Execution Failed:", e

        print "Done setting up input files"

    def setup_output_paths(self, mvars):
        output_path = "%(output_dir)s/%(sample_source_id)s" % mvars
        if path.exists(output_path):
            rmtree(output_path)
        try:
            makedirs(output_path)
        except OSError, e:
            print >>sys.stderr, "Cannot make directory:", e

    def setup_svn_info(self, mvars):
    	try:
    	    p = subprocess.Popen(
    	        ['svn', 'info'],
    	        stdout=subprocess.PIPE,
    	        cwd=mvars["minidir"])
    	    p.wait()
    	    svn_info = p.communicate()[0]
    	    mvars["svn_url"] = re.search("URL: (\S*)", svn_info).group(1)
    	    mvars["svn_revision"] = int(re.search("Revision: (\d*)", svn_info).group(1))
    	except: 
    	    print "WARNING: Unable to get svn info"
    	    mvars["svn_url"]="UNKNOWN"
    	    mvars["svn_revision"]=-1
        return mvars

    def apply_mvars(self, mvars, in_fname, out_fname):
        contents = file(in_fname).read()
        contents = contents % mvars
        f = file(out_fname, 'w')
        f.write(contents)
        f.close()
    
    def setup_job(
        self,
        mvars,
        sample_source_id,
        sample_source_description):

        mvars['date_code'] = date.today().strftime("%y%m%d")
        mvars['sample_source_path'] = getcwd()
        mvars['sample_source_id'] = sample_source_id % mvars
        mvars['sample_source_description'] = sample_source_description % mvars

        self.setup_input_data()
        self.setup_output_paths(mvars)
        mvars = self.setup_svn_info(mvars)

        print "mvars:"
        for key, value in mvars.iteritems():
            print "\t%s\t'%s'" % (key.rjust(26), value)

        copy("input/rosetta_inputs.db3", "%(output_dir)s/%(sample_source_id)s/rosetta_inputs.db3" % mvars)
        self.apply_mvars(mvars, "features.xml.TEMPLATE", "%(output_dir)s/%(sample_source_id)s/features.xml" % mvars)
        self.apply_mvars(mvars, "flags.TEMPLATE", "%(output_dir)s/%(sample_source_id)s/flags" % mvars)
        self.apply_mvars(mvars, "submit_local_job.sh.TEMPLATE", "%(output_dir)s/%(sample_source_id)s/submit_local_job.sh" % mvars)
        self.apply_mvars(mvars, "condor_submit_script.TEMPLATE", "%(output_dir)s/%(sample_source_id)s/condor_submit_script" % mvars)

    def submit_job(self, submit_type, mvars):
        if submit_type == "local":
            print "Submit local job..."
            p = subprocess.Popen(["sh", "%(output_dir)s/%(sample_source_id)s/submit_local_job.sh" % mvars])
            p.wait()
        elif submit_type == "lsf":
            print "Submit lsf job to the %(lsf_queue_name)s queue requesting %(num_cores)s cores..." % mvars 
            self.apply_mvars(mvars, "submit_lsf_job.sh.TEMPLATE", "%(output_dir)s/%(sample_source_id)s/submit_lsf_job.sh" % mvars)
            p = subprocess.Popen(["sh", "%(output_dir)s/%(sample_source_id)s/submit_lsf_job.sh" % mvars])
            p.wait()
        else:
            print "Submit condor job..."
            p = subprocess.Popen(['condor_submit', '%(output_dir)s/%(sample_source_id)s/condor_submit_script' % mvars])
            p.wait()

    def submit(
        self,
        mvars,
        submit_type,
        sample_source_id,
        sample_source_description,
        setup_only,
        execute_only):

        if(not execute_only):
            self.setup_job(
                mvars,
                sample_source_id,
                sample_source_description)
            
        if(not setup_only):
            self.submit_job(submit_type, mvars)

def main(argv):

    parser = OptionParser(usage="Usage: %prog [OPTIONS]")

    parser.add_option("--sample_source_id",
      default="top8000_r%(svn_revision)s_%(date_code)s",
      help="The output file is <sample_source_id>.db3, as it is substituted into features.xml.TEMPLATE. [default %default]")

    parser.add_option("--sample_source_description",
      default="Top 8000 chains in PDB from the Richardson Lab",
      help="This is included into the features database, 'SELECT description from sample_sources;'. [default %default]")

    parser.add_option("--local",
      default=False, action="store_true",
      help="Run this script locally rather submitting the job to run on a cluster.  This can be useful for debugging or one-off jobs. [default %default]")

    parser.add_option("--lsf",
      default=False, action="store_true",
      help="Generate features on an lsf MPI based cluster. [default %default]")

    parser.add_option("--arguments",
      default="../../_arguments.py",
      help="A file that contains a python dictionary of variables initialized by the cluster.py scientific benchmark driver script. [default %default]")

    parser.add_option("--setup_only",
      action="store_true", default=False,
      help="Only setup the submit jobs but do not execute them. [default %default]")

    parser.add_option("--submit_only",
      action="store_true", default=False,
      help="Assume the runs have been setup and submit them directly.  Note, this option can be used with the --local option. [default %default][5~")

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
      options.setup_only,
      options.submit_only)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
