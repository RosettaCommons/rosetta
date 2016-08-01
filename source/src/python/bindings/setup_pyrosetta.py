#!/usr/bin/env python
# Author: Jared Adolf-Bryfogle (jadolfbr@gmail.com)

# Purpose:
#  This file adds a file to your home directory, .pyrosetta_env, "
#  and by default adds a line to your bash or zsh profile (if not already present) to source this file.
#
#  After running this script, you will no longer need to run anything before using PyRosetta or add anything to your profile.

import os
import sys
import time
import re
from collections import defaultdict

from argparse import ArgumentParser





def get_profile_fname():
    """
    Gets the filename of the profile, which is different depending on operating system and shell.
    PLEASE edit this if you have a way to make this more robust.  It feels like the wild west getting this info!
    :rtype: str
    """
    home = os.environ['HOME']
    possible_names = defaultdict(list)
    possible_names['bash'] = [".bashrc", ".bash_profile", ".profile"]
    possible_names['zsh'] = [".zshrc", ".profile"]

    shell_name = os.path.basename(os.environ['SHELL'])

    #Use a common name first with the shell + rc in case of funky shells.
    common_name = shell_name+"rc"
    common_path = os.path.join(home, common_name)
    if os.path.exists(common_path):
        return common_path

    #Go through our common list
    for fname in possible_names[shell_name]:
        possible_path = os.path.join(os.environ['HOME'], fname)
        if os.path.exists(possible_path):
            print "Profile found: "+possible_path
            return possible_path

    #Exit and warn user to use the command-line option
    sys.exit("Could not find a profile path using common paths.  Please set one using the commandline using the option profile_path")


def ask_append_rc():
    append_rc = True

    v = raw_input("Append 'source "+outfile_path+"' to "+rc+"?   y / n; Default y  ")
    #print repr(v)
    if v:
        v = v.upper()
        if v == "Y" or v == "YES":
            append_rc = True
        elif v == "N" or v == "NO":
            append_rc = False
        else:
            sys.exit("Unrecognized input: "+v)

    return append_rc




if __name__ == "__main__":

    parser = ArgumentParser("This file adds a file to your home directory, "
                            ".pyrosetta_env, "
                            "and by default adds a line to your bash or zsh profile (if not already present) to source this file. "
                            "  After running this script, you will no longer need to run anything before using PyRosetta or add anything to your profile.")

    parser.add_argument("--source", "-s",
                        help = "The file to source in order to set the environment.  Default is SetPyRosettaEnvironment.sh",
                        default="SetPyRosettaEnvironment.sh")

    parser.add_argument("--outname", "-o",
                        help ="Name of the output file we will eventually source.  Default is .pyrosetta_env",
                        default = ".pyrosetta_env")

    parser.add_argument("--outdir", "-d",
                        help = "Output directory of the env file.  Default is $HOME directory")

    parser.add_argument("--shell-profile-path", "-p",
                        help = "If your shell profile path is different than the common ones, please set the path here.")




    options = parser.parse_args()

    source_environment_fname = options.source
    source_environment_dir = os.path.dirname(os.path.abspath(__file__))
    source_environment_path = os.path.join(source_environment_dir, source_environment_fname)

    outdir = os.environ["HOME"]

    if options.outdir:
        outdir = options.outdir

    if not os.path.exists(outdir):
        sys.exit("Outdir " + outdir + " does not exist! ")

    if not os.path.exists(source_environment_path):
        sys.exit("Source file does not exist! "+source_environment_path)


    #Allow the simple shell script to do its thing.  Then pull the PyRosetta path from it.
    os.system("source "+source_environment_path)
    pyrosetta = source_environment_dir
    print "PyRosetta at "+ pyrosetta
    outfile_path = os.path.join(outdir, options.outname)
    print "\nCreating "+ outfile_path
    OUTFILE = open(outfile_path, 'w')

    OUTFILE.write("# PyRosetta Paths.  Source this within your profile. \n")
    OUTFILE.write("# Date Created: "+ time.strftime("%m/%d/%Y")+"\n\n\n")
    OUTFILE.write("export PYROSETTA="+pyrosetta+"\n")
    OUTFILE.write("export PYTHONPATH=$PYROSETTA${PYTHONPATH+:$PYTHONPATH}\n")
    OUTFILE.write("export DYLD_LIBRARY_PATH=$PYROSETTA:$PYROSETTA/rosetta${DYLD_LIBRARY_PATH+:$DYLD_LIBRARY_PATH}\n")
    OUTFILE.write("export LD_LIBRARY_PATH=$PYROSETTA:$PYROSETTA/rosetta${LD_LIBRARY_PATH+:$LD_LIBRARY_PATH}\n")
    OUTFILE.write("export PYROSETTA_DATABASE=$PYROSETTA/database\n")
    OUTFILE.write("\n")
    OUTFILE.close()

    ##Now we figure out the bashrc or equivalent and then see if we already source it.  Append it

    rc = get_profile_fname()

    ##Append the file.
    if options.shell_profile_path:
        rc_file_path = options.shell_profile_path
    else:
        rc_file_path = os.path.join( os.environ["HOME"], rc)


    print "Reading "+rc_file_path

    READ_RCFILE = open(rc_file_path, 'r')
    RCFILE = open(rc_file_path, 'a')

    lines = READ_RCFILE.readlines()
    READ_RCFILE.close()

    source_line = "source "+outfile_path
    if re.search(source_line, " ".join(lines)):
        print "\nSource line already present in RC.  Complete!"
        print "  Re-Run this if you ever change paths to PyRosetta."
    else:
        append_rc = ask_append_rc()
        if append_rc:
            RCFILE.write("\n## PyRosetta ##\n")
            RCFILE.write(source_line+"\n")

            print "Complete - RC edited! "
            print "  Re-Run this if you ever move the PyRosetta directory..."
        else:
            print "Nothing to be done.  Please add source "+outfile_path +" to your profile"

    RCFILE.close()
    print "\nDone"


