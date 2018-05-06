# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.

from os.path import exists, expandvars, abspath, basename
from glob import glob
import subprocess
import sys
import os
import shutil
import string
import os.path
import re
import math
import time
import multiprocessing

import imp
file_path = os.path.split( os.path.abspath(__file__) ) [0]
imp.load_source('measure_params', file_path + '/measure_params.py')
from measure_params import compute_torsion, compute_dist, compute_squared_dist, cross, dot, vec_diff

#####################################################
class ErraserError(Exception):
    pass
#####################################################
def error_exit(message = ""):
    sys.stdout.flush()
    sys.stderr.flush()
    raise ErraserError(message)
#####################################################
def rosetta_bin_path(exe_file, rosetta_folder = "") :
    """
    Return the absolute path to the given executable file in Rosetta.
    If no input is given, return the Rosetta bin folder path.
    Only the prefix should be given. Rosetta will figure out the rest.
    User can provides either the Rosetta root path or the binary path as optional input,
    or the function will search for $ROSETTA environment variable.
    """
    if rosetta_folder == "" :
        rosetta_folder = expandvars("$ROSETTA")
        if rosetta_folder == "$ROSETTA" :
            error_exit("USER need to set environmental variable $ROSETTA and pointed it to the Rosetta folder!")

    exe_folder = rosetta_folder + "/main/source/bin/" #Default Rosetta folder structure
    if not exists(exe_folder) : #Otherwise, assume the input folder name is bin path
        exe_folder = rosetta_folder
    check_path_exist(exe_folder)
    name_extensions = [ ".python.linuxgccrelease", ".linuxgccrelease", ".linuxclangrelease", "", ".macosgccrelease", ".macosclangrelease",
                       "failed_to_find_Rosetta_path"] #this makes a better error message if pathing fails
    exe_path = ""
    for name in name_extensions :
        exe_path = exe_folder + exe_file + name
        if exists(exe_path) :
            break

    if not exists(exe_path):
        # last chance: glob anything in the build directory.
        bins = glob(exe_folder + exe_file + '.*')[0]
        if len(bins) > 0:
            exe_path = bins[0]
    check_path_exist(exe_path)
    return exe_path

def initialize_rna_tools(rosetta_folder = ""):
    if rosetta_folder == "" :
        rosetta_folder = expandvars("$ROSETTA")
        if rosetta_folder == "$ROSETTA" :
            error_exit("USER need to set environmental variable $ROSETTA and pointed it to the Rosetta folder!")
    exe_folder = rosetta_folder + "/tools/rna_tools/"
    print ['source', exe_folder+'INSTALL']
    subprocess.Popen(['source', exe_folder+'INSTALL'], shell=False)
    return

def get_fasta(pdb_file, out_file, rosetta_folder = ""):
    """
    Gets the FASTA contents using pdb2fasta.py
    Note that initialize_rna_tools has been called already
    """

    out_list = subprocess.check_output(['/usr/bin/python', expandvars("$ROSETTA")+'/tools/rna_tools/bin/pdb2fasta.py', pdb_file], shell=False).split('\n')
    with open(out_file, 'w') as f:
        for l in out_list:
            f.write("{}\n".format(l))
    return out_list

#####################################################
def rosetta_database_path(rosetta_folder = "") :
    """
    Return the absolute path to the Rosetta database.
    If no input is given, return the Rosetta bin folder path.
    User can provides either the Rosetta root path or the database path as optional input,
    or the function will search for $ROSETTA environment variable.
    """
    if rosetta_folder == "" :
        rosetta_folder = expandvars("$ROSETTA")
        if rosetta_folder == "$ROSETTA" :
            error_exit("USER needs to set environmental variable $ROSETTA and point it to the Rosetta folder!")

    database_folder = rosetta_folder + "/main/database/" #Default Rosetta folder structure
    if  not exists(database_folder) : #Otherwise, assume the input folder name is database path
        database_folder = rosetta_folder
    check_path_exist(database_folder)
    return database_folder

###################################################
def subprocess_call(command, out = sys.stdout, err = sys.stderr, is_append_file = False):
    """
    Execute the given commands in /usr/bin/sh.
    Error exit if failed.
    """
    print "COMMAND:", command #" ".join(command)
    #print "initial LD_LIBRARY_PATH is ", expandvars("$LD_LIBRARY_PATH")
    #os.environ["LD_LIBRARY_PATH"] = "/net/cci/amw579/local_gcc_build/bin/usr/local/lib64/:/net/cci/amw579/local_gcc_build/bin/usr/local/lib/:" + os.environ["LD_LIBRARY_PATH"]
    #print "final LD_LIBRARY_PATH is ", expandvars("$LD_LIBRARY_PATH")
    out_channel = sys.stdout
    if type(out) is str :
        if is_append_file :
            out_channel = open(out, 'a')
        else :
            out_channel = open(out, 'w')
    elif type(out) is file :
        out_channel = out

    err_channel = sys.stderr
    if type(err) is str :
        if is_append_file :
            err_channel = open(err, 'a')
        else :
            err_channel = open(err, 'w')
    elif type(err) is file :
        err_channel = err

    out_channel.flush()
    err_channel.flush()
    try:
      subprocess.check_call(command, shell=False, stdout = out_channel, stderr = err_channel, env=os.environ)
    except OSError, e:
      print e
    out_channel.flush()
    err_channel.flush()
###################################################

def subprocess_out(command, err = sys.stderr):
    """
    Execute the given commands in /usr/bin/sh.
    Error exit if failed.
    Return the output as a list of each lines.
    """

    err_channel = sys.stderr
    if type(err) is str :
        err_channel = open(err, 'w')
    elif type(err) is file :
        err_channel = err

    err_channel.flush()
    out_list = subprocess.check_output(command, shell=True).split('\n')
    if out_list[-1] == '' :
        out_list.pop(-1)
    err_channel.flush()
    if err_channel != sys.stderr :
        os.fsync(err_channel)
    return out_list
###################################################

def grep(pattern, input_file) :
    """
    Search for the given pattern in the input file.
    Return a list of lines containing the pattern.
    """
    check_path_exist(input_file)
    out_lines = [ line[:-1] for line in open(input_file) if re.search( pattern, line ) ]

    return out_lines
###################################################

def remove(pattern, ignore_errors=True) :
    """
    Remove files/folder that agree with the pattern.
    """
    for path in glob(pattern):
        if os.path.isfile(path):
            os.remove(path)
        elif os.path.isdir(path):
            shutil.rmtree(path, ignore_errors=ignore_errors)
        else:
            error_exit("Not a file nor a dir!!!")
###################################################
def copy(file1, dst) :
    """
    Copy files to new destination. Replace with new one if the two overlaped.
    """
    check_path_exist(file1)
    if os.path.isfile(dst) :
        os.remove(dst)
    shutil.copy(file1, dst)
###################################################
def move(file1, dst) :
    """
    Move files to new destination.
    """
    check_path_exist(file1)
    if os.path.isfile(dst) :
        os.remove(dst)
    shutil.move(file1, dst)
###################################################
def get_total_res(pdbname):
    """
    Return the total number of residues ina pdb file.
    """
    netpdbname = pdbname
    check_path_exist(netpdbname)

    # AMW: be alert to the possibility that ACE/NME residue mergers will require
    # compensation: Rosetta will think fewer residues than this...
    oldresnum = ''
    count = 0;
    for line in open(netpdbname):
        if len(line) > 4 and line[0:4] == 'ATOM':
            resnum = line[21:27]
            if resnum != oldresnum :
                count = count + 1
                oldresnum = resnum

    return count

def find_nearby_res(input_pdb, input_res, dist_cutoff, reload = True):
    """
    Find nearby residues to the input residues by a distance_cutoff.
    All the residues should be in the same chain with continous numbering starting from 1.
    """
    check_path_exist(input_pdb)

    # If input_res is a list of int (as would be if called by SWA_rebuild)
    # then we need old load_pdb_coord. Pray for homogeneity.
    try :
        coord_all
        coord_C1
        if reload:
            if not isinstance(input_res[0], str):#input_res[0]isinstance(input_res[0], (int,long)):
                [coord_all, coord_C1] = load_pdb_coord_old(input_pdb)
            else:
                [coord_all, coord_C1] = load_pdb_coord(input_pdb)
    except :
        #if isinstance(input_res[0], (int,long)):
        if not isinstance(input_res[0], str):#input_res[0]isinstance(input_res[0], (int,long)):
            [coord_all, coord_C1] = load_pdb_coord_old(input_pdb)
        else:
            [coord_all, coord_C1] = load_pdb_coord(input_pdb)

    res_list = []
    if isinstance(coord_all, dict):
        for i in input_res:
            if not i in coord_all.keys(): #range(1, len(coord_all) + 1) :
                print i
                print coord_all.keys()
                error_exit("Input residues outside the range of pdb residues!")
            for j in coord_all.keys(): #range(1, len(coord_all) + 1) :
                if (j in input_res or j in res_list) : continue
                dist_C1 = compute_dist( coord_C1[i], coord_C1[j] )
                if dist_C1 > dist_cutoff + 8:
                    continue
                for coord_target_atom in coord_all[i] :
                    found_qualifying_atom = False
                    for coord_all_atom in coord_all[j] :
                        dist = compute_dist( coord_target_atom, coord_all_atom)
                        if dist < dist_cutoff:
                            res_list.append(j)
                            found_qualifying_atom = True
                            break
                    if found_qualifying_atom:
                        break

        res_list.sort()
        return res_list
    else:
        for i in input_res:
            if not i in range(1, len(coord_all) + 1) :
                print i
                print coord_all.keys()
                error_exit("Input residues outside the range of pdb residues!")
            for j in range(1, len(coord_all) + 1) :
                if (j in input_res or j in res_list) : continue
                dist_C1 = compute_dist( coord_C1[i-1], coord_C1[j-1] )
                if dist_C1 > dist_cutoff + 8:
                    continue
                for coord_target_atom in coord_all[i-1] :
                    found_qualifying_atom = False
                    for coord_all_atom in coord_all[j-1] :
                        dist = compute_dist( coord_target_atom, coord_all_atom)
                        if dist < dist_cutoff:
                            res_list.append(j)
                            found_qualifying_atom = True
                            break
                    if found_qualifying_atom:
                        break

        res_list.sort()
        return res_list


################################################################
def parse_options( argv, tag, default):
    """
    Parse various input options.
    """
    tag_argv = '-%s' % tag
    if tag_argv in argv :
        pos = argv.index( tag_argv )   ###Position of the option name

        if pos == ( len( argv ) - 1 ) or argv[ pos+1 ][0] == '-' :
            error_exit("Invalid parse_option input")

        if default == "False" or default == "True" : # Boolean
            if argv[ pos + 1 ] == "True" or argv[ pos + 1 ] == "true" :
                value = True
            elif argv[ pos + 1 ] == "False" or argv[ pos + 1 ] == "false" :
                value = False
            else:
                error_exit('(%s != "True") and  (%s != "False")' % (tag, tag))
        elif isinstance( default, str ) :
            value = argv[ pos + 1 ]
        elif isinstance( default, int ) :
            value = int( argv[ pos + 1 ] )
        elif isinstance( default, float ) :
            value = float( argv[ pos + 1 ] )
        else:
            error_exit("Invalid parse_option default")

        print "%s = %s" %(tag, argv[ pos + 1 ])
        return value

    else: #Return the default value
        print "%s = %s" %(tag, default)
        if default=="False" or default=="True" :
            if( default == "True"):
                return True
            elif( default == "False"):
                return False
            else:
                error_exit('(%s != "True") and  (%s != "False")' % (tag, tag))
        else:
            return default
###############################

def parse_option_int_list( argv, tag ):
    """
    Parse input cmdline option with the following format:
    ... -tag 1 5 7-20 40 45-60 ...
    Return a list of expanded numbers
    """
    list_load = []
    tag_argv = '-%s' % tag

    if not tag_argv in argv:
        print "%s = %s" % (tag, list_load)
        return list_load

    pos = argv.index(tag_argv)
    for i in range( pos + 1, len(argv) ) :
        input_string = argv[i]
        if input_string[0] == '-' :
            break
        split_string = input_string.split('-')

        for num in split_string :
            try :
                int(num)
            except :
                error_exit("Incorrect input for -%s: instance %s" % (tag, input_string) )

        if len(split_string) == 1 :
            list_load.append(int(split_string[0]))
        elif len(split_string) == 2 :
            start_num = int(split_string[0])
            end_num = int(split_string[1])
            for j in range(start_num, end_num + 1) :
                list_load.append(j)
        else:
            error_exit("Incorrect input for -%s: instance %s" % (tag, input_string) )

    print "%s = %s" % (tag, list_load)
    return list_load

###############################
def parse_option_string_list ( argv, tag ) :
    """
    Parse input cmdline option with the following format:
    ... -tag ss bbb fef ...
    Return a list of input strings
    """
    list_load = []
    tag_argv = '-%s' % tag

    if tag_argv in argv :
        pos = argv.index(tag_argv)
        print pos
        for i in range( pos + 1, len(argv) ) :
            input_string = argv[i]
            if input_string[0] == '-' :
                break
            list_load.append(input_string)
    print "%s = %s" % (tag, list_load)
    return list_load
##########################################

def parse_option_chain_res_list ( argv, tag ) :
    """
    Parse input cmdline option with the following format:
    ... -tag A10 B20-22...
    Return a list of strings with number expanded:
    ['A10', 'B20', 'B21', 'B22', ... ]
    AMW: we will have to expand this so that it takes seg ID at the end
    """
    print argv
    list_load = []
    tag_argv = '-%s' % tag

    if not tag_argv in argv :
        print "%s = %s" % (tag, list_load)
        return list_load

    pos = argv.index(tag_argv)
    #print pos
    for i in range( pos + 1, len(argv) ) :
        input_string = argv[i]
        if input_string[0] == '-' :
            break
        input_string = input_string.upper()

        if ':' in input_string:
            chain_id = input_string[0]
            input_string = input_string[2:]
        else:
            chain_id = input_string[0]
            input_string = input_string[1:]
        split_string = input_string.split('-')

        for num in split_string :
            try :
                int(num)
            except :
                error_exit("Incorrect input for -%s: instance %s" % (tag, argv[i]) )

        if len(split_string) == 1 :
            list_load.append('%s:%s' % (chain_id, split_string[0]) )
        elif len(split_string) == 2 :
            start_num = int(split_string[0])
            end_num = int(split_string[1])
            for j in range(start_num, end_num + 1) :
                list_load.append('%s:%s' % (chain_id, j))
        else :
            error_exit("Incorrect input for -%s: instance %s" % (tag, argv[i]) )

    print "%s = %s" % (tag, list_load)
    return list_load
#####################################################
def rna_rosetta_ready_set( input_pdb, out_name, option, rosetta_bin = "", rosetta_database = "", rna_prot_erraser = False ) :
    """
    Call Rosetta to read in a pdb and output the model right away.
    Can be used to ensure the file have the rosetta format for atom name, ordering and phosphate OP1/OP2.
    """
    print "Looking to see if {} exists...".format(input_pdb)
    check_path_exist(input_pdb)
    print "It must!"

    command = [rosetta_bin_path("erraser_minimizer", rosetta_bin)]
    command.extend(["-database", rosetta_database_path(rosetta_database)])
    command.extend(["-s", input_pdb])
    if rna_prot_erraser:
        command.extend(["-rna:rna_prot_erraser", "true", "-rna:corrected_geo", "true"])
    command.extend(["-ready_set_only", "true"])
    command.extend(["-ignore_unrecognized_res", "-inout:skip_connect_info", "true"]) # -ignore_waters"
    command.extend(["-in:guarantee_no_DNA", str(option.guarantee_no_DNA).lower()])
    command.extend(["-out:file:write_pdb_link_records", "true"])

    # calebgeniesse: output virtual phosphates here
    command.extend(["-output_virtual", "true"])
    print "######Start submitting the Rosetta command for rna_rosetta_ready_set########"
    subprocess_call( command, sys.stdout, sys.stderr )
    move( input_pdb.replace(".pdb", "_0001.pdb") , out_name )
    print "######Rosetta section completed#############################################"
    return True
#####################################################
def extract_pdb( silent_file, output_folder_name, rosetta_bin = "",
                 rosetta_database = "", extract_first_only = False, output_virtual = False,
                 rna_prot_erraser = False ):
    """
    Extract pdb's from Rosetta silent files.
    """

    check_path_exist( silent_file )
    silent_file_abs = abspath( silent_file )

    remove_variant_types = "false"

    ##########################
    tags_string = ""

    for line in open(silent_file):
        if 'SCORE' in line and (not 'description' in line) :
            tag = line.split()[-1]
            tags_string += ' ' + tag
            if extract_first_only :
                break

    if tags_string == '' :#Empty silent file
        sys.stderr.write("Extract pdb: Empty silent file or illegal format! Skip the extraction.\n")
        return

    #########################################################

    base_dir = os.getcwd()

    if not exists( output_folder_name ) :
        os.mkdir(output_folder_name)

    os.chdir( output_folder_name )
    other_pdbs = glob('*.pdb')

    command = rosetta_bin_path("rna_extract", rosetta_bin)

    command += " -tags " + tags_string
    command += " -in:file:silent " + silent_file_abs
    command += " -in:file:silent_struct_type rna"
    command += " -database %s" % rosetta_database_path(rosetta_database)
    command += " -remove_variant_cutpoint_atoms " + remove_variant_types
    command += " -output_virtual " + str(output_virtual).lower()
    if( rna_prot_erraser ) :
        command += " -rna:rna_prot_erraser true -rna:corrected_geo true "

    print "######Start submitting the Rosetta command for extract_pdb##################"
    subprocess_call( command.split(), sys.stdout, sys.stderr )
    print "######Rosetta section completed#############################################"

    ############## Reformat pdb files withoriginal naming conventions ############
    pdbs = [pdb for pdb in glob('S_*.pdb') if pdb not in other_pdbs]
    for pdb in pdbs:
        idx = pdb.replace('S_','').replace('.pdb','')
        if not idx.isdigit() or len(idx) >= 6:
            continue
        move(pdb, pdb.replace(idx, idx.zfill(6)))
    ##############################################################################

    os.chdir( base_dir )
    return True

#####################################################
def phenix_release_tag():
    output = subprocess_out('phenix.version')
    for line in output:
        if 'Release tag:' in line:
            return int( line.split()[-1] )
    return None

#####################################################
def phenix_rna_validate(input_pdb, outliers_only = True):
    """
    Parse output of phenix.rna_validate.
    """

    def is_number(x):
        try:
            int(x)
            return True
        except ValueError:
            return False

    check_path_exist(input_pdb)
    command = "phenix.rna_validate " + input_pdb
    if outliers_only is False:
        command += " outliers_only=False"
    output = subprocess_out(command)
    data_type = None
    data_types = ["pucker", "bond", "angle", "suite"]
    data_headers = {
        # "<header>" : "<type>"
        "Sugar pucker" : "pucker",
        "Backbone bond lengths" : "bond",
        "Backbone bond angles" : "angle",
        "Backbone torsion suites" : "suite",
        # legacy format used prior to phenix release 1703
        "Pucker Outliers:" : "pucker",
        "Bond Length Outliers:" : "bond",
        "Angle Outliers:" : "angle",
        "Suite Outliers:" : "suite",
        "Suite Validation:" : "suite"
    }
    data = dict([(type, []) for type in data_types])

    output = filter(None, output)
    for line_idx, line in enumerate(output):
        if any (header in line for header in data_headers):
            data_header = filter(line.count, data_headers.keys())[0]
            data_type = data_headers[data_header]
            continue
        if line.isalpha():
            continue
        if data_type and line.startswith("   "):
            if "pucker" in data_type:
                if "yes" not in line:
                    continue
            cols = line.strip().replace(':','  ').split()

            # This is intended to capture "if the third column is the residue name"
            # but that's not necessarily alpha or 1 long.
            # but the res num is always numeric, and it is only second
            # in the same circumstance
            #if cols[2].isalpha() and len(cols[2]) == 1:
            if is_number(cols[1]):
                cols.insert(0, cols.pop(2))

            # This can also happen in the first element of cols
            if len(cols[0]) > 4:
                c = cols.pop(0)
                chain, res = c[0], c[1:]
                cols.insert(1, chain)
                cols.insert(2, res)

            if len(cols[1]) > 4:
                c = cols.pop(1)
                chain, res = c[0], c[1:]
                cols.insert(1, chain)
                cols.insert(2, res)
            data[data_type].append( cols )

    return data
#####################################################
def find_error_res(input_pdb):
    """
    Return a list of error resdiue in the given pdb file (RNA only).
    Use phenix.rna_validate.
    """

    data = phenix_rna_validate( input_pdb )
    error_res = []

    for error_type, output in data.iteritems():
        for cols in output:
            chn = cols[1]
            res = int( cols[2] )
            if "suite" in error_type:
                suitename = cols[3]
                suiteness = float( cols[4] )
                if suitename == "__" or not suiteness < 0.1:
                    continue
                if res > 1:
                    error_res.append( "%s:%s" % ( chn, res - 1 ) )
            error_res.append( "%s:%s" % ( chn, res ) )

    error_res = list(set(sorted(error_res)))
    return error_res
#####################################################
def pdb2fasta(input_pdb, fasta_out, using_protein = False):
    """
    Extract fasta info from pdb.
    """
    check_path_exist(input_pdb)

    print "making the fasta", fasta_out

    longer_names={' rA': 'a', ' rC': 'c', ' rG': 'g', ' rU': 'u',
                  'ADE': 'a', 'CYT': 'c', 'GUA': 'g', 'URI': 'u',
                  'A  ': 'a', 'C  ': 'c', 'G  ': 'g', 'U  ': 'u',
                  '  A': 'a', '  C': 'c', '  G': 'g', '  U': 'u'}

    if (using_protein):
        longer_names.update( { 'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
                               'CYS': 'C', 'GLU': 'E', 'GLN': 'Q', 'GLY': 'G',
                               'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
                               'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
                               'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'} )

    output = open(fasta_out, 'w')
    output.write( ">%s\n" % os.path.basename(input_pdb) )
    oldresnum = ''
    for line in open(input_pdb) :
        if len(line) > 20 :
            if line[0:3] == 'TER':
                output.write('\n')
            if line[0:4] == 'ATOM':
                resnum = line[21:27]
                if not resnum == oldresnum:
                    longname = line[17:20]
                    if longer_names.has_key(longname):
                        output.write( longer_names[longname] )
                    else:
                        output.write( "X[{0}]".format(longname) )
                oldresnum = resnum
    output.write('\n')
    output.close()
    return True

##################################################
def pdb_slice(input_pdb, out_name, segment) :
    """
    Slice the input pdb file with given segment.
    Segment format example: '1-3 8 10-30 45 46 40'
    Or input a python int list as the segment.
    Return list of cut residue
    """
    check_path_exist(input_pdb)

    print "in pdb_slice", segment, input_pdb, out_name
    kept_res = []
    if isinstance(segment, list) :
        kept_res = segment
    elif isinstance(segment, str) :
        # AMW: we will want to make this also string-residue-compatible
        for elem in segment.split() :
            if '-' in elem :
                res_range = elem.split('-')
                for i in range(int(res_range[0]), int(res_range[1]) + 1) :
                    kept_res.append(i)
            else :
                kept_res.append( int(elem) )

    if len(kept_res) == 0 :
        error_exit("Invalid input of 'segment' option!")

    output = open(out_name, 'w')
    previous_res = -1
    previous_chn = ""
    old_res = 0
    new_res = 0
    new_atom = 0

    # The algorithm for figuring out what lines need to be output is hard:
    # the LINK lines are chosen based on resnum match to the input, but
    # they need to match the output resnum.
    # A normal PDB would have the LINK lines first, but let's test if that
    # is actually necessary. I don't think it is.
    # 1. Determine what atom lines should be output.
    # 2. Learn old_res <-> kept_res correspondence
    # 3. Output corresponding LINK records.
    old_kept_res = []
    atomlines = [ line for line in open(input_pdb) if len(line) > 20 and line[0:4] == 'ATOM' ]
    for line in atomlines:
        current_chn = line[21]
        current_res = int(line[22:26])
        if current_res != previous_res or current_chn != previous_chn :
            #old_res += 1
            previous_res = current_res
            previous_chn = current_chn
            if isinstance(kept_res[0], str):
                old_res = "%s:%d" % (previous_chn, previous_res)
            else:
                old_res = previous_res
            if old_res in kept_res :
                old_kept_res.append( current_res )
                new_res += 1

        if old_res in kept_res :
            new_atom += 1
            output.write('%s%7d%s%4d%s' % (line[0:4], new_atom, line[11:22], new_res, line[26:]) )

    output.close()
    return kept_res
#####################################
def check_path_exist(path_name) :
    if not exists(path_name) :
        error_exit("Path %s does not exist!" % path_name)
#####################################
def load_pdb_coord_old(input_pdb) :
    """
    Load in the pdb and return the coordinates for each heavy atom in each residue.
    Also return a list of C1' coordinates
    AMW: old version that supports SWA_rebuild better.
    """
    check_path_exist(input_pdb)

    current_chn = ''
    current_res = ''
    coord_all = []
    coord_res = []
    coord_C1 = []
    for line in open(input_pdb) :
        if len(line) < 22 :
            continue
        if line[0:4] != 'ATOM' :
            continue
        if line[13] == 'H' or line[12] == 'H' or line[77] == 'H':
            continue
        chn = line[21]
        res = line[22:27]
        crs = int(res.strip())

        if res != current_res or chn != current_chn:
            if current_res != '' :
                coord_res.append(coord_cur)
                coord_all.append(coord_res)
            current_res = res
            current_chn = chn
            coord_res = []

        coord_cur = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        coord_res.append(coord_cur)
        if line[13:16] == "C1'" or line[13:16] == 'CA ' :
            coord_C1.append(coord_cur)
    coord_all.append(coord_res)
    if len(coord_C1) != len(coord_all) :
        error_exit("Number of residues != number of C1'!!!")

    return [coord_all, coord_C1]
#####################################
def load_pdb_coord(input_pdb) :
    """
    Load in the pdb and return the coordinates for each heavy atom in each residue.
    Also return a list of C1' coordinates
    AMW: now returns a dict from name (i.e. A:1) to coords.
    """
    check_path_exist(input_pdb)

    current_chn = ''
    current_res = ''
    coord_all = {}
    coord_res = []
    coord_C1 = {}
    for line in open(input_pdb) :
        if len(line) < 22 :
            continue
        if line[0:4] != 'ATOM' :
            continue
        if line[13] == 'H' or line[12] == 'H' or line[77] == 'H':
            continue
        chn = line[21]
        res = line[22:27]
        crs = "%s:%s" % ( chn, res.strip() )

        if res != current_res or chn != current_chn:
            if current_res != '' :
                coord_all["%s:%s" % ( current_chn, current_res.strip() )] = coord_res
            current_res = res
            current_chn = chn
            coord_res = []

        coord_cur = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
        coord_res.append(coord_cur)
        if line[13:16] == "C1'" or line[13:16] == 'CA ' :
            coord_C1[crs] = coord_cur
    coord_all["%s:%s" % ( current_chn, current_res.strip() ) ] = coord_res
    if len(coord_C1) != len(coord_all) :
        error_exit("Number of residues != number of C1'!!!")

    return [coord_all, coord_C1]
#############################################################
def regularize_pdb(input_pdb, out_name) :
    """
    Regularize the residue and atom naming of input pdb file (for RNA part).
    """
    check_path_exist(input_pdb)

    rna_types = {' rA': '  A', ' rC': '  C', ' rG': '  G', ' rU': '  U',
                 ' Ar': '  A', ' Cr': '  C', ' Gr': '  G', ' Ur': '  U',
                 'ADE': '  A', 'CYT': '  C', 'GUA': '  G', 'URI': '  U',
                 'A  ': '  A', 'C  ': '  C', 'G  ': '  G', 'U  ': '  U',
                 '  A': '  A', '  C': '  C', '  G': '  G', '  U': '  U'}
    atom_name_convert = { " O1P":" OP1", " O2P":" OP2", " O3P":" OP3"}

    output_pdb = open(out_name, 'w')

    for line in open(input_pdb) :
        if len(line) > 6 and ( line[0:4] == 'ATOM' or line[0:6] == 'ANISOU' or line[0:6] == 'HETATM' ) :
            res_name = line[17:20]
            atom_name = line[12:16]
            if not res_name in rna_types:
                output_pdb.write(line)
            else:
                res_name = rna_types[res_name]
                if atom_name in atom_name_convert :
                    atom_name = atom_name_convert[atom_name]
                output_pdb.write(line[:12] + atom_name + ' ' + res_name + line[20:])
        else:
            output_pdb.write(line)
#############################################################
def regularize_OP1_OP2(input_pdb, out_name) :
    """
    Regularize the chirality for OP1/OP2 to match PDB standard.
    """
    check_path_exist(input_pdb)

    def output_cached_lines( out, cached_lines,
                             P_coord, OP1_coord, OP2_coord, O5prime_coord,
                             OP1_line_atom, OP1_line_anisou, OP2_line_atom, OP2_line_anisou ) :

        if len(cached_lines) == 0 :
            return
        is_switch_OPs = False
        if P_coord != [] and OP1_coord != [] and OP2_coord != [] and O5prime_coord != [] :
            dot_product = dot( cross( vec_diff(OP1_coord, P_coord), vec_diff(OP2_coord, P_coord) ),
                               vec_diff(O5prime_coord, P_coord) )
            if dot_product > 0 :
                is_switch_OPs = True

        for line1 in cached_lines :
            new_line = line1
            if is_switch_OPs and len(line) > 26 :
                if line1[12:16] == ' OP1' :
                    new_line = line1[0:26]
                    if line1[0:4] == 'ATOM' or line1[0:6] == 'HETATM' :
                        new_line += OP2_line_atom[26:]
                    elif line1[0:6] == 'ANISOU' :
                        new_line += OP2_line_anisou[26:]
                if line1[12:16] == ' OP2' :
                    new_line = line1[0:26]
                    if line1[0:4] == 'ATOM' or line1[0:6] == 'HETATM' :
                        new_line += OP1_line_atom[26:]
                    elif line1[0:6] == 'ANISOU' :
                        new_line += OP1_line_anisou[26:]
            out.write(new_line)

    out = open(out_name, 'w')

    cached_lines = []
    P_coord = []
    O5prime_coord = []
    OP1_coord = []
    OP2_coord = []
    old_res_id = []
    rna_types = [ '  A', '  C', '  G', '  U']
    OP1_line_atom = ''
    OP2_line_anisou = ''
    OP2_line_atom = ''
    OP1_line_anisou = ''

    for line in open(input_pdb) :
        is_output_cache = False
        if len(line) > 26 and ( line[0:4] == 'ATOM' or line[0:6] == 'ANISOU' or line[0:6] == 'HETATM' ) :
            res_name = line[17:20]
            atom_name = line[12:16]
            res_id = line[21:26]
            if (not res_name in rna_types) or res_id != old_res_id :
                is_output_cache = True
                old_res_id = res_id
        else :
            is_output_cache = True

        if is_output_cache and len(cached_lines) != 0 :
            output_cached_lines( out, cached_lines,
                                 P_coord, OP1_coord, OP2_coord, O5prime_coord,
                                 OP1_line_atom, OP1_line_anisou, OP2_line_atom, OP2_line_anisou )
            P_coord = []
            O5prime_coord = []
            OP1_coord = []
            OP2_coord = []
            cached_lines = []

        if len(line) > 26 and ( line[0:4] == 'ATOM' or line[0:6] == 'ANISOU' or line[0:6] == 'HETATM' ) :
            res_name = line[17:20]
            atom_name = line[12:16]
            res_id = line[21:26]
            if res_name in rna_types :
                if atom_name == ' OP1' :
                    if line[0:4] == 'ATOM' or line[0:6] == 'HETATM' :
                        OP1_line_atom = line
                        OP1_coord.append( float( line[30:38] ) )
                        OP1_coord.append( float( line[38:46] ) )
                        OP1_coord.append( float( line[46:54] ) )
                    elif line[0:6] == 'ANISOU' :
                        OP1_line_anisou = line
                if atom_name == ' OP2' :
                    if line[0:4] == 'ATOM' or line[0:6] == 'HETATM' :
                        OP2_line_atom = line
                        OP2_coord.append( float( line[30:38] ) )
                        OP2_coord.append( float( line[38:46] ) )
                        OP2_coord.append( float( line[46:54] ) )
                    elif line[0:6] == 'ANISOU' :
                        OP2_line_anisou = line
                if atom_name == ' P  ' and (line[0:4] == 'ATOM' or line[0:6] == 'HETATM') :
                    P_coord.append( float( line[30:38] ) )
                    P_coord.append( float( line[38:46] ) )
                    P_coord.append( float( line[46:54] ) )
                if atom_name == " O5'" and (line[0:4] == 'ATOM' or line[0:6] == 'HETATM') :
                    O5prime_coord.append( float( line[30:38] ) )
                    O5prime_coord.append( float( line[38:46] ) )
                    O5prime_coord.append( float( line[46:54] ) )
        cached_lines.append(line)

    output_cached_lines( out, cached_lines,
                         P_coord, OP1_coord, OP2_coord, O5prime_coord,
                         OP1_line_atom, OP1_line_anisou, OP2_line_atom, OP2_line_anisou )

#############################################################
def rosetta2phenix_merge_back(orig_pdb, rs_pdb, out_name) :
    """
    Merge the rosetta pdb back to the phenix refined pdb.
    """
    check_path_exist(orig_pdb)
    check_path_exist(rs_pdb)

    rna_res = ["  G", "G  ", " rG", "GUA",
               "  A", "A  ", " rA", "ADE",
               "  U", "U  ", " rU", "URI",
               "  C", "C  ", " rC", "CYT"]


    grep = subprocess.Popen(["grep", "-r", "IO_STRING",
                         "{0}chemical/residue_type_sets/fa_standard/residue_types/nucleic/rna_phenix/".format(
                             rosetta_database_path())], stdout=subprocess.PIPE)
    rna_res.append(subprocess.check_output(["awk", "{print $2}"], stdin=grep.stdout).split("\n"))
    grep = subprocess.Popen(["grep", "-r", "IO_STRING",
                         "{0}chemical/residue_type_sets/fa_standard/residue_types/nucleic/rna_nonnatural/".format(
                             rosetta_database_path())], stdout=subprocess.PIPE)
    rna_res.append(subprocess.check_output(("awk", "{print $2}"), stdin=grep.stdout).split("\n"))


    def is_line_rna (line) :
        res_name = line[17:20]
        return (res_name in rna_res)

    ###################################
    output_pdb = open(out_name, 'w')
    lines_rs = open(rs_pdb).readlines()
    lines_orig = open(orig_pdb).readlines()
    atom_name_convert = { " O1P":" OP1", " O2P":" OP2", " O3P":" OP3",
                          "1H5'":" H5'", "2H5'":"H5''", "1H2'":" H2'",
                          "1H2 ":" H21", "2H2 ":" H22", "2HO'":"HO2'",
                          "3HO'":"HO3'", "5HO'":"HO5'", "1H4 ":" H41",
                          "2H4 ":" H42", "1H6 ":" H61", "2H6 ":" H62" }

    rna_types = {' rA': '  A', ' rC': '  C', ' rG': '  G', ' rU': '  U',
                 'ADE': '  A', 'CYT': '  C', 'GUA': '  G', 'URI': '  U',
                 'A  ': '  A', 'C  ': '  C', 'G  ': '  G', 'U  ': '  U',
                 '  A': '  A', '  C': '  C', '  G': '  G', '  U': '  U'}

    current_res = ''
    current_chain = ''
    res_rs = ''
    atom_index = 0
    for line_orig in lines_orig :
        header = line_orig[0:6]
        if header == 'ATOM  ' and is_line_rna(line_orig) :
            res_name_orig = line_orig[17:20]
            res_orig = line_orig[21:27]
            atom_name_orig = line_orig[12:16]
            if res_orig != current_res :
                lines_rs_temp = lines_rs[:]
                for line_rs in lines_rs_temp :
                    if (line_rs[0:4] != 'ATOM' or (not is_line_rna(line_rs) ) ) :
                        lines_rs.remove(line_rs)
                        continue
                    res = line_rs[21:27]
                    if res_rs != res :
                        res_rs = res
                        break
                    else :
                        if atom_index < 99999:
                            atom_index += 1
                        atom_name = line_rs[12:16].replace('*', "'")
                        if atom_name_convert.has_key(atom_name) :
                            atom_name = atom_name_convert[atom_name]
                        output_line = (line_rs[0:6] + str(atom_index).rjust(5) + ' ' + atom_name +
                        line_rs[16:21] + current_chain + current_res + line_rs[27:])
                        output_line = output_line.replace('rG',' G')
                        output_line = output_line.replace('rA',' A')
                        output_line = output_line.replace('rU',' U')
                        output_line = output_line.replace('rC',' C')
                        output_pdb.write(output_line)
                        lines_rs.remove(line_rs)
                current_res = res_orig

            if atom_index < 99999:
                atom_index += 1
            output_line = ""
            if not res_name_orig in rna_types:
                # could be a nonnatural rna
                output_line = (line_orig[0:6] + str(atom_index).rjust(5) + line_orig[11:17] +
                           res_name_orig + line_orig[20:])
                print output_line
            else:
                output_line = (line_orig[0:6] + str(atom_index).rjust(5) + line_orig[11:17] +
                           rna_types[res_name_orig] + line_orig[20:])

            lines_rs_temp = lines_rs[:]
            for line_rs in lines_rs_temp :
                if (line_rs[0:4] != 'ATOM' or (not is_line_rna(line_rs) ) ) :
                    lines_rs.remove(line_rs)
                    continue
                res = line_rs[21:27]
                res_name = line_rs[17:20]

                atom_name = line_rs[12:16].replace('*', "'")
                if atom_name_convert.has_key(atom_name) :
                    atom_name = atom_name_convert[atom_name]

                #if (rna_types[res_name] != rna_types[res_name_orig]) : continue
                if (atom_name != atom_name_orig) : continue
                if (res != res_rs) : break
                output_line = (output_line[0:27] + line_rs[27:55] + output_line[55:])
                lines_rs.remove(line_rs)
                break

            output_pdb.write(output_line)

        elif header == 'ANISOU' :
            output_line = ''
            res_name_orig = line_orig[17:20]
            if is_line_rna(line_orig) :
                output_line = ""
                if not res_name_orig in rna_types:
                    # could be a nonnatural rna
                    output_line = (line_orig[0:6] + str(atom_index).rjust(5) + line_orig[11:17] +
                                   res_name_orig + line_orig[20:])
                    print output_line
                else:
                    output_line = (line_orig[0:6] + str(atom_index).rjust(5) + line_orig[11:17] +
                                   rna_types[res_name_orig] + line_orig[20:])
            else :
                output_line = (line_orig[0:6] + str(atom_index).rjust(5) + line_orig[11:])
            output_pdb.write(output_line)
        elif (header == 'HETATM' or
            (header == 'ATOM  ' and (not is_line_rna(line_orig) ) ) ):
            res_name_orig = line_orig[17:20]
            if atom_index < 99999:
                atom_index += 1
            output_line = (line_orig[0:6] + str(atom_index).rjust(5) + line_orig[11:])
            output_pdb.write(output_line)
        elif header == 'LINK':
            output_pdb.write(line_orig)
        elif (header != 'MASTER' and header[0:3] != 'END' and header != 'ANISOU') :
            output_pdb.write(line_orig)

    output_pdb.write('END                                ')
    return True

###############################################
def rosetta2std_pdb (input_pdb, out_name, cryst1_line = "") :
    """
    Convert the rosetta pdb file to standard pdb format. (RNA only)
    Prepend the CRYST1 line if found in the pdb file or given in input.
    """
    check_path_exist(input_pdb)

    atom_name_convert = { " O1P":" OP1", " O2P":" OP2", " O3P":" OP3",
                          "1H5'":" H5'", "2H5'":"H5''", "1H2'":" H2'",
                          "1H2 ":" H21", "2H2 ":" H22", "2HO'":"HO2'",
                          "3HO'":"HO3'", "5HO'":"HO5'", "1H4 ":" H41",
                          "2H4 ":" H42", "1H6 ":" H61", "2H6 ":" H62" }

    res_name_convert = {" rG":"  G", " rA":"  A", " rU":"  U", " rC":"  C"}
    output = open(out_name, 'w')

    searched_lines = grep( "CRYST1", input_pdb)
    if len(searched_lines) != 0 :
        cryst1_line = searched_lines[0]

    if cryst1_line != "" :
        output.write("%s\n" % cryst1_line)
    for line in open(input_pdb) :
        output.write(line)
        #if len(line) > 2 and (line[0:3] == 'END' or line[0:3] == 'TER') :
        #    output.write(line)
        #elif line[0:4] == 'LINK':
        #    output.write(line)
        #elif len(line) > 5 and (line[0:6] == 'HETATM' or line[0:4] == 'ATOM') :
        #    output.write(line)
        # AMW: principle -- let Rosetta handle naming
        #elif len(line) > 3 and line[0:4] == 'ATOM' :
        #    new_line = line.replace("*", "'")
        #
        #    atom_name = new_line[12:16]
        #    if atom_name_convert.has_key(atom_name) :
        #        atom_name = atom_name_convert[atom_name]
        #
        #    res_name = line[17:20]
        #    if res_name_convert.has_key(res_name) :
        #        res_name = res_name_convert[res_name]
        #
        #    atom_type = line[13]
        #    if line[12] == 'H' :
        #        atom_type = 'H'
        #
        #    new_line_list = list(new_line)
        #    new_line_list[12:16] = atom_name
        #    new_line_list[17:20] = res_name
        #    new_line_list[77] = atom_type
        #
        #    output_line = string.join(new_line_list, '')
        #    output.write(output_line)

    output.close()
    return True
####################################
def pdb2rosetta (input_pdb, out_name, alter_conform = 'A', PO_dist_cutoff = 2.0, use_rs_atom_res_name = False, using_protein = False, renumbering=True) :
    """
    Convert regular pdb file to rosetta format.
    Return a list for converting original residues # into new ones,
    a list of residue convalently bonded to removed HETATM so they can be fixed in further refinement,
    a list of cutpoints (terminal residues) in the structure,
    and the CRYST1 line in the model.
    """

    check_path_exist(input_pdb)

    atom_name_convert = { " OP1":" O1P", " OP2":" O2P", " OP3":" O3P",
                          " H5'":"1H5'", "H5''":"2H5'", " H2'":"1H2'",
                          " H21":"1H2 ", " H22":"2H2 ", "HO2'":"2HO'",
                          "HO3'":"3HO'", "HO5'":"5HO'", " H41":"1H4 ",
                          " H42":"2H4 ", " H61":"1H6 ", " H62":"2H6 " }

    res_name_convert = { "  G":" rG", "G  ":" rG", "GUA":" rG", " rG":" rG", "rG ":" rG",
                         "  A":" rA", "A  ":" rA", "ADE":" rA", " rA":" rA", "rA ":" rA",
                         "  U":" rU", "U  ":" rU", "URI":" rU", " rU":" rU", "rU ":" rU",
                         "  C":" rC", "C  ":" rC", "CYT":" rC", " rC":" rC", "rC ":" rC",}

    rna_names = []
    grep = subprocess.Popen( ["grep", "-r", "IO_STRING", "{0}chemical/residue_type_sets/fa_standard/residue_types/nucleic/rna_phenix/".format(rosetta_database_path())], stdout=subprocess.PIPE )
    rna_names = subprocess.check_output(["awk", "{print $2}"], stdin=grep.stdout ).split("\n")
    grep = subprocess.Popen( ["grep", "-r", "IO_STRING", "{0}chemical/residue_type_sets/fa_standard/residue_types/nucleic/rna_nonnatural/".format(rosetta_database_path())], stdout=subprocess.PIPE )
    rna_names += subprocess.check_output(("awk", "{print $2}"), stdin=grep.stdout ).split("\n")

    #print rna_names

    protein_res_names = ['ALA', 'ARG', 'ASN', 'ASP',
                         'CYS', 'GLU', 'GLN', 'GLY',
                         'HIS', 'ILE', 'LEU', 'LYS',
                         'MET', 'PHE', 'PRO', 'SER',
                         'THR', 'TRP', 'TYR', 'VAL']

    # append contents of l-ncaa; ignore peptoids for now
    grep = subprocess.Popen( ["grep", "-r", "IO_STRING", "{0}chemical/residue_type_sets/fa_standard/residue_types/l-ncaa/".format(rosetta_database_path())], stdout=subprocess.PIPE )
    #protein_res_names.append( subprocess.check_output(("awk", "{print $2}"), stdin=grep.stdout ).split("\n") )
    protein_res_names += subprocess.check_output(("awk", "{print $2}"), stdin=grep.stdout).split("\n")

    #print protein_res_names

    res_conversion_list = []
    fixed_res_list = []
    cutpoint_list = []
    CRYST1_line = ''
    O3prime_coord_pre = []
    O3prime_coord_cur = []
    P_coord_cur = []
    is_previous_het = False
    is_current_het = False
    previous_res = ""
    res_no = 0
    atm_no = 0

    output = open(out_name, 'w')
    for line in open(input_pdb) :
        if len(line) <= 21 :
            continue
        if line[0:6] == 'CRYST1' :
            CRYST1_line = line[:-1]
        elif line[0:4] == 'LINK' :
            output.write(line)
        elif line[0:6] == 'ATOM  ' or line[0:6] == 'HETATM' :
            res_name = line[17:20]
            if line[0:6] == 'ATOM  ' :
                if res_name_convert.has_key(res_name) :
                    if use_rs_atom_res_name :
                        res_name = res_name_convert[res_name]
                elif using_protein and (res_name in protein_res_names):
                    pass
                elif not res_name in rna_names:
                    continue

            current_res = line[21:27]

            if current_res != previous_res:
                previous_res = current_res
                is_previous_het = is_current_het
                if (not is_previous_het) and (not is_current_het) :
                    if len(O3prime_coord_cur) == 0 or len(P_coord_cur) == 0 :
                        if res_no - 1 > 0 :
                            cutpoint_list.append(res_no - 1)

                is_current_het = ( line[0:6] == 'HETATM' )
                if len(O3prime_coord_cur) != 0 and len(O3prime_coord_cur) != 3 :
                    error_exit("More than one O3' in a residue!!")
                else :
                    O3prime_coord_pre = O3prime_coord_cur
                O3prime_coord_cur = []
                P_coord_cur = []
                # also increment for NC residues.
                if not is_current_het or ( res_name in protein_res_names or res_name in rna_names ):
                    res_no += 1
                    if line[71:75].strip() != '':
                        orig_res = '%s:%s:%s' % (line[21], line[71:75].strip(), line[22:27])
                    else:
                        orig_res = '%s:%s' % (line[21], line[22:27])
                    orig_res = orig_res.replace(' ', '')
                    res_conversion_list.append(orig_res)

            if line[16] != ' ' and line[16] != alter_conform :
                continue

            atom_name = line[12:16]
            if use_rs_atom_res_name :
                if atom_name_convert.has_key(atom_name) :
                    atom_name = atom_name_convert[atom_name]
                atom_name = atom_name.replace("'", "*")

            if atom_name == " O3*" or atom_name == " O3'" :
                O3prime_coord_cur = [float(line[30:38]), float(line[38:46]), float(line[46:54])]

            if atom_name == " P  " :
                P_coord_cur = [float(line[30:38]), float(line[38:46]), float(line[46:54])]

                if (not is_current_het) and is_previous_het :
                    if res_no - 1 > 0 :
                        cutpoint_list.append( res_no - 1 )

                if len(O3prime_coord_pre) == 3 :
                    PO_dist = compute_dist( P_coord_cur, O3prime_coord_pre )
                    if PO_dist > PO_dist_cutoff :
                        if (not is_previous_het) and (not is_current_het) :
                            if res_no - 1 > 0 :
                                cutpoint_list.append( res_no - 1 )
                    else :
                        if is_current_het and ( not is_previous_het ) :
                            if res_no > 0 :
                                fixed_res_list.append( res_no )

                        if ( not is_current_het ) and is_previous_het :
                            if res_no > 0 :
                                fixed_res_list.append( res_no )

            # HETATM records that contain protein or RNA residues must be okay.
            #print rna_names
            #print protein_res_names
            if not is_current_het or ( res_name in protein_res_names or res_name in rna_names ):
                # Can't support atom numbers bigger than 99999 in pdb format.
                if atm_no < 99999:
                    atm_no += 1
                if len(line) < 80 :
                    line = line[:-1]
                    for i in xrange( len(line), 81 ) :
                        line += ' '
                    line += '\n'
                line_list = list(line)
                line_list[6:11] = str(atm_no).rjust(5)
                line_list[12:16] = atom_name
                line_list[16] = ' '
                line_list[17:20] = res_name
                if renumbering:
                    line_list[21] = 'A'
                    line_list[22:26] = str(res_no).rjust(4)
                line_list[26] = ' '
                line_list[55:60] = " 1.00"
                output.write( string.join(line_list, '') )

    output.close()
    return [res_conversion_list, fixed_res_list, cutpoint_list, CRYST1_line]
##################################################
def res_num_convert(res_conversion_list, res_pdb):
    """
    Convert pdb residue number (A15, chainID + resnum) to Rosetta pdb numbering.
    """
    #print res_pdb
    #print res_conversion_list
    if res_pdb not in res_conversion_list:
        error_exit("The residue '%s' was not found in the model.  Please make sure " % res_pdb +
        "that you include both the chain ID and residue number, for example 'A15'." )
    else:
        return res_conversion_list.index(res_pdb) + 1
##################################################
def res_wise_rmsd(pdb1, pdb2) :
    """
    Compute the rmsd of each residue from two pdb file. No extra alignment is performed.
    Return a list of rmsd of each residue.
    Assuming the pdb files have the same ordering of res name / atom name.
    """
    check_path_exist(pdb1)
    check_path_exist(pdb2)

    coord_pdb1 = load_pdb_coord(pdb1) [0]
    coord_pdb2 = load_pdb_coord(pdb2) [0]

    if len(coord_pdb1.items()) != len(coord_pdb2.items()) :
        error_exit("Two pdbs have different # of residues!!!")

    res_rmsd_list = []
    for res, coords in coord_pdb1.iteritems():
        if len(coord_pdb1[res]) != len(coord_pdb2[res]) :
            error_exit("Residue %s have different # of atoms in the two pdbs!!!" % (i+1) )
        res_rmsd = 0.0
        for j in range(0, len(coord_pdb1[res])) :
            res_rmsd += compute_squared_dist(coord_pdb1[res][j], coord_pdb2[res][j])
        res_rmsd = math.sqrt(res_rmsd / len(coord_pdb1[res]))
        res_rmsd_list.append([res, res_rmsd])

    return res_rmsd_list

##################################################
def pdb_slice_into_chunks(input_pdb, n_chunk) :
    """
    Slice the pdb into smaller chunks.
    Return a list of res_list for each chunk.
    """

    if n_chunk == 1 :
        return []
    check_path_exist(input_pdb)

    def fill_gaps_and_remove_isolated_res(res_list, total_res) :
        res_list.sort()

        res_list.sort()
        res_remove = []
        for i in range(0, len(res_list)) :
            if i == 0 and res_list[i+1] - res_list[i] != 1 :
                res_remove.append(res_list[i])
            elif i == len(res_list) - 1 and res_list[i] - res_list[i-1] != 1 :
                res_remove.append(res_list[i])
            elif res_list[i] - res_list[i-1] != 1 and res_list[i+1] - res_list[i] != 1 :
                res_remove.append(res_list[i])

        for res in res_remove :
            res_list.remove(res)

        new_res = []
        for i in range(0, len(res_list)) :
            if i == len(res_list) - 1 :
                if total_res - res_list[i] > 0 and  total_res - res_list[i] <= 4 :
                    new_res += range(res_list[i] + 1, total_res + 1)
            elif res_list[i+1] - res_list[i] > 1 and  res_list[i+1] - res_list[i] <= 4 :
                new_res += range(res_list[i] + 1, res_list[i+1])
        res_list += new_res

        return res_remove

    total_res = get_total_res(input_pdb)
    print "Found %s residues in %s" % (total_res, input_pdb)
    res_list_sliced = []
    sliced_list_final = []
    res_list_unsliced = range(1, total_res + 1)
    res_list_current = []
    chunk_size = 0
    n_chunk_left = n_chunk
    while len(res_list_unsliced) != 0 :
        if len(res_list_current) == 0 :
            res_list_current.append( res_list_unsliced.pop(0) )
            chunk_size = int( len(res_list_unsliced) / n_chunk_left)
            n_chunk_left -= 1

        res_list_new = find_nearby_res(input_pdb, res_list_current, 3.5, False)
        res_remove = []
        for res in res_list_new :
            if res in res_list_sliced :
                res_remove.append(res)
        for res in res_remove :
            res_list_new.remove(res)

        if len(res_list_new) == 0 and (not len(res_list_current) >= chunk_size * 0.7) :
            while True :
                res = res_list_unsliced.pop(0)
                if not res in res_list_current :
                    res_list_new.append(res)
                    break

        res_list_current += res_list_new

        if len(res_list_current) >= chunk_size or len(res_list_new) == 0 :
            res_list_current.sort()
            fill_gaps_and_remove_isolated_res(res_list_current, total_res)

            res_remove = []
            for res in res_list_current :
                if res in res_list_unsliced :
                    res_remove.append(res)
            for res in res_remove :
                res_list_unsliced.remove(res)

            res_list_current = sorted(res_list_current)
            sliced_list_final.append( res_list_current )
            res_list_sliced += res_list_current
            if n_chunk_left == 1 :
                res_list_current = res_list_unsliced
                removed_res = fill_gaps_and_remove_isolated_res(res_list_current, total_res)
                res_list_current = sorted(res_list_current)
                sliced_list_final.append( res_list_current )
                while len(removed_res) != 0 :
                    res_remove = []
                    for res in removed_res :
                        res_list_near = find_nearby_res(input_pdb, [res], 2.0, False)
                        for res_list in sliced_list_final :
                            is_break = False
                            for res_near in res_list_near :
                                if res_near in res_list:
                                    res_list.append(res)
                                    res_list.sort()
                                    res_remove.append(res)
                                    is_break = True
                                    break
                            if is_break :
                                break

                    for res in res_remove :
                        removed_res.remove(res)
                break
            else :
                res_list_current = []
                res_list_current.append( res_list_unsliced.pop(0) )
                chunk_size = int( len(res_list_unsliced) / n_chunk_left * 1.1)
                n_chunk_left -= 1

    return sliced_list_final
####################################
def pdb_slice_with_patching( input_pdb, out_name, slice_res_list ) :
    """
    Slice the pdb file with given input residue list.
    Patch the nearby residues to the sliced pdb.
    Return a sorted new residue number list and a list of patched res (in new order).
    """
    check_path_exist( input_pdb )

    dist_cutoff = 5.0

    print 'SLICE_RES_LIST: ', slice_res_list
    patched_res = find_nearby_res( input_pdb, slice_res_list, dist_cutoff )
    print 'PATCHED_RES: ', patched_res
    patched_res.sort()
    all_res = slice_res_list + patched_res
    all_res.sort()
    print 'ALL_RES ', all_res
    patched_res_new = []
    for res in patched_res :
        patched_res_new.append( all_res.index(res) + 1 )
    print 'PATCHED_RES_NEW ', patched_res_new

    pdb_slice( input_pdb, out_name, all_res )
    return [all_res, patched_res_new]
####################################
def sliced2orig_merge_back( orig_pdb, new_pdb, out_name, res_list ) :
    """
    Merge processed sliced segment back to original pdb (Rosetta format).
    """
    check_path_exist( orig_pdb )
    check_path_exist( new_pdb )

    out = open(out_name, 'w')
    new_pdb_read = open(new_pdb)
    new_pdb_line = new_pdb_read.readline()
    res_new_pre = -10
    res_orig_pre = -10
    is_residue_done = False
    atom_num = 0
    for line in open(orig_pdb) :
        if len(line) < 4 :
            out.write(line)
            continue
        if line[0:4] == 'ATOM' :
            chain = line[21]
            res_num = int(line[22:26])
            if res_num != res_orig_pre :
                res_orig_pre = res_num
                is_residue_done = False
            if is_residue_done :
                continue
            if not res_num in res_list :
                atom_num += 1
                out.write('%s%7d%s' %(line[0:4], atom_num, line[11:]) )
            else :
                if len(new_pdb_line) > 4 and new_pdb_line[0:4] == 'ATOM' :
                    if new_pdb_line[21].isspace():
                        new_pdb_line = new_pdb_line[:21] + chain + new_pdb_line[22:]
                    atom_num += 1
                    out.write('%s%7d%s%4d%s' % (new_pdb_line[0:4], atom_num, new_pdb_line[11:22], res_num, new_pdb_line[26:]) )
                    res_new_pre = int( new_pdb_line[22:26] )

                while True:
                    new_pdb_line = new_pdb_read.readline()
                    if new_pdb_line == "" :
                        break
                    if len(new_pdb_line) > 4 and new_pdb_line[0:4] == 'ATOM' :
                        res_new = int( new_pdb_line[22:26] )
                        if res_new == res_new_pre :
                            if new_pdb_line[21].isspace():
                                new_pdb_line = new_pdb_line[:21] + chain + new_pdb_line[22:]
                            atom_num += 1
                            out.write('%s%7d%s%4d%s' % (new_pdb_line[0:4], atom_num, new_pdb_line[11:22], res_num, new_pdb_line[26:]) )
                            res_new_pre = int( new_pdb_line[22:26] )
                        else :
                            break
                is_residue_done = True
    out.close()
    return True
####################################
def find_chi_angle( input_pdb, res ) :
    """
    Find the chi angle value for given residue of the input pdb.
    """
    def coord_from_atm_name( atm_name, atm_coords_list ) :
        print atm_coords_list
        for i in atm_coords_list :
            if atm_name == i[0] :
                return i[1]
        error_exit("No proper atom name: %s is found!!" % atm_name)

    atm_coords_list = []
    is_purine = True
    print input_pdb
    with open(input_pdb) as foo:
        print foo.readlines()
    #atomlines = ( line for line in open(input_pdb) if len(line) > 6 and line[0:4] == 'ATOM' ) #and int( line[22:26] ) == int(res[2:]) and line[21] == res[0] )
    atomlines = [ line for line in open(input_pdb) if len(line) > 6 and line[0:4] == 'ATOM' ] #and int( line[22:26] ) == int(res[2:]) and line[21] == res[0] )
    for line in atomlines:
        print "|%s|" % line[22:26]
        print "|%s|" % line[21]
        atm_name = line[12:16].replace(' ', '')
        res_name = line[19]
        is_purine = res_name in 'GA'
        coord = [ float( line[30:38] ), float( line[38:46] ), float( line[46:54] ) ]
        atm_coords_list.append( [atm_name, coord] )
    atom1 = coord_from_atm_name("O4'", atm_coords_list)
    atom2 = coord_from_atm_name("C1'", atm_coords_list)
    if is_purine :
        atom3 = coord_from_atm_name('N9', atm_coords_list)
        atom4 = coord_from_atm_name('C4', atm_coords_list)
    else :
        atom3 = coord_from_atm_name('N1', atm_coords_list)
        atom4 = coord_from_atm_name('C2', atm_coords_list)
    return compute_torsion(atom1, atom2, atom3, atom4)


