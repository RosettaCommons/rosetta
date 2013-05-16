#!/usr/bin/python
# :noTabs=true:

# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

## @file   basic/options/options.py
## @brief  Program options generation script that is run to generate the options
## @author Sergey Lyskov (Sergey.Lyskov@jhu.edu)

# Script for generating option.gen.cc.hh, keys/OptionKeys.hh.gen.hh, keys/OptionKeys.cc.gen.hh files
# Use 'python options.py -Wiki' to generate a Wiki table

import sys

import options_class, options_rosetta
import os.path, os
#from difflib import Differ

Options = options_rosetta.Options
james_debug = 1

class KeepSameFile(object):
    def __init__(self,fname,opts):
        self.fname = fname
        self.opts = opts
        self.body = ""
    def write(self,s):
        self.body += s
    def close(self):
        ischanged = False
        try :
            with open(self.fname) as existing_file:
                if existing_file.read() != self.body:
                    ischanged = True
        except IOError :
            ischanged = True
        if ischanged:
            print "file",self.fname,"being updated"
            with open(self.fname,self.opts) as out:
                out.write(self.body)
            return 1
        return 0


def main(args):
    num_changed_files = 0
    if len(args) <= 1:  # no option give - just generating C++ files.

        # code below is for if we ever want to split options into separate groups
        # to speed up compilation times. totally untested!
        if james_debug:

            header1 = '// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-\n' + '// vi: set ts=2 noet:\n' + '// (c) Copyright Rosetta Commons Member Institutions.\n' + '// (c) This file is part of the Rosetta software suite and is made available under license.\n' + '// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.\n' + '// (c) For more information, see http://www.rosettacommons.org. Questions about this can be\n' + '// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.\n\n'
            header3 = '\n/// @brief  basic::options::OptionKeys collection\n' + '/// @author Stuart G. Mentzer (Stuart_Mentzer@objexx.com)\n' + '/// @author James M. Thompson (tex@u.washington.edu)\n' + '\n'
            header5 = '\n// Unit headers\n#include <basic/options/keys/OptionKeys.hh>\n\n' + 'namespace basic {\n' + 'namespace options {\n' + 'namespace OptionKeys {\n' + '\n'
            footer = '\n} // namespace OptionKeys\n' + '} // namespace options\n' + '} // namespace basic\n\n#endif\n'
            header_gen_hh = '#ifndef OPTION_CC_GEN_HH\n' +\
                '#define OPTION_CC_GEN_HH\n'+\
                '#include <basic/options/option.hh>\n' + \
                '#include <basic/options/option.cc.include.gen.hh>\n' + \
                '#include <utility/options/OptionCollection.hh>\n'

            footer_gen_hh = '\ninline void add_all_rosetta_options( utility::options::OptionCollection &option ) { add_rosetta_options_0(option); add_rosetta_options_1(option); add_rosetta_options_2(option); add_rosetta_options_3(option); }\n#endif\n'


            #def callback(arg, directory, files):
            #	for file in files:
            #		print os.path.join(directory, file), repr(arg)
            #os.path.walk(".", callback, "secret")

            output = {}
            output[ 'option.cc.gen.hh' ] = []
            output[ 'keys/OptionKeys.gen.hh' ] = {}
            output[ 'keys/OptionKeys.cc.gen' ] = []

            # Why is this here?  It seems to be unused.
            ##if not os.path.isdir('def'):
            #	os.mkdir( 'def' )


            for opt in Options:
                ns = str( opt.get_namespace(0) )

                # create new lists if necessary
                if not output[ 'keys/OptionKeys.gen.hh' ].has_key( ns ):
                    output[ 'keys/OptionKeys.gen.hh' ][ ns ] = []

                output[ 'option.cc.gen.hh' ].append( opt.getOptionCC() )
                output[ 'keys/OptionKeys.gen.hh' ][ ns ].append( opt.getOptionKeysHH() )
                output[ 'keys/OptionKeys.cc.gen' ].append( opt.getOptionKeysCC() )


            gen_hh_files = []
            gen_cc_files = []
            for file_prefix in output.keys():
                (dirname,filename) = os.path.split( file_prefix )
                if file_prefix == 'option.cc.gen.hh':
                    outfile = file_prefix
                    #print outfile
                    f = KeepSameFile(outfile, 'wb')
                    f.write( header_gen_hh )

                    split_len = len( output[ file_prefix ] ) / 4 + 1  # for now we generate 4 functions instead of 1
                    groups = [ output[ file_prefix ][i: i+split_len] for i in range(0, len(output[ file_prefix ]), split_len) ]
                    for i,g in enumerate(groups):
                            lines = 'inline void add_rosetta_options_%s( utility::options::OptionCollection &option ) {' % i
                            lines += "".join(g) + '\n}\n'
                            f.write( lines )
                    #lines = output[ file_prefix ]
                    f.write( footer_gen_hh )
                    num_changed_files += f.close()
                elif file_prefix == 'keys/OptionKeys.cc.gen':
                    split_len = len( output[ file_prefix ] ) / 4 + 1  # for now we split .cc just in four files
                    groups = [ output[ file_prefix ][i: i+split_len] for i in range(0, len(output[ file_prefix ]), split_len) ]
                    for i,g in enumerate(groups):
                        outfile = file_prefix + '%s.hh' % i
                        #print outfile
                        f = KeepSameFile(outfile, 'wb')
                        f.write( "".join(g) )
                        num_changed_files += f.close()
                else:
                    for ns in output[ file_prefix ].keys():
                        new_filename = ".".join( [ns, filename] )
                        outfile = os.path.join( dirname, new_filename )
                        header2  = '/// @file   basic/options/' + outfile
                        full_fn  = 'basic/options/keys/' + ns + '_OptionKeys_gen_HH'
                        inc_path = full_fn.replace('/','_')
                        inc_symb = 'INCLUDED_' + inc_path
                        header4  = '#ifndef ' + inc_symb + '\n#define ' + inc_symb + '\n'
                        gen_hh_files.append( new_filename )
                        header = [ header1, header2, header3, header4, header5 ]
                        output[ file_prefix ][ ns ].append( footer )
                        output[ file_prefix ][ ns ] = header + output[ file_prefix ][ ns ]
                        lines = "".join( output[ file_prefix ][ ns ] )
                        #print outfile
                        f = KeepSameFile(outfile, 'wb')
                        f.write( "".join(lines) )
                        num_changed_files += f.close()

            f = KeepSameFile('option.cc.include.gen.hh', 'wb')
            for include_file in gen_hh_files:
                f.write( '#include <basic/options/keys/' + include_file + '>\n' )
            num_changed_files += f.close()

        else:
            options_class.writeToFile(Options, 'option.cc.gen.hh', options_class.Option.getOptionCC)
            options_class.writeToFile(Options, 'keys/OptionKeys.gen.hh', options_class.Option.getOptionKeysHH)
            options_class.writeToFile(Options, 'keys/OptionKeys.cc.gen.hh', options_class.Option.getOptionKeysCC)

        # Generating Doxygen docs
        print "Generating Doxygen docs...",
        f = KeepSameFile("./../../../doc/options.dox", 'wb')
        f.write( options_class.getDoxygenPage(Options) )
        num_changed_files += f.close()
        print " Done!"
        print "number of files updated:",num_changed_files
        print "Total %s options." % len(Options)

    elif args[1] == '-Wiki':
        print options_class.printWikiTable(Options)


if __name__ == "__main__":  main(sys.argv)
