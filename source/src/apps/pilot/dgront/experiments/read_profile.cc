// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief Simplest test on reading binary checkpoint file generated from psiblast

#include <basic/options/option_macros.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
//#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <devel/init.hh>
#include <core/types.hh>
#include <core/sequence/SequenceProfile.hh>

#include <utility/excn/Exceptions.hh>

OPT_KEY( String, chk_file )

OPT_1GRP_KEY( File, out, pdb )


void register_options() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  NEW_OPT( chk_file, "input file - binary checkpoint from psi-blast","aaa");
  NEW_OPT( out::pdb, "provides a file name for the output PDB file","out_pose.pdb" );
//  option.add_relevant(chk_file);
}

int main( int argc, char * argv [] ) {

    try {
	using namespace core;
        using namespace basic::options;
        using namespace basic::options::OptionKeys;

	register_options();
	devel::init(argc, argv);

//------------- Read the native pose  ----------
    if ( option[ chk_file ].user() ) {
                std::cout << "reading " << option[chk_file]()<< std::endl;
                core::sequence::SequenceProfile q_prof;
                q_prof.read_from_binary_chk(option[chk_file]());
    }


    } catch ( utility::excn::EXCN_Base const & e ) {
              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
    }
    return 0;
}

