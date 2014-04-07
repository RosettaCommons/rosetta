// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief

// libRosetta headers
#include <core/fragment/ConstantLengthFragSet.hh>

#include <core/options/option_macros.hh>
#include <core/options/option.hh>
#include <core/options/keys/OptionKeys.hh>
#include <core/options/keys/in.OptionKeys.gen.hh>
#include <core/options/keys/out.OptionKeys.gen.hh>

#include <devel/init.hh>
#include <core/types.hh>

#include <utility/excn/Exceptions.hh>

// option key includes


void register_options() {
  using namespace core::options;
  using namespace core::options::OptionKeys;

  OPT(in::file::silent);
  OPT(in::file::s);
  OPT(in::file::residue_type_set);
  OPT(out::nooutput);
}

int
main( int argc, char * argv [] ) {
    try {
        using namespace core::options;
        using namespace core::fragment;
        using namespace core::options::OptionKeys;

	devel::init(argc, argv);
	register_options();

	ConstantLengthFragSetOP frags = new ConstantLengthFragSet;
	std::string filename = "frags.9mers";
	frags->read_fragment_file( filename );
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
        return 0;

}
