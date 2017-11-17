// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief

// libRosetta headers
#include <protocols/domain_assembly/CombineChainsMover.hh>
#include <protocols/jd2/JobDistributor.hh>

#include <basic/options/option_macros.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>
#include <protocols/moves/Mover.hh>
#include <devel/init.hh>
#include <core/types.hh>

#include <utility/excn/Exceptions.hh>

// option key includes


void register_options() {
  using namespace basic::options;
  using namespace basic::options::OptionKeys;

  OPT(in::file::silent);
  OPT(in::file::s);
  OPT(in::file::residue_type_set);
  OPT(out::nooutput);
}

int
main( int argc, char * argv [] ) {
    try {
        using namespace basic::options;
        using namespace basic::options::OptionKeys;

	devel::init(argc, argv);
	register_options();

	protocols::domain_assembly::CombineChainsMoverOP  mover = new protocols::domain_assembly::CombineChainsMover;

	// execution
	protocols::jd2::JobDistributor::get_instance()->go(mover);
    } catch (utility::excn::Exception const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
        return 0;
    }
