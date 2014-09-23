// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @ r_dock_tempered.cc
/// @ author Oliver Lange


// Rosetta headers
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/util.hh>
#include <devel/coupled_sidechains/CoupledSidechainProtocol2.hh>

#include <devel/init.hh>
#include <core/types.hh>

// Utility headers
#include <utility/vector1.hh>

// option key includes
#include <basic/options/option.hh>
#include <utility/excn/Exceptions.hh>


int
main( int argc, char * argv [] )
{
	try{
  using namespace basic::options;
	using namespace devel::coupled_sidechains;
  using namespace protocols::jd2;

	protocols::jd2::register_options();
	devel::coupled_sidechains::CoupledSidechainProtocol::register_options();

  // initialize core
  devel::init(argc, argv);
  //	core::init::init_random_generators(3,numeric::random::_RND_TestRun_, "mt19937"); //JQX from Sergery

	JobDistributor::get_instance()->go( protocols::moves::MoverOP( new CoupledSidechainProtocol ) );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

