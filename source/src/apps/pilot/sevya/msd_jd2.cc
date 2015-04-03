// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/apps/pilot/sevya/msd_jd2.cc
/// @brief  Experimental multi state design implementation. Written to encourage convergence in protein designs occurring
/// simultaneously, rather than enforce that they have identical sequences.
/// @author Alex Sevy


//core library


#include <devel/init.hh>
#include <basic/Tracer.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/moves/Mover.hh>
#include <utility/excn/EXCN_Base.hh>
//option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/option.hh>

static basic::Tracer TR("msd_jd2.main");

int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		devel::init( argc, argv );
		//Make JobDistributorFactory create a MSDJobDistributor instance
		option[ OptionKeys::run::msd_job_dist ].value( true );
		//Create dummy mover for job distributor
		protocols::moves::MoverOP mover;
		protocols::jd2::JobDistributor::get_instance()->go( mover );
		return 1;
	} catch ( utility::excn::EXCN_Base const & e ) {
		utility_exit_with_message("caught exception " + e.msg());
	}
}
