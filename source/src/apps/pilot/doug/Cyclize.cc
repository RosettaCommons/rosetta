// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file /src/apps/pilot/doug/cyclize_peptoid_peptide/Cyclization.cc
/// @brief Cyclization application
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )

// unit headers
#include <protocols/simple_moves/CyclizationMover.hh>

// protocols header
#include <protocols/jd2/JobDistributor.hh>

#include <protocols/moves/MoverContainer.hh>

// devel headers
#include <devel/init.hh>

// core headers


// basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/cyclization.OptionKeys.gen.hh>

// utility headers
#include <utility/excn/Exceptions.hh>

static basic::Tracer TR( "Cyclization" );

int
main( int argc, char * argv [] )
{
	try {

		using namespace basic::options;
		using namespace basic::options::OptionKeys::cyclization;
		using namespace protocols::simple_moves;
		using namespace protocols::moves;

		// init
		devel::init( argc, argv );

		// setup sequence mover
		SequenceMoverOP sm( new SequenceMover() );

		// setup the cyclization mover(s)
		if ( option[chains_to_cyclize].user() ) {
			core::Size num_cyclic_chains( option[chains_to_cyclize].value().size() );

			for ( core::Size i(1); i <= num_cyclic_chains; ++i ) {
				sm->add_mover( utility::pointer::make_shared< CyclizationMover >( option[chains_to_cyclize].value()[i], option[add_constraints].value(), true, option[num_min_rebuild].value() ) );
			}
		}

		// go go go
		protocols::jd2::JobDistributor::get_instance()->go( sm );

		TR << "\n+-----------------------------------------------------------------+\n"
			<<   "|                              DONE                               |\n"
			<<   "+-----------------------------------------------------------------+" << std::endl;

	} catch (utility::excn::Exception const & e ) {
		e.display();
		return -1;
	}

	return 0;
}
