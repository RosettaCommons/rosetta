// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/pilot/doug/cyclize_peptoid_peptide/Cyclization.cc
/// @brief Cyclization application
/// @author P. Douglas Renfrew ( renfrew@nyu.edu )


// protocols header
#include <protocols/jd2/JobDistributor.hh>

#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/PyMolMover.hh>

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/ScoreMover.hh>
#include <protocols/simple_moves/CyclizationMover.hh>
#include <protocols/simple_moves/RandomTorsionMover.hh>

// core headers
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/kinematics/MoveMap.hh>

// devel headers
#include <devel/init.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/cyclization.OptionKeys.gen.hh>

static thread_local basic::Tracer TR( "Cyclization" );

int
main( int argc, char * argv [] )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys::cyclization;
	using namespace protocols::simple_moves;

	// init
  devel::init( argc, argv );

	// setup the cyclization mover
	CyclizationMoverOP cyc_mvr( new CyclizationMover( option[chain_to_cyclize].value(), true, true, option[num_min_rebuild].value() ) );

	// setup the move map
	core::kinematics::MoveMapOP move_map( new core::kinematics::MoveMap() );
	move_map->set_bb( true );
	move_map->set_chi( true );

	// setup the random torsion mover
	RandomTorsionMoverOP ran_tor_mvr( new RandomTorsionMover( move_map, 5.0, 100 ) );

	// setup the sequnce mover
	protocols::moves::SequenceMoverOP seq_mvr( new protocols::moves::SequenceMover() );
	seq_mvr->add_mover( cyc_mvr );
	seq_mvr->add_mover( ran_tor_mvr );
	seq_mvr->add_mover( cyc_mvr );

	/*

		WRITE going to need a rebuild at residue connection mover (ctor with position)
		WRITE going to need a omega flipper mover that resists down chain propagation

		CREATE a sequence mover
		CREATE a rotamer trails mover
		CREATE a pack rotamers mover
		CREATE a min mover (or use the cyclization mover)
		CREATE a monte carlo mover

		ADD copper constraints

	*/

	// go go go
	protocols::jd2::JobDistributor::get_instance()->go( seq_mvr );

  TR << "\n+-----------------------------------------------------------------+\n"
     <<   "|                              DONE                               |\n"
     <<   "+-----------------------------------------------------------------+" << std::endl;

  return 0;
}
