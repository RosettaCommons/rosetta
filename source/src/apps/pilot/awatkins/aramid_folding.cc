// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Rhiju Das

#include <utility/json_utilities.hh>

#if defined(ZEROMQ)  and  defined(_NLOHMANN_JSON_ENABLED_)
#include <protocols/network/hal.hh>
#include <protocols/network/util.hh>
#include <json.hpp>
#include <protocols/network/ui_mover.hh>
#include <core/import_pose/import_pose.hh>
#endif

// libRosetta headers
#include <core/types.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <core/pose/Pose.hh>
#include <devel/init.hh>
#include <utility/vector1.hh>

#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/id/TorsionID.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/util.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperation.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>

#include <protocols/moves/TrialMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/RepeatMover.hh>
#include <protocols/minimization_packing/MinMover.hh>
#include <protocols/minimization_packing/RotamerTrialsMinMover.hh>
#include <protocols/simple_moves/RandomTorsionMover.hh>

//#include <protocols/moves/PyMOLMover.hh>

// C++ headers
#include <iostream>
#include <string>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/pose/annotated_sequence.hh>
#include <utility/excn/Exceptions.hh>

#include <utility/options/keys/OptionKeyList.hh>

using namespace core;
using namespace core::chemical;
using namespace core::scoring;
using namespace core::pose;
using namespace protocols;
using namespace basic::options::OptionKeys;
using utility::vector1;

#if defined(ZEROMQ)  and  defined(_NLOHMANN_JSON_ENABLED_)
using namespace utility;
using namespace protocols::network;
#endif

static basic::Tracer TR( "apps.pilot.awatkins.aramid_folding" );


///////////////////////////////////////////////////////////////////////////////
#if defined(ZEROMQ)  and  defined(_NLOHMANN_JSON_ENABLED_)
core::pose::Pose
#else
void
#endif
aramid_main()
{
	using namespace basic::options;
	using namespace core::pose;
	using namespace protocols::simple_moves;
	using namespace core::id;

	Pose pose;

	ScoreFunctionOP sfxn = core::scoring::get_score_function();

	core::chemical::ResidueTypeSetCOP rts = core::chemical::ChemicalManager::get_instance()->residue_type_set( chemical::FA_STANDARD );

	// set up by sequence
	make_pose_from_sequence( pose, "X[POLYARAMID_GLU]X[POLYARAMID_SER]X[POLYARAMID_GLU]X[POLYARAMID_SER]X[POLYARAMID_PHE]X[POLYARAMID_MET]X[POLYARAMID_MET]X[POLYARAMID_PHE]X[POLYARAMID_PHE]X[POLYARAMID_MET]X[POLYARAMID_PHE]X[POLYARAMID_PHE]X[POLYARAMID_PHE]X[POLYARAMID_PHE]X[POLYARAMID_GLU]X[POLYARAMID_SER]X[POLYARAMID_GLU]X[POLYARAMID_SER]", rts );
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		pose.set_torsion( TorsionID( ii, core::id::BB, 1 ), 180 );
		pose.set_torsion( TorsionID( ii, core::id::BB, 4 ), 180 );
		pose.set_torsion( TorsionID( ii, core::id::BB, 5 ), 180 );
	}

#if defined(ZEROMQ)  and  defined(_NLOHMANN_JSON_ENABLED_)
	protocols::network::AddUIObserver( pose );
#endif

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 600, 600 );

	core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
	movemap->set_chi( true );
	for ( Size ii = 1; ii <= pose.size(); ++ii ) {
		movemap->set_bb( ii, 1, true );
		movemap->set_bb( ii, 2, false );
		movemap->set_bb( ii, 3, false );
		movemap->set_bb( ii, 4, true );
		movemap->set_bb( ii, 5, false );
	}

	using namespace core::pack::task;
	using namespace core::pack::task::operation;
	TaskFactoryOP pert_tf( new TaskFactory() );
	pert_tf->push_back( TaskOperationCOP( new core::pack::task::operation::InitializeFromCommandline ) );
	operation::RestrictToRepackingOP pert_rtrp( new operation::RestrictToRepacking() );
	pert_tf->push_back( pert_rtrp );

	protocols::minimization_packing::MinMoverOP minmover( new protocols::minimization_packing::MinMover( movemap, sfxn, "lbfgs_armijo_nonmonotone", 0.001, true ) );
	protocols::minimization_packing::RotamerTrialsMinMoverOP rtminmover( new protocols::minimization_packing::RotamerTrialsMinMover( sfxn, pert_tf ) );

	// sample pose
	//
	moves::MonteCarloOP mc( new moves::MonteCarlo( pose, *sfxn, 10.0 ) );

	RandomTorsionMoverOP rtm( new RandomTorsionMover( movemap, 60, 2*pose.size() ) );
	moves::SequenceMoverOP inner_seq( new moves::SequenceMover( rtm, rtminmover ) );
	moves::RepeatMoverOP repeat_seq( new moves::RepeatMover( inner_seq, 10 ) );
	moves::SequenceMoverOP full_seq( new moves::SequenceMover( repeat_seq, minmover ) );

	moves::TrialMoverOP pert_trial( new moves::TrialMover( full_seq, mc ) );

	for ( Size jj = 1; jj <= 10; ++jj ) {
		TR << "Stage " << jj << "." << std::endl;
		for ( Size ii = 1; ii <= 1000; ++ii ) {
			pert_trial->apply( pose );
			if ( ii % 100 == 0 ) {
				mc->show_counters();
			}
		}
		mc->reset_counters();
		mc->set_temperature( 11 - jj );
	}

	pose.dump_pdb( "done.pdb" );



#if defined(ZEROMQ)  and  defined(_NLOHMANN_JSON_ENABLED_)
	return pose;
#endif
}

#if defined(ZEROMQ)  and  defined(_NLOHMANN_JSON_ENABLED_)
///////////////////////////////////////////////////////////////
core::pose::Pose
my_main()
{
	using namespace basic::options;

	return aramid_main();
	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}

#else
///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	using namespace basic::options;

	aramid_main();
	protocols::viewer::clear_conformation_viewers();
	exit( 0 );
}
#endif


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;

		std::cout << std::endl << "Basic usage:  " << argv[0] << "  -fasta <fasta file with sequence>  [ -native <native pdb file> ] " << std::endl;
		std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

		option.add_relevant( in::file::fasta );
		option.add_relevant( in::file::native );
		option.add_relevant( in::file::input_res );
		option.add_relevant( out::file::silent );
		option.add_relevant( out::nstruct );
		option.add_relevant( out::overwrite );

		////////////////////////////////////////////////////////////////////////////
		// setup
		////////////////////////////////////////////////////////////////////////////
		devel::init(argc, argv);

		////////////////////////////////////////////////////////////////////////////
		// end of setup
		////////////////////////////////////////////////////////////////////////////
#if defined(ZEROMQ)  and  defined(_NLOHMANN_JSON_ENABLED_)
		{ // creating dummy pose object to trigger database load so later we can create Pose immeditaly
			core::pose::Pose p;
			core::import_pose::pose_from_pdbstring(p, "ATOM     17  N   ILE A   1      16.327  47.509  23.466  1.00  0.00\n");
		}

		//protocols::network::hal(specification, hal_executioner, protocols::network::CommandLineArguments{argc, argv} );

		hal({
			{"main", {my_main, {
			}
			} } },
			protocols::network::CommandLineArguments{argc, argv} );

		return 0;


#else
		protocols::viewer::viewer_main( my_main );
#endif

	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

