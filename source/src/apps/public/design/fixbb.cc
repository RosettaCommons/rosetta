// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file /src/apps/public/design/fixbb.cc
/// @brief  Fixed backbone design.  Can do side chain minimization after PackRotamers by using the flag -minimize_sidechains.  This is SLOW.

//core library
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// TEMP
#include <core/scoring/hbonds/HBondOptions.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>

#include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/util.hh>

#include <core/kinematics/MoveMap.hh>

//protocols library (Movers)
#include <protocols/simple_moves/MinPackMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>

//utilities
#include <protocols/jd2/JobDistributor.hh>
#include <devel/init.hh>

// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>


//local options
namespace basic{ namespace options{ namespace OptionKeys{
basic::options::BooleanOptionKey const minimize_sidechains("minimize_sidechains");
basic::options::BooleanOptionKey const min_pack("min_pack");
basic::options::BooleanOptionKey const stochastic_pack("stochastic_pack");
}}}//basic::options::OptionKeys

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	option.add( minimize_sidechains, "Do minimization of side chains after rotamer packing").def(false);
	option.add( min_pack, "Pack and minimize sidechains simultaneously").def(false);
	option.add( stochastic_pack, "Pack using a continuous sidechains rotamer library").def(false);

	devel::init(argc, argv);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// end of setup
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//create a task factory: this will create a new PackerTask for each input pose
	core::pack::task::TaskFactoryOP main_task_factory = new core::pack::task::TaskFactory;
	main_task_factory->push_back( new core::pack::task::operation::InitializeFromCommandline );

	/// As of 2010/07/16, the ReadResfile operation is a no-op unless a resfile has been
	/// supplied on the command line, through the ResourceManager, or programmatically.
	/// Therefore, it is safe to add it without first checking to see if something has been
	/// provided on the command line.  If that check were here, then a resfile provided
	/// through the ResourceManager would not get read.
	main_task_factory->push_back( new core::pack::task::operation::ReadResfile );

	//create a ScoreFunction from commandline options
	core::scoring::ScoreFunctionOP score_fxn = core::scoring::getScoreFunction();

	/// TEMP!
	if ( false ) {
		using namespace core;
		scoring::methods::EnergyMethodOptionsOP emopts( new scoring::methods::EnergyMethodOptions( score_fxn->energy_method_options() ));
		emopts->hbond_options().use_hb_env_dep( false );
		emopts->hbond_options().decompose_bb_hb_into_pair_energies( true );
		score_fxn->set_energy_method_options( *emopts );
		score_fxn->set_weight( scoring::fa_pair, 0.0 );
	}


	//create the PackRotamersMover which will do the packing
	protocols::simple_moves::PackRotamersMoverOP pack_mover = new protocols::simple_moves::PackRotamersMover;

	// Use the symmetric packer if necessary
	if ( option[ symmetry::symmetry_definition ].user() ) {
		pack_mover = new protocols::simple_moves::symmetry::SymPackRotamersMover;
	}

	pack_mover->task_factory( main_task_factory );
	pack_mover->score_function( score_fxn );

	//This sequence mover will contain packing for sure, and may contain minimization
	protocols::moves::SequenceMoverOP seq_mover = new protocols::moves::SequenceMover;

	// make symmetric pose if necessary
	if ( option[ symmetry::symmetry_definition ].user() )  {
	    seq_mover->add_mover( new protocols::simple_moves::symmetry::SetupForSymmetryMover );
	}

	if ( option[ min_pack ] || option[ stochastic_pack ] ) {
		protocols::simple_moves::MinPackMoverOP minpack_mover = new protocols::simple_moves::MinPackMover;
		minpack_mover->task_factory( main_task_factory );
		minpack_mover->score_function( score_fxn );
		if ( option[ stochastic_pack ] ) minpack_mover->stochastic_pack( true );
		seq_mover->add_mover( minpack_mover );
	} else {
		seq_mover->add_mover( pack_mover );
	}

	//If sidechain minimization is requested, include that too
	if ( option[ minimize_sidechains ] ) {
            core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
            protocols::simple_moves::MinMoverOP min_mover;
            if (option [ symmetry::symmetry_definition ].user() ){
                min_mover = new protocols::simple_moves::symmetry::SymMinMover(
			movemap,
			score_fxn,
			basic::options::option[ basic::options::OptionKeys::run::min_type ].value(),
			0.01,
			true
		);
            }
            else {
		min_mover = new protocols::simple_moves::MinMover(
			movemap,
			score_fxn,
			basic::options::option[ basic::options::OptionKeys::run::min_type ].value(),
			0.01,
			true
		);
            }
	    protocols::simple_moves::TaskAwareMinMoverOP TAmin_mover = new protocols::simple_moves::TaskAwareMinMover(min_mover, main_task_factory);
            seq_mover->add_mover( TAmin_mover );
	} // end optional side chain minimization

	protocols::jd2::JobDistributor::get_instance()->go(seq_mover);
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	} 
}
