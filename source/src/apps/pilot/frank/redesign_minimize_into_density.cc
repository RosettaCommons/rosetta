// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//core library
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/pack/task/operation/TaskOperations.hh>
#include <core/pack/task/TaskFactory.hh>

#include <basic/options/option.hh>

#include <core/kinematics/MoveMap.hh>

//protocols library (Movers)
#include <protocols/simple_moves/PackRotamersMover.hh>
#include <protocols/simple_moves/symmetry/SymPackRotamersMover.hh>
#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/TaskAwareMinMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/electron_density/util.hh>

#include <core/scoring/electron_density/util.hh>
#include <protocols/simple_moves/symmetry/SymMinMover.hh>

//utilities
#include <protocols/jd2/JobDistributor.hh>
#include <devel/init.hh>

#include <utility/excn/Exceptions.hh>


// option key includes
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

int
main( int argc, char * argv [] )
{
    try {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	devel::init(argc, argv);

	// task factory
	core::pack::task::TaskFactoryOP main_task_factory = new core::pack::task::TaskFactory;
	main_task_factory->push_back( new core::pack::task::operation::InitializeFromCommandline );
	if ( option[ packing::resfile ].user() ) {
		main_task_factory->push_back( new core::pack::task::operation::ReadResfile );
	}

	//create a ScoreFunction from commandline options
	core::scoring::ScoreFunctionOP score_fxn = core::scoring::get_score_function();

	// now add density scores from cmd line
	if ( option[ edensity::mapfile ].user() ) {
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *score_fxn );
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

	// set pose for density scoring if a map was input
	//   + (potentially) dock map into density
	if ( option[ edensity::mapfile ].user() ) {
		seq_mover->add_mover( new protocols::electron_density::SetupForDensityScoringMover );
	}

	seq_mover->add_mover( pack_mover );

	// add constraints to the starting structure
	pose.remove_constraints();
	Size nres = pose.size();
	Real const cst_width( option[ OptionKeys::relax::coord_cst_width ].user()? option[ OptionKeys::relax::coord_cst_width ]() : 0.5 );
	Real const coord_sdev( option[ OptionKeys::relax::coord_cst_width ].user()? option[ OptionKeys::relax::coord_cst_stdev ]() : 0.5 );
	for ( Size ires = 1; ires<= nres - 1; ++ires ) {
		if (!pose.residue(ires).is_protein()) continue;
		Residue const & nat_i_rsd( pose.residue(i) );
		for ( Size iatom = 1; iatom<= nat_i_rsd.last_backbone_atom(); ++iatom ) {
			pose.add_constraint( new CoordinateConstraint(
							  AtomID(iatom,ires), AtomID(1,nres), nat_i_rsd.xyz( iatom ),
							  new BoundFunc( 0, cst_width, coord_sdev, "xyz" )) );
		}
	}

	if ( !option[ symmetry::symmetry_definition ].user() )  {
		core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
		movemap->set_bb(true);
		movemap->set_chi(true);
		movemap->set_jump(true);
		protocols::simple_moves::MinMoverOP min_mover = new protocols::simple_moves::MinMover(
			movemap, score_fxn, "lbfgs_armijo_nonmonotone", 0.01, true );
		seq_mover->add_mover( min_mover );
	} else {
		core::kinematics::MoveMapOP movemap = new core::kinematics::MoveMap;
		movemap->set_bb(true);
		movemap->set_chi(true);
		movemap->set_jump(true);
		protocols::simple_moves::symmetry::SymMinMoverOP min_mover = new protocols::simple_moves::symmetry::SymMinMover(
			movemap, score_fxn, "lbfgs_armijo_nonmonotone", 0.01, true );
		seq_mover->add_mover( min_mover );
	}

	protocols::jd2::JobDistributor::get_instance()->go(seq_mover);
    } catch ( utility::excn::EXCN_Base const & e ) {
                              std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
                                  }
        return 0;
}
