// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file 		src/apps/pilot/ralford/pymol_viewer_test.cc
///
/// @brief 		Pilot Application for Testing the PyMolViewer Additions
/// @details	Add a membrane to the pose and turn on the capability to view the membranes
///				in pymol as defined by the normal/center pair defined in the framework. 
///
/// @author 	Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified: 7/22/14
/// thanks to Andrew's strand_shifter app for some quick repack/min code

// Unit Headers
#include <devel/init.hh>

// Package Headers
#include <protocols/moves/Mover.hh>

#include <protocols/jd2/JobDistributor.hh> 
#include <protocols/jd2/util.hh>

#include <protocols/membrane/AddMembraneMover.hh> 
#include <protocols/membrane/InitialMembranePositionMover.hh> 

// Project Headers
#include <core/scoring/ScoreFunction.hh> 
#include <core/scoring/ScoreFunctionFactory.hh> 

#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/FoldTree.hh> 

#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>


#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/PackRotamersMover.hh>

#include <core/pose/Pose.hh> 
#include <core/types.hh> 

// Utility Headers
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>

#include <utility/pointer/owning_ptr.hh> 

#include <basic/Tracer.hh> 

using namespace protocols::moves;

static basic::Tracer TR( "apps.pilot.ralford.pymol_viewer_test" );

/// @brief View membrane planes based on normal/center
class MembraneViewMover : public Mover {

public: 

	MembraneViewMover() {}

	/// @brief Get Mover Name
	std::string get_name() const { return "MembraneViewMover"; }

	/// @brief Apply Membrane Relax
	void apply( Pose & pose ) {

		using namespace core::scoring; 
		using namespace protocols::membrane; 
		using namespace protocols::simple_moves; 
		using namespace core::pack::task; 
		using namespace core::optimization; 
		using namespace core::scoring; 
		using namespace core::kinematics;
		// Add a Membrane
		AddMembraneMoverOP add_memb = new AddMembraneMover(); 
		add_memb->apply( pose ); 

		// Set the Initial Position of the Membrane
		InitialMembranePositionMoverOP init_memb = new InitialMembranePositionMover();
		init_memb->apply( pose ); 

		ScoreFunctionOP sfxn = ScoreFunctionFactory::create_score_function("fa_menv_smooth_2014");

		TR << "Setting up a new movemap" << std::endl;
		// Setup a MoveMap Object
		MoveMapOP mm = new MoveMap(); 
		mm->set_bb( false ); 
		mm->set_chi( false ); 
		mm->set_jump( true ); 

		MinMoverOP stage1_min = new MinMover( mm, sfxn, "lbfgs_armijo_nonmonotone", 0.1, false );
		MinMoverOP stage2_min = new MinMover( mm, sfxn, "lbfgs_armijo_nonmonotone", 0.01, false );
		MinMoverOP stage3_min = new MinMover( mm, sfxn, "lbfgs_armijo_nonmonotone", 0.001, false );

		stage1_min->apply( pose ); 
		stage2_min->apply( pose );  
		stage3_min->apply( pose ); 

		TR << "Score: " << (*sfxn)(pose) << std::endl;

		//for ( Size i = 1; i <= 1; ++i ) {
		//	new_fastrelax_cycle( pose, sfxn ); 	
		//}

	}

	void
	new_fastrelax_cycle( Pose & pose, core::scoring::ScoreFunctionOP sfxn ) {

		using namespace protocols::simple_moves; 
		using namespace core::pack::task; 
		using namespace core::optimization; 
		using namespace core::scoring; 
		using namespace core::kinematics;

		TR << "Setting up packer task" << std::endl;
		PackerTaskOP repack_task = core::pack::task::TaskFactory::create_packer_task( pose ); 
		repack_task->restrict_to_repacking(); 
		repack_task->or_include_current( true ); 

		TR << "Setting up a new movemap" << std::endl;
		// Setup a MoveMap Object
		MoveMapOP mm = new MoveMap(); 
		mm->set_bb_true_range( 1, 183 ); 
		mm->set_chi( true ); 
		mm->set_jump( true ); 

		// Reweight fa_rep to 0.2, repack, and minimize @ 0.1 tol
		TR << "Running cycle 1 - low repulsive weight, pretty good tolerance..." << std::endl;
		sfxn->set_weight( fa_rep, 0.2 ); 
		PackRotamersMoverOP stage1_pack = new PackRotamersMover( sfxn, repack_task );
		MinMoverOP stage1_min = new MinMover( mm, sfxn, "lbfgs_armijo_nonmonotone", 0.1, false ); 
		stage1_pack->apply( pose ); 
		stage1_min->apply( pose ); 

		// Reweight fa_rep to 0.2, repack, and minimize @ 0.01 tol
		sfxn->set_weight( fa_rep, 0.25 ); 
		PackRotamersMoverOP stage2_pack = new PackRotamersMover( sfxn, repack_task );
		MinMoverOP stage2_min = new MinMover( mm, sfxn, "lbfgs_armijo_nonmonotone", 0.01, false ); 
		stage2_pack->apply( pose ); 
		stage2_min->apply( pose ); 

		// Reweight fa_rep to 0.2, repack, and minimize @ 0.001 tol
		sfxn->set_weight( fa_rep, 0.55 ); 
		PackRotamersMoverOP stage3_pack = new PackRotamersMover( sfxn, repack_task );
		MinMoverOP stage3_min = new MinMover( mm, sfxn, "lbfgs_armijo_nonmonotone", 0.001, false ); 
		stage3_pack->apply( pose ); 
		stage3_min->apply( pose ); 

		// Reweight fa_rep to 0.2, repack, and minimize @ 0.0001 tol
		sfxn->set_weight( fa_rep, 1.0 ); 
		PackRotamersMoverOP stage4_pack = new PackRotamersMover( sfxn, repack_task );
		MinMoverOP stage4_min = new MinMover( mm, sfxn, "lbfgs_armijo_nonmonotone", 0.0001, false ); 
		stage4_pack->apply( pose ); 
		stage4_min->apply( pose );

	}
};

typedef utility::pointer::owning_ptr< MembraneViewMover > MembraneViewMoverOP;
typedef utility::pointer::owning_ptr< MembraneViewMover const > MembraneViewMoverCOP; 

/// @brief Main method
int
main( int argc, char * argv [] )
{
	try {

		// Devel init factories
		devel::init(argc, argv);

		// Register JD2 options
		protocols::jd2::register_options();

		// Minimize with mp
		MembraneViewMoverOP mpview = new MembraneViewMover();
		protocols::jd2::JobDistributor::get_instance()->go( mpview );

		return 0; 

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}
