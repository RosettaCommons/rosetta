// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/loop_mover/refine/ShearMinCCDTrial.cc
/// @brief Concrete class derived from LoopRefineInnerCycle to implement the CCD min trial flavor of inner cycle refinement.
/// @details
///
/// @author Michael Pacella (mpacella88@gmail.com)

// Unit headers
#include <protocols/loops/loop_mover/refine/ShearMinCCDTrial.hh>
#include <protocols/loops/loop_mover/refine/ShearMinCCDTrialCreator.hh>

// Package headers

#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

// Project headers
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/kinematics/MoveMap.hh>
#include <protocols/simple_moves/BackboneMover.hh>
#include <core/optimization/Minimizer.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/pose/symmetry/util.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.loops.loop_mover.refine.ShearMinCCDTrial" );
using namespace core;

namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

///////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// BOILER PLATE CODE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

/// @brief default constructor
ShearMinCCDTrial::ShearMinCCDTrial() : LoopRefineInnerCycle()
{
	init();
}


/// @brief copy constructor
ShearMinCCDTrial::ShearMinCCDTrial( ShearMinCCDTrial const & rhs ) : LoopRefineInnerCycle(rhs)
{
	init_for_equal_operator_and_copy_constructor( *this, rhs );
}


/// @brief assignment operator
ShearMinCCDTrial & ShearMinCCDTrial::operator=( ShearMinCCDTrial const & rhs ){
	//abort self-assignment
	if ( this == &rhs ) return *this;
	LoopRefineInnerCycle::operator=( rhs );
	init_for_equal_operator_and_copy_constructor( *this, rhs );
	return *this;
}

//destructor
ShearMinCCDTrial::~ShearMinCCDTrial() {}

/// @brief Each derived class must specify its name.
std::string ShearMinCCDTrial::get_name() const
{
	return type();
}

moves::MoverOP ShearMinCCDTrial::clone() const
{
	return moves::MoverOP( new ShearMinCCDTrial( *this ) );
}

moves::MoverOP ShearMinCCDTrial::fresh_instance() const
{
	return moves::MoverOP( new ShearMinCCDTrial() );
}

/// @brief This mover retains state such that a fresh version is needed if the input Pose is about to change
bool ShearMinCCDTrial::reinitialize_for_new_input() const
{
	return true;
}

//void ShearMinCCDTrial::register_options()
//{
///  PUT THE LIST OF OPTIONS THAT ARE USED HERE  ///

///  RECURSIVELY CALL REGISTER OPTIONS ON ALL MOVERS THAT THIS CLASS HAS AN OWNING_PTR TO  ///
//}
///////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// END OF BOILER PLATE CODE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
ShearMinCCDTrial::ShearMinCCDTrial(
	LoopMover_Refine_CCDAP loop_mover,
	moves::MonteCarloOP mc,
	core::scoring::ScoreFunctionOP scorefxn,
	core::pack::task::TaskFactoryOP tf
) : LoopRefineInnerCycle( loop_mover, mc, scorefxn, tf )
{
	init();
}

void ShearMinCCDTrial::apply( Pose & pose )
{
	setup_objects( pose );

	pack::task::PackerTaskOP task_before_bb_perturbation = task_factory()->create_task_and_apply_taskoperations( pose );
	task_before_bb_perturbation->set_bump_check( true );

	LoopMover_Refine_CCDOP loop_mover_op( loop_mover() ); // lock AP
	Loops::const_iterator it( loop_mover_op->loops()->one_random_loop() );
	Loops one_loop;
	one_loop.add_loop( it );

	// set up movemap properly
	kinematics::MoveMapOP mm_one_loop( new kinematics::MoveMap() );
	loop_mover_op->setup_movemap( pose, one_loop, task_before_bb_perturbation->repacking_residues(), mm_one_loop );

	protocols::simple_moves::ShearMover shear_moves( mm_one_loop, mc()->temperature(), nmoves_ );
	shear_moves.apply( pose );

	if ( ! it->is_terminal( pose ) ) ccd_close_loops( pose, one_loop, *mm_one_loop);

	pack::task::PackerTaskOP task_after_bb_perturbation = task_factory()->create_task_and_apply_taskoperations( pose );
	task_after_bb_perturbation->set_bump_check( true );
	core::pack::rotamer_trials( pose, *scorefxn(), task_after_bb_perturbation );
	(*scorefxn())(pose); // update 10A nbr graph, silly way to do this

	kinematics::MoveMapOP all_loops_movemap = movemap();
	loop_mover_op->setup_movemap( pose, *loop_mover_op->loops(), task_after_bb_perturbation->repacking_residues(), all_loops_movemap );
	if ( loop_mover_op->flank_residue_min() ) {
		add_loop_flank_residues_bb_to_movemap(*loop_mover_op->loops(), *all_loops_movemap);
	} // added by JQX


	minimizer( pose )->run( pose, *all_loops_movemap, *scorefxn(), *min_options_ );
	std::string move_type = "shear_ccd_min";
	mc()->boltzmann( pose, move_type );
	mc()->show_scores();
}


void ShearMinCCDTrial::init()
{
	type( "ShearMinCCDTrial" );
	nmoves_ = 1;
	min_options_ = core::optimization::MinimizerOptionsOP( new core::optimization::MinimizerOptions("lbfgs_armijo_nonmonotone", 0.001, true /*use_nblist*/, false /*deriv_check*/ ) );
	init_options();
}

void ShearMinCCDTrial::init_for_equal_operator_and_copy_constructor(ShearMinCCDTrial & lhs, ShearMinCCDTrial const & rhs)
{
	// copy all data members from rhs to lhs
	lhs.nmoves_ = rhs.nmoves_;
	lhs.min_options_ = rhs.min_options_;
}

void ShearMinCCDTrial::init_options()
{
	//using namespace basic::options;
	// Set options here.
}

core::optimization::AtomTreeMinimizerOP ShearMinCCDTrial::minimizer( core::pose::Pose const & pose ) const
{
	// minimizer
	if ( ! minimizer_ ) {
		if ( core::pose::symmetry::is_symmetric( pose ) ) {
			// minimizer_ = dynamic_cast<core::optimization::AtomTreeMinimizer*>( new core::optimization::symmetry::SymAtomTreeMinimizer );
			minimizer_ = core::optimization::AtomTreeMinimizerOP( new core::optimization::symmetry::SymAtomTreeMinimizer );
		} else {
			minimizer_ = core::optimization::AtomTreeMinimizerOP( new core::optimization::AtomTreeMinimizer );
		}
	}
	return minimizer_;
}
void ShearMinCCDTrial::show( std::ostream & out ) const
{
	out << *this;
}

std::ostream & operator<<(std::ostream& out, ShearMinCCDTrial const & loop_refine_shear_CCD_min_trial_inner_cycle )
{
	out << loop_refine_shear_CCD_min_trial_inner_cycle.get_name() << "Concrete class derived from LoopRefineInnerCycle to implement the CCD min trial flavor of inner cycle refinement." << std::endl;
	return out;
}

void ShearMinCCDTrial::setup_objects( Pose const & pose )
{
	// TR << "Setting up data for " + get_name() + "." << std::endl;

	LoopRefineInnerCycle::setup_objects( pose );
}

ShearMinCCDTrialCreator::~ShearMinCCDTrialCreator() {}

moves::MoverOP ShearMinCCDTrialCreator::create_mover() const {
	return moves::MoverOP( new ShearMinCCDTrial() );
}

std::string ShearMinCCDTrialCreator::keyname() const {
	return "ShearMinCCDTrial";
}

} // namespace refine
} // namespace loop_mover
} // namespace loops
} // namespace protocols
