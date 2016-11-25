// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loops/loop_mover/refine/LoopRefineInnerCycle.cc
/// @brief Perform a small move followed CCD closure, packing and minimization
/// @details
///
/// @author Brian D. Weitzner ( brian.weitzner@gmail.com )

// Unit headers
#include <protocols/loops/loop_mover/refine/SmallMinCCDTrial.hh>
#include <protocols/loops/loop_mover/refine/SmallMinCCDTrialCreator.hh>

// Package headers
#include <protocols/loops/loops_main.hh>
#include <protocols/loops/loop_mover/LoopMover.hh>
#include <protocols/loops/loop_mover/refine/LoopMover_CCD.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

// Project headers
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/symmetry/SymAtomTreeMinimizer.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pose/Pose.hh>
#include <core/pose/symmetry/util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/simple_moves/BackboneMover.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// ObjexxFCL headers
#include <ObjexxFCL/format.hh>

// C++ headers
#include <fstream>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.loops.loop_mover.refine.SmallMinCCDTrial" );

namespace protocols {
namespace loops {
namespace loop_mover {
namespace refine {

///////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// BOILER PLATE CODE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

/// @brief default constructor
SmallMinCCDTrial::SmallMinCCDTrial() : LoopRefineInnerCycle()
{
	init();
}

/// @brief copy constructor
SmallMinCCDTrial::SmallMinCCDTrial( SmallMinCCDTrial const & rhs ) : LoopRefineInnerCycle(rhs)
{
	init_for_equal_operator_and_copy_constructor( *this, rhs );
}

/// @brief assignment operator
SmallMinCCDTrial & SmallMinCCDTrial::operator=( SmallMinCCDTrial const & rhs ){
	//abort self-assignment
	if ( this == &rhs ) return *this;
	LoopRefineInnerCycle::operator=( rhs );
	init_for_equal_operator_and_copy_constructor( *this, rhs );
	return *this;
}

//destructor
SmallMinCCDTrial::~SmallMinCCDTrial() {}

/// @brief Each derived class must specify its name.
// XRW TEMP std::string SmallMinCCDTrial::get_name() const
// XRW TEMP {
// XRW TEMP  return type();
// XRW TEMP }

//@brief clone operator, calls the copy constructor
protocols::moves::MoverOP
SmallMinCCDTrial::clone() const
{
	return protocols::moves::MoverOP( new SmallMinCCDTrial( *this ) );
}

/// @brief fresh_instance returns a default-constructed object for JD2
protocols::moves::MoverOP
SmallMinCCDTrial::fresh_instance() const
{
	return protocols::moves::MoverOP( new SmallMinCCDTrial() );
}

/// @brief This mover retains state such that a fresh version is needed if the input Pose is about to change
bool SmallMinCCDTrial::reinitialize_for_new_input() const
{
	return true;
}

void SmallMinCCDTrial::register_options()
{
	///  PUT THE LIST OF OPTIONS THAT ARE USED HERE  ///

	///  RECURSIVELY CALL REGISTER OPTIONS ON ALL MOVERS THAT THIS CLASS HAS AN OWNING_PTR TO  ///
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// END OF BOILER PLATE CODE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

// constructor with arguments
SmallMinCCDTrial::SmallMinCCDTrial(
	LoopMover_Refine_CCDAP loop_mover,
	moves::MonteCarloOP mc,
	core::scoring::ScoreFunctionOP scorefxn,
	core::pack::task::TaskFactoryOP tf
) : LoopRefineInnerCycle( loop_mover, mc, scorefxn, tf )
{
	init();
}

void SmallMinCCDTrial::apply( Pose & pose )
{
	using core::kinematics::MoveMapOP;
	using core::pack::rotamer_trials;
	using core::pack::task::PackerTaskOP;

	// TR << "Beginning apply function of " + get_name() + "." << std::endl;

	setup_objects( pose );

	// show( TR );

	// TODO: Determine if set_bumb_check can be done at the TaskFactory level and do it there if possible.
	// TODO: If not, make these two lines a function so I can't mess this up.
	PackerTaskOP task_before_bb_perturbation = task_factory()->create_task_and_apply_taskoperations( pose );
	task_before_bb_perturbation->set_bump_check( true );

	LoopMover_Refine_CCDOP loop_mover_op( loop_mover() ); // lock AP
	LoopsOP all_loops = loop_mover_op->loops();
	Loops one_loop = get_one_random_loop();

	// set up movemap for one loop
	MoveMapOP one_loop_movemap = movemap();
	loop_mover_op->setup_movemap( pose, one_loop, task_before_bb_perturbation->repacking_residues(), one_loop_movemap );

	debug_zero( pose );

	simple_moves::SmallMover small_moves( one_loop_movemap, mc()->temperature(), nmoves_ );
	small_moves.apply( pose );

	debug_one( pose );


	if ( ! one_loop[ one_loop.size() ].is_terminal( pose ) ) ccd_close_loops( pose, one_loop, *one_loop_movemap);

	debug_two( pose );

	// TODO: Determine if set_bumb_check can be done at the TaskFactory level and do it there if possible.
	// TODO: If not, make these two lines a function so I can't mess this up.
	PackerTaskOP task_after_bb_perturbation = task_factory()->create_task_and_apply_taskoperations( pose );
	task_after_bb_perturbation->set_bump_check( true );
	rotamer_trials( pose, *scorefxn(), task_after_bb_perturbation );
	debug_three( pose );


	MoveMapOP all_loops_movemap = movemap();
	loop_mover_op->setup_movemap( pose, *all_loops, task_after_bb_perturbation->repacking_residues(), all_loops_movemap );

	if ( loop_mover_op->flank_residue_min() ) {
		add_loop_flank_residues_bb_to_movemap(*all_loops, *all_loops_movemap);
	}

	minimizer( pose )->run( pose, *all_loops_movemap, *scorefxn(), *minimizer_options_ );

	debug_four( pose );

	std::string move_type = "small_ccd_min";
	mc()->boltzmann( pose, move_type );

	debug_five( pose );

	mc()->show_scores();
}


void SmallMinCCDTrial::setup_objects( Pose const & pose )
{
	// TR << "Setting up data for " + get_name() + "." << std::endl;

	LoopRefineInnerCycle::setup_objects( pose );
}

void SmallMinCCDTrial::init()
{
	using core::optimization::MinimizerOptions;
	type( "SmallMinCCDTrial" );

	nmoves_ = 1;
	minimizer_options_ = core::optimization::MinimizerOptionsOP( new MinimizerOptions( "lbfgs_armijo_nonmonotone", 0.001, true /*use_nblist*/, false /*deriv_check*/ ) );
	init_options();
}

void SmallMinCCDTrial::init_for_equal_operator_and_copy_constructor(
	SmallMinCCDTrial & lhs,
	SmallMinCCDTrial const & rhs
)
{
	// copy all data members from rhs to lhs
	lhs.minimizer_ = rhs.minimizer_;
}

void SmallMinCCDTrial::init_options()
{
	/* UNCOMMENT WHEN THERE ARE ACTUALLY OPTIONS TO PROCESS
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	*/
	// Set options here.
}

core::optimization::AtomTreeMinimizerOP SmallMinCCDTrial::minimizer( core::pose::Pose const & pose ) const
{
	// minimizer
	if ( ! minimizer_ ) {
		if ( core::pose::symmetry::is_symmetric( pose ) ) {
			// minimizer_ = dynamic_cast<core::optimization::AtomTreeMinimizer>( new core::optimization::symmetry::SymAtomTreeMinimizer );
			minimizer_ = core::optimization::AtomTreeMinimizerOP( new core::optimization::symmetry::SymAtomTreeMinimizer );
		} else {
			minimizer_ = core::optimization::AtomTreeMinimizerOP( new core::optimization::AtomTreeMinimizer );
		}
	}
	return minimizer_;
}

core::Size SmallMinCCDTrial::number_of_moves() const
{
	return nmoves_;
}

void SmallMinCCDTrial::set_number_of_moves( core::Size nmoves )
{
	nmoves_ = nmoves;
}

core::optimization::MinimizerOptionsOP SmallMinCCDTrial::minimizer_options() const
{
	return minimizer_options_;
}

void SmallMinCCDTrial::set_minimizer_options( core::optimization::MinimizerOptionsOP minimizer_options )
{
	if ( minimizer_options ) minimizer_options_ = minimizer_options;
}

void
SmallMinCCDTrial::show( std::ostream & out ) const
{
	out << *this;
}

std::ostream & operator<<(std::ostream& out, SmallMinCCDTrial const & small_min_ccd_trial )
{
	out << small_min_ccd_trial.get_name() << " is an awesome class." << std::endl;
	return out;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// EXCESSIVE DEBUGGING OUTPUT ////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////
void SmallMinCCDTrial::debug_zero( Pose & pose )
{
	if ( debug() ) {
		LoopMover_Refine_CCDOP loop_mover_op( loop_mover() ); // lock AP
		TR << "chutmp-debug small_move-0: " << "  " << (*scorefxn())(pose) << std::endl;
		TR << "small_move-0: " << pose.energies().total_energies().weighted_string_of( scorefxn()->weights() )
			<< " rmsd: " << ObjexxFCL::format::F(9,3,loop_rmsd( pose, *get_native_pose(), *loop_mover_op->loops() )) << std::endl;
		pose.dump_pdb("small_move-0.pdb");
	}

}

void SmallMinCCDTrial::debug_one( Pose & pose )
{
	if ( debug() ) {
		LoopMover_Refine_CCDOP loop_mover_op( loop_mover() ); // lock AP
		TR << "chutmp-debug small_move-1: " << "  " << (*scorefxn())(pose) << std::endl;
		TR << "small_move-1: " << pose.energies().total_energies().weighted_string_of( scorefxn()->weights() )
			<< " rmsd: " << ObjexxFCL::format::F(9,3,loop_rmsd( pose, *get_native_pose(), *(loop_mover_op->loops()) )) << std::endl;
		pose.dump_pdb("small_move-1.pdb");
		std::ofstream out("score.small_move_1");
		out << "scoring of input_pose " << (*scorefxn())(pose) << std::endl;
		scorefxn()->show( out );
		out << pose.energies();
	}
}

void SmallMinCCDTrial::debug_two( Pose & pose )
{
	if ( debug() ) {
		LoopMover_Refine_CCDOP loop_mover_op( loop_mover() ); // lock AP
		TR << "chutmp-debug small_move-2: " << "  " << (*scorefxn())(pose) << std::endl;
		TR << "small_move-2: " << pose.energies().total_energies().weighted_string_of( scorefxn()->weights() )
			<< " rmsd: " << ObjexxFCL::format::F(9,3,loop_rmsd( pose, *get_native_pose(), *(loop_mover_op->loops()) )) << std::endl;
		pose.dump_pdb("small_move-2.pdb");
		std::ofstream out("score.small_move_2");
		out << "scoring of input_pose " << (*scorefxn())(pose) << std::endl;
		scorefxn()->show( out );
		out << pose.energies();
	}
}

void SmallMinCCDTrial::debug_three( Pose & pose )
{
	if ( debug() ) {
		LoopMover_Refine_CCDOP loop_mover_op( loop_mover() ); // lock AP
		TR << "chutmp-debug small_move-3: " << "  " << (*scorefxn())(pose) << std::endl;
		TR << "small_move-3: " << pose.energies().total_energies().weighted_string_of( scorefxn()->weights() )
			<< " rmsd: " << ObjexxFCL::format::F(9,3,loop_rmsd( pose, *get_native_pose(), *(loop_mover_op->loops()) )) << std::endl;
		pose.dump_pdb("small_move-3.pdb");
	}
}

void SmallMinCCDTrial::debug_four( Pose & pose )
{
	if ( debug() ) {
		LoopMover_Refine_CCDOP loop_mover_op( loop_mover() ); // lock AP
		TR << "chutmp-debug small_move-4: " << "  " << (*scorefxn())(pose) << std::endl;
		TR << "small_move-4: " << pose.energies().total_energies().weighted_string_of( scorefxn()->weights() )
			<< " rmsd: " << ObjexxFCL::format::F(9,3,loop_rmsd( pose, *get_native_pose(), *(loop_mover_op->loops()) )) << std::endl;
		pose.dump_pdb("small_move-4.pdb");
	}
}

void SmallMinCCDTrial::debug_five( Pose & pose )
{
	if ( debug() ) {
		LoopMover_Refine_CCDOP loop_mover_op( loop_mover() ); // lock AP
		TR << "chutmp-debug small_move-5: " << "  " << (*scorefxn())(pose) << std::endl;
		TR << "small_move-5: " << pose.energies().total_energies().weighted_string_of( scorefxn()->weights() )
			<< " rmsd: " << ObjexxFCL::format::F(9,3,loop_rmsd( pose, *get_native_pose(), *(loop_mover_op->loops()) )) << std::endl;
		pose.dump_pdb("small_move-5.pdb");
	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////// END OF EXCESSIVE DEBUG OUTPUT //////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

// XRW TEMP SmallMinCCDTrialCreator::~SmallMinCCDTrialCreator() {}

// XRW TEMP moves::MoverOP SmallMinCCDTrialCreator::create_mover() const {
// XRW TEMP  return moves::MoverOP( new SmallMinCCDTrial() );
// XRW TEMP }

// XRW TEMP std::string SmallMinCCDTrialCreator::keyname() const {
// XRW TEMP  return "SmallMinCCDTrial";
// XRW TEMP }

std::string SmallMinCCDTrial::get_name() const {
	return mover_name();
}

std::string SmallMinCCDTrial::mover_name() {
	return "SmallMinCCDTrial";
}

void SmallMinCCDTrial::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	protocols::moves::xsd_type_definition_w_attributes( xsd, mover_name(),
		"Perform a small move followed CCD closure, packing and minimization. Takes only command line options.", attlist );
}

std::string SmallMinCCDTrialCreator::keyname() const {
	return SmallMinCCDTrial::mover_name();
}

protocols::moves::MoverOP
SmallMinCCDTrialCreator::create_mover() const {
	return protocols::moves::MoverOP( new SmallMinCCDTrial );
}

void SmallMinCCDTrialCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SmallMinCCDTrial::provide_xml_schema( xsd );
}


} // namespace refine
} // namespace loop_mover
} // namespace loops
} // namespace protocols
