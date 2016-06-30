// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/antibody/SnugDock.cc
/// @brief Dock and antigen to an antibody while optimizing the rigid body orientation of the VH and VL chains and performing CDR loop minimization.
/// @details
///
///
/// @author Jianqing Xu ( xubest@gmail.com )
/// @author Brian D. Weitzner ( brian.weitzner@gmail.com )
/// @author Jeliazko Jeliazkov ( jeliazkov@jhu.edu )

// Unit headers
#include <protocols/antibody/snugdock/SnugDock.hh>

// Package headers
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/CDRsMinPackMin.hh>
#include <protocols/antibody/RefineOneCDRLoop.hh>

// Project headers
#include <protocols/docking/DockingProtocol.hh>
#include <protocols/docking/DockMCMCycle.hh>
#include <protocols/docking/DockTaskFactory.hh>
#include <protocols/moves/ChangeFoldTreeMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/simple_moves/RepackSidechainsMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/types.hh>

// basic
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/Tracer.hh>


static THREAD_LOCAL basic::Tracer TR( "protocols.antibody.SnugDock" );
using namespace core;

namespace protocols {
namespace antibody {
namespace snugdock {

///////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// BOILER PLATE CODE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

/// @brief default constructor
SnugDock::SnugDock() :
	docking::DockingHighRes(),
	loop_refinement_method_( "refine_kic" ),
	h3_filter_( false ),
	debug_( false ),
	h3_filter_tolerance_( 20 ),
	high_res_kink_constraint_( false ),
	number_of_high_resolution_cycles_( 50 )
{
	init();
}

/// @brief copy constructor
SnugDock::SnugDock( SnugDock const & rhs ) : docking::DockingHighRes(rhs) {
	init_for_equal_operator_and_copy_constructor( *this, rhs );
}

/// @brief assignment operator
SnugDock & SnugDock::operator=( SnugDock const & rhs ) {
	//abort self-assignment
	if ( this == &rhs ) return *this;
	Mover::operator=( rhs );
	init_for_equal_operator_and_copy_constructor( *this, rhs );
	return *this;
}

//destructor
SnugDock::~SnugDock() {}

/// @brief Each derived class must specify its name.
std::string SnugDock::get_name() const {
	return type();
}

//@brief clone operator, calls the copy constructor
protocols::moves::MoverOP
SnugDock::clone() const {
	return protocols::moves::MoverOP( new SnugDock( *this ) );
}

/// @brief fresh_instance returns a default-constructed object for JD2
protocols::moves::MoverOP
SnugDock::fresh_instance() const {
	return protocols::moves::MoverOP( new SnugDock() );
}

/// @brief This mover retains state such that a fresh version is needed if the input Pose is about to change
bool SnugDock::reinitialize_for_new_input() const {
	return true;
}

void SnugDock::debug() {
	debug_ = true;
	number_of_high_resolution_cycles_ = 1;
}

void SnugDock::register_options() {
	moves::RandomMover::register_options();
	moves::SequenceMover::register_options();
	moves::ChangeFoldTreeMover::register_options();
	moves::TrialMover::register_options();
	docking::DockMCMCycle::register_options();
}
///////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////// END OF BOILER PLATE CODE //////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////

void SnugDock::init() {
	Mover::type( "SnugDock" );

	set_default();
	init_from_options();

}

void SnugDock::init_for_equal_operator_and_copy_constructor(SnugDock & lhs, SnugDock const & rhs) {
	// copy all data members from rhs to lhs
	lhs.antibody_info_ = rhs.antibody_info_;
	lhs.mc_ = rhs.mc_;

	// Movers
	lhs.high_resolution_step_ = rhs.high_resolution_step_;
	lhs.loop_refinement_method_ = rhs.loop_refinement_method_;
	lhs.pre_minimization_ = rhs.pre_minimization_;

	// H3 filter options
	lhs.h3_filter_ = rhs.h3_filter_;
	lhs.debug_ = rhs.debug_;
	lhs.h3_filter_tolerance_ = rhs.h3_filter_tolerance_;

	lhs.number_of_high_resolution_cycles_ = rhs.number_of_high_resolution_cycles_;
}
void SnugDock::set_default() {

	loop_refinement_method_ = "refine_kic";
	h3_filter_ = false;
	h3_filter_tolerance_= 20;
	number_of_high_resolution_cycles_ = 50;
	high_res_kink_constraint_ = false;

}

void SnugDock::init_from_options() {
	using basic::options::option;
	using namespace basic::options::OptionKeys;

	if ( option[ basic::options::OptionKeys::antibody::refine ].user() ) {
		loop_refinement_method_  = option[ basic::options::OptionKeys::antibody::refine ]() ;
	}
	/// Allow h3_filter to be turned on at expense of extra loop modeling
	if ( option[ basic::options::OptionKeys::antibody::h3_filter ].user() ) {
		h3_filter_  = option[ basic::options::OptionKeys::antibody::h3_filter ]() ;
	}
	if ( option[ basic::options::OptionKeys::antibody::h3_filter_tolerance ].user() ) {
		h3_filter_tolerance_  = option[ basic::options::OptionKeys::antibody::h3_filter_tolerance ]() ;
	}

	if ( option[ run::test_cycles ].user() ) {
		/// Ideally we would test a larger number of cycles because we are using a random mover in the apply,
		/// but because each submove can be quite long, this would take far too long.
		/// TODO: Create a scientific test for SnugDock that is run regularly.
		number_of_high_resolution_cycles( 5 );
	}
}

void SnugDock::apply( Pose & pose ) {
	TR << "Beginning apply function of " + get_name() + "." << std::endl;
	show( TR );

	if ( ! high_resolution_step_ ) setup_objects( pose );

	/// minimize the CDRs before move to full-atom SnugDock cycle. Remove clashes which may dissociate L-H
	pre_minimization_->apply(pose);

	TR << "Reinitializing the shared MC object before applying the high resolution phase of " + get_name() + "."
		<< std::endl;

	( * scorefxn() )( pose );
	mc_->reset( pose );

	for ( core::Size i = 0; i < number_of_high_resolution_cycles_; ++i ) {
		high_resolution_step_->apply( pose );
		TR << "SnugDock Moves Acceptance Rate" << std::endl;
		mc_->show_counters();
	}

	TR << "Setting the structure to the state with the best score observed during the simulation" << std::endl;
	mc_->recover_low( pose );

	/// Set the pose's foldtree to Ab-Ag docking (LH_A) to ensure the correct interface.
	pose.fold_tree(antibody_info_->get_FoldTree_LH_A(pose));
}

Size SnugDock::number_of_high_resolution_cycles() const {
	return number_of_high_resolution_cycles_;
}

void SnugDock::number_of_high_resolution_cycles( Size const number_of_high_resolution_cycles ) {
	number_of_high_resolution_cycles_ = number_of_high_resolution_cycles;
}

void SnugDock::set_antibody_info( AntibodyInfoOP antibody_info ) {
	antibody_info_ = antibody_info;
}

void SnugDock::setup_objects( Pose const & pose ) {
	using core::scoring::ScoreFunctionOP;
	using docking::DockMCMCycle;
	using docking::DockMCMCycleOP;
	using moves::ChangeFoldTreeMover;
	using moves::ChangeFoldTreeMoverOP;
	using moves::MonteCarloOP;
	using moves::RandomMover;
	using moves::RandomMoverOP;
	using moves::SequenceMover;
	using moves::SequenceMoverOP;
	using moves::TrialMover;
	using moves::TrialMoverOP;

	TR << "Setting up data for " + get_name() + "." << std::endl;

	/// AntibodyInfo is used to store information about the Ab-Ag complex and to generate useful helper objects based on
	/// that information (e.g. the various FoldTrees that are needed for SnugDock).
	if ( ! antibody_info_ ) antibody_info_ = AntibodyInfoOP( new AntibodyInfo( pose ) );


	pre_minimization_ = antibody::CDRsMinPackMinOP( new CDRsMinPackMin( antibody_info_ ) );


	/// A vanilla DockMCMCycle can be used because AntibodyInfo will always make the first jump in the FoldTree dockable.
	DockMCMCycleOP standard_dock_cycle( new DockMCMCycle );

	/// Set a default docking task factory to the DockMCMCycle.
	/// The DockingHighRes base class provides a mechanism to do this.
	tf2()->create_and_attach_task_factory( this, pose );
	standard_dock_cycle->set_task_factory( task_factory() );

	/// This MonteCarlo instance uses 'standard' with the docking patch and a temperature factor of 0.8.
	/// All movers in the high resolution step will share this MonteCarlo instance to provide consistent results.
	mc_ = standard_dock_cycle->get_mc();

	ChangeFoldTreeMoverOP set_foldtree_for_ab_ag_docking( new ChangeFoldTreeMover(
		antibody_info_->get_FoldTree_LH_A( pose )
		) );
	ChangeFoldTreeMoverOP set_foldtree_for_vH_vL_docking( new ChangeFoldTreeMover(
		antibody_info_->get_FoldTree_L_HA( pose )
		) );

	SequenceMoverOP antibody_antigen_dock_cycle( new SequenceMover(
		set_foldtree_for_ab_ag_docking,
		standard_dock_cycle
		) );

	SequenceMoverOP vH_vL_dock_cycle( new SequenceMover(
		set_foldtree_for_vH_vL_docking,
		standard_dock_cycle
		) );

	/// TODO: Does CDRsMinPackMin need a TaskFactory to be set?  Does it get this from AntibodyInfo?
	CDRsMinPackMinOP minimize_all_cdr_loops_base( new CDRsMinPackMin( antibody_info_ ) );
	TrialMoverOP minimize_all_cdr_loops( new TrialMover( minimize_all_cdr_loops_base, mc_ ) );

	/// FIXME: The chain break weight configuration and constraint weight should be handled by RefineOneCDRLoop.
	ScoreFunctionOP high_res_loop_refinement_scorefxn = scorefxn_pack()->clone();
	high_res_loop_refinement_scorefxn->set_weight( scoring::chainbreak, 1.0 );
	high_res_loop_refinement_scorefxn->set_weight( scoring::overlap_chainbreak, 10./3. );
	high_res_loop_refinement_scorefxn->set_weight( scoring::atom_pair_constraint, 100 );

	// update high-res sfxn with kink constraint, passed on from SnugDockProtocol
	if ( high_res_kink_constraint() ) {
		high_res_loop_refinement_scorefxn->set_weight( scoring::dihedral_constraint, 1.0 );
		high_res_loop_refinement_scorefxn->set_weight( scoring::angle_constraint, 1.0 );
	}

	RefineOneCDRLoopOP refine_cdr_h2_base( new RefineOneCDRLoop( antibody_info_, h2, loop_refinement_method_, high_res_loop_refinement_scorefxn ) );
	refine_cdr_h2_base->set_h3_filter( false );
	TrialMoverOP refine_cdr_h2( new TrialMover( refine_cdr_h2_base, mc_ ) );

	RefineOneCDRLoopOP refine_cdr_h3_base( new RefineOneCDRLoop( antibody_info_, h3, loop_refinement_method_, high_res_loop_refinement_scorefxn ) );
	refine_cdr_h3_base->set_h3_filter( h3_filter_ );
	refine_cdr_h3_base->set_num_filter_tries( h3_filter_tolerance_ );
	TrialMoverOP refine_cdr_h3( new TrialMover( refine_cdr_h3_base, mc_ ) );


	/// This is a very succinct description of what this mover does.  For a description in words, see the implementation
	/// of the streaming operator.
	if ( debug_ ) {
		high_resolution_step_ = moves::MoverContainerOP( new SequenceMover );
		high_resolution_step_->add_mover( antibody_antigen_dock_cycle );
		high_resolution_step_->add_mover( vH_vL_dock_cycle);
		high_resolution_step_->add_mover( minimize_all_cdr_loops);
		high_resolution_step_->add_mover( refine_cdr_h2 );
		high_resolution_step_->add_mover( refine_cdr_h3 );
	} else {
		high_resolution_step_ = moves::MoverContainerOP( new RandomMover );
		high_resolution_step_->add_mover( antibody_antigen_dock_cycle, 0.4 );
		high_resolution_step_->add_mover( vH_vL_dock_cycle, 0.4 );
		high_resolution_step_->add_mover( minimize_all_cdr_loops, 0.1 );
		high_resolution_step_->add_mover( refine_cdr_h2, 0.05 );
		high_resolution_step_->add_mover( refine_cdr_h3, 0.05 );
	}
}

void
SnugDock::show( std::ostream & out ) const {
	out << *this;
}

std::ostream & operator<<(std::ostream& out, SnugDock const & ) {
	out << "//////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
	out << "/// The following description borrows heavily from Fig. 1 of:" << std::endl;
	out << "/// Sircar A, Gray JJ (2010) SnugDock: Paratope Structural Optimization during" << std::endl;
	out << "/// Antibody-Antigen Docking Compensates for Errors in Antibody Homology Models." << std::endl;
	out << "/// PLoS Comput Biol 6(1): e1000644. doi:10.1371/journal.pcbi.1000644" << std::endl;
	out << "//////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
	out << "//////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
	out << "///" << std::endl;
	out << "/// Each iteration of this loop will perform one of five different perturbations followed by" << std::endl;
	out << "/// packing and minimization of the relevant regions of the pose and then with a Monte Carlo" << std::endl;
	out << "/// test of the Boltzmann criterion in the indicated frequencies:" << std::endl;
	out << "/// 1. 40% Perturb rigid body position of antibody relative to antigen" << std::endl;
	out << "///        a. Set the pose's FoldTree for Ab-Ag docking (FoldTree provided by AntibodyInfo)" << std::endl;
	out << "///        b. Run one cycle of standard high resolution docking (DockMCMCycle)" << std::endl;
	out << "/// 2. 40% Perturb rigid body position of vL relative to vH" << std::endl;
	out << "///        a. Set the pose's FoldTree for vL-vHAg (FoldTree provided by AntibodyInfo)" << std::endl;
	out << "///        b. Run one cycle of standard high resolution docking (DockMCMCycle)" << std::endl;
	out << "/// 3. 10% Select all CDRs for minimization" << std::endl;
	out << "///        a. Apply CDRsMinPackMin" << std::endl;
	out << "///        b. Monte Carlo accept or reject" << std::endl;
	out << "/// 4.  5% Perturb CDR H2 loop by small, shear and CCD moves followed by minimization" << std::endl;
	out << "///        a. Apply RefineOneCDRLoop (AntibodyInfo will setup the FoldTree for H2 refinement)" << std::endl;
	out << "///        b. Monte Carlo accept or reject" << std::endl;
	out << "/// 5.  5% Perturb CDR H3 loop by small, shear and CCD moves followed by minimization" << std::endl;
	out << "///        a. Apply RefineOneCDRLoop (AntibodyInfo will setup the FoldTree for H3 refinement)" << std::endl;
	out << "///        b. Monte Carlo accept or reject" << std::endl;
	out << "///" << std::endl;
	out << "//////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
	out << "//////////////////////////////////////////////////////////////////////////////////////////////" << std::endl;
	return out;
}

} // namespace snugdock
} // namespace antibody
} // namespace protocols
