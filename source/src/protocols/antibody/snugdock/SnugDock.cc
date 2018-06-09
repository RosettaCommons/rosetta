// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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
#include <protocols/antibody/snugdock/SnugDockCreator.hh>

// Package headers
#include <protocols/loops/Loops.hh>
#include <protocols/loops/loops_main.hh>
#include <protocols/antibody/AntibodyInfo.hh>
#include <protocols/antibody/CDRsMinPackMin.hh>
#include <protocols/antibody/RefineOneCDRLoop.hh>
#include <protocols/antibody/design/util.hh>

// Project headers
#include <protocols/docking/DockingProtocol.hh>
#include <protocols/docking/DockMCMCycle.hh>
#include <protocols/docking/DockTaskFactory.hh>
#include <protocols/moves/ChangeFoldTreeMover.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/moves/TrialMover.hh>
#include <protocols/minimization_packing/RepackSidechainsMover.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>
#include <protocols/simple_moves/DumpPdb.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/types.hh>

// basic
#include <basic/options/option.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/antibody.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/moves/mover_schemas.hh>


static basic::Tracer TR( "protocols.antibody.SnugDock" );
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
SnugDock::~SnugDock() = default;

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

	( * scorefxn() )( pose ); // necessary because setup object uses loops:select_loops_residues()
	if ( ! high_resolution_step_ ) setup_objects( pose );

	/// minimize the CDRs before move to full-atom SnugDock cycle. Remove clashes which may dissociate L-H
	pre_minimization_->apply(pose);

	TR << "Reinitializing the shared MC object before applying the high resolution phase of " + get_name() + "."
		<< std::endl;

	TR << "The Monte Carlo object!!" << std::endl;
	TR << *mc_ << std::endl;
	mc_->reset( pose );

	for ( core::Size i = 0; i < number_of_high_resolution_cycles_; ++i ) {
		high_resolution_step_->apply( pose );
		TR << "SnugDock Moves Acceptance Rate" << std::endl;
		mc_->show_counters();
	}

	TR << "Setting the structure to the state with the best score observed during the simulation" << std::endl;
	mc_->recover_low( pose );

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

void SnugDock::set_CDR_loops( std::map<CDRNameEnum,loops::Loop> cdr_loops ) {
	CDR_loops_ = cdr_loops;
}

void SnugDock::set_ab_has_light_chain( bool ab_has_light_chain ) {
	ab_has_light_chain_ = ab_has_light_chain;
}

void SnugDock::set_vh_vl_jump( bool vh_vl_jump ) {
	vh_vl_jump_ = vh_vl_jump;
}

void SnugDock::set_task_factory( core::pack::task::TaskFactoryCOP tf ) {
	tf_ = tf->clone();
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

	// create LoopsOP of all loops
	loops::LoopsOP CDR_loopsop ( new loops::Loops() );
	// add in the heavy chain CDRs
	CDR_loopsop->add_loop( CDR_loops_[h1] );
	CDR_loopsop->add_loop( CDR_loops_[h2] );
	CDR_loopsop->add_loop( CDR_loops_[h3] );
	// add in the light chain CDRs
	if ( ab_has_light_chain_ ) {
		CDR_loopsop->add_loop( CDR_loops_[l1] );
		CDR_loopsop->add_loop( CDR_loops_[l2] );
		CDR_loopsop->add_loop( CDR_loops_[l3] );
	}

	// set up movemap for minimization from loops
	// when the universal FT is compatible with AntibodyInfo
	// this should be removed and handled that way!
	kinematics::MoveMapOP CDR_movemap ( new kinematics::MoveMap() );
	CDR_movemap->clear();
	CDR_movemap->set_chi( false );
	CDR_movemap->set_bb( false );
	CDR_movemap->set_jump( false );

	utility::vector1< bool> bb_is_flexible( pose.size(), false );
	utility::vector1< bool> sc_is_flexible( pose.size(), false );

	TR << "Here are my loops:\n" << *CDR_loopsop << std::endl;

	loops::select_loop_residues( pose, *CDR_loopsop, false, bb_is_flexible, 10.0);
	CDR_movemap->set_bb( bb_is_flexible );

	loops::select_loop_residues( pose, *CDR_loopsop, true, sc_is_flexible, 10.0);
	CDR_movemap->set_chi( sc_is_flexible );

	// use the non-AntibodyInfo dependent CDRsMinPackMin (set TF later)
	pre_minimization_ = antibody::CDRsMinPackMinOP( new CDRsMinPackMin() );
	pre_minimization_->set_alter_foldtree( false );
	pre_minimization_->set_move_map( CDR_movemap );

	// Alter this because in the universal fold tree jump 1 is Ab-Ag, whereas jump 3 is H-L if it exists
	//DockMCMCycleOP standard_dock_cycle( new DockMCMCycle );
	DockMCMCycleOP ab_ag_dock_cycle( new DockMCMCycle );
	DockMCMCycleOP vh_vl_dock_cycle;

	if ( ab_has_light_chain_ ) {
		// when specifying a jump for docking, the sfxn has to be specified
		// for both docking packing -- no idea why, so I'm using the defaults
		vh_vl_dock_cycle = DockMCMCycleOP( new DockMCMCycle(vh_vl_jump_, core::scoring::ScoreFunctionFactory::create_score_function("docking", "docking_min"), core::scoring::get_score_function_legacy(core::scoring::PRE_TALARIS_2013_STANDARD_WTS)) );
		vh_vl_dock_cycle->set_task_factory( tf_ );
	}

	/// Set a default docking task factory to the DockMCMCycle.
	/// The DockingHighRes base class provides a mechanism to do this.
	ab_ag_dock_cycle->set_task_factory( tf_ );

	// Set the TF for CDRsMinPackMin
	pre_minimization_->set_task_factory( tf_ );

	/// This MonteCarlo instance uses 'standard' with the docking patch and a temperature factor of 0.8.
	/// All movers in the high resolution step will share this MonteCarlo instance to provide consistent results.
	mc_ = ab_ag_dock_cycle->get_mc();

	// wrap MinPackMin in TrialMover
	TrialMoverOP minimize_all_cdr_loops( new TrialMover( pre_minimization_, mc_ ) );

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

	RefineOneCDRLoopOP refine_cdr_h2_base( new RefineOneCDRLoop( CDR_loops_[h2], loop_refinement_method_, high_res_loop_refinement_scorefxn ) );
	TrialMoverOP refine_cdr_h2( new TrialMover( refine_cdr_h2_base, mc_ ) );

	RefineOneCDRLoopOP refine_cdr_h3_base( new RefineOneCDRLoop( CDR_loops_[h3], loop_refinement_method_, high_res_loop_refinement_scorefxn ) );
	TrialMoverOP refine_cdr_h3( new TrialMover( refine_cdr_h3_base, mc_ ) );


	/// This is a very succinct description of what this mover does.  For a description in words, see the implementation
	/// of the streaming operator.
	if ( debug_ ) {
		high_resolution_step_ = moves::MoverContainerOP( new SequenceMover );
		//high_resolution_step_->add_mover( moves::MoverOP(new simple_moves::DumpPdb("pre_highres.pdb")) );
		high_resolution_step_->add_mover( ab_ag_dock_cycle );
		//high_resolution_step_->add_mover( moves::MoverOP(new simple_moves::DumpPdb("ab_ag_dock.pdb")) );
		if ( ab_has_light_chain_ ) {
			high_resolution_step_->add_mover( vh_vl_dock_cycle );
			//high_resolution_step_->add_mover( moves::MoverOP(new simple_moves::DumpPdb("vh_vl_dock.pdb")) );
		}
		high_resolution_step_->add_mover( minimize_all_cdr_loops );
		//high_resolution_step_->add_mover( moves::MoverOP(new simple_moves::DumpPdb("min_all_cdrs.pdb")) );
		high_resolution_step_->add_mover( refine_cdr_h2 );
		//high_resolution_step_->add_mover( moves::MoverOP(new simple_moves::DumpPdb("refine_h2.pdb")) );
		high_resolution_step_->add_mover( refine_cdr_h3 );
		//high_resolution_step_->add_mover( moves::MoverOP(new simple_moves::DumpPdb("refine_h3.pdb")) );
	} else {
		high_resolution_step_ = moves::MoverContainerOP( new RandomMover );
		high_resolution_step_->add_mover( ab_ag_dock_cycle, 0.4 );
		if ( ab_has_light_chain_ ) { high_resolution_step_->add_mover( vh_vl_dock_cycle, 0.4 ); }
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

std::string SnugDock::get_name() const {
	return mover_name();
}

std::string SnugDock::mover_name() {
	return "SnugDock";
}

void SnugDock::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	//This is a DockingHighRes mover and doesn't actually read anything from the tag
	protocols::moves::xsd_type_definition_w_attributes(
		xsd, mover_name(),
		"Dock and antigen to an antibody while optimizing the rigid body orientation of the VH and VL chains and performing CDR loop minimization.",
		attlist );
}

std::string SnugDockCreator::keyname() const {
	return SnugDock::mover_name();
}

protocols::moves::MoverOP
SnugDockCreator::create_mover() const {
	return protocols::moves::MoverOP( new SnugDock );
}

void SnugDockCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	SnugDock::provide_xml_schema( xsd );
}


} // namespace snugdock
} // namespace antibody
} // namespace protocols
