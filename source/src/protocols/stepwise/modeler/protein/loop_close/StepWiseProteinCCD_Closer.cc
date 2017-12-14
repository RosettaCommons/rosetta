// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file StepWiseProteinCCD_Closer
/// @brief Makes a list of (phi, psi, omega) at moving_residues that could be useful for full-atom packing
/// @details
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/modeler/protein/loop_close/StepWiseProteinCCD_Closer.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/util.hh>
#include <protocols/simple_moves/TorsionSetMover.hh>
#include <protocols/loops/loop_closure/ccd/CCDLoopClosureMover.hh>
#include <protocols/loops/loop_closure/ccd/RamaCheck.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerSized.hh>
#include <protocols/loops/Loop.hh>
#include <core/types.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/kinematics/Jump.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/string.functions.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>

#ifdef WIN32
#include <time.h>
#endif

#include <algorithm>
#include <string>


using namespace core;
using core::Real;
using namespace protocols::stepwise::modeler;

//////////////////////////////////////////////////////////////////////////
// CCD loop closure. Applied to loops whose N-terminal
// and C-terminal segment have been modeled in, e.g., different silent
// files.
//
// Samples omega recursively for pre-prolines.
//
// Assumes that 'bridge_res' are contiguous and encompass cutpoint_closed.
//
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.stepwise.modeler.protein.loop_close.StepWiseProteinCCD_Closer" );

namespace protocols {
namespace stepwise {
namespace modeler {
namespace protein {
namespace loop_close {

//////////////////////////////////////////////////////////////////////////
//constructor!
StepWiseProteinCCD_Closer::StepWiseProteinCCD_Closer( working_parameters::StepWiseWorkingParametersCOP working_parameters ):
	ccd_loop_closure_mover_( /* NULL */ ),
	working_bridge_res_( working_parameters->working_bridge_res() ),
	moving_residues_( working_parameters->working_moving_res_list() ),
	is_pre_proline_( working_parameters->is_pre_proline() ),
	ccd_close_res_( 0 ),
	mm_( core::kinematics::MoveMapOP( new kinematics::MoveMap ) ),
	closed_loop_( false ),
	ntries_( 0 )
{
	Mover::type("StepWiseProteinCCD_Closer");
}

//////////////////////////////////////////////////////////////////////////
//destructor
StepWiseProteinCCD_Closer::~StepWiseProteinCCD_Closer() = default;

/////////////////////
std::string
StepWiseProteinCCD_Closer::get_name() const {
	return "StepWiseProteinCCD_Closer";
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinCCD_Closer::get_closure_solution( core::pose::Pose & pose )
{
	closed_loop_ = false;
	save_phi_psi_omega_over_loop_residues( pose ); // trying to avoid memory effects, where solutions will depend on order of when processed.
	CCD_loop_close_sample_omega_recursively( pose, 0 /*offset*/ );
	restore_phi_psi_omega_over_loop_residues( pose );
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinCCD_Closer::apply( core::pose::Pose & pose )
{
	get_closure_solution( pose );

	if ( main_chain_torsion_set().size() == 0 ) return;
	simple_moves::TorsionSetMover torsion_set_mover( which_torsions(), main_chain_torsion_set() );
	torsion_set_mover.apply( pose );
}

///////////////////////////////////////////////////////////////////////////
void
StepWiseProteinCCD_Closer::init( core::pose::Pose & pose ) {
	figure_out_loop( pose );
	figure_out_movemap( pose );
	setup_torsions( pose );
	fix_jump_atoms_at_loop_boundaries( pose );
	ntries_ = 0;
}

///////////////////////////////////////////////////////////////////////////
bool
StepWiseProteinCCD_Closer::CCD_loop_close( core::pose::Pose & pose )
{
	// Lazy initialization of ccd_loop_closure_mover_
	if ( ccd_loop_closure_mover_.get() == nullptr ) {
		ccd_loop_closure_mover_ = loops::loop_closure::ccd::CCDLoopClosureMoverOP( new loops::loop_closure::ccd::CCDLoopClosureMover() );
		ccd_loop_closure_mover_->max_cycles( 1000 );
		ccd_loop_closure_mover_->tolerance( 0.001 );
		ccd_loop_closure_mover_->check_rama_scores( false );
		ccd_loop_closure_mover_->set_verbose( false );
		ccd_loop_closure_mover_->rama()->max_rama_score_increase( 100.0 );

		// FYI: I think anything over 180 is meaningless below. ~Labonte
		ccd_loop_closure_mover_->max_total_torsion_delta_per_residue( 100.0, 500.0, 3600.0 );
	}

	// Does a typical use of this have changing Loop and/or MoveMap definitions?
	ccd_loop_closure_mover_->loop(  loop_ );
	ccd_loop_closure_mover_->movemap( mm_ );

	// cutoff for acceptance
	Real const rmsd_acceptance_cutoff =  1.5; // can be a little relaxed, since there will be a minimize afterwards.

	ccd_loop_closure_mover_->apply( pose );

	ntries_++;

	// FYI: CCDLoopClosureMover has a success() method. ~Labonte
	bool const loop_closed = ( ccd_loop_closure_mover_->deviation() <= rmsd_acceptance_cutoff );

	if ( TR.Debug.visible() ) {
		TR.Debug << "CCD dev: " << ccd_loop_closure_mover_->deviation() <<
			"  torsion_delta: " << ccd_loop_closure_mover_->torsion_delta() <<
			" rama_delta: " << ccd_loop_closure_mover_->rama_delta() <<
			"  number of cycles: " << ccd_loop_closure_mover_->actual_cycles() <<
			"  loop_closed:  " << loop_closed << std::endl;
	}

	if ( !loop_closed ) return false;
	grab_main_chain_torsion_set_list( pose );

	return true;
}


///////////////////////////////////////////////////////////////////////////
void
StepWiseProteinCCD_Closer::CCD_loop_close_sample_omega_recursively( core::pose::Pose & pose, int const offset ) {

	restore_phi_psi_over_loop_residues( pose ); // don't restore omega, since we'll be switching it around.

	//loop starts at position before cutpoint, ends after.
	Size const omega_pos = loop_.start() + offset;

	if ( omega_pos == loop_.stop() ) {
		closed_loop_ = CCD_loop_close( pose );
	} else {
		pose.set_omega( omega_pos, 180.0 );
		CCD_loop_close_sample_omega_recursively( pose, offset + 1 );

		if ( is_pre_proline_[ omega_pos ] ) {
			pose.set_omega( omega_pos, 0.0 );
			TR.Debug << "trying CIS omega! " << omega_pos << std::endl;
			CCD_loop_close_sample_omega_recursively( pose, offset + 1 );
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////
// note 'loop' includes takeoff and end points
// movemap will ensure that we only sample psi of N-terminal takeoff, and phi of C-terminal landing.
void
StepWiseProteinCCD_Closer::figure_out_loop( core::pose::Pose const & pose ) {

	using namespace protocols::loops;

	Size cutpoint = ccd_close_res_; // could be user input
	if ( cutpoint == 0 ) cutpoint = check_for_unique_cutpoint_flanked_by_bridge_res( pose );
	if ( cutpoint == 0 ) cutpoint = check_for_unique_cutpoint( pose );
	runtime_assert( cutpoint != 0 );
	runtime_assert( is_cutpoint_closed( pose, cutpoint ) );

	Size loop_start( cutpoint );
	while ( loop_start > 1 && working_bridge_res_.has_value( loop_start ) ) loop_start--;

	Size loop_stop( cutpoint + 1 );
	while ( loop_stop < pose.size() && working_bridge_res_.has_value( loop_stop ) ) loop_stop++;

	runtime_assert( ( int(loop_stop) - int(loop_start) - 1) < 4 ); // sanity check: can't have more than 3 bridge residues.
	loop_ = Loop( loop_start, loop_stop, cutpoint );
}

//////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseProteinCCD_Closer::figure_out_movemap( core::pose::Pose const & pose ) {
	using namespace core::id;
	mm_->clear();
	for ( Size k = 2; k < pose.residue( loop_.start() ).mainchain_torsions().size(); k++ ) {
		mm_->set( TorsionID( loop_.start(), id::BB, k ),  true );
	}
	for ( Size n = loop_.start()+1; n <= loop_.stop()-1; n++ ) {
		// Set non-omega true
		for ( Size k = 1; k < pose.residue( n ).mainchain_torsions().size(); k++ ) {
			mm_->set( TorsionID( n, id::BB, k ),  true );
		}
	}
	mm_->set( TorsionID( loop_.stop(), id::BB, 1 ),  true );
}

//////////////////////////////////////////////////////////////////////////////////////////////
void
StepWiseProteinCCD_Closer::setup_torsions( core::pose::Pose const & pose ) {

	using namespace core::id;
	which_torsions_.clear();

	//save psi,omega at 'takeoff' residue.
	for ( Size k = 2; k <= pose.residue( loop_.start() ).mainchain_torsions().size(); k++ ) which_torsions_.push_back( TorsionID( loop_.start(), BB, k ) );

	//save all mc torsions (for alphas, phi,psi,omega) inside loop
	for ( Size n = loop_.start() + 1; n < loop_.stop(); n++ ) {
		for ( Size k = 1; k <= pose.residue( loop_.start() ).mainchain_torsions().size(); k++ ) which_torsions_.push_back( TorsionID( n, BB, k ) );
	}

	//save phi at 'landing' residue.
	which_torsions_.push_back( TorsionID( loop_.stop(), BB, 1 ) );
}


//////////////////////////////////////////////////////////////////////////
utility::vector1< core::Real >
StepWiseProteinCCD_Closer::grab_main_chain_torsion_set_list( pose::Pose const & pose ){
	main_chain_torsion_set_.resize( which_torsions_.size() );
	std::transform(
		which_torsions_.begin(), which_torsions_.end(), main_chain_torsion_set_.begin(),
		[&]( core::id::TorsionID const & torsion ) {
			return pose.torsion( torsion );
		} );
	return main_chain_torsion_set_;
}

//////////////////////////////////////////////////////////////////////////
utility::vector1< core::Real > const &
StepWiseProteinCCD_Closer::main_chain_torsion_set() const
{
	return main_chain_torsion_set_;
}

//////////////////////////////////////////////////////////////////////////
utility::vector1< core::Real > const &
StepWiseProteinCCD_Closer::main_chain_torsion_set_save() const
{
	return main_chain_torsion_set_save_;
}

//////////////////////////////////////////////////////////////////////////
utility::vector1< core::id::TorsionID > const &
StepWiseProteinCCD_Closer::which_torsions() const
{
	return which_torsions_;
}


//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinCCD_Closer::restore_phi_psi_omega_over_loop_residues( pose::Pose & pose )
{
	for ( Size n = 1; n <= which_torsions_.size(); n++ )  {
		pose.set_torsion( which_torsions_[ n ], main_chain_torsion_set_save_[n] );
	}
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinCCD_Closer::restore_phi_psi_over_loop_residues( pose::Pose & pose )
{
	for ( Size n = 1; n <= which_torsions_.size(); n++ )  {
		if ( which_torsions_[n].torsion() > 2 ) continue; // phi, psi, but not omega
		pose.set_torsion( which_torsions_[ n ], main_chain_torsion_set_save_[n] );
	}
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinCCD_Closer::save_phi_psi_omega_over_loop_residues( pose::Pose const & pose ){
	main_chain_torsion_set_save_.clear();
	//  TR << ccd_close_res_ << " ";
	for ( auto const torsion : which_torsions_ )  {
		main_chain_torsion_set_save_.push_back(  pose.torsion( torsion ) );
		//   TR << "   " << which_torsions_[n] << " " << pose.torsion( torsion );
	}
	//  TR << std::endl;
}

//////////////////////////////////////////////////////////////////////////
void
StepWiseProteinCCD_Closer::fix_jump_atoms_at_loop_boundaries( pose::Pose & pose ){
	fix_protein_jump_atom( pose, loop_.start(), " N  " );
	fix_protein_jump_atom( pose, loop_.stop(),  " C  " );
}

//////////////////////////////////////////////////////////////////////////
Size
StepWiseProteinCCD_Closer::check_for_unique_cutpoint_flanked_by_bridge_res( pose::Pose const & pose ){
	Size cutpoint( 0 );
	Size nres( pose.size() );
	for ( Size n = 1; n <= nres; n++ ) {
		if ( ! is_cutpoint_closed( pose, n ) ) continue;

		if ( ( n > 1    && working_bridge_res_.has_value( n-1 ) ) ||
				( n < nres && working_bridge_res_.has_value( n+1 ) ) ) {
			runtime_assert( cutpoint == 0 ); // uniqueness check.
			cutpoint = n;
		}
	}
	return cutpoint;
}

//////////////////////////////////////////////////////////////////////////
Size
StepWiseProteinCCD_Closer::check_for_unique_cutpoint( pose::Pose const & pose ){
	Size cutpoint( 0 );
	Size nres( pose.size() );
	for ( Size n = 1; n <= nres; n++ ) {
		if  ( !is_cutpoint_closed( pose, n ) ) continue;

		runtime_assert( cutpoint == 0 ); // uniqueness check.
		cutpoint = n;
	}
	return cutpoint;
}

} //loop_close
} //protein
} //modeler
} //stepwise
} //protocols
