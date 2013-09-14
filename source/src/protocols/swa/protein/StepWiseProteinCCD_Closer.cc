// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseProteinCCD_Closer
/// @brief Makes a list of (phi, psi, omega) at moving_residues that
///              could be useful for full-atom packing
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/swa/protein/StepWiseProteinCCD_Closer.hh>
// AUTO-REMOVED #include <protocols/swa/protein/StepWiseProteinUtil.hh>
#include <protocols/swa/StepWiseJobParameters.hh>
// AUTO-REMOVED #include <protocols/swa/InputStreamWithResidueInfo.hh>

//////////////////////////////////
#include <core/types.hh>
// AUTO-REMOVED #include <core/chemical/util.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
// AUTO-REMOVED #include <core/io/silent/ProteinSilentStruct.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.fwd.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/scoring/Energies.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.tmpl.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/scoring/Ramachandran.hh>
#include <basic/Tracer.hh>
#include <core/kinematics/Jump.hh>

#include <protocols/loops/Loop.hh>
#include <protocols/loops/loop_closure/ccd/ccd_closure.hh>

// AUTO-REMOVED #include <ObjexxFCL/format.hh>
// AUTO-REMOVED #include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/FArray1D.hh>

#include <utility/exit.hh>

#ifdef WIN32
#include <time.h>
#endif

#include <string>

//Auto Headers
#include <core/chemical/ResidueType.hh>
#include <protocols/swa/sample_generators/StepWisePoseSampleGenerator.hh>
#include <utility/vector1.hh>


using namespace core;
using core::Real;

//////////////////////////////////////////////////////////////////////////
// CCD loop closure. Applied to loops whose N-terminal
// and C-terminal segment have been modeled in, e.g., different silent
// files.
//
// Later, should sample omega recursively for pre-prolines.
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.swa.protein.stepwise_protein_ccd_closer" );

namespace protocols {
namespace swa {
namespace protein {

  //////////////////////////////////////////////////////////////////////////
  //constructor!
  StepWiseProteinCCD_Closer::StepWiseProteinCCD_Closer( StepWisePoseSampleGeneratorOP sample_generator,
																												  StepWiseJobParametersOP job_parameters ):
		Mover(),
		sample_generator_( sample_generator ),
		working_bridge_res_( job_parameters->working_bridge_res() ),
		moving_residues_( job_parameters->working_moving_res_list() ),
		is_pre_proline_( job_parameters->is_pre_proline() ),
		ccd_close_res_( 0 ),
		verbose_( true )
  {
		Mover::type("StepWiseProteinCCD_Closer");
		if ( working_bridge_res_.size() >= 3 ) utility_exit_with_message( "Cannot supply more than 2 bridge residues");
  }

  //////////////////////////////////////////////////////////////////////////
  //destructor
  StepWiseProteinCCD_Closer::~StepWiseProteinCCD_Closer()
  {}

	/////////////////////
	std::string
	StepWiseProteinCCD_Closer::get_name() const {
		return "StepWiseProteinCCD_Closer";
	}

  //////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////
  void
  StepWiseProteinCCD_Closer::apply( core::pose::Pose & pose )
	{

		clock_t const time_start( clock() );

		Pose pose_save = pose;

		setup_torsions( pose );

		figure_out_loop( pose );

		sample_generator_->reset();

		pose_count_ = 0;
		while ( sample_generator_->has_another_sample() ){

			sample_generator_->get_next_sample( pose );

			//			if ( pose_count_ == 0 ) pose.dump_pdb( "first_sample_ccd.pdb" );

			save_phi_psi_omega_over_loop_residues( pose ); // trying to avoid memory effects, where solutions will depend on order of when processed.
			CCD_loop_close_sample_omega_recursively( pose, 0 /*offset*/ );
			restore_phi_psi_omega_over_loop_residues( pose );

		}

		Size const num_successes = main_chain_torsion_sets_for_moving_residues_.size();
		std::cout << "CCD closed: " << num_successes << " out of " << pose_count_ << " tries." << std::endl;

		// Kind of silly -- should have at least one output.
		if ( num_successes == 0  ) grab_main_chain_torsion_set_list( pose );

		std::cout << "Total time in StepWiseProteinCCD_Closer: " <<
			static_cast<Real>(clock() - time_start) / CLOCKS_PER_SEC << std::endl;

		pose = pose_save;

	}



	///////////////////////////////////////////////////////////////////////////
	bool
  StepWiseProteinCCD_Closer::CCD_loop_close( core::pose::Pose & pose ) {

		// param for ccd_closure
		int  const ccd_cycles = { 1000 }; // num of cycles of ccd_moves
		Real const ccd_tol = { 0.001 }; // criterion for a closed loop
		bool const rama_check = { true };
		Real const max_rama_score_increase = { 100.0 }; // dummy number when rama_check is false
		Real const max_total_delta_helix   = { 100.0 }; // max overall angle changes for a helical residue
		Real const max_total_delta_strand  = { 500.0 }; // ... for a residue in strand
		Real const max_total_delta_loop    = { 3600.0 }; // ... for a residue in loop
		// cutoff for acceptance
		Real const rmsd_acceptance_cutoff =  1.5; // can be a little relaxed, since there will be a minimize afterwards.
		// output for ccd_closure
		Real forward_deviation, backward_deviation; // actually loop closure msd, both dirs
		Real torsion_delta, rama_delta; // actually torsion and rama score changes, averaged by loop_size

		Size const n_cycles = protocols::loops::loop_closure::ccd::fast_ccd_loop_closure(
																						pose, mm_, loop_.start() , loop_.stop(), loop_.cut(), ccd_cycles,
																						ccd_tol, rama_check, max_rama_score_increase, max_total_delta_helix,
																						max_total_delta_strand, max_total_delta_loop, forward_deviation,
																						backward_deviation, torsion_delta, rama_delta
																						);

		pose_count_++;


		TR << "CCD forward_dev: " << forward_deviation << "  backward_dev: " << backward_deviation << "  torsion_delta: " << torsion_delta << " rama_delta: " << rama_delta << "  number of cycles: " << n_cycles << std::endl;
		if ( forward_deviation > rmsd_acceptance_cutoff || backward_deviation > rmsd_acceptance_cutoff ) return false;

		grab_main_chain_torsion_set_list( pose );

		return true;

	}


	///////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinCCD_Closer::CCD_loop_close_sample_omega_recursively( core::pose::Pose & pose, int const offset ) {

	 	restore_phi_psi_over_loop_residues( pose ); // don't restore omega, since we'll be switching it around.

		//loop starts at position before cutpoint, ends after.
		Size const omega_pos = loop_.start() + offset;

	 	if ( omega_pos == loop_.stop() ){

			CCD_loop_close( pose );

		} else {

			pose.set_omega( omega_pos, 180.0 );
			CCD_loop_close_sample_omega_recursively( pose, offset + 1 );

			if ( is_pre_proline_[ omega_pos ] ){
				pose.set_omega( omega_pos, 0.0 );
				TR << "trying CIS omega! " << omega_pos << std::endl;
				CCD_loop_close_sample_omega_recursively( pose, offset + 1 );
			}

		}
	}



	//////////////////////////////////////////////////////////////////////////////////////////////
	// note 'loop' includes takeoff and end points
	// movemap will ensure that we only sample psi of N-terminal takeoff, and phi of C-terminal landing.
	void
	StepWiseProteinCCD_Closer::figure_out_loop( core::pose::Pose const & pose ){

		using namespace protocols::loops;
		using namespace core::chemical;
		using namespace core::id;

		Size loop_start( 0 ), loop_end( 0 ), cutpoint( 0 );

		Size const num_bridge_res = working_bridge_res_.size();

		if ( ccd_close_res_ > 0 ) {

			cutpoint = ccd_close_res_;
			if ( !pose.residue_type( cutpoint ).has_variant_type( CUTPOINT_LOWER ) ||
					 !pose.residue_type( cutpoint+1 ).has_variant_type( CUTPOINT_UPPER ) ) utility_exit_with_message( "ccd_close_res needs to be specified as a cutpoint_closed!" );

		} else {

			// find the cutpoint...
			for ( Size n = 1; n <= pose.total_residue() ; n++ ){

				if ( num_bridge_res > 0  &&  ( (n < working_bridge_res_[1]-1) || (n > working_bridge_res_[num_bridge_res] ) ) ) continue;

				if ( pose.residue_type( n ).has_variant_type( CUTPOINT_LOWER ) &&
						 pose.residue_type( n+1 ).has_variant_type( CUTPOINT_UPPER ) ){
					if ( cutpoint > 0 ) utility_exit_with_message( "Cannot specify more than one closed cutpoint inside the bridge_residues in StepWiseProteinCCD_Closer! Consider using -ccd_close_res to specify, or else stick to one cutpoint_closed." );
					cutpoint = n;
				}
			}

		}

		if ( cutpoint == 0 ) utility_exit_with_message( "Must specify one closed cutpoint within bridge_res in StepWiseProteinCCD_Closer!" );

		if ( num_bridge_res  == 0 ) {
			// could be just twiddling torsions at cutpoint.
			loop_start = cutpoint;
			loop_end   = cutpoint+1;
		} else {
			loop_start = working_bridge_res_[1]-1;
			loop_end   = working_bridge_res_[num_bridge_res]+1;
		}

		if ( cutpoint < loop_start || cutpoint >= loop_end ) utility_exit_with_message( "cutpoint must lie within loop" );

		loop_ = Loop( loop_start, loop_end, cutpoint );

		// set up movemap.
		mm_.clear();
		mm_.set( TorsionID( loop_start, id::BB, 2 ),  true );
		for ( Size n = loop_start+1; n <= loop_end-1; n++ ){
			mm_.set( TorsionID( n, id::BB, 1 ),  true );
			mm_.set( TorsionID( n, id::BB, 2 ),  true );
		}
		mm_.set( TorsionID( loop_end, id::BB, 1 ),  true );

	}



	//////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinCCD_Closer::setup_torsions( pose::Pose const & pose ){

		// Need to fill torsions for moving and bridge residues.

		using namespace core::id;

		utility::vector1< bool > is_moving_or_bridge_res( pose.total_residue(), false );

		for ( Size n = 1; n <= moving_residues_.size(); n++ ) 	{
			is_moving_or_bridge_res[ moving_residues_[n] ] = true;
		}

		for ( Size n = 1; n <= working_bridge_res_.size(); n++ ) 	{
			is_moving_or_bridge_res[ working_bridge_res_[n] ] = true;
		}

		which_torsions_.clear();
		for ( Size n = 1; n <= pose.total_residue(); n++ ){

			if ( is_moving_or_bridge_res[ n ]  ) {

				//save phi,psi,omega inside loop
				for ( Size k = 1; k <= 3; k++ )	which_torsions_.push_back( TorsionID( n, BB, k ) );

			} else if ( n < pose.total_residue() && is_moving_or_bridge_res[ n+1 ] ) {

				//save psi,omega at 'takeoff' residue.
				for ( Size k = 2; k <= 3; k++ )	which_torsions_.push_back( TorsionID( n, BB, k ) );

			} else if ( n > 1 && is_moving_or_bridge_res[ n-1 ] ) {

				//save phi at 'landing' residue.
				for ( Size k = 1; k <= 1; k++ )	which_torsions_.push_back( TorsionID( n, BB, k ) );

			}

		}

	}


  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinCCD_Closer::grab_main_chain_torsion_set_list( pose::Pose const & pose ){

		utility::vector1< Real > main_chain_torsion_set_for_moving_residues;
		for ( Size n = 1; n <= which_torsions_.size(); n++ ) 	{
			main_chain_torsion_set_for_moving_residues.push_back(  pose.torsion( which_torsions_[ n ] ) );
		}
		main_chain_torsion_sets_for_moving_residues_.push_back( main_chain_torsion_set_for_moving_residues );

	}

  //////////////////////////////////////////////////////////////////////////
	utility::vector1< utility::vector1< core::Real > > const &
	StepWiseProteinCCD_Closer::main_chain_torsion_set_lists() const
	{
		return main_chain_torsion_sets_for_moving_residues_;
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
		for ( Size n = 1; n <= which_torsions_.size(); n++ ) 	{
			pose.set_torsion( which_torsions_[ n ], main_chain_torsion_set_for_moving_residues_save_[n] );
		}
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinCCD_Closer::restore_phi_psi_over_loop_residues( pose::Pose & pose )
	{
		for ( Size n = 1; n <= which_torsions_.size(); n++ ) 	{
			if ( which_torsions_[n].torsion() > 2 ) continue; // phi, psi, but not omega
			pose.set_torsion( which_torsions_[ n ], main_chain_torsion_set_for_moving_residues_save_[n] );
		}
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseProteinCCD_Closer::save_phi_psi_omega_over_loop_residues( pose::Pose const & pose ){
		main_chain_torsion_set_for_moving_residues_save_.clear();
		for ( Size n = 1; n <= which_torsions_.size(); n++ ) 	{
			main_chain_torsion_set_for_moving_residues_save_.push_back(  pose.torsion( which_torsions_[ n ] ) );
		}
	}



}
}
}
