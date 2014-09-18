// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/string.functions.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.modeler.working_parameters.StepWiseWorkingParameters" );

namespace protocols {
namespace stepwise {
namespace modeler {
namespace working_parameters {

	//Constructor
	StepWiseWorkingParameters::StepWiseWorkingParameters():
		// RNA
		output_extra_RMSDs_( false ),
		is_simple_full_length_job_params_( false ),
		five_prime_chain_break_res_( 0 ),
		add_virt_res_as_root_( false ),
		floating_base_( false ),
		floating_base_anchor_res_( 0 ),
		rebuild_bulge_mode_( false ),
		sample_both_sugar_base_rotamer_( false )
		// protein -- no variables.
	{}

	//Destructor
	StepWiseWorkingParameters::~StepWiseWorkingParameters()
	{}

	//////////////////////////////////////////////////////////////////////////////////////////
	// RNA stuff
	//////////////////////////////////////////////////////////////////////////////////////////
	Size
	StepWiseWorkingParameters::working_reference_res() const{ //the last static_residue that this attach to the moving residues

		if ( floating_base_anchor_res_ ) return working_floating_base_anchor_res();

		Size const num_nucleotides = working_moving_res_list_.size(); // this is number of working nucleotides
		//check that moving_res_ and working_moving_res_list list are intialized (note this is not foolproof)
		runtime_assert( num_nucleotides > 0 );
		runtime_assert( working_moving_res_ > 0 );

		Size const working_reference_res_ = ( is_prepend_ ) ? working_moving_res_ + num_nucleotides : working_moving_res_ - num_nucleotides;

		return working_reference_res_;
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	Size
	StepWiseWorkingParameters::gap_size_to_anchor() const { //the last static residue that this attach to the moving residues -- total sequence distance.
		Size const working_reference_res_ = working_reference_res();
		runtime_assert( sub_to_full_.find( working_reference_res_ ) != sub_to_full_.end() );
		Size const & reference_res = sub_to_full_.find( working_reference_res_ )->second;

		// check if these are really separated by a user-defined chain break ('cutpoint open');
		Size const check_cut_start = std::min( moving_res_, reference_res );
		Size const check_cut_stop = std::max( moving_res_, reference_res )-1;
		for ( Size i = check_cut_start; i <= check_cut_stop; i++ ) if ( cutpoint_open_list_.has_value( i ) ) return GAP_SIZE_DUMMY;

		int separation = std::abs( int( reference_res ) - int( moving_res_ ) );
		int gap_size_to_anchor = separation - 1;

		runtime_assert( gap_size_to_anchor >= 0 );
		return static_cast<Size>( gap_size_to_anchor );
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	core::kinematics::FoldTree const & StepWiseWorkingParameters::fold_tree() const{
		//		if ( fold_tree_.size() == 0 ) utility_exit_with_message( "fold_tree_.size() == 0" ); //Thie is the number of edge. simple_tree have 1 edge!
		return fold_tree_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void StepWiseWorkingParameters::set_fold_tree( core::kinematics::FoldTree const & setting ){
		fold_tree_ = setting;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	Size
	StepWiseWorkingParameters::actually_moving_res() const{
		return working_moving_res_;
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseWorkingParameters::set_global_sample_res_list( utility::vector1< core::Size > const & setting ){
		global_sample_res_list_ = setting;
		working_global_sample_res_list_ = apply_full_to_sub_mapping( global_sample_res_list_ );
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseWorkingParameters::set_force_syn_chi_res_list( utility::vector1< core::Size > const & setting ){
		force_syn_chi_res_list_ = setting;
		working_force_syn_chi_res_list_ = apply_full_to_sub_mapping( force_syn_chi_res_list_ );
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseWorkingParameters::set_force_north_sugar_list( utility::vector1< core::Size > const & setting ){
		force_north_sugar_list_ = setting;
		working_force_north_sugar_list_ = apply_full_to_sub_mapping( force_north_sugar_list_ );

		for ( Size n = 1; n <= force_north_sugar_list_.size(); n++ ){
			if ( force_south_sugar_list_.has_value( force_north_sugar_list_[n] ) ){
				utility_exit_with_message( "seq_num = " + ObjexxFCL::string_of( force_north_sugar_list_[n] ) + " is in both force_north_sugar_list_ and force_south_sugar_list_! " );
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseWorkingParameters::set_force_south_sugar_list( utility::vector1< core::Size > const & setting ){
		force_south_sugar_list_ = setting;
		working_force_south_sugar_list_ = apply_full_to_sub_mapping( force_south_sugar_list_ );

		for ( Size n = 1; n <= force_north_sugar_list_.size(); n++ ){
			if ( force_south_sugar_list_.has_value( force_north_sugar_list_[n] ) ){
				utility_exit_with_message( "seq_num = " + ObjexxFCL::string_of( force_north_sugar_list_[n] ) + " is in both force_north_sugar_list_ and force_south_sugar_list_! " );
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseWorkingParameters::set_protonated_H1_adenosine_list( utility::vector1< core::Size > const & setting ){
		protonated_H1_adenosine_list_ = setting;
		working_protonated_H1_adenosine_list_ = apply_full_to_sub_mapping( protonated_H1_adenosine_list_ );
	}
	//////////////////////////////////////////////////////////////////////////////////////////
	Size StepWiseWorkingParameters::working_floating_base_anchor_res() const{
		if ( floating_base_anchor_res_ == 0 ) return 0;
		if ( !is_working_res_[ floating_base_anchor_res_ ] ) return 0;
		return full_to_sub_.find( floating_base_anchor_res_ )->second;
	}


	//////////////////////////////////////////////////////////////////////////////////////////
	// Protein stuff
	//////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< bool > const
	StepWiseWorkingParameters::is_pre_proline() const {
		utility::vector1< bool > is_pre_proline;
		utility::vector1< Size > const working_res( working_res_list() );
		std::string const sequence( full_sequence() );
		for ( Size i = 1; i <= working_res.size(); i++ ){
			Size const & full_seq_pos = working_res[ i ];
			if ( full_seq_pos == sequence.size() ) {
				is_pre_proline.push_back( false );
			} else {
				is_pre_proline.push_back( sequence[ full_seq_pos ] == 'P' ); /*note offset by one -- this is the *next* sequence position*/
			}
		}
		return is_pre_proline;
	}



} //working_parameters
} //modeler
} //stepwise
} //protocols
