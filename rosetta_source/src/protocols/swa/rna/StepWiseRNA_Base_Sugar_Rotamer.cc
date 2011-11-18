// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_Base_Sugar_Rotamer
/// @brief Parameters to be passed between different modules of stepwise RNA building.
/// @detailed
/// @author Parin Sripakdeevong


//////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_Base_Sugar_Rotamer.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rna/RNA_FittedTorsionInfo.hh>
#include <core/scoring/rna/RNA_Util.hh>
#include <basic/Tracer.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <string>

using namespace core;
using core::Real;
using ObjexxFCL::fmt::F;

static basic::Tracer TR( "protocols.swa.rna.stepwise_rna_base_sugar_rotamer" );

namespace protocols {
namespace swa {
namespace rna {

	//////////////////////////////////////////////////////////////////////////////////////////
	// Constructor
	StepWiseRNA_Base_Sugar_Rotamer::StepWiseRNA_Base_Sugar_Rotamer(
											BaseState const & base_state, 
											PuckerState const & pucker_state, 
											core::scoring::rna::RNA_FittedTorsionInfo const & rna_fitted_torsion_info,
											core::Size const bin_size, 
											core::Size const bins4): //This is the bin for chi and epsilon, these two torsion angles vary from -20+mean to 20+mean
		base_state_(base_state),
		pucker_state_(pucker_state),
		rna_fitted_torsion_info_(rna_fitted_torsion_info),
		bin_size_( bin_size ), // must be 20, 10, or 5
		num_base_std_ID_( bins4 ) 
	{
		reset();

		if(base_state_== BOTH){
			num_base_ID_= 2; 
		}else if(base_state_== ANTI){
			num_base_ID_= 1;
		}else if(base_state == NONE){
			num_base_ID_= 0;
		}else{
			utility_exit_with_message( "Invalid pucker_state_" );
		}

		pucker_state_list_.clear();
		if(pucker_state_ == ALL){
			pucker_state_list_.push_back(NORTH);
			pucker_state_list_.push_back(SOUTH);
		}else if(pucker_state_ == NORTH){
			pucker_state_list_.push_back(NORTH);
		}else if(pucker_state_ == SOUTH ){
			pucker_state_list_.push_back(SOUTH);
		}else{
			utility_exit_with_message( "Invalid pucker_state_" );
		}

		if(num_base_std_ID_ != (1 + 40/bin_size_) ) utility_exit_with_message( "chi_bin_ != (1 + 40/bin_size_)" );
	}

  //////////////////////////////////////////////////////////////////////////
  //destructor
  StepWiseRNA_Base_Sugar_Rotamer::~StepWiseRNA_Base_Sugar_Rotamer()
  {}

  //////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_Base_Sugar_Rotamer::get_next_rotamer(){

		if(pucker_ID_ > pucker_state_list_.size()) return false;

		//This pucker_ID_old is always consistent the stored chi_ value
		pucker_ID_old_=pucker_ID_;

		PuckerState const & curr_pucker_state=current_pucker_state();

		if (curr_pucker_state == NORTH) {
			delta_ = rna_fitted_torsion_info_.ideal_delta_north();
			chi_ = rna_fitted_torsion_info_.gaussian_parameter_set_chi_north()[ base_ID_ ].center + bin_size_*(base_std_ID_-1) - 20;
			nu2_ = rna_fitted_torsion_info_.ideal_nu2_north();
			nu1_ = rna_fitted_torsion_info_.ideal_nu1_north();
		}	else if(curr_pucker_state == SOUTH) {
			delta_ = rna_fitted_torsion_info_.ideal_delta_south();
			chi_ = rna_fitted_torsion_info_.gaussian_parameter_set_chi_south()[ base_ID_ ].center + bin_size_*(base_std_ID_-1) - 20;
			nu2_ = rna_fitted_torsion_info_.ideal_nu2_south();
			nu1_ = rna_fitted_torsion_info_.ideal_nu1_south();
		}else{
			utility_exit_with_message( "Invalid current_pucker_state!" );
		}

		base_std_ID_++;

		if(base_std_ID_ > num_base_std_ID_){
			base_std_ID_=1;
			base_ID_++;
		}

		if(base_ID_ > num_base_ID_){
			base_ID_ = 1;
			pucker_ID_++;
		}

		return true;
	}

  //////////////////////////////////////////////////////////////////////////
	PuckerState const & 
	StepWiseRNA_Base_Sugar_Rotamer::current_pucker_state() const { 

//		std::cout << "pucker_ID_old_= " << pucker_ID_old_ << std::endl;		
//		std::cout << "pucker_state_list_.size()= " << pucker_state_list_.size() << std::endl;		

		PuckerState const & pucker_state=pucker_state_list_[pucker_ID_old_];
		if(pucker_state==ALL){
			utility_exit_with_message( "pucker_state should not equal ALL!" );
		}
		return pucker_state;
	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_Base_Sugar_Rotamer::reset() { 

		pucker_ID_=1;
		pucker_ID_old_=pucker_ID_;
		base_ID_=1;
		base_std_ID_=1;

		chi_=0.0;
		delta_=0.0;
		nu2_=0.0;
		nu1_=0.0;

	}


}
}
}
