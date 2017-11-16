// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////
///
/// @file protocols/scoring/methods/pcs2/PcsGridSearchParameter.cc
///
/// @brief
///
/// @details
///
/// @param
///
/// @return
///
/// @remarks
///
/// @references
///
/// @authorv Christophe Schmitz
///
////////////////////////////////////////////////


// Unit Headers
#include <protocols/scoring/methods/pcs2/PcsGridSearchParameter.hh>

// Package Headers

// Project Headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh>

// ObjexxFCL Headers

// C++ headers
#include <iostream>

#include <utility/vector1.hh>


static basic::Tracer TR_PcsGridSearchParameter( "protocols.scoring.methods.pcs.PcsGridSearchParameter" );


namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {

using namespace protocols;

PcsGridSearchParameter::~PcsGridSearchParameter(){
}

PcsGridSearchParameter &
PcsGridSearchParameter::operator=(PcsGridSearchParameter const & other){
	if ( this != &other ) {
		include_only_start_stage1_ = other.include_only_start_stage1_;
		include_only_start_stage2_ = other.include_only_start_stage2_;
		include_only_start_stage3_ = other.include_only_start_stage3_;
		include_only_start_stage4_ = other.include_only_start_stage4_;

		include_only_end_stage1_ = other.include_only_end_stage1_;
		include_only_end_stage2_ = other.include_only_end_stage2_;
		include_only_end_stage3_ = other.include_only_end_stage3_;
		include_only_end_stage4_ = other.include_only_end_stage4_;

		n_trial_min_stage1_ = other.n_trial_min_stage1_;
		n_trial_min_stage2_ = other.n_trial_min_stage2_;
		n_trial_min_stage3_ = other.n_trial_min_stage3_;
		n_trial_min_stage4_ = other.n_trial_min_stage4_;

		pcs_weight_stage1_ = other.pcs_weight_stage1_;
		pcs_weight_stage2_ = other.pcs_weight_stage2_;
		pcs_weight_stage3_ = other.pcs_weight_stage3_;
		pcs_weight_stage4_ = other.pcs_weight_stage4_;

		individual_scale_stage1_ = other.individual_scale_stage1_;
		individual_scale_stage2_ = other.individual_scale_stage2_;
		individual_scale_stage3_ = other.individual_scale_stage3_;
		individual_scale_stage4_ = other.individual_scale_stage4_;

		filenames_ = other.filenames_;
		individual_weights_ = other.individual_weights_;

	}
	return *this;
}


PcsGridSearchParameter::PcsGridSearchParameter(PcsGridSearchParameter const & other){

	include_only_start_stage1_ = other.include_only_start_stage1_;
	include_only_start_stage2_ = other.include_only_start_stage2_;
	include_only_start_stage3_ = other.include_only_start_stage3_;
	include_only_start_stage4_ = other.include_only_start_stage4_;

	include_only_end_stage1_ = other.include_only_end_stage1_;
	include_only_end_stage2_ = other.include_only_end_stage2_;
	include_only_end_stage3_ = other.include_only_end_stage3_;
	include_only_end_stage4_ = other.include_only_end_stage4_;

	n_trial_min_stage1_ = other.n_trial_min_stage1_;
	n_trial_min_stage2_ = other.n_trial_min_stage2_;
	n_trial_min_stage3_ = other.n_trial_min_stage3_;
	n_trial_min_stage4_ = other.n_trial_min_stage4_;

	pcs_weight_stage1_ = other.pcs_weight_stage1_;
	pcs_weight_stage2_ = other.pcs_weight_stage2_;
	pcs_weight_stage3_ = other.pcs_weight_stage3_;
	pcs_weight_stage4_ = other.pcs_weight_stage4_;

	individual_scale_stage1_ = other.individual_scale_stage1_;
	individual_scale_stage2_ = other.individual_scale_stage2_;
	individual_scale_stage3_ = other.individual_scale_stage3_;
	individual_scale_stage4_ = other.individual_scale_stage4_;


	filenames_ = other.filenames_;
	individual_weights_ = other.individual_weights_;
}


PcsGridSearchParameter::PcsGridSearchParameter(){

	include_only_start_stage1_ = 0;
	include_only_start_stage2_ = 0;
	include_only_start_stage3_ = 0;
	include_only_start_stage4_ = 0;

	include_only_end_stage1_ = 0;
	include_only_end_stage2_ = 0;
	include_only_end_stage3_ = 0;
	include_only_end_stage4_ = 0;

	n_trial_min_stage1_ = 0;
	n_trial_min_stage2_ = 0;
	n_trial_min_stage3_ = 0;
	n_trial_min_stage4_ = 0;

	pcs_weight_stage1_ = 1.0;
	pcs_weight_stage2_ = 1.0;
	pcs_weight_stage3_ = 1.0;
	pcs_weight_stage4_ = 1.0;

	individual_scale_stage1_ = -1.0;
	individual_scale_stage2_ = -1.0;
	individual_scale_stage3_ = -1.0;
	individual_scale_stage4_ = -1.0;

}


std::ostream &
operator << ( std::ostream& out, const PcsGridSearchParameter &me ){
	out << "pcs_weight_stage[1->4] " <<
		me.pcs_weight_stage1_ << " " <<
		me.pcs_weight_stage2_ << " " <<
		me.pcs_weight_stage3_ << " " <<
		me.pcs_weight_stage4_ << " " <<std::endl;

	out << "include_only_start_stage[1->4] " <<
		me.include_only_start_stage1_ << " " <<
		me.include_only_start_stage2_ << " " <<
		me.include_only_start_stage3_ << " " <<
		me.include_only_start_stage4_ << " " <<std::endl;

	out << "include_only_end_stage[1->4] " <<
		me.include_only_end_stage1_ << " " <<
		me.include_only_end_stage2_ << " " <<
		me.include_only_end_stage3_ << " " <<
		me.include_only_end_stage4_ << " " <<std::endl;

	out << "n_trial_min_stage[1->4] " <<
		me.n_trial_min_stage1_ << " " <<
		me.n_trial_min_stage2_ << " " <<
		me.n_trial_min_stage3_ << " " <<
		me.n_trial_min_stage4_ << " " <<std::endl;

	out << "individual_scale[1->4] " <<
		me.individual_scale_stage1_<< " " <<
		me.individual_scale_stage2_<< " " <<
		me.individual_scale_stage3_<< " " <<
		me.individual_scale_stage4_<< " " <<std::endl;

	out <<"filenames / individual_weights" << std::endl;
	core::Size i;
	for ( i = 1; i <= me.filenames_.size(); ++i ) {
		out <<"   "<< me.filenames_[i] << " " << me.individual_weights_[i]  << std::endl;
	}

	return out;
}

void
output_error_message(core::Size id){
	std::cerr << "In the definition of the grid search parameters for the PseudocontactShift, you should ensure that for all stages:" << std::endl;
	std::cerr << "1. the PCS_WEIGHT is positif" << std::endl;
	std::cerr << "2. the weight of each lanthanide is positif" << std::endl;
	std::cerr << "3. If the INCLUDE_ONLY tag is used: 0 < residue_start < residue_end" << std::endl;
	std::cerr << "4. If the INDIVIDUAL_SCALE tag is used: its value must be positive" << std::endl;
	std::cerr << "Looks like the problem is " << id << std::endl;
	utility_exit_with_message("Please, review your PseudocontactShift parameter file and correct it");
}

bool
valide_include_only_stage(core::Size start, core::Size end){
	if ( (start == 0) && (end == 0) ) {
		return true;
	}
	if ( start == 0 ) {
		return false;
	}
	if ( start >= end ) {
		return false;
	}
	return true;
}

void
PcsGridSearchParameter::control_grid_param(){

	core::Size i;
	if (
			(pcs_weight_stage1_ < 0)||
			(pcs_weight_stage2_ < 0)||
			(pcs_weight_stage3_ < 0)||
			(pcs_weight_stage4_ < 0)
			) {
		output_error_message(1);
	}
	if ( valide_include_only_stage(include_only_start_stage1_, include_only_end_stage1_) == false ) {
		output_error_message(3);
	}
	if ( valide_include_only_stage(include_only_start_stage2_, include_only_end_stage2_) == false ) {
		output_error_message(3);
	}
	if ( valide_include_only_stage(include_only_start_stage3_, include_only_end_stage3_) == false ) {
		output_error_message(3);
	}
	if ( valide_include_only_stage(include_only_start_stage4_, include_only_end_stage4_) == false ) {
		output_error_message(3);
	}
	if (
			(individual_scale_stage1_ < 0) ||
			(individual_scale_stage2_ < 0) ||
			(individual_scale_stage3_ < 0) ||
			(individual_scale_stage4_ < 0)
			) {
		if (
				(individual_scale_stage1_ != -1.0) &&
				(individual_scale_stage2_ != -1.0) &&
				(individual_scale_stage3_ != -1.0) &&
				(individual_scale_stage4_ != -1.0)
				) {
			output_error_message(4);
		}
	}

	for ( i = 1; i <= individual_weights_.size(); i++ ) {
		if ( individual_weights_[i] < 0 ) {
			output_error_message(2);
		}
	}

}

}//namespace pcs2
}//namespace methods
}//namespace scoring
}//namespace protocols
