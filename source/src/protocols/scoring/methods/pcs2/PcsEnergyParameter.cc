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
/// @file protocols/scoring/methods/pcs2/PcsEnergyParameter.cc
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


// Unit headers
#include <protocols/scoring/methods/pcs2/PcsEnergyParameter.hh>

// Package headers

// Project headers
#include <basic/Tracer.hh>

// Utility headers

// Numeric headers

// Objexx headers

// C++ headers

#include <utility/vector1.hh>


namespace protocols {
namespace scoring {
namespace methods {
namespace pcs2 {

static THREAD_LOCAL basic::Tracer TR_PcsEnergyParameter( "protocols.scoring.methods.pcs.PcsEnergyParameter" );

PcsEnergyParameter::PcsEnergyParameter(){

}

PcsEnergyParameter::~PcsEnergyParameter(){
}

PcsEnergyParameter &
PcsEnergyParameter::operator=(PcsEnergyParameter const & other){
	if ( this != &other ) {
		include_only_start_ = other.include_only_start_;
		include_only_end_ = other.include_only_end_;
		n_trial_min_ = other.n_trial_min_;
		pcs_weight_ = other.pcs_weight_;
		individual_scale_ = other.individual_scale_;
		vec_filename_ = other.vec_filename_;
		vec_individual_weight_ = other.vec_individual_weight_;
	}
	return *this;
}

PcsEnergyParameter::PcsEnergyParameter(PcsEnergyParameter const & other){

	include_only_start_ = other.include_only_start_;
	include_only_end_ = other.include_only_end_;
	n_trial_min_ = other.n_trial_min_;
	pcs_weight_ = other.pcs_weight_;
	individual_scale_ = other.individual_scale_;
	vec_filename_ = other.vec_filename_;
	vec_individual_weight_ = other.vec_individual_weight_;
}

void
PcsEnergyParameter::set_vector_name_and_weight(utility::vector1<std::string> const vec_filename,
	utility::vector1<core::Real> const vec_individual_weight){

	vec_filename_ = vec_filename;
	vec_individual_weight_ = vec_individual_weight;
}

void
PcsEnergyParameter::set_grid_param(
	core::Size const include_only_start,
	core::Size const include_only_end,
	core::Size const n_trial_min,
	core::Real const pcs_weight,
	core::Real const individual_scale
){

	include_only_start_ = include_only_start;
	include_only_end_ = include_only_end;
	n_trial_min_ = n_trial_min;
	pcs_weight_ = pcs_weight;
	individual_scale_ = individual_scale;
}

std::ostream &
operator<<(std::ostream& out, const PcsEnergyParameter &me){

	out <<"include_only_start " <<
		me.include_only_start_ << " include_only_end " <<
		me.include_only_end_ <<" n_trial_min " <<
		me.n_trial_min_ <<" pcs_weight " <<
		me.pcs_weight_<< std::endl;
	return out;
}


core::Size
PcsEnergyParameter::get_include_only_start() const{
	return include_only_start_;
}


core::Size
PcsEnergyParameter::get_include_only_end() const{
	return include_only_end_;
}

core::Size
PcsEnergyParameter::get_n_trial_min() const{
	return n_trial_min_;
}

core::Real
PcsEnergyParameter::get_pcs_weight() const{
	return pcs_weight_;
}


core::Real
PcsEnergyParameter::get_individual_scale() const{
	return individual_scale_;
}


utility::vector1<std::string> const &
PcsEnergyParameter::get_vector_filename() const{
	return vec_filename_;
}

utility::vector1<core::Real> const &
PcsEnergyParameter::get_vector_weight() const{
	return vec_individual_weight_;
}

} // PCS
} // methods
} // scoring
} // core
