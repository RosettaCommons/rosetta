// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/scoring/func/USOGFunc.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Package headers
#include <core/scoring/func/USOGFunc.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/constraints/util.hh>

// Project headers
#include <core/types.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>

// C/C++ headers
#include <algorithm>
#include <cmath>
#include <iostream>

#define SQRT_2PI 2.50662721600161

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace func {

Real USOGFunc::background_prob = exp(-10.); // the maximum constraint penalty == -log(background_prob)

USOGFunc::USOGFunc(const utility::vector1<core::Real>& means,
	const utility::vector1<core::Real>& std_devs,
	const utility::vector1<core::Real>& weights)
: means_(means), std_devs_(std_devs), weights_(weights) {
	if ( means_.size() != std_devs_.size() || means_.size() != weights_.size() ) {
		utility_exit_with_message("Unequal number of means, std_devs, weights");
	}
}

USOGFunc::USOGFunc(core::Real mean, core::Real std_dev, core::Real weight) {
	means_.push_back( mean );
	std_devs_.push_back( std_dev );
	weights_.push_back( weight );
}

FuncOP USOGFunc::clone() const {
	return FuncOP( new USOGFunc(*this) );
}

bool USOGFunc::operator == ( Func const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	USOGFunc const & other_downcast( static_cast< USOGFunc const & > (other) );
	if ( means_    != other_downcast.means_    ) return false;
	if ( std_devs_ != other_downcast.std_devs_ ) return false;
	if ( weights_  != other_downcast.weights_  ) return false;

	return true;
}

bool USOGFunc::same_type_as_me( Func const & other ) const
{
	return dynamic_cast< USOGFunc const * > ( &other );
}


core::Real USOGFunc::func(const core::Real x) const {
	Real score = 0;
	for ( core::Size i = 1; i <= numGaussians(); ++i ) {
		Real Z = (x - means_[i])/std_devs_[i];
		score += ( weights_[i] / (std_devs_[i] * SQRT_2PI)) * exp( -0.5 * Z * Z );
	}
	return -std::log(score + background_prob);
}

core::Real USOGFunc::dfunc(const core::Real x) const {
	Real score = 0, dscore = 0;
	for ( core::Size i = 1; i <= numGaussians(); ++i ) {
		Real Z = (x - means_[i])/std_devs_[i];
		Real score_i = ( weights_[i] / (std_devs_[i] * SQRT_2PI)) * exp( -0.5 * Z * Z );
		score += score_i;
		dscore += score * (-Z/std_devs_[i]);
	}
	return (-dscore / (score + background_prob));
}

core::Size USOGFunc::numGaussians() const {
	return means_.size();
}

void USOGFunc::show_definition(std::ostream& out) const {
	out << "USOGFUNC " << numGaussians();
	for ( core::Size i = 1; i <= numGaussians(); ++i ) {
		out << " " << means_[i] << " " << std_devs_[i] << " " << weights_[i];
	}
}

void USOGFunc::resetInstance() {
	means_.clear();
	weights_.clear();
	std_devs_.clear();
}

// Format: USOGFunc <num_gaussians> <mean_1> <std_dev1> <weight_1> ... <mean_n> <std_dev_n> <weight_n>
void USOGFunc::read_data(std::istream& in) {
	resetInstance();

	core::Size num_gaussians = 0;
	in >> num_gaussians;

	for ( core::Size i = 1; i <= num_gaussians; ++i ) {
		core::Real m = readValueOrDie(in);
		core::Real s = readValueOrDie(in);
		core::Real w = readValueOrDie(in);

		debug_assert(s > 0);
		debug_assert(w > 0);

		means_.push_back(m);
		weights_.push_back(w);
		std_devs_.push_back(s);
	}
}

// Utility functions
core::Real readValueOrDie(std::istream& in) {
	core::Real x;
	in >> x;

	if ( in.fail() ) {
		utility_exit_with_message("Failed to read floating point value ");
	}

	return x;
}

}  // namespace constraints
}  // namespace scoring
}  // namespace core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::func::USOGFunc::save( Archive & arc ) const {
	arc( cereal::base_class< Func >( this ) );
	arc( CEREAL_NVP( means_ ) ); // utility::vector1<core::Real>
	arc( CEREAL_NVP( std_devs_ ) ); // utility::vector1<core::Real>
	arc( CEREAL_NVP( weights_ ) ); // utility::vector1<core::Real>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::func::USOGFunc::load( Archive & arc ) {
	arc( cereal::base_class< Func >( this ) );
	arc( means_ ); // utility::vector1<core::Real>
	arc( std_devs_ ); // utility::vector1<core::Real>
	arc( weights_ ); // utility::vector1<core::Real>
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::func::USOGFunc );
CEREAL_REGISTER_TYPE( core::scoring::func::USOGFunc )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_func_USOGFunc )
#endif // SERIALIZATION
