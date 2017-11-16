// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerValueRange.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerValueRange.hh>
#include <basic/Tracer.hh>
#include <numeric/NumericTraits.hh>

using namespace core;

static basic::Tracer TR( "protocols.sampler.rigid_body.RigidBodyStepWiseSamplerValueRange" );
static Real const RADS_PER_DEG = numeric::NumericTraits < Real > ::pi() / 180.;

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rigid_body {

//Constructor
RigidBodyStepWiseSamplerValueRange::RigidBodyStepWiseSamplerValueRange():
	centroid_bin_min_( 0 ),
	centroid_bin_max_( 0 ),
	centroid_bin_size_( 0.0 ),
	euler_angle_bin_min_( 0 ),
	euler_angle_bin_max_( 0 ),
	euler_angle_bin_size_( 0.0 ),
	euler_z_bin_min_( 0 ),
	euler_z_bin_max_( 0 ),
	euler_z_bin_size_( 0.0 )
{}

//Destructor
RigidBodyStepWiseSamplerValueRange::~RigidBodyStepWiseSamplerValueRange()
{}

////////////////////////////////////////////////////////////////////////////////////////
void
RigidBodyStepWiseSamplerValueRange::init(){
	Real const max_angle_rounded = int( 180 / euler_angle_bin_size_ ) * euler_angle_bin_size_;
	set_euler_alpha_values( -max_angle_rounded + 0.5 * Real ( euler_angle_bin_size_ ),
		+max_angle_rounded - 0.5 * Real ( euler_angle_bin_size_ ),
		euler_angle_bin_size_ );

	Real const max_euler_z_rounded = int( 1.0 / euler_z_bin_size_ ) * euler_z_bin_size_;
	set_euler_z_values( -max_euler_z_rounded,
		+max_euler_z_rounded,
		euler_z_bin_size_ );
	set_euler_gamma_values( -max_angle_rounded + 0.5 * Real ( euler_angle_bin_size_ ),
		+max_angle_rounded - 0.5 * Real ( euler_angle_bin_size_ ),
		euler_angle_bin_size_ );

	////////////////////////////////////////////////////////////////////////////////
	Real max_distance_rounded = int( max_distance_/centroid_bin_size_ ) * centroid_bin_size_;
	set_x_values( -max_distance_rounded + 0.5 * Real(centroid_bin_size_),
		+max_distance_rounded - 0.5 * Real(centroid_bin_size_),
		centroid_bin_size_);
	set_y_values( -max_distance_rounded + 0.5 * Real(centroid_bin_size_),
		+max_distance_rounded - 0.5 * Real(centroid_bin_size_),
		centroid_bin_size_);

	// weird offset -- just matching Parin's original
	set_z_values( -max_distance_rounded,
		+max_distance_rounded - Real(centroid_bin_size_),
		centroid_bin_size_);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RigidBodyStepWiseSamplerValueRange::set_sampler_values( Real const & val_min, Real const & val_max, Real const & val_bin, utility::vector1< Real > & values ){
	values.clear();
	runtime_assert( val_min < val_max );
	for ( Real val = val_min; val <= val_max + 1.0e-6 /*floating point error*/; val += val_bin ) values.push_back( val );
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RigidBodyStepWiseSamplerValueRange::set_x_values( Real const centroid_x_min, Real const centroid_x_max, Real const centroid_x_bin ){
	set_sampler_values( centroid_x_min, centroid_x_max, centroid_x_bin, x_values_ );
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RigidBodyStepWiseSamplerValueRange::set_y_values( Real const centroid_y_min, Real const centroid_y_max, Real const centroid_y_bin ){
	set_sampler_values( centroid_y_min, centroid_y_max, centroid_y_bin, y_values_ );
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RigidBodyStepWiseSamplerValueRange::set_z_values( Real const centroid_z_min, Real const centroid_z_max, Real const centroid_z_bin ){
	set_sampler_values( centroid_z_min, centroid_z_max, centroid_z_bin, z_values_ );
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RigidBodyStepWiseSamplerValueRange::set_euler_alpha_values( Real const centroid_euler_alpha_min, Real const centroid_euler_alpha_max, Real const centroid_euler_alpha_bin ){
	set_sampler_values( RADS_PER_DEG * centroid_euler_alpha_min, RADS_PER_DEG * centroid_euler_alpha_max, RADS_PER_DEG * centroid_euler_alpha_bin, euler_alpha_values_ );
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RigidBodyStepWiseSamplerValueRange::set_euler_z_values( Real const centroid_euler_z_min, Real const centroid_euler_z_max, Real const centroid_euler_z_bin ){
	set_sampler_values( centroid_euler_z_min, centroid_euler_z_max, centroid_euler_z_bin, euler_z_values_ );
	// need to really get these in bounds...
	for ( Size n = 1; n <= euler_z_values_.size(); n++ ) {
		euler_z_values_[n] = std::max( euler_z_values_[n], -1.0 );
		euler_z_values_[n] = std::min( euler_z_values_[n], +1.0 );
	}
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RigidBodyStepWiseSamplerValueRange::set_euler_gamma_values( Real const centroid_euler_gamma_min, Real const centroid_euler_gamma_max, Real const centroid_euler_gamma_bin ){
	set_sampler_values( RADS_PER_DEG * centroid_euler_gamma_min, RADS_PER_DEG * centroid_euler_gamma_max, RADS_PER_DEG * centroid_euler_gamma_bin, euler_gamma_values_ );
}

} //rigid_body
} //sampler
} //stepwise
} //protocols
