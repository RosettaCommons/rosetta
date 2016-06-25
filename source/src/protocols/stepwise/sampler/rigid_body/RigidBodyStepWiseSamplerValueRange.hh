// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerValueRange.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_sampler_rigid_body_RigidBodyStepWiseSamplerValueRange_HH
#define INCLUDED_protocols_sampler_rigid_body_RigidBodyStepWiseSamplerValueRange_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerOneValue.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerValueRange.fwd.hh>
#include <core/types.hh>

// To Author(s) of this code: our coding convention explicitly forbid of using ‘using namespace ...’ in header files outside class or function body, please make sure to refactor this out!
using namespace core;

namespace protocols {
namespace stepwise {
namespace sampler {
namespace rigid_body {

class RigidBodyStepWiseSamplerValueRange: public utility::pointer::ReferenceCount {

public:

	//constructor
	RigidBodyStepWiseSamplerValueRange();

	//destructor
	~RigidBodyStepWiseSamplerValueRange();

public:

	void init();

	void
	set_sampler_values( Real const & val_min, Real const & val_max, Real const & val_bin, ValueList & values );

	void set_x_values( Real const centroid_x_min, Real const centroid_x_max, Real const centroid_x_bin );
	void set_y_values( Real const centroid_y_min, Real const centroid_y_max, Real const centroid_y_bin );
	void set_z_values( Real const centroid_z_min, Real const centroid_z_max, Real const centroid_z_bin );
	void set_euler_alpha_values( Real const centroid_euler_alpha_min, Real const centroid_euler_alpha_max, Real const centroid_euler_alpha_bin );
	void set_euler_z_values( Real const centroid_euler_z_min, Real const centroid_euler_z_max, Real const centroid_euler_z_bin );
	void set_euler_gamma_values( Real const centroid_euler_gamma_min, Real const centroid_euler_gamma_max, Real const centroid_euler_gamma_bin );

	void set_x_values( ValueList const & setting ){ x_values_ = setting; }
	ValueList x_values() const{ return x_values_; }

	void set_y_values( ValueList const & setting ){ y_values_ = setting; }
	ValueList y_values() const{ return y_values_; }

	void set_z_values( ValueList const & setting ){ z_values_ = setting; }
	ValueList z_values() const{ return z_values_; }

	void set_euler_alpha_values( ValueList const & setting ){ euler_alpha_values_ = setting; }
	ValueList euler_alpha_values() const{ return euler_alpha_values_; }

	void set_euler_z_values( ValueList const & setting ){ euler_z_values_ = setting; }
	ValueList euler_z_values() const{ return euler_z_values_; }

	void set_euler_gamma_values( ValueList const & setting ){ euler_gamma_values_ = setting; }
	ValueList euler_gamma_values() const{ return euler_gamma_values_; }

	void set_centroid_bin_min( int const & setting ){ centroid_bin_min_ = setting; }
	int centroid_bin_min() const{ return centroid_bin_min_; }

	void set_centroid_bin_max( int const & setting ){ centroid_bin_max_ = setting; }
	int centroid_bin_max() const{ return centroid_bin_max_; }

	void set_euler_angle_bin_min( int const & setting ){ euler_angle_bin_min_ = setting; }
	int euler_angle_bin_min() const{ return euler_angle_bin_min_; }

	void set_euler_angle_bin_max( int const & setting ){ euler_angle_bin_max_ = setting; }
	int euler_angle_bin_max() const{ return euler_angle_bin_max_; }

	void set_euler_z_bin_min( int const & setting ){ euler_z_bin_min_ = setting; }
	int euler_z_bin_min() const{ return euler_z_bin_min_; }

	void set_euler_z_bin_max( int const & setting ){ euler_z_bin_max_ = setting; }
	int euler_z_bin_max() const{ return euler_z_bin_max_; }

	void set_centroid_bin_size( core::Real const & setting ){ centroid_bin_size_ = setting; }
	core::Real centroid_bin_size() const{ return centroid_bin_size_; }

	void set_euler_angle_bin_size( core::Real const & setting ){ euler_angle_bin_size_ = setting; }
	core::Real euler_angle_bin_size() const{ return euler_angle_bin_size_; }

	void set_euler_z_bin_size( core::Real const & setting ){ euler_z_bin_size_ = setting; }
	core::Real euler_z_bin_size() const{ return euler_z_bin_size_; }

	void set_max_distance( core::Distance const & setting ){ max_distance_ = setting; }
	core::Distance max_distance() const{ return max_distance_; }

private:

	int centroid_bin_min_;
	int centroid_bin_max_;
	core::Real centroid_bin_size_;
	int euler_angle_bin_min_;
	int euler_angle_bin_max_;
	core::Real euler_angle_bin_size_;
	int euler_z_bin_min_;
	int euler_z_bin_max_;
	core::Real euler_z_bin_size_;

	core::Distance max_distance_;

	ValueList x_values_, y_values_, z_values_;
	ValueList euler_alpha_values_, euler_z_values_, euler_gamma_values_;

};

} //rigid_body
} //sampler
} //stepwise
} //protocols

#endif
