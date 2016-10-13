// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/helical_bundle/BundleGridSamplerHelper.cc
/// @brief  A class that stores the various parameter perturbations that will be sampled by the
/// BundleGridSampler mover.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit header
#include <protocols/helical_bundle/BundleGridSamplerHelper.hh>

// Package headers

// Project headers

// Basic headers
#include <basic/basic.hh>
#include <basic/Tracer.hh>

// Numeric headers
#include <numeric/random/random.hh>

// Utility Headers
#include <utility/assert.hh>
#include <utility/py/PyAssert.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <iostream>
#include <cstdio>


namespace protocols {
namespace helical_bundle {

static THREAD_LOCAL basic::Tracer TR( "protocols.helical_bundle.BundleGridSamplerHelper" );

/// @brief Constructor.
///
BundleGridSamplerHelper::BundleGridSamplerHelper() :
	nDoFs_(0),
	allowed_dof_types_(),
	allowed_dof_helix_indices_(),
	dof_samples_(),
	dof_lower_vals_(),
	dof_upper_vals_(),
	dof_sample_vals_(),
	cur_indices_()
{
}

BundleGridSamplerHelper::BundleGridSamplerHelper( BundleGridSamplerHelper const & src ) :
	utility::pointer::ReferenceCount(),
	utility::pointer::enable_shared_from_this< BundleGridSamplerHelper >(),
	nDoFs_( src.nDoFs_ ),
	allowed_dof_types_( src.allowed_dof_types_ ),
	allowed_dof_helix_indices_( src.allowed_dof_helix_indices_ ),
	dof_samples_( src.dof_samples_ ),
	dof_lower_vals_( src.dof_lower_vals_ ),
	dof_upper_vals_( src.dof_upper_vals_ ),
	dof_sample_vals_( src.dof_sample_vals_ ),
	cur_indices_( src.cur_indices_ )
{
}

BundleGridSamplerHelper::~BundleGridSamplerHelper() = default;


/// @brief make a copy of this BundleGridSamplerHelper object (allocate actual memory for it)
///
BundleGridSamplerHelperOP
BundleGridSamplerHelper::clone() const
{
	return BundleGridSamplerHelperOP( new BundleGridSamplerHelper( *this ) );
}

/// @brief Reset this BundleGridSamplerHelper object.
/// @details Clears all internal data.
void BundleGridSamplerHelper::reset() {
	nDoFs_=0;
	allowed_dof_types_.clear();
	allowed_dof_helix_indices_.clear();
	dof_samples_.clear();
	dof_lower_vals_.clear();
	dof_upper_vals_.clear();
	dof_sample_vals_.clear();
	cur_indices_.clear();
	return;
}

/// @brief Perform the pre-calculation that sets up the lists of parameter values to be sampled.
///
void BundleGridSamplerHelper::initialize_samples() {
	if ( TR.visible() ) TR << "Initializing samples." << std::endl;
	dof_sample_vals_.clear();
	cur_indices_.clear();

	for ( core::Size i=1, imax=nDoFs(); i<=imax; ++i ) {

		utility::vector1 < core::Real > curvect;

		if ( TR.Debug.visible() ) { TR.Debug << "lower=" << dof_lower_vals_[i] << " upper=" << dof_upper_vals_[i] << std::endl; TR.Debug.flush();}

		runtime_assert_string_msg( dof_lower_vals_[i] < dof_upper_vals_[i],
			"Internal error in BundleGridSamplerHelper::initialize_samples().  Lower values in parameter ranges must be less than upper values." );
		runtime_assert_string_msg( dof_samples_[i] > 1,
			"Internal error in BundleGridSamplerHelper::initialize_samples().  Sampled parameters must have at least two samples." );
		core::Real const sample_increment = (dof_upper_vals_[i]-dof_lower_vals_[i]) / ( static_cast<core::Real>( dof_samples_[i] - 1 ) );
		core::Real curval = dof_lower_vals_[i];

		for ( core::Size j=1; j<=dof_samples_[i]; ++j ) {
			curvect.push_back(curval);
			curval+=sample_increment;
		}

		dof_sample_vals_.push_back( curvect );
		cur_indices_.push_back(1); //Initialize this vector to (1,1,1,...)
	} //Looping through all DoFs.

	if ( TR.visible() ) {
		TR << "The following parameter values will be sampled:" << std::endl;
		char outstring[1024];
		sprintf(outstring, "Helix\tDoF\tParameter values" );
		TR << outstring << std::endl;

		for ( core::Size i=1, imax=nDoFs(); i<=imax; ++i ) {
			sprintf(outstring, "%lu\t%s\t", (unsigned long) allowed_dof_helix_indices_[i], DoF_name(allowed_dof_types_[i]).c_str() );
			TR << outstring;
			for ( core::Size j=1, jmax=dof_sample_vals_[i].size(); j<=jmax; ++j ) {
				sprintf(outstring, "%.4f\t", dof_sample_vals_[i][j]);
				TR << outstring;
			}
			TR << std::endl;
		}
		TR << std::endl;
		TR.flush();
	}

	return;
}

/// @brief RECURSIVE function that increments the current sample index.
/// @details This function adds 1 to the last DoF index.  If the last DoF
/// index value exceeds the number of samples for that DoF, it resets that
/// DoF index and increments the second-last DoF index by 1 by recursively
/// calling itself.  This function is overloaded; the public version by
/// default tries to increment the last index, while the private version
/// can be called recursively.
void BundleGridSamplerHelper::increment_cur_indices()
{
	increment_cur_indices( nDoFs() ); //Call the private, recursive function.
	return;
}

/// @brief RECURSIVE function that increments the current sample index.
/// @details This function adds 1 to the last DoF index.  If the last DoF
/// index value exceeds the number of samples for that DoF, it resets that
/// DoF index and increments the second-last DoF index by 1 by recursively
/// calling itself.
void BundleGridSamplerHelper::increment_cur_indices( core::Size const index_to_increment )
{
	runtime_assert_string_msg( index_to_increment > 0, "In BundleGridSamplerHelper::increment_cur_indices(): all possible permutations have already been sampled!" );

	++(cur_indices_[index_to_increment]);
	if ( cur_indices_[index_to_increment] > dof_sample_vals_[index_to_increment].size() ) {
		cur_indices_[index_to_increment] = 1;
		increment_cur_indices( index_to_increment - 1 );
	}

	return;
}

/// @brief Return the name of a DoF type given its enum.
///
std::string BundleGridSamplerHelper::DoF_name( DoFType const &type ) const
{
	if ( type == bgsh_r0 ) return "r0";
	else if ( type == bgsh_omega0 ) return "omega0";
	else if ( type == bgsh_delta_omega0 ) return "delta_omega0";
	else if ( type == bgsh_delta_omega1 ) return "delta_omega1";
	else if ( type == bgsh_delta_t ) return "delta_t";
	else if ( type == bgsh_z1_offset ) return "z1_offset";
	else if ( type == bgsh_z0_offset ) return "z0_offset";
	else if ( type == bgsh_epsilon ) return "epsilon";

	return "UNKNOWN";
}

} // namespace helical_bundle
} // namespace protocols

