// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/helical_bundle/PerturbBundleOptions.cc
/// @brief  A class for options for the PerturbBundle mover.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit header
#include <protocols/helical_bundle/PerturbBundleOptions.hh>

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


namespace protocols {
namespace helical_bundle {

static THREAD_LOCAL basic::Tracer TR( "protocols.helical_bundle.PerturbBundleOptions" );

/// @brief Constructor.
///
PerturbBundleOptions::PerturbBundleOptions() :
	helix_index_(0),
	perturbation_magnitude_(0.0),
	perturbation_type_(pt_gaussian),
	perturbable_(false),
	being_set_(false),
	use_defaults_(true),
	use_value_from_other_helix_(0),
	omega0_copies_pitch_instead_(false),
	samples_(0),
	default_value_(0),
	lower_value_(0),
	upper_value_(0)
{
}

PerturbBundleOptions::PerturbBundleOptions( PerturbBundleOptions const & src ) :
	utility::pointer::ReferenceCount(),
	utility::pointer::enable_shared_from_this< PerturbBundleOptions >() ,
	helix_index_(src.helix_index_),
	perturbation_magnitude_(src.perturbation_magnitude_),
	perturbation_type_(src.perturbation_type_),
	perturbable_(src.perturbable_),
	being_set_(src.being_set_),
	use_defaults_(src.use_defaults_),
	use_value_from_other_helix_(src.use_value_from_other_helix_),
	omega0_copies_pitch_instead_(src.omega0_copies_pitch_instead_),
	samples_(src.samples_),
	default_value_(src.default_value_),
	lower_value_(src.lower_value_),
	upper_value_(src.upper_value_)
{
}

PerturbBundleOptions::~PerturbBundleOptions() = default;


/// @brief make a copy of this PerturbBundleOptions object (allocate actual memory for it)
///
PerturbBundleOptionsOP
PerturbBundleOptions::clone() const
{
	return PerturbBundleOptionsOP( new PerturbBundleOptions( *this ) );
}

/// @brief Function to generate a random value, based on the perturbation type and perturbation magnitude.
/// @details Typically, one would then add this value to the current parameter value to perturb it.
core::Real PerturbBundleOptions::delta() const {
	if ( perturbation_type()==pt_uniform ) {
		return ((numeric::random::rg().uniform()*2.0)-1.0) * perturbation_magnitude();
	} else if ( perturbation_type()==pt_gaussian ) {
		return numeric::random::rg().gaussian() * perturbation_magnitude();
	} else {
		if ( TR.Warning.visible() ) TR.Warning << "An unknown perturber was specified.  Returning 0.0 from protocols::helical_bundle::PerturbBundleOptions::delta() function." << std::endl;
	}
	return 0.0;
}

} // namespace helical_bundle
} // namespace protocols

