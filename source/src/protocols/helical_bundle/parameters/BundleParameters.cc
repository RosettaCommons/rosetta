// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/helical_bundle/parameters/BundleParameters.cc
/// @brief  A class for holding parameters for parametric helical bundle generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit header
#include <protocols/helical_bundle/parameters/BundleParameters.hh>

// Package headers
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.hh>

// Project headers

// Basic headers
#include <basic/basic.hh>
#include <basic/Tracer.hh>

// Numeric headers

// Utility Headers
#include <utility/assert.hh>
#include <utility/py/PyAssert.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <iostream>


namespace protocols {
	namespace helical_bundle {
		namespace parameters {

			static thread_local basic::Tracer TR( "protocols.helical_bundle.parameters.BundleParameters" );

			/// @brief Constructor.
			///
			BundleParameters::BundleParameters() :
				r0_(0.0),
				omega0_(0.0),
				delta_omega0_(0.0),
				r1_(),
				omega1_(0.0),
				delta_omega1_all_(0.0),
				z1_(0.0),
				delta_omega1_(),
				delta_z1_(),
				invert_helix_(false),
				delta_t_(0.0),
				allow_dihedrals_(true),
				allow_bondangles_(false),
				allow_bondlengths_(false)
			{
			}

			BundleParameters::BundleParameters( BundleParameters const & src ) :
				core::conformation::parametric::Parameters( src ),
				r0_(src.r0()),
				omega0_(src.omega0()),
				delta_omega0_(src.delta_omega0()),
				r1_(src.r1_),
				omega1_(src.omega1()),
				delta_omega1_all_(src.delta_omega1_all()),
				z1_(src.z1()),
				delta_omega1_(src.delta_omega1_),
				delta_z1_(src.delta_z1_),
				invert_helix_(src.invert_helix()),
				delta_t_(src.delta_t()),
				allow_dihedrals_(src.allow_dihedrals()),
				allow_bondangles_(src.allow_bondangles()),
				allow_bondlengths_(src.allow_bondlengths())
			{
			}

			BundleParameters::~BundleParameters() {}


			///@brief make a copy of this residue( allocate actual memory for it )
			///
			core::conformation::parametric::ParametersOP
			BundleParameters::clone() const
			{
				return ParametersOP( new BundleParameters( *this ) );
			}

		} // namespace parameters
	} // namespace helical_bundle
} // namespace protocols

