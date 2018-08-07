// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/helical_bundle/parameters/BundleParameters.cc
/// @brief  A class for holding parameters for parametric helical bundle generation.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

// Unit header
#include <protocols/helical_bundle/parameters/BundleParameters.hh>

// Package headers
#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/ParametersSet.hh>

// Project headers
#include <protocols/helical_bundle/BundleParametrizationCalculator.hh>
#include <core/conformation/parametric/BooleanValuedParameter.hh>
#include <core/conformation/parametric/RealValuedParameter.hh>
#include <core/conformation/parametric/SizeValuedParameter.hh>
#include <core/conformation/parametric/RealVectorValuedParameter.hh>
#include <core/conformation/parametric/SizeVectorValuedParameter.hh>

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


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace helical_bundle {
namespace parameters {

static basic::Tracer TR( "protocols.helical_bundle.parameters.BundleParameters" );

/// @brief Constructor.
///
BundleParameters::BundleParameters() :
	core::conformation::parametric::Parameters()
{}

BundleParameters::BundleParameters( BundleParameters const & src ) :
	core::conformation::parametric::Parameters( src )
{}

BundleParameters::~BundleParameters() = default;


/// @brief make a copy of this residue( allocate actual memory for it )
///
core::conformation::parametric::ParametersOP
BundleParameters::clone() const
{
	return ParametersOP( new BundleParameters( *this ) );
}

/// @brief Get a summary of this ParametersSet object, for output to remark lines of a PDB file.
///
void BundleParameters::get_pdb_remark(std::stringstream &remark) const {
	using namespace core::conformation::parametric;
	using namespace protocols::helical_bundle;

	remark.setf( std::ios::fixed, std::ios::floatfield );
	remark.precision(8);
	remark << " PARAMETERS THAT ARE TYPICALLY SAMPLED:" << std::endl;
	for ( core::Size i(1); i < static_cast<core::Size>( BPC_end_of_list ); ++i ) {
		ParameterCOP curparam( parameter_cop(i) );
		remark << "   " << curparam->short_parameter_description() << " (" << curparam->parameter_name() << "," << curparam->parameter_units() << "): ";
		//Determine type:
		ParameterType const paramtype( curparam->parameter_type() );
		RealValuedParameterCOP realparam( utility::pointer::static_pointer_cast< RealValuedParameter const>( curparam ) );
		BooleanValuedParameterCOP boolparam( utility::pointer::static_pointer_cast< BooleanValuedParameter const>( curparam ) );
		SizeValuedParameterCOP sizeparam( utility::pointer::static_pointer_cast< SizeValuedParameter const>( curparam ) );
		SizeVectorValuedParameterCOP sizevectparam( utility::pointer::static_pointer_cast< SizeVectorValuedParameter const>( curparam ) );
		RealVectorValuedParameterCOP realvectparam( utility::pointer::static_pointer_cast< RealVectorValuedParameter const>( curparam ) );
		switch( paramtype ) {
		case PT_generic_real:
		case PT_generic_nonnegative_valued_real:
		case PT_generic_positive_valued_real:
		case PT_angle :
			remark << realparam->value();
			break;
		case PT_boolean :
			remark << ( boolparam->value() ? "TRUE" : "FALSE" );
			break;
		case PT_generic_integer:
		case PT_generic_whole_number:
		case PT_generic_natural_number :
			remark << sizeparam->value();
			break;
		case PT_generic_integer_vector:
		case PT_generic_whole_number_vector:
		case PT_generic_natural_number_vector :
			remark << std::endl << "     " << sizevectparam->value();
			break;
		case PT_generic_real_vector:
		case PT_generic_nonnegative_valued_real_vector:
		case PT_generic_positive_valued_real_vector:
		case PT_angle_vector :
			remark << std::endl << "     " << realvectparam->value();
			break;
		case PT_invalid_type :
			remark << std::endl;
			utility_exit_with_message( "Error in BundleParameters::get_pdb_remark(): Could not determine type of parameter " + curparam->parameter_name() + "!" );
		}

		remark << std::endl;
		if ( i == static_cast<core::Size>( BPC_last_parameter_to_be_sampled ) ) {
			remark << " PARAMETERS THAT ARE NOT TYPICALLY SAMPLED:" << std::endl;
		}
	}
}

} // namespace parameters
} // namespace helical_bundle
} // namespace protocols


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::helical_bundle::parameters::BundleParameters::save( Archive & arc ) const {
	arc( cereal::base_class< core::conformation::parametric::Parameters >( this ) );
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::helical_bundle::parameters::BundleParameters::load( Archive & arc ) {
	arc( cereal::base_class< core::conformation::parametric::Parameters >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::helical_bundle::parameters::BundleParameters );
CEREAL_REGISTER_TYPE( protocols::helical_bundle::parameters::BundleParameters )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_helical_bundle_parameters_BundleParameters )
#endif // SERIALIZATION
