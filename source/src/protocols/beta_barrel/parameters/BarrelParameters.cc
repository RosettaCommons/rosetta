// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/beta_barrel/parameters/BarrelParameters.cc
/// @brief  Parameters for a single strand in a parametric beta-barrel.
/// @author Andy Watkins

#include <protocols/beta_barrel/parameters/BarrelParameters.hh>

#include <core/conformation/parametric/Parameters.hh>
#include <core/conformation/parametric/BooleanValuedParameter.hh>
#include <core/conformation/parametric/RealValuedParameter.hh>
#include <core/conformation/parametric/SizeValuedParameter.hh>
#include <core/conformation/parametric/RealVectorValuedParameter.hh>
#include <core/conformation/parametric/SizeVectorValuedParameter.hh>

#include <basic/Tracer.hh>
#include <utility/stream_util.hh>
#include <utility/vector1.hh>

#ifdef    SERIALIZATION
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>
#include <cereal/types/polymorphic.hpp>
#endif

namespace protocols {
namespace beta_barrel {
namespace parameters {

static basic::Tracer TR( "protocols.beta_barrel.parameters.BarrelParameters" );

BarrelParameters::BarrelParameters() :
	core::conformation::parametric::Parameters()
{}

BarrelParameters::BarrelParameters( BarrelParameters const & src ) :
	core::conformation::parametric::Parameters( src )
{}

BarrelParameters::~BarrelParameters() = default;

core::conformation::parametric::ParametersOP
BarrelParameters::clone() const
{
	return utility::pointer::make_shared< BarrelParameters >( *this );
}

void BarrelParameters::get_pdb_remark( std::stringstream & remark ) const {
	using namespace core::conformation::parametric;

	remark.setf( std::ios::fixed, std::ios::floatfield );
	remark.precision( 8 );
	for ( core::Size i(1), imax(num_parameters()); i <= imax; ++i ) {
		ParameterCOP curparam( parameter_cop(i) );
		remark << "   " << curparam->short_parameter_description() << " (" << curparam->parameter_name() << "," << curparam->parameter_units() << "): ";
		ParameterType const paramtype( curparam->parameter_type() );
		switch( paramtype ) {
		case PT_generic_real:
		case PT_generic_nonnegative_valued_real:
		case PT_generic_positive_valued_real:
		case PT_angle:
			{
			RealValuedParameterCOP realparam( utility::pointer::static_pointer_cast< RealValuedParameter const >( curparam ) );
			remark << realparam->value();
			break;
		}
		case PT_boolean:
			{
			BooleanValuedParameterCOP boolparam( utility::pointer::static_pointer_cast< BooleanValuedParameter const >( curparam ) );
			remark << ( boolparam->value() ? "TRUE" : "FALSE" );
			break;
		}
		case PT_generic_integer:
		case PT_generic_whole_number:
		case PT_generic_natural_number:
			{
			SizeValuedParameterCOP sizeparam( utility::pointer::static_pointer_cast< SizeValuedParameter const >( curparam ) );
			remark << sizeparam->value();
			break;
		}
		case PT_generic_integer_vector:
		case PT_generic_whole_number_vector:
		case PT_generic_natural_number_vector:
			{
			SizeVectorValuedParameterCOP svparam( utility::pointer::static_pointer_cast< SizeVectorValuedParameter const >( curparam ) );
			remark << std::endl << "     " << svparam->value();
			break;
		}
		case PT_generic_real_vector:
		case PT_generic_nonnegative_valued_real_vector:
		case PT_generic_positive_valued_real_vector:
		case PT_angle_vector:
			{
			RealVectorValuedParameterCOP rvparam( utility::pointer::static_pointer_cast< RealVectorValuedParameter const >( curparam ) );
			remark << std::endl << "     " << rvparam->value();
			break;
		}
		case PT_invalid_type:
			remark << std::endl;
			utility_exit_with_message( "Error in BarrelParameters::get_pdb_remark(): Could not determine type of parameter " + curparam->parameter_name() + "!" );
		}
		remark << std::endl;
	}
}

} // namespace parameters
} // namespace beta_barrel
} // namespace protocols


#ifdef    SERIALIZATION

template< class Archive >
void
protocols::beta_barrel::parameters::BarrelParameters::save( Archive & arc ) const {
	arc( cereal::base_class< core::conformation::parametric::Parameters >( this ) );
}

template< class Archive >
void
protocols::beta_barrel::parameters::BarrelParameters::load( Archive & arc ) {
	arc( cereal::base_class< core::conformation::parametric::Parameters >( this ) );
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::beta_barrel::parameters::BarrelParameters );
CEREAL_REGISTER_TYPE( protocols::beta_barrel::parameters::BarrelParameters )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_beta_barrel_parameters_BarrelParameters )
#endif
