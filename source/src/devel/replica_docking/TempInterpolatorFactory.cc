// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Zhe Zhang

#include <devel/replica_docking/TempInterpolatorFactory.hh>
#include <devel/replica_docking/TempInterpolator.hh>

#include <protocols/canonical_sampling/TemperatureController.hh>

// Basic headers
#include <basic/Tracer.hh>

// Utility headers
#include <utility/exit.hh> // runtime_assert, throw utility::excn::EXCN_RosettaScriptsOption
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace devel {
namespace replica_docking {

static basic::Tracer TR( "devel.replica_docking.TempInterpolatorFactory" );

TempInterpolatorFactory::TempInterpolatorFactory()= default;

TempInterpolatorFactory::~TempInterpolatorFactory() = default;

/// @brief return new TempInterpolator by key lookup in mover_prototype_map_ (new TempInterpolator parses Tag if provided)
TempInterpolatorBaseOP
TempInterpolatorFactory::new_tempInterpolator( utility::tag::TagCOP tag, core::Size n_levels )
{
	if ( !tag->hasOption( "curve" ) ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "Error: interpolation curve type required !" );
	}
	std::string curve = tag->getOption< std::string >( "curve" );
	if ( curve == "const" ) {
		return TempInterpolatorBaseOP( new TempFixValue( tag->getOption< core::Real >( "value" ) ) );
	} else {
		if ( tag->hasOption( "start" ) && tag->hasOption( "end" ) ) {
			auto start = tag->getOption< core::Real >( "start" );
			auto end = tag->getOption< core::Real >( "end" );
			return TempInterpolatorBaseOP( new TempInterpolator( n_levels, start, end, curve ) );
		} else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError,  "Error: start and end value must be given for linear and expotential interpolation !" );
		}
	}
}

void
TempInterpolatorFactory::attributes_for_interpolators(
	utility::tag::XMLSchemaDefinition & xsd,
	utility::tag::AttributeList & interp_attributes
)
{
	using namespace utility::tag;
	using Attr = utility::tag::XMLSchemaAttribute;

	XMLSchemaRestriction temp_interpolator_type;
	temp_interpolator_type.name( "temp_interpolator_type" );
	temp_interpolator_type.base_type( xs_string );
	temp_interpolator_type.add_restriction( xsr_enumeration, "const" );
	temp_interpolator_type.add_restriction( xsr_enumeration, "exponential" );
	temp_interpolator_type.add_restriction( xsr_enumeration, "linear" );
	xsd.add_top_level_element( temp_interpolator_type );

	std::string const nonconst_warning("Only used if curve = 'linear' or 'exponential'.");

	interp_attributes
		+ Attr::required_attribute( "curve", "temp_interpolator_type", "Does this Interpolator use a 'const', 'linear', or 'exponential' curve?" )
		+ Attr( "start", xsct_real, "start temperature to interpolate from. " + nonconst_warning ) // not required because it is not used if curve = const
		+ Attr( "end", xsct_real, "end temperature to interpolate to; " + nonconst_warning ) // not requried because it is not used if curve = const
		+ Attr( "value", xsct_real, "temperature value, only used if curve='const' (in other words, interpolate from 'value' to 'value'" ); // not requried because it is not used if curve != const
}


} //namespace devel
} //namespace replica_docking
