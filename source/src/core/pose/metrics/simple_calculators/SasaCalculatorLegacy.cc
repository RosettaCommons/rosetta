// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/metrics/SasaCalculatorLegacy.cc
/// @brief  SasaCalculatorLegacy class
/// @author John Karanicolas

// Unit headers
#include <core/pose/metrics/simple_calculators/SasaCalculatorLegacy.hh>

//#include <core/conformation/Residue.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
#include <basic/MetricValue.hh>
#include <core/scoring/sasa.hh>

#include <utility/assert.hh>

#include <utility/vector1.hh>


using namespace core;
using namespace core::pose;
using namespace core::pose::metrics;

#ifdef    SERIALIZATION
// Project serialization headers
#include <core/id/AtomID_Map.srlz.hh>

// Utility serialization headers
#include <utility/vector1.srlz.hh>
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace core {
namespace pose {
namespace metrics {
namespace simple_calculators {

void SasaCalculatorLegacy::lookup( std::string const & key, basic::MetricValueBase * valptr ) const {

	if ( key == "total_sasa" ) {
		basic::check_cast( valptr, &total_sasa_, "total_sasa expects to return a Real" );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( total_sasa_ );

	} else if ( key == "atom_sasa" ) {
		basic::check_cast( valptr, &atom_sasa_, "atom_sasa expects to return a id::AtomID_Map< Real >" );
		(static_cast<basic::MetricValue<id::AtomID_Map< Real > > *>(valptr))->set( atom_sasa_ );

	} else if ( key == "residue_sasa" ) {
		basic::check_cast( valptr, &residue_sasa_, "residue_sasa expects to return a utility::vector1< Real >" );
		(static_cast<basic::MetricValue<utility::vector1< Real > > *>(valptr))->set( residue_sasa_ );
		//} else if ( key == "residue_hsasa") {
		// basic::check_cast( valptr, &residue_sasa_h_, "residue_sasa_h expects to return a utility::vector1< Real >");
		// (static_cast<basic::MetricValue<utility::vector1< Real > > *>(valptr))->set( residue_sasa_h_);
	} else {
		basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
		utility_exit();
	}

}


std::string SasaCalculatorLegacy::print( std::string const & key ) const {

	if ( key == "total_sasa" ) {
		return utility::to_string( total_sasa_ );
	} else if ( key == "atom_sasa" ) {
		basic::Error() << "id::AtomID_Map< Real > has no output operator, for metric " << key << std::endl;
		utility_exit();
	} else if ( key == "residue_sasa" ) {
		return utility::to_string( residue_sasa_ );
		//} else if ( key == "residue_hsasa") {
		// return utility::to_string( residue_hsasa_);
	}

	basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
	utility_exit();
	return "";

}


void SasaCalculatorLegacy::recompute( Pose const & this_pose ) {
	total_sasa_ = core::scoring::calc_per_atom_sasa( this_pose, atom_sasa_, residue_sasa_, probe_radius_);
}


} // simple_calculators
} // metrics
} // pose
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::metrics::simple_calculators::SasaCalculatorLegacy::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::metrics::StructureDependentCalculator >( this ) );
	arc( CEREAL_NVP( total_sasa_ ) ); // core::Real
	arc( CEREAL_NVP( atom_sasa_ ) ); // core::id::AtomID_Map<core::Real>
	arc( CEREAL_NVP( residue_sasa_ ) ); // utility::vector1<core::Real>
	arc( CEREAL_NVP( probe_radius_ ) ); // core::Real
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::metrics::simple_calculators::SasaCalculatorLegacy::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::metrics::StructureDependentCalculator >( this ) );
	arc( total_sasa_ ); // core::Real
	arc( atom_sasa_ ); // core::id::AtomID_Map<core::Real>
	arc( residue_sasa_ ); // utility::vector1<core::Real>
	arc( probe_radius_ ); // core::Real
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::metrics::simple_calculators::SasaCalculatorLegacy );
CEREAL_REGISTER_TYPE( core::pose::metrics::simple_calculators::SasaCalculatorLegacy )

CEREAL_REGISTER_DYNAMIC_INIT( core_pose_metrics_simple_calculators_SasaCalculatorLegacy )
#endif // SERIALIZATION
