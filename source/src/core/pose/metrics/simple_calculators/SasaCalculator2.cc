// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pose/metrics/SasaCalculator2.cc
/// @brief  SasaCalculator2 class
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com).  Based on SasaCalculatorLegacy

// Unit headers
#include <core/pose/metrics/simple_calculators/SasaCalculator2.hh>

#include <core/scoring/sasa/SasaCalc.hh>

// Utility headers
#include <basic/Tracer.hh>
#include <utility/exit.hh>
#include <utility/stream_util.hh>
#include <utility/string_util.hh>
#include <basic/MetricValue.hh>
#

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
using utility::vector1;
using core::Real;

SasaCalculator2::SasaCalculator2(){

	sasa_calc_ = core::scoring::sasa::SasaCalcOP( new core::scoring::sasa::SasaCalc() );
}

SasaCalculator2::SasaCalculator2(core::Real probe_r)

{
	sasa_calc_ = core::scoring::sasa::SasaCalcOP( new core::scoring::sasa::SasaCalc() );
	sasa_calc_->set_probe_radius(probe_r);

}
void SasaCalculator2::lookup( std::string const & key, basic::MetricValueBase * valptr ) const {

	if ( key == "total_sasa" ) {
		basic::check_cast( valptr, &total_sasa_, "total_sasa expects to return a Real" );
		(static_cast<basic::MetricValue<Real> *>(valptr))->set( total_sasa_);

	} else if ( key == "total_sasa_sc" ) {
		basic::check_cast( valptr, &total_sasa_sc_, "total_sasa_sc expects to return a Real");
		(static_cast<basic::MetricValue<Real> *>(valptr))->set(total_sasa_sc_);

	} else if ( key == "total_hsasa" ) {
		basic::check_cast( valptr, &total_hsasa_, "total_hsasa expects to return a Real");
		(static_cast<basic::MetricValue<Real> *>(valptr))->set(total_hsasa_);

	} else if ( key == "total_hsasa_sc" ) {
		basic::check_cast( valptr, &total_hsasa_sc_, "total_hsasa_sc expects to return a Real");
		(static_cast<basic::MetricValue<Real> *>(valptr))->set(total_hsasa_sc_);

	} else if ( key == "total_rel_hsasa" ) {
		basic::check_cast( valptr, &total_rel_hsasa_, "total_rel_hsasa expects to return a Real");
		(static_cast<basic::MetricValue<Real> *>(valptr))->set(total_rel_hsasa_);

	} else if ( key == "atom_sasa" ) {
		basic::check_cast( valptr, &atom_sasa_ , "atom_sasa expects to return a id::AtomID_Map< Real >" );
		(static_cast<basic::MetricValue<id::AtomID_Map< Real > > *>(valptr))->set( atom_sasa_ );

	} else if ( key == "residue_sasa" ) {
		basic::check_cast( valptr, &residue_sasa_, "residue_sasa expects to return a utility::vector1< Real >" );
		(static_cast<basic::MetricValue<utility::vector1< Real > > *>(valptr))->set( residue_sasa_);

	} else if ( key == "residue_sasa_sc" ) {
		basic::check_cast( valptr, &residue_sasa_sc_, "residue_sasa_sc expects to return a utility::vector1< Real >" );
		(static_cast<basic::MetricValue<utility::vector1< Real > > *>(valptr))->set( residue_sasa_sc_);

	} else if ( key == "residue_hsasa" ) {
		basic::check_cast( valptr, &residue_hsasa_, "residue_hsasa expects to return a utility::vector1< Real >");
		(static_cast<basic::MetricValue<utility::vector1< Real > > *>(valptr))->set(residue_hsasa_);

	} else if ( key == "residue_hsasa_sc" ) {
		basic::check_cast( valptr, &residue_hsasa_sc_, "residue_hsasa_sc expects to return a utility::vector1< Real >");
		(static_cast<basic::MetricValue<utility::vector1< Real > > *>(valptr))->set(residue_hsasa_sc_);

	} else if ( key == "residue_rel_hsasa" ) {
		basic::check_cast( valptr, &residue_rel_hsasa_, "residue_rel_hsasa expects to return a utility::vector1<Real>");
		(static_cast<basic::MetricValue<utility::vector1< Real > > *>(valptr))->set(residue_rel_hsasa_);
	} else {
		basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
		utility_exit();
	}

}


std::string SasaCalculator2::print( std::string const & key ) const {

	if ( key == "total_sasa" ) {
		return utility::to_string( total_sasa_ );
	} else if ( key == "total_sasa_sc" ) {
		return utility::to_string( total_sasa_sc_);
	} else if ( key == "total_hsasa" ) {
		return utility::to_string( total_hsasa_);
	} else if ( key == "total_hsasa_sc" ) {
		return utility::to_string( total_hsasa_sc_);
	} else if ( key == "total_rel_hsasa" ) {
		return utility::to_string(total_rel_hsasa_);
	} else if ( key == "atom_sasa" ) {
		basic::Error() << "id::AtomID_Map< Real > has no output operator, for metric " << key << std::endl;
		utility_exit();
	} else if ( key == "residue_sasa" ) {
		return utility::to_string( residue_sasa_ );
	} else if ( key == "residue_sasa_sc" ) {
		return utility::to_string( residue_sasa_sc_ );
	} else if ( key == "residue_hsasa" ) {
		return utility::to_string( residue_hsasa_ );
	} else if ( key == "residue_hsasa_sc" ) {
		return utility::to_string( residue_hsasa_sc_ );
	} else if ( key == "residue_rel_hsasa" ) {
		return utility::to_string(residue_rel_hsasa_);
	}

	basic::Error() << "This Calculator cannot compute metric " << key << std::endl;
	utility_exit();
	return "";

}


void SasaCalculator2::recompute( Pose const & this_pose ) {
	total_sasa_ = sasa_calc_->calculate(this_pose);
	sasa_calc_->fill_data(total_hsasa_, total_rel_hsasa_, atom_sasa_, residue_sasa_, residue_hsasa_, residue_rel_hsasa_);

	total_sasa_sc_ = sasa_calc_->get_total_sasa_sc();
	total_hsasa_sc_ = sasa_calc_->get_total_hsasa_sc();

	residue_sasa_sc_ = sasa_calc_->get_residue_sasa_sc();
	residue_hsasa_sc_ = sasa_calc_->get_residue_hsasa_sc();

}


} // simple_calculators
} // metrics
} // pose
} // core


#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::pose::metrics::simple_calculators::SasaCalculator2::save( Archive & arc ) const {
	arc( cereal::base_class< core::pose::metrics::StructureDependentCalculator >( this ) );
	arc( CEREAL_NVP( sasa_calc_ ) ); // core::scoring::sasa::SasaCalcOP
	arc( CEREAL_NVP( total_sasa_ ) ); // core::Real
	arc( CEREAL_NVP( total_hsasa_ ) ); // core::Real
	arc( CEREAL_NVP( total_rel_hsasa_ ) ); // core::Real
	arc( CEREAL_NVP( atom_sasa_ ) ); // core::id::AtomID_Map<core::Real>
	arc( CEREAL_NVP( residue_sasa_ ) ); // vector1<core::Real>
	arc( CEREAL_NVP( residue_hsasa_ ) ); // vector1<core::Real>
	arc( CEREAL_NVP( residue_rel_hsasa_ ) ); // vector1<core::Real>
	arc( CEREAL_NVP( total_sasa_sc_ ) ); // core::Real
	arc( CEREAL_NVP( total_hsasa_sc_ ) ); // core::Real
	arc( CEREAL_NVP( residue_sasa_sc_ ) ); // vector1<core::Real>
	arc( CEREAL_NVP( residue_hsasa_sc_ ) ); // vector1<core::Real>
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::pose::metrics::simple_calculators::SasaCalculator2::load( Archive & arc ) {
	arc( cereal::base_class< core::pose::metrics::StructureDependentCalculator >( this ) );
	arc( sasa_calc_ ); // core::scoring::sasa::SasaCalcOP
	arc( total_sasa_ ); // core::Real
	arc( total_hsasa_ ); // core::Real
	arc( total_rel_hsasa_ ); // core::Real
	arc( atom_sasa_ ); // core::id::AtomID_Map<core::Real>
	arc( residue_sasa_ ); // vector1<core::Real>
	arc( residue_hsasa_ ); // vector1<core::Real>
	arc( residue_rel_hsasa_ ); // vector1<core::Real>
	arc( total_sasa_sc_ ); // core::Real
	arc( total_hsasa_sc_ ); // core::Real
	arc( residue_sasa_sc_ ); // vector1<core::Real>
	arc( residue_hsasa_sc_ ); // vector1<core::Real>
}

SAVE_AND_LOAD_SERIALIZABLE( core::pose::metrics::simple_calculators::SasaCalculator2 );
CEREAL_REGISTER_TYPE( core::pose::metrics::simple_calculators::SasaCalculator2 )

CEREAL_REGISTER_DYNAMIC_INIT( core_pose_metrics_simple_calculators_SasaCalculator2 )
#endif // SERIALIZATION
