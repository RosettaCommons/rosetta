// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/multistage_rosetta_scripts/MRSJobSummary.cc
/// @brief  The implementation for class protocols::multistage_rosetta_scripts::MRSJobSummary's methods
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Unit headers
#include <protocols/multistage_rosetta_scripts/MRSJobSummary.hh>
#include <protocols/multistage_rosetta_scripts/cluster/ClusterMetric.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace multistage_rosetta_scripts {

MRSJobSummary::MRSJobSummary() :
	jd3::standard::EnergyJobSummary(),
	cluster_metric_( 0 )
{}

MRSJobSummary::MRSJobSummary( core::Real energy ) :
	jd3::standard::EnergyJobSummary( energy ),
	cluster_metric_( 0 )
{}

MRSJobSummary::MRSJobSummary(
	core::Real energy,
	cluster::ClusterMetricOP cluster_metric
) :
	jd3::standard::EnergyJobSummary( energy ),
	cluster_metric_( cluster_metric )
{}

MRSJobSummary::~MRSJobSummary(){}

} // namespace multistage_rosetta_scripts
} // namespace protocols

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::multistage_rosetta_scripts::MRSJobSummary::save( Archive & arc ) const {
	arc( cereal::base_class< protocols::jd3::standard::EnergyJobSummary >( this ) );
	arc( CEREAL_NVP( cluster_metric_ ) ); // cluster::ClusterMetricOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::multistage_rosetta_scripts::MRSJobSummary::load( Archive & arc ) {
	arc( cereal::base_class< protocols::jd3::standard::EnergyJobSummary >( this ) );
	arc( cluster_metric_ ); // cluster::ClusterMetricOP
}

SAVE_AND_LOAD_SERIALIZABLE( protocols::multistage_rosetta_scripts::MRSJobSummary );
CEREAL_REGISTER_TYPE( protocols::multistage_rosetta_scripts::MRSJobSummary )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_multistage_rosetta_scripts_MRSJobSummary )
#endif // SERIALIZATION
