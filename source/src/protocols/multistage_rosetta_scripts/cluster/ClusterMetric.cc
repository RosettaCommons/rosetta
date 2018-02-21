
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW CoMotion, email: license@uw.edu.

/// @file   protocols/multistage_rosetta_scripts/cluster/ClusterMetric.cc
/// @brief  Auto-generated serialization template functions
/// @author Jack Maguire, jackmaguire1444@gmail.com

#include <protocols/multistage_rosetta_scripts/cluster/ClusterMetric.hh>

#ifdef SERIALIZATION
#include <utility/serialization/serialization.hh>
#include <cereal/types/polymorphic.hpp>
#endif

protocols::multistage_rosetta_scripts::cluster::ClusterMetric::ClusterMetric(){}
protocols::multistage_rosetta_scripts::cluster::ClusterMetric::~ClusterMetric(){}

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
protocols::multistage_rosetta_scripts::cluster::ClusterMetric::save( Archive & ) const {}

/// @brief Automatically generated deserialization method
template< class Archive >
void
protocols::multistage_rosetta_scripts::cluster::ClusterMetric::load( Archive & ) {}

SAVE_AND_LOAD_SERIALIZABLE( protocols::multistage_rosetta_scripts::cluster::ClusterMetric );
CEREAL_REGISTER_TYPE( protocols::multistage_rosetta_scripts::cluster::ClusterMetric )

CEREAL_REGISTER_DYNAMIC_INIT( protocols_multistage_rosetta_scripts_cluster_ClusterMetric )
#endif // SERIALIZATION
