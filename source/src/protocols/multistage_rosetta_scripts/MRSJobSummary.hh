// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/multistage_rosetta_scripts/MRSJobSummary.hh
/// @brief  The definition for class protocols::multistage_rosetta_scripts::MRSJobSummary
/// @details This summary build off of the EnergyJobSummary.
/// So far, it just includes an optional ClusterMetricOP.
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_protocols_multistage_rosetta_scripts_MRSJobSummary_HH
#define INCLUDED_protocols_multistage_rosetta_scripts_MRSJobSummary_HH

// Unit headers
#include <protocols/multistage_rosetta_scripts/MRSJobSummary.fwd.hh>
#include <protocols/multistage_rosetta_scripts/cluster/ClusterMetric.fwd.hh>

// Package headers
#include <protocols/jd3/standard/MoverAndPoseJob.hh>

#ifdef    SERIALIZATION
// Cereal headers
#include <cereal/types/polymorphic.fwd.hpp>
#endif // SERIALIZATION

namespace protocols {
namespace multistage_rosetta_scripts {

class MRSJobSummary : public jd3::standard::EnergyJobSummary
{
public:
	MRSJobSummary();

	MRSJobSummary( core::Real energy );

	MRSJobSummary( core::Real energy, cluster::ClusterMetricOP cluster_metric );

	~MRSJobSummary();

public://Setters and Getters
	void set_cluster_metric( cluster::ClusterMetricOP setting ) {
		cluster_metric_ = std::move( setting );
	}

	cluster::ClusterMetricOP cluster_metric() {
		return cluster_metric_;
	}

	cluster::ClusterMetricCOP cluster_metric() const {
		return cluster_metric_;
	}

private:
	cluster::ClusterMetricOP cluster_metric_;
#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif // SERIALIZATION

};

} // namespace multistage_rosetta_scripts
} // namespace protocols



#ifdef    SERIALIZATION
CEREAL_FORCE_DYNAMIC_INIT( protocols_multistage_rosetta_scripts_MRSJobSummary )
#endif // SERIALIZATION


#endif //INCLUDED
