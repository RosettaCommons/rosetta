// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/multistage_rosetta_scripts/cluster/metrics/JumpMetricCreator.hh
/// @brief
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_protocols_multistage_rosetta_scripts_cluster_metric_JumpMetricCreator_hh
#define INCLUDED_protocols_multistage_rosetta_scripts_cluster_metric_JumpMetricCreator_hh

// unit headers
//#include <protocols/multistage_rosetta_scripts/cluster/metrics/JumpMetric.hh>
#include <protocols/multistage_rosetta_scripts/cluster/ClusterMetricCreator.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

namespace protocols {
namespace multistage_rosetta_scripts {
namespace cluster {
namespace metrics {

class JumpMetricCreator : public ClusterMetricCreator {
public:
	ClusterMetricOP create_metric() const override;

	std::string keyname() const override {
		return "Jump";
	}

	void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const override;
};

} // namespace metrics
} // namespace cluster
} // namespace multistage_rosetta_scripts
} // namespace protocols

#endif
