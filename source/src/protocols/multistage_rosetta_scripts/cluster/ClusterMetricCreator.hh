// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/multistage_rosetta_scripts/cluster/ClusterMetricCreator.hh
/// @brief
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_protocols_multistage_rosetta_scripts_cluster_ClusterMetricCreator_HH
#define INCLUDED_protocols_multistage_rosetta_scripts_cluster_ClusterMetricCreator_HH

#include <utility/pointer/owning_ptr.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>
#include <protocols/multistage_rosetta_scripts/cluster/ClusterMetric.fwd.hh>

namespace protocols {
namespace multistage_rosetta_scripts {
namespace cluster {

class ClusterMetricCreator : public utility::pointer::ReferenceCount
{
public:
	ClusterMetricCreator(){}
	virtual ~ClusterMetricCreator(){}

	/// @brief Return a new metric.
	virtual ClusterMetricOP create_metric() const = 0;

	/// @brief Return the tag name associated with this factory.
	virtual std::string keyname() const = 0;

	/// @brief Describe the schema for the Cluster Metric that this Creator is responsible for
	virtual void provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const = 0;

};

using ClusterMetricCreatorOP = utility::pointer::shared_ptr< ClusterMetricCreator >;
using ClusterMetricCreatorCOP = utility::pointer::shared_ptr< ClusterMetricCreator const >;

} // cluster
} // multistage_rosetta_scripts
} // protocols

#endif
