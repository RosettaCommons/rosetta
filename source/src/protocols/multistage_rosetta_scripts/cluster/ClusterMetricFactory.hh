// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/multistage_rosetta_scripts/cluster/ClusterMetricFactory.hh
/// @brief Class for instantiating arbitrary Cluster Metrics from a string
/// @author Jack Maguire, jackmaguire1444@gmail.com

#ifndef INCLUDED_protocols_multistage_rosetta_scripts_cluster_ClusterMetricFactory_HH
#define INCLUDED_protocols_multistage_rosetta_scripts_cluster_ClusterMetricFactory_HH

// Package headers
#include <protocols/multistage_rosetta_scripts/cluster/ClusterMetric.hh>
#include <protocols/multistage_rosetta_scripts/cluster/ClusterMetricCreator.hh>

// Basic headers
#include <basic/datacache/DataMap.fwd.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/tag/XMLSchemaGeneration.fwd.hh>

// C++ headers
#include <map>
#include <string>

namespace protocols {
namespace multistage_rosetta_scripts {
namespace cluster {

class ClusterMetricFactory : public utility::SingletonBase< ClusterMetricFactory > {
public:
	friend class utility::SingletonBase< ClusterMetricFactory >;

	void factory_register( ClusterMetricCreatorOP creator );
	bool has_type( std::string const & ) const;

	void provide_xml_schema(
		std::string const & selector_name,
		utility::tag::XMLSchemaDefinition & xsd
	) const;

	ClusterMetricOP new_cluster_metric(
		std::string const & selector_name,
		core::pose::Pose const & pose,
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & datamap
	) const;

	/// @brief Should the Factory throw an exception or call utility::exit when it encounters the
	/// second of two ClusterMetricCreators with the same keyname?  It's default behavior is to
	/// call utility::exit, but this method allows you to set it so that it will throw an
	/// exception instead (which is unit testable).
	void set_throw_on_double_registration();

	void define_cluster_metric_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const;

	/// @brief Read access to the map of creator names to creators -- for unit testing purposes only
	std::map< std::string, ClusterMetricCreatorOP > const & creator_map() const;

	static std::string cluster_metric_xml_schema_group_name();

private:
	ClusterMetricFactory();

	// Unimplemented -- uncopyable
	ClusterMetricFactory( ClusterMetricFactory const & ) = delete;
	ClusterMetricFactory const & operator = ( ClusterMetricFactory const & ) = delete;

private:
	std::map< std::string, ClusterMetricCreatorOP > creator_map_;
	bool throw_on_double_registration_;
};


} //namespace cluster
} //namespace multistage_rosetta_scripts
} //namespace protocols


#endif
