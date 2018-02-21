// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/multistage_rosetta_scripts/cluster/ClusterMetricFactory.cc
/// @brief  Class for instantiating arbitrary Cluster Metrics from a string
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Unit headers
#include <protocols/multistage_rosetta_scripts/cluster/ClusterMetricFactory.hh>

// Package headers
#include <protocols/multistage_rosetta_scripts/cluster/ClusterMetric.hh>
#include <protocols/multistage_rosetta_scripts/cluster/ClusterMetricCreator.hh>
#include <protocols/multistage_rosetta_scripts/cluster/util.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/xml_schema_group_initialization.hh>

// Boost headers
#include <boost/bind.hpp>

namespace protocols {
namespace multistage_rosetta_scripts {
namespace cluster {

void
ClusterMetricFactory::factory_register( ClusterMetricCreatorOP creator )
{
	if ( creator_map_.find( creator->keyname() ) != creator_map_.end() ) {
		std::string const err_msg =
			"Factory Name Conflict: Two or more ClusterMetricCreators registered with the name " + creator->keyname();
		if ( throw_on_double_registration_ ) {
			throw CREATE_EXCEPTION(utility::excn::Exception,  err_msg );
		} else {
			utility_exit_with_message( err_msg );
		}
	}
	creator_map_[ creator->keyname() ] = creator;
}

bool ClusterMetricFactory::has_type( std::string const & selector_type ) const
{
	return creator_map_.find( selector_type ) != creator_map_.end();
}

void
ClusterMetricFactory::provide_xml_schema(
	std::string const &selector_name,
	utility::tag::XMLSchemaDefinition & xsd
) const {
	if ( ! has_type( selector_name ) ) {
		std::string err_msg =  "No ClusterMetricCreator with the name '" + selector_name + "' has been registered with the ClusterMetricFactory";
		throw CREATE_EXCEPTION(utility::excn::Exception,  err_msg );
	}
	auto iter = creator_map_.find( selector_name );
	iter->second->provide_xml_schema( xsd );
}

ClusterMetricOP ClusterMetricFactory::new_cluster_metric(
	std::string const & selector_name,
	core::pose::Pose const & pose,
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap
) const {
	if ( ! has_type( selector_name ) ) {
		std::string err_msg =  "No ClusterMetricCreator with the name '" + selector_name + "' has been registered with the ClusterMetricFactory";
		throw CREATE_EXCEPTION(utility::excn::Exception,  err_msg );
	}
	auto iter = creator_map_.find( selector_name );
	ClusterMetricOP new_selector = iter->second->create_metric();
	new_selector->parse_my_tag( pose, tag, datamap );
	return new_selector;
}

/// @details By convention, the named assigned to each of the complexTypes for ClusterMetrics should be
/// what is returned by the function "complex_type_name_for_cluster_metric" (declared in
/// protocols/multistage_rosetta_scripts/cluster/util.hh) when given the argument returned by that ClusterMetric's
/// ClusterMetricCreator's keyname() function. So long as the writing of XML schema for your residue
/// selector is accomplished by the calling the functions in protocols/multistage_rosetta_scripts/cluster/util.hh, then
/// this should happen automatically.
void ClusterMetricFactory::define_cluster_metric_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	try {
		utility::tag::define_xml_schema_group (
			creator_map_,
			cluster_metric_xml_schema_group_name(),
			& complex_type_name_for_cluster_metric,
			xsd );
	} catch ( utility::excn::Exception const & e ) {
		throw CREATE_EXCEPTION(utility::excn::Exception,
			"Could not generate an XML Schema for ClusterMetrics from ClusterMetricFactory; offending class"
			" must call protocols::multistage_rosetta_scripts::cluster::complex_type_name_for_cluster_metric when defining"
			" its XML Schema\n" + e.msg() );
	}

}

std::map< std::string, ClusterMetricCreatorOP > const &
ClusterMetricFactory::creator_map() const {
	return creator_map_;
}

std::string ClusterMetricFactory::cluster_metric_xml_schema_group_name()
{
	return "cluster_metric";
}


void ClusterMetricFactory::set_throw_on_double_registration()
{
	throw_on_double_registration_ = true;
}


ClusterMetricFactory::ClusterMetricFactory() :
	throw_on_double_registration_( false )
{}


} //namespace cluster
} //namespace multistage_rosetta_scripts
} //namespace protocols
