// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/simple_metrics/util.cc
/// @brief Util files for SimpleMetrics.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <core/simple_metrics/util.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <basic/Tracer.hh>

#include <core/pose/Pose.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/simple_metrics/SimpleMetric.hh>
#include <core/simple_metrics/SimpleMetricFactory.hh>
#include <core/select/util.hh>

#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
#include <utility/tag/util.hh>
#include <utility/string_util.hh>
#include <basic/datacache/DataMap.hh>

static basic::Tracer TR( "core.simple_metrics.util" );


namespace core {
namespace simple_metrics {

using namespace core::select;
using namespace core::select::residue_selector;

std::string
complex_type_name_for_simple_metric( std::string const & simple_metric_name){
	return "simple_metric_" + simple_metric_name + "_complex_type";
}

void
xsd_simple_metric_type_definition_w_attributes(
	utility::tag::XMLSchemaDefinition & xsd,
	std::string const & rs_type,
	std::string const & description,
	utility::tag::AttributeList const & attributes
)
{
	utility::tag::XMLSchemaComplexTypeGenerator ct_gen;
	ct_gen.complex_type_naming_func( & complex_type_name_for_simple_metric )
		.element_name( rs_type )
		.description( description )
		.add_attributes( attributes )
		.add_optional_name_attribute()
		.write_complex_type_to_schema( xsd );
}

std::map< core::Size, core::Size >
get_residue_mapping_from_selectors(
	select::residue_selector::ResidueSelectorCOP residue_selector,
	select::residue_selector::ResidueSelectorCOP residue_selector_ref,
	core::pose::Pose const & pose,
	core::pose::Pose const & ref_pose
){

	using namespace utility;

	vector1< core::Size > selector_res;
	vector1< core::Size > reference_res;

	if ( residue_selector ) {
		selector_res = get_residues_from_subset(residue_selector->apply( pose ));
	}
	if ( residue_selector_ref ) {
		reference_res = get_residues_from_subset(residue_selector_ref->apply( ref_pose ));
	}

	//Fail if residue selector selections do not match.
	if ( residue_selector && residue_selector_ref ) {
		if ( selector_res.size() != reference_res.size() ) {
			utility_exit_with_status("SimpleMetrics: Both set residue selectors must select the same number of residues in order to run RMSD calculation!");
		}
	}

	//Fail if only reference selector was set.
	if ( residue_selector_ref && (! residue_selector) ) {
		utility_exit_with_message("SimpleMetrics: Cannot only set the reference residue selector! If they are both the same, please set the main selector!");
	}

	//Setup the main residue mapping we will use depending on what is set.
	std::map< core::Size, core::Size > residue_map;
	if ( residue_selector && residue_selector_ref ) {
		for ( core::Size i = 1; i <= selector_res.size(); ++i ) {
			residue_map[ selector_res[ i] ] = reference_res[ i ];
		}
	} else if ( residue_selector ) {
		for ( core::Size res : selector_res ) {
			residue_map[ res ] = res;
		}
	} else {
		for ( core::Size res = 1; res <= pose.size(); ++res ) {
			residue_map[ res ] = res;
		}
	}
	return residue_map;
}

utility::vector1< SimpleMetricCOP >
get_metrics_from_datamap_and_subtags(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	std::string tag_name)
{
	utility::vector1< SimpleMetricCOP > metrics;

	if ( tag->hasOption(tag_name ) ) {
		std::string cst_gen_cslist = tag->getOption< std::string >( tag_name );
		utility::vector1< std::string > cst_gen_vector = utility::string_split( cst_gen_cslist, ',' );
		for ( std::string gen: cst_gen_vector ) {
			//Retrieve from the data map
			if ( datamap.has( "SimpleMetric", gen ) ) {
				SimpleMetricCOP metric( datamap.get_ptr< SimpleMetric >( "SimpleMetric", gen ) );
				metrics.push_back( metric );
				TR << "Added simple metric " << metric->name() << "." << std::endl;
			} else {
				throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "SimpleMetric " + gen + " not found in basic::datacache::DataMap.");
			}
		}
	}
	//Note that this mover is a little unusual in that it adds any constraint generators defined as subtags to the DataMap.
	for ( auto subtag=tag->getTags().begin(); subtag!=tag->getTags().end(); ++subtag ) {
		SimpleMetricOP metric = SimpleMetricFactory::get_instance()->new_simple_metric( (*subtag)->getName(), *subtag, datamap );
		metrics.push_back( metric );
		TR << "Added simple metric " << metric->name() << "." << std::endl;
	}

	return metrics;
}

SimpleMetricCOP
get_metric_from_datamap_and_subtags(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap,
	std::string tag_name)
{
	if ( tag->hasOption(tag_name ) ) {
		std::string gen = tag->getOption< std::string >( tag_name );
		if ( datamap.has( "SimpleMetric", gen ) ) {
			SimpleMetricCOP metric( datamap.get_ptr< SimpleMetric >( "SimpleMetric", gen ) );
			return metric;
		} else {
			throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "SimpleMetric " + gen + " not found in basic::datacache::DataMap.");
		}

	}

	//If multiple metrics are passed in, we fail.
	utility::vector1< SimpleMetricOP > metrics;
	for ( auto subtag=tag->getTags().begin(); subtag!=tag->getTags().end(); ++subtag ) {
		SimpleMetricOP metric = SimpleMetricFactory::get_instance()->new_simple_metric( (*subtag)->getName(), *subtag, datamap );
		metrics.push_back( metric );
	}

	if ( metrics.size() > 1 ) {
		throw CREATE_EXCEPTION(utility::excn::RosettaScriptsOptionError, "This class only accepts a single SimpleMetric as a subtag.");
	} else {
		return metrics[1];
	}
}


} //core
} //simple_metrics


