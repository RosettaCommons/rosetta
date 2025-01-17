// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/drug_design/RDKitMetricFilter.cc
/// @brief
/// @author Rocco Moretti (rmorettiase@gmail.com)

//Unit Headers
#include <protocols/drug_design/RDKitMetricFilter.hh>
#include <protocols/drug_design/RDKitMetricFilterCreator.hh>

//Project Headers
#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>
#include <utility/tag/Tag.hh>
#include <basic/datacache/DataMap.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/rdkit/RestypeToRDMol.hh>
#include <core/chemical/rdkit/util.hh>

#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/util.hh>
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>

#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace drug_design {

static basic::Tracer TR( "protocols.drug_design.RDKitMetricFilter" );

protocols::filters::FilterOP
RDKitMetricFilterCreator::create_filter() const { return protocols::filters::FilterOP( new RDKitMetricFilter ); }

std::string
RDKitMetricFilterCreator::keyname() const { return RDKitMetricFilter::class_name(); }

void RDKitMetricFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	RDKitMetricFilter::provide_xml_schema( xsd );
}

std::string RDKitMetricFilter::class_name() {
	return "RDKitMetric";
}

void RDKitMetricFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	attlist
		+ XMLSchemaAttribute::required_attribute( "residue", xsct_refpose_enabled_residue_number, "Residue to compute the metrics on" )
		+ XMLSchemaAttribute::required_attribute( "metric", xs_string, "Metric to calculate" )
		+ XMLSchemaAttribute( "lower_threshold", xsct_real,
		"Filter is false if the metric is less than this value" )
		+ XMLSchemaAttribute( "upper_threshold", xsct_real,
		"Filter is false if the metric is greater than this value" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(),
		"Test the number of heteroatom-heteroatom (non C/H) bonds in a residue.",
		attlist );
}

void
RDKitMetricFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & )
{
	residue_ = tag->getOption< std::string >( "residue" );
	metric( tag->getOption<std::string>( "metric" ) ); // For value checking.
	lower_threshold_ = tag->getOption<core::Real>( "lower_threshold", lower_threshold_ );
	upper_threshold_ = tag->getOption<core::Real>( "upper_threshold", upper_threshold_ );
	TR << "RDKitMetric filter for residue " << residue_ << " metric " << metric_ << " between " << lower_threshold_ << " and " << upper_threshold_ << std::endl;
}

void
RDKitMetricFilter::metric( std::string const & setting ) {
	if ( setting.size() != 0 ) {
		std::map< std::string, std::string > const metrics( core::chemical::rdkit::get_metric_names() );
		if ( metrics.count( setting ) == 0 ) {
			TR.Error << "Metric " << setting << " not recognized by Rosetta as a valid RDKit metric. Valid metrics are:\n\n";
			for ( std::map< std::string, std::string >::const_iterator itr( metrics.begin() ), itr_end( metrics.end());
					itr != itr_end; ++itr ) {
				TR.Error << itr->first << "\t-\t" << itr->second << "\n";
			}
			TR << std::endl;
			utility_exit_with_message("Rosetta doesn't understand '"+setting+"' as an RDKit metric.");
		}
	}
	metric_ = setting;
}

bool
RDKitMetricFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const value( compute( pose ) );
	if ( value < lower_threshold_ ) {
		TR << "Failing RDKitMetric filter on residue " << residue_ << " with " << metric_ << " value " << value << " which is below the lower threshold of " << lower_threshold_ << std::endl;
		return false;
	} else if ( value > upper_threshold_ ) {
		TR << "Failing RDKitMetric filter on residue " << residue_ << " with " << metric_ << " value " << value << " which is above the upper threshold of " << upper_threshold_ << std::endl;
		return false;
	} else {
		return true;
	}
}

void
RDKitMetricFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const value( compute( pose ) );
	out << "RDKit metric '" << metric_ << "' for residue " << residue_ << " is " << value << std::endl;
}

core::Real
RDKitMetricFilter::report_sm( core::pose::Pose const & pose ) const {
	return( compute( pose ) );
}

core::Real
RDKitMetricFilter::compute( core::pose::Pose const & pose ) const {
	core::Size resnum( core::pose::parse_resnum( residue_, pose ) );
	if ( resnum == 0 || resnum > pose.total_residue() ) {
		TR.Error << "Attempted to access residue " << residue_ << " in pose with " <<  pose.total_residue() << " residues. Failing filter. " << std::endl;
		utility_exit_with_message("Cannot apply RDKitMetric filter on non-existant residue!");
	}
	if ( metric_.size() == 0 ) {
		TR.Error << "Metric must be set for RDKitMetric" << std::endl;
		utility_exit_with_message("Must set metric in RDKitMetricFilter.");
	}

	core::chemical::MutableResidueTypeOP restype( utility::pointer::make_shared< core::chemical::MutableResidueType >(pose.residue(resnum).type()) );
	TR << "Calculating RDKit metric '" << metric_ << "' value for residue " << residue_ << ", of type " << restype->name() << std::endl;

	core::Real retval(0);

	// Use a neutral, hydrogen free molecule for calculating the metric values - should be the default
	core::chemical::rdkit::RestypeToRDMol converter(*restype);
	::RDKit::RWMolOP rdmol(converter.Mol() );
	retval = core::chemical::rdkit::rdkit_metric( *rdmol, metric_ );

	return retval;
}

}
}
