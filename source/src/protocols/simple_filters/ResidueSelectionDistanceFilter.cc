// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_filters/ResidueSelectionDistanceFilter.cc
/// @details loops over the residues in a residue selector and reports the average distance
/// @author TJ Brunette (tjbrunette@gmail.com)

#include <protocols/simple_filters/ResidueSelectionDistanceFilter.hh>
#include <protocols/simple_filters/ResidueSelectionDistanceFilterCreator.hh>

#include <utility/tag/Tag.hh>
#include <utility/vector1.hh>
#include <basic/datacache/DataMap.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>
#include <core/select/residue_selector/ResidueSelector.hh>
#include <core/select/residue_selector/util.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <iostream>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>


namespace protocols {
namespace simple_filters {

static basic::Tracer TR( "protocols.simple_filters.ResidueSelectionDistanceFilter" );

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP ResidueSelectionDistanceFilterCreator::create_filter() const { return protocols::filters::FilterOP( new ResidueSelectionDistanceFilter ); }

// XRW TEMP std::string
// XRW TEMP ResidueSelectionDistanceFilterCreator::keyname() const { return "ResidueSelectionDistance"; }

ResidueSelectionDistanceFilter::~ResidueSelectionDistanceFilter(){}


bool
ResidueSelectionDistanceFilter::apply( core::pose::Pose const & pose ) const {
	core::Real const distance( compute( pose ) );
	std::stringstream resids;
	core::select::residue_selector::ResidueSubset subset;
	subset =selector_->apply( pose );
	for ( core::Size resid=1; resid<=pose.total_residue(); ++resid ) {
		if ( subset[ resid ] ) {
			resids << resid << " ";
		}
	}
	TR << "ResidueSelectionDistanceFilter with distance threshold of "<<distance_threshold_<<" between residues "<< resids.str() << " is " << distance << std::endl;
	return( distance<=distance_threshold_ );
}

void
ResidueSelectionDistanceFilter::report( std::ostream & out, core::pose::Pose const & pose ) const {
	core::Real const distance( compute( pose ) );
	std::stringstream resids;
	core::select::residue_selector::ResidueSubset subset;
	subset =selector_->apply( pose );
	for ( core::Size resid=1; resid<=pose.total_residue(); ++resid ) {
		if ( subset[ resid ] ) {
			resids << resid << " ";
		}
	}
	out <<  "ResidueSelectionDistanceFilter with distance threshold of "<<distance_threshold_<<" between residues "<< resids.str() << " is " << distance << std::endl;
}

core::Real
ResidueSelectionDistanceFilter::report_sm( core::pose::Pose const & pose ) const {
	core::Real const distance( compute( pose ) );
	return( distance );
}
core::Real
ResidueSelectionDistanceFilter::compute( core::pose::Pose const & pose ) const {
	core::select::residue_selector::ResidueSubset subset;
	subset =selector_->apply( pose );
	utility::vector1<Size> valid_res;
	for ( core::Size resid=1; resid<=pose.total_residue(); ++resid ) {
		if ( subset[ resid ] ) {
			valid_res.push_back( resid );
		}
	}
	if ( valid_res.size()<=1 ) { //can not measure distance if < 1 residue
		return(0);
	}
	core::Size ct=0;
	core::Real total_distance=0;
	for ( core::Size ii=1; ii<=valid_res.size(); ++ii ) {
		for ( core::Size jj=ii+1; jj<=valid_res.size(); ++jj ) {
			core::conformation::Residue const res_res1( pose.conformation().residue( valid_res[ii] ) );
			core::conformation::Residue const res_res2( pose.conformation().residue( valid_res[jj] ) );
			if ( atom_to_measure_=="nbr" ) {
				core::Real const distance( res_res1.xyz( res_res1.nbr_atom() ).distance( res_res2.xyz( res_res2.nbr_atom() ) ) );
				total_distance += distance;
				ct+=1;
			} else {
				core::Real const distance( res_res1.xyz( atom_to_measure_ ).distance( res_res2.xyz( atom_to_measure_ ) ) );
				total_distance += distance;
				ct+=1;
			}
		}
	}
	return( total_distance/core::Real(ct) );
}

void
ResidueSelectionDistanceFilter::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & data,
	filters::Filters_map const &,
	moves::Movers_map const &,
	core::pose::Pose const &
) {
	atom_to_measure_ = "nbr";
	if ( tag->hasOption("atom") ) {
		atom_to_measure_ = tag->getOption< std::string >("atom");
	}
	distance_threshold_ = tag->getOption<core::Real>( "distance", 99.0 );
	selector_ = protocols::rosetta_scripts::parse_residue_selector( tag, data );
}

std::string ResidueSelectionDistanceFilter::name() const {
	return class_name();
}

std::string ResidueSelectionDistanceFilter::class_name() {
	return "ResidueSelectionDistance";
}

void ResidueSelectionDistanceFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;
	attlist + XMLSchemaAttribute("atom", xs_string, "Atom of interest")
		+ XMLSchemaAttribute::attribute_w_default("distance", xsct_real, "Distance threshold", "99.0");

	core::select::residue_selector::attributes_for_parse_residue_selector( attlist, "residue_selector", "parses residue selector" );

	protocols::filters::xsd_type_definition_w_attributes( xsd, class_name(), "XRW TO DO", attlist );
}

std::string ResidueSelectionDistanceFilterCreator::keyname() const {
	return ResidueSelectionDistanceFilter::class_name();
}

protocols::filters::FilterOP
ResidueSelectionDistanceFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new ResidueSelectionDistanceFilter );
}

void ResidueSelectionDistanceFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	ResidueSelectionDistanceFilter::provide_xml_schema( xsd );
}



}
}
