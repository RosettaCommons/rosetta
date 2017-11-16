// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_filters/NeighborTypeFilter.cc
/// @brief
/// @author Sarel Fleishman (sarelf@u.washington.edu), Jacob Corn (jecorn@u.washington.edu)

#include <protocols/simple_filters/NeighborTypeFilter.hh>
#include <protocols/simple_filters/NeighborTypeFilterCreator.hh>

#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/Conformation.hh>
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <protocols/rosetta_scripts/util.hh>
#include <core/pose/selection.hh>
// XSD XRW Includes
#include <utility/tag/XMLSchemaGeneration.hh>
#include <protocols/filters/filter_schemas.hh>

namespace protocols {
namespace simple_filters {

static basic::Tracer neighbor_type_filter_tracer( "protocols.simple_filters.NeighborTypeFilter" );

// XRW TEMP protocols::filters::FilterOP
// XRW TEMP NeighborTypeFilterCreator::create_filter() const { return protocols::filters::FilterOP( new NeighborTypeFilter ); }

// XRW TEMP std::string
// XRW TEMP NeighborTypeFilterCreator::keyname() const { return "NeighborType"; }

NeighborTypeFilter::~NeighborTypeFilter() = default;

void
NeighborTypeFilter::parse_my_tag( utility::tag::TagCOP tag, basic::datacache::DataMap &, filters::Filters_map const &, moves::Movers_map const &, core::pose::Pose const & )
{
	residue_types_.assign( core::chemical::num_canonical_aas, false );
	utility::vector0< utility::tag::TagCOP > const neighbor_type_tags( tag->getTags() );
	for ( auto nt_tag_ptr : neighbor_type_tags ) {
		if ( nt_tag_ptr->getName() == "Neighbor" ) {
			std::string const type( nt_tag_ptr->getOption<std::string>( "type" ) );
			residue_types_[ core::chemical::aa_from_name( type ) ] = true;
		}
	}
	target_residue_ = core::pose::get_resnum_string( tag );
	distance_threshold_ = tag->getOption<core::Real>( "distance", 8.0 );

	neighbor_type_filter_tracer<<"NeighborTypeFilter with distance threshold of "<<distance_threshold_<<" around residue "<<target_residue_<<std::endl;
}

bool
NeighborTypeFilter::apply( core::pose::Pose const & pose ) const
{
	std::vector< core::Size > neighbors = compute( pose );
	if ( neighbors.size() == 0 ) return false;
	neighbor_type_filter_tracer<<"neighbours of residue "<<target_residue_<<": ";
	for ( std::vector< core::Size >::const_iterator n_it=neighbors.begin(); n_it!=neighbors.end(); ++n_it ) {
		neighbor_type_filter_tracer<<pose.residue( *n_it ).name3()<<*n_it<<" ";
	}
	neighbor_type_filter_tracer<<std::endl;
	return true;
}

void
NeighborTypeFilter::report( std::ostream & out, core::pose::Pose const & pose ) const
{
	std::vector< core::Size > neighbors = compute( pose );
	out<<"neighbours of residue "<<target_residue_<<": ";
	for ( std::vector< core::Size >::const_iterator n_it=neighbors.begin(); n_it!=neighbors.end(); ++n_it ) {
		out<<pose.residue( *n_it ).name3()<<*n_it<<" ";
	}
	out<<'\n';
}

core::Real
NeighborTypeFilter::report_sm( core::pose::Pose const & pose ) const
{
	std::vector< core::Size > neighbors = compute( pose );
	return( neighbors.size() != 0 );
}

std::string NeighborTypeFilter::name() const {
	return class_name();
}

std::string NeighborTypeFilter::class_name() {
	return "NeighborType";
}

std::string NeighborType_subelement_ct_name( std::string const & name ) {
	return "NeighborType_subelement_" + name + "Type";
}

void NeighborTypeFilter::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd )
{
	using namespace utility::tag;
	AttributeList attlist;

	AttributeList NeighborType_subtag_attributes;
	NeighborType_subtag_attributes
		+ XMLSchemaAttribute( "type", xs_string, "XRW TO DO");

	XMLSchemaSimpleSubelementList subelements;
	subelements.complex_type_naming_func( & NeighborType_subelement_ct_name );
	subelements.add_simple_subelement( "Neighbors", NeighborType_subtag_attributes, "XRW TO DO");

	core::pose::attributes_for_get_resnum_string( attlist, "" );

	attlist + XMLSchemaAttribute::attribute_w_default( "distance", xsct_real, "XRW TO DO", "8.0");
	attlist + XMLSchemaAttribute::attribute_w_default( "confidence", xsct_real, "Probability that this pose will be filtered out if the filter fails", "1.0");

	XMLSchemaComplexTypeGenerator complex_type_generator;
	complex_type_generator
		.element_name( class_name() )
		.description( "Filter for poses that place a neighbour of the types specified around a target residue in the partner protein." )
		.complex_type_naming_func( & protocols::filters::complex_type_name_for_filter )
		.add_attributes( attlist )
		.set_subelements_repeatable( subelements )
		.write_complex_type_to_schema( xsd );
}

std::string NeighborTypeFilterCreator::keyname() const {
	return NeighborTypeFilter::class_name();
}

protocols::filters::FilterOP
NeighborTypeFilterCreator::create_filter() const {
	return protocols::filters::FilterOP( new NeighborTypeFilter );
}

void NeighborTypeFilterCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const
{
	NeighborTypeFilter::provide_xml_schema( xsd );
}


std::vector< core::Size >
NeighborTypeFilter::compute( core::pose::Pose const & pose ) const
{
	std::vector< core::Size > neighbors;

	core::Size const chain2begin( pose.conformation().chain_begin( 2 ) );
	core::Size const residue_num( pose.size() );

	core::Size const target_residue( core::pose::parse_resnum( target_residue_, pose ) );
	core::Size const start( target_residue < chain2begin ? chain2begin : 1 );
	core::Size const end( target_residue < chain2begin ? residue_num : chain2begin-1 );
	core::conformation::Residue const res_target( pose.residue( target_residue ) );

	runtime_assert( target_residue <= residue_num );
	for ( core::Size res=start; res<=end; ++res ) {
		core::conformation::Residue const resi( pose.residue( res ) );
		if ( !residue_types_[ resi.aa() ] ) continue;

		core::Real const distance( res_target.xyz( res_target.nbr_atom() ).distance( resi.xyz( resi.nbr_atom() )) );
		if ( distance <= distance_threshold_ ) {
			neighbors.push_back( res );
		}
	}
	return neighbors;
}

}
}
