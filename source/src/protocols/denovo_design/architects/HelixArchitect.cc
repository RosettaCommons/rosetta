// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/architects/HelixArchitect.cc
/// @brief Architect for helices
/// @author Tom Linsky (tlinsky@gmail.com)

// Unit headers
#include <protocols/denovo_design/architects/HelixArchitect.hh>
#include <protocols/denovo_design/architects/HelixArchitectCreator.hh>
#include <protocols/denovo_design/architects/DeNovoArchitectFactory.hh>
// Protocol headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/util.hh>

// Basic/Utility headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

static basic::Tracer TR( "protocols.denovo_design.architects.HelixArchitect" );

namespace protocols {
namespace denovo_design {
namespace architects {

HelixArchitect::HelixArchitect( std::string const & id_value ):
	DeNovoArchitect( id_value ),
	motifs_(),
	lengths_()
{
}

HelixArchitect::~HelixArchitect() = default;

DeNovoArchitectOP
HelixArchitect::clone() const
{
	return DeNovoArchitectOP( new HelixArchitect( *this ) );
}

std::string
HelixArchitect::type() const
{
	return HelixArchitect::class_name();
}

void
HelixArchitectCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	HelixArchitect::provide_xml_schema( xsd );
}

void
HelixArchitect::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;
	XMLSchemaRestriction length_string;
	length_string.name( "length_string" );
	length_string.base_type( xs_string );
	length_string.add_restriction( xsr_pattern, "[0-9]+(:[0-9]+)?(,[0-9]+(:[0-9]+)?)*" );
	xsd.add_top_level_element( length_string );

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "length", "length_string", "Comma-separated list of single integers and colon-separated ranges to specify all possible helix lengths" );

	DeNovoArchitect::add_common_denovo_architect_attributes( attlist );
	DeNovoArchitectFactory::xsd_architect_type_definition_w_attributes( xsd, class_name(), "Design helical segments", attlist );

}



void
HelixArchitect::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & )
{
	set_lengths( tag->getOption< std::string >( "length" ) );
}

DeNovoArchitect::StructureDataOP
HelixArchitect::design( core::pose::Pose const &, core::Real & random ) const
{
	if ( motifs_.empty() ) {
		utility_exit_with_message( type() + "Architect requires one or more motifs or lengths to be specified" );
	}
	core::Size const motif_idx = extract_int( random, 1, motifs_.size() );

	components::Segment const & chosen_motif = *motifs_[ motif_idx ];
	TR << "Designing segment: " << chosen_motif << std::endl;
	StructureDataOP sd( new components::StructureData( id() ) );
	sd->add_segment( chosen_motif );

	return sd;
}

void
HelixArchitect::set_lengths( std::string const & lengths_str )
{
	set_lengths( parse_length_str< core::Size >( lengths_str ) );
}

void
HelixArchitect::set_lengths( Lengths const & lengths )
{
	lengths_ = lengths;
	// add motifs
	motifs_.clear();
	for ( unsigned long length : lengths ) {
		std::stringstream ss;
		ss << 'L' << std::string( length, 'H' ) << 'L';
		std::stringstream abego;
		abego << 'X' << std::string( length, 'A' ) << 'X';

		components::SegmentOP motif( new components::Segment( id() ) );
		motif->extend( ss.str(), abego.str() );
		motifs_.push_back( motif );
	}
}

components::SegmentCOPs::const_iterator
HelixArchitect::motifs_begin() const
{
	return motifs_.begin();
}

components::SegmentCOPs::const_iterator
HelixArchitect::motifs_end() const
{
	return motifs_.end();
}

///////////////////////////////////////////////////////////////////////////////
/// Creator
///////////////////////////////////////////////////////////////////////////////
DeNovoArchitectOP
HelixArchitectCreator::create_architect( std::string const & architect_id ) const
{
	return DeNovoArchitectOP( new HelixArchitect( architect_id ) );
}

std::string
HelixArchitectCreator::keyname() const
{
	return HelixArchitect::class_name();
}

} //protocols
} //denovo_design
} //architects
