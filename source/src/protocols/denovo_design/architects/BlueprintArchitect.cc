// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/architects/BlueprintArchitect.cc
/// @brief Designs a structure using a Blueprint file
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/architects/BlueprintArchitect.hh>
#include <protocols/denovo_design/architects/BlueprintArchitectCreator.hh>
#include <protocols/denovo_design/architects/DeNovoArchitectFactory.hh>
// Protocol headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/util.hh>
#include <protocols/jd2/parser/BluePrint.hh>

// Basic/Utililty headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.architects.BlueprintArchitect" );

namespace protocols {
namespace denovo_design {
namespace architects {

BlueprintArchitect::BlueprintArchitect( std::string const & id_value ):
	DeNovoArchitect( id_value ),
	blueprint_()
{
}

BlueprintArchitect::~BlueprintArchitect()
{}

DeNovoArchitectOP
BlueprintArchitect::clone() const
{
	return DeNovoArchitectOP( new BlueprintArchitect( *this ) );
}

std::string
BlueprintArchitect::type() const
{
	return BlueprintArchitect::class_name();
}

void
BlueprintArchitect::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & )
{
	std::string const bp_file = tag->getOption< std::string >( "blueprint", "" );
	if ( bp_file.empty() ) {
		std::stringstream msg;
		msg << "BlueprintArchitect: No blueprint file specified!  "
			<< "You must specify a blueprint file using the \"blueprint\" option" << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
	}
	protocols::jd2::parser::BluePrint bp( bp_file );
	set_blueprint( bp );
}
void
BlueprintArchitect::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ){
	using namespace utility::tag;
	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute::required_attribute( "blueprint", xs_string, "Path to blueprint file" );
	DeNovoArchitect::add_common_denovo_architect_attributes( attlist );
	DeNovoArchitectFactory::xsd_architect_type_definition_w_attributes( xsd, class_name(), "Architect that constructs pose using a blueprint", attlist );

}
void
BlueprintArchitectCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const{
	BlueprintArchitect::provide_xml_schema( xsd );
}

DeNovoArchitect::StructureDataOP
BlueprintArchitect::design( core::pose::Pose const & pose, core::Real & ) const
{
	if ( !blueprint_ ) {
		std::stringstream msg;
		msg << "BluePrintArchitect::design(): No blueprint specified!" << std::endl;
		utility_exit_with_message( msg.str() );
	}

	StructureDataOP sd( new StructureData( id() ) );
	components::SegmentCounts counts;
	components::add_segments_for_chain( id(), *sd, blueprint_->secstruct(), abego_str( blueprint_->abego() ), counts );

	set_template_segments( *sd, pose );

	// find pairings
	set_helix_pairings( *sd );
	//std::string const strand_pairings = blueprint_->strand_pairings();
	//std::string const hss_triplets = blueprint_->hss_triplets();


	return sd;
}

void
BlueprintArchitect::set_blueprint( protocols::jd2::parser::BluePrint const & bp )
{
	blueprint_ = protocols::jd2::parser::BluePrintCOP( new protocols::jd2::parser::BluePrint( bp ) );
}

/// @brief Adds helix pairings to the given SD using the HHPAIR line of the blueprint
/// @param[in,out] sd StructureData object to be modified
void
BlueprintArchitect::set_helix_pairings( StructureData & ) const
{
	std::string const helix_pairings = blueprint_->helix_pairings();
	if ( helix_pairings.empty() ) return;

	core::Size h1, h2;
	char orientation;
	get_helix_pairings( helix_pairings, h1, h2, orientation );
	TR << "Found paring h1=" << h1 << " h2=" << h2 << " or=" << orientation << std::endl;
}

void
invalid_hhpair( std::string const & hhpair_str )
{
	std::stringstream msg;
	msg << "BlueprintArchitect::get_helix_pairings(): invalid HHPAIR string in blueprint: "
		<< hhpair_str << std::endl;
	utility_exit_with_message( msg.str() );
}

/// @brief Converts HHPAIR data in the blueprint into usable information
/// @param[in]  hhpair_str  Blueprint HHPAIR string
/// @param[out] h1          Helix number (N-->C ordering) for first paired helix
/// @param[out] h2          Helix number (N-->C ordering) for first paired helix
/// @param[out] orientation Orientation of the helices (P or A)
void
BlueprintArchitect::get_helix_pairings(
	std::string const & hhpair_str,
	core::Size & h1,
	core::Size & h2,
	char & orientation ) const
{
	typedef utility::vector1< std::string > Strings;
	Strings const fields = utility::string_split( hhpair_str, '.' );

	if ( fields.size() != 2 ) {
		invalid_hhpair( hhpair_str );
	}

	if ( fields.rbegin()->size() != 1 ) {
		invalid_hhpair( hhpair_str );
	}

	Strings const helices = utility::string_split( *fields.begin(), '-' );

	if ( helices.size() != 2 ) {
		invalid_hhpair( hhpair_str );
	}

	h1 = boost::lexical_cast< core::Size >( *helices.begin() );
	h2 = boost::lexical_cast< core::Size >( *helices.rbegin() );
	orientation = *fields.rbegin()->begin();
}

/// @brief add templated segments from pose. These are designated by a residue number
/// in the blueprint instead of 0
/// @param[in,out] sd    StructureData to be modified
/// @param[in]     pose  Pose containing template residues
void
BlueprintArchitect::set_template_segments( StructureData & sd, core::pose::Pose const & pose ) const
{
	SegmentNames const templates = get_template_segments( sd );
	for ( SegmentNames::const_iterator s=templates.begin(); s!=templates.end(); ++s ) {
		TR << "Adding template for " << *s << std::endl;
		sd.set_template_pose( *s, pose, blueprint_->resnum( sd.segment( *s ).lower() ), blueprint_->resnum( sd.segment( *s ).upper() ) );
	}
}

/// @brief gets names of templated segments from pose. These are designated by a residue number
/// in the blueprint instead of 0.  All residues in the segment must be sequential and not built
/// denovo
/// @param[in] sd    StructureData to be scanned
/// @returns SegmentNameSet of segments with templates
SegmentNames
BlueprintArchitect::get_template_segments( StructureData const & sd ) const
{
	SegmentNames has_template;
	for ( SegmentNameList::const_iterator s=sd.segments_begin(); s!=sd.segments_end(); ++s ) {
		core::Size template_residues = 0;
		core::Size prev_pose_residue = 0;
		for ( core::Size resid=sd.segment( *s ).lower(); resid<=sd.segment( *s ).upper(); ++resid ) {
			core::Size const pose_residue = blueprint_->resnum( resid );
			if ( pose_residue == 0 ) {
				prev_pose_residue = 0;
				continue;
			}

			if ( ( prev_pose_residue > 0 ) && ( prev_pose_residue + 1 != pose_residue ) ) {
				std::stringstream msg;
				msg << class_name() << "::get_template_segments(): Pose residue numbering is not sequential! "
					<< "Previous residue: " << prev_pose_residue << " Current residue: " << pose_residue
					<< std::endl;
				utility_exit_with_message( msg.str() );
			}

			++template_residues;
		}

		if ( template_residues == 0 ) continue;

		if ( template_residues != sd.segment( *s ).length() ) {
			bool contains_fixed = false;
			for ( core::Size resid=sd.segment( *s ).lower(); resid<=sd.segment( *s ).upper(); ++resid ) {
				if ( blueprint_->buildtype( resid ) != '.' ) continue;
				contains_fixed = true;
				break;
			}

			// this is not a template segments if it doesn't contain fixed residues
			// we assume it can be built denovo safely
			if ( !contains_fixed ) continue;

			std::stringstream msg;
			msg << class_name() << "::get_template_segments(): A secondary structure element ("
				<< *s << ") is not either fully taken from the input pose or fully de-novo. "
				<< "It contains " << template_residues << " residues from the input pose and "
				<< sd.segment( *s ).length() << " residues total." << std::endl;
			utility_exit_with_message( msg.str() );
		}
		has_template.push_back( *s );
	}
	return has_template;
}

///////////////////////////////////////////////////////////////////////////////
/// Creator
///////////////////////////////////////////////////////////////////////////////
DeNovoArchitectOP
BlueprintArchitectCreator::create_architect( std::string const & architect_id ) const
{
	return DeNovoArchitectOP( new BlueprintArchitect( architect_id ) );
}

std::string
BlueprintArchitectCreator::keyname() const
{
	return BlueprintArchitect::class_name();
}

} //protocols
} //denovo_design
} //architects
