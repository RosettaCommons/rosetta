// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/architects/StrandArchitect.cc
/// @brief Architect that creates a beta strand
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/architects/StrandArchitect.hh>
#include <protocols/denovo_design/architects/StrandArchitectCreator.hh>
#include <protocols/denovo_design/architects/DeNovoArchitectFactory.hh>
// Protocol headers
#include <protocols/denovo_design/architects/StructureArchitect.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/util.hh>

// Basic/Utililty headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>
#include <utility/tag/XMLSchemaGeneration.hh>
// Boost headers
#include <boost/assign.hpp>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.architects.StrandArchitect" );

namespace protocols {
namespace denovo_design {
namespace architects {

StrandArchitect::StrandArchitect( std::string const & id_value ):
	DeNovoArchitect( id_value ),
	motifs_(),
	lengths_(),
	bulges_(),
	updated_( false )
{
}

StrandArchitect::~StrandArchitect()
{}

StrandArchitect::DeNovoArchitectOP
StrandArchitect::clone() const
{
	return DeNovoArchitectOP( new StrandArchitect( *this ) );
}

std::string
StrandArchitect::type() const
{
	return StrandArchitect::class_name();
}

void
StrandArchitectCreator::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) const {
	StrandArchitect::provide_xml_schema( xsd );
}

void
StrandArchitect::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & )
{
	std::string const lengths_str = tag->getOption< std::string >( "length", "" );
	if ( !lengths_str.empty() ) set_length( lengths_str );

	std::string const bulge_str = tag->getOption< std::string >( "bulge", "" );
	if ( !bulge_str.empty() ) set_bulges( bulge_str );

	if ( ! updated_ ) enumerate_permutations();
}

void
StrandArchitect::provide_xml_schema( utility::tag::XMLSchemaDefinition & xsd ) {
	using namespace utility::tag;

	//Define the length restriction and the bulge restriction
	XMLSchemaRestriction bulge_string;
	bulge_string.name( "bulge_string" );
	bulge_string.base_type( xs_string );
	bulge_string.add_restriction( xsr_pattern, "[0-9]+([,;][0-9])+" );
	xsd.add_top_level_element( bulge_string );

	XMLSchemaRestriction register_shift_string;
	register_shift_string.name( "register_shift_string" );
	register_shift_string.base_type( xs_string );
	register_shift_string.add_restriction( xsr_pattern, "-?[0-9]+([,:]-?[0-9])+" );
	xsd.add_top_level_element( register_shift_string );

	XMLSchemaRestriction orientation_string;
	orientation_string.name( "orientation_string" );
	orientation_string.base_type( xs_string );
	orientation_string.add_restriction( xsr_pattern, "[U,D](,[U,D])*(:[U,D](,[U,D])*)*" );
	xsd.add_top_level_element( orientation_string );

	XMLSchemaRestriction length_string;
	length_string.name( "strand_length_string" );
	length_string.base_type( xs_string );
	length_string.add_restriction( xsr_pattern, "[0-9]+(:[0-9]+)?(,[0-9]+(:[0-9]+)?)*" );
	xsd.add_top_level_element( length_string );

	AttributeList attlist;
	attlist
		+ XMLSchemaAttribute( "length", "strand_length_string", "Comma-separated list of single integers and hyphen-separated ranges to specify all possible strand lengths" )
		+ XMLSchemaAttribute( "bulge", "bulge_string", "Specifies where bulges occur in a strand" )
		+ XMLSchemaAttribute( "register_shift", "register_shift_string", "Specifies what register shifts are to be used.  This option is only used in the context of the BetaSheetArchitect." )
		+ XMLSchemaAttribute( "orientation", "orientation_string", "Specifies the orientation of strands in a beta sheet.  This option is only used in the context of the BetaSheetArchitect." );

	DeNovoArchitect::add_common_denovo_architect_attributes( attlist );
	DeNovoArchitectFactory::xsd_architect_type_definition_w_attributes( xsd, class_name(), "Architect to construct a beta strand", attlist );

}

StrandArchitect::StructureDataOP
StrandArchitect::design( core::pose::Pose const &, core::Real & random ) const
{
	check_updated();
	if ( motifs_.size() == 0 ) {
		std::stringstream msg;
		msg << class_name() << "::design(): List of possible strand permutations is empty!"
			<< std::endl;
		utility_exit_with_message( msg.str() );
	}

	core::Size const idx = extract_int( random, 1, motifs_.size() );
	return StructureDataOP( new StructureData( *motifs_[ idx ] ) );
}

/// @brief Given a list of bulges, build a motif
/// @param[in] bulges    List of bulge positions to place onto the strand
/// @param[in] secstruct Secondary structure for the new segment
/// @param[in] abego     Abego for the new segment, without bulges placed
StrandArchitect::StructureDataOP
StrandArchitect::create_motif(
	StrandBulges const & bulges,
	std::string const & secstruct,
	std::string const & abego ) const
{
	StructureDataOP sd( new StructureData( id() ) );
	std::string abego_str = abego;
	for ( auto const & b : bulges ) {
		if ( b ) abego_str[ b ] = 'A';
	}
	sd->add_segment( components::Segment( id(), secstruct, abego_str, false, false ) );
	store_bulges( *sd, bulges );
	return sd;
}

components::StructureDataCOPs
StrandArchitect::compute_permutations() const
{
	components::StructureDataCOPs motifs;
	for ( core::Size const & l : lengths_ ) {
		std::stringstream secstruct;
		secstruct << 'L' << std::string( l, 'E' ) << 'L';
		std::stringstream abego;
		abego << 'X' << std::string( l, 'B' ) << 'X';

		if ( bulges_.empty() ) {
			motifs.push_back( create_motif( StrandBulges(), secstruct.str(), abego.str() ) );
		} else {
			for ( StrandBulges const & bulgelist : bulges_ ) {
				motifs.push_back( create_motif( bulgelist, secstruct.str(), abego.str() ) );
			}
		}
	}
	if ( motifs.empty() ) {
		std::stringstream msg;
		msg << "StrandArchitect: no strand permutations could be generated with the given user input." << std::endl;
		msg << "Lengths: " << lengths_ << std::endl;
		msg << "Bulges: " << bulges_ << std::endl;
		utility_exit_with_message( msg.str() );
	}
	return motifs;
}

void
StrandArchitect::enumerate_permutations()
{
	motifs_ = compute_permutations();
	updated_ = true;
}

/// @brief If architect is updated, simply exit without error.  If not, exit with an error message
void
StrandArchitect::check_updated() const
{
	if ( !updated_ ) {
		std::stringstream msg;
		msg << "StrandArchitect: Motif list needs updating, but it is being accessed. You probably need to call enumerate_permutations()" << std::endl;
		utility_exit_with_message( msg.str() );
	}
}

components::StructureDataCOPs::const_iterator
StrandArchitect::motifs_begin() const
{
	check_updated();
	return motifs_.begin();
}

components::StructureDataCOPs::const_iterator
StrandArchitect::motifs_end() const
{
	check_updated();
	return motifs_.end();
}

StrandArchitect::StrandOrientation
StrandArchitect::int_to_orientation( int const orient )
{
	if ( orient == 0 ) {
		std::stringstream msg;
		msg << class_name() << ": int_to_orientation(): Found orientation inside StructureData with value ("
			<< orient << "), but orientations must be between 1 and " << components::ORIENTATIONS_END - 1 << "." << std::endl;
		utility_exit_with_message( msg.str() );
	}
	if ( orient >= components::ORIENTATIONS_END ) {
		std::stringstream msg;
		msg << class_name() << ": int_to_orientation(): Found orientation inside StructureData with value ("
			<< orient << "), but there are only " << components::ORIENTATIONS_END - 1 << " possible orientations." << std::endl;
		utility_exit_with_message( msg.str() );
	}
	return static_cast< StrandOrientation >( orient );
}

StrandBulges
StrandArchitect::retrieve_bulges( StructureData const & sd ) const
{
	return retrieve_bulges_s( sd, id() );
}

StrandBulges
StrandArchitect::retrieve_bulges_s( StructureData const & sd, std::string const & segment_id )
{
	std::string const value = sd.get_data_str( segment_id, bulge_keyname() );
	utility::vector1< std::string > const fields = utility::string_split( value, ';' );
	StrandBulges bulges;
	for ( std::string const & position : fields ) {
		if ( position.empty() ) continue;
		bulges.push_back( boost::lexical_cast< SegmentResid >( position ) );
	}
	return bulges;
}

void
StrandArchitect::store_bulges( StructureData & sd, StrandBulges const & bulges ) const
{
	std::stringstream storage;
	if ( bulges.empty() ) {
		sd.set_data_str( id(), bulge_keyname(), "" );
	} else {
		for ( SegmentResid const r : bulges ) {
			if ( !storage.str().empty() ) storage << ';';
			storage << r;
		}
		sd.set_data_str( id(), bulge_keyname(), storage.str() );
	}
}

std::string const
StrandArchitect::bulge_keyname()
{
	return "bulge";
}

void
StrandArchitect::set_length( std::string const & length_str )
{
	set_length( parse_length_str< core::Size >( length_str ) );
}

void
StrandArchitect::set_length( Lengths const & lengths_val )
{
	lengths_ = lengths_val;
	needs_update();
}

void
StrandArchitect::set_bulges( std::string const & bulges_str )
{
	// bulges string if of format LENGTHS_B1;LENGTHS_B2;...;LENGTHS_BN
	// LENGTHS_B1 .. LENGTHS_BN are lengths strings
	// e.g. "2,3,4;5:8;10" means there are three bulges.  The first is at 2, 3 or 4.
	//      the second is at 5, 6, 7 or 8.  The third is at 10.
	AllowedStrandBulges allowed_bulges;
	utility::vector1< std::string > const fields = utility::string_split( bulges_str, ';' );
	for ( std::string const & bulge_positions : fields ) {
		allowed_bulges.push_back( parse_length_str< SegmentResid >( bulge_positions ) );
	}

	set_bulges( allowed_bulges );
}

void
StrandArchitect::set_bulges( AllowedStrandBulges const & bulges )
{
	bulges_ = bulges;
	needs_update();
}

void
StrandArchitect::needs_update()
{
	updated_ = false;
}

///////////////////////////////////////////////////////////////////////////////
/// Creator
///////////////////////////////////////////////////////////////////////////////
std::string
StrandArchitectCreator::keyname() const
{
	return StrandArchitect::class_name();
}

DeNovoArchitectOP
StrandArchitectCreator::create_architect( std::string const & architect_id ) const
{
	return DeNovoArchitectOP( new StrandArchitect( architect_id ) );
}

} //protocols
} //denovo_design
} //architects
