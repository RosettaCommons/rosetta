// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/architects/StrandArchitect.cc
/// @brief Architect that creates a beta strand
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/architects/StrandArchitect.hh>
#include <protocols/denovo_design/architects/StrandArchitectCreator.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>

// Basic/Utililty headers
#include <basic/Tracer.hh>
#include <utility/string_util.hh>
#include <utility/tag/Tag.hh>

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
	orientations_( 1, UP ),
	register_shifts_( 1, 0 ),
	paired_strands_( "", "" ),
	bulges_( 1, 0 ),
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
StrandArchitect::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & )
{
	std::string const lengths_str = tag->getOption< std::string >( "length", "" );
	if ( !lengths_str.empty() ) set_length( lengths_str );

	std::string const orientation_str = tag->getOption< std::string >( "orientation", "" );
	if ( !orientation_str.empty() ) set_orientation( orientation_str );

	std::string const register_shift_str = tag->getOption< std::string >( "register_shift", "" );
	if ( !register_shift_str.empty() ) set_register_shift( register_shift_str );

	std::string const bulge_str = tag->getOption< std::string >( "bulge", "" );
	if ( !bulge_str.empty() ) set_bulge( bulge_str );

	if ( ! updated_ ) enumerate_permutations();
}

StrandArchitect::StructureDataOP
StrandArchitect::design( core::pose::Pose const &, core::Real & ) const
{
	return StructureDataOP();
}

void
StrandArchitect::set_paired_strands( PairedStrandNames const & strands )
{
	paired_strands_ = strands;
	needs_update();
}

components::StructureDataCOPs
StrandArchitect::compute_permutations() const
{
	components::StructureDataCOPs motifs;
	for ( Lengths::const_iterator l=lengths_.begin(); l!=lengths_.end(); ++l ) {
		std::stringstream secstruct;
		secstruct << 'L' << std::string( *l, 'E' ) << 'L';
		std::stringstream abego;
		abego << 'X' << std::string( *l, 'B' ) << 'X';

		for ( StrandOrientations::const_iterator o=orientations_.begin(); o!=orientations_.end(); ++o ) {
			for ( RegisterShifts::const_iterator s=register_shifts_.begin(); s!=register_shifts_.end(); ++s ) {
				for ( StrandBulges::const_iterator b=bulges_.begin(); b!=bulges_.end(); ++b ) {
					StructureDataOP sd( new StructureData( id() ) );
					std::string abego_str = abego.str();
					if ( *b != 0 ) {
						abego_str[ *b ] = 'A';
					}
					sd->add_segment( id(), components::Segment( secstruct.str(), abego_str, false, false ) );
					store_register_shift( *sd, *s );
					store_bulge( *sd, *b );
					store_orientation( *sd, *o );
					store_paired_strands( *sd, paired_strands_ );
					motifs.push_back( sd );
				}
			}
		}
	}
	if ( motifs.empty() ) {
		std::stringstream msg;
		msg << "StrandArchitect: no strand permutations could be generated with the given user input." << std::endl;
		msg << "Lengths: " << lengths_ << std::endl;
		msg << "Orientations: " << orientations_ << std::endl;
		msg << "Shifts: " << register_shifts_ << std::endl;
		msg << "Bulges: " << bulges_ << std::endl;
		msg << "Paired Strands: " << paired_strands_ << std::endl;
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

components::StructureDataCOPs::const_iterator
StrandArchitect::motifs_begin() const
{
	if ( !updated_ ) {
		std::stringstream msg;
		msg << "StrandArchitect: Motif list needs updating, but it is being accessed. You probably need to call enumerate_permutations()" << std::endl;
		utility_exit_with_message( msg.str() );
	}
	return motifs_.begin();
}

components::StructureDataCOPs::const_iterator
StrandArchitect::motifs_end() const
{
	if ( !updated_ ) {
		std::stringstream msg;
		msg << "StrandArchitect: Motif list needs updating, but it is being accessed. You probably need to call enumerate_permutations()" << std::endl;
		utility_exit_with_message( msg.str() );
	}
	return motifs_.end();
}

/// @brief retrieve register shift
RegisterShift
StrandArchitect::retrieve_register_shift( StructureData const & sd ) const
{
	return sd.get_data_int( id(), register_shift_keyname() );
}

StrandOrientation
StrandArchitect::int_to_orientation( int const orient )
{
	if ( orient >= ORIENTATIONS_END ) {
		std::stringstream msg;
		msg << class_name() << ": int_to_orientation(): Found orientation inside StructureData with value ("
			<< orient << "), but there are only " << ORIENTATIONS_END - 1 << " possible orientations." << std::endl;
		utility_exit_with_message( msg.str() );
	}
	return static_cast< StrandOrientation >( orient );
}

/// @brief retrieve orientation
StrandOrientation
StrandArchitect::retrieve_orientation( StructureData const & sd ) const
{
	return int_to_orientation( sd.get_data_int( id(), orientation_keyname() ) );
}

StrandBulge
StrandArchitect::retrieve_bulge( StructureData const & sd ) const
{
	return static_cast< StrandBulge >( sd.get_data_int( id(), bulge_keyname() ) );
}

void
StrandArchitect::store_register_shift( StructureData & sd, RegisterShift const shift ) const
{
	sd.set_data_int( id(), register_shift_keyname(), shift );
}

void
StrandArchitect::store_orientation( StructureData & sd, StrandOrientation const orient ) const
{
	sd.set_data_int( id(), orientation_keyname(), orient );
}

void
StrandArchitect::store_bulge( StructureData & sd, StrandBulge const bulge ) const
{
	sd.set_data_int( id(), bulge_keyname(), bulge );
}

std::string const
StrandArchitect::register_shift_keyname()
{
	return "register_shift";
}

std::string const
StrandArchitect::orientation_keyname()
{
	return "orientation";
}


std::string const
StrandArchitect::bulge_keyname()
{
	return "bulge";
}

std::string const
StrandArchitect::paired_strands_keyname()
{
	return "paired_strands";
}

PairedStrandNames
StrandArchitect::retrieve_paired_strands( StructureData const & sd ) const
{
	std::string const & pair_str = sd.get_data_str( id(), paired_strands_keyname() );
	return str_to_paired_strands( pair_str );
}

PairedStrandNames
StrandArchitect::str_to_paired_strands( std::string const & pair_str )
{
	utility::vector1< std::string > const paired_strands = utility::string_split( pair_str, ',' );
	if ( paired_strands.size() != 2 ) {
		std::stringstream msg;
		msg << "StrandArchitect: paired strands list must contain exactly two comma-separated segment names."
			<< pair_str << " is the paired_strands string found." << std::endl;
		utility_exit_with_message( msg.str() );
	}
	return PairedStrandNames( paired_strands[1], paired_strands[2] );
}

void
StrandArchitect::store_paired_strands( StructureData & sd, PairedStrandNames const & strands ) const
{
	std::stringstream pair_str;
	pair_str << strands;
	sd.set_data_str( id(), paired_strands_keyname(), pair_str.str() );
}

/// @brief set register shift
void
StrandArchitect::set_register_shift( RegisterShifts const & shifts )
{
	register_shifts_ = shifts;
	needs_update();
}

/// @brief set allowed register shifts from a string
void
StrandArchitect::set_register_shift( std::string const & val )
{
	RegisterShifts retval;
	utility::vector1< std::string > const str_shifts( utility::string_split( val, ',' ) );
	for ( utility::vector1< std::string >::const_iterator s=str_shifts.begin(); s!=str_shifts.end(); ++s ) {
		TR.Debug << *s << " " << val << std::endl;
		if ( s->empty() ) continue;
		utility::vector1< std::string > const ranges( utility::string_split( *s, ':' ) );
		if ( ranges.size() == 1 ) {
			retval.push_back( boost::lexical_cast< RegisterShift >( ranges[1] ) );
		} else if ( ranges.size() == 2 ) {
			RegisterShift const start( boost::lexical_cast< RegisterShift >( ranges[1] ) );
			RegisterShift const end( boost::lexical_cast< RegisterShift >( ranges[2] ) );
			for ( RegisterShift i=start; i<=end; ++i ) {
				retval.push_back( i );
			}
		} else {
			utility_exit_with_message( "Invalid register shift input: " + val );
		}
	}
	set_register_shift( retval );
}

/// @brief set allowed orientations from a string
void
StrandArchitect::set_orientation( std::string const & val )
{
	StrandOrientations retval;
	utility::vector1< std::string > const str_orients( utility::string_split( val, ',' ) );
	for ( utility::vector1< std::string >::const_iterator s=str_orients.begin(); s!=str_orients.end(); ++s ) {
		TR.Debug << *s << " " << val << std::endl;
		if ( s->empty() ) continue;
		utility::vector1< std::string > const ranges( utility::string_split( *s, ':' ) );
		// check input string to make sure only "A", "P" are specified
		for ( core::Size i=1; i<=ranges.size(); ++i ) {
			if ( ( ranges[i] != "U" ) && ( ranges[i] != "D" ) ) {
				utility_exit_with_message( "Invalid orientation character: " + ranges[i] );
			}
		}
		if ( ranges.size() == 1 ) {
			if ( ranges[1] == "U" ) {
				retval.push_back( UP );
			} else {
				retval.push_back( DOWN );
			}
		} else if ( ranges.size() == 2 ) {
			retval.push_back( UP );
			retval.push_back( DOWN );
		} else {
			utility_exit_with_message( "Invalid orientation input: " + val );
		}
	}
	set_orientation( retval );
}

void
StrandArchitect::set_orientation( StrandOrientations const & orientations )
{
	orientations_ = orientations;
	needs_update();
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
StrandArchitect::set_bulge( std::string const & bulges_str )
{
	set_bulge( parse_length_str< StrandBulge >( bulges_str ) );
}

void
StrandArchitect::set_bulge( StrandBulges const & bulges )
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

std::ostream &
operator<<( std::ostream & os, PairedStrandNames const & paired_strands )
{
	os << paired_strands.segment1 << ',' << paired_strands.segment2;
	return os;
}

///////////////////////////////////////////////////////////////////////////////
/// PairedStrandNames
///////////////////////////////////////////////////////////////////////////////
bool
PairedStrandNames::operator<( PairedStrandNames const & other ) const
{
	if ( segment1 < other.segment1 ) return true;
	if ( segment1 == other.segment1 ) {
		return ( segment2 < other.segment2 );
	}
	return false;
}

} //protocols
} //denovo_design
} //architects
