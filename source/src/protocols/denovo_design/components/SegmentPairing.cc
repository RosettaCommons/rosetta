// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/components/SegmentPairing.cc
/// @brief Handles user-specified pairing between/among segments
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit header
#include <protocols/denovo_design/components/SegmentPairing.hh>

// Protocol headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// Boost headers
#include <boost/assign.hpp>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.components.SegmentPairing" );

namespace protocols {
namespace denovo_design {
namespace components {

std::string
	SegmentPairing::TAG_NAME = "Pairing";

std::string
SegmentPairing::type_to_str( PairingType const & type )
{
	if ( type == HELIX ) return HelixPairing::class_name();
	if ( type == STRAND ) return StrandPairing::class_name();
	return "UNKNOWN";
}

SegmentPairingOP
create_segment_pairing( std::string const & type_name )
{
	if ( type_name == StrandPairing::class_name() ) return StrandPairingOP( new StrandPairing );
	if ( type_name == HelixPairing::class_name() ) return HelixPairingOP( new HelixPairing );
	std::stringstream msg;
	msg << "Unknown pairing type: " << type_name << std::endl;
	utility_exit_with_message( msg.str() );
	return SegmentPairingOP();
}

SegmentPairing::SegmentPairing():
	utility::pointer::ReferenceCount(),
	segments_()
{
}

SegmentPairing::SegmentPairing( SegmentNames const & segment_names ):
	utility::pointer::ReferenceCount(),
	segments_( segment_names )
{
}

void
SegmentPairing::parse_my_tag( utility::tag::Tag const & tag )
{
	std::string const segment_str = tag.getOption< std::string >( "segments", "" );
	if ( !segment_str.empty() ) set_segments( segment_str );
	parse_tag( tag );
}

SegmentNames const &
SegmentPairing::segments() const
{
	return segments_;
}

void
SegmentPairing::set_segments( std::string const & segment_str )
{
	utility::vector1< std::string > const segments = utility::string_split( segment_str, ',' );
	set_segments( SegmentNames( segments.begin(), segments.end() ) );
}

void
SegmentPairing::set_segments( SegmentNames const & segments )
{
	segments_ = segments;
}

std::ostream &
operator<<( std::ostream & os, SegmentPairing const & pairing )
{
	utility::tag::Tag tag;
	tag.setName( SegmentPairing::TAG_NAME );

	std::string const & type_str = SegmentPairing::type_to_str( pairing.type() );
	tag.setOption< std::string >( "type", type_str );

	std::stringstream seg_str;
	for ( SegmentNames::const_iterator s=pairing.segments().begin(); s!=pairing.segments().end(); ++s ) {
		if ( s != pairing.segments().begin() ) seg_str << ',';
		seg_str << *s;
	}
	tag.setOption< std::string >( "segments", seg_str.str() );

	pairing.to_xml( tag );
	os << tag;
	return os;
}

HelixPairing::HelixPairing():
	SegmentPairing(),
	parallel_( true )
{
}

HelixPairing::HelixPairing( SegmentName const & h1, SegmentName const & h2, bool const is_parallel ):
	SegmentPairing( boost::assign::list_of (h1)(h2).convert_to_container< SegmentNames >() ),
	parallel_( is_parallel )
{
}

SegmentPairingOP
HelixPairing::clone() const
{
	return SegmentPairingOP( new HelixPairing( *this ) );
}

void
HelixPairing::parse_tag( utility::tag::Tag const & tag )
{
	// check segments -- there should be exactly two
	if ( segments().size() != 2 ) {
		std::stringstream msg;
		msg << "You must specify at least two segments to HelixPairing.  Specified segments = "
			<< segments() << std::endl;
		throw utility::excn::EXCN_RosettaScriptsOption( msg.str() );
	}
	parallel_ = tag.getOption< bool >( "parallel" );
}

std::string
HelixPairing::pairing_string( StructureData const & sd ) const
{
	using protocols::fldsgn::topology::SS_Info2;
	SS_Info2 ss_info( sd.ss() );

	std::string const & segment1 = *segments().begin();
	std::string const & segment2 = *(++segments().begin());
	core::Size const h1 = ss_info.helix_id( sd.segment( segment1 ).safe() );
	core::Size const h2 = ss_info.helix_id( sd.segment( segment2 ).safe() );

	std::stringstream pair_str;
	if ( h1 < h2 ) pair_str << h1 << '-' << h2 << '.';
	else pair_str << h2 << '-' << h1 << '.';

	if ( parallel_ ) pair_str << 'P';
	else pair_str << 'A';

	return pair_str.str();
}

void
HelixPairing::to_xml( utility::tag::Tag & tag ) const
{
	tag.setOption< bool >( "parallel", parallel_ );
}

StrandPairing::StrandPairing():
	SegmentPairing(),
	orient1_( architects::UP ),
	orient2_( architects::DOWN ),
	shift_( 0 )
{
}

StrandPairing::StrandPairing(
	SegmentName const & s1,
	SegmentName const & s2,
	architects::StrandOrientation const & orient1,
	architects::StrandOrientation const & orient2,
	architects::RegisterShift const & shift ):
	SegmentPairing( boost::assign::list_of (s1)(s2).convert_to_container< SegmentNames >() ),
	orient1_( orient1 ),
	orient2_( orient2 ),
	shift_( shift )
{
}

SegmentPairingOP
StrandPairing::clone() const
{
	return SegmentPairingOP( new StrandPairing( *this ) );
}

void
StrandPairing::parse_tag( utility::tag::Tag const & tag )
{
	orient1_ = architects::StrandArchitect::int_to_orientation( tag.getOption< core::Size >( "orient1", orient1_ ) );
	orient2_ = architects::StrandArchitect::int_to_orientation( tag.getOption< core::Size >( "orient2", orient2_ ) );
	shift_ = tag.getOption< architects::RegisterShift >( "shift", shift_ );
}

void
StrandPairing::to_xml( utility::tag::Tag & tag ) const
{
	tag.setOption< architects::StrandOrientation >( "orient1", orient1_ );
	tag.setOption< architects::StrandOrientation >( "orient2", orient2_ );
	tag.setOption< architects::RegisterShift >( "shift", shift_ );
}

bool
StrandPairing::parallel() const
{
	return ( orient1_ == orient2_ );
}

std::string
StrandPairing::pairing_string( StructureData const & ) const
{
	return "IMPLEMENT_ME";
}

} //protocols
} //denovo_design
} //components
