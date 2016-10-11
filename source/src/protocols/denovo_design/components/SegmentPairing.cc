// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/denovo_design/components/SegmentPairing.cc
/// @brief Handles user-specified pairing between/among segments
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit header
#include <protocols/denovo_design/components/SegmentPairing.hh>

// Protocol headers
#include <protocols/denovo_design/architects/StrandArchitect.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>

// Core headers
#include <core/select/residue_selector/ResidueVector.hh>

// Basic headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// Boost headers
#include <boost/assign.hpp>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.components.SegmentPairing" );

namespace protocols {
namespace denovo_design {
namespace components {

std::string SegmentPairing::TAG_NAME = "Pairing";

std::string
SegmentPairing::type_to_str( PairingType const & type )
{
	if ( type == HELIX ) return HelixPairing::class_name();
	if ( type == STRAND ) return StrandPairing::class_name();
	if ( type == HELIX_SHEET ) return HelixSheetPairing::class_name();
	return "UNKNOWN";
}

SegmentPairingOP
SegmentPairing::create( std::string const & type_name )
{
	if ( type_name == StrandPairing::class_name() ) return StrandPairingOP( new StrandPairing );
	if ( type_name == HelixPairing::class_name() ) return HelixPairingOP( new HelixPairing );
	if ( type_name == HelixSheetPairing::class_name() ) return HelixSheetPairingOP( new HelixSheetPairing );
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

bool
SegmentPairing::has_segment( std::string const & segment ) const
{
	for ( std::string const & s : segments_ ) {
		if ( s == segment ) return true;
	}
	return false;
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
	tag.set_quote_options( true );
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

core::select::residue_selector::ResidueVector
get_strand_residues(
	StrandOrientation const & orient,
	core::Size const start,
	core::Size const stop,
	core::select::residue_selector::ResidueVector const & bulges )
{
	TR.Debug << "Bulges are at " << bulges << std::endl;
	std::set< core::Size > const bulge_set( bulges.begin(), bulges.end() );

	core::select::residue_selector::ResidueVector resids;
	if ( orient == UP ) {
		for ( core::Size resid=start; resid<=stop; ++resid ) {
			if ( bulge_set.find( resid ) != bulge_set.end() ) continue;
			resids.push_back( resid );
		}
	} else if ( orient == DOWN ) {
		for ( core::Size resid=stop; resid>=start; --resid ) {
			if ( bulge_set.find( resid ) != bulge_set.end() ) continue;
			resids.push_back( resid );
		}
	} else {
		std::stringstream msg;
		msg << "mark_paired_residues(): Invalid orientation found for strand ["
			<< start << "," << stop << "] : "  << orient << std::endl;
		utility_exit_with_message( msg.str() );
	}

	return resids;
}

core::select::residue_selector::ResidueVector
get_bulges( core::Size const start, core::Size const stop, std::string const & abego )
{
	core::select::residue_selector::ResidueVector bulges;
	for ( core::Size resid=start; resid<=stop; ++resid ) {
		if ( abego[resid-1] == 'A' )  bulges.push_back( resid );
	}
	return bulges;
}

void
add_paired_residues(
	StructureData const & sd,
	StrandPairing const & p,
	protocols::fldsgn::topology::SS_Info2 const & ss_info,
	ResiduePairs & pairs )
{
	using core::select::residue_selector::ResidueVector;

	protocols::fldsgn::topology::Strands const & strands = ss_info.strands();

	debug_assert( p.segments().size() == 2 );
	std::string const seg1_name = *p.segments().begin();
	std::string const seg2_name = *p.segments().rbegin();
	std::string const abego = sd.abego();

	core::Size const strand1 = ss_info.strand_id( sd.segment( seg1_name ).safe() );
	core::Size const strand2 = ss_info.strand_id( sd.segment( seg2_name ).safe() );

	StrandOrientation const o1 = p.orient1();
	StrandOrientation const o2 = p.orient2();
	RegisterShift const shift = p.shift();

	ResidueVector const bulges1 = get_bulges( strands[strand1]->begin(), strands[strand1]->end(), abego );
	ResidueVector const bulges2 = get_bulges( strands[strand2]->begin(), strands[strand2]->end(), abego );
	ResidueVector const s1_resids = get_strand_residues( o1, strands[strand1]->begin(), strands[strand1]->end(), bulges1 );
	ResidueVector const s2_resids = get_strand_residues( o2, strands[strand2]->begin(), strands[strand2]->end(), bulges2 );

	for ( core::Size s1_idx=1; s1_idx<=s1_resids.size(); ++s1_idx ) {
		core::Size const s2_idx = s1_idx - shift;
		if ( s2_idx < 1 ) continue;
		if ( s2_idx > s2_resids.size() ) continue;
		// we have a pair
		core::Size const res1 = s1_resids[ s1_idx ];
		core::Size const res2 = s2_resids[ s2_idx ];
		TR.Debug << "Paired: " << res1 << " <--> " << res2 << std::endl;
		pairs.push_back( ResiduePair( res1, res2 ) );
	}
}

ResiduePairs
SegmentPairing::get_strand_residue_pairs( StructureData const & sd )
{
	protocols::fldsgn::topology::SS_Info2 const ss_info( sd.ss() );

	ResiduePairs pairs;
	SegmentPairingCOPs const pairings = get_pairings( sd, STRAND );
	for ( SegmentPairingCOP const & pairing : pairings ) {
		debug_assert( utility::pointer::dynamic_pointer_cast< StrandPairing const >( pairing ) );
		StrandPairingCOP spairing = utility::pointer::static_pointer_cast< StrandPairing const >( pairing );
		add_paired_residues( sd, *spairing, ss_info, pairs );
	}
	return pairs;
}

SegmentPairingCOPs
SegmentPairing::get_pairings( StructureData const & sd, PairingType const & type )
{
	SegmentPairingCOPs my_type;
	for ( SegmentPairingCOPs::const_iterator p=sd.pairings_begin(); p!=sd.pairings_end(); ++p ) {
		if ( (*p)->type() != type ) continue;
		my_type.push_back( *p );
	}
	return my_type;
}

std::string
SegmentPairing::get_pairing_str( StructureData const & sd, PairingType const & type )
{
	std::stringstream pair_str;
	SegmentPairingCOPs const my_type = get_pairings( sd, type );
	for ( SegmentPairingCOPs::const_iterator p=my_type.begin(); p!=my_type.end(); ++p ) {
		if ( !pair_str.str().empty() ) pair_str << ';';
		pair_str << (*p)->pairing_string( sd );
	}
	return pair_str.str();
}

std::string
SegmentPairing::get_strand_pairings( StructureData const & sd )
{
	return get_pairing_str( sd, STRAND );
}

std::string
SegmentPairing::get_helix_pairings( StructureData const & sd )
{
	return get_pairing_str( sd, HELIX );
}

std::string
SegmentPairing::get_hss_triplets( StructureData const & sd )
{
	return get_pairing_str( sd, HELIX_SHEET );
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
	orient1_( UP ),
	orient2_( DOWN ),
	shift_( 0 )
{
}

StrandPairing::StrandPairing(
	SegmentName const & s1,
	SegmentName const & s2,
	StrandOrientation const & orient1,
	StrandOrientation const & orient2,
	RegisterShift const & shift ):
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
	shift_ = tag.getOption< RegisterShift >( "shift", shift_ );
}

void
StrandPairing::to_xml( utility::tag::Tag & tag ) const
{
	tag.setOption< StrandOrientation >( "orient1", orient1_ );
	tag.setOption< StrandOrientation >( "orient2", orient2_ );
	tag.setOption< RegisterShift >( "shift", shift_ );
}

bool
StrandPairing::parallel() const
{
	return ( orient1_ == orient2_ );
}

std::string
StrandPairing::pairing_string( StructureData const & sd ) const
{
	using protocols::fldsgn::topology::SS_Info2;
	std::stringstream pair_str;

	// collect strand numbers
	SS_Info2 const ss_info( sd.ss() );
	core::Size const res_in_s1 = sd.segment( *segments().begin() ).safe();
	core::Size const s1 = ss_info.strand_id( res_in_s1 );
	core::Size const res_in_s2 = sd.segment( *( segments().begin() + 1 ) ).safe();
	core::Size const s2 = ss_info.strand_id( res_in_s2 );

	debug_assert( s1 != s2 );
	if ( s1 < s2 ) pair_str << s1 << '-' << s2;
	else pair_str << s2 << '-' << s1;

	pair_str << '.';
	if ( parallel() ) pair_str << 'P';
	else pair_str << 'A';

	pair_str << '.';

	// find the pairing in n-->c order
	if ( s1 < s2 ) pair_str << nobu_register_shift( sd, *ss_info.strand( s1 ), *ss_info.strand( s2 ), shift_, orient1_ );
	else pair_str << nobu_register_shift( sd, *ss_info.strand( s2 ), *ss_info.strand( s1 ), -shift_, orient2_ );

	return pair_str.str();
}

StrandOrientation
StrandPairing::orient1() const
{
	return orient1_;
}

StrandOrientation
StrandPairing::orient2() const
{
	return orient2_;
}

RegisterShift
StrandPairing::shift() const
{
	return shift_;
}

core::Size
count_bulges(
	protocols::fldsgn::topology::Strand const & strand,
	std::string const & abego )
{
	core::Size bulge_count = 0;
	for ( core::Size resid=strand.begin(); resid<=strand.end(); ++resid ) {
		if ( abego[ resid - 1 ] != 'B' ) ++bulge_count;
	}
	return bulge_count;
}

RegisterShift
StrandPairing::nobu_register_shift(
	StructureData const & sd,
	protocols::fldsgn::topology::Strand const & s1,
	protocols::fldsgn::topology::Strand const & s2,
	RegisterShift const nc_order_shift,
	StrandOrientation const & nc_order_orient ) const
{
	core::Size const bulges1 = count_bulges( s1, sd.abego() );
	core::Size const bulges2 = count_bulges( s2, sd.abego() );

	RegisterShift preliminary_shift;
	// E1 = E2 + X + Y, X = shift_
	if ( nc_order_orient == UP ) {
		// X is all that matters
		preliminary_shift = nc_order_shift;
	} else if ( nc_order_orient == DOWN ) {
		// Y is all that matters
		RegisterShift const len1 = s1.length() - bulges1;
		RegisterShift const len2 = s2.length() - bulges2;
		RegisterShift const flipped_shift = len1 - ( len2 + nc_order_shift );
		preliminary_shift = flipped_shift;
	}  else {
		std::stringstream msg;
		msg << "StrandPairings:: invalid orientation: " << nc_order_orient << std::endl;
		utility_exit_with_message( msg.str() );
	}

	return preliminary_shift;
}

HelixSheetPairing::HelixSheetPairing():
	SegmentPairing()
{}

HelixSheetPairing::HelixSheetPairing(
	SegmentName const & helix,
	SegmentName const & s1,
	SegmentName const & s2 ):
	SegmentPairing( boost::assign::list_of(helix)(s1)(s2).convert_to_container< SegmentNames >() )
{}

SegmentPairingOP
HelixSheetPairing::clone() const
{
	return SegmentPairingOP( new HelixSheetPairing( *this ) );
}

void
HelixSheetPairing::parse_tag( utility::tag::Tag const & )
{
}

void
HelixSheetPairing::to_xml( utility::tag::Tag & ) const
{
}

std::string
HelixSheetPairing::pairing_string( StructureData const & sd ) const
{
	using protocols::fldsgn::topology::SS_Info2;

	// FORMAT: H,S1-S2
	// H  : helix number
	// S1 : strand 1 number
	// S2 : strand 2 number

	SS_Info2 const ss_info( sd.ss() );

	// get helix id number
	SegmentNames::const_iterator segment_id = segments().begin();
	if ( segment_id == segments().end() ) {
		std::stringstream msg;
		msg << class_name() << ": no helix name was specified! Segments="
			<< segments() << std::endl << "SD = " << sd << std::endl;
		utility_exit_with_message( msg.str() );
	}
	core::Size const res_in_helix = sd.segment( *segment_id ).safe();
	core::Size const helix = ss_info.helix_id( res_in_helix );

	// now get strand1 id number
	++segment_id;
	if ( segment_id == segments().end() ) {
		std::stringstream msg;
		msg << class_name() << ": no strand 1 name was specified! Segments="
			<< segments() << std::endl << "SD = " << sd << std::endl;
		utility_exit_with_message( msg.str() );
	}
	core::Size const res_in_s1 = sd.segment( *segment_id ).safe();
	core::Size const s1 = ss_info.strand_id( res_in_s1 );

	// get strand2 id number
	++segment_id;
	if ( segment_id == segments().end() ) {
		std::stringstream msg;
		msg << class_name() << ": no strand 2 name was specified! Segments="
			<< segments() << std::endl << "SD = " << sd << std::endl;
		utility_exit_with_message( msg.str() );
	}
	core::Size const res_in_s2 = sd.segment( *segment_id ).safe();
	core::Size const s2 = ss_info.strand_id( res_in_s2 );

	std::stringstream pair_str;
	pair_str << helix << ',' << s1 << '-' << s2;
	return pair_str.str();
}

} //protocols
} //denovo_design
} //components
