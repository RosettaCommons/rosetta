// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/denovo_design/connection/ConnectionArchitect.cc
/// @brief Architect for covalently joining two segments of a pose
/// @author Tom Linsky (tlinsky@uw.edu)

// Unit headers
#include <protocols/denovo_design/connection/ConnectionArchitect.hh>

// Protocol headers
#include <protocols/denovo_design/components/IdealAbegoGenerator.hh>
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>
#include <protocols/denovo_design/components/StructureDataFactory.hh>
#include <protocols/denovo_design/util.hh>

// Core headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

// Basic/Utililty headers
#include <basic/Tracer.hh>
#include <utility/tag/Tag.hh>

// Boost headers
#include <boost/algorithm/string/predicate.hpp>
#include <boost/assign.hpp>
#include <numeric/random/random.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.denovo_design.connection.ConnectionArchitect" );

namespace protocols {
namespace denovo_design {
namespace connection {

ConnectionArchitect::ConnectionArchitect( std::string const & id_value ):
	protocols::denovo_design::architects::StructureArchitect( id_value ),
	ideal_abego_(),
	motifs_(),
	segment1_ids_(),
	segment2_ids_(),
	chain1_( 0 ),
	chain2_( 0 )
{
}

ConnectionArchitect::~ConnectionArchitect()
{}

ConnectionArchitectOP
ConnectionArchitect::clone() const
{
	return ConnectionArchitectOP( new ConnectionArchitect( *this ) );
}

std::string
ConnectionArchitect::type() const
{
	return ConnectionArchitect::class_name();
}

void
ConnectionArchitect::parse_tag( utility::tag::TagCOP tag, basic::datacache::DataMap & )
{
	set_user_chain1( tag->getOption< core::Size >( "chain1", user_chain1() ) );
	set_user_chain2( tag->getOption< core::Size >( "chain2", user_chain2() ) );

	std::string const seg1_ids = tag->getOption< std::string >( "segment1", "" );
	if ( !seg1_ids.empty() ) set_segment1_ids( seg1_ids );

	std::string const seg2_ids = tag->getOption< std::string >( "segment2", "" );
	if ( !seg2_ids.empty() ) set_segment2_ids( seg2_ids );

	std::string const motif_str = tag->getOption< std::string >( "motif", "" );
	std::string const cutpoint_str = tag->getOption< std::string >( "cutpoint", "" );
	if ( !motif_str.empty() ) set_motifs( motif_str, cutpoint_str );
}

void
ConnectionArchitect::apply( components::StructureData & sd ) const
{
	core::Real random = numeric::random::rg().uniform();
	apply( sd, random );
}

void
ConnectionArchitect::apply( components::StructureData & sd, core::Real & random ) const
{
	AreConnectablePredicate connectable( false );
	MotifOPs const connection_candidates = compute_connection_candidates( sd, connectable );
	if ( connection_candidates.empty() ) {
		std::stringstream ss;
		ss << id() << ": No valid free connection points could be found. User-specified segment1 ids: " << segment1_ids_
			<< " and segment2 ids: " << segment2_ids_ << std::endl;
		ss << sd << std::endl;
		throw EXCN_ConnectionSetupFailed( ss.str() );
	}

	MotifOP motif = choose_motif( connection_candidates, random );
	if ( motif ) {
		connect( sd, *motif );
	} else {
		std::stringstream ss;
		ss << id() << ": Failed to choose a motif to use for the connection. User-specified segment1 ids: " << segment1_ids_
			<< " and segment2 ids: " << segment2_ids_ << std::endl;
		ss << sd << std::endl;
		throw EXCN_ConnectionSetupFailed( ss.str() );
	}
}

/// @brief returns list of allowed segment 1 ids
SegmentNames const &
ConnectionArchitect::segment1_ids() const
{
	return segment1_ids_;
}

/// @brief returns list of allowed segment2 ids
SegmentNames const &
ConnectionArchitect::segment2_ids() const
{
	return segment2_ids_;
}

/// @brief sets list of segment1 ids from string
void
ConnectionArchitect::set_segment1_ids( std::string const & segment1_str )
{
	set_segment1_ids( parse_segment_names( segment1_str ) );
}

/// @brief sets list of segment1 ids from list
void
ConnectionArchitect::set_segment1_ids( SegmentNames const & segments )
{
	segment1_ids_ = segments;
}

/// @brief sets list of segment2 ids from string
void
ConnectionArchitect::set_segment2_ids( std::string const & segment2_str )
{
	set_segment2_ids( parse_segment_names( segment2_str ) );
}

/// @brief sets list of segment1 ids from list
void
ConnectionArchitect::set_segment2_ids( SegmentNames const & segments )
{
	segment2_ids_ = segments;
}

/// @brief gets user-specified chain number for the lower chain to be connected. 0 if not specified
core::Size
ConnectionArchitect::user_chain1() const
{
	return chain1_;
}

/// @brief sets user-specified chain number for the lower chain to be connected. 0 if unspecified
void
ConnectionArchitect::set_user_chain1( core::Size const chain )
{
	chain1_ = chain;
}

/// @brief gets user-specified chain number for the upper chain to be connected. 0 if not specified
core::Size
ConnectionArchitect::user_chain2() const
{
	return chain2_;
}

/// @brief sets user-specified chain number for the upper chain to be connected. 0 if unspecified
void
ConnectionArchitect::set_user_chain2( core::Size const chain )
{
	chain2_ = chain;
}

/// @brief sets motifs using a motif string and a string of cutpoints
void
ConnectionArchitect::set_motifs( std::string const & motif_str, std::string const & cutpoints_str )
{
	MotifCOPs const str_motifs = parse_motif_string( motif_str );
	Lengths const lengths = parse_length_string( cutpoints_str );

	if ( lengths.empty() ) {
		motifs_ = str_motifs;
		TR.Debug << "set_motifs(): Built " << str_motifs.size() << " motifs from "
			<< motif_str << std::endl;
		return;
	}

	MotifCOPs motifs;
	for ( MotifCOPs::const_iterator m_it=str_motifs.begin(); m_it!=str_motifs.end(); ++m_it ) {
		debug_assert( *m_it );
		MotifCOP m = *m_it;
		for ( Lengths::const_iterator l=lengths.begin(); l!=lengths.end(); ++l ) {
			if ( *l > m->elem_length() ) continue;
			MotifOP newmotif( new Motif( *m ) );
			newmotif->set_cutpoint( *l );
			motifs.push_back( newmotif );
		}
	}

	TR.Debug << "set_motifs(): Built " << motifs.size() << " motifs from "
		<< motif_str << " and " << cutpoints_str << std::endl;

	motifs_ = motifs;
}

/// @brief sets motifs via a vector
void
ConnectionArchitect::set_motifs( MotifCOPs const & motifs )
{
	motifs_ = motifs;
}

MotifCOPs
ConnectionArchitect::parse_motif_string( std::string const & motif_str ) const
{
	MotifCOPs motifs;
	utility::vector1< std::string > const motif_strs = utility::string_split( motif_str, ',' );
	for ( utility::vector1< std::string >::const_iterator mstr=motif_strs.begin(); mstr!=motif_strs.end(); ++mstr ) {
		MotifOP newmotif( new Motif() );
		newmotif->parse_motif( *mstr );
		motifs.push_back( newmotif );
	}
	return motifs;
}

MotifOPs
ConnectionArchitect::compute_connection_candidates(
	components::StructureData const & sd,
	AreConnectablePredicate const & connectable ) const
{
	SegmentPairs const seg_pairs = segment_pairs( sd );

	debug_assert( seg_pairs.size() > 0 );

	// set of valid lengths for use when trying to compute idealized abegos
	LengthSet const length_set = lengths();

	MotifOPs candidates;
	for ( SegmentPairs::const_iterator pair=seg_pairs.begin(); pair!=seg_pairs.end(); ++pair ) {
		//if ( ! pair_allowed( pair->first, pair->second ) ) {
		//  TR.Debug << "ConnectionArchitect: Connection to segments " << pair->first
		//   << ", " << pair->second << " DISALLOWED by user setting." << std::endl;
		//  continue;
		//}

		MotifOPs const motifs = motifs_for_pair( *pair, sd, length_set );

		// if nothing is specified, try using an empty motif
		debug_assert( !motifs.empty() );
		for ( MotifOPs::const_iterator m=motifs.begin(); m!=motifs.end(); ++m ) {
			if ( connectable( sd, **m ) ) {
				TR.Debug << "c1, c2, len : connectable " << pair->first << " <--> "
					<< pair->second << ", " << **m << std::endl;
				candidates.push_back( *m );
			}
		}
	}

	return candidates;
}

ConnectionArchitect::SegmentPairs
ConnectionArchitect::segment_pairs( components::StructureData const & sd ) const
{
	// if a segment exists that was built by a connection of the same name, maintain the same termini as before
	SegmentNames local_comp1_ids;
	SegmentNames local_comp2_ids;

	// determine alternate id -- strip off all parent names
	std::string alt_id = id();
	for ( SegmentNameList::const_iterator c=sd.segments_begin(); c!=sd.segments_end(); ++c ) {
		if ( boost::ends_with( *c, id() ) ) {
			alt_id = *c;
		}
	}

	SegmentNames const available_uppers = available_upper_termini( sd );
	if ( chain1_ ) {
		debug_assert( chain1_ <= available_uppers.size() );
		local_comp1_ids.push_back( available_uppers[ chain1_ ] );
	} else {
		local_comp1_ids = available_uppers;
	}

	SegmentNames const available_lowers = available_lower_termini( sd );
	if ( chain2_ ) {
		debug_assert( chain2_ <= available_lowers.size() );
		local_comp2_ids.push_back( available_lowers[ chain2_ ] );
	} else {
		local_comp2_ids = available_lowers;
	}

	if ( local_comp1_ids.empty() ) {
		std::stringstream err;
		err << "Connection " << id() << ": " << " no available segment1 upper termini were found matching the user's input.";
		err << "Input ids: " << segment1_ids_ << " User chain: " << chain1_ << " Perm: " << std::endl;
		err << sd << std::endl;
		throw utility::excn::EXCN_Msg_Exception( err.str() );
	}

	if ( local_comp2_ids.empty() ) {
		std::stringstream err;
		err << "Connection " << id() << ": " << " no available segment2 lower termini were found matching the user's input.";
		err << "Input ids: " << segment2_ids_ << " User chain: " << chain2_ << " Perm: " << std::endl;
		err << sd << std::endl;
		throw utility::excn::EXCN_Msg_Exception( err.str() );
	}

	return combine_segment_names( local_comp1_ids, local_comp2_ids );
	/*
	bool const repeat = ( sd.has_segment( alt_id ) ||
	( sd.has_data_int( id(), "length" ) && ( sd.get_data_int( id(), "length" ) == 0 ) ) );

	if ( repeat ) {
	// if we are re-calling a connection, assume we want to rebuild it
	// copy conection points and delete the old loop
	local_comp1_ids.push_back( lower_segment_id( perm ) );
	local_comp2_ids.push_back( upper_segment_id( perm ) );
	for ( SegmentNameList::const_iterator l=segment_names().begin(); l!=segment_names().end(); ++l ) {
	if ( perm.has_segment( *l ) ) {
	perm.delete_segment( *l );
	} else {
	perm.disconnect_segments( lower_segment_id( perm ), upper_segment_id( perm ) );
	}
	}
	} else {
	if ( chain1_ ) {
	utility::vector1< std::string > const available_uppers = find_available_upper_termini( perm );
	debug_assert( chain1_ <= available_uppers.size() );
	local_comp1_ids.push_back( available_uppers[ chain1_ ] );
	} else {
	local_comp1_ids = find_available_upper_termini( perm );
	}
	if ( chain2_ ) {
	utility::vector1< std::string > const available_lowers = find_available_lower_termini( perm );
	debug_assert( chain2_ <= available_lowers.size() );
	local_comp2_ids.push_back( available_lowers[ chain2_ ] );
	} else {
	local_comp2_ids = find_available_lower_termini( perm );
	}
	} */
}

ConnectionArchitect::SegmentPairs
ConnectionArchitect::combine_segment_names(
	SegmentNames const & seg1s,
	SegmentNames const & seg2s ) const
{
	SegmentPairs pairs;

	for ( SegmentNames::const_iterator s1=seg1s.begin(); s1!=seg1s.end(); ++s1 ) {
		for ( SegmentNames::const_iterator s2=seg2s.begin(); s2!=seg2s.end(); ++s2 ) {
			pairs.push_back( std::make_pair( *s1, *s2 ) );
		}
	}
	return pairs;
}

SegmentNames
ConnectionArchitect::available_upper_termini( components::StructureData const & sd ) const
{
	SegmentNames local_ids;
	if ( segment1_ids_.empty() ) {
		local_ids = sd.available_upper_termini();
	} else {
		for ( SegmentNames::const_iterator s=segment1_ids_.begin(); s!=segment1_ids_.end(); ++s ) {
			if ( sd.has_free_upper_terminus( *s ) ) {
				local_ids.push_back( *s );
			}
		}
	}

	TR.Debug << "Found " << local_ids.size() << " available upper termini: "
		<< local_ids << std::endl;
	return local_ids;
}

/// @brief finds usable/available upper termini (i.e. those for comp1)
SegmentNames
ConnectionArchitect::available_lower_termini( components::StructureData const & sd ) const
{
	SegmentNames local_ids;
	if ( segment2_ids_.empty() ) {
		local_ids = sd.available_lower_termini();
	} else {
		for ( SegmentNames::const_iterator s=segment2_ids_.begin(); s!=segment2_ids_.end(); ++s ) {
			if ( sd.has_free_lower_terminus( *s ) ) {
				local_ids.push_back( *s );
			}
		}
	}

	TR.Debug << "Found " << local_ids.size() << " available lower termini: "
		<< local_ids << std::endl;
	return local_ids;
}

MotifOP
ConnectionArchitect::choose_motif( MotifOPs const & motifs, core::Real & random ) const
{
	// now we just choose a valid combo from the list
	core::Size const info_idx = extract_int( random, 1, motifs.size() );
	return motifs[ info_idx ];
}

ConnectionArchitect::LengthSet
ConnectionArchitect::lengths() const
{
	LengthSet length_set;
	for ( MotifCOPs::const_iterator m=motifs_.begin(); m!=motifs_.end(); ++m ) {
		debug_assert( *m );
		length_set.insert( (*m)->length() );
	}
	return length_set;
}

MotifOPs
ConnectionArchitect::motifs_for_pair(
	SegmentPair const & pair,
	components::StructureData const & sd,
	LengthSet const & length_set ) const
{
	MotifOPs all_motifs;
	if ( ideal_abego_ ) {
		core::Size const res1 = sd.segment( pair.first ).stop();
		core::Size const res2 = sd.segment( pair.second ).start();
		all_motifs = ideal_abego_->generate( sd.abego( res1 ), sd.abego( res2 ), length_set );
	} else {
		for ( MotifCOPs::const_iterator m=motifs_.begin(); m!=motifs_.end(); ++m ) {
			all_motifs.push_back( (*m)->clone() );
		}
	}

	// if list is empty, generate a 0-length motif
	if ( all_motifs.empty() ) {
		all_motifs.push_back( MotifOP( new Motif() ) );
	}

	// attach lower and upper terminal segments
	for ( MotifOPs::const_iterator m=all_motifs.begin(); m!=all_motifs.end(); ++m ) {
		(*m)->set_lower_segment( pair.first );
		(*m)->set_upper_segment( pair.second );
	}

	return all_motifs;
}

/// @brief creates StructureData from given ConnectionInfo object
void
ConnectionArchitect::connect( components::StructureData & sd, Motif & motif ) const
{
	// add segment and remove upper/lower segment names that were added
	// when enumerating possibilities
	// Segment names can be retrieved from the string cache in StructureData
	std::string const lower = motif.lower_segment();
	std::string const upper = motif.upper_segment();
	sd.set_data_str( id(), "segment1", lower );
	sd.set_data_str( id(), "segment2", upper );

	// special case when motif length == 0
	if ( motif.length() == 0 ) {
		sd.move_segment( lower, upper );
		sd.delete_trailing_residues( lower );
		sd.delete_leading_residues( upper );
		sd.connect_segments( lower, upper );
		return;
	}

	motif.set_lower_segment( "" );
	motif.set_upper_segment( "" );
	sd.add_segment( id(), motif );

	sd.move_segment( lower, id() );
	sd.delete_trailing_residues( lower );
	sd.delete_leading_residues( id() );
	sd.connect_segments( lower, id() );

	sd.move_segment( id(), upper );
	sd.delete_trailing_residues( id() );
	sd.delete_leading_residues( upper );
	sd.connect_segments( id(), upper );

	// loop can move relative to other segments, so put it in its own movable group
	sd.set_movable_group( id(), sd.choose_new_movable_group() );
}

SegmentNames
parse_segment_names( std::string const & segment_name_str )
{
	if ( segment_name_str.empty() ) {
		return SegmentNames();
	}
	utility::vector1< std::string > const names = utility::string_split( segment_name_str, ',' );
	return SegmentNames( names.begin(), names.end() );
}

core::Real
calc_approx_loop_length( std::string const & abego )
{
	// these values are based on SIN(ANGLECHANGE/2)*3.8
	static std::map< std::string, core::Real > const distmap =
		boost::assign::map_list_of
		("AA",2.18)
		("AB",1.90)
		("AE",3.79)
		("AG",3.67)
		("BA",3.44)
		("BB",3.57)
		("BE",0.98)
		("BG",0.33)
		("EA",0.33)
		("EB",0.66)
		("EE",2.91)
		("EG",3.57)
		("GA",3.74)
		("GB",3.67)
		("GE",2.91)
		("GG",1.90);
	core::Real max_dist = 0.0;
	for ( core::Size i=1, endi=abego.size(); i<endi; ++i ) {
		std::string const dyad = abego.substr(i-1,2);
		std::map< std::string, core::Real >::const_iterator d = distmap.find(dyad);
		if ( d == distmap.end() ) {
			max_dist += 3.8;
		} else {
			max_dist += d->second;
		}
	}
	return max_dist;
}

AreConnectablePredicate::AreConnectablePredicate( bool const allow_cyclic ):
	allow_cyclic_( allow_cyclic )
{}

/// @brief checks whether two segments can be connected
bool
AreConnectablePredicate::operator()(
	components::StructureData const & sd,
	Motif const & motif ) const
{
	std::string const & segment1 = motif.lower_segment();
	std::string const & segment2 = motif.upper_segment();
	debug_assert( !segment1.empty() );
	debug_assert( !segment2.empty() );

	// if one id or the other is already connected to something, this isn't connectable
	if ( ! sd.segment(segment1).has_free_upper_terminus() ) {
		TR.Debug << segment1 << " and " << segment2 << " not connectable due to no free c anchor." << std::endl;
		TR.Debug << " available upper termini= " << sd.available_upper_termini() << std::endl;
		return false;
	}

	if ( ! sd.segment(segment2).has_free_lower_terminus() ) {
		TR.Debug << segment1 << " and " << segment2 << " not connectable due to no free n anchor." << std::endl;
		TR.Debug << " available lower termini= " << sd.available_lower_termini() << std::endl;
		return false;
	}

	// ensure the two segments are on different chains
	// connecting a segment to itself is illegal unless allow_cyclic is true
	std::pair< std::string, std::string > const termini1 = sd.termini( segment1 );
	std::pair< std::string, std::string > const termini2 = sd.termini( segment2 );
	if ( !allow_cyclic_ && ( termini1 == termini2 ) ) {
		TR.Debug << segment1 << " and " << segment2 << " not connectable because it would create a cyclic peptide." << std::endl;
		return false;
	}

	// if there is a cutpoint, the two segments should be in the same group
	// if not, the two segments should be in different movable groups
	if ( ! check_movable_groups( sd, motif ) ) {
		return false;
	}

	return check_distance( sd, motif );
}

AreConnectablePredicate::MovableGroupSet
AreConnectablePredicate::connected_movable_groups( components::StructureData const & sd, std::string const & seg_name ) const
{
	SegmentNameList const & segments = sd.connected_segments( seg_name, true );
	MovableGroupSet mgs;
	for ( SegmentNameList::const_iterator s=segments.begin(); s!=segments.end(); ++s ) {
		mgs.insert( sd.segment( *s ).movable_group() );
	}
	return mgs;
}

bool
AreConnectablePredicate::check_movable_groups(
	components::StructureData const & sd,
	Motif const & motif ) const
{
	MovableGroupSet const mgs1 = connected_movable_groups( sd, motif.lower_segment() );
	MovableGroupSet const mgs2 = connected_movable_groups( sd, motif.upper_segment() );

	// no duplicate movable groups should be present
	MovableGroupSet intersection;
	std::set_intersection( mgs1.begin(), mgs1.end(), mgs2.begin(), mgs2.end(),
		std::inserter( intersection, intersection.begin() ) );

	if ( ! motif.cutpoint() ) {
		if ( !intersection.empty() ) {
			TR.Debug << motif.lower_segment() << " and " << motif.upper_segment()
				<< " are not connectable because there is not cutpoint, and both are connected to "
				<< "segments in the same movable group. Movable groups connected to both segments: "
				<< intersection << std::endl;
			return false;
		}
		return true;
	} else {
		if ( intersection.empty() ) {
			TR.Debug << motif.lower_segment() << " and " << motif.upper_segment()
				<< " are not connectable because there is a cutpoint, but neither is connected to "
				<< "any segments with common movable groups. Motif=" << motif << " SD=" << sd << std::endl;
			return false;
		}
		return true;
	}
}

bool
AreConnectablePredicate::check_distance(
	components::StructureData const & sd,
	Motif const & motif ) const
{
	// max 3.8 angstroms per residue, plus ~1.5 angstroms for N-C bond
	static core::Real const bond_dist = 1.5;
	static core::Real const max_dist_per_res = 3.8;

	if ( !motif.cutpoint() ) return true;

	std::string const & segment1 = motif.lower_segment();
	std::string const & segment2 = motif.upper_segment();

	core::pose::PoseCOP template1 = sd.segment( segment1 ).template_pose();
	if ( !template1 ) {
		TR.Debug << "skipping distance calculation because " << segment1 << " has no template pose." << std::endl;
		return true;
	}

	core::pose::PoseCOP template2 = sd.segment( segment2 ).template_pose();
	if ( !template2 ) {
		TR.Debug << "skipping distance calculation because " << segment2 << " has no template pose." << std::endl;
		return true;
	}

	core::Real avg_dist = calc_approx_loop_length(
		sd.abego( sd.segment( segment1 ).stop() ) +
		motif.abego() +
		sd.abego( sd.segment( segment2 ).start() ) );
	avg_dist /= static_cast< core::Real >( motif.length() );

	// get residue and atom from segment 1
	core::Size const res1 = template1->total_residue();
	core::chemical::ResidueType const & rtype = template1->residue( res1 ).residue_type_set()->name_map( template1->residue(res1).name3() );
	std::string const & aname = rtype.atom_name( rtype.upper_connect_atom() );
	core::Vector const xyz1 = template1->residue(res1).xyz( aname );

	// get residue and atom from segment 2
	core::Size const res2 = 1;
	core::chemical::ResidueType const & rtype2 = template2->residue( res2 ).residue_type_set()->name_map( template2->residue(res2).name3() );
	std::string const & aname2 = rtype2.atom_name( rtype2.lower_connect_atom() );
	core::Vector const xyz2 = template2->residue(res2).xyz( aname2 );

	// compute distance
	core::Real const dist = xyz1.distance( xyz2 );

	// if the distance is > the fully extended Ca-Ca distance of an nres residue insert, this connection is physically impossible
	if ( dist > (max_dist_per_res*motif.elem_length() + bond_dist) ) {
		TR.Debug << segment1 << " and " << segment2 << " with motif " << motif
			<< "\tnot connectable due to nres="
			<< motif.length() << " distance=" << dist << std::endl;
		return false;
	}

	return true;
}

} //protocols
} //denovo_design
} //connection
