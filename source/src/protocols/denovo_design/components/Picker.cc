// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/denovo_design/components/Picker.cc
/// @brief Auto-caching fragment picker
/// @details
/// @author Tom Linsky


//Unit Headers
#include <protocols/denovo_design/components/Picker.hh>

//Project Headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>

//Protocol Headers
#include <protocols/denovo_design/util.hh>
#include <protocols/forge/build/Interval.hh>
#include <protocols/loops/Loop.hh>
#include <protocols/loops/Loops.hh>

//Core Headers
#include <core/conformation/Conformation.hh>
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/sequence/ABEGOManager.hh>

//Basic/Utility/Numeric Headers
#include <basic/Tracer.hh>

#include <utility/excn/Exceptions.hh>

// Boost/ObjexxFCL Headers
#include <boost/lexical_cast.hpp>

//C++ Headers

static THREAD_LOCAL basic::Tracer TR("protocols.denovo_design.components.Picker");

////////////////////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace denovo_design {
namespace components {
////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Autocaching fragment picker
//////////////////////////////////////////////////////////////////////////////////////////////////////

Picker::Picker() :
	utility::SingletonBase< Picker >(),
	vlb_(),
	fragcache_(),
	n_frags_( 250 )
{
}

Picker::~Picker() {}

Picker *
Picker::create_singleton_instance()
{
	return new Picker();
}

void
Picker::set_nfrags( core::Size const nfrags )
{
	n_frags_ = nfrags;
}

/// @brief Simply picks fragments for the given position, caching as needed
/// @param[in] chain_aa The complete amino acid string, typically from a Pose;
///                     can be empty.  If empty, sequence bias is not used to pick fragments.
/// @param[in] chain_ss The complete secondary structure string, typically from a Pose.
/// @param[in] chain_abego The complete abego string, typically from a setter, set_abego
/// @param[in] resid The starting residue to pick fragments from; Pose
///  numbering (i.e. 1-based indexing).
/// @param[in] frag_length The desired length of the fragments
/// @param[in] n_frags The number of fragments to pick per position.
core::fragment::FrameOP
Picker::pick_fragments_for_resid(
	std::string const & chain_aa,
	std::string const & chain_ss,
	utility::vector1< std::string > const & chain_abego,
	core::Size const resid,
	core::Size const frag_length )
{
	// Get fragment key
	std::string const key = ss_key( chain_aa, chain_ss, chain_abego, resid, frag_length );
	TR.Debug << "Fragment SS Key=" << key << std::endl;

	FragmentMap::iterator cached = fragcache_.find( key );
	if ( cached == fragcache_.end() ) {
		core::fragment::FrameList framelist = vlb_.pick_fragments_public(
			chain_ss,                                          // secondary structure of entire pose
			chain_aa,                                          // amino acid sequence
			chain_abego,                                       // abegos for entire pose
			protocols::forge::build::Interval( resid, resid ), // start/end of region for picking
			frag_length,                                       // fragment length
			n_frags_ );                                        // # fragments per position
		core::fragment::FrameOP frame = *framelist.begin();
		cached = fragcache_.insert( std::make_pair( key, frame ) ).first;
		return frame->clone_with_frags();
	}

	debug_assert( cached != fragcache_.end() );
	debug_assert( cached->second );
	core::fragment::FrameOP frame = cached->second->clone_with_frags();
	debug_assert( frame );
	frame->shift_to( resid );
	TR.Debug << "Retrieved " << key << " from cache" << std::endl;
	return frame;
}

/// @brief pick fragments of a given length, padding when necessary -- CACHES fragments for all positions
/// @param[in] complete_aa The complete amino acid string, typically from a Pose;
///            can be empty.  If empty, sequence bias is not used to pick fragments.
/// @param[in] complete_ss The complete secondary structure string, typically from a Pose.
/// @param[in] complete_abego The complete abego string, typically from a setter, set_abego
/// @param[in] interval The interval [left, right] to pick fragments from; Pose
///  numbering (i.e. 1-based indexing).
/// @param[in] frag_length The desired length of the fragments
/// @param[in] n_frags The number of fragments to pick per position.
core::fragment::FrameList
Picker::get_framelist(
	std::string const & complete_aa,
	std::string const & complete_ss,
	utility::vector1< std::string > const & complete_abego,
	utility::vector1< core::Size > const & chain_endings,
	core::Size const start_res,
	core::Size const end_res,
	core::Size const frag_length )
{
	assert( start_res >= 1 );
	assert( end_res >= 1 );
	assert( start_res <= complete_ss.size() );
	assert( end_res <= complete_ss.size() );
	assert( (!complete_abego.size()) || ( start_res <= complete_abego.size() ) );
	assert( (!complete_abego.size()) || ( end_res <= complete_abego.size() ) );

	core::Size const closest_chain_ending = get_closest_chain_ending( chain_endings, complete_ss.size(), end_res );

	// cut ss and abego to stay in current chain
	TR.Debug << "Closest chain ending is " << closest_chain_ending << std::endl;
	std::string const chain_ss = complete_ss.substr( 0, closest_chain_ending );
	utility::vector1< std::string > const chain_abego = truncate_abego( complete_abego, closest_chain_ending );
	std::string const chain_aa = complete_aa.substr( 0, closest_chain_ending );

	core::fragment::FrameList framelist;
	for ( core::Size resid=start_res; resid<=end_res; ++resid ) {
		framelist.push_back( pick_fragments_for_resid( chain_aa, chain_ss, chain_abego, resid, frag_length ) );
	}

	TR << "Obtained fragments at " << framelist.size() << " positions from " << start_res
		<< " to " << end_res << "." << std::endl;
	return framelist;
}

utility::vector1< std::string >
truncate_abego( utility::vector1< std::string > const & complete_abego, core::Size const closest_chain_ending )
{
	utility::vector1< std::string > chain_abego;
	core::Size cur_resid = 1;
	for ( utility::vector1< std::string >::const_iterator ab=complete_abego.begin(); ab!=complete_abego.end(); ++ab, ++cur_resid ) {
		if ( cur_resid > closest_chain_ending ) break;
		chain_abego.push_back( *ab );
	}

	return chain_abego;
}

core::Size
get_closest_chain_ending( utility::vector1< core::Size > const & chain_endings, core::Size total_residue, core::Size const end_res )
{
	for ( utility::vector1< core::Size >::const_iterator r=chain_endings.begin(); r!=chain_endings.end(); ++r ) {
		if ( ( *r > end_res ) && ( *r < total_residue ) ) {
			total_residue = *r;
		}
	}
	return total_residue;
}

core::fragment::ConstantLengthFragSetOP
Picker::pick_and_cache_fragments(
	std::string const & complete_ss,
	utility::vector1< std::string > const & complete_abego,
	utility::vector1< core::Size > const & chain_endings,
	core::Size const start_res,
	core::Size const end_res,
	core::Size const frag_length )
{
	core::fragment::ConstantLengthFragSetOP fragset =
		core::fragment::ConstantLengthFragSetOP( new core::fragment::ConstantLengthFragSet( frag_length ) );

	fragset->add( get_framelist( "", complete_ss, complete_abego, chain_endings, start_res, end_res, frag_length ) );
	return fragset;
}

core::fragment::ConstantLengthFragSetOP
Picker::pick_and_cache_fragments(
	std::string const & complete_aa,
	std::string const & complete_ss,
	utility::vector1< std::string > const & complete_abego,
	utility::vector1< core::Size > const & chain_endings,
	core::Size const start_res,
	core::Size const end_res,
	core::Size const frag_length )
{
	core::fragment::ConstantLengthFragSetOP fragset =
		core::fragment::ConstantLengthFragSetOP( new core::fragment::ConstantLengthFragSet( frag_length ) );

	fragset->add( get_framelist( complete_aa, complete_ss, complete_abego, chain_endings, start_res, end_res, frag_length ) );
	return fragset;
}

/// @brief picks and caches fragments for the listed components with size frag_size
/// this does cache and should be the primary function called by objects needing fragments
core::fragment::ConstantLengthFragSetOP
Picker::pick_and_cache_fragments(
	std::string const & ss,
	std::string const & abego,
	protocols::loops::Loops const & loops,
	utility::vector1< core::Size > const & chain_endings,
	core::Size const frag_length )
{
	utility::vector1< std::string > const abegovec = abego_vector( abego );

	core::fragment::ConstantLengthFragSetOP fragset =
		core::fragment::ConstantLengthFragSetOP( new core::fragment::ConstantLengthFragSet( frag_length ) );

	for ( protocols::loops::Loops::const_iterator l=loops.begin(); l!=loops.end(); ++l ) {
		fragset->add( get_framelist( "", ss, abegovec, chain_endings, l->start(), l->stop(), frag_length ) );
	}
	return fragset;
}

/// @brief picks and caches fragments for the listed components with size frag_size
/// this does cache and should be the primary function called by objects needing fragments
core::fragment::ConstantLengthFragSetOP
Picker::fragments_for_permutation(
	StructureData const & perm,
	SegmentNameList const & comp_ids,
	core::Size const frag_length )
{
	core::fragment::ConstantLengthFragSetOP fragset =
		core::fragment::ConstantLengthFragSetOP( new core::fragment::ConstantLengthFragSet( frag_length ) );

	utility::vector1< core::Size > chain_endings;
	for ( SegmentNameList::const_iterator s=perm.segments_begin(); s!=perm.segments_end(); ++s ) {
		if ( perm.segment( *s ).has_free_upper_terminus() ) {
			chain_endings.push_back( perm.segment( *s ).upper() );
		}
	}

	utility::vector1< std::string > const abegovec = abego_vector( perm.abego() );
	for ( SegmentNameList::const_iterator s=comp_ids.begin(); s!=comp_ids.end(); ++s ) {
		fragset->add( get_framelist( "", perm.ss(), abegovec, chain_endings, perm.segment(*s).lower(), perm.segment(*s).upper(), frag_length ) );
	}

	return fragset;
}

/// @brief picks and caches fragments for the listed components with size frag_size
/// this does cache and should be the primary function called by objects needing fragments
/// replaces X in the desired abego with whatever is in the pose
core::fragment::ConstantLengthFragSetOP
Picker::fragments_for_permutation_take_X_from_pose(
	StructureData const & perm,
	core::pose::Pose const & pose,
	SegmentNames const & comp_ids,
	core::Size const frag_length )
{
	core::fragment::ConstantLengthFragSetOP fragset =
		core::fragment::ConstantLengthFragSetOP( new core::fragment::ConstantLengthFragSet( frag_length ) );

	//utility::vector1< std::string > complete_abego = core::sequence::ABEGOManager().get_symbols( *(perm.pose()), 1 );

	for ( SegmentNames::const_iterator c=comp_ids.begin(); c!=comp_ids.end(); ++c ) {
		Segment residues = perm.segment( *c );
		assert( residues.stop() >= residues.start() );
		assert( residues.start() >= 1 );
		utility::vector1< std::string > poseabego = core::sequence::ABEGOManager().get_symbols( pose, 1 );
		utility::vector1< std::string > newabego = abego_vector( perm.abego() );
		for ( core::Size i=1; i<=newabego.size(); ++i ) {
			if ( newabego[i] == "X" ) {
				newabego[i] = poseabego[i];
			}
		}
		fragset->add( get_framelist( "", perm.ss(), newabego, pose.conformation().chain_endings(), residues.lower(), residues.upper(), frag_length ) );
	}
	return fragset;
}

/// @brief generates a key based on secondary structure to be used in fragcache
std::string
Picker::ss_key(
	std::string const & aa,
	std::string const & ss,
	utility::vector1< std::string > const & abego,
	core::Size const resid,
	core::Size const fragsize ) const
{
	debug_assert( resid >= 1 );
	debug_assert( resid <= ss.size() );
	debug_assert( (!abego.size()) || ( resid <= abego.size() ) );

	core::Size const stop = resid + fragsize;
	TR.Debug << "Generating key from resid " << resid << " ss=" << ss << std::endl;

	std::stringstream key;
	key << fragsize;
	for ( core::Size n=resid; n<stop; ++n ) {
		if ( n >= ss.size() ) key << 'D';
		else key << ss[ n - 1 ];
	}
	for ( core::Size n=resid; n<stop; ++n ) {
		if ( n > abego.size() ) key << 'X';
		else key << abego[n];
	}
	if ( !aa.empty() ) {
		for ( core::Size n=resid; n<stop; ++n ) {
			if ( n > aa.size() ) key << 'Z';
			else key << aa[ n - 1 ];
		}
	}

	return key.str();
}

} // components
} // denovo_design
} // protocols

// singleton stuff
namespace utility {

using protocols::denovo_design::components::Picker;

#if defined MULTI_THREADED && defined CXX11
template<> std::mutex utility::SingletonBase< Picker >::singleton_mutex_{};
template<> std::atomic< Picker * > utility::SingletonBase< Picker >::instance_( NULL );
#else
template<> Picker * utility::SingletonBase< Picker >::instance_( NULL );
#endif

}
