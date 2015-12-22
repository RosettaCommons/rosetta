// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/denovo_design/components/Picker.cc
/// @brief Auto-caching fragment picker
/// @detailed
/// @author Tom Linsky


//Unit Headers
#include <protocols/denovo_design/components/Picker.hh>

//Project Headers
#include <protocols/denovo_design/components/Segment.hh>
#include <protocols/denovo_design/components/StructureData.hh>

//Protocol Headers
#include <protocols/denovo_design/util.hh>
#include <protocols/forge/build/Interval.hh>

//Core Headers
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
	protocols::forge::components::VarLengthBuild(),
	n_frags_( 200 )
{
	fragcache_.clear();
}

Picker::Picker( core::Size const n_frags ) :
	protocols::forge::components::VarLengthBuild(),
	n_frags_( n_frags )
{
	fragcache_.clear();
}

Picker::~Picker() {}

/// @brief pick fragments of a given length, padding when necessary -- DOES NOT CACHE
/// @param[in] complete_ss The complete secondary structure string, typically from a Pose.
/// @param[in] complete_aa The complete amino acid string, typically from a Pose;
///            can be empty.  If empty, sequence bias is not used to pick fragments.
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

	// find ss key
	std::string const key = ss_key( complete_aa, complete_ss, complete_abego, start_res, end_res, frag_length );
	TR.Debug << "Fragment SS Key=" << key << std::endl;

	// look for fragments in cache
	std::map< std::string, core::fragment::FrameList >::iterator it = fragcache_.find(key);
	if ( it == fragcache_.end() ) {
		// not found
		core::fragment::FrameList framelist = pick_fragments(
			complete_ss,                                               // secondary structure of entire pose
			complete_aa,                                               // amino acid sequence
			complete_abego,                                            // abegos for entire pose
			protocols::forge::build::Interval( start_res, end_res ),   // start/end of region for picking
			frag_length,                                               // fragment length
			n_frags_ );                                                // # fragments per position
		fragcache_[key] = framelist;
		TR << "Saved " << key << " in cache" << std::endl;
		assert( framelist.size() );
		return framelist;
	} else {
		// fragments exist
		core::Size res = start_res;
		core::Size prevframe = 0;
		core::fragment::FrameList newlist;
		for ( core::fragment::FrameList::const_iterator f = it->second.begin(); f != it->second.end(); ++f ) {
			debug_assert( *f );
			if ( prevframe && ( prevframe + 1 != (*f)->start() ) ) {
				std::stringstream msg;
				msg << "Picker: frames are not in order. Prev=" << prevframe << " Cur=" << (*f)->start() << std::endl;
				throw utility::excn::EXCN_Msg_Exception( msg.str() );
			}
			core::fragment::FrameOP frame = (*f)->clone_with_frags();
			debug_assert( frame );
			frame->shift_to(res);
			newlist.push_back( frame );
			TR.Debug << "Added frame with nr_frags=" << frame->nr_frags() << std::endl;
			++res;
			prevframe = (*f)->start();
		}
		TR << "Retreiving " << key << " from cache" << std::endl;
		assert( newlist.size() );
		return newlist;
	}
}

core::fragment::ConstantLengthFragSetOP
Picker::pick_and_cache_fragments(
	std::string const & complete_ss,
	utility::vector1< std::string > const & complete_abego,
	core::Size const start_res,
	core::Size const end_res,
	core::Size const frag_length )
{
	core::fragment::ConstantLengthFragSetOP fragset =
		core::fragment::ConstantLengthFragSetOP( new core::fragment::ConstantLengthFragSet( frag_length ) );

	fragset->add( get_framelist( "", complete_ss, complete_abego, start_res, end_res, frag_length ) );
	return fragset;
}

core::fragment::ConstantLengthFragSetOP
Picker::pick_and_cache_fragments(
	std::string const & complete_aa,
	std::string const & complete_ss,
	utility::vector1< std::string > const & complete_abego,
	core::Size const start_res,
	core::Size const end_res,
	core::Size const frag_length )
{
	core::fragment::ConstantLengthFragSetOP fragset =
		core::fragment::ConstantLengthFragSetOP( new core::fragment::ConstantLengthFragSet( frag_length ) );

	fragset->add( get_framelist( complete_aa, complete_ss, complete_abego, start_res, end_res, frag_length ) );
	return fragset;
}

/// @brief picks and caches fragments for the listed components with size frag_size
/// this does cache and should be the primary function called by objects needing fragments
core::fragment::ConstantLengthFragSetOP
Picker::fragments_for_permutation(
	StructureData const & perm,
	StringList const & comp_ids,
	core::Size const frag_length )
{
	core::fragment::ConstantLengthFragSetOP fragset =
		core::fragment::ConstantLengthFragSetOP( new core::fragment::ConstantLengthFragSet( frag_length ) );

	//utility::vector1< std::string > complete_abego = core::sequence::ABEGOManager().get_symbols( *(perm.pose()), 1 );

	for ( StringList::const_iterator s = comp_ids.begin(); s != comp_ids.end(); ++s ) {
		Segment const & residues = perm.segment( *s );
		debug_assert( residues.stop() >= residues.start() );
		debug_assert( residues.start() >= 1 );

		std::string ss = perm.ss();
		utility::vector1< std::string > abego = perm.abego();
		if ( residues.upper_segment().empty() ) {
			ss = ss.substr( 0, residues.cterm_resi() );
			abego = utility::vector1< std::string >( abego.begin(), abego.begin() + residues.cterm_resi() );
		}
		fragset->add( get_framelist( "", ss, abego, residues.nterm_resi(), residues.cterm_resi(), frag_length ) );
	}
	return fragset;
}

/// @brief picks and caches fragments for the listed components with size frag_size
/// this does cache and should be the primary function called by objects needing fragments
/// replaces X in the desired abego with whatever is in the pose
core::fragment::ConstantLengthFragSetOP
Picker::fragments_for_permutation_take_X_from_pose(
	StructureData const & perm,
	utility::vector1< std::string > const & comp_ids,
	core::Size const frag_length )
{
	core::fragment::ConstantLengthFragSetOP fragset =
		core::fragment::ConstantLengthFragSetOP( new core::fragment::ConstantLengthFragSet( frag_length ) );

	//utility::vector1< std::string > complete_abego = core::sequence::ABEGOManager().get_symbols( *(perm.pose()), 1 );

	for ( core::Size i=1; i<=comp_ids.size(); ++i ) {
		Segment residues = perm.segment( comp_ids[i] );
		assert( residues.stop() >= residues.start() );
		assert( residues.start() >= 1 );
		utility::vector1< std::string > poseabego = core::sequence::ABEGOManager().get_symbols( *(perm.pose()), 1 );
		utility::vector1< std::string > newabego = perm.abego();
		for ( core::Size i=1; i<=newabego.size(); ++i ) {
			if ( newabego[i] == "X" ) {
				newabego[i] = poseabego[i];
			}
		}
		fragset->add( get_framelist( "", perm.ss(), newabego, residues.nterm_resi(), residues.cterm_resi(), frag_length ) );
	}
	return fragset;
}

/// @brief generates a key based on secondary structure to be used in fragcache
std::string
Picker::ss_key(
	std::string const & aa,
	std::string const & ss,
	utility::vector1< std::string > const & abego,
	core::Size const start,
	core::Size const end,
	core::Size const fragsize ) const
{
	assert( start <= end );
	assert( start >= 1 );
	assert( start <= ss.size() );
	assert( end <= ss.size() );
	assert( (!abego.size()) || ( start <= abego.size() ) );
	assert( (!abego.size()) || ( end <= abego.size() ) );
	core::Size const key_len = end - start + fragsize;
	std::string key = boost::lexical_cast< std::string >( fragsize ) + ss.substr(start-1, key_len);
	for ( core::Size i=start, c=1; i<=abego.size() && c<=key_len; ++i ) {
		assert( abego[i].size() == 1 );
		key += abego[i];
		++c;
	}
	if ( aa.size() ) {
		key += aa.substr(start-1, key_len);
	}
	return key;
}

} // components
} // denovo_design
} // protocols

