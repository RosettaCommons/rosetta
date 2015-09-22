// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/inputter/SilentFileInputter.cc
/// @brief An inputter that takes a list of silent files
/// @author Ken Jung

// Unit Headers
#include <protocols/inputter/SilentFileInputter.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/pose/Pose.hh>
#include <utility/lua/LuaIterator.hh>

#include <core/pose/util.hh>
#include <utility/string_util.hh>
// tracer
#include <basic/Tracer.hh>

namespace protocols {
namespace inputter {

static THREAD_LOCAL basic::Tracer TR( "protocols.inputter.SilentFileInputter" );

#ifdef USELUA
void lregister_SilentFileInputter( lua_State * lstate ) {
	lregister_Inputter( lstate );
	luabind::module(lstate, "protocols")
	[
		luabind::namespace_("inputter")
		[
			luabind::class_<SilentFileInputter, Inputter>("SilentFileInputter")
		]
	];
}
#endif

SilentFileInputter::~SilentFileInputter(){}

core::pose::PoseSP SilentFileInputter::get_nth_pose( int n ) {
	offset_ = true;
	core::chemical::ResidueTypeSetCOP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	core::pose::PoseSP tmppose(new core::pose::Pose());
	if ( multiply_over_all_ ) {
		// we go through all the filenames once, then more times according to multiplier
		for ( int i = 1; i<n; i++ ) {
			tags_[curr_idx_].first++;
			if ( tags_[curr_idx_].first >= multiplier_ )  {
				// this only can happen when curr_idx_ == 0
				tags_.pop_front();
			} else {
				curr_idx_++;
				curr_idx_ = curr_idx_ == tags_.size() ? 0 : curr_idx_;
			}
		}
		sfd_[ tags_[curr_idx_].second ]->fill_pose( *tmppose, *residue_set );
		tags_[curr_idx_].first++;
		core::pose::add_comment( *tmppose, "inputfile", tags_[curr_idx_].second );
		core::pose::add_comment( *tmppose, "filemultiplier", utility::to_string(tags_[curr_idx_].first) );
		if ( tags_[curr_idx_].first >= multiplier_ )  {
			// this only can happen when curr_idx_ == 0
			tags_.pop_front();
		} else {
			curr_idx_++;
			curr_idx_ = curr_idx_ == tags_.size()  ? 0 : curr_idx_;
		}
		return tmppose;
	} else {
		// we go through each file name MULTIPLIER times, then pop it
		for ( int i = 1; i<n; i++ ) {
			tags_.front().first++;
			if ( tags_.front().first >= multiplier_ ) {
				tags_.pop_front();
			}
		}
		sfd_[ tags_.front().second ]->fill_pose( *tmppose, *residue_set );
		tags_.front().first++;
		core::pose::add_comment( *tmppose, "inputfile", tags_.front().second );
		core::pose::add_comment( *tmppose, "file_multiplier", utility::to_string(tags_.front().first) );
		if ( tags_.front().first >= multiplier_ ) {
			tags_.pop_front();
		}
		return tmppose;
	}
}

bool SilentFileInputter::has_nth_pose( int n ) {
	int sum = 0;
	for ( core::Size i = 0; i < tags_.size(); i++ ) {
		sum += (multiplier_ - tags_[i].first );
		if ( sum >= n ) {
			return true;
		}
	}
	return false;
}

InputterSP SilentFileInputter::create() {
	return InputterSP( new SilentFileInputter() );
}

#ifdef USELUA
void SilentFileInputter::parse_def( utility::lua::LuaObject const & def,
				utility::lua::LuaObject const & tasks,
				utility::lua::LuaObject & inputters ) {
	if( def["multiplier"] ) multiplier_ = def["multiplier"].to<int>();	
	if( def["multiply_over_all"] ) multiply_over_all_ = def["multiply_over_all"].to<bool>();	
	for (utility::lua::LuaIterator i=def["filelist"].begin(), end; i != end; ++i) {
		sfd_.read_file( (*i).to<std::string>() );
	}
	utility::vector1< std::string > tags = sfd_.tags();
	for( utility::vector1< std::string >::iterator itr = tags.begin(), end = tags.end(); itr != end; ++itr ) {
		tags_.push_back( std::pair<int, std::string> ( 0, *itr) );
	}
}

void SilentFileInputter::lregister( lua_State * lstate ) {
	lregister_SilentFileInputter( lstate );
}
#endif

} // inputter
} // protocols
