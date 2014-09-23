// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/inputter/PDBInputter.cc
/// @brief An inputter that takes a list of pdbs
/// @author Ken Jung

// Unit Headers
#include <protocols/inputter/PDBInputter.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <utility/lua/LuaIterator.hh>

#include <core/pose/util.hh>
#include <utility/string_util.hh>
// tracer
#include <basic/Tracer.hh>

namespace protocols {
namespace inputter {

static thread_local basic::Tracer TR( "protocols.inputter.PDBInputter" );

#ifdef USELUA
void lregister_PDBInputter( lua_State * lstate ) {
	lregister_Inputter( lstate );
	luabind::module(lstate, "protocols")
	[
		luabind::namespace_("inputter")
		[
			luabind::class_<PDBInputter, Inputter>("PDBInputter")
		]
	];
}
#endif

PDBInputter::~PDBInputter(){}

core::pose::PoseSP PDBInputter::get_nth_pose( int n ) {
	offset_ = true;
	core::chemical::ResidueTypeSetCOP residue_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
	core::pose::PoseSP tmppose(new core::pose::Pose());
	if( multiply_over_all_ ) {
		// we go through all the filenames once, then more times according to multiplier
		for( int i = 1; i<n; i++){
			file_names_[curr_idx_].first++;
			if( file_names_[curr_idx_].first >= multiplier_ )  {
				// this only can happen when curr_idx_ == 0
				file_names_.pop_front();
			} else {
				curr_idx_++;
				curr_idx_ = curr_idx_ == file_names_.size() ? 0 : curr_idx_;
			}
		}
		core::import_pose::pose_from_pdb( *tmppose, *residue_set, file_names_[curr_idx_].second );
		file_names_[curr_idx_].first++;
		core::pose::add_comment( *tmppose, "inputfile", file_names_[curr_idx_].second );
		core::pose::add_comment( *tmppose, "filemultiplier", utility::to_string(file_names_[curr_idx_].first) );
		if( file_names_[curr_idx_].first >= multiplier_ )  {
			// this only can happen when curr_idx_ == 0
			file_names_.pop_front();
		} else {
			curr_idx_++;
			curr_idx_ = curr_idx_ == file_names_.size()  ? 0 : curr_idx_;
		}
		return tmppose;
	} else {
		// we go through each file name MULTIPLIER times, then pop it
		for( int i = 1; i<n; i++){
			file_names_.front().first++;
			if( file_names_.front().first >= multiplier_ ) 
				file_names_.pop_front();
		}
		core::import_pose::pose_from_pdb( *tmppose, *residue_set, file_names_.front().second );
		file_names_.front().first++;
		core::pose::add_comment( *tmppose, "inputfile", file_names_.front().second );
		core::pose::add_comment( *tmppose, "file_multiplier", utility::to_string(file_names_.front().first) );
		if( file_names_.front().first >= multiplier_ ) 
			file_names_.pop_front();
		return tmppose;
	}
}

bool PDBInputter::has_nth_pose( int n ) {
	int sum = 0;
	for( core::Size i = 0; i < file_names_.size(); i++ ) {
		sum += (multiplier_ - file_names_[i].first	); 
		if (sum >= n )
			return true;
	}
	return false;
}

InputterSP PDBInputter::create() {
	return InputterSP( new PDBInputter() );
}

#ifdef USELUA
void PDBInputter::parse_def( utility::lua::LuaObject const & def,
				utility::lua::LuaObject const & tasks,
				utility::lua::LuaObject & inputters ) {
	if( def["multiplier"] ) multiplier_ = def["multiplier"].to<int>();	
	if( def["multiply_over_all"] ) multiply_over_all_ = def["multiply_over_all"].to<bool>();	
	for (utility::lua::LuaIterator i=def["filelist"].begin(), end; i != end; ++i) {
		file_names_.push_back( std::pair<int, std::string> ( 0, (*i).to<std::string>()) );
	}
}

void PDBInputter::lregister( lua_State * lstate ) {
	lregister_PDBInputter( lstate );
}
#endif
} // inputter
} // protocols
