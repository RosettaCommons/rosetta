// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/picking_old/FragmentLibraryManager.cc
/// @brief  singleton class for accessing fragment libraries
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <core/fragment/picking_old/FragmentLibraryManager.hh>

// package headers
#include <core/fragment/picking_old/vall/VallLibrary.hh>
#include <core/fragment/picking_old/vall/vall_io.hh>

// project headers
#include <basic/database/open.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <utility/vector1.hh>
#include <utility/thread/threadsafe_creation.hh>

// Boost headers
#include <boost/bind.hpp>
#include <boost/function.hpp>

namespace core {
namespace fragment {
namespace picking_old {

/// @brief default constructor
FragmentLibraryManager::FragmentLibraryManager() :
	vall_( NULL )
{}

/// @brief return instance of standard Vall library
/// This should be called in the constructor, and the Vall should
/// not be cleared from memory after it's loaded.
vall::VallLibrary const & FragmentLibraryManager::get_Vall() const {
	using namespace basic::options::OptionKeys;
	using basic::options::option;
	using basic::database::full_name;

	// TODO: we need proper locking support here for multi-threaded access
	if ( vall_ == NULL ) {
		vall_ = new vall::VallLibrary();
		vall::vall_library_from_file( full_name( option[ in::file::vall ][1] ), *vall_, 62471 );
	}

	return *vall_;
}


/// @brief clear standard Vall library from memory
/// WARNING WARNING WARNING! THREAD UNSAFE!
void FragmentLibraryManager::clear_Vall() {
	if ( vall_ != NULL ) {
		delete vall_;
		vall_ = NULL;
	}
}


} // picking_old
} // fragment
} // core
