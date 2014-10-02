// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragment/picking_old/FragmentLibraryManager.hh
/// @brief  singleton class for accessing fragment libraries
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_FragmentLibraryManager_hh
#define INCLUDED_core_fragment_picking_old_FragmentLibraryManager_hh

// unit headers
#include <core/fragment/picking_old/FragmentLibraryManager.fwd.hh>

// package headers
#include <core/fragment/picking_old/vall/VallLibrary.fwd.hh>

// Utility headers
#include <utility/SingletonBase.hh>

namespace core {
namespace fragment {
namespace picking_old {


/// @brief  singleton class for accessing fragment libraries
class FragmentLibraryManager : public utility::SingletonBase< FragmentLibraryManager >
{
public:
	friend class utility::SingletonBase< FragmentLibraryManager >;

private: // constructor


	/// @brief default constructor
	FragmentLibraryManager();

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static FragmentLibraryManager * create_singleton_instance();

public: // access

	/// @brief return instance of standard Vall library
	vall::VallLibrary const & get_Vall() const;


public: // memory management


	/// @brief clear standard Vall library from memory
	/// WARNING WARNING WARNING. THREAD UNSAFE!
	void clear_Vall();

private: // data

	// *** WARNING -- pointers for all libraries below must be
	// initialized to NULL in constructor ***

	/// @brief standard Vall library
	mutable vall::VallLibrary * vall_;


};


} // picking_old
} // fragment
} // core


#endif /* INCLUDED_core_fragment_picking_old_FragmentLibraryManager_HH */
