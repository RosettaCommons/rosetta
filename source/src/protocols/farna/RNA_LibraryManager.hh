// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/farna/RNA_LibraryManager.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_farna_RNA_LibraryManager_HH
#define INCLUDED_protocols_farna_RNA_LibraryManager_HH

#include <protocols/farna/RNA_LibraryManager.fwd.hh>
#include <protocols/farna/RNA_Fragments.fwd.hh>
#include <protocols/farna/RNA_JumpLibrary.fwd.hh>
#include <protocols/farna/BasePairStepLibrary.fwd.hh>
#include <utility/SingletonBase.hh>
#include <map>
#include <string>

namespace protocols {
namespace farna {

/// @brief Holds all libraries relevant to FARFAR as CONST copies:
///
///           JumpLibrary
///           FullAtomRNA_Fragments
///           BasePairStepLibrary, ...
///
///  However, those const libraries include some mutable data to allow for efficient, lazy loading --
///  need to put mutexes around those functions that update those data.

class RNA_LibraryManager : public utility::SingletonBase< RNA_LibraryManager >
{

public:
	friend class utility::SingletonBase< RNA_LibraryManager >;

private:

	/// @brief private constructor
	RNA_LibraryManager(){};

	/// @brief private singleton creation function to be used with
	/// utility::thread::threadsafe_singleton
	static RNA_LibraryManager * create_singleton_instance();

public:

	RNA_Fragments const &
	rna_fragment_library( std::string const & tag );

	RNA_JumpLibrary const &
	rna_jump_library( std::string const & tag );

	RNA_JumpLibraryCOP const &
	rna_jump_library_cop( std::string const & tag );

	BasePairStepLibrary const &
	canonical_base_pair_step_library();

	BasePairStepLibrary const &
	general_base_pair_step_library();

private:

	std::map< std::string, RNA_FragmentsCOP > rna_fragment_libraries_;
	std::map< std::string, RNA_JumpLibraryCOP > rna_jump_libraries_;
	BasePairStepLibraryCOP canonical_base_pair_step_library_;
	BasePairStepLibraryCOP general_base_pair_step_library_;

};

} //farna
} //protocols

#endif
