// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief
/// MergeBehaviorManager

#ifndef INCLUDED_core_chemical_MergeBehaviorManager_hh
#define INCLUDED_core_chemical_MergeBehaviorManager_hh


// Unit headers
#include <core/chemical/MergeBehaviorManager.fwd.hh>

// STL headers
#include <list>

#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <map>

namespace core {
namespace chemical {

class MergeBehaviorManager {
public:
	typedef std::map< std::string, std::string > AtomRenamingMap;
	typedef std::pair< merge_residue_behavior, AtomRenamingMap > ResidueMergeInstructions;
	typedef std::map< std::string, ResidueMergeInstructions > MergeBehaviorMap;

	MergeBehaviorManager( std::string const & database_directory );
	~MergeBehaviorManager();

	ResidueMergeInstructions const &
	merge_behavior_for_name3( std::string const & name3 ) const;

private:

	MergeBehaviorMap
	read_merge_behaviors_from_database_file( std::string const & filename );

private:
	MergeBehaviorMap merge_behaviors_;

	// no_behavior_ is a pair of do not merge with an empty map
	ResidueMergeInstructions no_behavior_;
};

}
}

#endif

