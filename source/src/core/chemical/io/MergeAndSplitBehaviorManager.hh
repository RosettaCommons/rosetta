// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/io/MergeAndSplitBehaviorManager.hh
/// @brief   Declarations for MergeAndSplitBehaviorManager.
/// @author  Andy Watkins
/// @author  Labonte <JWLabonte@jhu.edu>


#ifndef INCLUDED_core_chemical_io_MergeAndSplitBehaviorManager_HH
#define INCLUDED_core_chemical_io_MergeAndSplitBehaviorManager_HH


// Unit headers
#include <core/chemical/io/MergeAndSplitBehaviorManager.fwd.hh>
#include <core/chemical/io/merge_and_split_behaviors_io.hh>  // typedefs for custom data structures

// Utility header
#include <utility/VirtualBase.hh>


namespace core {
namespace chemical {
namespace io {

class MergeAndSplitBehaviorManager : public utility::VirtualBase {
public:  // Standard methods ///////////////////////////////////////////////////
	/// @brief  Default constructor
	MergeAndSplitBehaviorManager();

	/// @brief  Standard constructor
	MergeAndSplitBehaviorManager( std::string const & database_directory );

	// Destructor
	~MergeAndSplitBehaviorManager() override;


public:  // Accessors //////////////////////////////////////////////////////////
	/// @brief  What is the merge behavior for this residue by PDB 3-letter code?
	ResidueMergeInstructions const & merge_behavior_for_name3( std::string const & name3 ) const;

	/// @brief  What is the split behavior for this residue by PDB 3-letter code?
	SplitBehaviors const & split_behavior_for_name3( std::string const & name3 ) const;


private:  // Private data //////////////////////////////////////////////////////
	MergeBehaviorMap merge_behaviors_;
	ResidueMergeInstructions const NO_MERGE_BEHAVIOR_ = ResidueMergeInstructions( mrb_do_not_merge, AtomRenamingMap() );

	SplitBehaviorsMap split_behaviors_;
	SplitBehaviors const NO_SPLIT_BEHAVIORS_ = SplitBehaviors();

#ifdef    SERIALIZATION
public:
	template< class Archive > void save( Archive & arc ) const;
	template< class Archive > void load( Archive & arc );
#endif  // SERIALIZATION

};  // class MergeAndSplitBehaviorManager

}  // namespace io
}  // namespace chemical
}  // namespace core

#endif  // INCLUDED_core_chemical_io_MergeAndSplitBehaviorManager_HH
