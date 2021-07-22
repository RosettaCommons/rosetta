// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/io/MergeAndSplitBehaviorManager.cc
/// @brief   Declarations for MergeAndSplitBehaviorManager.
/// @author  Andy Watkins
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit header
#include <core/chemical/io/MergeAndSplitBehaviorManager.hh>

// Utility header
#include <utility/vector1.hh>

#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>
#include <utility/vector1.srlz.hh>

// Cereal headers
#include <cereal/types/string.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/types/vector.hpp>
#endif // SERIALIZATION

// Basic headers
#include <basic/Tracer.hh>


// Construct tracer.
static basic::Tracer TR( "core.chemical.io.MergeAndSplitBehaviorManager" );


namespace core {
namespace chemical {
namespace io {

// Public methods //////////////////////////////////////////////////////////////
// Standard methods ////////////////////////////////////////////////////////////
// Default constructor
MergeAndSplitBehaviorManager::MergeAndSplitBehaviorManager()
{}

// Standard constructor
MergeAndSplitBehaviorManager::MergeAndSplitBehaviorManager( std::string const & database_directory )
{
	using namespace std;
	using namespace io;

	string const & merge_filename( database_directory + "merge_residue_behaviors.txt" );
	MergeBehaviorMap merge_behaviors( read_merge_behaviors_from_database_file( merge_filename ) );
	merge_behaviors_.insert( merge_behaviors.begin(), merge_behaviors.end() );

	string const & split_filename( database_directory + "split_residue_behaviors.txt" );
	SplitBehaviorsMap split_behaviors( read_split_behaviors_from_database_file( split_filename ) );
	split_behaviors_.insert( split_behaviors.begin(), split_behaviors.end() );
}

// Destructor
MergeAndSplitBehaviorManager::~MergeAndSplitBehaviorManager() = default;


// Accessors ///////////////////////////////////////////////////////////////////
/// @return  If key not found, returns mrb_do_not_merge setting & an empty AtomRenamingMap.
ResidueMergeInstructions const &
MergeAndSplitBehaviorManager::merge_behavior_for_name3( std::string const & name3 ) const
{
	if ( merge_behaviors_.count( name3 ) ) {
		return merge_behaviors_.find( name3 )->second;
	}
	return NO_MERGE_BEHAVIOR_;
}

/// @return  If key not found, returns an empty SplitBehaviors type.
SplitBehaviors const &
MergeAndSplitBehaviorManager::split_behavior_for_name3( std::string const & name3 ) const
{
	if ( split_behaviors_.count( name3 ) ) {
		return split_behaviors_.find( name3 )->second;
	}
	return NO_SPLIT_BEHAVIORS_;
}

}  // namespace io
}  // namespace chemical
}  // namespace core


#ifdef    SERIALIZATION
/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::io::MergeAndSplitBehaviorManager::save( Archive & arc ) const
{
	arc( CEREAL_NVP( merge_behaviors_ ) );  // MergeBehaviorMap
	// EXEMPT NO_MERGE_BEHAVIOR_
	arc( CEREAL_NVP( split_behaviors_ ) );  // SplitBehaviorsMap
	// EXEMPT NO_SPLIT_BEHAVIORS_
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::io::MergeAndSplitBehaviorManager::load( Archive & arc )
{
	arc( merge_behaviors_ );  // MergeBehaviorMap
	// EXEMPT NO_MERGE_BEHAVIOR_
	arc( split_behaviors_ );  // SplitBehaviorsMap
	// EXEMPT NO_SPLIT_BEHAVIORS_
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::io::MergeAndSplitBehaviorManager );
#endif // SERIALIZATION
