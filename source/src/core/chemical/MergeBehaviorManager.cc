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

// Rosetta headers
#include <core/chemical/MergeBehaviorManager.hh>
#include <core/types.hh>

// Basic headers
#include <basic/database/open.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/file/FileName.hh>
#include <utility/io/izstream.hh>

// C++ headers
#include <fstream>
#include <string>
//#include <sstream>
#include <set>
#include <algorithm>

#include <utility/vector1.hh>
#include <utility/file/file_sys_util.hh>
#include <utility/io/util.hh>


#ifdef    SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/string.hpp>
#include <cereal/types/map.hpp>
#include <cereal/types/utility.hpp>
#endif // SERIALIZATION

namespace core {
namespace chemical {

static basic::Tracer TR( "core.chemical.MergeBehaviorManager" );

MergeBehaviorManager::AtomRenamingMap
mrb_map_from_correspondence( std::string const & correspondence ) {
	MergeBehaviorManager::AtomRenamingMap mrb_map;
	std::istringstream line_word_by_word( correspondence );

	while ( !line_word_by_word.fail() ) {
		std::string token;
		std::getline(line_word_by_word, token, ',');

		// split by -
		std::string key, value;
		std::istringstream tokstream( token );
		std::getline(tokstream, key, '-');
		std::getline(tokstream, value, '-');

		mrb_map[ key ] = value;
	}
	return mrb_map;
}

merge_residue_behavior
mrb_from_name( std::string const & mrb ) {
	if ( mrb == "do_not_merge" ) { return mrb_do_not_merge; }
	if ( mrb == "merge_w_prev" ) { return mrb_merge_w_prev; }
	if ( mrb == "merge_w_next" ) { return mrb_merge_w_next; }
	utility_exit_with_message( "Unable to convert string \"" + mrb + "\" to a merge residue behavior" );
	return mrb_do_not_merge;
}

MergeBehaviorManager::MergeBehaviorManager() :
	no_behavior_( std::make_pair( mrb_do_not_merge, AtomRenamingMap() ) )
{}

MergeBehaviorManager::MergeBehaviorManager( std::string const & database_directory ) :
	no_behavior_( std::make_pair( mrb_do_not_merge, AtomRenamingMap() ) )
{
	//AtomRenamingMap mp;
	//no_behavior_ = std::make_pair( mrb_do_not_merge, mp );

	//std::string const & filename( basic::database::full_name( database_directory + "merge_residue_behaviors.txt" ) );
	std::string const & filename( database_directory + "merge_residue_behaviors.txt" );
	MergeBehaviorMap merge_behaviors( read_merge_behaviors_from_database_file( filename ) );
	merge_behaviors_.insert( merge_behaviors.begin(), merge_behaviors.end() );
}

MergeBehaviorManager::~MergeBehaviorManager() = default;

MergeBehaviorManager::ResidueMergeInstructions const &
MergeBehaviorManager::merge_behavior_for_name3( std::string const & name3 ) const {
	if ( merge_behaviors_.count( name3 ) ) { return merge_behaviors_.find( name3 )->second; }
	return no_behavior_;
}


MergeBehaviorManager::MergeBehaviorMap
MergeBehaviorManager::read_merge_behaviors_from_database_file( std::string const & filename )
{
	using namespace std;
	using namespace core;
	using namespace utility;

	if ( ! utility::file::file_exists( filename ) ) return MergeBehaviorMap();

	vector1< string > const lines( utility::io::get_lines_from_file_data( filename ) );
	MergeBehaviorMap behavior_map;

	Size const n_lines( lines.size() );
	for ( uint i( 1 ); i <= n_lines; ++i ) {

		if ( lines[ i ][0] == '#' ) continue;

		istringstream line_word_by_word( lines[ i ] );
		string key;  // The map key is a PDB 3-letter code.
		string behavior;  // This is the behavior: do_not_merge, merge_w_prev, merge_w_next
		string correspondence;  // This is atom-by-atom correspondence

		std::getline(line_word_by_word, key, '\t');
		std::getline(line_word_by_word, behavior, '\t');
		std::getline(line_word_by_word, correspondence, '\t');

		behavior_map[ key ] = make_pair( mrb_from_name( behavior ), mrb_map_from_correspondence( correspondence ) );
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "Read " << behavior_map.size() << " merge behaviors from " << filename << '.' << endl;
	}

	return behavior_map;
}

}
}

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::chemical::MergeBehaviorManager::save( Archive & arc ) const {
	arc( CEREAL_NVP( merge_behaviors_ ) ); // MergeBehaviorMap
	arc( CEREAL_NVP( no_behavior_ ) ); // ResidueMergeInstructions
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::chemical::MergeBehaviorManager::load( Archive & arc ) {
	arc( merge_behaviors_ ); // MergeBehaviorMap
	arc( no_behavior_ ); // ResidueMergeInstructions
}

SAVE_AND_LOAD_SERIALIZABLE( core::chemical::MergeBehaviorManager );
#endif // SERIALIZATION
