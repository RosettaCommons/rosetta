// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    core/chemical/io/merge_and_split_behaviors_io.cc
/// @brief   Database input/output function definitions for residue-splitting behaviors.
/// @author  Labonte <JWLabonte@jhu.edu>


// Project header
#include <core/types.hh>

// Unit header
#include <core/chemical/io/merge_and_split_behaviors_io.hh>  // also includes typedefs for custom data structures

// Utility headers
#include <utility/io/util.hh>
#include <utility/file/file_sys_util.hh>

// Basic header
#include <basic/Tracer.hh>

// C++ header
#include <sstream>


// Construct tracer.
static basic::Tracer TR( "core.chemical.io.merge_and_split_behaviors_io" );


namespace core {
namespace chemical {
namespace io {

// Parse a string of comma-delimited pairs of old PDB atom names to the new Rosetta atom names,
// separated by a hyphen.
// For example, " C3 - C1 , C4 - C2 " means that atom " C3 " is renamed to " C1 " and " C4 " to " C2 ".
AtomRenamingMap
get_atom_renamings( std::string const & instructions ) {
	using namespace std;

	AtomRenamingMap atom_renamings;
	istringstream renamings_atoms_by_atoms( instructions );

	string token;
	getline( renamings_atoms_by_atoms, token, ',' );
	while ( ! renamings_atoms_by_atoms.fail() ) {
		// Split by -
		string key, value;
		istringstream tokstream( token );
		getline( tokstream, key, '-' );
		getline( tokstream, value, '-' );

		atom_renamings[ key ] = value;

		getline( renamings_atoms_by_atoms, token, ',' );
	}
	return atom_renamings;
}

// Convert behavior read from file as string to enum value.
merge_residue_behavior
mrb_from_name( std::string const & mrb ) {
	if ( mrb == "do_not_merge" ) { return mrb_do_not_merge; }
	if ( mrb == "merge_w_prev" ) { return mrb_merge_w_prev; }
	if ( mrb == "merge_w_next" ) { return mrb_merge_w_next; }
	utility_exit_with_message( "Unable to convert string \"" + mrb + "\" to a merge residue behavior" );
}

// Database files for merge behaviors are stored in the Rosetta database in the same folder as ResidueTypeSets.
MergeBehaviorMap
read_merge_behaviors_from_database_file( std::string const & filename )
{
	using namespace std;
	using namespace core;
	using namespace utility;

	if ( ! file::file_exists( filename ) ) {
		TR.Debug << "File " << filename << " does not exist." << endl;
		return MergeBehaviorMap();
	}

	vector1< string > const lines( utility::io::get_lines_from_file_data( filename ) );
	MergeBehaviorMap behavior_map;

	for ( auto line : lines ) {
		istringstream line_word_by_word( line );
		string key;  // The map key is a PDB 3-letter code.
		string behavior;  // This is the behavior: do_not_merge, merge_w_prev, merge_w_next
		string correspondence;  // This is atom-by-atom correspondence

		getline( line_word_by_word, key, '\t') ;
		getline( line_word_by_word, behavior, '\t' );
		getline( line_word_by_word, correspondence, '\t' );

		behavior_map[ key ] = make_pair( mrb_from_name( behavior ), get_atom_renamings( correspondence ) );
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "Read " << behavior_map.size() << " merge behaviors from " << filename << '.' << endl;
	}

	return behavior_map;
}


// Parse groups of atom-renaming instructions.
utility::vector1< AtomRenamingMap >
get_groups_of_atom_renamings( utility::vector1< std::string > const & instructions )
{
	using namespace std;
	using namespace utility;

	Size const n_instructions( instructions.size() );
	vector1< AtomRenamingMap > atom_renamings_groups( n_instructions );

	for ( core::uint i( 1 ); i <= n_instructions; ++i ) {
		atom_renamings_groups[ i ] = get_atom_renamings( instructions[ i ] );
	}

	return atom_renamings_groups;
}

// Parse instructions for splitting PDB residues.
// <residue_names> is a semicolon-delimited list of base residue types,
// given as a comma-separated pair of Rosetta 3-letter code and residue name,
// (NOT 3-letter code,) i.e., what one could enter in a Rosetta-formatted HETNAM record,
// for example: "Glc,->3)-alpha-D-Glcp;Gal,->2)-beta-D-Galp".
// <behavior> is a semicolon-delimited grouping of sets of atom renaming instructions.
// The sets are ordered in the same manner as the names in residue_names;
// that is, the first residue named has its atoms renamed first.
SplitBehaviors
get_SplitBehaviors( std::string const & residue_names, std::string const & behavior )
{
	using namespace std;
	using namespace utility;

	istringstream names_pair_by_pair( residue_names );
	istringstream instructions_res_by_res( behavior );
	string name_pair, instruction;
	vector1< string > names, instructions;

	getline( names_pair_by_pair, name_pair, ';' );
	while ( ! names_pair_by_pair.fail() ) {
		names.push_back( name_pair );
		getline( names_pair_by_pair, name_pair, ';' );
	}

	getline( instructions_res_by_res, instruction, ';' );
	while ( ! instructions_res_by_res.fail() ) {
		instructions.push_back( instruction );
		getline( instructions_res_by_res, instruction, ';' );
	}

	vector1< pair< string, string > > residue_pairs;
	for ( auto residue_pair : names ) {
		string code, base_name;
		istringstream residues_name_by_name( residue_pair );
		getline( residues_name_by_name, code, ',' );
		getline( residues_name_by_name, base_name, ',' );
		residue_pairs.push_back( make_pair( code, base_name ) );
	}

	return make_pair( residue_pairs, get_groups_of_atom_renamings( instructions ) );
}

// Return a mapping of PDB 3-letter codes to a set of splitting instructions.
// Database files for split behaviors are stored in the Rosetta database in the same folder as ResidueTypeSets.
SplitBehaviorsMap
read_split_behaviors_from_database_file( std::string const & filename )
{
	using namespace std;
	using namespace utility;
	using namespace utility::io;

	if ( ! file::file_exists( filename ) ) {
		TR.Debug << "File " << filename << " does not exist." << endl;
		return SplitBehaviorsMap();
	}

	vector1< string > const lines( get_lines_from_file_data( filename ) );
	SplitBehaviorsMap behaviors_map;

	for ( auto line : lines ) {
		istringstream line_column_by_column( line );

		string key;  // The map key is a PDB 3-letter code (old residue).
		string residue_names;  // comma-delimited list of new base residue types, by Rosetta residue name
		string behavior;  // semicolon-delimited list of atom-renaming instructions

		// The three columns in the file are tab-delimited.
		getline( line_column_by_column, key, '\t' );
		getline( line_column_by_column, residue_names, '\t' );
		getline( line_column_by_column, behavior, '\t' );

		behaviors_map[ key ] = get_SplitBehaviors( residue_names, behavior );
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "Read " << behaviors_map.size() << " residue-splitting behaviors from " << filename << '.' << endl;
	}

	return behaviors_map;
}

}  // namespace io
}  // namespace chemical
}  // namespace core
