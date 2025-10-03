// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   utility/gemmi_util.hh
/// @brief  Utilities for working with Gemmi CIF file data
///
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_utility_gemmi_util_HH
#define INCLUDED_utility_gemmi_util_HH

// Unit headers
#include <utility/gemmi_util.fwd.hh>
#include <gemmi/cif.hpp>

#include <string>
#include <vector>
#include <initializer_list>

namespace utility {

///////////////  READING UTILS

/// As the default as_char() is not robust to empty strings
inline
char
as_char(std::string const & value, char null) {
	if ( value.size() == 0 ) {
		return null;
	}
	if ( value.size() == 2 && (value == "''"  || value == "\"\"" ) ) {
		return null;
	}
	return gemmi::cif::as_char(value, null);
}

/// @brief find the index for the given column name in the table.
/// If it can't be found, return a negative number
inline int find_gemmi_column(gemmi::cif::Table & table, std::string const & name) {
	if ( table.width() == 0 ) { return -1; } // No columns
	gemmi::cif::Table::Row const & tags = table.tags();
	for ( int ii = 0; ii < int(tags.size()); ++ii ) {
		std::string const & tag = tags[ii];
		if ( name == tag || name == tag.substr( table.prefix_length ) ) {
			return ii;
		}
	}
	return -1;
}

///////////////  WRITING UTILS

inline
void
normalize_table_name(std::string & table_name) {
	if ( table_name.size() == 0 ) {
		table_name = "_TABLE.";
	}
	if ( table_name[0] != '_' ) {
		table_name = '_' + table_name;
	}
	if ( table_name[ table_name.size()-1 ] != '.' ) {
		table_name = table_name + '.';
	}
}

/// @brief Gets a table with the given name (cif category) and column names for writing
/// If the table does not exist, create it.
/// If it does already exist, make sure that all the provided column names are present.
///
/// (The returned table is simply a view to the underlying Block.)
///
/// Rows can be added with gemmi_append_row() below
inline
gemmi::cif::Table
gemmi_get_table(gemmi::cif::Block & block, std::string table_name, std::vector<std::string> const & columns ) {
	normalize_table_name(table_name);

	if ( ! block.has_mmcif_category(table_name) ) {
		block.init_loop(table_name, columns);
	} else {
		// Already has the table, make sure we have the columns.
		block.find_mmcif_category(table_name).ensure_loop();
		gemmi::cif::Loop * loop = block.find_mmcif_category(table_name).get_loop();
		runtime_assert( loop != nullptr );
		std::vector<std::string> new_tags;
		for ( std::string const & tag: columns ) {
			if ( ! loop->has_tag( table_name + tag ) ) {
				new_tags.push_back( table_name + tag );
			}
		}
		if ( ! new_tags.empty() ) {
			loop->add_columns( new_tags, "?" );
		}
	}

	return block.find( table_name, columns );
}

template < class Iterable >
void
gemmi_append_row(gemmi::cif::Table & table, Iterable const & values ) {
	std::vector< std::string > quoted;
	for ( auto iter(values.begin()); iter != values.end(); ++iter ) {
		if ( *iter == "?" || *iter == "." ) {
			// Assume we actually want it as a null, rather than a quoted question mark
			quoted.push_back( *iter );
		} else {
			quoted.push_back( gemmi::cif::quote( *iter ) );
		}
	}
	table.append_row(quoted);
}

inline
void
gemmi_append_row(gemmi::cif::Table & table, std::initializer_list<std::string> const & init_list ) {
	gemmi_append_row< std::initializer_list<std::string> >(table, init_list);
}

/// @brief Adds a new table (actually a 'Loop' object) to the given block, with the given column names
/// Returns a reference to the newly added loop object (which can be augmented with the `gemmi_add_row()` function)
inline
gemmi::cif::Loop &
gemmi_add_table(gemmi::cif::Block & block, std::string table_name, std::vector<std::string> const & columns) {
	if ( table_name.size() == 0 ) {
		table_name = "_TABLE.";
	}
	if ( table_name[0] != '_' ) {
		table_name = '_' + table_name;
	}
	if ( table_name[ table_name.size()-1 ] != '.' ) {
		table_name = table_name + '.';
	}

	return block.init_loop(table_name, columns);
}

/// @brief Adds a row to the table. Takes care of quoting the entries properly
template< class Iterable >
void
gemmi_add_row(gemmi::cif::Loop & loop, Iterable const & value ) {
	std::vector< std::string > quoted;
	for ( auto iter(value.begin()); iter != value.end(); ++iter ) {
		if ( *iter == "?" || *iter == "." ) {
			// Assume we actually want it as a null, rather than a quoted question mark
			quoted.push_back( *iter );
		} else {
			quoted.push_back( gemmi::cif::quote( *iter ) );
		}
	}
	loop.add_row( quoted );
}

inline
void
gemmi_add_row(gemmi::cif::Loop & loop, std::initializer_list<std::string> const & init_list ) {
	gemmi_add_row< std::initializer_list<std::string> >(loop, init_list);
}


} // namespace utility


#endif // INCLUDED_utility_gemmi_util_HH
