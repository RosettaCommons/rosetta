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
char
as_char(std::string const & value, char null);

/// @brief find the index for the given column name in the table.
/// If it can't be found, return a negative number
int find_gemmi_column(gemmi::cif::Table & table, std::string const & name);

///////////////  WRITING UTILS

void
normalize_table_name(std::string & table_name);

/// @brief Gets a table with the given name (cif category) and column names for writing
/// If the table does not exist, create it.
/// If it does already exist, make sure that all the provided column names are present.
///
/// (The returned table is simply a view to the underlying Block.)
///
/// Rows can be added with gemmi_append_row() below
gemmi::cif::Table
gemmi_get_table(gemmi::cif::Block & block, std::string table_name, std::vector<std::string> const & columns );

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
gemmi::cif::Loop &
gemmi_add_table(gemmi::cif::Block & block, std::string table_name, std::vector<std::string> const & columns);

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
