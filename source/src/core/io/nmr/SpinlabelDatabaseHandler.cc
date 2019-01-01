// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/io/nmr/SpinlabelDatabaseHandler.cc
/// @brief   Implementation of class SpinlabelDatabaseHandler
/// @details last Modified: 10/07/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/io/nmr/SpinlabelDatabaseHandler.hh>

// Project headers

// Package headers
#include <core/types.hh>

// Utility headers
#include <utility/io/util.hh>
#include <utility/vector1.hh>
#include <utility/exit.hh>

// Basic header
#include <basic/Tracer.hh>
#include <basic/database/open.hh>

// C++ headers
#include <string>
#include <map>
#include <sstream>

namespace core {
namespace io {
namespace nmr {

static basic::Tracer TR("core.io.nmr.SpinlabelDatabaseHandler");

/// @brief Construct from full spinlabel name, 3-letter code,
///        1-letter code and radical atom name.
SpinlabelDatabaseEntry::SpinlabelDatabaseEntry(
	std::string const & fullname,
	std::string const & threelettercode,
	char const & onelettercode,
	std::string const & radicalatom
) :
	fullname_(fullname),
	three_letter_code_(threelettercode),
	one_letter_code_(onelettercode),
	radical_atom_(radicalatom),
	distance_potential_histogram_(""),
	ensemble_conformers_("")
{ }

/// @brief Destructor
SpinlabelDatabaseEntry::~SpinlabelDatabaseEntry() { }

/// @brief Empty default constructor
SpinlabelDatabaseHandler::SpinlabelDatabaseHandler() {
	read_in_database_file(basic::database::full_name( "scoring/nmr/spinlabel/spinlabel_properties.txt"), spinlabel_data_table_ );
}

/// Public member methods
/// @brief return map with paramagnetic ion data
SpinlabelDatabaseHandler::SpinlabelDatabaseMap const &
SpinlabelDatabaseHandler::get_spinlabel_data_table() const {
	return get_instance()->spinlabel_data_table_;
}

/// @brief return data for one paramagnetic ion
SpinlabelDatabaseEntry const &
SpinlabelDatabaseHandler::get_spinlabel_data(std::string const & spinlabel) const {
	SpinlabelDatabaseMap::const_iterator iter, iter_end;
	iter = get_instance()->spinlabel_data_table_.find(spinlabel);
	iter_end = get_instance()->spinlabel_data_table_.end();
	if ( iter == iter_end ) {
		utility_exit_with_message("ERROR: Entry for spinlabel " + spinlabel + " does not exist in spinlabel database.");
	}
	return (*iter).second;
}

/// @brief Some utility function to read in database file
void read_in_database_file(std::string const & filename, SpinlabelDatabaseHandler::SpinlabelDatabaseMap & table) {
	TR.Debug << "Opening spinlabel database file " << filename << " ... " << std::endl;
	utility::vector1< std::string > const lines( utility::io::get_lines_from_file_data( filename ) );
	Size n_lines(lines.size());
	for ( Size i = 1; i <= n_lines; ++i ) {
		std::istringstream line_words( lines[i] );
		std::string key, fullname, threelettercode, radicalatom, histogram_file, ensemble_conformers_file;
		char onelettercode;
		if ( !( line_words >> key >> fullname >> threelettercode >> onelettercode >> radicalatom >> histogram_file >> ensemble_conformers_file ) ) {
			TR.Warning << "Ignoring line " << i << " from file " << filename << "." << std::endl;
			continue;
		}
		if ( histogram_file == "N/A" || histogram_file == "n/a" ) {
			TR.Debug << "No histogram found for spinlabel with key " << key << ". A spline distance potential for scoring of atom pair distances may not be used." << std::endl;
		}
		if ( ensemble_conformers_file == "N/A" || ensemble_conformers_file == "n/a" ) {
			TR.Debug << "No conformers file found for spinlabel with key " << key << ". It cannot be used for PRE scoring in centroid mode." << std::endl;

		}
		SpinlabelDatabaseEntry entry(fullname, threelettercode, onelettercode, radicalatom);
		if ( histogram_file != "N/A" && histogram_file != "n/a" ) {
			entry.set_path_to_distance_potential_histogram( basic::database::full_name(histogram_file) );
		}
		if ( ensemble_conformers_file != "N/A" && ensemble_conformers_file != "n/a" ) {
			entry.set_path_to_ensemble_conformers( basic::database::full_name(ensemble_conformers_file) );
		}
		table.insert( std::make_pair(key, entry) );
	}
	TR.Debug << "Read " << table.size() << " entries from database file " << filename << "." << std::endl;
}

} // namespace nmr
} // namespace io
} // namespace core
