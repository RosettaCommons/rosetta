// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/io/nmr/ParamagneticDatabaseHandler.cc
/// @brief   Function definitions of ParamagneticDatabaseHandler class
/// @details last Modified: 05/11/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

// Unit headers
#include <core/io/nmr/ParamagneticDatabaseHandler.hh>

// Package headers
#include <core/io/nmr/ParaIon.hh>

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
#include <cmath>


namespace core {
namespace io {
namespace nmr {

static basic::Tracer TR("core.io.nmr.ParamagneticDatabaseHandler");

/// @brief Empty constructor
ParamagneticDatabaseHandler::ParamagneticDatabaseHandler() :
	ion_data_table_( read_in_database_file( basic::database::full_name( "chemical/element_sets/default/para_ion_properties.txt" ) ) )
{}

/// Public member methods
/// @brief return map with paramagnetic ion data
std::map< std::string, ParaIon > const &
ParamagneticDatabaseHandler::get_ion_data_table() const {
	return get_instance()->ion_data_table_;
}

/// @brief return data for one paramagnetic ion
ParaIon
ParamagneticDatabaseHandler::get_ion_data(std::string const & ion) {
	return get_instance()->ion_data_table_[ion];
}

/// @brief Some utility function to read in database file
std::map< std::string, ParaIon >
read_in_database_file(std::string const & filename) {
	utility::vector1< std::string > const lines( utility::io::get_lines_from_file_data( filename ) );
	std::map< std::string, ParaIon > data_table;
	Size n_lines(lines.size());
	for ( Size i = 1; i <= n_lines; ++i ) {
		std::istringstream line_words( lines[i] );
		std::string label;
		core::Real chg, S, L, J, gJ, tau_e;
		if ( !(line_words >> label >> chg >> S >> L >> J >> gJ >> tau_e) ) {
			TR.Warning << "Ignoring line " << i << " from file " << filename << "." << std::endl;
			continue;
		}
		data_table.insert(std::pair<std::string, ParaIon>(label, ParaIon(label, chg, S, L, J, gJ, tau_e)));
	}
	if ( TR.Debug.visible() ) {
		TR << "Read " << data_table.size() << " entries from database file " << filename << "." << std::endl;
	}
	return data_table;
}

} // namespace nmr
} // namespace io
} // namespace core
