// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/io/nmr/ParamagneticDatabaseHandler.hh
/// @brief   classes to handle ion properties database file
/// @details last Modified: 05/11/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_io_nmr_ParamagneticDatabaseHandler_HH
#define INCLUDED_core_io_nmr_ParamagneticDatabaseHandler_HH

// Unit headers
#include <core/io/nmr/ParamagneticDatabaseHandler.fwd.hh>

// Package headers
#include <core/io/nmr/ParaIon.hh>

// Package headers
#include <core/types.hh>

// Utility headers
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>
#include <map>

namespace core {
namespace io {
namespace nmr {

class ParamagneticDatabaseHandler : public utility::SingletonBase< ParamagneticDatabaseHandler > {
	friend class utility::SingletonBase< ParamagneticDatabaseHandler >;

private: // Methods

	/// @brief Empty default constructor
	ParamagneticDatabaseHandler();

private: // Data

	/// @brief map of ion label and ion properties
	std::map< std::string, ParaIon > ion_data_table_;

public: // Methods

	/// @brief return map with paramagnetic ion data
	std::map< std::string, ParaIon > const & get_ion_data_table() const;

	/// @brief return data for one paramagnetic ion
	ParaIon get_ion_data(std::string const & ion);
};

/// @brief Some utility function to read in database file
std::map< std::string, ParaIon >
read_in_database_file(std::string const & filename);

} // namespace nmr
} // namespace io
} // namespace core

#endif // INCLUDED_core_io_nmr_ParamagneticDatabaseHandler_HH
