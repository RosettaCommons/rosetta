// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/io/nmr/ParamagneticDatabaseHandler.hh
/// @brief   classes to handle spinlabel properties database file
/// @details last Modified: 10/07/16
/// @author  Georg Kuenze (georg.kuenze@vanderbilt.edu)

#ifndef INCLUDED_core_io_nmr_SpinlabelDatabaseHandler_HH
#define INCLUDED_core_io_nmr_SpinlabelDatabaseHandler_HH

// Unit headers
#include <core/io/nmr/SpinlabelDatabaseHandler.fwd.hh>

// Project headers

// Package headers
#include <core/types.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <utility/SingletonBase.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>
#include <map>

namespace core {
namespace io {
namespace nmr {

class SpinlabelDatabaseEntry : public utility::pointer::ReferenceCount {

private: // Data

	/// @brief Empty default constructor
	SpinlabelDatabaseEntry();

public: // Methods

	/// @brief Construct from full spinlabel name, 3-letter code,
	///        1-letter code and radical atom name.
	SpinlabelDatabaseEntry(
		std::string const & fullname,
		std::string const & threelettercode,
		char const & onelettercode,
		std::string const & radicalatom
	);

	/// @brief Destructor
	~SpinlabelDatabaseEntry();

	std::string const & fullname() const { return fullname_; }
	std::string const & three_letter_code() const { return three_letter_code_; }
	char const & one_letter_code() const { return one_letter_code_; }
	std::string const & radical_atom() const { return radical_atom_; }
	std::string const & distance_potential_histogram() const { return distance_potential_histogram_; }
	std::string const & ensemble_conformers() const { return ensemble_conformers_; }

	void set_path_to_distance_potential_histogram(std::string const & filename) { distance_potential_histogram_ = filename; }
	void set_path_to_ensemble_conformers(std::string const & filename) { ensemble_conformers_ = filename; }

private: // Data

	std::string fullname_;
	std::string three_letter_code_;
	char one_letter_code_;
	std::string radical_atom_;
	// Paths to database files
	std::string distance_potential_histogram_;
	std::string ensemble_conformers_;

};


class SpinlabelDatabaseHandler : public utility::SingletonBase< SpinlabelDatabaseHandler > {
	friend class utility::SingletonBase< SpinlabelDatabaseHandler >;

public:
	typedef std::map< std::string, SpinlabelDatabaseEntry > SpinlabelDatabaseMap;

private: // Methods

	/// @brief Empty default constructor
	SpinlabelDatabaseHandler();

private: // Data

	/// @brief map of ion label and ion properties
	SpinlabelDatabaseMap spinlabel_data_table_;

public: // Methods

	/// @brief return map with spinlabel database entries
	SpinlabelDatabaseMap const & get_spinlabel_data_table() const;

	/// @brief return data for spinlabel
	SpinlabelDatabaseEntry const & get_spinlabel_data(std::string const & spinlabel) const;
};

/// @brief Some utility function to read in database file
void
read_in_database_file(
	std::string const & filename,
	SpinlabelDatabaseHandler::SpinlabelDatabaseMap & table
);

} // namespace nmr
} // namespace io
} // namespace core

#endif // INCLUDED_core_io_nmr_SpinlabelDatabaseHandler_HH
