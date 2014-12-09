// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/io/silent/EnergyNames.hh
///
/// @brief class representing a set of energy names used in silent-file input/output
/// @author James Thompson

#ifndef INCLUDED_core_io_silent_EnergyNames_hh
#define INCLUDED_core_io_silent_EnergyNames_hh

// mini headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <core/io/silent/SharedSilentData.hh>

// C++ Headers
#include <string>
#include <sstream>

#include <utility/vector1_bool.hh>


namespace core {
namespace io {
namespace silent {


class EnergyNames : public SharedSilentData {

public:
	EnergyNames() {}
	EnergyNames( std::string const & line ) {
		std::istringstream score_line_stream( line );
		std::string tag;
		//tr.debug << "reading score names from " << line << std::endl;
		score_line_stream >> tag; // score:
		if ( score_line_stream.fail() || tag != "SCORE:" ) {
			//tr.error << "bad format in score line of silent file" << std::endl;
			//tr.error << "tag = "  << tag << std::endl;
			//tr.error << "line = " << line << std::endl;
		}

		score_line_stream >> tag; // first score name
		while ( !score_line_stream.fail() ) {
			energy_names_.push_back( tag );
			score_line_stream >> tag;
		}
	} // EnergyNames( std::string const & line )

	void energy_names( utility::vector1< std::string > new_e_names ) {
		 energy_names_ = new_e_names;
	}

	utility::vector1< std::string > energy_names() const {
		return energy_names_;
	}

	bool operator == (const EnergyNames & other) const {
		using std::string;
		using utility::vector1;

		vector1< string > other_names = other.energy_names();
		if ( energy_names_.size() != other_names.size() ) {
			return false;
		}

		for ( vector1< string >::const_iterator it1 = energy_names_.begin(), it2 = other_names.begin(),
					end1 = energy_names_.end(), end2 = other_names.end();
					it1 != end1 && it2 != end2; ++it1, ++it2
		) {
			if ( *it1 != *it2 ) return false;
		}

		return true;
	} // operator ==


private:
	utility::vector1< std::string > energy_names_;
}; // EnergyNames


} // namespace silent
} // namespace io
} // namespace core

#endif
