// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/quota/ABEGO_SS_Map.cc
/// @brief map ABEGO vs. SS combinations
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/quota/ABEGO_SS_Map.hh>

// utility headers
#include <utility/vector1.hh>
#include <core/types.hh>

#include <sstream>

#include <utility/options/BooleanVectorOption.hh>


namespace protocols {
namespace frag_picker {
namespace quota {

/// @details Auto-generated virtual destructor
ABEGO_SS_Map::~ABEGO_SS_Map() {}

char ABEGO_SS_Map::all_abego_[6] = {'X','A','B','E','G','O'};
char ABEGO_SS_Map::all_ss_[4] = {'X','H','E','L'};

ABEGO_SS_Map::ABEGO_SS_Map(utility::vector1< std::pair<Size,Size> > ss_abego_types) {

	for ( Size i=1; i<=3; i++ ) {
		utility::vector1<bool> row;
		for ( Size j=1; j<=5; j++ ) {
			row.push_back(false);
		}
		ss_abego_types_.push_back(row);
	}
	for ( Size i=1; i<=ss_abego_types.size(); i++ ) {
		set_status(ss_abego_types[i],true);
	}
}

std::string ABEGO_SS_Map::show_valid() {

	std::ostringstream ss;
	for ( Size i=1; i<=3; i++ ) {
		for ( Size j=1; j<=5; j++ ) {
			if ( ! ss_abego_types_[i][j] ) continue;
			ss<< "("<<all_ss_[i]<<","<<all_abego_[j]<<") ";
		}
	}

	return ss.str();
}

Size torsion2big_bin_id(core::Real const phi,  core::Real const psi,  core::Real const omega) {

	if ( std::abs( omega ) < 90 ) {
		return 5; // cis-omega
	} else if ( phi >= 0.0 ) {
		if ( -100 < psi && psi <= 100 ) {
			return 4; // alpha-L
		} else {
			return 3; // E
		}
	} else {
		if ( -125 < psi && psi <= 50 ) {
			return 1; // helical
		} else {
			return 2; // extended
		}
	}
	return 0;
}

Size abego_index(char abego_class) {

	switch(abego_class) {
	case 'A' :
	case 'a' :
		return 1;
	case 'B' :
	case 'b' :
		return 2;
	case 'E' :
	case 'e' :
		return 3;
	case 'G' :
	case 'g' :
		return 4;
	case 'O' :
	case 'o' :
		return 5;
	default :
		std::stringstream ss (std::stringstream::out);
		ss << "[ERROR] Unrecognized abego class id: "<<abego_class<<"\n";
		utility_exit_with_message(ss.str());
		return 0;
	}
}

Size ss_index(char ss_class) {

	switch(ss_class) {
	case 'H' :
	case 'h' :
		return 1;
	case 'E' :
	case 'e' :
		return 2;
	case 'C' :
	case 'c' :
	case 'L' :
	case 'l' :
		return 3;
	default :
		std::stringstream ss (std::stringstream::out);
		ss << "[ERROR] Unrecognized ss class id: "<<ss_class<<"\n";
		utility_exit_with_message(ss.str());
		return 0;
	}
}

} // quota
} // frag_picker
} // protocols

