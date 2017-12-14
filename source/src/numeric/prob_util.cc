// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file numeric/prob_util.cc
/// @author Christopher Miles (cmiles@uw.edu)

// Unit header
#include <numeric/prob_util.hh>

// C/C++ headers
#include <iostream>
#include <fstream>
#include <string>

// Utility headers
#include <utility/exit.hh>
#include <utility/vector1.hh>

namespace numeric {

void read_probabilities_or_die(const std::string& filename, utility::vector1<double>* probs) {
	using namespace std;
	assert(probs);
	probs->clear();

	ifstream in(filename.c_str());
	if ( !in.is_open() ) {
		utility_exit_with_message("Error reading probabilities from " + filename);
	}

	double p;
	in >> p;
	while ( in.good() ) {
		probs->push_back(p);
		in >> p;
	}
	in.close();
}

void print_probabilities(const utility::vector1<double>& probs, std::ostream& out) {
	for ( unsigned i = 1; i <= probs.size(); ++i ) {
		out << "P(" << i << ") = " << probs[i] << std::endl;
	}
}

}  // namespace numeric
