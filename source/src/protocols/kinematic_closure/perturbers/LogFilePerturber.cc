// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Unit headers
#include <protocols/kinematic_closure/types.hh>
#include <protocols/kinematic_closure/ClosureProblem.hh>
#include <protocols/kinematic_closure/perturbers/LogFilePerturber.hh>

// Core headers
#include <core/pose/Pose.hh>

// Protocol headers
#include <protocols/loops/Loop.hh>

// Numeric headers
#include <numeric/conversions.hh>

// Boost headers
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>

// C++ headers
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <iterator>

#define foreach BOOST_FOREACH

namespace protocols {
namespace kinematic_closure {
namespace perturbers {

using namespace std;

LogFilePerturber::LogFilePerturber(string path) {
	using numeric::conversions::to_radians;

	string line;
	ifstream file (path.c_str());

	if ( !file.is_open() ) {
		utility_exit_with_message("LogFilePerturber can't open file: " + path);
	}

	while ( getline(file, line) ) {
		string header;
		list<string> fields;

		boost::trim(line);
		boost::split(fields, line, boost::is_any_of("\t "));
		header = fields.front();
		fields.pop_front();

		ParameterList dofs;

		transform(
			fields.begin(),
			fields.end(),
			back_inserter(dofs),
			boost::lexical_cast<double, string>);

		if ( header == "TORSIONS" ) {
			for ( Size i = 1; i <= dofs.size(); i++ ) to_radians(dofs[i]);
			torsion_angles_.push_back(dofs);
		}
		if ( header == "ANGLES" ) {
			for ( Size i = 1; i <= dofs.size(); i++ ) to_radians(dofs[i]);
			bond_angles_.push_back(dofs);
		}
		if ( header == "LENGTHS" ) {
			bond_lengths_.push_back(dofs);
		}
	}

	file.close();
	iteration_ = 1;
}

void LogFilePerturber::perturb_subset(
	Pose const &, IndexList const &, ClosureProblemOP problem) {

	if ( iteration_ > torsion_angles_.size() ) {
		utility_exit_with_message("LogFilePerturber ran out of frames!");
	}

	copy(
		torsion_angles_[iteration_].begin(),
		torsion_angles_[iteration_].end(),
		problem->perturb_torsions().begin());

	copy(
		bond_angles_[iteration_].begin(),
		bond_angles_[iteration_].end(),
		problem->perturb_angles().begin());

	copy(
		bond_lengths_[iteration_].begin(),
		bond_lengths_[iteration_].end(),
		problem->perturb_lengths().begin());

	iteration_ += 1;
}

void LogFilePerturber::log_torsions(ostream & out, ParameterList torsions) {
	out << "TORSIONS ";
	foreach ( Real x, torsions ) { out << x << " "; }
	out << endl;
}

void LogFilePerturber::log_angles(ostream & out, ParameterList angles) {
	out << "ANGLES ";
	foreach ( Real x, angles ) { out << x << " "; }
	out << endl;
}

void LogFilePerturber::log_lengths(ostream & out, ParameterList lengths) {
	out << "LENGTHS ";
	foreach ( Real x, lengths ) { out << x << " "; }
	out << endl;
}

}
}
}
