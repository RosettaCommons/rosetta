// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

// Unit Headers
#include <core/scoring/constraints/Constraint.hh>

// Project Headers

#include <core/id/AtomID.hh>

// Utility header
#include <utility/vector1.hh>

// C++ Headers

#include <algorithm>

namespace core {
namespace scoring {
namespace constraints {

/// @details Auto-generated virtual destructor
Constraint::~Constraint() {}

Size
Constraint::show_violations(
	std::ostream & out,
	pose::Pose const &,
	Size,
	Real threshold
) const {
	out << "Constraint_show_violation stubbed out!\n" ;
	threshold = 1; //to make compile happy
	return 0;
}

utility::vector1< core::Size >
Constraint::residues() const {
	utility::vector1< int > pos_list;
	for ( Size i=1; i<= natoms(); ++i ) {
		int const seqpos( atom(i).rsd() );
		// seqpos already in list?
		if ( std::find( pos_list.begin(), pos_list.end(), seqpos ) == pos_list.end() ) {
			pos_list.push_back( seqpos );
		}
	}
	return pos_list;
}

}
}
}

