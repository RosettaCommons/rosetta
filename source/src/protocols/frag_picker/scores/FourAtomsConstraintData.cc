// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/scores/FourAtomsConstraintData.hh
/// @brief  provides a holder for data necessary to create a single constraint based on four atoms
/// @author Dominik Gront (dgront@chem.uw.edu.pl)

#include <protocols/frag_picker/scores/FourAtomsConstraintData.hh>


namespace protocols {
namespace frag_picker {
namespace scores {

FourAtomsConstraintData::~FourAtomsConstraintData() {}

/// @brief makes a new object
FourAtomsConstraintData::FourAtomsConstraintData(core::scoring::func::FuncOP function,
	core::Size first_atom, core::Size second_offset, core::Size second_atom,
	core::Size third_offset, core::Size third_atom, core::Size fourth_offset,
	core::Size fourth_atom) {
	func_ = function;
	first_atom_ = first_atom;
	second_atom_ = second_atom;
	third_atom_ = third_atom;
	fourth_atom_ = fourth_atom;
	second_offset_ = second_offset;
	third_offset_ = third_offset;
	fourth_offset_ = fourth_offset;
}

}
}
}
