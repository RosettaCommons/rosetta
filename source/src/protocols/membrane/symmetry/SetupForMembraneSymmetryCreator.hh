// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file		protocols/membrane/symmetry/SetupForMembraneSymmetryCreator.hh
///
/// @brief		Setup a Symmetric Membrane Protein Using the Membrane Framework
/// @details	The setup for membrane symmetry class first adds a membrane residue
///				to the asymmetric unit of a protein, creates a symmetric complex, and
///				then adds the remainder of the membrane framework to capture the entire 
///				symmetric system. This should work for both relax and docking protocols
///				currently in Rosetta. 
///
///				Last Modified: 10/3/14
///
/// @author		Rebecca Alford (rfalford12@gmail.com)

#ifndef INCLUDED_protocols_membrane_symmetry_SetupForMembraneSymmetryCreator_hh
#define INCLUDED_protocols_membrane_symmetry_SetupForMembraneSymmetryCreator_hh

// Utility headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace membrane {
namespace symmetry {

/// @brief Mover Creator
class SetupForMembraneSymmetryCreator : public protocols::moves::MoverCreator {
	
public:
	
	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();
	
};

} // symmetry
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_symmetry_SetupForMembraneSymmetryCreator_hh
