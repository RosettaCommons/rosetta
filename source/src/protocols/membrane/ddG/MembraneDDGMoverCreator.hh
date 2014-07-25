// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file       protocols/membrane/ddG/MembraneDDGMoverCreator.hh
///
/// @brief      Compute ddG Scores for Membrane Protein (Rosetta Scripts Creator Hook)
/// @details	Initialize a membrane pose, compute an initial membrane position,
///				compute a native score, make mutation, repack sidechains, and score new
///				structure. Uses a user-provided resfile for repacking/mutations
///
/// @author     Rebecca Alford (rfalford12@gmail.com)
/// @note       Last Modified (7/7/14)

#ifndef INCLUDED_protocols_membrane_ddG_MembraneDDGMoverCreator_hh
#define INCLUDED_protocols_membrane_ddG_MembraneDDGMoverCreator_hh

// Utility headers
#include <protocols/moves/MoverCreator.hh>

namespace protocols {
namespace membrane {
namespace ddG {

/// @brief Mover Creator
class MembraneDDGMoverCreator : public protocols::moves::MoverCreator {
	
public:
	
	virtual protocols::moves::MoverOP create_mover() const;
	virtual std::string keyname() const;
	static std::string mover_name();
	
};

} // ddG
} // membrane
} // protocols

#endif // INCLUDED_protocols_membrane_ddG_MembraneDDGMoverCreator_hh
