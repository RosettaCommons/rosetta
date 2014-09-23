// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/BumpGrid.fwd.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_BumpGrid_fwd_hh
#define INCLUDED_protocols_match_BumpGrid_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace match {

class Bool3DGrid;
typedef utility::pointer::shared_ptr< Bool3DGrid > Bool3DGridOP;
typedef utility::pointer::shared_ptr< Bool3DGrid const > Bool3DGridCOP;

/// The different radii types used in collision detection.
/// These values are taken from Probe / Reduce from the Richardson lab.
/// Code inside BumpGrid.cc assumes the radii listed here being in non-decreasing order,
/// and that ZERO is the 0th element of this enumeration.
enum ProbeRadius
{
	ZERO = 0,   // radius for virtual atoms or atoms that should not be bump-checked.
	H_ARO,      // radius 1.0  -- Aromatic hydrogen
	H_ALA,      // radius 1.17 -- Alaphatic hydrogen
	OXY,        // radius 1.4  -- Oxygen
	NIT,        // radius 1.55 -- Nitrogen
	C_CAR,      // radius 1.65 -- Carbonyl carbon
	C_ALA,      // radius 1.75 -- Alaphatic and aromatic carbon
	SULPH,      // radius 1.8 -- Sulfur and Phosphorus

	n_probe_radii = SULPH // keep this guy last
};

class BumpGrid;

typedef utility::pointer::shared_ptr< BumpGrid > BumpGridOP;
typedef utility::pointer::shared_ptr< BumpGrid const > BumpGridCOP;

}
}

#endif
