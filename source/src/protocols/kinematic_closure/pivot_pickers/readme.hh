// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// #error Don't include; for documentation only.

/// @namespace protocols::kinematic_closure::pivot_pickers
///
/// @brief Algorithms for choosing pivot residues.
///
/// @details The pivot torsions are used by the kinematic closure algorithm to
/// ensure that sampling doesn't cause the breaks in the backbone.  The choice
/// of pivots can be significant, because the pivot torsions are the most
/// difficult degrees of freedom to control.  In other words, since these DOFs
/// are used to ensure that the backbone is closed, it is hard to impose other
/// restraints on them.  This is in contrast with the pivot torsions, which can
/// be set arbitrarily.  So if you're studying a loop and you know that some of
/// the residues have really restrained torsions, it would be better to not
/// represent these residues as pivots.  This level of customization may
/// require you to write a custom PivotPicker subclass, but it's not too hard.
