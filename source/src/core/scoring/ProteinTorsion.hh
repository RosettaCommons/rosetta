// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ProteinTorsion.hh
/// @brief  Ramachandran potential class delcaration
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_scoring_ProteinTorsion_hh
#define INCLUDED_core_scoring_ProteinTorsion_hh


namespace core {
namespace scoring {

enum ProteinTorsion {
	PHI = 1,
	PSI,
	OMEGA,
	CHI1,
	CHI2,
	CHI3,
	CHI4,
	protein_torsion_end
};


}
}

#endif
