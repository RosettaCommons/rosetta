// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/disulfides/DisulfideAtomIndices.fwd.hh
/// @brief  Disulfide Atom Indices class forward declaration
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_disulfides_DisulfideAtomIndices_fwd_hh
#define INCLUDED_core_scoring_disulfides_DisulfideAtomIndices_fwd_hh

namespace core {
namespace scoring {
namespace disulfides {


enum DisulfideDerivativeAtom {
	NO_DERIVATIVES_FOR_ATOM = 0,
	CYS_C_ALPHA,
	CYS_C_BETA,
	CYS_S_GAMMA,
	CYS_S_DELTA,
	CYS_CEN
};


class DisulfideAtomIndices;

}
}
}

#endif
