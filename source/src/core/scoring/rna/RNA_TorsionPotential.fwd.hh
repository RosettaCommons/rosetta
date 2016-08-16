// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/RNA_TorsionPotential.fwd.hh
/// @brief  RNA_TorsionPotential potential class forward delcaration
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

#ifndef INCLUDED_core_scoring_rna_RNA_TorsionPotential_fwd_hh
#define INCLUDED_core_scoring_rna_RNA_TorsionPotential_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace rna {

class RNA_TorsionPotential;

typedef  utility::pointer::shared_ptr< RNA_TorsionPotential > RNA_TorsionPotentialOP;
typedef  utility::pointer::shared_ptr< RNA_TorsionPotential const > RNA_TorsionPotentialCOP;

} //rna
} //scoring
} //core

#endif
