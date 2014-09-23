// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/OmegaTetherEnergy.fwd.hh
/// @brief  OmegaTether energy method class forward declaration
/// @author Phil Bradley
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_scoring_rna_RNA_SugarCloseEnergy_FWD_HH
#define INCLUDED_core_scoring_rna_RNA_SugarCloseEnergy_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {
namespace rna {

class RNA_SugarCloseEnergy;

typedef utility::pointer::shared_ptr< RNA_SugarCloseEnergy > RNA_SugarCloseEnergyOP;
typedef utility::pointer::shared_ptr< RNA_SugarCloseEnergy const > RNA_SugarCloseEnergyCOP;

} //rna
} //scoring
} //core


#endif
