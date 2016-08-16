// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/ResidueNeighborList.fwd.hh
/// @brief  Forward declaration of a container class for use by the Etable and FA_Elec
///         classes for storing residue-pair level atom-neighbor information
/// @author Andrew Leaver-Fay

#ifndef INCLUDED_core_scoring_ResidueNeighborList_fwd_hh
#define INCLUDED_core_scoring_ResidueNeighborList_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace scoring {

class ResidueNblistData;
typedef utility::pointer::shared_ptr< ResidueNblistData > ResidueNblistDataOP;
typedef utility::pointer::shared_ptr< ResidueNblistData const > ResidueNblistDataCOP;

class ResiduePairNeighborList;

typedef utility::pointer::shared_ptr< ResiduePairNeighborList > ResiduePairNeighborListOP;
typedef utility::pointer::shared_ptr< ResiduePairNeighborList const > ResiduePairNeighborListCOP;

}
}

#endif
