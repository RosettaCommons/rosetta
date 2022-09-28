// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/rdf/RDFBase.fwd.hh
///
/// @brief Forward headers for ChemistryBase
/// @author Steven Combs


#ifndef INCLUDED_core_chemical_modifications_ChemistryBase_fwd_hh
#define INCLUDED_core_chemical_modifications_ChemistryBase_fwd_hh


#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace chemical {
namespace modifications {

class ChemistryBase;

typedef utility::pointer::shared_ptr<ChemistryBase> ChemistryBaseOP;
typedef utility::pointer::shared_ptr<ChemistryBase const> ChemistryBaseCOP;


}
}
}



#endif
