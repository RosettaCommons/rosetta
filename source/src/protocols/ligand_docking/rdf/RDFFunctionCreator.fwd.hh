// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/ligand_docking/rdf/RDFFunctionCreator.fwd.hh
/// @brief  Forward Header for base class for RDFFunctionCreator
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_ligand_docking_rdf_RDFFunctionCreator_fwd_hh
#define INCLUDED_protocols_ligand_docking_rdf_RDFFunctionCreator_fwd_hh

// Utility Headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace ligand_docking {
namespace rdf {

/// @brief Abstract base class for a FeaturesReporter factory; the Creator class is responsible for
/// creating a particular FeaturesReporter class.
class RDFFunctionCreator;

typedef utility::pointer::owning_ptr< RDFFunctionCreator > RDFFunctionCreatorOP;
typedef utility::pointer::owning_ptr< RDFFunctionCreator const > RDFFunctionCreatorCOP;

}
} //namespace
} //namespace

#endif
