// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/drug_design/substituion_support.fwd.hh
/// @brief use RDKit to substitute items based on matched templates
/// @author Rocco Moretti (rmorettiase@gmail.com)

#ifndef INCLUDED_protocols_drug_design_substitution_support_fwd_hh
#define INCLUDED_protocols_drug_design_substitution_support_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace drug_design {

class AtomSubstitution;
typedef utility::pointer::shared_ptr< AtomSubstitution > AtomSubstitutionOP;
typedef utility::pointer::shared_ptr< AtomSubstitution const > AtomSubstitutionCOP;

class MoleculeSubstitution;
typedef utility::pointer::shared_ptr< MoleculeSubstitution > MoleculeSubstitutionOP;
typedef utility::pointer::shared_ptr< MoleculeSubstitution const > MoleculeSubstitutionCOP;

}
}

#endif // Include guard
