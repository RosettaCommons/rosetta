// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/constraint_generator/MetalContactsConstraintGenerator.fwd.hh
/// @brief This constraint generator takes residue selectors for a residue/residues containing metal ion(s) and for residue(s) for which to set up contacts. It allows users to specify which base atoms will be used to define angles/dihedrals to constrain; ideal values for angles/dihedrals/distances; and an option to constrain to native values.
/// @author guffysl (guffy@email.unc.edu)

#ifndef INCLUDED_protocols_constraint_generator_MetalContactsConstraintGenerator_fwd_hh
#define INCLUDED_protocols_constraint_generator_MetalContactsConstraintGenerator_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace protocols {
namespace constraint_generator {

class MetalContactsConstraintGenerator;

typedef utility::pointer::shared_ptr< MetalContactsConstraintGenerator > MetalContactsConstraintGeneratorOP;
typedef utility::pointer::shared_ptr< MetalContactsConstraintGenerator const > MetalContactsConstraintGeneratorCOP;

} //protocols
} //constraint_generator

#endif //INCLUDED_protocols_constraint_generator_MetalContactsConstraintGenerator_fwd_hh
