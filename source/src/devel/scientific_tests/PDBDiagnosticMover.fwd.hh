// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file devel/scientific_tests/PDBDiagnosticMover.fwd.hh
/// @brief does some very lightweight modeling on a PDB.  Meant to be run against the whole PDB as a scientific test
/// @author Steven Lewis (smlewi@gmail.com)

#ifndef INCLUDED_devel_scientific_tests_PDBDiagnosticMover_fwd_hh
#define INCLUDED_devel_scientific_tests_PDBDiagnosticMover_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace devel {
namespace scientific_tests {

class PDBDiagnosticMover;

typedef utility::pointer::shared_ptr< PDBDiagnosticMover > PDBDiagnosticMoverOP;
typedef utility::pointer::shared_ptr< PDBDiagnosticMover const > PDBDiagnosticMoverCOP;

} //devel
} //scientific_tests

#endif //INCLUDED_devel_scientific_tests_PDBDiagnosticMover_fwd_hh
