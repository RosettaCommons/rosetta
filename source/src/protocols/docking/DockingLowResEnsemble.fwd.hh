// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/docking/DockinLowResEnsemble.fwd.hh
/// @brief forward declarations for low-res ensemble docking
/// @author Daisuke Kuroda

#ifndef INCLUDED_protocols_docking_DockingLowRes_Ensemble_fwd_hh
#define INCLUDED_protocols_docking_DockingLowRes_Ensemble_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace docking {

/// @brief this mover does the low-resolution centroid mode phase of the EnsembleDock algorithm
class DockingLowResEnsemble;
typedef utility::pointer::shared_ptr< DockingLowResEnsemble > DockingLowResEnsembleOP;
typedef utility::pointer::shared_ptr< DockingLowResEnsemble const > DockingLowResEnsembleCOP;

}
}

#endif
