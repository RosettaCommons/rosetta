// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/docking/DockingInitialPerturbation.fwd.hh
/// @brief
/// @author Evan Baugh (ebaugh1@jhu.edu)

#ifndef INCLUDED_protocols_docking_DockingInitialPerturbation_fwd_hh
#define INCLUDED_protocols_docking_DockingInitialPerturbation_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace docking {

class DockingInitialPerturbation;
typedef utility::pointer::shared_ptr< DockingInitialPerturbation > DockingInitialPerturbationOP;
typedef utility::pointer::shared_ptr< DockingInitialPerturbation const > DockingInitialPerturbationCOP;

class DockingSlideIntoContact; // fwd declaration
typedef utility::pointer::shared_ptr< DockingSlideIntoContact > DockingSlideIntoContactOP;
typedef utility::pointer::shared_ptr< DockingSlideIntoContact const > DockingSlideIntoContactCOP;

class FaDockingSlideIntoContact; // fwd declaration
typedef utility::pointer::shared_ptr< FaDockingSlideIntoContact > FaDockingSlideIntoContactOP;
typedef utility::pointer::shared_ptr< FaDockingSlideIntoContact const > FaDockingSlideIntoContactCOP;

class SlideIntoContact; // fwd declaration
typedef utility::pointer::shared_ptr< SlideIntoContact > SlideIntoContactOP;
typedef utility::pointer::shared_ptr< SlideIntoContact const > SlideIntoContactCOP;

} // docking
} // protocols

#endif
