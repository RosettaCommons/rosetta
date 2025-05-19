// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/chemical/modifications/Reprotonate.fwd.hh
/// @brief  Heuristically reprotonate a molecule to be in neutral pH
/// @author Rocco Moretti (rmorettiase@gmail.com)


#ifndef INCLUDED_core_chemical_modifications_Reprotonate_fwd_hh
#define INCLUDED_core_chemical_modifications_Reprotonate_fwd_hh


#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace chemical {
namespace modifications {

class Reprotonate;

typedef utility::pointer::shared_ptr< Reprotonate > ReprotonateOP;


}
}
}
#endif
