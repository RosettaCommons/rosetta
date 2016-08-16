// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/conformation/RotamerSetBase.hh
/// @brief  Abstract base class for a class that holds disembodied Residues (ie not contained within a Pose) as Rotamers
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)


#ifndef INCLUDED_core_conformation_RotamerSetBase_fwd_hh
#define INCLUDED_core_conformation_RotamerSetBase_fwd_hh


// Utility headers
#include <utility/pointer/owning_ptr.hh>

//Auto Headers
namespace core {
namespace conformation {

class RotamerSetBase;

typedef utility::pointer::shared_ptr< RotamerSetBase > RotamerSetBaseOP;
typedef utility::pointer::shared_ptr< RotamerSetBase const > RotamerSetBaseCOP;

} // namespace conformation
} // namespace core


#endif
