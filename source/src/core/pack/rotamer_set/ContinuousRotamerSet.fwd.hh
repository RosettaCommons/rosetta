// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/RotamerSet/ContinuousRotamerSet.fwd.hh
/// @brief  Forward declaration of the ContinuousRotamerSet class
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com)


#ifndef INCLUDED_core_pack_rotamer_set_ContinuousRotamerSet_FWD_HH
#define INCLUDED_core_pack_rotamer_set_ContinuousRotamerSet_FWD_HH

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace pack {
namespace rotamer_set {

class ContinuousRotamerSet;
class ContinuousRotamerSets;

typedef utility::pointer::shared_ptr< ContinuousRotamerSet > ContinuousRotamerSetOP;
typedef utility::pointer::shared_ptr< ContinuousRotamerSet const > ContinuousRotamerSetCOP;

typedef utility::pointer::shared_ptr< ContinuousRotamerSets > ContinuousRotamerSetsOP;
typedef utility::pointer::shared_ptr< ContinuousRotamerSets const > ContinuousRotamerSetsCOP;

} // namespace rotamer_set
} // namespace pack
} // namespace core


#endif
