// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/util/ABEGOManager.fwd.hh
/// @brief foward header file for class of ABEGO
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

#ifndef INCLUDED_core_sequence_ABEGOManager_FWD_hh
#define INCLUDED_core_sequence_ABEGOManager_FWD_hh

// utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace sequence {

class ABEGOManager;

typedef utility::pointer::shared_ptr<ABEGOManager> ABEGOManagerOP;
typedef utility::pointer::shared_ptr<ABEGOManager const> ABEGOManagerCOP;

} // namespace sequence
} // namespace core


#endif /* INCLUDED_core_sequence_ABEGOManager_FWD_hh */

