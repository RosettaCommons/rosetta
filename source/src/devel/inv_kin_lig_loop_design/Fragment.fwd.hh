// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   devel/InvKinLigLoopDesign/Fragment.fwd.hh
///
/// @brief
/// @author

#ifndef DEVEL_INVKINLIGLOOPDESIGN_FRAGMENT_FWD_HH
#define DEVEL_INVKINLIGLOOPDESIGN_FRAGMENT_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace devel {
namespace inv_kin_lig_loop_design {

namespace Fragment {

class ResEntry;
class Entry;
class File;

typedef utility::pointer::shared_ptr< File > FileOP;
typedef utility::pointer::shared_ptr< File const > FileCOP;

} // namespace Fragment

struct Librarian;

} // namespace inv_kin_lig_loop_design
} // namespace devel

#endif
