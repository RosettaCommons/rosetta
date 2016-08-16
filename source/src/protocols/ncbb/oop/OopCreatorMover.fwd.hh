// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/ncbb/oop/OopCreatorMover.fwd.hh
///
/// @brief
/// @author Andy Watkins


#ifndef INCLUDED_protocols_ncbb_oop_OopCreatorMover_fwd_hh
#define INCLUDED_protocols_ncbb_oop_OopCreatorMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace ncbb {
namespace oop {


class OopCreatorMover; // fwd declaration
typedef utility::pointer::shared_ptr< OopCreatorMover > OopCreatorMoverOP;
typedef utility::pointer::shared_ptr< OopCreatorMover const > OopCreatorMoverCOP;


} // namespace oop
} // namespace ncbb
} // namespace protocols

#endif // INCLUDED_protocols_ncbb_oop_OopCreatorMover_FWD_HH
