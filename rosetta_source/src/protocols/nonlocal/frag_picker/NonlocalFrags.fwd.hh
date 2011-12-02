// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/nonlocal/frag_picker/NonlocalFrags.fwd.hh
/// @author David Kim (dekim@u.washington.edu)

#ifndef PROTOCOLS_NONLOCAL_FRAG_PICKER_NONLOCALFRAGS_FWD_HH_
#define PROTOCOLS_NONLOCAL_FRAG_PICKER_NONLOCALFRAGS_FWD_HH_

#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace nonlocal {
namespace frag_picker {

class NonlocalFrags;
typedef utility::pointer::owning_ptr<NonlocalFrags> NonlocalFragsOP;
typedef utility::pointer::owning_ptr<NonlocalFrags const> NonlocalFragsCOP;

}  // namespace nonlocal
}  // namespace frag_picker
}  // namespace protocols

#endif  // PROTOCOLS_FRAG_PICKER_NONLOCAL_NONLOCALFRAGS_FWD_HH_
