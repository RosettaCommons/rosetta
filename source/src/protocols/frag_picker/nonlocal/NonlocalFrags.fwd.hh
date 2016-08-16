// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/frag_picker/nonlocal/NonlocalFrags.fwd.hh
/// @author David Kim (dekim@u.washington.edu)

#ifndef INCLUDED_PROTOCOLS_FRAG_PICKER_NONLOCAL_NONLOCALFRAGS_FWD_HH
#define INCLUDED_PROTOCOLS_FRAG_PICKER_NONLOCAL_NONLOCALFRAGS_FWD_HH

#include <utility/pointer/owning_ptr.fwd.hh>

namespace protocols {
namespace frag_picker {
namespace nonlocal {

class NonlocalFrags;
typedef utility::pointer::shared_ptr<NonlocalFrags> NonlocalFragsOP;
typedef utility::pointer::shared_ptr<NonlocalFrags const> NonlocalFragsCOP;

}  // namespace nonlocal
}  // namespace frag_picker
}  // namespace protocols

#endif  // PROTOCOLS_FRAG_PICKER_NONLOCAL_NONLOCALFRAGS_FWD_HH_
