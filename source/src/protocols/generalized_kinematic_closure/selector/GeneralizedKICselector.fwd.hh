// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/generalized_kinematic_closure/selector/GeneralizedKICselector.fwd.hh
/// @brief  Defines owning pointers for GeneralizedKICselector class.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_generalized_kinematic_closure_selector_GeneralizedKICselector_fwd_hh
#define INCLUDED_protocols_generalized_kinematic_closure_selector_GeneralizedKICselector_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace generalized_kinematic_closure {
namespace selector {

class GeneralizedKICselector; // fwd declaration
typedef utility::pointer::shared_ptr< GeneralizedKICselector > GeneralizedKICselectorOP;
typedef utility::pointer::shared_ptr< GeneralizedKICselector const > GeneralizedKICselectorCOP;
typedef utility::vector1<GeneralizedKICselectorOP> GeneralizedKICselectorOPs;
typedef utility::vector1<GeneralizedKICselectorCOP> GeneralizedKICselectorCOPs;

} // namespare selector
} // namespace generalized_kinematic_closure
} // namespace protocols

#endif // INCLUDED_protocols_generalized_kinematic_closure_selector_GeneralizedKICselector_fwd_hh
