// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/generalized_kinematic_closure/perturber/GeneralizedKICperturber.fwd.hh
/// @brief  Defines owning pointers for GeneralizedKICperturber class.
/// @author Vikram K. Mulligan (vmullig@uw.edu)

#ifndef INCLUDED_protocols_generalized_kinematic_closure_perturber_GeneralizedKICperturber_fwd_hh
#define INCLUDED_protocols_generalized_kinematic_closure_perturber_GeneralizedKICperturber_fwd_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace generalized_kinematic_closure {
namespace perturber {

class GeneralizedKICperturber; // fwd declaration
typedef utility::pointer::shared_ptr< GeneralizedKICperturber > GeneralizedKICperturberOP;
typedef utility::pointer::shared_ptr< GeneralizedKICperturber const > GeneralizedKICperturberCOP;
typedef utility::vector1<GeneralizedKICperturberOP> GeneralizedKICperturberOPs;
typedef utility::vector1<GeneralizedKICperturberCOP> GeneralizedKICperturberCOPs;

} // namespare perturber
} // namespace generalized_kinematic_closure
} // namespace protocols

#endif // INCLUDED_protocols_generalized_kinematic_closure_perturber_GeneralizedKICperturber_fwd_hh
