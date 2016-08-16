// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/antibody/grafting/scs_functor.fwd.hh
/// @brief Structural Predicates for Structural Component Selector (SCS) layer
/// @author Sergey Lyskov


#ifdef CXX11

#ifndef INCLUDED_protocols_antibody_grafting_scs_functor_fwd_hh
#define INCLUDED_protocols_antibody_grafting_scs_functor_fwd_hh


#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace antibody {
namespace grafting {

class SCS_Functor;
typedef utility::pointer::shared_ptr< SCS_Functor > SCS_FunctorOP;
typedef utility::pointer::shared_ptr< SCS_Functor const > SCS_FunctorCOP;

} // namespace grafting
} // namespace antibody
} // namespace protocols

#endif // INCLUDED_protocols_antibody_grafting_scs_functor_fwd_hh

#endif // CXX11
