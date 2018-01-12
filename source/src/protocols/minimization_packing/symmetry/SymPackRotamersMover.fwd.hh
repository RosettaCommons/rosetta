// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Ingemar Andre

#ifndef INCLUDED_protocols_minimization_packing_symmetry_SymPackRotamersMover_fwd_hh
#define INCLUDED_protocols_minimization_packing_symmetry_SymPackRotamersMover_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace minimization_packing {
namespace symmetry {

class SymPackRotamersMover;
typedef utility::pointer::shared_ptr< SymPackRotamersMover > SymPackRotamersMoverOP;
typedef utility::pointer::shared_ptr< SymPackRotamersMover const > SymPackRotamersMoverCOP;

} // symmetry
} // moves
} // protocols

#endif //INCLUDED_protocols_minimization_packing_symmetry_SymPackRotamersMover_fwd
