// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/simple_moves/CartesianMD.fwd.hh
/// @brief Forward declarations for CartesianMD


#ifndef INCLUDED_protocols_simple_moves_CartesianMD_fwd_hh
#define INCLUDED_protocols_simple_moves_CartesianMD_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace md {

class CartesianMD;

typedef utility::pointer::shared_ptr< CartesianMD > CartesianMDOP;
typedef utility::pointer::shared_ptr< CartesianMD const > CartesianMDCOP;

} //movers
} //protocols

#endif /*INCLUDED_protocols_protein_interface_design_movers_CartesianMD_FWD_HH*/

