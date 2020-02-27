// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/kinematics/jacobian/ModuleType1.fwd.hh
/// @brief Class that contains the linearized Cartesian motion/force space of a 1,2,3-amino acid module
/// @author teunhoevenaars (teunhoevenaars@gmail.com)

#ifndef INCLUDED_core_kinematics_jacobian_ModuleType1_fwd_hh
#define INCLUDED_core_kinematics_jacobian_ModuleType1_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>


// Forward
namespace core {
namespace kinematics {
namespace jacobian {

class ModuleType1;

using ModuleType1OP = utility::pointer::shared_ptr< ModuleType1 >;
using ModuleType1COP = utility::pointer::shared_ptr< ModuleType1 const >;

} //jacobian
} //kinematics
} //core

#endif //INCLUDED_core_kinematics_jacobian_ModuleType1_fwd_hh
