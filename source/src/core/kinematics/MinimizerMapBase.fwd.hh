// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/kinematics/MinimizerMapBase.fwd.hh
/// @brief  Forward declaration of the minimizer-map base class for communicating between
///         an AtomTree and the class that will use the AtomTree to perform some gradient-based
///         minimization.
/// @author Phil Bradley


#ifndef INCLUDED_core_kinematics_MinimizerMapBase_fwd_hh
#define INCLUDED_core_kinematics_MinimizerMapBase_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace core {
namespace kinematics {


class MinimizerMapBase;

typedef utility::pointer::shared_ptr< MinimizerMapBase > MinimizerMapBaseOP;
typedef utility::pointer::shared_ptr< MinimizerMapBase const > MinimizerMapBaseCOP;


} // namespace kinematics
} // namespace core


#endif
