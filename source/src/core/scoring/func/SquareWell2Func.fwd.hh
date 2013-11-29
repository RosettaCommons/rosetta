// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/func/SquareWell2Func.fwd.hh
/// @brief forward declaration for Mixture functions
/// @author Rhiju Das
/// @author Jianqing Xu


#ifndef INCLUDED_core_scoring_constraints_SquareWell2Func_fwd_hh
#define INCLUDED_core_scoring_constraints_SquareWell2Func_fwd_hh

#include <utility/pointer/owning_ptr.fwd.hh>

// C++ Headers

namespace core {
namespace scoring {
namespace constraints {

class SquareWell2Func;
typedef utility::pointer::owning_ptr< SquareWell2Func > SquareWell2FuncOP;
typedef utility::pointer::owning_ptr< SquareWell2Func const > SquareWell2FuncCOP;

//class CircularSquareWell2Func;
//typedef utility::pointer::owning_ptr< CircularSquareWell2Func > CircularSquareWell2FuncOP;
//typedef utility::pointer::owning_ptr< CircularSquareWell2Func const > CircularSquareWell2FuncCOP;

}
}
}

#endif
