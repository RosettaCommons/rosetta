// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

//////////////////////////////////////////////
///
/// @file basic/svd/SVD_Solver.fwd.hh
///
/// @brief SVD solver class
///
/// @authorv Christophe Schmitz & Srivatsan Raman
///
////////////////////////////////////////////////

#ifndef INCLUDED_basic_svd_SVD_Solver_fwd_hh
#define INCLUDED_basic_svd_SVD_Solver_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace basic {
namespace svd {

class SVD_Solver;
typedef utility::pointer::shared_ptr< SVD_Solver > SVD_SolverOP;
typedef utility::pointer::shared_ptr< SVD_Solver const > SVD_SolverCOP;

} //SVD
} // basic

#endif
