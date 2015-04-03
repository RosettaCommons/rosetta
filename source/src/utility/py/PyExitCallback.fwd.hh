// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/utility/PyExitCallback.fwd.hh
/// @brief  Class to specify utility_exit callback in PyRosetta.
/// @author Sergey Lyskov


#ifndef INCLUDED_utility_PyExitCallback_fwd_hh
#define INCLUDED_utility_PyExitCallback_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace utility { namespace py {

class PyExitCallback;

typedef utility::pointer::shared_ptr< PyExitCallback > PyExitCallbackOP;
typedef utility::pointer::shared_ptr< PyExitCallback const > PyExitCallbackCOP;

} } // namespace py namespace utility

#endif // INCLUDED_utility_PyExitCallback_fwd_hh
