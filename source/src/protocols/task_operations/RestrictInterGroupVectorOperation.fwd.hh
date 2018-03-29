// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/task_operations/RestrictInterGroupVectorOperation.fwd.hh
/// @brief
/// @author  Ben Stranges stranges@unc.edu

#ifndef INCLUDED_protocols_task_operations_RestrictInterGroupVectorOperation_fwd_hh
#define INCLUDED_protocols_task_operations_RestrictInterGroupVectorOperation_fwd_hh

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace task_operations {

class RestrictInterGroupVectorOperation;

typedef utility::pointer::shared_ptr< RestrictInterGroupVectorOperation > RestrictInterGroupVectorOperationOP;
typedef utility::pointer::shared_ptr< RestrictInterGroupVectorOperation const > RestrictInterGroupVectorOperationCOP;


} // task_operations
} // protocols


#endif
