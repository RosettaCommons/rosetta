// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/thread_manager/RosettaThreadAllocation.hh
/// @author Jack Maguire

#ifndef INCLUDED_basic_thread_manager_RosettaThreadAllocation_hh
#define INCLUDED_basic_thread_manager_RosettaThreadAllocation_hh

#include <basic/thread_manager/RosettaThreadAllocation.fwd.hh>

// Platform headers
#include <platform/types.hh>
#include <utility/vector1.hh>

namespace basic {
namespace thread_manager {

struct RosettaThreadAllocation {

	RosettaThreadAllocation() = default;
	RosettaThreadAllocation( RosettaThreadAllocation && ) = default;

	RosettaThreadAllocation( RosettaThreadAllocation const & ) = delete;

	~RosettaThreadAllocation();

	//threads_.size() will be 1 thread smaller than you expect
	//  because the "master" thread will also run
	utility::vector1< platform::Size > thread_ids;
};


} //thread_manager
} //basic

#endif //INCLUDED_basic_thread_manager_RosettaThreadAllocation_hh
