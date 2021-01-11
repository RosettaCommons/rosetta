// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file basic/thread_manager/RosettaThreadAllocation.cc
/// @author Jack Maguire, jackmaguire1444@gmail.com

#include <basic/thread_manager/RosettaThreadAllocation.hh>
#include <basic/thread_manager/RosettaThreadManager.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "basic.thread_manager.RosettaThreadAllocation" );

namespace basic {
namespace thread_manager {

RosettaThreadAllocation::~RosettaThreadAllocation(){
#ifndef NDEBUG
	TR.Debug << "RosettaThreadAllocation is destructing with " << thread_ids.size() << " elements" << std::endl;
#endif
	RosettaThreadManager::get_instance()->release_threads( * this );
}

} //thread_manager
} //basic

