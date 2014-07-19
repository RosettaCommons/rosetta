// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    protocols/wum2/EndPoint.cc
/// @brief   Non MPI version of EndPoint
/// @details This class is required because SingleNode role needs to use an EndPoint that is not MPI dependent (ie just a wrapper for 2 queues)
/// @author  Ken Jung

#include <protocols/wum2/EndPoint.hh>
#include <protocols/wum2/WorkUnit.hh>

namespace protocols {
namespace wum2 {

#ifndef __native_client__
EndPoint::EndPoint( boost::function < boost::uint64_t () > role_available_mem  ) : role_available_mem_( role_available_mem) {}
#endif

} // wum2
} // protocols
