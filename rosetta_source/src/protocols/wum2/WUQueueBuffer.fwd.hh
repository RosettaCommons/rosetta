// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/wum2/WUQueueBuffer.fwd.hh
/// @brief  Foward decls for WUQueueBuffer, memory aware structure that links irecv and isend
/// 				with their mpi::request status, so buffers are maintained until messages are sent/recvd
/// @author Ken Jung

#ifndef INCLUDED_protocols_wum2_WUQueueBuffer_fwd_hh
#define INCLUDED_protocols_wum2_WUQueueBuffer_fwd_hh

#include <boost/shared_ptr.hpp>

namespace protocols {
namespace wum2 {

class WUQueueBuffer;
typedef boost::shared_ptr< WUQueueBuffer > WUQueueBufferSP;


} //namespace wum2
} //namespace protocols

#endif

