// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/wum2/EndPoint.fwd.hh
/// @brief  Non MPI version of EndPoint
///  This class is required because SingleNode role needs to use an EndPoint that is not MPI dependent (ie just a wrapper for 2 queues)
/// @author Ken Jung

#ifndef INCLUDED_protocols_wum2_EndPoint_fwd_hh
#define INCLUDED_protocols_wum2_EndPoint_fwd_hh

#include <boost/shared_ptr.hpp>

namespace protocols {
namespace wum2 {

class EndPoint;
typedef boost::shared_ptr< EndPoint > EndPointSP;


} //namespace wum2
} //namespace protocols

#endif

