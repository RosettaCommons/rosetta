// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/wum2/MPI_EndPoint.fwd.hh
/// @brief  Handles communication between different roles in mpi
// Every role in mpi needs one of these for each different role it communicates with
// for example, slave only communicates with master so slave needs one
// master communicates with slave and pool, so master needs two

// role_available_mem is a functor to the role's available memory function
//
// this allows a role with multiple endpoints have its endpoints be aware of 
// other endpoint memory usage and act accordingly
/// @author Ken Jung

#ifndef INCLUDED_protocols_wum2_MPI_EndPoint_fwd_hh
#define INCLUDED_protocols_wum2_MPI_EndPoint_fwd_hh

#include <boost/shared_ptr.hpp>

namespace protocols {
namespace wum2 {

class MPI_EndPoint;
typedef boost::shared_ptr< MPI_EndPoint > MPI_EndPointSP;


} //namespace wum2
} //namespace protocols

#endif

