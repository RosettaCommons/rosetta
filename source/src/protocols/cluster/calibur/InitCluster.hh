// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/kalngyk/InitCluster.hh
/// @author SC Li & YK Ng (kalngyk@gmail.com)

#ifndef external_calibur_InitCluster_HH
#define external_calibur_InitCluster_HH

#ifndef PYROSETTA
#include <time.h> // DO NOT AUTO-REMOVE
#endif
#if !defined(__WIN32__) && !defined(WIN32)
#include <sys/resource.h> // DO NOT AUTO-REMOVE
#endif
#include <core/types.hh>

namespace protocols {
namespace cluster {
namespace calibur {

//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

#if !defined(__WIN32__) && !defined(WIN32)
#ifndef PYROSETTA
core::Real
__timeval_difference(struct timeval * x, struct timeval * y);

core::Real
_get_elapsed(int set_start);
#endif
#endif

}
}
}

#endif
