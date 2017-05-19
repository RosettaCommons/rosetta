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

#include <vector>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <math.h>
#include <assert.h>
#ifndef PYROSETTA
#include <time.h>
#endif
#ifndef __WIN32__
#include <sys/resource.h>
#endif

namespace protocols {
namespace cluster {
namespace calibur {

//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

#ifndef __WIN32__
#ifndef PYROSETTA
double
__timeval_difference(struct timeval * x, struct timeval * y);

double
_get_elapsed(int set_start);
#endif
#endif

}
}
}

#endif
