// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file external/calibur/SimPDB.fwd.hh
/// @author SC Li & YK Ng (kalngyk@gmail.com)

#ifndef external_calibur_SimPDB_FWD_HH
#define external_calibur_SimPDB_FWD_HH

#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace cluster {
namespace calibur {

class SimPDB;
typedef utility::pointer::shared_ptr< SimPDB > SimPDBOP;


}
}
}

#endif
