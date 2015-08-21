// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/TrajectoryReportToDB.fwd.hh
/// @brief  Report features data to database multiple times per structure, creating a trajectory
/// @author Kyle Barlow (kb@kylebarlow.com)

#ifndef INCLUDED_protocols_features_TrajectoryReportToDB_fwd_hh
#define INCLUDED_protocols_features_TrajectoryReportToDB_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace features {

class TrajectoryReportToDB;
typedef utility::pointer::shared_ptr< TrajectoryReportToDB > TrajectoryReportToDBOP;
typedef utility::pointer::shared_ptr< TrajectoryReportToDB const > TrajectoryReportToDBCOP;

}
}

#endif //INCLUDED_protocols_features_TrajectoryReportToDB_FWD_HH
