// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/ReportToDB.fwd.hh
/// @brief  report data to sqlite database
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_ReportToDB_fwd_hh
#define INCLUDED_protocols_features_ReportToDB_fwd_hh

// Utility headers
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace features {

class ReportToDB;
typedef utility::pointer::shared_ptr< ReportToDB > ReportToDBOP;
typedef utility::pointer::shared_ptr< ReportToDB const > ReportToDBCOP;

}
}

#endif //INCLUDED_protocols_features_ReportToDB_FWD_HH
