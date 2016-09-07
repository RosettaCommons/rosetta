// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/Report.hh
/// @brief  report feature data to database
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_Report_hh
#define INCLUDED_protocols_features_Report_hh

// Unit Headers
#include <protocols/features/Report.fwd.hh>

// Package Headers
#ifdef WIN32
#include <protocols/features/FeaturesReporter.hh>
#endif

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>

//Auto Headers
#include <protocols/features/FeaturesReporter.fwd.hh>

namespace protocols {
namespace features {

class Report : public utility::pointer::ReferenceCount {

public:
	typedef utility::vector1< FeaturesReporterOP > FeaturesReporters;

public:
	Report();

	Report(Report const & src);

	~Report() override;

	virtual
	core::Size
	version() = 0;

	virtual
	void
	apply(
		core::pose::Pose const & pose,
		utility::sql_database::sessionOP db_sesion) = 0;


protected:

	FeaturesReporters reporters_;

};

}//features
}//protocols

#endif
