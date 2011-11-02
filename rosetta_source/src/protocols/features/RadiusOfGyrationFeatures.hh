// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/RadiusOfGyrationFeatures.hh
/// @brief  report the radius of gyration to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_RadiusOfGyrationFeatures_hh
#define INCLUDED_protocols_features_RadiusOfGyrationFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/RadiusOfGyrationFeatures.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// Utility Headers
// AUTO-REMOVED #include <utility/sql_database/DatabaseSessionManager.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


namespace protocols{
namespace features{

class RadiusOfGyrationFeatures : public protocols::features::FeaturesReporter {
public:
	RadiusOfGyrationFeatures();

	RadiusOfGyrationFeatures( RadiusOfGyrationFeatures const & src );

	virtual ~RadiusOfGyrationFeatures();

	///@brief return sql statements that setup the right tables
	std::string
	schema() const;

	///@brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		core::Size stuct_id,
		utility::sql_database::sessionOP db_interface
	);

};

} // features namespace
} // protocols namespace

#endif // include guard
