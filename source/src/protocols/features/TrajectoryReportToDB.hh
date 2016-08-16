// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/TrajectoryReportToDB.hh
/// @brief  Report features data to database multiple times per structure, creating a trajectory
/// @author Kyle Barlow (kb@kylebarlow.com)

#ifndef INCLUDED_protocols_features_TrajectoryReportToDB_hh
#define INCLUDED_protocols_features_TrajectoryReportToDB_hh

// Unit Headers
#include <protocols/features/TrajectoryReportToDB.fwd.hh>
#include <protocols/features/ReportToDB.hh>

#include <protocols/features/TrajectoryMapFeatures.fwd.hh>

// // Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

// // Utility Headers
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.hh>

// // Boost Headers

// // C++ Headers
// #include <string>


namespace protocols {
namespace features {

class TrajectoryReportToDB : public ReportToDB {

public:
	TrajectoryReportToDB();

	TrajectoryReportToDB(
		core::Size stride
	);

	TrajectoryReportToDB(
		std::string const & name
	);

	TrajectoryReportToDB(
		utility::sql_database::sessionOP db_session,
		std::string const & batch_name,
		std::string const & batch_description,
		bool use_transactions=true,
		core::Size cache_size=2000
	);

	TrajectoryReportToDB(TrajectoryReportToDB const & src);

	virtual ~TrajectoryReportToDB();

	virtual moves::MoverOP fresh_instance() const;

	virtual moves::MoverOP clone() const;

	void
	parse_my_tag(
		utility::tag::TagCOP tag,
		basic::datacache::DataMap & data,
		Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & pose );

	void
	parse_stride_tag_item(
		utility::tag::TagCOP tag);

	void
	apply(
		Pose& pose
	);

	// Getters/Setters
	virtual std::string name() { return "TrajectoryReportToDB"; }
	virtual std::string get_name() const { return "TrajectoryReportToDB"; }

	void
	set_stride(
		core::Size stride
	);

	core::Size
	get_stride () const;

	/// @details Probably only needed for unit testing
	std::map<std::string, core::Size>
	get_cycle_counts() const;

private:
	/// @brief Add trajectory_map_features_reporter_ to list of features reporters
	void initialize_trajectory_reporter();

	/// @brief Stores the number of cycles each output tag has been applied
	std::map<std::string, core::Size> cycle_counts_;

	/// @brief Features will be reported this number of cycles
	core::Size stride_;

	/// @brief Special pointer to TrajectoryMapFeatures reporter that we add
	/// so that we can update the current_cycle variable in that reporter
	TrajectoryMapFeaturesOP trajectory_map_features_reporter_;
};

} // namespace
} // namespace

#endif //include guard
