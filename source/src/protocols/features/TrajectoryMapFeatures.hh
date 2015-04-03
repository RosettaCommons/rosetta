// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/TrajectoryMapFeatures.hh
/// @brief  Map trajectory structure_ids to cycle number
/// @details Will only work properly when specially created by TrajectoryReportToDB
/// @author Kyle Barlow (kb@kylebarlow.com)

#ifndef INCLUDED_protocols_features_TrajectoryMapFeatures_hh
#define INCLUDED_protocols_features_TrajectoryMapFeatures_hh

// Unit Headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/TrajectoryMapFeatures.fwd.hh>

//External

// Project Headers
#include <core/types.hh>
#include <utility/vector1.hh>

// C++ Headers
#include <string>

namespace protocols{
namespace features{

class TrajectoryMapFeatures : public protocols::features::FeaturesReporter {
public:
	TrajectoryMapFeatures();

	TrajectoryMapFeatures(TrajectoryMapFeatures const & src);

	virtual ~TrajectoryMapFeatures();

	/// @brief return string with class name
	std::string
	type_name() const;

	/// @brief generate the table schemas and write them to the database
	virtual void
	write_schema_to_db(utility::sql_database::sessionOP db_session) const;

	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const;

	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1< bool > const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session);

	void
	delete_record(
		StructureID struct_id,
		utility::sql_database::sessionOP db_sesion);

	void
	set_current_cycle(
		core::Size stride
	);

	core::Size
	get_current_cycle () const;

private:

	/// @brief Stores the cycle count for the current pose being saved
	/// @details Rather than storing the cycle count in pose (where it might be
	/// arbitrarily cleared), we are counting cycles/steps by mapping to output tag.
	/// Since the features reporter does not know what tag it is at apply time,
	/// but the ReportToDB subclassed mover does, this TrajectoryMapFeatures
	/// reporter is set up to query this variable to get the current cycle count,
	/// which should be up to date because the original apply() call is made to
	/// the TrajectoryMapFeatures mover.
	/// This variable being up-to-date and meaningful depends on it being set by
	/// the TrajectoryMapFeatures Mover (or an unimplemented at the time of this writing
	/// mover)
	core::Size current_cycle_;

};

} // namespace
} // namespace

#endif // include guard
