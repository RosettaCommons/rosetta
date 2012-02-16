// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/features/FeaturesReporter.hh
/// @brief  Base class to report geometry and scores to features Statistics Scientific Benchmark
/// @author Matthew O'Meara

#ifndef INCLUDED_protocols_features_FeaturesReporter_hh
#define INCLUDED_protocols_features_FeaturesReporter_hh

// Unit Headers
#include <utility/pointer/ReferenceCount.hh>
#include <protocols/features/FeaturesReporter.fwd.hh>

// Project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <protocols/filters/Filter.fwd.hh>
#include <protocols/moves/Mover.fwd.hh>
#include <protocols/moves/DataMap.fwd.hh>
#include <utility/tag/Tag.fwd.hh>
#include <utility/vector1.fwd.hh>

// C++ Headers
#include <string>
#include <set>

#include <utility/vector1.hh>

#ifdef WIN32
#include <utility/sql_database/DatabaseSessionManager.hh>
#else
#include <utility/sql_database/DatabaseSessionManager.fwd.hh>
#endif


namespace protocols{
namespace features{

class FeaturesReporter : public utility::pointer::ReferenceCount {
public:

	virtual
	std::string
	type_name() const  {
		return "Unknown_FeaturesReporter";
	}

	///@brief return sql statements that sets up the appropriate tables
	///to contain the features.
	virtual
	std::string
	schema() const { return "";}

	///@brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	virtual
	utility::vector1<std::string>
	features_reporter_dependencies() const {
		utility::vector1<std::string> dependencies;
		return dependencies;
	}


	///@brief convience function to write the schema to a database.
	///There is no need to overload this in the derived class.
	void
	write_schema_to_db(
		utility::sql_database::sessionOP db_session
	) const;

	///@brief collect all the feature data for the pose.
	core::Size
	report_features(
		core::pose::Pose const & /*pose*/,
		core::Size /*parent_id*/,
		utility::sql_database::sessionOP /*db_session*/
	);

	///@brief collect all the feature data for the pose.
	virtual
	core::Size
	report_features(
		core::pose::Pose const & /*pose*/,
		utility::vector1< bool > const & /*relevant_residues*/,
		core::Size /*parent_id*/,
		utility::sql_database::sessionOP /*db_session*/
	);

	virtual
	void
	parse_my_tag(
		utility::tag::TagPtr const tag,
		protocols::moves::DataMap & /*data*/,
		protocols::filters::Filters_map const & /*filters*/,
		protocols::moves::Movers_map const & /*movers*/,
		core::pose::Pose const & /*pose*/);

	virtual
	void
	load_into_pose(
		utility::sql_database::sessionOP db_session,
		core::Size struct_id,
		core::pose::Pose & pose) {}

	virtual
	void
	delete_record(
		core::Size struct_id,
		utility::sql_database::sessionOP db_session) {}

protected:

	std::string
	find_tag(
		core::pose::Pose const & pose
	) const;

	///@brief a helper function for deleting data associated with a given structure from  feature database
	///WARNING table_name must be sanitized!
	void
	delete_records_from_table(
		std::string const & table_name,
		core::Size struct_id,
		utility::sql_database::sessionOP db_session);

};

} // namespace
} // namespace

#endif // include guard
