// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   src/protocols/features/JobDataFeatures.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_features_JobDataFeatures_hh
#define INCLUDED_protocols_features_JobDataFeatures_hh

//unit headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/JobDataFeatures.fwd.hh>

//External

//platform headers
#include <protocols/jd2/Job.fwd.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace features {

class JobDataFeatures :public protocols::features::FeaturesReporter {
public:
	JobDataFeatures();
	JobDataFeatures(JobDataFeatures const & src);
	~JobDataFeatures() override;

	/// @brief return string with class name
	std::string type_name() const override;

	/// @brief generate the table schemas and write them to the database
	void
	write_schema_to_db(utility::sql_database::sessionOP db_session) const override;

	/// @brief return the set of features reporters that are required to
	///also already be extracted by the time this one is used.
	utility::vector1<std::string>
	features_reporter_dependencies() const override;

	/// @brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1<bool> const & relevant_residues,
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	) override;

	void
	load_into_pose(
		utility::sql_database::sessionOP db_session,
		StructureID struct_id,
		core::pose::Pose & pose
	) override;

	void delete_record(
		StructureID struct_id,
		utility::sql_database::sessionOP db_session
	) override;
private:
	void insert_string_rows(StructureID struct_id, utility::sql_database::sessionOP db_session, protocols::jd2::JobCOP job) const;

	void insert_string_string_rows(StructureID struct_id, utility::sql_database::sessionOP db_session, protocols::jd2::JobCOP job) const;

	void insert_string_real_rows(StructureID struct_id, utility::sql_database::sessionOP db_session, protocols::jd2::JobCOP job) const;

	void load_string_data(utility::sql_database::sessionOP  db_session, StructureID struct_id, core::pose::Pose & pose);

	void load_string_string_data(utility::sql_database::sessionOP  db_session, StructureID struct_id, core::pose::Pose & pose);

	void load_string_real_data(utility::sql_database::sessionOP  db_session, StructureID struct_id, core::pose::Pose & pose);
};

}
}


#endif /* JOBDATAFEATURES_HH_ */
