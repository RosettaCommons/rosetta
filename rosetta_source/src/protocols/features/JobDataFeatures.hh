// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/protocols/features/JobDataFeatures.hh
/// @author Sam DeLuca

#ifndef INCLUDED_protocols_features_JobDataFeatures_hh_
#define INCLUDED_protocols_features_JobDataFeatures_hh_

//unit headers
#include <protocols/features/FeaturesReporter.hh>
#include <protocols/features/JobDataFeatures.fwd.hh>

//platform headers
#include <protocols/jd2/Job.fwd.hh>


namespace protocols {
namespace features {

class JobDataFeatures :public protocols::features::FeaturesReporter {
public:
	JobDataFeatures();
	JobDataFeatures(JobDataFeatures const & src);
	virtual ~JobDataFeatures();

	///@brief return string with class name
	std::string type_name() const;

	///@brief return sql statements that set up the right tables
	std::string schema() const;

	///@brief collect all the feature data for the pose
	core::Size
	report_features(
		core::pose::Pose const & pose,
		utility::vector1<bool> const & relevant_residues,
		core::Size struct_id,
		utility::sql_database::sessionOP db_session
	);

	void delete_record(
		core::Size struct_id,
		utility::sql_database::sessionOP db_session
	);
private:
	void insert_string_rows(core::Size struct_id, utility::sql_database::sessionOP db_session, protocols::jd2::JobCOP job) const;

	void insert_string_string_rows(core::Size struct_id, utility::sql_database::sessionOP db_session, protocols::jd2::JobCOP job) const;

	void insert_string_real_rows(core::Size struct_id, utility::sql_database::sessionOP db_session, protocols::jd2::JobCOP job) const;
};

}
}



#endif /* JOBDATAFEATURES_HH_ */
