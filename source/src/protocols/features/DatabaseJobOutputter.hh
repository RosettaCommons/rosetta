// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/features/DatabaseJobOutputter.hh
/// @brief  header file for DatabaseJobOutputter class
/// @author Matthew O'Meara (mattjomeara@gmail.com)


#ifndef INCLUDED_protocols_features_DatabaseJobOutputter_hh
#define INCLUDED_protocols_features_DatabaseJobOutputter_hh

// unit Headers
#include <protocols/features/DatabaseJobOutputter.fwd.hh>
#include <protocols/jd2/FileJobOutputter.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.fwd.hh>

// project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

// utility Headers
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>

#include <protocols/features/ProteinSilentReport.fwd.hh>


namespace protocols {
namespace features {

/// @details this is a implementation of JobOutputter for database-based output.
class DatabaseJobOutputter : public protocols::jd2::FileJobOutputter
{
public:

	typedef protocols::jd2::FileJobOutputter parent;

	DatabaseJobOutputter();

	~DatabaseJobOutputter() override;

	static void register_options();

	/// @brief load options from option sytem
	void
	load_options_from_option_system();

	/// @brief Set database name
	void
	set_database_name(std::string const & database_name);

	/// @brief Get database name
	std::string
	get_database_name() const;

	/// @brief Set database postgreSQL schema
	void
	set_database_pq_schema(std::string const & database_pq_schema);

	/// @brief Get database postgresQL schema
	std::string
	get_database_pq_schema() const;


	/// @brief see parent class for explanation
	void flush() override;

	/// @brief this function outputs the final result of a job.

	void final_pose( protocols::jd2::JobOP job, core::pose::Pose const & pose, std::string const & tag ) override;

	/// @brief this function is intended for saving mid-protocol poses;
	/// for example the final centroid structure in a combined
	/// centroid/fullatom protocol.

	void other_pose(
		protocols::jd2::JobOP job,
		core::pose::Pose const & pose,
		std::string const & tag,
		int copy_count = -1,
		bool score_only = false
	) override;

	/// @brief this function is not used for output, but it belongs here
	/// since it needs to check the same output locations as the class
	/// normally writes to.  This class checks wherever output goes to
	/// see if the job's expected output already exists (on disk or
	/// whatever).  This is the most basic form of checkpointing.

	bool job_has_completed( protocols::jd2::JobCOP job ) override;

public: // accessors

	/// @brief this is the master function for determining the
	/// unique output identifier for a job

	std::string output_name( protocols::jd2::JobCOP job ) override;

private: // members

	protocols::features::ProteinSilentReportOP protein_silent_report_;
	std::string database_name_;
	std::string database_pq_schema_;
	std::string path_;

}; // DatabaseJobOutputter

} // namespace features
} // namespace protocols

#endif //INCLUDED_protocols_features_DatabaseJobOutputter_hh
