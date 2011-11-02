// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/DatabaseJobOutputter.hh
/// @brief  header file for DatabaseJobOutputter class
/// @author Matthew O'Meara (mattjomeara@gmail.com)


#ifndef INCLUDED_protocols_jd2_DatabaseJobOutputter_hh
#define INCLUDED_protocols_jd2_DatabaseJobOutputter_hh

// unit Headers
#include <protocols/jd2/DatabaseJobOutputter.fwd.hh>
#include <protocols/jd2/FileJobOutputter.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.fwd.hh>

// project Headers
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <protocols/features/ProteinSilentReport.hh>

// utility Headers
#include <utility/file/FileName.hh>
#include <utility/vector1.hh>

// C++ headers
#include <string>

#include <protocols/features/ProteinSilentReport.fwd.hh>


namespace protocols {
namespace jd2 {

///@details this is a implementation of JobOutputter for database-based output.
class DatabaseJobOutputter : public protocols::jd2::FileJobOutputter
{
public:

  typedef protocols::jd2::FileJobOutputter parent;

  DatabaseJobOutputter();
  virtual ~DatabaseJobOutputter();

	static void register_options();

	/// @brief load options from option sytem
	void
	load_options_from_option_system();

	/// @brief Set database file name
	void
	set_database_fname(std::string const & database_fname);

	/// @brief Get database file name
	std::string
	get_database_fname() const;

	///@brief see parent class for explanation
	virtual void flush();

  ///@brief this function outputs the final result of a job.
  virtual
  void final_pose( JobCOP job, core::pose::Pose const & pose );

	/// @brief this function is intended for saving mid-protocol poses;
	/// for example the final centroid structure in a combined
	/// centroid/fullatom protocol.
  virtual
  void other_pose(
		JobCOP job,
		core::pose::Pose const & pose,
		std::string const & tag );

	/// @brief this function is not used for output, but it belongs here
	/// since it needs to check the same output locations as the class
	/// normally writes to.  This class checks wherever output goes to
	/// see if the job's expected output already exists (on disk or
	/// whatever).  This is the most basic form of checkpointing.
  virtual
  bool job_has_completed( JobCOP job );

public: // accessors

	/// @brief this is the master function for determining the
	/// unique output identifier for a job
  virtual
  std::string output_name( JobCOP job );

private: // members

	protocols::features::ProteinSilentReportOP protein_silent_report_;
	std::string database_fname_;

}; // DatabaseJobOutputter

} // namespace
} // namespace

#endif //include guard
