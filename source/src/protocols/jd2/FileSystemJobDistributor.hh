// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/FileSystemJobDistributor.hh
/// @brief  August 2008 job distributor as planned at RosettaCon08 - Simplest class FileSystemJobDistributor
/// @author Andrew Leaver-Fay
/// @author Steven Lewis smlewi@gmail.com

#ifndef INCLUDED_protocols_jd2_FileSystemJobDistributor_hh
#define INCLUDED_protocols_jd2_FileSystemJobDistributor_hh

// Unit headers
#include <protocols/jd2/FileSystemJobDistributor.fwd.hh>

// Package headers
#include <protocols/jd2/JobDistributor.hh>
//#include <protocols/jd2/JobInputter.fwd.hh>
//#include <protocols/jd2/JobOutputter.fwd.hh>
#include <protocols/jd2/Job.fwd.hh>
// AUTO-REMOVED #include <protocols/jd2/JobDistributorFactory.fwd.hh>

//#include <protocols/moves/Mover.fwd.hh>

// Utility headers
#include <core/types.hh>

// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {

class FileSystemJobDistributor : public JobDistributor
{
protected:
	FileSystemJobDistributor();

public:
	///WARNING WARNING!  SINGLETONS' DESTRUCTORS ARE NEVER CALLED IN MINI!  DO NOT TRY TO PUT THINGS IN THIS FUNCTION!
	///here's a nice link explaining why: http://www.research.ibm.com/designpatterns/pubs/ph-jun96.txt
	virtual ~FileSystemJobDistributor();

	virtual void restart();

	virtual
	core::Size
	get_new_job_id();

	virtual
	void
	mark_current_job_id_for_repetition();

	virtual
	void
	remove_bad_inputs_from_job_list();

	virtual
	void
	current_job_finished();

	virtual
	void
	go( protocols::moves::MoverOP mover );

	friend class JobDistributorFactory; // calls protected ctor
protected:

	///@brief This function is called upon a successful job completion; it has been virtualized so BOINC and MPI can delay/protect output
	///base implementation is just a call to the job outputter
	virtual
	void
	job_succeeded( core::pose::Pose & pose, core::Real run_time, std::string const & tag );

	///@brief This function is called when we five up on the job;  it has been virtualized so BOINC and MPI can delay/protect output
	///base implementation is just a call to the job outputter
	virtual
	void
	job_failed( core::pose::Pose & pose, bool will_retry );

	virtual void handle_interrupt();

private:
	/// @brief delete all temporary files for this job
	void delete_in_progress_files();

	///@brief tacks an extension (.in_progress) on to filename - for multiple process/one directory behaviour
	std::string const temporary_file_name( JobCOP job ) const;

	std::string extension_;
    std::string path_;

	//	core::Size next_job_to_try_assigning_; //don't it is redundant with current_job_id_ in the base-class

	core::Size retry_count_;
};

}//jd2
}//protocols

#endif //INCLUDED_protocols_jd2_FileSystemJobDistributor_HH
