// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/NoOutputJobOutputter.hh
/// @brief  header file for NoOutputJobOutputter class, part of August 2008 job distributor as planned at RosettaCon08
/// @author Steven Lewis smlewi@gmail.com


#ifndef INCLUDED_protocols_jd2_NoOutputJobOutputter_hh
#define INCLUDED_protocols_jd2_NoOutputJobOutputter_hh

//unit headers
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

//utility headers
//C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {

/// @details the NoOutputJobOutputter class is responsible for NOT dealing with output
class NoOutputJobOutputter : public JobOutputter {
public:

	//constructor -- reads cmd-line to initialize evaluators
	NoOutputJobOutputter() {};

	//////////////////////////////creating output functions/////////////////////////////////////////

	/// @brief this function takes a string and writes it to disk (separately from Tracer output).
	///use some sort of extention option system - default .dat?  .data?
	virtual
	void file( JobCOP, std::string const & )  {};

	/// @brief this function outputs the final result of a job.
	virtual
	void final_pose( JobOP, core::pose::Pose const &, std::string const & ) {};

	/// @brief this function is intended for saving mid-protocol poses; for example the final centroid structure in a combined centroid/fullatom protocol.
	virtual
	void other_pose( JobOP, core::pose::Pose const & , std::string const &, int /*copy_count = -1*/, bool /*score_only = false*/  ) {};

	/////////////////////////////////state of output functions/////////////////////////////////

	/// @brief this function is not used for output, but it belongs here since it needs to check the same output locations as the class normally writes to.  This class checks wherever output goes to see if the job's expected output already exists (on disk or whatever).  This is the most basic form of checkpointing.
	virtual
	bool job_has_completed( JobCOP ) { return false; };

	virtual std::string output_name( JobCOP job ) {
		return affixed_numbered_name( job );
	}

}; // NoOutputJobOutputter

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_NoOutputJobOutputter_HH
