// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/NoOutputJobOutputter.hh
/// @brief  header file for ScoreOnlyJobOutputter class, allows the user to output only score files
/// @author Steven Lewis smlewi@gmail.com


#ifndef INCLUDED_protocols_jd2_ScoreOnlyJobOutputter_hh
#define INCLUDED_protocols_jd2_ScoreOnlyJobOutputter_hh

//unit headers
#include <protocols/jd2/FileJobOutputter.hh>
#include <protocols/jd2/ScoreOnlyJobOutputter.fwd.hh>
#include <protocols/jd2/Job.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

//utility headers
//C++ headers
#include <string>

namespace protocols {
namespace jd2 {

///@details Joboutputter that only outputs score files, useful for screening protocols
class ScoreOnlyJobOutputter : public FileJobOutputter {
public:

	//constructor -- reads cmd-line to initialize evaluators
	ScoreOnlyJobOutputter();

	//////////////////////////////creating output functions/////////////////////////////////////////

	///@brief this function takes a string and writes it to disk (separately from Tracer output).
	///use some sort of extention option system - default .dat?  .data?
	//virtual
	//void file( JobCOP, std::string const & );

	///@brief this function outputs the final result of a job.
	virtual
	void final_pose( JobCOP job, core::pose::Pose const &  pose);

	///@brief this function is intended for saving mid-protocol poses; for example the final centroid structure in a combined centroid/fullatom protocol.
	virtual
	void other_pose( JobCOP job, core::pose::Pose const & pose, std::string const &  tag);

	/////////////////////////////////state of output functions/////////////////////////////////

	///@brief this function is not used for output, but it belongs here since it needs to check the same output locations as the class normally writes to.  This class checks wherever output goes to see if the job's expected output already exists (on disk or whatever).  This is the most basic form of checkpointing.
	virtual
	bool job_has_completed( JobCOP ) ;

	virtual std::string output_name( JobCOP job );

}; // NoOutputJobOutputter

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_NoOutputJobOutputter_HH
