// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/NoOutputJobOutputter.hh
/// @brief  header file for ScoreOnlyJobOutputter class, allows the user to output only score files
/// @author Steven Lewis smlewi@gmail.com


#ifndef INCLUDED_protocols_jd2_ScoreOnlyJobOutputter_hh
#define INCLUDED_protocols_jd2_ScoreOnlyJobOutputter_hh

//unit headers
#include <protocols/jd2/FileJobOutputter.hh>
#include <protocols/jd2/ScoreOnlyJobOutputter.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

//utility headers
//C++ headers
#include <string>

#include <core/types.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {

/// @details Joboutputter that only outputs score files, useful for screening protocols
class ScoreOnlyJobOutputter : public FileJobOutputter {
public:

	//constructor -- reads cmd-line to initialize evaluators
	ScoreOnlyJobOutputter();

	//////////////////////////////creating output functions/////////////////////////////////////////

	/// @brief this function takes a string and writes it to disk (separately from Tracer output).
	///use some sort of extention option system - default .dat?  .data?
	//virtual
	//void file( JobCOP, std::string const & );

	/// @brief this function outputs the final result of a job.

	void final_pose( JobOP job, core::pose::Pose const & pose, std::string const & tag ) override;

	/// @brief this function is intended for saving mid-protocol poses; for example the final centroid structure in a combined centroid/fullatom protocol.

	void other_pose( JobOP job, core::pose::Pose const & pose, std::string const &  tag, int copy_count = -1, bool score_only = false) override;

	/////////////////////////////////state of output functions/////////////////////////////////

	/// @brief this function is not used for output, but it belongs here since it needs to check the same output locations as the class normally writes to.  This class checks wherever output goes to see if the job's expected output already exists (on disk or whatever).  This is the most basic form of checkpointing.

	bool job_has_completed( JobCOP job ) override ;

	std::string output_name( JobCOP job ) override;

private:

	void read_done_jobs();

	// list of tags already written
	utility::vector1< std::string > score_file_tags_;

}; // ScoreOnlyJobOutputter

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_NoOutputJobOutputter_HH
