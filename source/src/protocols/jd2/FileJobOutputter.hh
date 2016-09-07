// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/jd2/FileJobOutputter.hh
/// @brief  header file for FileJobOutputter class, part of August 2008 job distributor as planned at RosettaCon08
/// @author Steven Lewis smlewi@gmail.com


#ifndef INCLUDED_protocols_jd2_FileJobOutputter_hh
#define INCLUDED_protocols_jd2_FileJobOutputter_hh

//unit headers
#include <protocols/jd2/FileJobOutputter.fwd.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

//utility headers
#include <utility/file/FileName.hh>

//C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {

/// @details this is a middle-layer implementation of JobOutputter for file-based output.  It handles scorefile output, as scorefiles are common to file-outputting methods.
class FileJobOutputter : public protocols::jd2::JobOutputter
{
public:

	typedef protocols::jd2::JobOutputter parent;

	FileJobOutputter();

	~FileJobOutputter() override;

	void set_defaults();
	//////////////////////////////creating output functions/////////////////////////////////////////

	/// @brief this function takes a string and writes it to disk (separately from Tracer output).

	void file( JobCOP job, std::string const & data ) override;

	/// @brief this function outputs the final result of a job.

	void final_pose( JobOP job, core::pose::Pose const & pose, std::string const & tag ) override = 0;

	/// @brief this function is intended for saving mid-protocol poses; for example the final centroid structure in a combined centroid/fullatom protocol.

	void other_pose( JobOP job, core::pose::Pose const & pose, std::string const & tag, int copy_count = -1, bool score_only = false ) override = 0;

	/////////////////////////////////state of output functions/////////////////////////////////

	/// @brief this function is not used for output, but it belongs here since it needs to check the same output locations as the class normally writes to.  This class checks wherever output goes to see if the job's expected output already exists (on disk or whatever).  This is the most basic form of checkpointing.  The base implementation looks for a pdb with the job's name already in existence.

	bool job_has_completed( JobCOP job ) override = 0;

	/// @brief this is the master function for determining the unique output identifier for a job

	std::string output_name( JobCOP job ) override = 0;

	utility::file::FileName const & scorefile_name() {
		return scorefile_name_;
	}

	bool
	write_scorefile() const {
		return write_scorefile_; }

	//////////////////////////////////////scorefile functions/////////////////////////////////////
protected:
	/// @brief this function will handle the scorefile.  If you need to make it virtual do so.  Latter two arguments are for redirecting the output to a different scorefile for "other_pose"s.  Also adds StringReal job info to the score file.
	virtual
	void scorefile( JobCOP job, core::pose::Pose const & pose, std::string prefix_tag = "", std::string suffix_tag = "", std::string scorefile = "" );

	/// @brief this function will handle the scorefile for arbitrary poses.
	//void other_scorefile( std::string const & tag, core::pose::Pose const & pose, std::string const & o_scorefile );

private:
	bool write_scorefile_;
	utility::file::FileName scorefile_name_;

}; // FileJobOutputter

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_FileJobOutputter_HH
