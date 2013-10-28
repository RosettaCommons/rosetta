// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/AtomTreeDiffJobOutputter.hh
/// @brief  header file for AtomTreeDiffJobOutputter class
/// @author Gordon Lemmon (gordon.h.lemmon@vanderbilt.edu); Rocco Moretti (rmoretti@u.washington.edu)


#ifndef INCLUDED_protocols_jd2_AtomTreeDiffJobOutputter_hh
#define INCLUDED_protocols_jd2_AtomTreeDiffJobOutputter_hh

//unit headers
#include <protocols/jd2/AtomTreeDiffJobOutputter.fwd.hh>
#include <protocols/jd2/FileJobOutputter.hh>
#include <protocols/jd2/JobOutputter.hh>
#include <protocols/jd2/Job.fwd.hh>

//project headers
#include <core/pose/Pose.hh>
#include <core/import_pose/atom_tree_diffs/atom_tree_diff.hh>

//utility headers
#include <utility/file/FileName.hh>
#include <utility/io/ozstream.hh>

//C++ headers
#include <string>
#include <set>

#include <utility/vector1.hh>


namespace protocols {
namespace jd2 {

///@details this is a middle-layer implementation of JobOutputter for file-based output.  It handles scorefile output, as scorefiles are common to file-outputting methods.
class AtomTreeDiffJobOutputter : public protocols::jd2::FileJobOutputter
{
public:

	AtomTreeDiffJobOutputter();

	~AtomTreeDiffJobOutputter();

	//////////////////////////////creating output functions/////////////////////////////////////////

	///@brief this function outputs the final result of a job.
	void final_pose( JobOP job, core::pose::Pose const & pose );

	///@brief this function is intended for saving mid-protocol poses; for example the final centroid structure in a combined centroid/fullatom protocol.
	void other_pose( JobOP job, core::pose::Pose const & pose, std::string const & tag, int copy_count = -1, bool score_only = false );

	/////////////////////////////////state of output functions/////////////////////////////////

	///@brief this function is not used for output, but it belongs here since it needs to check the same output locations as the class normally writes to.  This class checks wherever output goes to see if the job's expected output already exists (on disk or whatever).  This is the most basic form of checkpointing.  The base implementation looks for a pdb with the job's name already in existence.
	bool job_has_completed( JobCOP job );

	///@brief this is the master function for determining the unique output identifier for a job
	std::string output_name( JobCOP job );

	///@brief what precision should the atom_tree_diff be ouput at?
	void set_precision(int bb_precision, int sc_precision, int bondlen_precision);

	///@brief use input as reference pose?
	void use_input_for_ref(bool use_input=true);

private:

	///@brief Appends pose to the silent file
	void
	dump_pose(
		std::string const & tag,
		core::pose::Pose const & pose,
		std::map< std::string, core::Real > scores,
		JobCOP job
	);

	utility::io::ozstream out_;
	std::string outfile_name_;
	std::set< std::string > used_tags_;
	std::string last_ref_tag_;
	core::pose::Pose last_ref_pose_;
	core::import_pose::atom_tree_diffs::AtomTreeDiff atom_tree_diff_;

	///@brief precision to output atom tree diff at.
	int bb_precision_, sc_precision_, bondlen_precision_;

	///@brief use input as reference pose? (default false)
	bool use_input_;

}; // AtomTreeDiffJobOutputter

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_AtomTreeDiffJobOutputter_HH
