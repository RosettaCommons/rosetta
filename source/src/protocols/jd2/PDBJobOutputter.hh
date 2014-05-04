// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/PDBJobOutputter.hh
/// @brief  header file for PDBJobOutputter class, part of August 2008 job distributor as planned at RosettaCon08
/// @author Steven Lewis smlewi@gmail.com


#ifndef INCLUDED_protocols_jd2_PDBJobOutputter_hh
#define INCLUDED_protocols_jd2_PDBJobOutputter_hh

//unit headers
#include <protocols/jd2/PDBJobOutputter.fwd.hh>
#include <protocols/jd2/FileJobOutputter.hh>
#include <protocols/jd2/Job.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

//utility headers
// AUTO-REMOVED #include <utility/io/ozstream.hh>

//C++ headers
#include <string>

#include <utility/vector1.hh>
#include <utility/io/ozstream.fwd.hh>
#include <iostream>


namespace protocols {
namespace jd2 {

///@details this simplest implementation of JobOutputter outputs raw PDBs and associated files, uncompressed.
class PDBJobOutputter : public protocols::jd2::FileJobOutputter
{
public:

	typedef protocols::jd2::FileJobOutputter parent;

	PDBJobOutputter();

	virtual ~PDBJobOutputter();

	//////////////////////////////creating output functions/////////////////////////////////////////

	///@brief this function takes a string and writes it to disk (separately from Tracer output).  This implementation writes a single file whose filename is based on the job and a user-specified extension (default .data)
	//	virtual --> moved to FileJobOutputter
	//	void file( JobCOP job, std::string const & data );

	///@brief this function outputs the final result of a job.  This implementation will write a PDB file (plus scores).
	virtual
	void final_pose( JobOP job, core::pose::Pose const & pose, std::string const & tag );

	///@brief this function is intended for saving mid-protocol poses; for example the final centroid structure in a combined centroid/fullatom protocol.  This implementation will write a PDB file (plus scores).
	virtual
	void other_pose( JobOP job, core::pose::Pose const & pose, std::string const & tag, int copy_count = -1, bool score_only = false );

	/////////////////////////////////state of output functions/////////////////////////////////

	///@brief this function is not used for output, but it belongs here since it needs to check the same output locations as the class normally writes to.  This class checks wherever output goes to see if the job's expected output already exists (on disk or whatever).  This is the most basic form of checkpointing.  The base implementation looks for a pdb with the job's name already in existence.
	virtual
	bool job_has_completed( JobCOP job );

	///@brief this is the master function for determining the unique output identifier for a job
	virtual
	std::string output_name( JobCOP job );

protected:
	///@brief this private function provides the extended name, not just the output name.  e.g output_name returns 1UBQ_0001, this returns 1UBQ_0001.pdb.  In this case the extension is .pdb
	virtual
	std::string extended_name( JobCOP job, std::string const suffix = "" );

	////////////////////////////////////////score-related functions///////////////////////////////////

	///@brief this function extracts the pose's scores for printing
	virtual
	void extract_scores( core::pose::Pose const & pose, utility::io::ozstream & out );

	//THIS FUNCTION WILL MOVE HIGHER IN THE HIERARCHY AT SOME POINT
	///@brief this function extracts the pose's scores for printing
	virtual
	void extract_data_from_Job( JobCOP job, utility::io::ozstream & out );

	//This function is deprecated for now - might return in the future
// 	///@brief this function extracts the pose's extra data/scores for printing
// 	virtual
// 	void extract_extra_scores( core::pose::Pose const & pose, utility::io::ozstream & out );

	//////////////////////////////////////protected PDB output/////////////////////////////////////
	///@brief handles ozstream output; shared by both pdb output functions
	virtual
	void dump_pose( JobCOP job, core::pose::Pose const & pose, utility::io::ozstream & out );

	////////////////////////////////////////data////////////////////////////////////////
private:
	std::string extension_;
    std::string path_;
}; // PDBJobOutputter

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_PDBJobOutputter_HH
