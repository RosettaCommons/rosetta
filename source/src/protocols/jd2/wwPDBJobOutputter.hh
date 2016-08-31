// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/wwPDBJobOutputter.hh
/// @brief  header file for wwPDBJobOutputter class, interstitial between FileJobOutputter and PDBJobOutputter or mmCIFJobOutputter to share functionality.
/// @author Steven Lewis smlewi@gmail.com


#ifndef INCLUDED_protocols_jd2_wwPDBJobOutputter_hh
#define INCLUDED_protocols_jd2_wwPDBJobOutputter_hh

//unit headers
#include <protocols/jd2/wwPDBJobOutputter.fwd.hh>
#include <protocols/jd2/FileJobOutputter.hh>
#include <protocols/jd2/Job.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>

//utility headers

//C++ headers
#include <string>

#include <utility/vector1.hh>
#include <utility/io/ozstream.fwd.hh>
#include <iostream>


namespace protocols {
namespace jd2 {

/// @details this simplest implementation of JobOutputter outputs raw wwPDBs and associated files, uncompressed.
class wwPDBJobOutputter : public protocols::jd2::FileJobOutputter
{
public:

	typedef protocols::jd2::FileJobOutputter parent;

	wwPDBJobOutputter();

	virtual ~wwPDBJobOutputter();

	//////////////////////////////creating output functions/////////////////////////////////////////

	/// @brief this function outputs the final result of a job.  This implementation will write a wwPDB-format file (plus scores).  It calls a pure virtual so that child classes can write PDB or mmCIF format.
	virtual
	void final_pose( JobOP job, core::pose::Pose const & pose, std::string const & tag );

	/// @brief this function is intended for saving mid-protocol poses; for example the final centroid structure in a combined centroid/fullatom protocol.  This implementation will write a wwPDB-format file (plus scores).  It calls a pure virtual so that child classes can write PDB or mmCIF format.
	virtual
	void other_pose( JobOP job, core::pose::Pose const & pose, std::string const & tag, int copy_count = -1, bool score_only = false );

	/////////////////////////////////state of output functions/////////////////////////////////

	/// @brief this function is not used for output, but it belongs here since it needs to check the same output locations as the class normally writes to.  This class checks wherever output goes to see if the job's expected output already exists (on disk or whatever).  This is the most basic form of checkpointing.  The base implementation looks for a pdb/cif with the job's name already in existence.
	virtual
	bool job_has_completed( JobCOP job );

	/// @brief this is the master function for determining the unique output identifier for a job
	virtual
	std::string output_name( JobCOP job );

protected:
	/// @brief this private function provides the extended name, not just the output name.  e.g output_name returns 1UBQ_0001, this returns 1UBQ_0001.pdb or 1UBQ_0001.cif.
	virtual
	std::string extended_name( JobCOP job, std::string const & suffix = "" );

	////////////////////////////////////////score-related functions///////////////////////////////////

	/// @brief this function extracts the pose's scores and outputs them as a string to be packaged in an output structure.
	/// @details Refactored in the 2016 Chemical XRW (eXtreme Rosetta Workshop) by Vikram K. Mulligan (vmullig@uw.edu).
	/// @param[in] job Const-access owning pointer to the job from which the data will be extracted.
	/// @param[out] data_out A string in which the data will be stored, that can later be passed to whatever container wants it.
	virtual
	std::string extract_data_from_Job( JobCOP job );

protected:
	/// @brief This pure virtual is implemented by child functions to actually write a file - this is where mmCIF and PDB output differ.
	virtual
	void dump_pose(
		JobCOP job,
		core::pose::Pose const & pose,
		utility::io::ozstream & out,
		std::string const & filename) = 0;

	///////////////////protected: child classes PDBJO and mmCIFJO need this datum.////////////////////////////////
protected:

	/// @brief setter for output file paths, in case child class needs to override -out:path:all with -out:path:[PDB/mmCIF]
	void set_path(std::string const & path);

	/// @brief getter for output file path
	std::string const & get_path();

	/// @brief setter for output file extensions, child class must set
	void set_extension(std::string const & extension);

	/// @brief getter for output file extension
	std::string const & get_extension();

	////////////////////////////////////////data////////////////////////////////////////
private:
	std::string path_;
	std::string extension_;

}; // wwPDBJobOutputter

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_wwPDBJobOutputter_HH
