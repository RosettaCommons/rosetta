// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/mmCIFJobOutputter.hh
/// @brief  header file for mmCIFJobOutputter class, part of JD2, implemented post 2016-chemXRW
/// @author Steven Lewis smlewi@gmail.com


#ifndef INCLUDED_protocols_jd2_mmCIFJobOutputter_hh
#define INCLUDED_protocols_jd2_mmCIFJobOutputter_hh

//unit headers
#include <protocols/jd2/mmCIFJobOutputter.fwd.hh>
#include <protocols/jd2/wwPDBJobOutputter.hh>
#include <protocols/jd2/Job.fwd.hh>

//project headers
#include <core/pose/Pose.fwd.hh>
#include <core/io/StructFileRep.fwd.hh>

//utility headers

//C++ headers
#include <string>

#include <utility/vector1.hh>
#include <utility/io/ozstream.fwd.hh>
#include <iostream>


namespace protocols {
namespace jd2 {

/// @details outputs mmCIFs and associated files, uncompressed.
class mmCIFJobOutputter : public protocols::jd2::wwPDBJobOutputter
{
public:

	typedef protocols::jd2::wwPDBJobOutputter parent;

	mmCIFJobOutputter();

	virtual ~mmCIFJobOutputter();

	//////////////////////////////creating output functions/////////////////////////////////////////

	/// @brief this function outputs the final result of a job.  This implementation will write a mmCIF file (plus scores).
	//Exists in parent wwPDBJO
	//virtual
	//void final_pose( JobOP job, core::pose::Pose const & pose, std::string const & tag );

	/// @brief this function is intended for saving mid-protocol poses; for example the final centroid structure in a combined centroid/fullatom protocol.  This implementation will write a mmCIF file (plus scores).
	//Exists in parent wwPDBJO
	//virtual
	//void other_pose( JobOP job, core::pose::Pose const & pose, std::string const & tag, int copy_count = -1, bool score_only = false );

	/////////////////////////////////state of output functions/////////////////////////////////

	/// @brief this function is not used for output, but it belongs here since it needs to check the same output locations as the class normally writes to.  This class checks wherever output goes to see if the job's expected output already exists (on disk or whatever).  This is the most basic form of checkpointing.  The base implementation looks for an mmCIF with the job's name already in existence.
	//Exists in parent wwPDBJO
	//virtual
	//bool job_has_completed( JobCOP job );

	/// @brief this is the master function for determining the unique output identifier for a job
	//Exists in parent wwPDBJO
	//virtual
	//std::string output_name( JobCOP job );

protected:
	/// @brief this private function provides the extended name, not just the output name.  e.g output_name returns 1UBQ_0001, this returns 1UBQ_0001.cif.  In this case the extension is .cif
	//Exists in parent wwPDBJO
	//virtual
	//std::string extended_name( JobCOP job, std::string const suffix = "" );

	////////////////////////////////////////score-related functions///////////////////////////////////

	/// @details Refactored (in PDBJobOutputter) in the 2016 Chemical XRW (eXtreme Rosetta Workshop) by Vikram K. Mulligan (vmullig@uw.edu).
	/// @param[in] job Const-access owning pointer to the job from which the data will be extracted.
	/// @param[out] data_out A string in which the data will be stored, that can later be passed to whatever container wants it.
	//Exists in parent wwPDBJO
	//virtual
	//std::string extract_data_from_Job( JobCOP job );

protected:
	//////////////////////////////////////protected mmCIF output/////////////////////////////////////
	/// @brief This is the function actually different between mmCIF and PDB output, and unshared by the wwPDB parent class.  Here it causes a cif file to be written.  Pure virtual in the base class.  filename is an optional label for the score data table, not an actual control.
	virtual
	void dump_pose(
		JobCOP job,
		core::pose::Pose const & pose,
		utility::io::ozstream & out,
		std::string const & filename = "" );

private:
	///@details this function takes "extra data" associated with the job and writes it to JOBNAME.extradata.  Filename is the location the parent mmCIF got written to.  It's pretty stupid to do that via this string instead of the Job object, given that that's the Job object's purpose, but otherwise we lose the "tag" ability specified several layers up in wwPDBJO.
	void dump_extra_data_file(
		JobCOP job,
		//core::pose::Pose const & pose,
		std::string const & parent_filename );

	///@details this function takes energies associated with the pose and writes it to JOBNAME.energies.  Filename is the location the parent mmCIF got written to.  It's pretty stupid to do that via this string instead of the Job object, given that that's the Job object's purpose, but otherwise we lose the "tag" ability specified several layers up in wwPDBJO.  tagging_filename is filename from dump_pose, passaged for extract_scores further along.
	void dump_energies_file(
		//JobCOP job,
		core::io::StructFileRepCOP sfr,
		std::string const & parent_filename );



	////////////////////////////////////////data////////////////////////////////////////
private:
	std::string extension_;
	//std::string path_; //exists in parent class
}; // mmCIFJobOutputter

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_mmCIFJobOutputter_HH
