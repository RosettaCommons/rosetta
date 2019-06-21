// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/mmtfJobOutputter.hh
/// @brief  header file for mmtfJobOutputter class, part of JD2, implemented post 2016-chemXRW
/// @author Danny Farrell danpf@uw.edu


#ifndef INCLUDED_protocols_jd2_mmtfJobOutputter_hh
#define INCLUDED_protocols_jd2_mmtfJobOutputter_hh

//unit headers
#include <protocols/jd2/mmtfJobOutputter.fwd.hh>
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

/// @details outputs mmtfs and associated files, uncompressed.
class mmtfJobOutputter : public protocols::jd2::wwPDBJobOutputter
{
public:

	typedef protocols::jd2::wwPDBJobOutputter parent;

	mmtfJobOutputter();

	~mmtfJobOutputter() override;


protected:
	//////////////////////////////////////protected mmtf output/////////////////////////////////////
	/// @brief This is the function actually different between mmtf and PDB output, and unshared by the wwPDB parent class.  Here it causes a cif file to be written.  Pure virtual in the base class.  filename is an optional label for the score data table, not an actual control.

	void dump_pose(
		JobCOP job,
		core::pose::Pose const & pose,
		utility::io::ozstream & out,
		std::string const & filename = "" ) override;

private:
	///@details this function takes "extra data" associated with the job and writes it to JOBNAME.extradata.  Filename is the location the parent mmtf got written to.  It's pretty stupid to do that via this string instead of the Job object, given that that's the Job object's purpose, but otherwise we lose the "tag" ability specified several layers up in wwPDBJO.
	void dump_extra_data_file(
		JobCOP job,
		//core::pose::Pose const & pose,
		std::string const & parent_filename );

	///@details this function takes energies associated with the pose and writes it to JOBNAME.energies.  Filename is the location the parent mmtf got written to.  It's pretty stupid to do that via this string instead of the Job object, given that that's the Job object's purpose, but otherwise we lose the "tag" ability specified several layers up in wwPDBJO.  tagging_filename is filename from dump_pose, passaged for extract_scores further along.
	void dump_energies_file(
		//JobCOP job,
		core::io::StructFileRepCOP sfr,
		std::string const & parent_filename );



	////////////////////////////////////////data////////////////////////////////////////
private:
	std::string extension_;
	//std::string path_; //exists in parent class
}; // mmtfJobOutputter

} // namespace jd2
} // namespace protocols

#endif //INCLUDED_protocols_jd2_mmtfJobOutputter_HH
