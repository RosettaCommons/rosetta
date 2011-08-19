// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/jd2/EnzdesJobOutputter.hh
/// @brief header file for EnzdesJobOutputter
/// @author Florian Richter (floric@u.washington.edu), september 2010


#ifndef INCLUDED_protocols_jd2_EnzdesJobOutputter_hh
#define INCLUDED_protocols_jd2_EnzdesJobOutputter_hh

//unit headers
#include <protocols/jd2/EnzdesJobOutputter.fwd.hh>
#include <protocols/jd2/PDBJobOutputter.hh>

//project headers
#include <core/io/silent/SilentFileData.fwd.hh>
#include <protocols/enzdes/EnzFilters.fwd.hh>
#include <protocols/jd2/SilentFileJobOutputter.fwd.hh>

namespace protocols {
namespace jd2 {

/// @brief for now this class only writes a different scorefile
/// than the default one written by the FileJobOutputter. the structure
/// output format is pdb
class EnzdesJobOutputter : public protocols::jd2::PDBJobOutputter
{
public: //constructor / destructor

  typedef protocols::jd2::PDBJobOutputter parent;

  EnzdesJobOutputter();

  ~EnzdesJobOutputter();

	void final_pose( JobCOP job, core::pose::Pose const & pose );

	bool job_has_completed( JobCOP job );

protected: //Job Outputter interface

  void scorefile(
    JobCOP job,
    core::pose::Pose const & pose,
    std::string tag = "",
    std::string scorefile = ""
 );

private: //data

	core::io::silent::SilentFileDataOP scorefile_writer_;
	protocols::enzdes::EnzdesScorefileFilterOP enz_scofile_;
	bool silent_output_;
	SilentFileJobOutputterOP silent_job_outputter_;
};

}//jd2
}//protocols



#endif //INCLUDED_protocols_jd2_EnzdesJobOutputter_FWD_HH
