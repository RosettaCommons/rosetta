// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

#ifndef INCLUDED_protocols_ddg_ddGData_hh
#define INCLUDED_protocols_ddg_ddGData_hh

// Unit header
#include <protocols/ddg/ddGData.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <fstream>

#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/vector1.hh>
#include <map>


namespace protocols {
namespace ddG {

class ddGData {
public:
	//constructors
	ddGData();

	ddGData(std::string filename);

	//virtual ~ddGData();

	//iterator functions
	bool end();
	void get_next_filenames();
	utility::vector1< core::pose::Pose > read_mut_data();
	utility::vector1< core::pose::Pose > read_wt_data();
	core::Real read_exp_data();

private:

	std::string file_to_read_;
	std::string curr_mut_filename;
	std::string curr_wt_filename;
	core::Real experimental;
	bool curr_mut_read;
	bool curr_wt_read;
	std::ifstream inputstream;
};//class ddGData


}//namespace ddG
} //namespace protocols

#endif
