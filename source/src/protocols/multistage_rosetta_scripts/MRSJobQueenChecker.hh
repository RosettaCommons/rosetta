// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/multistage_rosetta_scripts/MRSJobQueenChecker.hh
/// @brief This class is derived from the MRSJobQueen and is intended to be a "dry run" equivalent - used to check conditions of the run without doing heavy computation
/// @author Jack Maguire, jack@med.unc.edu


#ifndef INCLUDED_protocols_multistage_rosetta_scripts_MRSJobQueenChecker_HH
#define INCLUDED_protocols_multistage_rosetta_scripts_MRSJobQueenChecker_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/multistage_rosetta_scripts/MRSJobQueenChecker.fwd.hh>
#include <protocols/multistage_rosetta_scripts/MRSJobQueen.hh>

//JD3
#include <protocols/jd3/JobDigraph.fwd.hh>
#include <protocols/jd3/LarvalJob.fwd.hh>

namespace protocols {
namespace multistage_rosetta_scripts {

class MRSJobQueenChecker : public MRSJobQueen {

public:

	//constructor
	MRSJobQueenChecker();

	//destructor
	~MRSJobQueenChecker();

public:
	///@brief create serialized poses for every input job and use these to predict the amount of memory required for archiving results
	core::Size estimate_number_of_bytes_needed_for_archiving();

	static std::pair< core::Size, core::Size >
	fa_and_cen_sizes_for_archives( core::pose::PoseOP pose );

public://overrides

	jd3::JobDigraphOP
	initial_job_dag()
	override;

	std::list< jd3::LarvalJobOP > determine_job_list( core::Size job_dag_node_index, core::Size max_njobs ) override;

};

} //multistage_rosetta_scripts
} //protocols

#endif
