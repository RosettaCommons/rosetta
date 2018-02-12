// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/multistage_rosetta_scripts/MRSJobQueenChecker.cc
/// @brief
/// @detailed
/// @author Jack Maguire, jackmaguire1444@gmail.com


#include <protocols/multistage_rosetta_scripts/MRSJobQueenChecker.hh>

#include <basic/Tracer.hh>
#include <core/pose/Pose.hh>

#include <protocols/jd3/LarvalJob.hh>
#include <protocols/jd3/JobDigraph.hh>
#include <protocols/simple_moves/SwitchResidueTypeSetMover.hh>

#ifdef SERIALIZATION
#include <cereal/archives/binary.hpp>
#include <cereal/types/polymorphic.hpp>
#include <cereal/types/list.hpp>
#include <cereal/types/utility.hpp>
#include <cereal/details/helpers.hpp>
#endif

static basic::Tracer TR( "protocols.multistage_rosetta_scripts.MRSJobQueenChecker" );

namespace protocols {
namespace multistage_rosetta_scripts {

//Constructor
MRSJobQueenChecker::MRSJobQueenChecker() : MRSJobQueen()
{}

//Destructor
MRSJobQueenChecker::~MRSJobQueenChecker()
{}

core::Size
MRSJobQueenChecker::estimate_number_of_bytes_needed_for_archiving(){

#if defined(SERIALIZATION)
	TR << "Estimating total memory required for archiving" << std::endl;
#else
	utility_exit_with_message( "The memory estimation tool is only currently available with extras=serialization" );
#endif

	core::Size max_num_bytes = 0;

	core::Size max_fa_pose_size = 0;

	utility::vector1< jd3::standard::PreliminaryLarvalJob > const & input_jobs = preliminary_larval_jobs();

	for ( core::Size ii = 1; ii <= num_input_structs(); ++ii ) {
		jd3::standard::PreliminaryLarvalJob ljob_ii = input_jobs[ ii ];
		core::pose::PoseOP pose_ii = pose_for_inner_job( ljob_ii.inner_job, * options_for_job( * ljob_ii.inner_job ) );
		std::pair< core::Size, core::Size > fa_and_cen_sizes = fa_and_cen_sizes_for_archives( pose_ii );

		TR << "job #" << ii << "'s pose is " << fa_and_cen_sizes.first << " bytes in full-atom mode and " << fa_and_cen_sizes.second << " bytes in centroid mode" << std::endl;

		if ( max_fa_pose_size  < fa_and_cen_sizes.first  ) max_fa_pose_size  = fa_and_cen_sizes.first;
		//if( max_cen_pose_size < fa_and_cen_sizes.second ) max_cen_pose_size = fa_and_cen_sizes.second;
	}

	short unsigned int max_stage = 1;
	for ( short unsigned int i=1; i < num_stages(); ++i ) {
		core::Size const num_possible_bytes = max_fa_pose_size * i * num_results_to_keep_for_stage( i );
		//core::Size const num_possible_bytes = max_fa_pose_size * 2 * num_results_to_keep_for_stage( i );//TODO add this as an option
		if ( num_possible_bytes > max_num_bytes ) {
			max_stage = i;
			max_num_bytes = num_possible_bytes;
		}
	}

	TR << "Assuming all poses are in full-atom mode, the maximum amount of memory needed at once is after node " << max_stage << " : " << max_num_bytes << " bytes." << std::endl;

	return max_num_bytes;
}

inline core::Size
archive_size_for_pose( core::pose::PoseOP pose ){
#if /*defined(USEMPI) &&*/ defined(SERIALIZATION)
	std::ostringstream oss;
	{
		cereal::BinaryOutputArchive arc( oss );
		arc( pose );
	}
	return oss.str().size();
#else
	return pose->size();
#endif
}

std::pair< core::Size, core::Size >
MRSJobQueenChecker::fa_and_cen_sizes_for_archives( core::pose::PoseOP pose ){
	//TODO do some kind of scoring here?

	protocols::simple_moves::SwitchResidueTypeSetMover to_fa( "fa_standard" );
	to_fa.apply( *pose );
	core::Size const fa = archive_size_for_pose( pose );

	protocols::simple_moves::SwitchResidueTypeSetMover to_cen( "centroid" );
	to_cen.apply( *pose );
	core::Size const cen = archive_size_for_pose( pose );

	return std::make_pair( fa, cen );
}

jd3::JobDigraphOP
MRSJobQueenChecker::initial_job_dag() {
	MRSJobQueen::initial_job_dag();
	determine_validity_of_stage_tags();
	jd3::JobDigraphOP const dag ( new jd3::JobDigraph( 0 ) );
	return dag;
}

std::list< jd3::LarvalJobOP >
MRSJobQueenChecker::determine_job_list( core::Size, core::Size )  {
	std::list< jd3::LarvalJobOP > const jobs;
	return jobs;
}


} //multistage_rosetta_scripts
} //protocols
