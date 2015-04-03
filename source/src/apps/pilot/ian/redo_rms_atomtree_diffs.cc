// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   demo/ian_test/redo_rms_atomtree_diffs.cc
///
/// @brief  Re-calculate ligand RMS values assuming some member of this ensemble is the correct pose.
/// @author Ian Davis (ian.w.davis@gmail.com)

// must be here to avoid VC++ ambiguous symbol w/ ObjexxFCL::byte
// for boinc builds - dek
#include <protocols/jobdist/JobDistributors.hh>


#include <numeric/conversions.hh>
#include <numeric/random/random_permutation.hh>
#include <numeric/xyzVector.io.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray1.io.hh>
#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/io/ozstream.hh>
#include <utility/pointer/owning_ptr.hh>

#include <devel/init.hh>
#include <core/types.hh>
#include <core/conformation/Residue.hh>
#include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
#include <core/scoring/rms_util.hh>
#include <basic/Tracer.hh>

#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/ligand_docking/LigandBaseProtocol.hh>


#include <ctime>
#include <fstream>
#include <set>
#include <sstream>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>


int
main( int argc, char * argv [] )
{
	try {

	using core::io::atom_tree_diffs::AtomTreeDiff;
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	using namespace protocols;
	using namespace protocols::jobdist;
	using namespace protocols::moves;
	basic::Tracer TR( "redo_rms_atomtree_diffs.main" );

	// Parses command line options and inits RNG.
	// Doesn't seem to hurt to do it again if already done once (?)
	devel::init(argc, argv);

	time_t overall_start_time = time(NULL);
	utility::vector1< BasicJobOP > input_jobs = load_s_and_l();
	// Reduce read contention between processes by randomizing the order in which structures are processed
	numeric::random::random_permutation( input_jobs, numeric::random::rg() );
	PlainPdbJobDistributor< BasicJobOP > jobdist( input_jobs );

	// A "native" pose for the diff reference point.
	// This used to be required but will be rarely used now.
	//core::pose::PoseOP native_pose;
	//if( option[ in::file::native ].user() ) {
	//	native_pose = new core::pose::Pose();
	//	core::import_pose::pose_from_pdb( *native_pose, option[ in::file::native ]().name() );
	//}

	// Tags to be compared to.
	utility::vector1< std::string > ref_tags;
	if( option[ in::file::tags ].active() ) {
		ref_tags = option[ in::file::tags ]();
	}
	// Now -native is used as an external reference point for RMS values
	core::pose::PoseOP native_pose;
	if( option[ in::file::native ].user() ) {
		native_pose = new core::pose::Pose();
		core::import_pose::pose_from_pdb( *native_pose, option[ in::file::native ]().name() );
	}
	if( ref_tags.empty() && native_pose() == NULL ) utility_exit_with_message("Must provide list of tags to compare to!");

	utility::io::ozstream out( option[ out::file::silent ]() );
	BasicJobOP curr_job;
	int curr_nstruct, num_structures_processed = 0;
	jobdist.startup();
	while( jobdist.next_job(curr_job, curr_nstruct) ) {
		TR << "Reading silent file " << curr_job->input_tag() << " ... ";
		AtomTreeDiff atdiff( curr_job->input_tag() );
		AtomTreeDiff::ScoresPairList const & scores_list = atdiff.scores();
		TR << scores_list.size() << " structures" << std::endl;

		//std::set< std::string > desired_tags;
		//protocols::ligand_docking::select_best_poses(atdiff, desired_tags);
		//TR << "Keeping " << desired_tags.size() << " best structures" << std::endl;

		// Step 1:  load structures to compare to
		utility::vector1< core::pose::PoseOP > ref_poses;
		for(core::Size ii = 1; ii <= ref_tags.size(); ++ii) {
			core::pose::PoseOP ref_pose = new core::pose::Pose();
      //if( native_pose() == NULL )
				atdiff.read_pose(ref_tags[ii], *ref_pose);
      //else atdiff.read_pose(ref_tags[ii], *ref_pose, *native_pose);
			ref_poses.push_back( ref_pose );
		}
		if( native_pose() != NULL ) {
			ref_tags.push_back("native");
			ref_poses.push_back(native_pose);
		}

		// Step 2:  actually do the rms comparison
		// By iterating in the original order and testing against desired_tags,
		// we're always seeking forward on the disk, which I think should help performance.
		for(core::Size ii = 1; ii <= scores_list.size(); ++ii) {
			std::string output_tag = scores_list[ii].first;
			// If this tag was not one of the best, skip it.
			//if( desired_tags.count(output_tag) == 0 ) continue;
			std::map< std::string, core::Real > scores = scores_list[ii].second;

			core::pose::PoseOP the_pose = new core::pose::Pose();
			//if( native_pose() == NULL )
				atdiff.read_pose(output_tag, *the_pose);
			//else atdiff.read_pose(output_tag, *the_pose, *native_pose);

			for(core::Size i = 1; i <= ref_poses.size(); ++i) {
				std::string ref_tag( ref_tags[i] );
				core::pose::PoseOP ref_pose( ref_poses[i] );

				using namespace core::scoring;
				core::Size const last_rsd = ref_pose->total_residue();
				assert( !ref_pose->residue(last_rsd).is_polymer() );
				scores[ref_tag+"_ligand_auto_rms_with_super"] = automorphic_rmsd(ref_pose->residue(last_rsd), the_pose->residue(last_rsd), true /*superimpose*/);
				scores[ref_tag+"_ligand_auto_rms_no_super"] = automorphic_rmsd(ref_pose->residue(last_rsd), the_pose->residue(last_rsd), false /*don't superimpose*/);
			}

			// Just dump SCORE lines to output, conformation didn't change.
			core::io::atom_tree_diffs::dump_score_line(out, output_tag, scores);

			num_structures_processed += 1;
			TR << "Finished " << output_tag << std::endl;
		}
	} // loop over jobs and nstructs
	jobdist.shutdown();

	time_t overall_end_time = time(NULL);
	TR << "Finished all " << num_structures_processed << " structures in " << (overall_end_time - overall_start_time) << " seconds." << std::endl;
	if ( num_structures_processed == 0 )
		basic::Warning() << "No structures processed.  Existing output files may have been skipped, did you mean to delete them?" << std::endl;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}

