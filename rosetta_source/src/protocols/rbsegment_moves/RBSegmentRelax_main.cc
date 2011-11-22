// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Srivatsan Raman
/// @author Frank DiMaio

#include <protocols/jobdist/JobDistributors.hh>
#include <protocols/jobdist/Jobs.hh>
// AUTO-REMOVED #include <protocols/jobdist/standard_mains.hh>
#include <core/types.hh>

// AUTO-REMOVED #include <core/init.hh>

#include <core/kinematics/Jump.hh>
#include <core/fragment/FragSet.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/import_pose/import_pose.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
// AUTO-REMOVED #include <basic/options/after_opts.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>

#include <protocols/loops/loops_main.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMover_QuickCCD.hh>

// AUTO-REMOVED #include <numeric/random/random.hh>
// AUTO-REMOVED #include <numeric/random/random_permutation.hh>

#include <core/kinematics/FoldTree.hh>

// AUTO-REMOVED #include <protocols/rbsegment_Moves/RBSegmentMover.hh>
#include <protocols/rbsegment_moves/RBSegmentRelax.hh>
// AUTO-REMOVED #include <protocols/rbsegment_Moves/FragInsertAndAlignMover.hh>
// AUTO-REMOVED #include <protocols/loops/LoopBuild.hh>
// AUTO-REMOVED #include <protocols/viewer/viewers.hh>
// AUTO-REMOVED #include <protocols/relax_protocols.hh>
// AUTO-REMOVED #include <protocols/frags/TorsionFragment.hh>
// AUTO-REMOVED #include <protocols/evaluation/RmsdEvaluator.hh>
// AUTO-REMOVED #include <protocols/moves/PackRotamersMover.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/task/TaskFactory.hh>

#include <core/io/silent/SilentStructFactory.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
// AUTO-REMOVED #include <ctime>

//silly using/typedef

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/loops.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/RBSegmentRelax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

//Auto Headers
#include <utility/vector1.hh>
//Auto Headers
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/rbsegment_moves/RBSegment.hh>
#include <basic/options/option.hh>




using basic::T;
using basic::Error;
using basic::Warning;

basic::Tracer TRb("rbsegmove_main");

namespace protocols {
// namespace rbsegment_moves {

////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
int
RBSegmentRelax_main( bool boinc_mode ) {
	using namespace rbsegment_moves;
	using namespace jobdist;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::id;

	//////////////////
	// rbmover scorefn
	core::scoring::ScoreFunctionOP scorefxn_rb
	         = core::scoring::ScoreFunctionFactory::create_score_function( option[ OptionKeys::RBSegmentRelax::rb_scorefxn ]() );

	core::pose::Pose start_pose, pose, native_pose;
	core::chemical::ResidueTypeSetCAP rsd_set;

	std::string pdbfilename;
	if ( option[ OptionKeys::RBSegmentRelax::input_pdb ].user() )
		pdbfilename = option[ OptionKeys::RBSegmentRelax::input_pdb ]().name();
	else
		pdbfilename = option[ OptionKeys::in::file::s ]()[1];

	// if full-atom load starting structure as full-atom to recover sidechains later
	if ( option[ in::file::fullatom ]() ) {
		core::import_pose::pose_from_pdb( start_pose, pdbfilename );
	} else {
		core::import_pose::centroid_pose_from_pdb( start_pose, pdbfilename );
	}

	// roughly guess at secondary structure
	core::pose::set_ss_from_phipsi( start_pose );

	// native structure
	if ( option[ in::file::native ].user() ) {
		core::import_pose::pose_from_pdb( native_pose, option[ OptionKeys::in::file::native ]() );
    core::pose::set_ss_from_phipsi( native_pose );

#ifdef BOINC_GRAPHICS
    // set native for graphics
    boinc::Boinc::set_graphics_native_pose( native_pose );
#endif

		core::util::switch_to_residue_type_set( native_pose, core::chemical::CENTROID );
	}

	// Read RB segs, auto generate loops
	utility::vector1< protocols::rbsegment_moves::RBSegment > rbsegs;
	utility::vector1< int > cutpts;
	protocols::loops::Loops loops;
	std::string rbfilename( option[ OptionKeys::RBSegmentRelax::rb_file ]().name() );

	for (int i=1; i<=start_pose.fold_tree().num_cutpoint() ; ++i)
		cutpts.push_back( start_pose.fold_tree().cutpoint(i) );
	int last_peptide_res = start_pose.total_residue();
	while ( !start_pose.residue( last_peptide_res ).is_protein() )
		last_peptide_res--;
	protocols::rbsegment_moves::read_RBSegment_file( rbsegs, loops, rbfilename, true, last_peptide_res , cutpts  );

	////////////////////////
	////////////////////////
	// job distributor initialization
	utility::vector1< protocols::jobdist::BasicJobOP > input_jobs;
	int const nstruct_flag = option[ out::nstruct ];
	int const nstruct = std::max( 1, nstruct_flag );
	protocols::jobdist::BasicJobOP job = new protocols::jobdist::BasicJob("S", "rbseg", nstruct);
	input_jobs.push_back( job );
	protocols::jobdist::BaseJobDistributorOP jobdist;

	// output nonidealized silent file or PDBs?
	bool silent_output;
	if ( boinc_mode || option[ OptionKeys::out::file::silent ].user() ) {
		TRb.Debug << "Outputting silent file\n";
		jobdist = new protocols::jobdist::PlainSilentFileJobDistributor( input_jobs );
		silent_output = true;
	} else {
		TRb.Debug << "Outputting PDBs\n";
		jobdist = new protocols::jobdist::PlainPdbJobDistributor( input_jobs );
		silent_output = false;
	}

	protocols::jobdist::BasicJobOP prev_job, curr_job;
	int curr_nstruct;
	jobdist->startup();

	// read fragments
	utility::vector1< core::fragment::FragSetOP > frag_libs;
	if ( option[ OptionKeys::loops::frag_files ].user() )
		protocols::loops::read_loop_fragments( frag_libs );

	/////
	/////
	while ( jobdist->next_job(curr_job, curr_nstruct) ) { // loop over jobs
		std::string curr_job_tag = curr_job->output_tag( curr_nstruct );

		pose = start_pose;
//		if ( option[ in::file::fullatom ]() )
//			core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );

#ifdef BOINC_GRAPHICS
	// attach boinc graphics pose observer
		protocols::boinc::Boinc::attach_graphics_current_pose_observer( pose );
#endif

		// the rigid body movement mover
		RBSegmentRelax shaker( scorefxn_rb, rbsegs, loops );
		shaker.initialize( frag_libs );
		shaker.set_randomize( 2 );  //???
		shaker.apply( pose );

		////
		////  output
		if ( silent_output ) {
			PlainSilentFileJobDistributor *jd =
					 dynamic_cast< PlainSilentFileJobDistributor * > (jobdist());

			std::string silent_struct_type( "binary" );  // default to binary
			if ( option[ out::file::silent_struct_type ].user() ) {
				silent_struct_type = option[ OptionKeys::out::file::silent_struct_type ];
			}

			core::io::silent::SilentStructOP ss
				= core::io::silent::SilentStructFactory::get_instance()->get_silent_struct( silent_struct_type );

			ss->fill_struct( pose, curr_job_tag );

			jd->dump_silent( curr_nstruct, *ss );
		} else {
			jobdist->dump_pose_and_map( curr_job_tag, pose );    // output PDB
		}
	}
	jobdist->shutdown();
	return 0;
}

//}
} // namespace protocols

////////////////////////////////////////////////////////

