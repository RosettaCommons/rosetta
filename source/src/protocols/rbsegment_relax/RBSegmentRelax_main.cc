// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Srivatsan Raman
/// @author Frank DiMaio

//#include <protocols/jobdist/JobDistributors.hh>
//#include <protocols/jobdist/Jobs.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <core/types.hh>

#include <core/kinematics/Jump.hh>
#include <core/fragment/FragSet.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>

#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <protocols/loops/loops_main.hh>
#include <core/kinematics/FoldTree.hh>

#include <protocols/rbsegment_relax/RBSegmentRelax.hh>
#include <protocols/rbsegment_relax/RBSegmentRelax_main.hh>

#include <core/io/silent/SilentStructFactory.hh>

#include <protocols/boinc/boinc.hh>  // required for boinc build

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

// option key includes
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/RBSegmentRelax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

//Auto Headers
#include <utility/vector1.hh>
#include <core/chemical/ChemicalManager.fwd.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/util/SwitchResidueTypeSet.hh>
#include <protocols/rbsegment_relax/RBSegment.hh>
#include <basic/options/option.hh>

using basic::T;
using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TRb( "rbsegmove_main" );

namespace protocols {
// namespace rbsegment_relax {

////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

RBSegmentRelaxImpl::RBSegmentRelaxImpl(){
	//////////////////
	// rbmover scorefn
	scorefxn_rb_ = core::scoring::ScoreFunctionFactory::create_score_function(
		basic::options::option[ basic::options::OptionKeys::RBSegmentRelax::rb_scorefxn ]()
	);
	core::pose::Pose native_pose;

	// native structure
	if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
		core::import_pose::pose_from_file( native_pose, basic::options::option[ basic::options::OptionKeys::in::file::native ]() , core::import_pose::PDB_file);
		core::pose::set_ss_from_phipsi( native_pose );

#ifdef BOINC_GRAPHICS
		// set native for graphics
		boinc::Boinc::set_graphics_native_pose( native_pose );
#endif

		core::util::switch_to_residue_type_set( native_pose, core::chemical::CENTROID );
	}

}

RBSegmentRelaxImpl::~RBSegmentRelaxImpl(){}

void RBSegmentRelaxImpl::apply( core::pose::Pose & pose ){
	core::pose::set_ss_from_phipsi( pose ); // TODO is this line really necessary? is this not already in import_pose.cc

	// Read RB segs, auto generate loops
	utility::vector1< rbsegment_relax::RBSegment > rbsegs;
	utility::vector1< int > cutpts;
	protocols::loops::Loops loops;
	std::string rbfilename( basic::options::option[ basic::options::OptionKeys::RBSegmentRelax::rb_file ]().name() );

	for ( int i=1; i<=pose.fold_tree().num_cutpoint() ; ++i ) {
		cutpts.push_back( pose.fold_tree().cutpoint(i) );
	}
	int last_peptide_res = pose.total_residue();
	while ( !pose.residue( last_peptide_res ).is_protein() )
			last_peptide_res--;
	read_RBSegment_file( rbsegs, loops, rbfilename, true, last_peptide_res , cutpts  );

	// read fragments
	utility::vector1< core::fragment::FragSetOP > frag_libs;
	if ( basic::options::option[ basic::options::OptionKeys::loops::frag_files ].user() ) {
		protocols::loops::read_loop_fragments( frag_libs );
	}

#ifdef BOINC_GRAPHICS
	// attach boinc graphics pose observer
	protocols::boinc::Boinc::attach_graphics_current_pose_observer( pose );
#endif

	rbsegment_relax::RBSegmentRelax shaker( scorefxn_rb_, rbsegs, loops );
	shaker.initialize( frag_libs );
	shaker.set_randomize( 2 );  //???
	shaker.apply( pose );


}

std::string
RBSegmentRelaxImpl::get_name() const{
	return "RBSegmentRelaxImpl";
}

int
RBSegmentRelax_main() {
	// using namespace rbsegment_relax;
	// using namespace jobdist;
	// using namespace basic::options;
	// using namespace basic::options::OptionKeys;
	// using namespace core::scoring;
	// using namespace core::chemical;
	// using namespace core::id;


	//core::pose::Pose start_pose, pose;
	//core::chemical::ResidueTypeSetCAP rsd_set;

	//std::string pdbfilename;
	//if ( option[ OptionKeys::RBSegmentRelax::input_pdb ].user() )
	// pdbfilename = option[ OptionKeys::RBSegmentRelax::input_pdb ]().name();
	//else
	// pdbfilename = option[ OptionKeys::in::file::s ]()[1];

	// if full-atom load starting structure as full-atom to recover sidechains later
	//if ( option[ in::file::fullatom ]() ) {
	// core::import_pose::pose_from_file( start_pose, pdbfilename , core::import_pose::PDB_file);
	//} else {
	// core::import_pose::centroid_pose_from_pdb( start_pose, pdbfilename , core::import_pose::PDB_file);
	//}

	// roughly guess at secondary structure
	//core::pose::set_ss_from_phipsi( start_pose );


	////////////////////////
	////////////////////////
	// job distributor initialization
	//utility::vector1< protocols::jobdist::BasicJobOP > input_jobs;
	//int const nstruct_flag = option[ out::nstruct ];
	//int const nstruct = std::max( 1, nstruct_flag );
	//protocols::jobdist::BasicJobOP job = new protocols::jobdist::BasicJob("S", "rbseg", nstruct);
	//input_jobs.push_back( job );
	//protocols::jobdist::BaseJobDistributorOP jobdist;

	// output nonidealized silent file or PDBs?
	//bool silent_output;
	//if ( boinc_mode || option[ OptionKeys::out::file::silent ].user() ) {
	// TRb.Debug << "Outputting silent file\n";
	// jobdist = new protocols::jobdist::PlainSilentFileJobDistributor( input_jobs );
	// silent_output = true;
	//} else {
	// TRb.Debug << "Outputting PDBs\n";
	// jobdist = new protocols::jobdist::PlainPdbJobDistributor( input_jobs );
	// silent_output = false;
	//}

	//protocols::jobdist::BasicJobOP prev_job, curr_job;
	//int curr_nstruct;
	//jobdist->startup();


	/////
	/////
	//while ( jobdist->next_job(curr_job, curr_nstruct) ) { // loop over jobs
	// std::string curr_job_tag = curr_job->output_tag( curr_nstruct );

	// pose = start_pose;
	//  if ( option[ in::file::fullatom ]() )
	//   core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID );

	//#ifdef BOINC_GRAPHICS
	//  attach boinc graphics pose observer
	//  protocols::boinc::Boinc::attach_graphics_current_pose_observer( pose );
	//#endif

	// the rigid body movement mover
	//  RBSegmentRelax shaker( scorefxn_rb, rbsegs, loops );
	//  shaker.initialize( frag_libs );
	//  shaker.set_randomize( 2 );  //???
	//  shaker.apply( pose );

	////
	////  output
	//  if ( silent_output ) {
	//   PlainSilentFileJobDistributor *jd =
	//      dynamic_cast< PlainSilentFileJobDistributor * > (jobdist());

	//   std::string silent_struct_type( "binary" );  // default to binary
	//   if ( option[ out::file::silent_struct_type ].user() ) {
	//    silent_struct_type = option[ OptionKeys::out::file::silent_struct_type ];
	//   }

	//   core::io::silent::SilentStructOP ss
	//    = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct( silent_struct_type );

	//   ss->fill_struct( pose, curr_job_tag );

	//   jd->dump_silent( curr_nstruct, *ss );
	//  } else {
	//   jobdist->dump_pose_and_map( curr_job_tag, pose );    // output PDB
	//  }
	// }
	// jobdist->shutdown();

	RBSegmentRelaxImplOP rb_segment_relax_impl( new RBSegmentRelaxImpl() );
	protocols::jd2::JobDistributor::get_instance()->go( rb_segment_relax_impl );

	return 0;
}

//}
} // namespace protocols

////////////////////////////////////////////////////////

