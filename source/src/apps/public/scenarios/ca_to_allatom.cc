// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/apps/public/scenarios/ca_to_allatom.cc
/// @brief An app which constructs an all atom model from a CA-only input trace.
///        The inputs are the starting CA trace, a list of residue ranges corresponding
///        to secondary structure elements (as best as can be inferred from the initial trace)
///        and a set of backbone torsion fragments.
/// @author Frank DiMaio

// libRosetta headers


#include <protocols/jd2/JobDistributor.hh>

#include <core/types.hh>

#include <core/conformation/Residue.hh>

#include <core/chemical/ResidueTypeSet.fwd.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>

#include <core/scoring/electron_density/util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/util.hh>

#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/RBSegmentRelax.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

#include <core/io/silent/SilentStructFactory.hh>
#include <core/id/AtomID.hh>

#include <core/io/silent/silent.fwd.hh>

#include <protocols/rbsegment_relax/RBSegmentRelax.hh>

#include <protocols/loops/loops_main.hh>
#include <protocols/loops/Loops.hh>


#include <utility/options/OptionCollection.hh>
#include <protocols/simple_filters/RmsdEvaluator.hh>
#include <protocols/evaluation/EvaluatorFactory.hh>
#include <protocols/viewer/viewers.hh>

#include <utility/vector1.hh>

#include <basic/options/option_macros.hh>

// Numeric Headers
#include <numeric/random/random.hh>

// Platform Headers
#include <platform/types.hh>


#include <basic/Tracer.hh>
using basic::T;
using basic::Warning;
using basic::Error;

// tracer
static THREAD_LOCAL basic::Tracer TZ( "pilot_apps::ca_to_allatom" );

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

#include <core/fragment/FragSet.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/evaluation/util.hh>
#include <protocols/moves/MoverStatistics.hh>
#include <protocols/rbsegment_relax/RBSegment.hh>
#include <numeric/xyzVector.io.hh>

//Auto Headers
#include <core/kinematics/FoldTree.hh>

#include <utility/excn/Exceptions.hh>

// add options
OPT_1GRP_KEY( Real, ca_to_allatom, frag_randomness )
OPT_1GRP_KEY( Boolean, ca_to_allatom, no_lr )
OPT_1GRP_KEY( Boolean, ca_to_allatom, fix_ligands )

// recenters all non-randomized atoms in the pose
// applies the same trasformation to an alternate pose (if given, e.g. from constraints)
numeric::xyzVector< core::Real > recenter_with_missing( core::pose::Pose &pose) {
	numeric::xyzVector< core::Real > cog(0,0,0), x_i;

	core::Size count=0;
	for ( core::Size ir = 1; ir <= pose.total_residue(); ir++ ) {
		if ( !pose.residue(ir).is_protein() ) continue;
		x_i = pose.xyz( core::id::AtomID( 2, ir ) );
		if ( x_i.length() < 800 ) { // was this randomized??
			cog += x_i;
			count++;
		}
	}

	if ( count > 0 ) {
		cog /= (core::Real) count;
	}

	for ( core::Size ir = 1; ir <= pose.total_residue(); ir++ ) {
		for ( core::Size at = 1; at <= pose.residue( ir ).natoms(); at++ ) {
			pose.set_xyz(  core::id::AtomID( at, ir ),  pose.xyz( core::id::AtomID( at, ir )) - cog  );
		}
	}
	TZ << "Transforming by " << (-cog) << std::endl;

	return (-cog);
}

//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
class CaToAllAtom : public protocols::moves::Mover{
public:
	CaToAllAtom();
	virtual ~CaToAllAtom();
	void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;
private:
	core::scoring::ScoreFunctionOP scorefxn_cst_;
	core::pose::Pose native_pose_;

};

typedef utility::pointer::shared_ptr< CaToAllAtom > CaToAllAtomOP;

CaToAllAtom::CaToAllAtom(){
	// RB scoring function; get from command line
	scorefxn_cst_ = core::scoring::ScoreFunctionFactory::create_score_function(
		basic::options::option[ basic::options::OptionKeys::RBSegmentRelax::rb_scorefxn ]() );

	// grab edens scores from CL as well
	if ( basic::options::option[ basic::options::OptionKeys::edensity::mapfile ].user() ) {
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *scorefxn_cst_ );
	}

	// Native pose
	if ( basic::options::option[ basic::options::OptionKeys::in::file::native ].user() ) {
		core::import_pose::pose_from_file( native_pose_, basic::options::option[ basic::options::OptionKeys::in::file::native ]() , core::import_pose::PDB_file);
		core::pose::set_ss_from_phipsi( native_pose_ ); /// Is this necessary? Done by import pose?
	}
}

CaToAllAtom::~CaToAllAtom(){}

void CaToAllAtom::apply( core::pose::Pose & pose ){
	// load rbsegs
	utility::vector1< protocols::rbsegment_relax::RBSegment > rbsegs,rbsegs_remap;
	utility::vector1< int > cutpts;
	protocols::loops::Loops loops;
	//std::string filename( basic::options::option[ basic::options::OptionKeys::RBSegmentRelax::rb_file ]().name() );

	using protocols::evaluation::PoseEvaluatorOP;
	protocols::evaluation::MetaPoseEvaluatorOP evaluator( new protocols::evaluation::MetaPoseEvaluator );
	protocols::evaluation::EvaluatorFactory::get_instance()->add_all_evaluators(*evaluator);
	evaluator->add_evaluation( PoseEvaluatorOP( new protocols::simple_filters::SelectRmsdEvaluator( native_pose_, "_native" ) ) );

	utility::vector1< core::fragment::FragSetOP > frag_libs;
	bool hasLoopFile = basic::options::option[ basic::options::OptionKeys::loops::frag_files ].user();
	if ( hasLoopFile ) {
		protocols::loops::read_loop_fragments( frag_libs );
	}

	for ( int i=1; i<=pose.fold_tree().num_cutpoint() ; ++i ) {
		cutpts.push_back( pose.fold_tree().cutpoint(i) );
	}
	int last_peptide_res = pose.total_residue();
	while ( !pose.residue( last_peptide_res ).is_protein() ) last_peptide_res--;

	std::string rbfilename( basic::options::option[ basic::options::OptionKeys::RBSegmentRelax::rb_file ]().name() );
	protocols::rbsegment_relax::read_RBSegment_file( rbsegs, loops, rbfilename, true, last_peptide_res , cutpts  );

	protocols::rbsegment_relax::RBSegmentRelax shaker( scorefxn_cst_, rbsegs, loops );
	shaker.initialize( frag_libs , basic::options::option[ basic::options::OptionKeys::ca_to_allatom::frag_randomness ]() );
	shaker.set_bootstrap( true );
	shaker.set_skip_lr( basic::options::option[ basic::options::OptionKeys::ca_to_allatom::no_lr ] );
	shaker.set_fix_ligands( basic::options::option[ basic::options::OptionKeys::ca_to_allatom::fix_ligands ] );
	shaker.apply( pose );
}

std::string CaToAllAtom::get_name() const{
	return "CaToAllAtom";
}

void*
ca_to_allatom_main( void * )
{
	using namespace protocols;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility::file;

	// // load rbsegs
	// utility::vector1< protocols::rbsegment_relax::RBSegment > rbsegs,rbsegs_remap;
	// utility::vector1< int > cutpts;
	// protocols::loops::Loops loops;
	// std::string filename( option[ OptionKeys::RBSegmentRelax::rb_file ]().name() );
	//
	// evaluation::MetaPoseEvaluatorOP evaluator = new evaluation::MetaPoseEvaluator;
	// evaluation::EvaluatorFactory::get_instance()->add_all_evaluators(*evaluator);
	// evaluator->add_evaluation( new simple_filters::SelectRmsdEvaluator( native_pose_, "_native" ) );

	// utility::vector1< protocols::jobdist::BasicJobOP > input_jobs = protocols::jobdist::load_s_and_l();

	// protocols::jobdist::BaseJobDistributorOP jobdist;
	// bool const silent_output = option[ out::file::silent ].user();
	// if ( silent_output ) {
	//  TZ << "Silent Output Mode " << std::endl;
	//  jobdist = new protocols::jobdist::PlainSilentFileJobDistributor(input_jobs);
	// } else {
	//  TZ << "PDB Output Mode " << std::endl;
	//  jobdist = new protocols::jobdist::PlainPdbJobDistributor(input_jobs, "none");
	// }

	// if( option[ out::nooutput ]() ){
	//  jobdist->disable_output();
	//  jobdist->enable_ignorefinished();
	// }

	// read fragments
	// utility::vector1< core::fragment::FragSetOP > frag_libs;
	// bool hasLoopFile = basic::options::option[ basic::options::OptionKeys::loops::frag_files ].user();
	// if ( hasLoopFile )
	//  protocols::loops::read_loop_fragments( frag_libs );

	// protocols::jobdist::BasicJobOP curr_job, prev_job;
	// int curr_nstruct;
	// jobdist->startup();
	// while( jobdist->next_job(curr_job, curr_nstruct) ) {
	//  std::string curr_job_tag = curr_job->output_tag( curr_nstruct );
	//
	//  // read as-needed
	//  if ( !prev_job || curr_job->input_tag() != prev_job->input_tag() ) {
	//   core::import_pose::pose_from_file( start_pose, curr_job->input_tag() , core::import_pose::PDB_file);

	//   for (int i=1; i<=start_pose.fold_tree().num_cutpoint() ; ++i)
	//    cutpts.push_back( start_pose.fold_tree().cutpoint(i) );
	//   int last_peptide_res = start_pose.total_residue();
	//   while ( !start_pose.residue( last_peptide_res ).is_protein() ) last_peptide_res--;
	//
	//   std::string rbfilename( option[ OptionKeys::RBSegmentRelax::rb_file ]().name() );
	//   protocols::rbsegment_relax::read_RBSegment_file( rbsegs, loops, rbfilename, true, last_peptide_res , cutpts  );
	//  }
	//  pose = start_pose;

	// the rigid body movement mover
	//  protocols::rbsegment_relax::RBSegmentRelax shaker( scorefxn_cst, rbsegs, loops );
	//  shaker.initialize( frag_libs , option[ ca_to_allatom::frag_randomness ]() );
	//  shaker.set_bootstrap( true );
	//  shaker.set_skip_lr( option[ ca_to_allatom::no_lr ] );
	//  shaker.set_fix_ligands( option[ ca_to_allatom::fix_ligands ] );
	//  shaker.apply( pose );

	////
	////  output
	//  if ( silent_output ) {
	//   protocols::jobdist::PlainSilentFileJobDistributor *jd =
	//      dynamic_cast< protocols::jobdist::PlainSilentFileJobDistributor * > (jobdist());
	//
	//   std::string silent_struct_type( "binary" );  // default to binary
	//   if ( option[ out::file::silent_struct_type ].user() ) {
	//    silent_struct_type = option[ OptionKeys::out::file::silent_struct_type ];
	//   }
	//
	//   core::io::silent::SilentStructOP ss
	//    = core::io::silent::SilentStructFactory::get_instance()->get_silent_struct( silent_struct_type );
	//
	//   ss->fill_struct( pose, curr_job_tag );
	//
	//   jd->dump_silent( curr_nstruct, *ss );
	//  } else {
	//   jobdist->dump_pose_and_map( curr_job_tag, pose );    // output PDB
	//  }
	//
	//  prev_job = curr_job;
	// } // loop over jobs
	// jobdist->shutdown();

	CaToAllAtomOP ca_to_all_atom( new CaToAllAtom() );
	protocols::jd2::JobDistributor::get_instance()->go( ca_to_all_atom );

	return 0;
}


int
main( int argc, char * argv [] )
{
	try {
		// options, random initialization
		NEW_OPT( ca_to_allatom::frag_randomness, "fragment randomness", 0.0 );
		NEW_OPT( ca_to_allatom::no_lr,           "skip lr?", false );
		NEW_OPT( ca_to_allatom::fix_ligands,     "fix ligands?", false );

		devel::init( argc, argv );
		protocols::viewer::viewer_main( ca_to_allatom_main );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}
