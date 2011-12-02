// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/apps/public/scenarios/ca_to_allatom.cc
/// @brief An app which constructs an all atom model from a CA-only input trace.
///        The inputs are the starting CA trace, a list of residue ranges corresponding
///        to secondary structure elements (as best as can be inferred from the initial trace)
///        and a set of backbone torsion fragments.
/// @author Frank DiMaio

// libRosetta headers

#include <protocols/jobdist/JobDistributors.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>

#include <core/types.hh>

#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
// AUTO-REMOVED
// AUTO-REMOVED #include <core/conformation/util.hh>
// AUTO-REMOVED #include <core/chemical/ResidueSelector.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
// AUTO-REMOVED #include <core/optimization/AtomTreeMinimizer.hh>
// AUTO-REMOVED #include <core/optimization/MinimizerOptions.hh>

// AUTO-REMOVED #include <core/scoring/sasa.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>
//#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
// AUTO-REMOVED #include <core/scoring/constraints/AtomPairConstraint.hh>
// AUTO-REMOVED #include <core/scoring/constraints/Func.hh>
// AUTO-REMOVED #include <core/scoring/constraints/HarmonicFunc.hh>
#include <core/scoring/electron_density/util.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/pose/util.hh>

#include <devel/init.hh>

#include <basic/options/option.hh>
#include <basic/options/after_opts.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/RBSegmentRelax.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// AUTO-REMOVED #include <basic/basic.hh>
// AUTO-REMOVED #include <basic/database/open.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentStructFactory.hh>
#include <core/id/AtomID.hh>

#include <core/io/silent/silent.fwd.hh>
// AUTO-REMOVED #include <core/io/silent/BinaryProteinSilentStruct.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>
// AUTO-REMOVED #include <core/scoring/constraints/CoordinateConstraint.hh>

#include <protocols/rigid/RBSegmentRelax.hh>
//#include <protocols/rbsegment_Moves/RBSegmentMover.hh>
//#include <protocols/rbsegment_Moves/RBSegmentRelax.hh>
//#include <protocols/rbsegment_Moves/FragInsertAndAlignMover.hh>

#include <protocols/loops/loops_main.hh>
// AUTO-REMOVED #include <protocols/loops/ccd_closure.hh>
#include <protocols/loops/Loops.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMover.fwd.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMover.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMover_QuickCCD.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMover_QuickCCD_Moves.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMover_CCD.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMover_KIC.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMoverFactory.hh>
// AUTO-REMOVED #include <protocols/loops/LoopRelaxMover.hh>
// AUTO-REMOVED #include <protocols/loops/looprelax_protocols.hh>
// AUTO-REMOVED #include <protocols/loops/LoopBuild.hh>

// AUTO-REMOVED #include <protocols/relax_protocols.hh>
#include <utility/options/OptionCollection.hh>
#include <protocols/simple_filters/RmsdEvaluator.hh>
#include <protocols/evaluation/EvaluatorFactory.hh>
#include <protocols/viewer/viewers.hh>
//#include <protocols/moves/BackboneMover.hh>

#include <utility/vector1.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>

#include <basic/options/option_macros.hh>

// Numeric Headers
#include <numeric/random/random.hh>

// Platform Headers
#include <platform/types.hh>

static numeric::random::RandomGenerator RG(86759999);

#include <basic/Tracer.hh>
using basic::T;
using basic::Warning;
using basic::Error;

// tracer
static basic::Tracer TZ("pilot_apps::ca_to_allatom");

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

#include <core/fragment/FragSet.hh>
#include <core/import_pose/import_pose.hh>
#include <protocols/evaluation/util.hh>
#include <protocols/moves/MoverStatistics.hh>
#include <protocols/rigid/RBSegment.hh>
#include <numeric/xyzVector.io.hh>

//Auto Headers
#include <core/kinematics/FoldTree.hh>

// add options
OPT_1GRP_KEY( Real, ca_to_allatom, frag_randomness )
OPT_1GRP_KEY( Boolean, ca_to_allatom, no_lr )
OPT_1GRP_KEY( Boolean, ca_to_allatom, fix_ligands )

// recenters all non-randomized atoms in the pose
// applies the same trasformation to an alternate pose (if given, e.g. from constraints)
numeric::xyzVector< core::Real > recenter_with_missing( core::pose::Pose &pose) {
	numeric::xyzVector< core::Real > cog(0,0,0), x_i;

	core::Size count=0;
	for( core::Size ir = 1; ir <= pose.total_residue(); ir++){
		if (!pose.residue(ir).is_protein() ) continue;
		x_i = pose.xyz( core::id::AtomID( 2, ir ) );
		if (x_i.length() < 800 ) { // was this randomized??
			cog += x_i;
			count++;
		}
	}

	if (count > 0)
		cog /= (core::Real) count;

	for( core::Size ir = 1; ir <= pose.total_residue(); ir++){
		for( core::Size at = 1; at <= pose.residue( ir ).natoms(); at++){
			pose.set_xyz(  core::id::AtomID( at, ir ),  pose.xyz( core::id::AtomID( at, ir )) - cog  );
		}
	}
	TZ << "Transforming by " << (-cog) << std::endl;

	return (-cog);
}




//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////


void*
ca_to_allatom_main( void * )
{
	using namespace protocols;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace utility::file;

	// RB scoring function; get from command line
	core::scoring::ScoreFunctionOP scorefxn_cst
	         = core::scoring::ScoreFunctionFactory::create_score_function(  option[ RBSegmentRelax::rb_scorefxn ]() );

	// grab edens scores from CL as well
	if ( option[ edensity::mapfile ].user() ) {
		core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *scorefxn_cst );
	}

	// Native pose
	core::pose::Pose native_pose, pose, start_pose;
	if ( option[ in::file::native ].user() ) {
		core::import_pose::pose_from_pdb(
			native_pose, option[ in::file::native ]()
		);
		core::pose::set_ss_from_phipsi( native_pose );
	}

	// load rbsegs
	utility::vector1< protocols::rigid::RBSegment > rbsegs,rbsegs_remap;
	utility::vector1< int > cutpts;
	protocols::loops::Loops loops;
	std::string filename( option[ OptionKeys::RBSegmentRelax::rb_file ]().name() );

	evaluation::MetaPoseEvaluatorOP evaluator = new evaluation::MetaPoseEvaluator;
	evaluation::EvaluatorFactory::get_instance()->add_all_evaluators(*evaluator);
	evaluator->add_evaluation( new simple_filters::SelectRmsdEvaluator( native_pose, "_native" ) );

	utility::vector1< protocols::jobdist::BasicJobOP > input_jobs = protocols::jobdist::load_s_and_l();

	protocols::jobdist::BaseJobDistributorOP jobdist;
	bool const silent_output = option[ out::file::silent ].user();
	if ( silent_output ) {
		TZ << "Silent Output Mode " << std::endl;
		jobdist = new protocols::jobdist::PlainSilentFileJobDistributor(input_jobs);
	} else {
		TZ << "PDB Output Mode " << std::endl;
		jobdist = new protocols::jobdist::PlainPdbJobDistributor(input_jobs, "none");
	}

	if( option[ out::nooutput ]() ){
		jobdist->disable_output();
		jobdist->enable_ignorefinished();
	}

	// read fragments
	utility::vector1< core::fragment::FragSetOP > frag_libs;
	bool hasLoopFile = basic::options::option[ basic::options::OptionKeys::loops::frag_files ].user();
	if ( hasLoopFile )
		protocols::loops::read_loop_fragments( frag_libs );

	protocols::jobdist::BasicJobOP curr_job, prev_job;
	int curr_nstruct;
	jobdist->startup();
	while( jobdist->next_job(curr_job, curr_nstruct) ) {
		std::string curr_job_tag = curr_job->output_tag( curr_nstruct );

		// read as-needed
		if ( !prev_job || curr_job->input_tag() != prev_job->input_tag() ) {
			core::import_pose::pose_from_pdb( start_pose, curr_job->input_tag() );

			for (int i=1; i<=start_pose.fold_tree().num_cutpoint() ; ++i)
				cutpts.push_back( start_pose.fold_tree().cutpoint(i) );
			int last_peptide_res = start_pose.total_residue();
			while ( !start_pose.residue( last_peptide_res ).is_protein() ) last_peptide_res--;

			std::string rbfilename( option[ OptionKeys::RBSegmentRelax::rb_file ]().name() );
			protocols::rigid::read_RBSegment_file( rbsegs, loops, rbfilename, true, last_peptide_res , cutpts  );
		}
		pose = start_pose;

		// the rigid body movement mover
		protocols::rigid::RBSegmentRelax shaker( scorefxn_cst, rbsegs, loops );
		shaker.initialize( frag_libs , option[ ca_to_allatom::frag_randomness ]() );
		shaker.set_bootstrap( true );
		shaker.set_skip_lr( option[ ca_to_allatom::no_lr ] );
		shaker.set_fix_ligands( option[ ca_to_allatom::fix_ligands ] );
		shaker.apply( pose );

		////
		////  output
		if ( silent_output ) {
			protocols::jobdist::PlainSilentFileJobDistributor *jd =
					 dynamic_cast< protocols::jobdist::PlainSilentFileJobDistributor * > (jobdist());

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

		prev_job = curr_job;
	} // loop over jobs
	jobdist->shutdown();

	return 0;
}


int
main( int argc, char * argv [] )
{
	// options, random initialization
	NEW_OPT( ca_to_allatom::frag_randomness, "fragment randomness", 0.0 );
	NEW_OPT( ca_to_allatom::no_lr,           "skip lr?", false );
	NEW_OPT( ca_to_allatom::fix_ligands,     "fix ligands?", false );

	devel::init( argc, argv );
	protocols::viewer::viewer_main( ca_to_allatom_main );

	return 0;
}
