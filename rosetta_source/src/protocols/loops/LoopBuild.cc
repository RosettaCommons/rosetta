// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @brief demo program for implementing loop relax + FA relax
/// @author Srivatsan Raman
/// @author James Thompson
/// @author Mike Tyka
/// @author Daniel J. Mandell

// include these first for building on Visual Studio
#include <protocols/jobdist/JobDistributors.hh>
#include <protocols/jobdist/Jobs.hh>
// AUTO-REMOVED #include <protocols/jobdist/standard_mains.hh>
#include <protocols/moves/MoverStatus.hh>
#include <core/types.hh>
// AUTO-REMOVED #include <core/sequence/util.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>
#include <core/io/pdb/pose_io.hh>
// AUTO-REMOVED #include <basic/options/util.hh>
// AUTO-REMOVED #include <basic/options/after_opts.hh>
// AUTO-REMOVED #include <core/optimization/MinimizerOptions.hh>
// AUTO-REMOVED #include <core/optimization/AtomTreeMinimizer.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>
// AUTO-REMOVED #include <core/scoring/constraints/util.hh>
#include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/chemical/VariantType.hh>
// AUTO-REMOVED
// AUTO-REMOVED #include <core/scoring/rms_util.hh>
#include <core/pose/util.hh>
#include <protocols/evaluation/RmsdEvaluator.hh>
// AUTO-REMOVED #include <core/fragment/FragmentIO.hh>
// AUTO-REMOVED #include <numeric/random/random.hh>
// AUTO-REMOVED #include <numeric/random/random_permutation.hh>
#include <protocols/loops/loops_main.hh>
// AUTO-REMOVED #include <protocols/loops/ccd_closure.hh>
#include <protocols/loops/util.hh>
#include <protocols/loops/Loops.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMover.fwd.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMover.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMover_QuickCCD.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMover_QuickCCD_Moves.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMover_CCD.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMover_KIC.hh>
// AUTO-REMOVED #include <protocols/loops/LoopMoverFactory.hh>
// AUTO-REMOVED #include <protocols/loops/looprelax_protocols.hh>
#include <protocols/loops/LoopBuild.hh>
// AUTO-REMOVED #include <protocols/viewer/viewers.hh>
// AUTO-REMOVED #include <protocols/relax_protocols.hh>
// AUTO-REMOVED #include <protocols/moves/PackRotamersMover.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/task/TaskFactory.hh>

#include <core/io/silent/silent.fwd.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStructFactory.hh>

#include <protocols/moves/symmetry/SetupForSymmetryMover.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>

// AUTO-REMOVED #include <core/fragment/ConstantLengthFragSet.hh>
// AUTO-REMOVED #include <core/fragment/ConstantLengthFragSet.fwd.hh>

// AUTO-REMOVED #include <core/kinematics/visualize.hh>

#include <core/fragment/FragSet.hh>
// AUTO-REMOVED #include <protocols/abinitio/AbrelaxApplication.hh>
// AUTO-REMOVED #include <protocols/abinitio/ClassicAbinitio.hh>
#include <utility/exit.hh>

// AUTO-REMOVED #include <protocols/electron_density/util.hh>
#include <core/scoring/electron_density/util.hh>

#include <protocols/evaluation/PoseEvaluator.hh>
#include <protocols/evaluation/util.hh>

// AUTO-REMOVED #include <protocols/checkpoint/CheckPointer.hh>

// AUTO-REMOVED #include <protocols/comparative_modeling/util.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif

//silly using/typedef

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
// AUTO-REMOVED #include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <protocols/loops/LoopRelaxMover.hh>

#include <core/import_pose/import_pose.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <utility/vector1.hh>
#include <basic/options/option.hh>



using basic::T;
using basic::Error;
using basic::Warning;

namespace protocols {
namespace loops {

int
LoopRelax_main( bool boinc_mode ) {
	basic::Tracer TR("protocols::loopbuild");

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::scoring::constraints;
	using namespace core::chemical;
	using namespace core::id;
	using namespace jobdist;

	std::string remodel            ( option[ OptionKeys::loops::remodel ]() );
	std::string const intermedrelax( option[ OptionKeys::loops::intermedrelax ]() );
	std::string const refine       ( option[ OptionKeys::loops::refine ]() );
	std::string const relax        ( option[ OptionKeys::loops::relax ]() );
	bool const keep_time      ( option[ OptionKeys::loops::timer ]() );

	TR << "==== Loop protocol: ================================================="
		<< std::endl;
	TR << " remodel        " << remodel        << std::endl;
	TR << " intermedrelax  " << intermedrelax  << std::endl;
	TR << " refine         " << refine         << std::endl;
	TR << " relax          " << relax          << std::endl;

	core::pose::Pose start_pose, pose, init_pose_obj;
	core::chemical::ResidueTypeSetCAP rsd_set;

	// DJM: must first load the structure in fullatom mode. Otherwise
	// disulfide-bonded cysteines will not have their residue types set
	// correctly, which causes crashes upon switching to full-atom and repacking.
	if ( option[ in::file::fullatom ]() ) {
		// if full-atom load starting structure as full-atom to recover
		// sidechains later
		rsd_set
			= ChemicalManager::get_instance()->residue_type_set( "fa_standard" );
		core::import_pose::pose_from_pdb(
			start_pose, *rsd_set, option[ OptionKeys::loops::input_pdb ]().name()
		);
	} else { // no full-atom, create a centroid PDB
		rsd_set
			= ChemicalManager::get_instance()->residue_type_set( "centroid" );
		core::import_pose::pose_from_pdb(
			start_pose, *rsd_set, option[ OptionKeys::loops::input_pdb ]().name()
		);
	}

	bool psipred_ss2_ok = loops::set_secstruct_from_psipred_ss2( start_pose );
	if ( !psipred_ss2_ok ) {
		std::string dssp_name( option[ OptionKeys::in::file::dssp ]().name() );
		bool dssp_ok = loops::set_secstruct_from_dssp(start_pose, dssp_name);
		if ( !dssp_ok ) {
			core::pose::set_ss_from_phipsi( start_pose );
		}
	}

	// symmetrize start pose & loopfile
	if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
			protocols::moves::symmetry::SetupForSymmetryMover pre_mover;
			pre_mover.apply( start_pose );
	}

	// bit of a hack for looprelax-into-density
	// set pose for density scoring if a map was input
	// (potentially) dock map into density -- do this here so we only need to
	// dock the pose once
	if ( option[ edensity::mapfile ].user() ) {
		protocols::electron_density::SetupForDensityScoringMover pre_mover;
		pre_mover.apply( start_pose );
	}

	core::pose::Pose native_pose;
	if ( option[ in::file::native ].user() ) {
		core::import_pose::pose_from_pdb(
			native_pose, option[ in::file::native ]()
		);
		core::pose::set_ss_from_phipsi( native_pose );
	} else	{
		native_pose = start_pose;
	}

	// fragment initialization
	// is there any way to clean this up? This logic is very convoluted.
	utility::vector1< core::fragment::FragSetOP > frag_libs;
	if ( remodel == "perturb_ccd" || remodel == "quick_ccd" ||
			remodel == "quick_ccd_moves" || remodel == "old_loop_relax" ||
			remodel == "sdwindow" ||
			option[ OptionKeys::loops::build_initial ].value() ||
	   	( option[ OptionKeys::loops::frag_files ].user()
	      	&& (refine == "refine_ccd" || intermedrelax != "no" || relax != "no")
			)
	) {
		// these protocols optionally take a fragment set .. only load if
		// specified
		read_loop_fragments( frag_libs );
	}

	evaluation::MetaPoseEvaluatorOP evaluator = new evaluation::MetaPoseEvaluator;
	evaluation::read_common_evaluator_options(*evaluator);
	evaluator->add_evaluation(
		new evaluation::SelectRmsdEvaluator( native_pose, "_native" )
	);

	// job distributor initialization
	utility::vector1< BasicJobOP > input_jobs;
	core::Size const nstruct = static_cast< core::Size > (
		std::max( (int) 1, (int)option[ out::nstruct ] )
	);
	BasicJobOP job = new BasicJob("S", "lr", nstruct);
	input_jobs.push_back( job );
	BaseJobDistributorOP jobdist;

	// output nonidealized silent file or PDBs?
	bool silent_output;
	if ( boinc_mode || option[ OptionKeys::out::file::silent ].user() ) {
		TR.Debug << "Outputting non-idealized silent file\n";
		jobdist = new PlainSilentFileJobDistributor( input_jobs );
		silent_output = true;
	} else {
		TR.Debug << "Outputting PDBs\n";
		jobdist = new PlainPdbJobDistributor( input_jobs );
		silent_output = false;
	}

	BasicJobOP curr_job;
	int curr_nstruct;
	jobdist->startup();

	TR << "Annotated sequence of start_pose: "
		<< start_pose.annotated_sequence(true) << std::endl;

	// load loopfile
	protocols::loops::Loops my_loops = protocols::loops::get_loops_from_file();

	while ( jobdist->next_job(curr_job, curr_nstruct) ) { // loop over jobs
		std::string curr_job_tag = curr_job->output_tag( curr_nstruct );
		clock_t starttime = clock();

		pose = start_pose;

		LoopRelaxMover mover;
		mover.frag_libs( frag_libs );
		mover.loops( my_loops );
		mover.relax( relax );
		mover.refine( refine );
		mover.remodel( remodel );
		mover.intermedrelax( intermedrelax );
		mover.set_current_tag( curr_job_tag );

		// add density wts from cmd line to looprelax scorefunctions
		if ( option[ edensity::mapfile ].user() ) {
			ScoreFunctionOP lr_cen_scorefxn = get_cen_scorefxn();
			core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *lr_cen_scorefxn );
			ScoreFunctionOP lr_fa_scorefxn = get_fa_scorefxn();
			core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *lr_fa_scorefxn );
			mover.scorefxns( lr_cen_scorefxn, lr_fa_scorefxn );
		}

		mover.apply( pose );

		//////////////////////////////////////////////////////////////////////////////
		////
		////   Filter
		////
		core::Real final_score( 0.0 );
		using std::string;
		using core::pose::getPoseExtraScores;
		if ( core::pose::getPoseExtraScores(
				pose, std::string("final_looprelax_score"), final_score
			)
		) {
			if ( option[ OptionKeys::loops::final_score_filter ].user() &&
			      final_score > option[ OptionKeys::loops::final_score_filter ]()
			) {
				TR.Debug <<  "FailedFilter " << final_score << " > "
					<< option[  OptionKeys::loops::final_score_filter ]() << std::endl;
				continue;
			}
		}

		/////////////////////////////////////////////////////////////////////////////
		////
		////   Output
		////

		if ( silent_output ) {
			PlainSilentFileJobDistributor *jd =
				dynamic_cast< PlainSilentFileJobDistributor * >
				(jobdist());

			std::string silent_struct_type( "binary" );  // default to binary
			if ( option[ out::file::silent_struct_type ].user() ) {
				silent_struct_type = option[ out::file::silent_struct_type ]();
			}

			using namespace core::io::silent;
			SilentStructOP ss = SilentStructFactory::get_instance()->get_silent_struct(
				silent_struct_type
			);

			ss->fill_struct( pose, curr_job_tag );

			// run PoseEvaluators
			evaluator->apply( pose, curr_job_tag, *ss );

			jd->dump_silent( curr_nstruct, *ss );
		} else {
			if ( ! ( refine == "refine_kic" || remodel == "perturb_kic") ) {
				jobdist->dump_pose_and_map( curr_job_tag, pose );    // output PDB
			}
			else { // DJM: for KIC output, until jobdist matures
				// if closure failed don't output
				if ( remodel == "perturb_kic" ) {
					if ( mover.get_last_move_status() != protocols::moves::MS_SUCCESS ) {
						TR << "Initial kinematic closure failed. Not outputting."
							<< std::endl;
						continue;
					}
				}
				std::string outfile(
					option[ OptionKeys::out::path::path ]().name() + "/" +
					option[ OptionKeys::out::prefix ] +
					"_" + right_string_of(curr_nstruct,4,'0') + ".pdb"
				);

				std::ofstream out( outfile.c_str(), std::ios::out | std::ios::binary );
				core::io::pdb::dump_pdb( pose, out );
				if ( remodel != "no" ) {
					core::Real cen_looprms=0.0;
					getPoseExtraScores( pose, "cen_looprms", cen_looprms );
					out << "loop_cenrms: " << cen_looprms << std::endl;
					// mirror to tracer
					TR << "loop_cenrms: " << cen_looprms << std::endl;
				}
				if ( refine == "refine_kic" ) {
					core::Real final_looprms=0.0;
					core::Real final_score=0.0;
					core::Real final_chainbreak=0.0;
					getPoseExtraScores( pose, "looprms", final_looprms );
					getPoseExtraScores( pose, "final_looprelax_score", final_score );
					getPoseExtraScores( pose, "final_chainbreak", final_chainbreak );
					out << "loop_rms: " << final_looprms << std::endl;
					out << "total_energy: " << final_score << std::endl;
					out << "chainbreak: " << final_chainbreak << std::endl;
					// mirror to tracer
					TR << "loop_rms: " << final_looprms << std::endl;
					TR << "total_energy: " << final_score << std::endl;
					TR << "chainbreak: " << final_chainbreak << std::endl;
				}
				out.flush(); // make sure buffer is flushed before attempting to gzip
				if (option[ OptionKeys::out::pdb_gz ]()) {
					utility::file::gzip( outfile, true );
				}
			}
		}
		clock_t stoptime = clock();
		if ( keep_time ) {
			TR << "Job " << curr_nstruct << " took "<< ((double) stoptime - starttime )/CLOCKS_PER_SEC
			<< " seconds" << std::endl;
		}
	} // loop over jobs
	jobdist->shutdown();

	return 0;
} // Looprelax_main

} // namespace loops
} // namespace protocols
