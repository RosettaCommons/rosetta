// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   demo/ian_test/extract_atomtree_diffs.cc
///
/// @brief  Convert atomtree_diff silent files into normal PDBs.
/// @author Ian Davis (ian.w.davis@gmail.com)

// must be here to avoid VC++ ambiguous symbol w/ ObjexxFCL::byte
// for boinc builds - dek
#include <protocols/jobdist/JobDistributors.hh>


// AUTO-REMOVED #include <numeric/conversions.hh>
#include <numeric/random/random_permutation.hh>
// AUTO-REMOVED #include <numeric/xyzVector.io.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray1D.hh>
// AUTO-REMOVED #include <ObjexxFCL/FArray1.io.hh>
#include <utility/exit.hh>
#include <utility/file/FileName.hh>
#include <utility/io/ozstream.hh>
#include <utility/pointer/owning_ptr.hh>

#include <devel/init.hh>
#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
// AUTO-REMOVED #include <core/pack/task/PackerTask.hh>
// AUTO-REMOVED #include <core/pack/task/TaskFactory.hh>
// AUTO-REMOVED #include <core/scoring/rms_util.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/dunbrack/RotamerConstraint.hh>
// AUTO-REMOVED #include <basic/prof.hh>
#include <basic/Tracer.hh>

// AUTO-REMOVED #include <protocols/docking/DockingInitialPerturbation.hh>
#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh> //for addding constraints if demanded by user
// AUTO-REMOVED #include <protocols/rigid/RB_geometry.hh>
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/ligand_docking/LigandBaseProtocol.hh>
// AUTO-REMOVED #include <protocols/simple_moves/MinMover.hh>
// AUTO-REMOVED #include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/Mover.hh>
// AUTO-REMOVED #include <protocols/moves/MoverContainer.hh>
// AUTO-REMOVED #include <protocols/simple_moves/PackRotamersMover.hh>
// AUTO-REMOVED #include <protocols/rigid/RigidBodyMover.hh>
// AUTO-REMOVED #include <protocols/moves/TrialMover.hh>


// AUTO-REMOVED #include <ctime>
#include <fstream>
#include <set>
#include <sstream>

#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif

// option key includes

#include <basic/options/option_macros.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>
#include <utility/vector0.hh>
#include <utility/vector1.hh>




int
main( int argc, char * argv [] )
{

	try {

	OPT(in::path::database);
	OPT(in::file::extra_res_fa);
	OPT(packing::unboundrot);
	OPT(packing::ex1::ex1);
	OPT(packing::ex1aro::ex1aro);
	OPT(packing::ex2::ex2);
	OPT(packing::extrachi_cutoff);
	OPT(packing::no_optH);
	OPT(packing::flip_HNQ);
	OPT(docking::ligand::soft_rep);
	OPT(docking::ligand::old_estat);

	OPT(in::file::s);
	OPT(in::file::tags);
	OPT(out::path::pdb);
	OPT(enzdes::cstfile);

	using basic::options::option;
	using namespace basic::options::OptionKeys;
	using namespace protocols;
	using namespace protocols::jobdist;
	using namespace protocols::moves;
	basic::Tracer TR( "extract_atomtree_diffs.main" );

	// Parses command line options and inits RNG.
	// Doesn't seem to hurt to do it again if already done once (?)
	devel::init(argc, argv);

	time_t overall_start_time = time(NULL);
	utility::vector1< BasicJobOP > input_jobs = load_s_and_l();
	// Reduce read contention between processes by randomizing the order in which structures are processed
	numeric::random::random_permutation( input_jobs, numeric::random::rg() );
	PlainPdbJobDistributor jobdist( input_jobs );

	// A "native" pose for the diff reference point.
	// This used to be required but will be rarely used now.
	core::pose::PoseOP native_pose;
	if( option[ in::file::native ].user() ) {
		native_pose = core::pose::PoseOP( new core::pose::Pose() );
		core::import_pose::pose_from_pdb( *native_pose, option[ in::file::native ]().name() );
	}

	// Tags to be extracted.  If empty, extract top 5% by total score.
	std::set< std::string > desired_tags;
	if( option[ in::file::tags ].active() ) {
		desired_tags.insert( option[ in::file::tags ]().begin(), option[ in::file::tags ]().end() );
	}

	// Reading the enzdes constraints creates new residue types we may need, in the case of covalent constraints.
	// Actually, the residue types aren't created until we apply the constraints to the reference poses, below.
	protocols::toolbox::match_enzdes_util::EnzConstraintIOOP constraint_io = NULL;
	if( option[basic::options::OptionKeys::enzdes::cstfile].user() ){
		//we need the residue type set, assuming FA standard is used
		core::chemical::ResidueTypeSetCAP restype_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
		option[basic::options::OptionKeys::run::preserve_header ].value(true);
		constraint_io = protocols::toolbox::match_enzdes_util::EnzConstraintIOOP( new protocols::toolbox::match_enzdes_util::EnzConstraintIO( restype_set ) );
		constraint_io->read_enzyme_cstfile(basic::options::option[basic::options::OptionKeys::enzdes::cstfile]);
	}

	// Create a protocol object just to steal its score function.
	core::scoring::ScoreFunctionOP scorefxn = protocols::ligand_docking::LigandBaseProtocol().get_scorefxn();

	BasicJobOP curr_job;
	int curr_nstruct, num_structures_processed = 0;
	jobdist.startup();
	while( jobdist.next_job(curr_job, curr_nstruct) ) {
		TR << "Reading silent file " << curr_job->input_tag() << " ... ";
		core::import_pose::atom_tree_diffs::AtomTreeDiff atdiff( curr_job->input_tag() );
		core::import_pose::atom_tree_diffs::ScoresPairList const & scores_list = atdiff.scores();
		TR << scores_list.size() << " structures" << std::endl;

		if( option[basic::options::OptionKeys::enzdes::cstfile].user() ){
			utility::vector1< core::pose::PoseOP > const & ref_poses = atdiff.all_ref_poses();
			for( core::Size i = 1; i <= ref_poses.size(); ++i ) {
				constraint_io->add_constraints_to_pose( *ref_poses[i], scorefxn, false );
			}
		}

		//if( !option[ in::file::tags ].active() ) {
		//	desired_tags.clear();
		//	protocols::ligand_docking::select_best_poses(atdiff, desired_tags);
		//	TR << "Keeping " << desired_tags.size() << " best structures\n";
		//}

		// By iterating in the original order and testing against desired_tags,
		// we're always seeking forward on the disk, which I think should help performance.
		for(core::Size ii = 1; ii <= scores_list.size(); ++ii) {
			std::string output_tag = scores_list[ii].first;
			//std::map< std::string, core::Real > scores = scores_list[ii].second;

			// If this tag was not requested on cmd line, skip it.
			if( !desired_tags.empty() && desired_tags.count(output_tag) == 0 ) continue;

			core::pose::PoseOP the_pose( new core::pose::Pose() );
			if( native_pose == NULL ) atdiff.read_pose(output_tag, *the_pose);
			else atdiff.read_pose(output_tag, *the_pose, *native_pose);

			// Score new structure.  Cached energies (including *residue* energies)
			// must be up-to-date in order to get sensible output.  If you remove these
			// lines, you *must* insert equivalent logic at the end of all apply() methods
			// (or at least for all movers that might be passed to this function).
			core::pack::dunbrack::load_unboundrot( *the_pose ); // adds scoring bonuses for the "unbound" rotamers, if any
			(*scorefxn)( *the_pose );
			/// Now handled automatically.  scorefxn->accumulate_residue_total_energies( *the_pose );

			jobdist.dump_pose_and_map(output_tag, *the_pose);

			num_structures_processed += 1;
			TR << "Finished " << output_tag  << std::endl;
		}
	} // loop over jobs and nstructs
	jobdist.shutdown();

	time_t overall_end_time = time(NULL);
	TR << "Finished all " << num_structures_processed << " structures in " << (overall_end_time - overall_start_time) << " seconds." << std::endl;
	if ( num_structures_processed == 0 )
		basic::Warning() << "No structures processed.  Existing output files may have been skipped, did you mean to delete them?" << std::endl;
	return 0;

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

