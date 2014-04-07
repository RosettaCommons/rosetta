// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   apps/pilot/rmoretti/ligand_dock_jd1.cc
///
/// @brief Old, JD1 version of ligand_dock, kept in case someone really needs it.
/// @author Ian Davis (ian.w.davis@gmail.com)

// must be here to avoid VC++ ambiguous symbol w/ ObjexxFCL::byte
// for boinc builds - dek
#include <protocols/jobdist/JobDistributors.hh>


#include <core/types.hh>
// AUTO-REMOVED #include <core/io/pdb/pose_io.hh>
#include <basic/options/option.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>
// AUTO-REMOVED #include <basic/prof.hh>
#include <basic/Tracer.hh>
#include <basic/datacache/CacheableString.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/ScoreFunction.hh>

#include <protocols/toolbox/match_enzdes_util/EnzConstraintIO.hh> //for addding constraints if demanded by user
#include <protocols/jobdist/Jobs.hh>
#include <protocols/jobdist/standard_mains.hh>
#include <protocols/ligand_docking/LigandDockProtocol.hh>
#include <protocols/moves/Mover.hh>

#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>
#include <utility/file/FileName.hh>
#include <utility/io/ozstream.hh>

#if defined(WIN32) || defined(__CYGWIN__)
	#include <ctime>
#endif

#include <fstream>

// ObjexxFCL headers
#include <ObjexxFCL/string.functions.hh>


// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/enzdes.OptionKeys.gen.hh>

#include <core/import_pose/import_pose.hh>
#include <core/pose/Pose.hh>


using namespace ObjexxFCL;

//////////////////////////////////////////////////////////////////////////////
/// @details Assumes option system has already been initialized!
int
ligand_dock_main_jd1()
{
	using basic::options::option;
	using namespace basic::options::OptionKeys;
	using namespace protocols;
	using namespace protocols::jobdist;
	using namespace protocols::moves;
	using namespace protocols::ligand_docking;
	basic::Tracer TR("ligand_dock.main");

	// Build overall docking protocol Mover
	LigandDockProtocolOP dockingProtocol = new LigandDockProtocol();

	time_t overall_start_time = time(NULL);
	bool const use_silent_in = option[ in::file::silent ].active();
	utility::vector1< BasicJobOP > input_jobs;
	core::import_pose::atom_tree_diffs::AtomTreeDiffOP atdiff;
	if( use_silent_in ) {
		using core::import_pose::atom_tree_diffs::AtomTreeDiff;
		int const nstruct = std::max( 1, option[ out::nstruct ]() );
		atdiff = new AtomTreeDiff( *(option[ in::file::silent ]().begin()) );
		if( (option [ in::file::silent ]()).size() > 1 ) utility_exit_with_message("ligand_dock.main can only handle one input silent file at a time!");
		if( option[ in::file::tags ].active() ) {
			// Do only the ones specified in tags.  If it's a number, try it as an index too.
			utility::vector1< std::string > tags = option[ in::file::tags ]();
			for(core::Size i = 1; i <= tags.size(); ++i) {
				core::Size tag_as_int = 0;
				if( is_int(tags[i]) ) tag_as_int = int_of(tags[i]);
				if( atdiff->has_tag( tags[i] ) ) {
					input_jobs.push_back(new BasicJob( tags[i], tags[i], nstruct ));
				} else if( 0 < tag_as_int && tag_as_int <= atdiff()->scores().size() ) {
					input_jobs.push_back(new BasicJob( atdiff()->scores()[tag_as_int].first, atdiff()->scores()[tag_as_int].first, nstruct ));
				} else {
					utility_exit_with_message("Can't find tag "+tags[i]);
				}
			}
		} else {
			// Just do them all!
			for(core::Size i = 1; i <= atdiff()->scores().size(); ++i) {
				input_jobs.push_back(new BasicJob( atdiff()->scores()[i].first, atdiff()->scores()[i].first, nstruct ));
			}
		}
	} else {
		input_jobs = load_s_and_l();
		// Reduce read contention between processes by randomizing the order in which structures are processed
		// Don't want to do this with silent file input -- slows access and screws up reference structure tracking.
		numeric::random::random_permutation( input_jobs, numeric::random::RG );
	}
	std::string outfile = "silent.out";
	if( option[ out::file::silent].user() ) outfile = option[ out::file::silent];
	AtomTreeDiffJobDistributor jobdist( input_jobs, outfile );
	jobdist.set_precision(6,3,1); // makes silent file much smaller, ~3x vs. default 6,4,2

	// A "native" pose for calculating RMS (optional)
	core::pose::PoseOP rms_native_pose;
	if( option[ in::file::native ].user() ) {
		rms_native_pose = new core::pose::Pose();
		core::import_pose::pose_from_pdb( *rms_native_pose, option[ in::file::native ]().name() );
	}

	protocols::toolbox::match_enzdes_util::EnzConstraintIOOP constraint_io = NULL;
	if( option[basic::options::OptionKeys::enzdes::cstfile].user() ){
		//we need the residue type set, assuming FA standard is used
		core::chemical::ResidueTypeSetCAP restype_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
		option[basic::options::OptionKeys::run::preserve_header ].value(true);
		constraint_io = new protocols::toolbox::match_enzdes_util::EnzConstraintIO( restype_set );
		constraint_io->read_enzyme_cstfile(basic::options::option[basic::options::OptionKeys::enzdes::cstfile]);
		core::scoring::ScoreFunctionOP scorefunc = dockingProtocol->get_scorefxn();
		if (scorefunc->has_zero_weight( core::scoring::coordinate_constraint ) ) scorefunc->set_weight(core::scoring::coordinate_constraint, 1.0 ) ;
		if (scorefunc->has_zero_weight( core::scoring::atom_pair_constraint ) ) scorefunc->set_weight(core::scoring::atom_pair_constraint, 1.0 ) ;
		if (scorefunc->has_zero_weight( core::scoring::angle_constraint ) ) scorefunc->set_weight(core::scoring::angle_constraint, 1.0 ) ;
		if (scorefunc->has_zero_weight( core::scoring::dihedral_constraint ) ) scorefunc->set_weight(core::scoring::dihedral_constraint, 1.0 ) ;
	}

	BasicJobOP curr_job, prev_job;
	int curr_nstruct = 0, num_structures_processed = 0;
	core::pose::PoseOP input_pose; // starts NULL, coords *never* modified!
	jobdist.startup();
	while( jobdist.next_job(curr_job, curr_nstruct) ) {
		time_t pdb_start_time = time(NULL);
		TR << "Starting " << curr_job->output_tag(curr_nstruct) << " ..." << std::endl;

		// we read each PDB just once to save on disk I/O
		if( curr_job.get() != prev_job.get() || input_pose.get() == NULL ) {
			input_pose = new core::pose::Pose();
			if( use_silent_in ) atdiff->read_pose( curr_job->input_tag(), *input_pose );
			else core::import_pose::pose_from_pdb( *input_pose, curr_job->input_tag() );

			//if constraints are requested
			if( option[basic::options::OptionKeys::enzdes::cstfile].user() ){
				constraint_io->add_constraints_to_pose( *input_pose, dockingProtocol->get_scorefxn(), false );
			}
		}

		// Make a modifiable copy of the pose read from disk
		core::pose::PoseOP the_pose = new core::pose::Pose( *input_pose );
		the_pose->data().set(core::pose::datacache::CacheableDataType::JOBDIST_OUTPUT_TAG, new basic::datacache::CacheableString(curr_job->output_tag(curr_nstruct)));
		dockingProtocol->apply( *the_pose );

		// Score new structure and add to silent file
		std::map< std::string, core::Real > scores;
		core::import_pose::atom_tree_diffs::map_of_weighted_scores(*the_pose, *(dockingProtocol->get_scorefxn()), scores);
		if( rms_native_pose.get() == NULL ) { // input serves as native
			dockingProtocol->append_ligand_docking_scores(*input_pose, *the_pose,
				dockingProtocol->get_scorefxn(), scores, constraint_io);
		} else {
			dockingProtocol->append_ligand_docking_scores(*rms_native_pose, *the_pose,
				dockingProtocol->get_scorefxn(), scores, constraint_io);
		}
		// Want to recycle the reference poses from the input silent file, or else every entry becomes a new reference pose!
		if( use_silent_in ) {
			//std::cout << "Ref addr " << (atdiff->ref_pose_for( curr_job->input_tag() ))() << std::endl;
			jobdist.dump_pose( curr_job->output_tag(curr_nstruct), scores, *(atdiff->ref_pose_for( curr_job->input_tag() )), *the_pose );
		}
		else jobdist.dump_pose( curr_job->output_tag(curr_nstruct), scores, *input_pose, *the_pose );

		prev_job = curr_job; // pointer assignment, not a copy op
		num_structures_processed += 1;
		time_t pdb_end_time = time(NULL);
		TR << "Finished " << curr_job->output_tag(curr_nstruct) << " in " << (long)(pdb_end_time - pdb_start_time) << " seconds." << std::endl;
	} // loop over jobs and nstructs

	time_t overall_end_time = time(NULL);
	TR << "Finished all " << num_structures_processed << " structures in " << (long)(overall_end_time - overall_start_time) << " seconds." << std::endl;
	if ( num_structures_processed == 0 )
		basic::Warning() << "No structures processed.  Existing output files may have been skipped, did you mean to delete them?" << std::endl;
	jobdist.shutdown(); // under BOINC, this will cause program exit!
	return 0;
}

#include <devel/init.hh>
// AUTO-REMOVED #include <basic/options/option_macros.hh>
#include <basic/options/keys/packing.OptionKeys.gen.hh>
#include <basic/options/keys/docking.OptionKeys.gen.hh>

#include <core/chemical/ResidueType.hh>
#include <utility/vector0.hh>

/// This wrapper exists so David Kim's BOINC executable can call my real main() method.
int
main( int argc, char * argv [] )
{

	try {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	option.add_relevant(in::path::database);
	option.add_relevant(in::file::extra_res_fa);
	option.add_relevant(packing::unboundrot);
	option.add_relevant(packing::ex1::ex1);
	option.add_relevant(packing::ex1aro::ex1aro);
	option.add_relevant(packing::ex2::ex2);
	option.add_relevant(packing::extrachi_cutoff);
	option.add_relevant(packing::no_optH);
	option.add_relevant(packing::flip_HNQ);
	option.add_relevant(docking::ligand::soft_rep);
	option.add_relevant(docking::ligand::old_estat);

	option.add_relevant(in::file::s);
	option.add_relevant(in::file::native);
	option.add_relevant(out::nstruct);
	option.add_relevant(out::suffix);
	option.add_relevant(out::path::pdb);
	option.add_relevant(docking::uniform_trans);
	option.add_relevant(docking::randomize2);
	option.add_relevant(docking::dock_pert);
	option.add_relevant(docking::ligand::random_conformer);
	option.add_relevant(docking::ligand::improve_orientation);
	option.add_relevant(docking::ligand::minimize_ligand);
	option.add_relevant(docking::ligand::harmonic_torsions);
	option.add_relevant(docking::ligand::minimize_backbone);
	option.add_relevant(docking::ligand::harmonic_Calphas);
	option.add_relevant(docking::ligand::protocol);
	option.add_relevant(docking::ligand::start_from);
	option.add_relevant(docking::ligand::mutate_same_name3);
	option.add_relevant(enzdes::cstfile);

	// Parses command line options and inits RNG.
	// Doesn't seem to hurt to do it again if already done once (?)
	// Except in unit testing mode, where it wipes out e.g. -database
	devel::init(argc, argv);

	return ligand_dock_main_jd1();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

