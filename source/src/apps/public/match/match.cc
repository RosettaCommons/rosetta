// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   apps/public/match.cc
/// @brief
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini


#include <devel/init.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>

#include <core/id/AtomID.hh>

#include <basic/options/option.hh>

#include <basic/options/keys/match.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>

#include <core/pose/Pose.hh>
//#include <core/scoring/ScoreFunctionFactory.hh>

#include <basic/options/util.hh>

#include <protocols/match/Matcher.hh>
#include <protocols/match/MatcherMover.hh>
#include <protocols/match/MatcherTask.hh>

#include <protocols/match/output/ProcessorFactory.hh>
#include <protocols/match/output/MatchProcessor.hh>


#include <core/pack/rotamers/SingleLigandRotamerLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibraryFactory.hh>

#include <basic/Tracer.hh>

#include <utility/string_util.hh>

#include <core/import_pose/import_pose.hh>
#include <utility/vector1.hh>

#include <utility/excn/Exceptions.hh>

#if defined(WIN32) || defined(__CYGWIN__)
#include <ctime>
#endif


void
match_main();

/* moved to protocols/match/MatcherMover.cc as one step toward integrating match app with MatcherMover */
//void
//set_ligpose_rotamer( core::pose::Pose & ligpose );

int main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		devel::init( argc, argv );

		match_main();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
}

void
match_main()
{
	using namespace core;
	using namespace core::chemical;
	using namespace core::id;
	using namespace core::io;
	using namespace core::pose;
	using namespace protocols::match;
	using namespace protocols::match::downstream;
	using namespace protocols::match::output;
	using namespace protocols::match::upstream;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	basic::Tracer T( "apps.public.match.match" );

	//we need this for the output to be correct
	option[OptionKeys::run::preserve_header ].value(true);

	MatcherTaskOP mtask( new MatcherTask );

	utility::vector1< std::string > input_jobs = basic::options::start_files();
	if ( input_jobs.size() == 0 ) utility_exit_with_message("No input scaffold structures specified for matcher. Check for -s <pdbfile> in arguments.");


	pose::Pose scaffold;
	core::import_pose::pose_from_file( scaffold, input_jobs[ 1 ] , core::import_pose::PDB_file);
	scaffold.update_residue_neighbors();

	pose::Pose ligpose;
	core::conformation::ResidueOP ligres = core::conformation::ResidueFactory::create_residue(
		core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD )->name_map(
		option[ basic::options::OptionKeys::match::lig_name ] )
	);
	ligpose.append_residue_by_jump( *ligres, 1 );

	if ( option[ OptionKeys::match::ligand_rotamer_index ].user() ) {
		set_ligpose_rotamer( ligpose );
	}


	Size cent, nbr1, nbr2;
	ligres->select_orient_atoms( cent, nbr1, nbr2 );

	T << "Selecting orientation atoms:";
	T << " " << ligres->atom_name( cent );
	T << " " << ligres->atom_name( nbr1 );
	T << " " << ligres->atom_name( nbr2 ) << std::endl;

	mtask->set_upstream_pose( scaffold );

	utility::vector1< core::id::AtomID > oats( 3 );
	oats[ 1 ] = AtomID( nbr2, 1 ); oats[ 2 ] = AtomID( nbr1, 1 ); oats[ 3 ] = AtomID( cent, 1 );

	mtask->set_downstream_pose( ligpose,  oats );
	mtask->initialize_from_command_line();
	//task->consolidate_matches( false );

	time_t matcher_start_time = time(NULL);
	long processing_time(0);
	MatcherOP matcher( new Matcher );
	matcher->initialize_from_task( *mtask );
	MatchProcessorOP processor = ProcessorFactory::create_processor( matcher, mtask );

	time_t find_hits_end_time;
	if ( matcher->find_hits() ) {
		find_hits_end_time = time(NULL);
		time_t process_start_time( time(NULL) );
		matcher->process_matches( *processor );
		processing_time = (long) (time(NULL) - process_start_time);

	} else {
		find_hits_end_time = time(NULL);
	}
	long find_hits_time = (long)(find_hits_end_time - matcher_start_time );
	time_t matcher_end_time = time(NULL);
	T << "Matcher ran for " << (long)(matcher_end_time - matcher_start_time) << " seconds, where finding hits took " << find_hits_time << " seconds and processing the matches took " << processing_time << " seconds."  << std::endl;

}

