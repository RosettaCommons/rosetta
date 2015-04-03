// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

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
#include <protocols/match/MatcherTask.hh>

#include <protocols/match/output/ProcessorFactory.hh>
#include <protocols/match/output/MatchProcessor.hh>


#include <core/pack/dunbrack/SingleLigandRotamerLibrary.hh>

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

void
set_ligpose_rotamer( core::pose::Pose & ligpose );

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
	using namespace core::io::pdb;
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
	if( input_jobs.size() == 0 ) utility_exit_with_message("No input scaffold structures specified for matcher. Check for -s <pdbfile> in arguments.");


	pose::Pose scaffold;
	core::import_pose::pose_from_pdb( scaffold, input_jobs[ 1 ] );
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

	time_t matcher_start_time = time(NULL), find_hits_end_time( time( NULL ) );
	long processing_time(0);
	MatcherOP matcher( new Matcher );
	matcher->initialize_from_task( *mtask );
	MatchProcessorOP processor = ProcessorFactory::create_processor( matcher, mtask );

	if( matcher->find_hits() ){
		find_hits_end_time = time(NULL);
		time_t process_start_time( time(NULL) );
		matcher->process_matches( *processor );
		processing_time = (long) (time(NULL) - process_start_time);

	}
	else{
		find_hits_end_time = time(NULL);
	}
	long find_hits_time = (long)(find_hits_end_time - matcher_start_time );
	time_t matcher_end_time = time(NULL);
	T << "Matcher ran for " << (long)(matcher_end_time - matcher_start_time) << " seconds, where finding hits took " << find_hits_time << " seconds and processing the matches took " << processing_time << " seconds."  << std::endl;

}

void
set_ligpose_rotamer( core::pose::Pose & ligpose )
{

	// Retrieve the rotamer library for this ligand;
	// check that the requested ligand-rotamer-index is in-bounds.
	// Relplace-residue on the ligpose with the rotamer from the library.

	using namespace core;
	using namespace core::conformation;
	using namespace core::scoring;
	using namespace core::pack::dunbrack;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	runtime_assert( ligpose.total_residue() == 1 ); // we're expecting a one-residue pose.

	if ( option[ match::ligand_rotamer_index ] < 1 ) {
		utility_exit_with_message( "Illegal rotamer index given in command line flag match::ligand_rotamer_index ("
			+ utility::to_string( option[ match::ligand_rotamer_index ]() ) + ").  Must be greater than 0." );
	}

	RotamerLibrary const & rotlib( * core::pack::dunbrack::RotamerLibrary::get_instance() );
	SingleResidueRotamerLibraryCOP res_rotlib( rotlib.get_rsd_library( ligpose.residue_type( 1 ) ) );

	if ( res_rotlib != 0 ) {
		SingleLigandRotamerLibraryCOP lig_rotlib(
			utility::pointer::dynamic_pointer_cast< SingleLigandRotamerLibrary const >
			( res_rotlib ));

		if ( lig_rotlib == 0 ) {
			utility_exit_with_message( "Failed to retrieve a ligand rotamer library for "
				+ ligpose.residue_type(1).name() + " after finding the flag match::ligand_rotamer_index <int> on the command line");
		}

		Size const nligrots = lig_rotlib->get_rotamers().size();

		if ( (Size) option[ match::ligand_rotamer_index ] > nligrots ) {
			utility_exit_with_message( "Illegal rotamer index given in command line flag match::ligand_rotamer_index ("
				+ utility::to_string( option[ match::ligand_rotamer_index ]() ) + "). Index exceeds the number"
				" of ligand rotamers ( " + utility::to_string( nligrots ) + ")" );
		}

		ResidueCOP ligrot = lig_rotlib->get_rotamers()[ option[ match::ligand_rotamer_index ] ];
		ligpose.replace_residue( 1, *ligrot, false );
	} else {
		utility_exit_with_message( "Failed to find ligand rotamer library for " +
			ligpose.residue(1).name() + " after finding the flag -match::ligand_rotamer_index on the command line." );
	}

}

