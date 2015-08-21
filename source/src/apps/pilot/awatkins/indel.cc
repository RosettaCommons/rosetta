// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

// Project Headers
#include <utility/pointer/owning_ptr.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/jd2/Job.hh>

// Filter headers
#include <protocols/indel/IndelOptimizationMover.hh>

// Utility Headers
#include <devel/init.hh>
#include <basic/options/util.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/Tracer.hh>
#include <utility/exit.hh>

// C++ headers
#include <string>
#include <sstream>

//The original author used a lot of using declarations here.  This is a stylistic choice.
// Namespaces
using namespace core;
using namespace basic::options;
using namespace basic::options::OptionKeys;

// application specific options
namespace indel {
// pert options
IntegerOptionKey const start_res( "indel::start_res" );
IntegerOptionKey const end_res( "indel::end_res" );
IntegerOptionKey const loop_length( "indel::loop_length" );
IntegerOptionKey const num_to_dock( "indel::num_to_dock" );
BooleanOptionKey const dump_initial_results( "indel::dump_initial_results" );
}

int
main( int argc, char* argv[] )
{
	try {
		/*********************************************************************************************************************
		Common Setup
		***************************( *******************************************************************************************/

		// add application specific options to options system

		option.add( indel::start_res, "The first residue to delete" ).def( 1 );
		option.add( indel::end_res, "The last residue to delete (will be set to start_res if not specified)" ).def( 1 );
		option.add( indel::loop_length, "The number of residues behind and in front of the deleted residue(s) to remodel" ).def( 4 );
		option.add( indel::num_to_dock, "If docking is necessary, the number of models to build before docking" ).def( 10 );
		option.add( indel::dump_initial_results, "Dump the initial PDBs created pre-docking" ).def( false );

		devel::init(argc, argv);

		Size start_res = option[ indel::start_res ].value();
		Size end_res;
		if ( option[ indel::end_res ].user() ) {
			end_res = option[ indel::end_res ].value();
		} else {
			end_res = start_res;
		}

		Size loop_length = option[ indel::loop_length ].value();

		// Set up LoopRelaxMover
		std::string remodel            ( option[ OptionKeys::loops::remodel ]() );
		std::string const intermedrelax( option[ OptionKeys::loops::intermedrelax ]() );
		std::string const refine       ( option[ OptionKeys::loops::refine ]() );
		std::string const relax        ( option[ OptionKeys::loops::relax ]() );
		bool frag_files = option[ OptionKeys::loops::frag_files ].user();
		Size num_to_dock( option[ indel::num_to_dock ]() );
		bool dump_initial_results = option[ indel::dump_initial_results ].user();

		//create mover instance
		protocols::indel::IndelOptimizationMoverOP indel_mover(
			new protocols::indel::IndelOptimizationMover( start_res, end_res, loop_length,
			remodel, intermedrelax, refine, relax, frag_files, num_to_dock,
			dump_initial_results ) );

		//call job distributor
		protocols::jd2::JobDistributor::get_instance()->go( indel_mover );
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cerr << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}//main
