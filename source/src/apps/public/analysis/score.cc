// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file Mike Tyka
/// @brief


// libRosetta headers
#include <protocols/simple_moves/ScoreMover.hh>
#include <protocols/simple_moves/TailsScoreMover.hh>

#include <protocols/jobdist/standard_mains.hh>

#include <protocols/simple_moves/MinMover.hh>
#include <protocols/simple_moves/ProlineFixMover.hh>
#include <protocols/moves/MoverContainer.hh>
#include <protocols/simple_moves/ConstraintSetMover.hh>
#include <protocols/topology_broker/TopologyBroker.hh>
#include <protocols/topology_broker/util.hh>

#include <core/scoring/ScoreFunctionFactory.hh> // get_score_function
#include <core/scoring/constraints/util.hh>
#include <core/scoring/electron_density/util.hh>
#include <protocols/simple_moves/symmetry/SetupForSymmetryMover.hh>
#include <protocols/simple_moves/SuperimposeMover.hh>
#include <core/kinematics/MoveMap.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/exit.hh>

#include <basic/Tracer.hh>
#include <devel/init.hh>

// C++ headers
#include <iostream>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/rescore.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/symmetry.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/relax.OptionKeys.gen.hh>
#include <basic/options/keys/krassk.OptionKeys.gen.hh>
#include <basic/options/keys/broker.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <basic/options/option.hh>

using namespace core;
using namespace basic::options;

static THREAD_LOCAL basic::Tracer TR( "protocols.moves.ScoreMover" );

namespace score_app { BooleanOptionKey linmin( "score_app:linmin" );
BooleanOptionKey superimpose_to_native( "score_app:superimpose_to_native" ); }

int
main( int argc, char * argv [] )
{
	try {

		using namespace protocols;
		using namespace protocols::moves;
		using namespace core::scoring;
		using namespace basic::options;
		using namespace basic::options::OptionKeys;
		using namespace utility::file;

		simple_moves::ScoreMover::register_options();
		protocols::jobdist::register_options_universal_main();
		option.add_relevant( in::file::fullatom );
		option.add_relevant( relax::fast );
		option.add( score_app::linmin, "Do a quick round of linmin before reporting the score" );
		option.add_relevant( score_app::linmin );
		option.add_relevant( out::output                  );
		option.add_relevant( out::nooutput                );
		option.add_relevant( in::file::fullatom           );
		option.add_relevant( rescore::verbose             );
		option.add_relevant( in::file::repair_sidechains  );
		option.add( score_app::superimpose_to_native, "superimpose structure to native" );


		// scoring should by default not produce output files - that's so annoying
		// unless of coures the user insists.

		// initialize core
		devel::init(argc, argv);

		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << " Rosetta Tool:   score   -  rescores PDBs and silent files, extracts PDBs from silent files, assembles PDBs into silent files. " << std::endl;
		std::cout << " Usage:                                                                  " << std::endl;
		std::cout << "   PDB input:      -in:file:s *.pdb   or  " << std::endl;
		std::cout << "                   -in:file:l  list_of_pdbs  " << std::endl;
		std::cout << "                   -no_optH                                    Dont change positions of Hydrogen atoms! (default true, specify false if you want optH)" << std::endl;
		std::cout << "   Silent input:   -in:file:silent silent.out                  silent input filesname " << std::endl;
		std::cout << "                   -in:file:tags                               specify specific tags to be extracted, if left out all will be taken " << std::endl;
		std::cout << "                   -in:file:fullatom                           for full atom structures " << std::endl;
		std::cout << "                   -in:file:binary_silentfile                  for non-ideal structures (such as from looprelax) " << std::endl;
		std::cout << "                   -in:file:silent_optH                        Call optH when reading silent files (useful for HisD/HisE determination)" << std::endl;
		std::cout << "                   -score_app:linmin                           Run a quick linmin before scoring" << std::endl;
		std::cout << "   Native:         -in:file:native                             native PDB (rms, maxsub and gdtm scores will be calculated)" << std::endl;
		std::cout << "   Scorefunction:  -score:weights  weights                     weight set or weights file " << std::endl;
		std::cout << "                   -score:patch  patch                         patch set " << std::endl;
		std::cout << "                   -score:optH_weights                         Weights file for optH (default standard.wts w/ sc12 patch)" << std::endl;
		std::cout << "                   -score:optH_patch                           Weights patch file for optH" << std::endl;
		std::cout << "                   -rescore:verbose                            display score breakdown " << std::endl;
		std::cout << "   Output:         -out:nooutput                               don't print PDB structures (default now) " << std::endl;
		std::cout << "                   -out:output                                 force printing of PDB structures " << std::endl;
		std::cout << "                   -out:file:silent                            write silent-out file " << std::endl;
		std::cout << "                   -out:file:scorefile name                    write scorefile (default default.sc)" << std::endl;
		std::cout << "                   -out:prefix  myprefix                       prefix the output structures with a string " << std::endl;
		std::cout << "  Examples: " << std::endl;
		std::cout << "   score  -database ~/minirosetta_database -in:file:silent silent.out -in::file::binary_silentfile -in::file::fullatom -native 1a19.pdb " << std::endl;
		std::cout << "   Will rescore all structures in silent.out, in full atom mode and accounting for nonideal structure if present. Additionally " << std::endl;
		std::cout << "   it will print a PDB for every structure with -out:output flag " << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;

		//The following lines are to ensure one can rescore the pcs energy term (that uses TopologyClaimer)
		if ( option[ broker::setup ].user() ) {
			protocols::topology_broker::TopologyBrokerOP top_bro_OP( new  topology_broker::TopologyBroker() );
			try {
				add_cmdline_claims(*top_bro_OP, false /* do_I_need_fragments */);
			}
catch ( utility::excn::EXCN_Exception &excn )  {
	excn.show( TR.Error );
	utility_exit();
}
		}

		// do not output pdb by default, unless with -out:output flag
		if ( !option[ out::output ].user() ) {
			option[ out::nooutput ].value( true );
		}

		// get scorefxn and add constraints if defined
		core::scoring::ScoreFunctionOP sfxn = core::scoring::get_score_function();
		if ( option[ in::file::fullatom ]() ) {
			core::scoring::constraints::add_fa_constraints_from_cmdline_to_scorefxn( *sfxn );
		} else {
			core::scoring::constraints::add_constraints_from_cmdline_to_scorefxn( *sfxn );
		}

		// now add density scores from cmd line
		if ( option[ edensity::mapfile ].user() ) {
			core::scoring::electron_density::add_dens_scores_from_cmdline_to_scorefxn( *sfxn );
		}

		// create a ScoreMover

		simple_moves::ScoreMoverOP scoretmp;
		if ( option[ krassk::tail_mode] ) {
			scoretmp = simple_moves::ScoreMoverOP( new simple_moves::TailsScoreMover(sfxn) );
		} else {
			scoretmp = simple_moves::ScoreMoverOP( new simple_moves::ScoreMover(sfxn) );
		}

		if (  option[ rescore::verbose ] ) {
			scoretmp->set_verbose( true );
		} else {
			scoretmp->set_verbose( false );
		}

		// save it to a mover that will be passed to job_distributor
		MoverOP mover = scoretmp;

		// do sth more than just scoring
		if ( option[ score_app::linmin ]() || option[ in::file::repair_sidechains ]() ) {
			assert( sfxn );
			protocols::moves::SequenceMoverOP seqmov( new protocols::moves::SequenceMover );
			if ( option[ in::file::repair_sidechains ]() ) {
				protocols::simple_moves::ProlineFixMoverOP pfm( new protocols::simple_moves::ProlineFixMover );
				seqmov->add_mover( pfm );
			}
			if ( option[ score_app::linmin ]() ) {
				core::kinematics::MoveMapOP movemap( new core::kinematics::MoveMap );
				movemap->set_bb( true ); movemap->set_chi( true );
				protocols::simple_moves::MinMoverOP minmover( new protocols::simple_moves::MinMover(
					movemap, sfxn, "linmin", 1e-4,
					true /*use_nblist*/, true /*deriv_check*/, true /*verbose driv check*/ ) );
				seqmov->add_mover( minmover );
			}

			seqmov->add_mover( mover );
			mover = seqmov;
		} // if infile remediation necessary

		// add constraints from cmd line
		if ( option[ OptionKeys::constraints::cst_fa_file ].user() || option[ OptionKeys::constraints::cst_file ].user() ) {
			protocols::moves::SequenceMoverOP seqmov( new protocols::moves::SequenceMover );
			protocols::simple_moves::ConstraintSetMoverOP loadCsts( new protocols::simple_moves::ConstraintSetMover );
			if ( option[ OptionKeys::constraints::cst_fa_file ].user() ) {
				loadCsts->constraint_file( core::scoring::constraints::get_cst_fa_file_option() );
			} else {
				loadCsts->constraint_file( core::scoring::constraints::get_cst_file_option() );
			}
			seqmov->add_mover( loadCsts );
			seqmov->add_mover( mover );
			mover = seqmov;
		}

		// set pose for density scoring if a map was input
		//   + (potentially) dock map into density
		if ( option[ edensity::mapfile ].user() ) {
			protocols::moves::SequenceMoverOP seqmov( new protocols::moves::SequenceMover );
			seqmov->add_mover( MoverOP( new protocols::electron_density::SetupForDensityScoringMover ) );
			seqmov->add_mover( mover );
			mover = seqmov;
		}

		// set pose for symmetry
		if ( option[ OptionKeys::symmetry::symmetry_definition ].user() )  {
			protocols::moves::SequenceMoverOP seqmov( new protocols::moves::SequenceMover );
			seqmov->add_mover( MoverOP( new protocols::simple_moves::symmetry::SetupForSymmetryMover ) );
			seqmov->add_mover( mover );
			mover = seqmov;
		}

		if ( option[ score_app::superimpose_to_native ]() ) {
			if ( !option[ in::file::native ].user() ) {
				TR << "No native specified. Cannot align to native..." << '\n';
			} else {
				// read native structure
				core::pose::Pose native;
				core::import_pose::pose_from_pdb( native, option[ basic::options::OptionKeys::in::file::native ] );
				protocols::moves::SequenceMoverOP seqmov( new protocols::moves::SequenceMover );
				seqmov->add_mover( MoverOP( new protocols::simple_moves::SuperimposeMover( native ) ) );
				seqmov->add_mover( mover );
				mover = seqmov;
			}
		}

		// operate this mover and output pdbs/scorefile
		protocols::jobdist::universal_main( *mover );

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

