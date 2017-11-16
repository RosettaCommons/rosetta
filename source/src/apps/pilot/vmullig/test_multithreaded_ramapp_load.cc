// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/vmullig/test_multithreaded_ramapp_load.cc
/// @brief A pilot app to test the load of RamaPrePro data when many threads are trying to do this at once.
/// @author Vikram K. Mulligan (vmullig@u.washington.edu)

// devel headers
#include <devel/init.hh>

// core headers
#include <core/types.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/RamaPrePro.hh>
#include <core/scoring/Ramachandran2B.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/PairEPotential.hh>
#include <core/scoring/GenBornPotential.hh>
#include <core/scoring/HydroxylTorsionPotential.hh>
#include <core/scoring/facts/FACTSPotential.hh>
#include <core/scoring/dna/DNA_BasePotential.hh>
#include <core/pose/Pose.hh>
#include <core/import_pose/import_pose.hh>
#include <core/conformation/Conformation.hh>
#include <core/chemical/ResidueType.hh>

// protocol headers

// utility headers
#include <utility/excn/Exceptions.hh>

// basic headers
#include <basic/Tracer.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <utility/options/OptionCollection.hh>
#include <basic/options/option_macros.hh>

#include <boost/bind.hpp>

#ifdef MULTI_THREADED
#include <thread>
#include <mutex>
#include <atomic>
#endif

static basic::Tracer TR("test_multithreaded_ramapp_load");

OPT_KEY( String, testname )

void register_options() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	option.add_relevant( in::file::s );
	NEW_OPT ( testname, "The name of the test to run.  Options include \"pairE\", \"genborn\", \"hxl\", \"facts\", \"dnabase\", \"rama_prepro\", \"rama2b\", and \"rama\".  Defaults to \"rama_prepro\".", "rama_prepro" );
}

void
rama_prepro_test( core::Size const thread_index, core::pose::PoseOP pose ) {
	core::scoring::RamaPrePro const & ramapp( core::scoring::ScoringManager::get_instance()->get_RamaPrePro() );
	core::Real rama_score;
	utility::vector1< core::Real > grad;
	utility::vector1< core::Real > mainchain_tors(2);
	mainchain_tors[1] = pose->phi(5);
	mainchain_tors[2] = pose->psi(5);
	ramapp.eval_rpp_rama_score(pose->conformation(), pose->residue_type_ptr(5), pose->residue_type_ptr(6), mainchain_tors, rama_score, grad, false);
	TR << "Thread " << thread_index << " evaluated ramaprepro score of " << rama_score << "." << std::endl;
}

void
rama2b_test( core::Size const thread_index, core::pose::PoseOP pose ) {
	core::scoring::Ramachandran2B const &rama2b( core::scoring::ScoringManager::get_instance()->get_Ramachandran2B() );
	core::Real rama_score, drama_dphi, drama_dpsi;
	rama2b.eval_rama_score_residue( pose->residue_type(5).aa(), pose->phi(5), pose->psi(5), rama_score, drama_dphi, drama_dpsi );
	TR << "Thread " << thread_index << " evaluated rama2b score of " << rama_score << "." << std::endl;
}

void
rama_test( core::Size const thread_index, core::pose::PoseOP pose ) {
	core::scoring::Ramachandran const &rama( core::scoring::ScoringManager::get_instance()->get_Ramachandran() );
	core::Real const rama_score( rama.eval_rama_score_residue( pose->residue_type(5).aa(), pose->phi(5), pose->psi(5) ) );
	TR << "Thread " << thread_index << " evaluated rama score of " << rama_score << "." << std::endl;
}

void
pairE_test( core::Size const thread_index, core::pose::PoseOP pose ) {
	core::scoring::PairEPotential const &pairE( core::scoring::ScoringManager::get_instance()->get_PairEPotential() );
	core::Real const pairE_score( pairE.pair_term_energy( pose->residue(5), 10, pose->residue(6), 10 ) );
	TR << "Thread " << thread_index << " evaluated pairE score of " << pairE_score << "." << std::endl;
}

void
genborn_test( core::Size const thread_index, core::pose::PoseOP pose ) {
	core::scoring::GenBornPotential const &genborn( core::scoring::ScoringManager::get_instance()->get_GenBornPotential() );
	genborn.get_all_born_radii( *pose );
	TR << "Thread " << thread_index << " evaluated Born radii for pose." << std::endl;
}

void
hxl_test( core::Size const thread_index, core::pose::PoseOP pose ) {
	core::scoring::HydroxylTorsionPotential const &hxl( core::scoring::ScoringManager::get_instance()->get_HydroxylTorsionPotential() );
	core::Real const hxlE( hxl.eval_residue_energy( pose->residue(5) ) );
	TR << "Thread " << thread_index << " evaluated hydrxyl torsion energy of " << hxlE << "." << std::endl;
}

void
facts_test( core::Size const thread_index, core::pose::PoseOP pose ) {
	core::scoring::FACTSPotential const &facts( core::scoring::ScoringManager::get_instance()->get_FACTSPotential() );
	facts.setup_for_derivatives( *pose );
	TR << "Thread " << thread_index << " called FACTSPotential::setup_for_derivatives()." << std::endl;
}

void
dnabase_test( core::Size const thread_index, core::pose::PoseOP pose ) {
	core::scoring::dna::DNA_BasePotential const &dnabase( core::scoring::ScoringManager::get_instance()->get_DNA_BasePotential() );
	core::Real const dnaE( dnabase.base_pair_score( pose->residue(5), pose->residue(6) ) );
	TR << "Thread " << thread_index << " evaluated DNA base pair score of " << dnaE << "." << std::endl;
}

void
thread_fxn( core::Size const thread_index, core::pose::PoseOP pose ) {
	TR << "Launching thread " << thread_index << std::endl;

	std::string const test_name( basic::options::option[basic::options::OptionKeys::testname]() );
	if ( !test_name.compare( "rama_prepro" ) ) {
		rama_prepro_test(thread_index, pose);
	} else if ( !test_name.compare("rama2b") ) {
		rama2b_test( thread_index, pose );
	} else if ( !test_name.compare("rama") ) {
		rama_test( thread_index, pose );
	} else if ( !test_name.compare("pairE") ) {
		pairE_test( thread_index, pose );
	} else if ( !test_name.compare("genborn") ) {
		genborn_test( thread_index, pose );
	} else if ( !test_name.compare("hxl") ) {
		hxl_test( thread_index, pose );
	} else if ( !test_name.compare("facts") ) {
		facts_test( thread_index, pose );
	} else if ( !test_name.compare("dnabase") ) {
		dnabase_test( thread_index, pose );
	} else {
		TR.Warning << "Invalid test specified!  (Could not recognize \"" << test_name << "\"." << std::endl;
	}

	TR << "Thread " << thread_index << " terminated." << std::endl;
}


int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;
		using namespace basic::options::OptionKeys;

		register_options();
		devel::init( argc, argv );

		if ( ! option [ in::file::s ].user() ) {
			utility_exit_with_message("Please specify -in:file:s for the input PDB.");
		}

		utility::vector1< core::pose::PoseOP > poses;
		poses.push_back( core::import_pose::pose_from_file( basic::options::option[basic::options::OptionKeys::in::file::s]()[1], false, core::import_pose::PDB_file ) );
		for ( core::Size i=2; i<=6; ++i ) {
			poses.push_back(poses[1]->clone());
		}

		utility::vector1< std::thread > threads;
		for ( core::Size i(1); i<=6; ++i ) {
			threads.push_back( std::thread( boost::bind( &thread_fxn, i, poses[i] ) ) );
		}
		for ( core::Size i(1); i<=6; ++i ) {
			threads[i].join();
		}

		TR << "Execution completed." << std::endl;


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "Caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
}
