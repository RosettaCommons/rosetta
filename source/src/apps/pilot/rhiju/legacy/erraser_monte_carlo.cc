// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file swa_monte_Carlo.cc
/// @author Rhiju Das (rhiju@stanford.edu)

// libRosetta headers
#include <core/types.hh>

// do we need all these?
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/rms_util.hh>
#include <core/chemical/rna/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/id/types.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/import_pose/import_pose.hh>
#include <devel/init.hh>


//////////////////////////////////////////////////////////
#include <protocols/stepwise/modeler/rna/StepWiseRNA_Modeler.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/monte_carlo/rna/TransientCutpointHandler.hh>
#include <protocols/farna/util.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/MonteCarlo.hh>

//////////////////////////////////////////////////
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> // for option[ out::file::silent  ] and etc.
#include <basic/options/keys/in.OptionKeys.gen.hh> // for option[ in::file::tags ] and etc.
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option_macros.hh>
#include <core/pose/Pose.hh>
#include <utility/vector1.hh>


//////////////////////////////////////////////////////////
#include <numeric/random/random.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>


// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <ctime>
#include <list>

using namespace core;
using namespace protocols;
using namespace ObjexxFCL;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using utility::vector1;
using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;


OPT_KEY( IntegerVector, sample_res )
OPT_KEY( IntegerVector, design_res )
OPT_KEY( Integer, cycles )
OPT_KEY( Real, temperature )
OPT_KEY( Boolean, cart_min )
OPT_KEY( Boolean, minimize_single_res )
OPT_KEY( Integer, num_random_samples )

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// eventually could make this a class.
// writing as an app for quick testing.
// places to cleanup are marked below
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void
setup_design_res( pose::Pose & pose ){

	using namespace core::pose::full_model_info;
	utility::vector1< Size > design_res_ = option[ design_res ]();
	FullModelInfo & full_model_info = nonconst_full_model_info( pose );
	std::string full_sequence = full_model_info.full_sequence();

	for ( Size i = 1; i <= design_res_.size(); i++ ){
		full_sequence[ design_res_[ i ] - 1 ] = 'n';
	}

	full_model_info.set_full_sequence( full_sequence );
	//	std::cout << "NEW SEQUENCE: " << const_full_model_info( pose ).full_sequence();

}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
erraser_monte_carlo()
{

  using namespace core::pose;
  using namespace core::chemical;
  using namespace core::kinematics;
  using namespace core::scoring;
  using namespace core::io::silent;
  using namespace core::optimization;
  using namespace core::import_pose;
	using namespace protocols::stepwise;
	using namespace protocols::stepwise::modeler::rna;
	using namespace protocols::stepwise::monte_carlo::rna;
	using namespace protocols::moves;

	clock_t const time_start( clock() );

	// read starting pose(s) from disk
	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	pose::Pose pose;
	import_pose::pose_from_pdb( pose, *rsd_set, option[ in::file::s ][1] );
	protocols::farna::figure_out_reasonable_rna_fold_tree( pose );
	protocols::farna::virtualize_5prime_phosphates( pose );
	setup_design_res( pose );

	protocols::viewer::add_conformation_viewer ( pose.conformation(), "current", 400, 400 );

	PoseOP native_pose;
	if (option[ in::file::native ].user() ) {
		native_pose = PoseOP( new Pose );
		import_pose::pose_from_pdb( *native_pose, *rsd_set, option[ in::file::native ]() );
	}

	// what are residues that can move? Assume they are contiguous. Can later put in a check!
	if ( !option[ sample_res ].user() ) utility_exit_with_message( "Must specify sample_res" );
	utility::vector1< Size > sample_res_list = option[ sample_res ]();

	// scorefunction
	core::scoring::ScoreFunctionOP scorefxn = get_score_function();

	bool const cart_min_ = option[ cart_min ]();
	if ( cart_min_ && scorefxn->get_weight( cart_bonded ) == 0.0 ) scorefxn->set_weight( cart_bonded, 1.0 );

	// ned to package following into a function or Mover.
	CartesianMinimizer cartesian_minimizer;
	float const dummy_tol( 0.00000025);
	bool const use_nblist( true );
	MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false, false );
	options.nblist_auto_update( true );
	MoveMap move_map;
	figure_out_stepwise_rna_movemap( move_map, pose, sample_res_list );

	// encapsulate this:
	if ( cart_min_ ){
		std::cout << "BEFORE CARTESIAN MINIMIZE" << std::endl;
		scorefxn->show( std::cout, pose );
		cartesian_minimizer.run( pose, move_map, *scorefxn, options );
		std::cout << "AFTER CARTESIAN MINIMIZE" << std::endl;
		pose.dump_pdb( "after_cart_min.pdb" );
		scorefxn->show( std::cout, pose );
	}

	// monte carlo setup.
	MonteCarloOP monte_carlo_ = new MonteCarlo( pose, *scorefxn, option[ temperature]() );

	for ( Size n = 1; n <= Size( option[ cycles ]() ); n++ ){

		Size const erraser_res = numeric::random::rg().random_element( sample_res_list );
		bool const did_mutation = mutate_res_if_allowed( pose, erraser_res );


		TransientCutpointHandler cutpoint_handler( erraser_res );
		cutpoint_handler.set_minimize_res( sample_res_list ); // its needed to figure out fold tree.
		cutpoint_handler.put_in_cutpoints( pose );

		StepWiseRNA_Modeler erraser_modeler( erraser_res, scorefxn );
		StepWiseModelerOptionsOP options = new StepWiseModelerOptions;
		options->set_choose_random( true );
		options->set_force_centroid_interaction( true );
		options->set_kic_modeler_if_relevant( true );
		options->set_use_phenix_geo( option[ basic::options::OptionKeys::rna::corrected_geo ]() );
		options->set_num_random_samples( option[ num_random_samples ]() );
		options->set_num_pose_minimize( 1 );
		erraser_modeler.set_options( options );
		if ( !option[ minimize_single_res ]() )  erraser_modeler.set_minimize_res( sample_res_list );

		erraser_modeler.apply( pose );
		std::string move_type = erraser_modeler.get_num_sampled() ? "erraser" : "erraser_no_op";
		if ( did_mutation ) move_type += "-mut";

		cutpoint_handler.take_out_cutpoints( pose );

		std::cout << "After ERRASER move and formally remove cutpoint" << std::endl;
		scorefxn->show( std::cout, pose );

		if ( cart_min_ ){
			cartesian_minimizer.run( pose, move_map, *scorefxn, options );
			std::cout << "After Cartesian Minimizer" << std::endl;
			scorefxn->show( std::cout, pose );
		}

		monte_carlo_->boltzmann( pose, move_type );
		std::cout << "Monte Carlo accepted? " << monte_carlo_->mc_accepted() << std::endl;
		monte_carlo_->show_counters();
		pose.dump_pdb( "latest.pdb" );
	}

	monte_carlo_->recover_low( pose );
	monte_carlo_->show_counters();

	pose.dump_pdb( "final.pdb" );
	std::cout << "Total time for monte carlo: " << static_cast<Real>( clock() - time_start ) / CLOCKS_PER_SEC << " seconds." << std::endl;

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	clock_t const my_main_time_start( clock() );

	erraser_monte_carlo();

	protocols::viewer::clear_conformation_viewers();

	std::cout << "Total time to run " << static_cast<Real>( clock() - my_main_time_start ) / CLOCKS_PER_SEC << " seconds." << std::endl;

  exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

  using namespace basic::options;
	utility::vector1< Size > blank_size_vector;
	utility::vector1< std::string > blank_string_vector;

	NEW_OPT( sample_res, "residues to build", blank_size_vector );
	NEW_OPT( design_res, "residues to re-design", blank_size_vector );
	NEW_OPT( cycles, "Number of Monte Carlo cycles", 50 );
	NEW_OPT( temperature, "Monte Carlo temperature", 1.0 );
	NEW_OPT( cart_min, "Do Cartesian minimizations", false );
	NEW_OPT( minimize_single_res, "Minimize the residue that just got rebuilt, instead of all", false );
	NEW_OPT( num_random_samples, "Number of samples from swa residue sampler before minimizing best", 1 );

  ////////////////////////////////////////////////////////////////////////////
  // setup
  ////////////////////////////////////////////////////////////////////////////
  devel::init(argc, argv);


  ////////////////////////////////////////////////////////////////////////////
  // end of setup
  ////////////////////////////////////////////////////////////////////////////
  protocols::viewer::viewer_main( my_main );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}



