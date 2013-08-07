// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

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
#include <core/chemical/rna/RNA_Util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/id/types.hh>
#include <core/init/init.hh>
#include <core/optimization/CartesianMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/pose/util.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>
#include <core/import_pose/import_pose.hh>

//////////////////////////////////////////////////////////
#include <protocols/swa/rna/ERRASER_Modeler.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/rna/RNA_ProtocolUtil.hh>
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

static numeric::random::RandomGenerator RG(2391021);  // <- Magic number, do not change it!

OPT_KEY( IntegerVector, sample_res )
OPT_KEY( Integer, cycles )
OPT_KEY( Real, temperature )
OPT_KEY( Boolean, cart_min )

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// eventually could make this a class.
// writing as an app for quick testing.
// places to cleanup are marked below
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
prepare_fold_tree_for_erraser( pose::Pose & pose,
															 core::Size const erraser_res,
															 utility::vector1< Size > const sample_res_list ){

  using namespace core::kinematics;
  using namespace core::chemical::rna;

	FoldTree f = pose.fold_tree();

	Size start( sample_res_list[1]-1 ), end( sample_res_list[ sample_res_list.size() ]+1 );
	Size const cutpoint = erraser_res;

	std::cout << "ADDING CUTPOINT TO FOLD_TREE: " << start << " " << end << " " << cutpoint << std::endl;
	f.new_jump( start, end, cutpoint );

	Size const which_jump = f.jump_nr( start, end );

	runtime_assert( which_jump > 0 );
	runtime_assert( f.upstream_jump_residue( which_jump ) == start );
	runtime_assert( f.downstream_jump_residue( which_jump ) == end );

	// may need to be smarter about whether start/end are upstream/downstream
	f.set_jump_atoms( which_jump, default_jump_atom( pose.residue(start) ), default_jump_atom( pose.residue(end) ) );

	pose.fold_tree( f );
	std::cout << "NEW FOLD TREE " << pose.fold_tree() << std::endl;
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
	using namespace protocols::swa::rna;
	using namespace protocols::moves;

	clock_t const time_start( clock() );

	// read starting pose(s) from disk
	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	pose::Pose pose;
	import_pose::pose_from_pdb( pose, *rsd_set, option[ in::file::s ][1] );
	protocols::rna::figure_out_reasonable_rna_fold_tree( pose );
	protocols::rna::virtualize_5prime_phosphates( pose );

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
	core::scoring::ScoreFunctionOP scorefxn = getScoreFunction();

	bool const cart_min_ = option[ cart_min ]();
	if ( cart_min_ && scorefxn->get_weight( cart_bonded ) == 0.0 ) scorefxn->set_weight( cart_bonded, 1.0 );

	// ned to package following into a function or Mover.
	CartesianMinimizer cartesian_minimizer;
	float const dummy_tol( 0.00000025);
	bool const use_nblist( true );
	MinimizerOptions options( "dfpmin", dummy_tol, use_nblist, false, false );
	options.nblist_auto_update( true );
	MoveMap move_map;
	figure_out_swa_rna_movemap( move_map, pose, sample_res_list );

	// encapsulate this:
	if ( cart_min_ ){
		std::cout << "BEFORE CARTESIAN MINIMIZE" << std::endl;
		scorefxn->show( std::cout, pose );
		cartesian_minimizer.run( pose, move_map, *scorefxn, options );
		std::cout << "AFTER CARTESIAN MINIMIZE" << std::endl;
		pose.dump_pdb( "after_cart_min.pdb" );
		scorefxn->show( std::cout, pose );
	}

	utility::vector1< Size > fixed_res;
	for ( Size i = 1; i <= pose.total_residue(); i++ ) if ( !sample_res_list.has_value( i) ) fixed_res.push_back( i );

	// monte carlo setup.
	MonteCarloOP monte_carlo_ = new MonteCarlo( pose, *scorefxn, option[ temperature]() );

	for ( Size n = 1; n <= option[ cycles ](); n++ ){

		FoldTree const f_save = pose.fold_tree();

		// create reasonable fold tree
		Size const erraser_res = RG.random_element( sample_res_list );
		prepare_fold_tree_for_erraser( pose, erraser_res, sample_res_list );

		// add chainbreak variants -- put this in fold_tree?
		add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_LOWER, erraser_res   );
		add_variant_type_to_pose_residue( pose, chemical::CUTPOINT_UPPER, erraser_res+1 );

		// run ERRASER_Modeler
		ERRASER_ModelerOP erraser_modeler = new ERRASER_Modeler( erraser_res, scorefxn );
		erraser_modeler->set_choose_random( true );
		erraser_modeler->set_force_centroid_interaction( true );
		erraser_modeler->set_fixed_res( fixed_res );
		erraser_modeler->set_use_phenix_geo ( option[ basic::options::OptionKeys::rna::corrected_geo ]() );
		//		erraser_modeler->set_verbose( true );
		if ( n == 1 ) erraser_modeler->set_skip_sampling( true ); // just do minimize on first round.

		erraser_modeler->apply( pose );

		std::string const move_type = erraser_modeler->get_num_sampled() ? "erraser" : "erraser_no_op";

		// display not necessary?
		std::cout << "After ERRASER move" << std::endl;
		scorefxn->show( std::cout, pose );

		// remove chainbreak variants. along with fold_tree restorer, put into separate function.
		remove_variant_type_from_pose_residue( pose, chemical::CUTPOINT_LOWER, erraser_res   );
		remove_variant_type_from_pose_residue( pose, chemical::CUTPOINT_UPPER, erraser_res+1 );
		// return to simple fold tree
		pose.fold_tree( f_save );

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
	NEW_OPT( cycles, "Number of Monte Carlo cycles", 50 );
	NEW_OPT( temperature, "Monte Carlo temperature", 1.0 );
	NEW_OPT( cart_min, "Do Cartesian minimizations", false );

  ////////////////////////////////////////////////////////////////////////////
  // setup
  ////////////////////////////////////////////////////////////////////////////
  core::init::init(argc, argv);


  ////////////////////////////////////////////////////////////////////////////
  // end of setup
  ////////////////////////////////////////////////////////////////////////////
  protocols::viewer::viewer_main( my_main );


	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
	}

}



