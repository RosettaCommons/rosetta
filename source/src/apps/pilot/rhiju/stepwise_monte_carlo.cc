// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file swa_monte_carlo.cc
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
#include <core/chemical/rna/RNA_Util.hh>
#include <core/scoring/rms_util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/id/types.hh>
#include <core/init/init.hh>
#include <core/pose/util.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>
#include <core/import_pose/import_pose.hh>

//////////////////////////////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_Modeler.hh>
#include <protocols/swa/rna/StepWiseRNA_Modeler.hh>
#include <protocols/swa/StepWiseUtil.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/monte_carlo/RNA_SWA_MonteCarloUtil.hh>
#include <protocols/swa/monte_carlo/RNA_AddMover.hh>
#include <protocols/swa/monte_carlo/RNA_DeleteMover.hh>
#include <protocols/swa/monte_carlo/RNA_AddOrDeleteMover.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelInfoUtil.hh>
#include <protocols/swa/monte_carlo/types.hh>
#include <protocols/rna/RNA_ProtocolUtil.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <basic/datacache/BasicDataCache.hh>

//////////////////////////////////////////////////
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh> // for option[ out::file::silent  ] and etc.
#include <basic/options/keys/in.OptionKeys.gen.hh> // for option[ in::file::tags ] and etc.
#include <basic/options/keys/rna.OptionKeys.gen.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/option_macros.hh>
#include <basic/Tracer.hh>
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

static basic::Tracer TR( "apps.pilot.rhiju.stepwise_monte_carlo" );

OPT_KEY( IntegerVector, input_res )
OPT_KEY( IntegerVector, cutpoint_open )
OPT_KEY( Integer, cycles )
OPT_KEY( Real, minimize_single_res_frequency )
OPT_KEY( Boolean, allow_internal_moves )
OPT_KEY( Real, temperature )
OPT_KEY( Real, add_delete_frequency )
OPT_KEY( Real, just_min_after_mutation_frequency )
OPT_KEY( Integer, num_random_samples )
OPT_KEY( IntegerVector, sample_res )


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This should go to its own class soon!
// do it right after fang's erraser/swa unification & corrected_geo.
bool
apply_swa_mover( pose::Pose & pose,
								 core::scoring::ScoreFunctionOP scorefxn,
								 bool const minimize_single_res,
								 std::string & move_type ){

	using namespace protocols::swa;
	using namespace protocols::swa::rna;
	using namespace protocols::swa::monte_carlo;
	using namespace core::pose::full_model_info;

	utility::vector1< Size > const & moving_res = nonconst_full_model_info_from_pose( pose ).moving_res_list();
	if ( moving_res.size() == 0 ) return false;

	Size remodel_res( 0 );

	bool const allow_internal_moves_ = option[ allow_internal_moves ]();
	if ( allow_internal_moves_ ){
		remodel_res =	RG.random_element( moving_res );
	} else {
		utility::vector1< Size > possible_res;
		get_potential_resample_residues( pose, possible_res );
		if ( possible_res.size() == 0 ) return false;
		remodel_res =	RG.random_element( possible_res );
	}
	TR << "About to remodel residue " << remodel_res << " in " << pose.annotated_sequence() << std::endl;

	bool const did_mutation = mutate_res_if_allowed( pose, remodel_res ); // based on 'n' in full_model_info.full_sequence

	Real const just_min_after_mutation_frequency_ = option[ just_min_after_mutation_frequency ]();
	bool just_min_after_mutation_ = ( did_mutation && ( RG.uniform() < just_min_after_mutation_frequency_ ) );

	if ( is_at_terminus( pose, remodel_res ) ){
		swa::rna::StepWiseRNA_Modeler stepwise_rna_modeler( remodel_res, scorefxn );
		//	stepwise_rna_modeler->set_use_phenix_geo ( option[ basic::options::OptionKeys::rna::corrected_geo ]() );
		stepwise_rna_modeler.set_force_centroid_interaction( true );
		stepwise_rna_modeler.set_choose_random( true );
		stepwise_rna_modeler.set_num_random_samples( option[ num_random_samples ]() );
		stepwise_rna_modeler.set_num_pose_minimize( 1 );

		if ( just_min_after_mutation_ ) stepwise_rna_modeler.set_skip_sampling( true );
		if ( ! minimize_single_res ) stepwise_rna_modeler.set_minimize_res( moving_res );

		stepwise_rna_modeler.apply( pose );

		move_type = "swa";

	} else {

		runtime_assert( allow_internal_moves_ );

		// note that following looks almost exactly like stepwise_rna_modeler -- would be best to unify the classes!
		// we should be able to do this after Fang unifies StepWiseRNA_ResidueSampler and StepWiseRNA_AnalyticCloseSampler

		std::cout << "GOING TO RUN ERRASER " << remodel_res << " in " << pose.annotated_sequence() << std::endl;
		StepWiseRNA_Modeler erraser_modeler( remodel_res, scorefxn );
		erraser_modeler.set_use_phenix_geo ( option[ basic::options::OptionKeys::rna::corrected_geo ]() );
		erraser_modeler.set_force_centroid_interaction( true );
		erraser_modeler.set_kic_sampling( true );
		erraser_modeler.set_choose_random( true );
		erraser_modeler.set_num_random_samples( option[ num_random_samples ]() );
		erraser_modeler.set_num_pose_minimize( 1 );

		if ( just_min_after_mutation_ ) erraser_modeler.set_skip_sampling( true );
		if ( ! minimize_single_res ) erraser_modeler.set_minimize_res( moving_res );

		erraser_modeler.apply( pose );

		move_type = "erraser";

	}

	if ( did_mutation ) move_type += "-mut";
	if ( just_min_after_mutation_ ) move_type = "mut";

	return true;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
stepwise_monte_carlo()
{

  using namespace core::pose;
  using namespace core::chemical;
  using namespace core::kinematics;
  using namespace core::scoring;
  using namespace core::io::silent;
  using namespace core::optimization;
  using namespace core::import_pose;
  using namespace core::pose::full_model_info;
	using namespace protocols::swa;
	using namespace protocols::swa::rna;
	using namespace protocols::swa::monte_carlo;
	using namespace protocols::moves;

	clock_t const time_start( clock() );

	std::string const fasta_file = option[ in::file::fasta ]()[1];
	core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( fasta_file )[1];
	std::string const desired_sequence = fasta_sequence->sequence();

	// read starting pose(s) from disk
	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );
	pose::Pose pose;
	import_pose::pose_from_pdb( pose, *rsd_set, option[ in::file::s ][1] );
	protocols::rna::figure_out_reasonable_rna_fold_tree( pose );
	protocols::rna::virtualize_5prime_phosphates( pose );

	protocols::viewer::add_conformation_viewer ( pose.conformation(), "current", 800, 800 );

	// native pose
	PoseOP native_pose;
	if (option[ in::file::native ].user() ) {
		native_pose = PoseOP( new Pose );
		import_pose::pose_from_pdb( *native_pose, *rsd_set, option[ in::file::native ]() );
	}

	// output silent file.
	std::string const silent_file = option[ out::file::silent ]();
	SilentFileData silent_file_data;

	//SubToFull (minimal object needed for add/delete)
	utility::vector1< Size > input_res_list = option[ input_res ]();
	utility::vector1< Size > start_moving_res_list /*blank*/;
	utility::vector1< Size > cutpoint_open_in_full_model  = option[ cutpoint_open ]();
	if ( input_res_list.size() != pose.total_residue() ) utility_exit_with_message( "input_res size does not match pose size" );
	FullModelInfoOP full_model_info_op =	new FullModelInfo( input_res_list, start_moving_res_list, desired_sequence, cutpoint_open_in_full_model );
	pose.data().set( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO, full_model_info_op );
	update_pdb_info_from_full_model_info( pose ); // for output pdb or silent file -- residue numbering.

	//setup rmsd res as everything to be sampled.
	utility::vector1< Size > rmsd_res;
	for ( Size i = 1; i <= desired_sequence.size(); i++ ) if ( !input_res_list.has_value(i) ) rmsd_res.push_back( i );

	// scorefunction
	core::scoring::ScoreFunctionOP scorefxn;
	if ( option[ score::weights ].user() ) scorefxn = getScoreFunction();
	else  scorefxn = ScoreFunctionFactory::create_score_function( "rna/rna_res_level_energy.wts" );

	// mover setup
	RNA_DeleteMoverOP rna_delete_mover = new RNA_DeleteMover;
	rna_delete_mover->set_minimize_scorefxn( scorefxn );

	RNA_AddMoverOP rna_add_mover = new RNA_AddMover( rsd_set, scorefxn );
	rna_add_mover->set_start_added_residue_in_aform( false );
	rna_add_mover->set_presample_added_residue(  true );
	rna_add_mover->set_presample_by_swa(  true );
	//	rna_add_mover->set_minimize_all_rebuilt_res( ! option[ minimize_single_res ]() );
	rna_add_mover->set_num_random_samples( option[ num_random_samples ]() );

	RNA_AddOrDeleteMoverOP rna_add_or_delete_mover = new RNA_AddOrDeleteMover( rna_add_mover, rna_delete_mover );
	rna_add_or_delete_mover->set_sample_res( option[ sample_res ]() );

	// final setup
	Pose start_pose = pose;
	Size num_struct = option[ out::nstruct ]();
	std::string move_type;
	Real const add_delete_frequency_ = option[ add_delete_frequency ]();
	Real const minimize_single_res_frequency_ = option[ minimize_single_res_frequency ]();
	bool success;

	// main loop
	for ( Size n = 1; n <= num_struct; n++ ) {

		pose = start_pose;

		MonteCarloOP monte_carlo_ = new MonteCarlo( pose, *scorefxn, option[ temperature]() );

		Size k( 0 );

		scorefxn->show( TR, pose );

		while (  k <= option[ cycles ]() ){

			bool success( true );
			bool const minimize_single_res = ( RG.uniform() <= minimize_single_res_frequency_ );

			if ( RG.uniform() < add_delete_frequency_ ){
				rna_add_or_delete_mover->set_minimize_all_rebuilt_res(  ! minimize_single_res );
				rna_add_or_delete_mover->apply( pose, move_type );
			} else {
				// later make this an actual class!
				success = apply_swa_mover( pose, scorefxn, minimize_single_res, move_type );
			}

			if ( !success ) continue;

			k++;
			scorefxn->show( TR, pose );

			if ( minimize_single_res ) move_type += "-minsngl";
			monte_carlo_->boltzmann( pose, move_type );

			// following can be removed later.
			TR << "Monte Carlo accepted? " << monte_carlo_->mc_accepted() << std::endl;
			TR << "Score: " << (*scorefxn)( pose ) << "  Num missing residues: " << get_number_missing_residues( pose );
			if ( native_pose ) TR << "  RMSD: " <<  get_all_atom_rmsd( pose, *native_pose, rmsd_res );
			TR << std::endl;
			monte_carlo_->show_counters();

		}

		monte_carlo_->recover_low( pose );
		monte_carlo_->show_counters();

		// output to silent filter
		BinaryRNASilentStruct s( pose, "S_"+lead_zero_string_of( n, 6 ) );
		s.add_energy( "missing", get_number_missing_residues( pose ) );
		if ( native_pose ) 	s.add_energy( "rms", get_all_atom_rmsd( pose, *native_pose, rmsd_res ) );
		silent_file_data.write_silent_struct( s, silent_file, false /*score_only*/ );

	}


}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	clock_t const my_main_time_start( clock() );

	stepwise_monte_carlo();

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

	NEW_OPT( input_res, "input residue", blank_size_vector );
	NEW_OPT( cutpoint_open, "cutpoint open (defined in full model)", blank_size_vector );
	NEW_OPT( cycles, "Number of Monte Carlo cycles", 50 );
	NEW_OPT( temperature, "Monte Carlo temperature", 1.0 );
	NEW_OPT( add_delete_frequency, "Frequency of add/delete vs. resampling", 0.5 );
	NEW_OPT( just_min_after_mutation_frequency, "After a mutation, how often to just minimize (without further sampling the mutated residue)", 0.5 );
	NEW_OPT( minimize_single_res_frequency, "Frequency with which to minimize the residue that just got rebuilt, instead of all", 0.0 );
	NEW_OPT( allow_internal_moves, "Allow moves in which internal cutpoints are created to allow ERRASER rebuilds", false );
	NEW_OPT( num_random_samples, "Number of samples from swa residue sampler before minimizing best", 1 );
	NEW_OPT( sample_res, "specify particular residues that should be rebuild (as opposed to all missing in starting PDB)", blank_size_vector );

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



