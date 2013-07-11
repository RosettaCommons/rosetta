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
#include <core/scoring/rna/RNA_Util.hh>
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

OPT_KEY( IntegerVector, input_res )
OPT_KEY( Integer, cycles )
OPT_KEY( Boolean, minimize_single_res )
OPT_KEY( Real, temperature )


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
apply_swa_mover( pose::Pose & pose,
								 core::scoring::ScoreFunctionOP scorefxn,
								 std::string & move_type ){

	using namespace protocols::swa::monte_carlo;
	using namespace protocols::swa::rna;
	using namespace core::pose::full_model_info;

	move_type = "swa";

	utility::vector1< Size > possible_res;
	get_potential_resample_residues( pose, possible_res );

	if ( possible_res.size() == 0 ) return false;

	Size const remodel_res = RG.random_element( possible_res );

	swa::rna::StepWiseRNA_Modeler stepwise_rna_modeler( remodel_res, scorefxn );
	stepwise_rna_modeler.set_choose_random( true );
	stepwise_rna_modeler.set_force_centroid_interaction( true );
	//	stepwise_rna_modeler->set_use_phenix_geo ( option[ basic::options::OptionKeys::rna::corrected_geo ]() );

	if ( ! option[ minimize_single_res]() ) stepwise_rna_modeler.set_minimize_res( nonconst_full_model_info_from_pose( pose ).moving_res_list() );

	stepwise_rna_modeler.apply( pose );


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

	protocols::viewer::add_conformation_viewer ( pose.conformation(), "current", 400, 400 );

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
	utility::vector1< Size > cutpoint_open_in_full_model /*blank*/;
	FullModelInfoOP full_model_info_op =	new FullModelInfo(  input_res_list, start_moving_res_list,	desired_sequence, cutpoint_open_in_full_model );
	pose.data().set( core::pose::datacache::CacheableDataType::FULL_MODEL_INFO, full_model_info_op );
	update_pdb_info_from_sub_to_full( pose ); // for output pdb or silent file -- residue numbering.

	//setup rmsd res as everything to be sampled.
	std::list< Size > rmsd_res;
	for ( Size i = 1; i <= desired_sequence.size(); i++ ) if ( !Contain_seq_num(i, input_res_list ) ) rmsd_res.push_back( i );

	// scorefunction
	core::scoring::ScoreFunctionOP scorefxn = getScoreFunction();

	// mover setup
	RNA_DeleteMoverOP rna_delete_mover = new RNA_DeleteMover;

	RNA_AddMoverOP rna_add_mover = new RNA_AddMover( rsd_set, scorefxn );
	rna_add_mover->set_start_added_residue_in_aform( false );
	rna_add_mover->set_presample_added_residue(  true );
	rna_add_mover->set_presample_by_swa(  true );
	rna_add_mover->set_minimize_all_rebuilt_res( ! option[ minimize_single_res ]() );
	RNA_AddOrDeleteMoverOP rna_add_or_delete_mover = new RNA_AddOrDeleteMover( rna_add_mover, rna_delete_mover );

	// final setup
	Pose start_pose = pose;
	Size num_struct = option[ out::nstruct ]();
	std::string move_type;
	Real const add_delete_frequency( 0.2 );
	bool success;

	// main loop
	for ( Size n = 1; n <= num_struct; n++ ) {

		pose = start_pose;
		MonteCarloOP monte_carlo_ = new MonteCarlo( pose, *scorefxn, option[ temperature]() );

		Size k( 0 );

		while (  k <= option[ cycles ]() ){

			bool success( true );
			if ( RG.uniform() < add_delete_frequency ){
				rna_add_or_delete_mover->apply( pose, move_type );
			} else {
				// later make this an actual class!
				success = apply_swa_mover( pose, scorefxn, move_type );
			}

			if ( !success ) continue;
			k++;
			monte_carlo_->boltzmann( pose, move_type );

			// following can be removed later.
			std::cout << "Monte Carlo accepted? " << monte_carlo_->mc_accepted() << std::endl;
			monte_carlo_->show_counters();
			pose.dump_pdb( "latest.pdb" );

		}

		monte_carlo_->recover_low( pose );
		monte_carlo_->show_counters();

		// output to silent file
		BinaryRNASilentStruct s( pose, "S_"+lead_zero_string_of( n, 6 ) );
		// following is a quick hack -- need to actually take into account moving suites
		if ( native_pose ) 	s.add_energy( "rms", all_atom_rmsd( *native_pose, pose, rmsd_res ) );
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
	NEW_OPT( cycles, "Number of Monte Carlo cycles", 50 );
	NEW_OPT( temperature, "Monte Carlo temperature", 1.0 );
	NEW_OPT( minimize_single_res, "Minimize the residue that just got rebuilt, instead of all", false );

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



