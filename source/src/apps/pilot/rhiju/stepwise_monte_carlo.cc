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
#include <protocols/swa/monte_carlo/RNA_AddMover.hh>
#include <protocols/swa/monte_carlo/SubToFullInfo.hh>
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
	using namespace protocols::swa::rna;
	using namespace protocols::swa::monte_carlo;
	using namespace protocols::moves;

	clock_t const time_start( clock() );

	// read starting pose(s) from disk
	ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	std::string const fasta_file = option[ in::file::fasta ]()[1];
	core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( fasta_file )[1];
	std::string const desired_sequence = fasta_sequence->sequence();

	pose::Pose pose;
	import_pose::pose_from_pdb( pose, *rsd_set, option[ in::file::s ][1] );
	protocols::rna::figure_out_reasonable_rna_fold_tree( pose );
	protocols::rna::virtualize_5prime_phosphates( pose );
	utility::vector1< Size > input_res_list = option[ input_res ]();

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

	//SubToFull (minimal object needed for add/delete!)
	utility::vector1< Size > start_moving_res_list /*blank*/;
	std::map< Size, Size > sub_to_full;
	for ( Size i = 1; i <= input_res_list.size(); i++ ) sub_to_full[i] = input_res_list[i];
	utility::vector1< Size > cutpoints_in_full_pose /*blank*/;
	SubToFullInfoOP sub_to_full_info_op =	new SubToFullInfo(  sub_to_full, start_moving_res_list,	desired_sequence, cutpoints_in_full_pose );
	pose.data().set( core::pose::datacache::CacheableDataType::SUB_TO_FULL_INFO, sub_to_full_info_op );

	//setup rmsd res as everything sampled.
	std::list< Size > rmsd_res;
	for ( Size i = 1; i <= desired_sequence.size(); i++ ) if ( !Contain_seq_num(i, input_res_list ) ) rmsd_res.push_back( i );

	// scorefunction
	core::scoring::ScoreFunctionOP scorefxn = getScoreFunction();

	Pose start_pose = pose;

	Size num_struct = option[ out::nstruct ]();

	for ( Size n = 1; n <= num_struct; n++ ) {

		pose = start_pose;

		// following is hard-wired for now.
		// first build scratch residue A7.
		RNA_AddMover add_mover( rsd_set, scorefxn );
		add_mover.apply( pose, 6 /*this is NOT global numbering, but perhaps it should be*/, CHAIN_TERMINUS_5PRIME /*prepend*/);

		{
			// then sample residue A7
			StepWiseRNA_Modeler stepwise_rna_modeler( 6, scorefxn );
			stepwise_rna_modeler.set_choose_random( true );
			stepwise_rna_modeler.set_force_centroid_interaction( true );
			//	stepwise_rna_modeler->set_use_phenix_geo ( option[ basic::options::OptionKeys::rna::corrected_geo ]() );
			stepwise_rna_modeler.apply( pose );
		}

		// then build in C6
		add_mover.apply( pose, 6 /*this is NOT global numbering, but perhaps it should be*/, CHAIN_TERMINUS_5PRIME /*prepend*/);
		add_variant_type_to_pose_residue( pose, CUTPOINT_LOWER, 5 );
		add_variant_type_to_pose_residue( pose, CUTPOINT_UPPER, 6 );

		// then sample residue C6 (with CCD loop closure)
		{
			// then sample residue A7
			StepWiseRNA_Modeler stepwise_rna_modeler( 6, scorefxn );
			stepwise_rna_modeler.set_choose_random( true );
			stepwise_rna_modeler.set_force_centroid_interaction( true );
			//	stepwise_rna_modeler->set_use_phenix_geo ( option[ basic::options::OptionKeys::rna::corrected_geo ]() );
			stepwise_rna_modeler.apply( pose );
		}

		// output to silent file
		BinaryRNASilentStruct s( pose, "S_"+lead_zero_string_of( n, 6 ) );
		// following is a quick hack.
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



