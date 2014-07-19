// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// This file is part of the Rosetta software suite and is made available under license.
// The Rosetta software is developed by the contributing members of the Rosetta Commons consortium.
// (C) 199x-2009 Rosetta Commons participating institutions and developers.
// For more information, see http://www.rosettacommons.org/.

/// @file
/// @brief


// libRosetta headers
#include <core/scoring/rms_util.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/pose/annotated_sequence.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/stepwise/sampling/rna/helix/RNA_HelixAssembler.hh>

#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/init/init.hh>
#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>
#include <utility/exit.hh>

#include <numeric/xyzVector.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ headers
#include <iostream>
#include <string>

//silly using/typedef

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/rna.OptionKeys.gen.hh>

#include <utility/excn/Exceptions.hh>

//Auto Headers
using namespace core;
using namespace protocols;
using namespace basic::options::OptionKeys;

using utility::vector1;

using io::pdb::dump_pdb;

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( StringVector,  seq )
OPT_KEY( Boolean,  minimize_all )
OPT_KEY( Boolean,  dump )
OPT_KEY( String,  finish_weights )

/////////////////////////////////////////////////
void
rna_build_helix_test(){

 	using namespace core::scoring;
 	using namespace core::chemical;
 	using namespace core::sequence;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose;
	using namespace core::kinematics;
	using namespace core::io::silent;
	using namespace protocols::stepwise::sampling::rna::helix;

	std::string full_sequence;
	if ( option[ in::file::fasta ].user() ) {
		std::string const fasta_file = option[ in::file::fasta ]()[1];
		SequenceOP fasta_sequence = read_fasta_file( fasta_file )[1];
		full_sequence = fasta_sequence->sequence();
	} else if ( option[ seq ].user() ){
		utility::vector1< std::string > full_sequence_segments = option[ seq ]();
		full_sequence = "";
		for (Size n = 1; n <= full_sequence_segments.size(); n++ ) full_sequence += full_sequence_segments[ n ];
	} else {
		utility_exit_with_message( "Need to specify -fasta <fasta file> or -seq <sequence>." );
	}

	std::string silent_file = "";
	bool output_silent( false );
	if ( option[ out::file::silent ].user() ) {
		silent_file = option[ out::file::silent  ]();
		output_silent = true;
	}
	SilentFileData silent_file_data;

	bool const is_use_phenix_geo = option[ basic::options::OptionKeys::rna::corrected_geo] ();
	ResidueTypeSetCAP rsd_set;
	rsd_set = ChemicalManager::get_instance()->residue_type_set( FA_RNA );

	pose::Pose pose;
	std::string sequence_to_build;
	for (Size n = 1; n <= full_sequence.size(); n++ ) if ( !is_blank_seq( full_sequence[n-1]) ) sequence_to_build.push_back( full_sequence[ n-1 ] );
	pose::make_pose_from_sequence( pose, sequence_to_build,	*rsd_set );

	RNA_HelixAssembler rna_helix_assembler;
	//rna_helix_assembler.random_perturbation( true );
	rna_helix_assembler.set_minimize_all( option[ minimize_all ]() );
	rna_helix_assembler.set_dump( option[ dump ]() );
	rna_helix_assembler.use_phenix_geo( is_use_phenix_geo );
	if ( option[ basic::options::OptionKeys::score::weights ].user() ) {
		ScoreFunctionOP scorefxn = get_score_function();
		rna_helix_assembler.set_scorefxn ( scorefxn );
	}
	// if following is blank, there won't be a 'final' minimize.
	if ( option[ finish_weights ].user() ){
		ScoreFunctionOP finish_scorefxn = ScoreFunctionFactory::create_score_function( option[ finish_weights ]()  );
		rna_helix_assembler.set_finish_scorefxn( finish_scorefxn );
	}
	rna_helix_assembler.set_model_and_remove_capping_residues( true );

	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	Size const nstruct = option[ out::nstruct ];

	std::string outfile = full_sequence+".pdb";
	if ( option[ out::file::o ].user() ) outfile = option[ out::file::o ]();

	for (Size n = 1; n <= nstruct; n++ ) {
		rna_helix_assembler.apply( pose, full_sequence );

		if ( output_silent ) {
			std::string const tag( "S_"+ObjexxFCL::lead_zero_string_of(n, 3) );
			BinarySilentStruct s( pose, tag );
			silent_file_data.write_silent_struct( s, silent_file, false /*write score only*/ );
		}

		pose.dump_pdb( outfile );
	}

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace basic::options;

	rna_build_helix_test();
	protocols::viewer::clear_conformation_viewers();
	exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
try {
	using namespace basic::options;

	std::cout << std::endl << "Basic usage:  " << argv[0] << "  -seq <sequence of first strand> <sequence of second strand>   -o <name of output pdb file> " << std::endl;
	std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

	utility::vector1< std::string > blank_string_vector;
	NEW_OPT( seq, "Input sequence", blank_string_vector );
	NEW_OPT( minimize_all, "minimize all torsions in response to each base pair addition", false );
	NEW_OPT( dump, "dump intermediate pdbs", false );
	NEW_OPT( finish_weights, "[ optional ] score function to do a second minimize with after each base pair addition", "" );

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
	return -1;
}

}
