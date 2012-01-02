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
// AUTO-REMOVED #include <core/scoring/rna/RNA_Util.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintIO.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>
#include <core/chemical/ResidueTypeSet.hh>
// AUTO-REMOVED #include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>

// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunctionFactory.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Jump.hh>
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinaryRNASilentStruct.hh>
// AUTO-REMOVED #include <core/io/silent/RNA_SilentStruct.hh>
// AUTO-REMOVED #include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/pose/annotated_sequence.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/rna/RNA_HelixAssembler.hh>

#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/pose/util.hh>
//#include <basic/database/open.hh>
#include <devel/init.hh>
#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>
// AUTO-REMOVED #include <utility/io/ozstream.hh>
// AUTO-REMOVED #include <utility/io/izstream.hh>
#include <utility/exit.hh>

#include <numeric/xyzVector.hh>
// AUTO-REMOVED #include <numeric/conversions.hh>

// AUTO-REMOVED #include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ headers
//#include <cstdlib>
// AUTO-REMOVED #include <fstream>
#include <iostream>
#include <string>

//silly using/typedef

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

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
OPT_KEY( String,  seq )

/////////////////////////////////////////////////
void
rna_build_helix_test(){

 	using namespace core::scoring;
 	using namespace core::chemical;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose;
	using namespace core::kinematics;
	using namespace core::io::silent;
	using namespace protocols::rna;

	std::string full_sequence;
	if ( option[ in::file::fasta ].user() ) {
		std::string const fasta_file = option[ in::file::fasta ]()[1];
		core::sequence::SequenceOP fasta_sequence = core::sequence::read_fasta_file( fasta_file )[1];
		full_sequence = fasta_sequence->sequence();
	} else if ( option[ seq ].user() ){
		full_sequence = option[ seq ]();
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

	ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( RNA );

	pose::Pose pose;
	pose::make_pose_from_sequence( pose, full_sequence,	*rsd_set );

	RNA_HelixAssembler rna_helix_assembler;
	//rna_helix_assembler.random_perturbation( true );
	protocols::viewer::add_conformation_viewer( pose.conformation(), "current", 400, 400 );

	Size const nstruct = option[ out::nstruct ];

	for (Size n = 1; n <= nstruct; n++ ) {
		rna_helix_assembler.apply( pose, full_sequence );

		std::string const tag( "S_"+lead_zero_string_of(n, 3) );
		if ( output_silent ) {
			BinaryRNASilentStruct s( pose, tag );
			silent_file_data.write_silent_struct( s, silent_file, false /*write score only*/ );
		}
		pose.dump_pdb( full_sequence+".pdb" );
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
	using namespace basic::options;

	NEW_OPT( seq, "Input sequence", "" );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);


	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );


}
