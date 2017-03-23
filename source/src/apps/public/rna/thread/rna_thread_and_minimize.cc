// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief


// libRosetta headers
#include <protocols/rna/movers/RNAThreadAndMinimizeMover.hh>

#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/pose/rna/util.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <utility/string_util.hh>
#include <protocols/rna/denovo/util.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/viewer/viewers.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <core/pose/Pose.hh>
#include <protocols/simple_moves/AddConstraintsToCurrentConformationMover.hh>

#include <devel/init.hh>

#include <core/io/pdb/pdb_writer.hh>

#include <utility/vector1.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ headers
//#include <cstdlib>
#include <utility/io/ozstream.hh>
#include <fstream>
#include <iostream>
#include <string>

// option key includes

#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

#include <utility/excn/Exceptions.hh>

using namespace core;
using namespace basic;
using namespace protocols;
using namespace basic::options::OptionKeys;

using utility::vector1;

typedef  numeric::xyzMatrix< Real > Matrix;


//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( String, seq )
OPT_KEY( String, mutation_list )
OPT_KEY( String, insertion_list )
OPT_KEY( String, input_sequence_type )
OPT_KEY( Boolean, long_strategy )

///////////////////////////////////////////////////////////////////////
/*void
setup_mask( utility::vector1< sequence::SequenceOP > const & sequences_from_alignment, ObjexxFCL::FArray1D_bool & sequence_mask )
{
using namespace basic::options;
using namespace basic::options::OptionKeys;

utility::vector1< std::string > sequences;
for ( Size n=1; n <= sequences_from_alignment.size(); n++ )  sequences.push_back( sequences_from_alignment[n]->sequence() );

Size const alignment_length = sequences[1].size();
sequence_mask = ObjexxFCL::FArray1D_bool( alignment_length, false );

///////////////////////////////////////////////////////////////
// First pass -- rule out any part that is gapped in the alignment.
for ( Size i = 1; i <= alignment_length; i++ ) {

bool found_a_gap( false );
for ( Size n = 1; n <= sequences.size(); n++ ) {
if ( sequences[n][i-1] == '-' ) {
found_a_gap = true;
break;
}
}

if ( !found_a_gap ) {
sequence_mask( i ) = true;
}

}

// Save into a file.
Size count( 0 );
if ( option[ sequence_mask_file ].user() ) {
utility::io::ozstream out( option[ OptionKeys::sequence_mask_file ]() );
for ( Size i = 1; i <= alignment_length; i++ ) {
if ( sequences[1][i-1] == '-' ) continue; //Assume native numbering
count++;
out << sequence_mask(i);
}
out << std::endl;
out.close();
}
}*/

//////////////////////////////////////////////////////////////////////////////////////////////////////
inline std::string
remove_dashes( std::string const & s ){
	std::string s_out( "" );
	for ( Size i = 0; i < s.size(); i++ ) {
		if ( s[i] != '-' ) s_out += s[i];
	}
	return s_out;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
// Find which template to use, based on input pdb name, or, if that fails, the template sequence.
Size
figure_out_which_sequence_is_template(
	utility::vector1< sequence::SequenceOP > const & sequences_from_alignment,
	pose::Pose const & template_pose,
	std::string const & template_file_name ) {

	Size which_sequence( 0 );

	for ( Size n=1; n<=sequences_from_alignment.size(); n++ ) {
		if ( sequences_from_alignment[n]->id() == template_file_name ) {
			which_sequence = n;
			break;
		}
	}

	//////////////////////////////////////////////////////////////////
	// If could not find template based on name, try based on sequence...
	std::string const template_sequence = template_pose.sequence();
	if ( which_sequence == 0 ) {
		for ( Size n=1; n<=sequences_from_alignment.size(); n++ ) {
			if ( sequences_from_alignment[n]->ungapped_sequence() == template_sequence ) {
				which_sequence = n;
				break;
			}
		}
		if ( which_sequence > 0 )  std::cout  << "Using " << sequences_from_alignment[which_sequence]->id() << " from alignment FASTA file, as its sequence corresponds to input template PDB!" << std::endl;
	}

	if ( which_sequence == 0 )  utility_exit_with_message( "Could not figure out which sequence in fasta file was the template based on name or on sequence!!" );

	if ( which_sequence == 1 ) std::cout << "WARNING! WARNING! Assuming template corresponds to first sequence in the alignment FASTA file " <<  sequences_from_alignment[which_sequence]->id() <<  ". But that first sequence should be the target sequence!" << std::endl;

	runtime_assert( template_sequence == remove_dashes( sequences_from_alignment[which_sequence]->ungapped_sequence()  ) );

	return which_sequence;
}

/*void
change_sequence_accordingly(
std::string & target_sequence_from_alignment,
Size const seqpos,
std::string const & new_identity
) {
std::string temp_seq = "";
Size additional_counter = 0;
for ( Size ii = 0; ii < target_sequence_from_alignment.size(); ++ii ) {
if ( )temp_seq += target_sequence_from_alignment[ ii ];
}
}*/

///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	using namespace basic::options;

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::pose;
	using namespace core::sequence;

	////////////////////////////////////////////////
	//Read in template pdb.
	Pose pose;
	std::string template_file = option[ in::file::s ][1];
	core::chemical::ResidueTypeSetCOP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::FA_STANDARD );
	core::import_pose::pose_from_file( pose, *rsd_set, template_file , core::import_pose::PDB_file);
	core::pose::rna::figure_out_reasonable_rna_fold_tree( pose );

	////////////////////////////////////////////////
	//Read in fasta file or user-inputted sequence
	std::string template_sequence_from_alignment, target_sequence_from_alignment, mutlist;
	ObjexxFCL::FArray1D_bool sequence_mask;
	if ( option[ seq ].user() ) {

		runtime_assert( ! option[ in::file::fasta ].user() );
		template_sequence_from_alignment = pose.annotated_sequence();
		target_sequence_from_alignment = option[ seq ]();
		std::cout << "Input pose seq: " << template_sequence_from_alignment << std::endl;
		std::cout << "Input tgt seq:  " << target_sequence_from_alignment << std::endl;

		std::cout << "Clean pose seq: " << core::pose::rna::remove_bracketed(template_sequence_from_alignment) << std::endl;
		std::cout << "Clean tgt seq:  " << core::pose::rna::remove_bracketed(target_sequence_from_alignment) << std::endl;
		if ( core::pose::rna::remove_bracketed(target_sequence_from_alignment).size() != core::pose::rna::remove_bracketed(template_sequence_from_alignment).size() ) utility_exit_with_message( "Input sequence from -seq should match size of template pose specified by -s!" );

	} else if ( option[ mutation_list ].user() ) {

		template_sequence_from_alignment = pose.annotated_sequence();
		target_sequence_from_alignment = pose.annotated_sequence();
		mutlist = option[ mutation_list ].value();
		// Mutate target sequence with options in mutation list.
		/*utility::vector1< std::string > mutation_list = utility::string_split( option[ mutation_list ].value(), ' ' );
		for ( auto const & mutation : mutation_list ) {
		// expect: seqpos,identity
		std::istringstream iss( mutation );
		Size seqpos, std::string new_identity;
		std::getline( iss, seqpos, ',');
		iss >> new_identity;//std::getline( iss, identity, ',')
		//change_sequence_accordingly( target_sequence_from_alignment, seqpos, new_identity );
		}*/
	}

	protocols::rna::movers::RNAThreadAndMinimizeMoverOP rtm( new protocols::rna::movers::RNAThreadAndMinimizeMover(
		target_sequence_from_alignment,
		template_sequence_from_alignment,
		option[ long_strategy ](),
		option[ input_sequence_type ](),
		mutlist,
		option[ insertion_list ].value()
		) );

	protocols::jd2::JobDistributor::get_instance()->go( rtm );

	core::pose::rna::virtualize_5prime_phosphates( pose );

	std::string outfile( "threaded.pdb" );
	if ( option[ out::file::o ].user() ) outfile = option[ out::file::o ]();

	std::cout << "Outputting: " << outfile << std::endl;
	pose.dump_pdb( outfile );

	exit( 0 );
}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
	try {
		using namespace basic::options;

		NEW_OPT( seq, "target sequence (can include dashes). Must specify either this or -fasta. Length of this sequence must exactly equal length of template pose specified by -s", "");
		NEW_OPT( mutation_list, "Mutations to make in the original numbering of the template PDB that's been provided. Provide sequence position and new residue identity, in modomics or IUPAC or annotated format, i.e. 13,c 14,X[4SU]", "");
		NEW_OPT( insertion_list, "Space-separated list of insertions, specified by a residue number, a comma, and a sequence of residues in any format", "");
		NEW_OPT( input_sequence_type, "Format of input sequence: MODOMICS or IUPAC (or omit or use any other value to use annotated sequence format)", "" );
		NEW_OPT( long_strategy, "Make mutations one at a time, stabilizing after each", false );


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
