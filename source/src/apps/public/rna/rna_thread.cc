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
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <utility/string_util.hh>

//Mmmm.. constraints.
//#include <core/scoring/constraints/CoordinateConstraint.hh>

#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/tree/Atom.hh>
#include <core/id/AtomID_Map.hh>
#include <core/id/AtomID.hh>
#include <core/id/DOF_ID.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/MoveMap.hh>

#include <core/pose/PDBInfo.hh>

#include <protocols/farna/RNA_ProtocolUtil.hh>

#include <protocols/viewer/viewers.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>

#include <devel/init.hh>

#include <core/io/pdb/pose_io.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>
#include <utility/io/izstream.hh>

#include <numeric/xyzVector.hh>
#include <numeric/conversions.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ headers
//#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

#include <utility/excn/Exceptions.hh>

using namespace core;
using namespace basic;
using namespace protocols;
using namespace basic::options::OptionKeys;

using utility::vector1;

using io::pdb::dump_pdb;

typedef  numeric::xyzMatrix< Real > Matrix;

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Adapted from easy_target_test.cc (threading routine for CASP8 proteins, 2008) to RNA
// for RNA_puzzle GIR1 in 2012.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Definition of new OptionKeys
// these will be available in the top-level OptionKey namespace:
// i.e., OPT_KEY( Type, key ) -->  OptionKey::key
// to have them in a namespace use OPT_1GRP_KEY( Type, grp, key ) --> OptionKey::grp::key
OPT_KEY( String, seq )
OPT_KEY( Integer,seq_offset )
OPT_KEY( String, sequence_mask_file )


///////////////////////////////////////////////////////////////////////
void
setup_mask(
					 ObjexxFCL::FArray1D_bool & sequence_mask,
					 utility::vector1< sequence::SequenceOP > const & sequences_from_alignment )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	utility::vector1< std::string > sequences;
	for (Size n=1; n <= sequences_from_alignment.size(); n++ )		sequences.push_back( sequences_from_alignment[n]->sequence() );

	Size const alignment_length = sequences[1].size();
	sequence_mask = ObjexxFCL::FArray1D_bool( alignment_length, false );

	///////////////////////////////////////////////////////////////
	// First pass -- rule out any part that is gapped in the alignment.
	for (Size i = 1; i <= alignment_length; i++ ) {

		bool found_a_gap( false );
		for (Size n = 1; n <= sequences.size(); n++ ){
			if ( sequences[n][i-1] == '-' ) {
				found_a_gap = true;
				break;
			}
		}

		if (!found_a_gap) {
			sequence_mask( i ) = true;
		}

	}

	///////////////////////////////////////////////////////////////
	// Second pass -- look for ungapped part that are bracketed by gapped regions
	bool look_for_fishy_gaps( false ); //this was from protein stuff, and is not actually in use for RNA.

	if (look_for_fishy_gaps ){
		Size const look_for_gap( 4 );
		for (int i = 1; i <= int(alignment_length); i++ ) {

			if ( sequence_mask(i) == false ) continue; //Don't worry about it.

			bool found_gap_before( false ), found_gap_after( false );
			for (int offset = -1 * (int)look_for_gap; offset < 0; ++offset ) {
				if ( i+offset > 1 && sequence_mask( i+offset ) == false ) {
					found_gap_before = true;
					break;
				}
			}
			if (!found_gap_before) continue;

			for (int offset = 1; Size(offset) <= look_for_gap; offset ++ ) {
				if ( i+offset <= int(alignment_length) && sequence_mask( i+offset ) == false ) {
					found_gap_after = true;
					break;
				}
			}
			if (!found_gap_after) continue;

			std::cout << "MASK: Region that is not nominally gapped but looks fishy: " << i << std::endl;
			sequence_mask( i ) = false;
		}
	}

	// Save into a file.
	Size count( 0 );
	if ( option[ sequence_mask_file ].user() ) {
		utility::io::ozstream out( option[ OptionKeys::sequence_mask_file ]() );
		for (Size i = 1; i <= alignment_length; i++ ) {
			if (sequences[1][i-1] == '-' ) continue; //Assume native numbering
			count++;
			out << sequence_mask(i);
		}
		out << std::endl;
		out.close();
	}

}

////////////////////////////////////////////////////////////////////////////
void
setup_alignment_map( std::map< Size, Size > & mapping,
										 std::string const & sequence_from_alignment ){
	Size count( 0 );
	for (Size i = 1; i <= sequence_from_alignment.size(); i++ ) {
		if ( sequence_from_alignment[ i-1 ] != '-' ) {
			count++;
			mapping[ i ] = count;
		} else {
			mapping[ i ] = 0;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
std::string
remove_dashes( std::string const & s ){
	std::string s_out( "" );
	for ( Size i = 0; i < s.size(); i++ ){
		if (s[i] != '-') s_out += s[i];
	}
	return s_out;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void
prepare_threaded_model(
											 pose::Pose & pose,
											 std::string const & target_sequence_from_alignment,
											 std::string const & template_sequence_from_alignment,
											 ObjexxFCL::FArray1D_bool & sequence_mask /* should make this optional! */,
											 Size offset = 0
											 )
{

	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::pose;
	using namespace protocols::farna;
	using namespace utility;

	// following creates pose from scratch! Assumes that we need to add/delete, and it changes the
	// fold-tree -- both are probably bad ideas for a general threading routine.
	// Would be better to add or delete residues as needed, trying to be smart about fold_tree [how?].
	// And better to skip this part if add/delete not necessary -- i.e., just isolate the "mutate sequence" part.
	Pose const template_pose = pose;
	pose.clear();

	Size const alignment_length = target_sequence_from_alignment.size();
	runtime_assert( alignment_length ==  template_sequence_from_alignment.size() );

	if ( sequence_mask.size() == 0 ) sequence_mask = ObjexxFCL::FArray1D_bool( alignment_length, true );
	runtime_assert( alignment_length ==  sequence_mask.size() );

	std::map <Size, Size> template_alignment2sequence, target_alignment2sequence;
	// Note that the Sequence object already has a nice function for this called resnum() -- we should just use it.
	setup_alignment_map( template_alignment2sequence, template_sequence_from_alignment );
	setup_alignment_map( target_alignment2sequence,   target_sequence_from_alignment );

	std::string target_sequence;
	utility::vector1< Size > working_res;
	for (Size i = 1; i <= alignment_length; i++ ){

		if ( !sequence_mask( i ) )  continue;

		Size const pdb_number( template_alignment2sequence[i] );

		if ( i == 1 ||  ( sequence_mask( i-1 ) && pdb_number > 1 && !template_pose.fold_tree().is_cutpoint( pdb_number-1 ) )  ){
			pose.append_residue_by_bond( template_pose.residue( pdb_number ) );
		} else {
			pose.append_residue_by_jump( template_pose.residue( pdb_number ), pose.total_residue() );
		}

		//std::cout << "creating alignment ==>  template:" << pdb_number << "  target: " << target_alignment2sequence[i] << std::endl;
		//std::cout << target_sequence_from_alignment[i-1] << ' ' << template_pose.residue( pdb_number ).name1() << std::endl;

		target_sequence += target_sequence_from_alignment[i-1];
		working_res.push_back( target_alignment2sequence[i] + offset); // will be used for PDB numbering
	}

	// fix up numbering.
	PDBInfoOP pdb_info = new PDBInfo( pose );
	pdb_info->set_numbering( working_res );
	pose.pdb_info( pdb_info );

	std::string current_sequence = pose.sequence();
	std::cout << "TARGET: " << target_sequence << std::endl;
	std::cout << "CURRENT: " << current_sequence << std::endl;
	utility::vector1< Size > changed_pos;

	std::cout << "WORKING_RES " << make_tag( working_res ) << std::endl;
	std::cout << "WORKING_RES " << make_tag_with_dashes( working_res ) << std::endl;

	/////////////////////
	//Mutate sequence
	/////////////////////
	utility::vector1< int> changed_pos_working;
	for (Size i = 1; i <= target_sequence.size(); i++ ){

		char const new_seq = target_sequence[i-1];
		if ( mutate_position( pose, i, new_seq ) ){
			changed_pos.push_back( i );
			changed_pos_working.push_back( working_res[i] );
		}

	}

	std::cout << "Changed residues (without offset applied): " << make_tag_with_dashes( changed_pos ) << std::endl;
	std::cout << "Changed residues (in residue numbering): "    << make_tag_with_dashes( changed_pos_working ) << std::endl;

}


//////////////////////////////////////////////////////////////////////////////////////////////////////
// Find which template to use, based on input pdb name, or, if that fails, the template sequence.
Size
figure_out_which_sequence_is_template( utility::vector1< sequence::SequenceOP > const & sequences_from_alignment,
																			 pose::Pose const & template_pose,
																			 std::string const & template_file_name ){

	Size which_sequence( 0 );

	for( Size n=1; n<=sequences_from_alignment.size(); n++ ){
		if ( sequences_from_alignment[n]->id() == template_file_name ) {
			which_sequence = n;
			break;
		}
	}

	//////////////////////////////////////////////////////////////////
	// If could not find template based on name, try based on sequence...
	std::string const template_sequence = template_pose.sequence();
	if ( which_sequence == 0 ){
		for( Size n=1; n<=sequences_from_alignment.size(); n++ ){
			if ( sequences_from_alignment[n]->ungapped_sequence() == template_sequence ) {
				which_sequence = n;
				break;
			}
		}
		if ( which_sequence > 0 )  std::cout  << "Using " << sequences_from_alignment[which_sequence]->id() << " from alignment FASTA file, as its sequence corresponds to input template PDB!" << std::endl;
	}

	if ( which_sequence == 0 )		utility_exit_with_message( "Could not figure out which sequence in fasta file was the template based on name or on sequence!!" );

	if ( which_sequence == 1 ) std::cout << "WARNING! WARNING! Assuming template corresponds to first sequence in the alignment FASTA file " <<  sequences_from_alignment[which_sequence]->id() <<  ". But that first sequence should be the target sequence!" << std::endl;

	runtime_assert( template_sequence == remove_dashes( sequences_from_alignment[which_sequence]->ungapped_sequence()  ) );

	return which_sequence;

}


////////////////////////////////////////////////////////////////////////////
void
rna_thread_test(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::pose;
	using namespace core::sequence;

	////////////////////////////////////////////////
	//Read in template pdb.
	Pose pose;
	std::string template_file = option[ in::file::s ][1];
	core::chemical::ResidueTypeSetCAP rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );
	core::import_pose::pose_from_pdb( pose, *rsd_set, template_file );
	protocols::farna::figure_out_reasonable_rna_fold_tree( pose );

	////////////////////////////////////////////////
	//Read in fasta file or user-inputted sequence
	std::string template_sequence_from_alignment, target_sequence_from_alignment;
	ObjexxFCL::FArray1D_bool sequence_mask;
	if ( option[ seq ].user() ) {

		runtime_assert( ! option[ in::file::fasta ].user() );
		template_sequence_from_alignment = pose.sequence();
		target_sequence_from_alignment = option[ seq ]();
		if ( target_sequence_from_alignment.size() != template_sequence_from_alignment.size() ) utility_exit_with_message( "Input sequence from -seq should match size of template pose specified by -s!" );

	} else {

		std::string const fasta_file ( option[ in::file::fasta ]()[1] );
		utility::vector1< std::string > sequences, pdb_names;
		utility::vector1< SequenceOP > sequences_from_alignment = read_fasta_file( fasta_file );
		Size const which_sequence = figure_out_which_sequence_is_template( sequences_from_alignment, pose, template_file );
		template_sequence_from_alignment = sequences_from_alignment[ which_sequence ]->sequence();
		target_sequence_from_alignment  = sequences_from_alignment[ 1 ]->sequence();

		// We can/should replace above with standard Rosetta reader!
		// And if we use Sequence object, that will return useful info like ungapped_sequence()...
		// Also cool -- Sequence can include start numbering [= offset + 1].

		//Figure out a mask ... this is currently not doing anything if we just have a template and a target sequence.
		setup_mask( sequence_mask, sequences_from_alignment );

	}

	prepare_threaded_model( pose, target_sequence_from_alignment,
													template_sequence_from_alignment,
													sequence_mask, option[ seq_offset ]() );
	protocols::farna::virtualize_5prime_phosphates( pose );

	std::string outfile( "threaded.pdb" );
	if ( option[ out::file::o ].user() ) outfile = option[ out::file::o ]();

	std::cout << "Outputting: " << outfile << std::endl;
	pose.dump_pdb( outfile );

}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{

	using namespace basic::options;
	rna_thread_test();
	exit( 0 );

}


///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{
try {
	using namespace basic::options;

	//Uh, options?
	NEW_OPT( seq, "target sequence (can include dashes). Must specify either this or -fasta. Length of this sequence must exactly equal length of template pose specified by -s", "");
	NEW_OPT( seq_offset, "Integer to add to all residue numbers in output PDB", 0 );
	NEW_OPT( sequence_mask_file, "Output name for sequence mask file (not in use at the moment)", "sequence_mask.txt" );
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
}
}
