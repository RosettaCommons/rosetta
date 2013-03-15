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
#include <core/scoring/rms_util.hh>
#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueMatcher.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ResidueSelector.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/chemical/VariantType.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/sequence/util.hh>

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

#include <core/io/silent/SilentFileData.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/task/operation/TaskOperations.hh>

#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/optimization/MinimizerOptions.hh>

#include <protocols/rna/RNA_ProtocolUtil.hh>

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
OPT_KEY( Boolean, chi_stats )
OPT_KEY( Boolean, copy_native_chi )
OPT_KEY( Boolean, rhiju_fold_tree )
OPT_KEY( Boolean, use_native_CA )
OPT_KEY( Boolean, vary_geometry )
OPT_KEY( Boolean, mask_loop )
OPT_KEY( Boolean, cst_relax )
OPT_KEY( String, sequence_mask_file )
OPT_KEY( String, secstruct_file )
OPT_KEY( String, cst_file )
OPT_KEY( Integer, seq_offset )

//////////////////////////////////////////////////////////////////
// Following could probably be replaced with Rosetta's
// more general fasta reader
//////////////////////////////////////////////////////////////////
Size
read_alignment_fasta_file(
				utility::vector1< std::string >  & sequences,
				utility::vector1< std::string >  & pdb_names,
				std::string const & fasta_file
 )
{
	sequences.clear();
	pdb_names.clear();

	utility::io::izstream data_stream( fasta_file );

	if ( data_stream.fail() )		utility_exit_with_message( "Could not find fasta file: " + fasta_file  );

	std::string line;
	std::string sequence = "";

	while( getline(data_stream, line) 	) {

		while ( line.size() == 0 || line[0] == '>' || line[0] == ' '  ) {

			//Finished a sequence?
			if ( sequence.size() > 0 ) 	sequences.push_back( sequence );

			//Is this a tag with the pdb names?
			if ( line.size() > 0 && line[0] == '>' ) {
				pdb_names.push_back( line.substr( 1, line.size()-1 ) );
			}

			sequence = "";
			getline(data_stream, line);
			if ( data_stream.fail() ) break;
		}

		if ( data_stream.fail() ) break;

		//We're inside a sequence.
		sequence += line;
	}
	if (sequence.size() > 0 ) sequences.push_back( sequence );

	if ( sequences.size() == 0 )		utility_exit_with_message( "Could not find sequence inside file" + fasta_file  );

	Size const alignment_length = sequences[1].size();
	for (Size n = 1; n <= sequences.size(); n++ ){
		std::cout << pdb_names[n] << "  " << sequences[n] << std::endl;
		assert( sequences[n].size() == 	alignment_length );
	}

	return alignment_length;
}

///////////////////////////////////////////////////////////////////////
void
setup_mask(
					 ObjexxFCL::FArray1D_bool & sequence_mask,
					 utility::vector1< std::string > const & sequences )
{

	sequence_mask = false;

	Size const alignment_length = sequences[1].size();

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
	bool look_for_fishy_gaps( false ); //this was from protein stuff.

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

	////////////////////////////////////////////////
	// Debug output
	std::cout << "SETUP MASK ==>  exclude ";
	for (Size i = 1; i <= alignment_length; i++ ) {
		if ( sequence_mask(i) == false ) std::cout << i << " " ;
	}
	std::cout << std::endl;

	////////////////////////////////////////////////
	// Save into a file.
	Size count( 0 );
	utility::io::ozstream out( options::option[ options::OptionKeys::sequence_mask_file ]  );
	for (Size i = 1; i <= alignment_length; i++ ) {
		if (sequences[1][i-1] == '-' ) continue; //Assume native numbering
		count++;
		out << sequence_mask(i);
	}
	out << std::endl;
	out.close();

}

////////////////////////////////////////////////////////////////////////////
void
setup_alignment_map( utility::vector1< std::map< Size, Size > > & alignment2sequence, utility::vector1 < std::string > const & sequences )
{

	//Figure out mapping of alignment to each sequence.
	Size const alignment_length( sequences[1].size() );
	for (Size n = 1; n <= sequences.size(); n++ ){
		Size count( 0 );
		std::map< Size, Size > mapping;
		for (Size i = 1; i <= alignment_length; i++ ) {
			if ( sequences[n][i-1] != '-' ) {
				count++;
				mapping[ i ] = count;
			} else {
				mapping[ i ] = 0;
			}
		}
		alignment2sequence.push_back( mapping );
	}

}


/////////////////////////////////////////////////////////////////////////////////
// Would be useful utility to put somewhere else in Rosetta...
//  perhaps will need for silent output.
/////////////////////////////////////////////////////////////////////////////////
std::string
make_tag_with_dashes( utility::vector1< Size > working_res ){

	using namespace ObjexxFCL;
	std::string tag = "";

	utility::vector1< std::pair<Size,Size> > working_res_segments;
	Size start_segment = working_res[1];
	Size last_res = working_res[1];

	for (Size n = 2; n<= working_res.size(); n++ ){
		if ( working_res[n] != last_res+1 ){
			working_res_segments.push_back( std::make_pair( start_segment, last_res ) );
			start_segment = working_res[n];
		}
		last_res = working_res[n];
	}
	working_res_segments.push_back( std::make_pair( start_segment, last_res ) );

	for (Size n = 1; n <= working_res_segments.size(); n++ ){
		if ( n > 1 ) tag += " ";
		std::pair< Size, Size > const & segment = working_res_segments[n];
		if ( segment.first == segment.second ){
			tag += string_of( segment.first );
		} else{
			tag += string_of( segment.first )+"-"+string_of(segment.second);
		}
	}

	return tag;
}


/////////////////////////////////////////////////////////////////////////////////void
// Would be useful utility to put somewhere else in Rosetta...
//  perhaps will need for silent output.
/////////////////////////////////////////////////////////////////////////////////
std::string
make_tag( utility::vector1< Size > working_res ){

	using namespace ObjexxFCL;
	std::string tag = "";

	for (Size n = 1; n <= working_res.size(); n++ ){
		if ( n > 1 ) tag += " ";
		tag += string_of( working_res[n] );
	}

	return tag;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
void
prepare_full_length_start_model(
	 pose::Pose & template_pose,
	 pose::Pose & pose,
	 utility::vector1<std::string> const & sequences,
	 ObjexxFCL::FArray1D_bool const & sequence_mask,
	 utility::vector1< std::map< Size, Size > > & alignment2sequence,
	 utility::vector1< std::string > const & pdb_names,
	 std::string const & which_file
																)
{
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	pose.clear();
	Size const alignment_length( alignment2sequence[1].size() );

	// Find which alignment to use, based on input pdb name.
	Size which_sequence( 0 );
	for( Size n=1; n<=pdb_names.size(); n++ ){
		//		std::cout << pdb_names[n] << " " << which_file << std::endl;
		if ( pdb_names[n] == which_file ) {
			which_sequence = n;
			break;
		}
	}

	runtime_assert( which_sequence > 0 );

	template_pose.dump_pdb( "template.pdb" );

	ResidueTypeSet const & rsd_set( template_pose.residue(1).residue_type_set() );

	std::string const full_desired_sequence = sequences[1];
	std::string desired_sequence = "";
	utility::vector1< Size > working_res;

	for (Size i = 1; i <= alignment_length; i++ ){

		if ( !sequence_mask( i ) )  continue;
		Size const pdb_number( alignment2sequence[ which_sequence ][i] );

		if ( i == 1 || sequence_mask( i-1 ) ){
			pose.append_residue_by_bond( template_pose.residue( pdb_number ) );
		} else {
			pose.append_residue_by_jump( template_pose.residue( pdb_number ), pose.total_residue() );
		}
		desired_sequence += full_desired_sequence[i-1];

		std::cout << "creating alignment ==>  template:" << pdb_number << "  target: " << alignment2sequence[1][i] << std::endl;
		std::cout << full_desired_sequence[i-1] << ' ' << template_pose.residue( pdb_number ).name1() << std::endl;

		working_res.push_back( alignment2sequence[1][i] + option[ seq_offset ]() );
	}

	pose.dump_pdb( "before_mutations.pdb" );

	std::cout << "WORKING_RES " << make_tag( working_res ) << std::endl;


	std::string current_sequence = pose.sequence();
	std::cout << "DESIRED:" << desired_sequence << std::endl;
	std::cout << "CURRENT:" << current_sequence << std::endl;
	utility::vector1< Size > changed_pos;

	//Write over sequence?
	for (Size i = 1; i <= desired_sequence.size(); i++ ){

		char const new_seq = desired_sequence[i-1];
		if ( new_seq == current_sequence[i-1] ) continue;

		changed_pos.push_back( i );

		ResidueTypeCOP new_rsd_type( ResidueSelector().set_name1( new_seq ).exclude_variants().select( rsd_set )[1] );
		ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type, pose.residue( i ), pose.conformation() ) );
		Real const save_chi = pose.chi(i);
		pose.replace_residue( i, *new_rsd, false );
		pose.set_chi( i, save_chi );

	}


	std::cout << "Changed residues (without offset applied): " << make_tag_with_dashes( changed_pos ) << std::endl;

	///////////////////
	///////////////////
	return;
	///////////////////
	///////////////////

	// do the following later... this is copied from a protein modeling routine, and might be worth keeping in here?
	std::cout << "ABOUT TO DO LOOPS" << std::endl;
	//Now need to add in residues that don't exist already. This is tricky. Need to do it segment by segment.
	//Loop definition
	bool in_loop( false );
	Size start( 0 ), end( 0 );
	utility::vector1< std::pair< Size, Size > > loops;
	for (Size i = 1; i <= alignment_length; i++ ){
		Size const pdb_number( alignment2sequence[1][i] );
		if (pdb_number == 0 ) continue;
		//std::cout << "POS:" << i <<  " " << pdb_number << " " << sequence_mask(i) << " " << in_loop << std::endl;

		if ( !sequence_mask(i) && !in_loop) {
			in_loop = true;
			start = pdb_number;
		}
		if ( sequence_mask(i) && in_loop ) {
			in_loop = false;
			end = pdb_number-1;
			loops.push_back( std::make_pair( start, end ) );
		}
	}

	std::cout << "about to do termini" << std::endl;
	//1. Will need to handle special case for termini -- not coded yet!!
	//2. Probably will want to randomly move around cutpoint (by "prepending"!).
	for (Size n = 1; n <= loops.size(); n++ ) {

		std::cout << "LOOP " << n << " ==> " << loops[n].first << " " << loops[n].second << std::endl;

		for (Size i = loops[n].first; i <= loops[n].second; i++ ) {
			ResidueTypeCOP new_rsd_type( ResidueSelector().set_name1( desired_sequence[i-1] ).exclude_variants().select( rsd_set )[1] );
			ResidueOP new_rsd( ResidueFactory::create_residue( *new_rsd_type ) );
			pose.conformation().append_polymer_residue_after_seqpos( *new_rsd, i-1, true );
			pose.set_omega( i, 180.0 );
		}

	}

}

////////////////////////////////////////////////////////////////////////////
void
rna_thread_test(){

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::scoring;
	using namespace core::pose;

	////////////////////////////////////////////////
	//Read in fasta file.
	std::string const fasta_file ( option[ in::file::fasta ]()[1] );
	utility::vector1< std::string > sequences, pdb_names;
	Size const alignment_length = read_alignment_fasta_file( sequences, pdb_names, fasta_file );

	//Figure out a mask.
	ObjexxFCL::FArray1D_bool sequence_mask( alignment_length, false );
	setup_mask( sequence_mask, sequences );

	//Get from alignment to index inside each pdb.
	utility::vector1< std::map< Size, Size > > alignment2sequence;
	setup_alignment_map( alignment2sequence, sequences );

	////////////////////////////////////////////////
	//Read in native pdb
	Pose native_pose, template_pose, temp_pose;
	PoseOP pose_op( new Pose );
	Pose & pose( *pose_op );

	//std::string native_file = option[ in::file::native ];
	//	io::pdb::pose_from_pdb( native_pose, native_file );

	//Read in template pdb.
	std::string template_file = option[ in::file::s ][1];

	core::chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( "rna" );

	core::import_pose::pose_from_pdb( template_pose, *rsd_set, template_file );

	prepare_full_length_start_model( template_pose, pose, sequences, sequence_mask,  alignment2sequence, pdb_names, template_file );
	protocols::rna::virtualize_5prime_phosphates( pose );

	std::string outfile( "threaded.pdb" );
	if ( option[ out::file::o ].user() ) outfile = option[ out::file::o ]();

	std::cout << "Outputting: " << outfile << std::endl;
	pose.dump_pdb( outfile );

	// later put in erraser-style refinement of threaded mutated residues --
	// kind of an investment of time, though.

}


////////////////////////////////////////////////////////////////////////////////////////
void
initialize_sequence_mask( pose::Pose & pose, ObjexxFCL::FArray1D_bool & sequence_mask ) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	sequence_mask.dimension( pose.total_residue() );
	sequence_mask = true;

	if ( option[ sequence_mask_file ].user() ) {
		utility::io::izstream data_stream( option[ sequence_mask_file ] );
		std::string line;
		getline(data_stream, line);
		for (Size n = 1; n <= pose.total_residue(); n++ ) {
			sequence_mask( n ) = line[n-1];
		}
	}

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
	using namespace basic::options;

	//Uh, options?
	NEW_OPT( chi_stats, "", false );
	NEW_OPT( copy_native_chi, "", false );
	NEW_OPT( rhiju_fold_tree, "", false );
	NEW_OPT( use_native_CA, "", false );
	NEW_OPT( vary_geometry, "", false );
	NEW_OPT( mask_loop, "", false );
	NEW_OPT( cst_relax, "", false );
	NEW_OPT( sequence_mask_file, "", "sequence_mask.txt" );
	NEW_OPT( secstruct_file, "", "blah.secstruct" );
	NEW_OPT( cst_file, "", "cst.txt" );
	NEW_OPT( seq_offset, "", 0 );

	////////////////////////////////////////////////////////////////////////////
	// setup
	////////////////////////////////////////////////////////////////////////////
	devel::init(argc, argv);


	////////////////////////////////////////////////////////////////////////////
	// end of setup
	////////////////////////////////////////////////////////////////////////////

	protocols::viewer::viewer_main( my_main );

}
