// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>
#include <protocols/stepwise/modeler/util.hh> // for reroot
#include <protocols/stepwise/modeler/rna/util.hh> // for virtualize_free_rna_moieties
#include <protocols/stepwise/modeler/rna/checker/VDW_CachedRepScreenInfo.hh> // for fill_vdw_cached_rep_screen_info_from_command_line
#include <protocols/electron_density/SetupForDensityScoringMover.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/FullModelParameters.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/datacache/CacheableDataType.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/util.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/rna/util.hh>
#include <core/chemical/types.hh>
#include <core/chemical/VariantType.hh>
#include <core/import_pose/import_pose.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/sequence/Sequence.hh>
#include <core/sequence/util.hh>
#include <utility/stream_util.hh>
#include <utility/vector1.functions.hh>
#include <utility/tools/make_vector1.hh>

#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/magnesium.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/datacache/BasicDataCache.hh>
#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.setup.FullModelInfoSetupFromCommandLine" );

using namespace core;
using namespace core::pose;
using namespace core::pose::full_model_info;
using namespace core::id;
using utility::vector1;
using std::pair;
using std::make_pair;
using core::kinematics::FoldTree;

//////////////////////////////////////////////////////////////////////////////////////
//
// All setup functions, including PDB readin & FullModelInfo initialization for
//  stepwise application.
//
// Getting to be a big file, might be good to separate some functions out, and/or
//  turn into a class?
//
//                                -- rhiju, 2014
//
//////////////////////////////////////////////////////////////////////////////////////

namespace protocols {
namespace stepwise {
namespace setup {


//////////////////////////////////////////////////////////////////////////////////////
core::pose::PoseOP
get_pdb_with_full_model_info( std::string const & input_file,
	core::chemical::ResidueTypeSetCAP rsd_set ) {

	core::pose::PoseOP pose = get_pdb_and_cleanup( input_file, rsd_set );

	// a bit wasteful, since this checks command-line options over and over again, but hey this works.
	//  fill_full_model_info_from_command_line( *pose );

	return pose;
}

//////////////////////////////////////////////////////////////////////////////////////
// might be better to move these into core (e.g., core::pose::full_model_info ),
// or into a new protocols/full_model_setup/ directory.
core::pose::PoseOP
get_pdb_and_cleanup( std::string const & input_file,
	core::chemical::ResidueTypeSetCAP rsd_set )
{
	using namespace core::pose;
	PoseOP input_pose( new Pose );
	core::chemical::ResidueTypeSetCOP rsd_set_op( rsd_set );
	import_pose::pose_from_file( *input_pose, *rsd_set_op, input_file , core::import_pose::PDB_file);
	tag_into_pose( *input_pose, input_file );
	cleanup( *input_pose );
	make_sure_full_model_info_is_setup( *input_pose );
	return input_pose;
}


//////////////////////////////////////////////////////////////////////////////////////
// currently have stuff we need for RNA... put any protein cleanup here too.
void
cleanup( pose::Pose & pose, bool const force_cut_at_rna_chainbreak /* = false */  ) {
	rna::figure_out_reasonable_rna_fold_tree( pose, force_cut_at_rna_chainbreak);
	rna::virtualize_5prime_phosphates( pose );
	pose.conformation().detect_disulfides();
}

///////////////////////////////////////////////////////////////
void
initialize_native_and_align_pose( PoseOP & native_pose,
	PoseOP & align_pose,
	core::chemical::ResidueTypeSetCAP rsd_set,
	PoseCOP start_pose ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	if ( option[ in::file::native ].user() )  {
		align_pose = native_pose = get_pdb_with_full_model_info( option[ in::file::native ](), rsd_set );
	}
	if ( option[ OptionKeys::stepwise::align_pdb ].user() ) {
		align_pose = get_pdb_with_full_model_info(  option[ OptionKeys::stepwise::align_pdb ](), rsd_set );
	}
	if ( align_pose == 0 && option[ in::file::s ].user() ) {
		align_pose = start_pose->clone();
	}
	if ( option[ OptionKeys::stepwise::virtualize_free_moieties_in_native ]() ) { // could generalize to proteins
		if ( native_pose != 0 )  modeler::rna::virtualize_free_rna_moieties( *native_pose );
		if ( align_pose  != 0 ) modeler::rna::virtualize_free_rna_moieties( *align_pose );
	}

	// really? -- this results in rms and rms_fill calculation even without native -- confusing. rhiju, dec. 2016
	//	if ( native_pose == 0 && align_pose != 0 ) native_pose = align_pose;
}


///////////////////////////////////////////////////////////////////////////////////////
pose::PoseOP
initialize_pose_and_other_poses_from_command_line( core::chemical::ResidueTypeSetCAP rsd_set ) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::io::silent;

	utility::vector1< std::string > const & input_pdb_files    = option[ in::file::s ]();
	utility::vector1< std::string > const & input_silent_files = option[ in::file::silent ]();

	utility::vector1< pose::PoseOP > input_poses;
	for ( Size n = 1; n <= input_pdb_files.size(); n++ ) {
		input_poses.push_back( get_pdb_and_cleanup( input_pdb_files[ n ], rsd_set ) );
	}
	for ( Size n = 1; n <= input_silent_files.size(); n++ ) {
		PoseOP pose( new Pose );
		SilentFileData silent_file_data;
		core::chemical::ResidueTypeSetCOP rsd_set_op( rsd_set );
		silent_file_data.read_file( input_silent_files[n] );
		silent_file_data.begin()->fill_pose( *pose, *rsd_set_op );
		input_poses.push_back( pose );
	}

	pair< vector1< Size >, vector1< char > > const & input_resnum_and_chain = option[ in::file::input_res ].resnum_and_chain();
	vector1< Size > const & input_res_list = input_resnum_and_chain.first;
	if ( input_res_list.size() ) {
		vector1< char > input_chain_list = input_resnum_and_chain.second;
		Size input_res_count = 0;
		for ( Size n = 1; n <= input_poses.size(); n++ ) {
			Pose & pose = *input_poses[ n ];
			PDBInfoOP pdb_info( new PDBInfo( pose ) );
			vector1< Size > input_res_for_pose;
			vector1< char > input_chain_for_pose;
			for ( Size k = 1; k <= pose.size(); k++ ) {
				input_res_count++;
				runtime_assert( input_res_count <= input_res_list.size() );
				Size const & number_in_full_model = input_res_list[ input_res_count ];
				input_res_for_pose.push_back( number_in_full_model );
				input_chain_for_pose.push_back( input_chain_list[ input_res_count ] );
			}
			pdb_info->set_numbering( input_res_for_pose );
			pdb_info->set_chains(    input_chain_for_pose );
			pose.pdb_info( pdb_info );
		}
		runtime_assert( input_res_count == input_res_list.size() );
	}

	if ( input_poses.size() == 0 ) input_poses.push_back( core::pose::PoseOP( new Pose ) ); // just a blank pose for now.

	if ( option[ full_model::other_poses ].user() ) {
		get_other_poses( input_poses, option[ full_model::other_poses ](), rsd_set );
	}

	fill_full_model_info_from_command_line( input_poses );  //FullModelInfo (minimal object needed for add/delete)
	protocols::stepwise::modeler::rna::checker::fill_vdw_cached_rep_screen_info_from_command_line( *input_poses[1] );
	return input_poses[1];
}

///////////////////////////////////////////////////////////////////////////////////////
void
get_other_poses( utility::vector1< pose::PoseOP > & other_poses,
	utility::vector1< std::string > const & other_files,
	core::chemical::ResidueTypeSetCAP rsd_set ) {

	for ( Size n = 1; n <= other_files.size(); n++ ) {
		other_poses.push_back( get_pdb_and_cleanup( other_files[ n ], rsd_set ) );
	}
}

///////////////////////////////////////////////////////////////////////////////////////
void
fill_full_model_info_from_command_line( pose::Pose & pose ) {
	vector1< PoseOP > other_pose_ops; // dummy
	fill_full_model_info_from_command_line( pose, other_pose_ops );
}

///////////////////////////////////////////////////////////////////////////////////////
void
fill_full_model_info_from_command_line( vector1< PoseOP > & pose_ops ) {

	runtime_assert( pose_ops.size() > 0 );

	vector1< PoseOP > other_pose_ops;
	for ( Size n = 2; n <= pose_ops.size(); n++ ) other_pose_ops.push_back( pose_ops[n] );

	fill_full_model_info_from_command_line( *(pose_ops[1]), other_pose_ops );
}


///////////////////////////////////////////////////////////////////////////////////////
// this is needed so that first pose holds pointers to other poses.
void
fill_full_model_info_from_command_line( pose::Pose & pose, vector1< PoseOP > & other_pose_ops ) {

	vector1< Pose * > pose_pointers;
	pose_pointers.push_back( & pose );
	for ( Size n = 1; n <= other_pose_ops.size(); n++ ) pose_pointers.push_back( & (*other_pose_ops[ n ]) );

	fill_full_model_info_from_command_line( pose_pointers );

	nonconst_full_model_info( pose ).set_other_pose_list( other_pose_ops );

}


///////////////////////////////////////////////////////////////////////////////////////
vector1< Size >
get_cutpoints_from_numbering( vector1< core::sequence::SequenceCOP > const & fasta_sequences,
	vector1< char > const & conventional_chains,
	vector1< int  > const & conventional_numbering ) {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	vector1< Size > cutpoints;
	Size ntot( 0 );
	// explicit chain boundaries
	for ( Size n = 1; n <= fasta_sequences.size(); n++ ) {
		std::string const sequence = core::pose::rna::remove_bracketed( fasta_sequences[ n ]->sequence() );
		vector1< Size > const & spacer_pos = fasta_sequences[ n ]->spacer_positions();
		for ( Size q = 1; q <= spacer_pos.size(); q++ ) cutpoints.push_back( ntot + spacer_pos[ q ] );
		ntot += sequence.size();
		if ( n != fasta_sequences.size() /*very end is not 'cutpoint'*/ ) cutpoints.push_back( ntot );
	}
	runtime_assert( ntot == conventional_chains.size() );
	runtime_assert( ntot == conventional_numbering.size() );
	for ( Size k = 1; k < ntot; k++ ) {
		if ( cutpoints.has_value( k ) ) continue;
		if ( conventional_chains[ k ] != conventional_chains[ k+1 ] ) {
			cutpoints.push_back( k );
		} else { // break within chain.
			if ( conventional_numbering[ k+1 ] > conventional_numbering[ k ] + 1  ) {
				if ( !option[ full_model::allow_jump_in_numbering ]() ) {
					TR << TR.Red << "There appears to be a break in numbering at " << conventional_numbering[ k ] << " so adding to cutpoint list. If that is wrong, use flag -allow_jump_in_numbering." << TR.Reset << std::endl;
					cutpoints.push_back( k );
				}
			}
		}
	}

	std::sort( cutpoints.begin(), cutpoints.end() );
	return cutpoints;
}


/////////////////////////////////////////////////////////////////////////////////////////
// Kind of a hack. Would be better to somehow figure out if residue is not polymeric --
//  perhaps try to instantiate with ResidueTypeSet?
void
get_extra_cutpoints_from_names( Size const nres,
	vector1< Size > & cutpoints,
	std::map< Size, std::string > const & non_standard_residue_map )
{
	for ( Size n = 1; n <= nres; n++ ) {
		if ( cutpoints.has_value( n ) ) continue;
		if ( !modeler::stepwise_addable_residue( n, non_standard_residue_map ) ) {
			cutpoints.push_back( n );
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////
vector1< Size >
get_cutpoints( vector1< core::sequence::SequenceCOP > const & fasta_sequences,
	std::map< Size, std::string > const & non_standard_residue_map,
	vector1< char > const & conventional_chains,
	vector1< int  > const & conventional_numbering ) {
	vector1< Size > cutpoints = get_cutpoints_from_numbering( fasta_sequences, conventional_chains, conventional_numbering );
	get_extra_cutpoints_from_names( conventional_numbering.size(), cutpoints, non_standard_residue_map );
	return cutpoints;
}

///////////////////////////////////////////////////////////////////////////////////////
// looks for tab-delimited tags like 'chain:A' and 'res_num:5-20' in fasta IDs.
// move this to util.cc in core?
///////////////////////////////////////////////////////////////////////////////////////
void
get_conventional_chains_and_numbering( vector1< core::sequence::SequenceCOP > const & fasta_sequences,
	vector1< char > & conventional_chains,
	vector1< int  > & conventional_numbering ) {
	using utility::string_split;
	bool found_info_in_previous_sequence( false );
	Size count( 0 );
	for ( Size n = 1; n <= fasta_sequences.size(); n++ ) {
		utility::vector1< char > chains;
		vector1< int > resnum;
		bool found_info( false );
		std::string tag;
		std::stringstream ss( fasta_sequences[n]->id() );
		while ( ss.good() ) {
			ss >> tag;
			bool string_is_ok( false );
			pair< std::vector< int >, std::vector< char > > resnum_and_chain = utility::get_resnum_and_chain( tag, string_is_ok );
			if ( !string_is_ok ) continue;
			for ( Size n = 0; n < resnum_and_chain.first.size(); n++ ) {
				resnum.push_back( resnum_and_chain.first[n] );
				chains.push_back( resnum_and_chain.second[n] );
			}
			found_info = true;
		}
		if ( n > 1 ) runtime_assert( found_info == found_info_in_previous_sequence );

		Size const clean_len = core::pose::rna::remove_bracketed( fasta_sequences[n]->sequence() ).size();
		if ( !found_info || resnum.size() != clean_len ) { //fasta_sequences[n]->sequence().size() /*happens with stray numbers*/ ) {
			resnum.clear();
			for ( Size q = 1; q <= clean_len; ++q ) { //fasta_sequences[n]->sequence().size(); q++ ) {
				resnum.push_back( ++count );
				chains.push_back( ' ' ); // unknown chain
			}
		}
		std::string const sequence = fasta_sequences[n]->sequence();
		runtime_assert( clean_len == resnum.size() ); //sequence.size() == resnum.size() );
		for ( Size q = 1; q <= clean_len /*sequence.size()*/; q++ ) conventional_chains.push_back( chains[ q ] );
		for ( Size q = 1; q <= clean_len /*sequence.size()*/; q++ ) conventional_numbering.push_back( resnum[q] );

		found_info_in_previous_sequence  = found_info;
	}
}

///////////////////////////////////////////////////////////////////////////////////////
// Need a 'water bank' to hold up to 6 waters for each Mg(2+) during stepwise modeling
//  with water hydration. Some of these waters will get virtualized, but its important
//  that number of residues stays constant.
//
// Following is based on an assumption that the *only* waters to be modeled in stepwise
//  would be ones associated with Mg(2+)... obviously that is not very general.
//
// It might be smarter to look inside the PDBs (esp. any native PDBs) to count
//   waters associated with each Mg(2+), one by one, and supplement water bank to
//   ensure 6 waters per Mg(2+).
void
setup_water_bank_for_magnesiums( std::map< Size, std::string > & non_standard_residue_map,
	vector1< core::sequence::SequenceOP > & fasta_sequences ) {
	using namespace core::sequence;

	// how many magnesiums are there? how many waters are there?
	Size offset( 0 ), num_magnesiums( 0 ), num_waters( 0 );
	for ( Size n = 1; n <= fasta_sequences.size(); n++ ) {
		std::string const sequence = fasta_sequences[n]->sequence();
		for ( Size i = 1; i <= sequence.size(); i++ ) {
			if ( sequence[ i - 1 ] == 'Z' ) {
				std::map< Size, std::string >::const_iterator it = non_standard_residue_map.find( offset+i );
				if ( it != non_standard_residue_map.end() && it->second == "MG" )  num_magnesiums++;
			} else if ( sequence[ i - 1 ] == 'w' ) {
				std::map< Size, std::string >::const_iterator it = non_standard_residue_map.find( offset+i );
				if ( it != non_standard_residue_map.end() && it->second == "HOH" )  num_waters++;
			}
		}
		offset += sequence.size();
	}

	runtime_assert( num_waters <= 6 * num_magnesiums ); // currently do not model any non-Mg(2+) waters.

	// are there at least six waters for each magnesium? if not, let's put them in.
	// Use chain 'w' semi-arbitrarily, and numbering from 1001 onward. -- this could be a disaster.
	std::string sequence;
	std::stringstream id;
	Size extra_water_number( 1000 );
	for ( Size n = num_waters; n <= 6 * num_magnesiums; n++ ) {
		sequence += 'w';
		offset++; // number of residues.
		extra_water_number++; // 1001, 1002, ...
		non_standard_residue_map[ offset ] = "HOH";
		id << " w:" << extra_water_number;
	}

	fasta_sequences.push_back( SequenceOP( new Sequence( sequence, id.str() ) ) );
}

void
look_for_dna( vector1< core::sequence::SequenceOP > & fasta_sequences )
{
	using namespace core::sequence;
	std::map<char, std::string> rna_name = { {'a',"RAD"}, {'c',"RCY"}, {'g',"RGU"}, {'t',"5MU"} };
	for ( Size n = 1; n <= fasta_sequences.size(); n++ ) {
		auto fasta_sequence = fasta_sequences[ n ];
		auto id = fasta_sequences[ n ]->id();
		std::stringstream ss( fasta_sequences[n]->id() );
		bool is_DNA( false );
		while ( ss.good() ) {
			std::string tag;
			ss >> tag;
			if ( tag == "DNA" ) {
				is_DNA = true; break;
			}
		}
		auto sequence = fasta_sequences[ n ]->sequence();
		utility::vector1< std::string > fullname_list; // a vector of non-standard full names
		std::vector< Size > oneletter_to_fullname_index; // for each one-letter sequence, zero means no fullname given
		std::string one_letter_sequence;
		parse_sequence( sequence, fullname_list, oneletter_to_fullname_index, one_letter_sequence );
		if ( !is_DNA ) {
			for ( Size k = 1; k <= one_letter_sequence.size(); k++ ) {
				if ( one_letter_sequence[ k-1 ] == 't' ) {
					Size const pos = oneletter_to_fullname_index[ k - 1 ];
					if ( pos == 0 || fullname_list[ pos ].substr( 3, 14 ) != ":deoxy_O2prime" ) {
						utility_exit_with_message( "Seeing a 't' in sequence -- if you want this to be DNA, include the tag 'DNA' in the FASTA id for this sequence. if you want this to be RNA, use X[5MU] as the character." );
					}
				}
			}
		}
		if ( is_DNA ) {
			auto fasta_sequence_new = fasta_sequence;
			vector1< std::string > fullname_list_new;
			for ( Size k = 1; k <= oneletter_to_fullname_index.size(); k++ ) {
				Size const pos = oneletter_to_fullname_index[ k - 1 ];
				std::string fullname_new;
				if ( pos > 0 ) {
					fullname_new = fullname_list[ pos ];
				} else {
					char q = one_letter_sequence[ k - 1 ];
					if ( !rna_name.count( q ) ) utility_exit_with_message( "The character "+std::string(&q)+" is not a recognized DNA nucleotide. Use a,c,g,t." );
					fullname_new = rna_name[ q ];
				}
				if ( fullname_new.substr(3,14) != ":deoxy_O2prime" ) {
					fullname_new = fullname_new.substr(0,3) + ":deoxy_O2prime" + fullname_new.substr(3,std::string::npos);
				}
				fullname_list_new.push_back( fullname_new );
			}
			std::string sequence_new;
			for ( Size k = 1; k <= one_letter_sequence.size(); k++ ) {
				sequence_new.push_back( one_letter_sequence[ k-1 ] );
				sequence_new.push_back( '[' );
				sequence_new.append( fullname_list_new[ k ] );
				sequence_new.push_back( ']' );
			}
			fasta_sequences[ n ] = SequenceOP( new Sequence( sequence_new, id ) );
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////
FullModelParametersOP
get_sequence_information(
	std::string const & fasta_file,
	vector1< Size > & cutpoint_open_in_full_model )
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	vector1< core::sequence::SequenceOP > fasta_sequences = core::sequence::read_fasta_file( fasta_file );

	// calebgeniesse: setup for edensity scoring, if map is provided via cmd-line
	if ( option[ edensity::mapfile ].user() ) {
		// update fasta_sequences accordingly
		Size idx = fasta_sequences.size();
		fasta_sequences[idx]->append_char('X');
		fasta_sequences[idx]->id( fasta_sequences[idx]->id() + " z:1" );
	}

	look_for_dna( fasta_sequences );
	std::map< Size, std::string > non_standard_residue_map  = core::sequence::parse_out_non_standard_residues( fasta_sequences /*will reduce to one-letter*/ );
	if ( option[ magnesium::hydrate ]() ) setup_water_bank_for_magnesiums( non_standard_residue_map, fasta_sequences );
	std::string const desired_sequence           = core::sequence::get_concatenated_sequence( fasta_sequences );

	FullModelParametersOP full_model_parameters( new FullModelParameters( desired_sequence ) );
	vector1< char > conventional_chains;
	vector1< int  > conventional_numbering;
	get_conventional_chains_and_numbering( fasta_sequences, conventional_chains, conventional_numbering );

	full_model_parameters->set_conventional_numbering( conventional_numbering );
	full_model_parameters->set_conventional_chains( conventional_chains );
	full_model_parameters->set_non_standard_residue_map( non_standard_residue_map );
	cutpoint_open_in_full_model  = get_cutpoints( fasta_sequences, non_standard_residue_map,
		conventional_chains, conventional_numbering );
	return full_model_parameters;
}


///////////////////////////////////////////////////////////////////////////////////////
void
fill_full_model_info_from_command_line( vector1< Pose * > & pose_pointers ) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	if ( !option[ in::file::fasta ].user() ) {
		for ( Size n = 1; n <= pose_pointers.size(); n++ ) make_sure_full_model_info_is_setup( *pose_pointers[n] );
		return;
	}

	std::string const fasta_file = option[ in::file::fasta ]()[1];
	vector1< Size > cutpoint_open_in_full_model = option[ full_model::cutpoint_open ]();
	FullModelParametersOP full_model_parameters = get_sequence_information( fasta_file, cutpoint_open_in_full_model );
	std::string const & desired_sequence = full_model_parameters->full_sequence();

	// calebgeniesse: setup for edensity scoring, if map is provided via cmd-line
	if ( option[ edensity::mapfile ].user() ) {
		// update pose
		if ( pose_pointers[1]->total_residue() > 0 ) {
			setup_for_density_scoring( *pose_pointers[1] );
		}
	}

	if ( option[ full_model::cutpoint_open ].user() ) {
		cutpoint_open_in_full_model =
			full_model_parameters->conventional_to_full( option[ full_model::cutpoint_open ].resnum_and_chain() );
	}

	vector1< Size > extra_minimize_res =
		full_model_parameters->conventional_to_full( option[ full_model::extra_min_res ].resnum_and_chain() );
	vector1< Size > sample_res         =
		full_model_parameters->conventional_to_full( option[ full_model::sample_res ].resnum_and_chain() ); //stuff that can be resampled.
	vector1< Size > working_res        =
		full_model_parameters->conventional_to_full( option[ full_model::working_res ].resnum_and_chain() ); //all working stuff
	vector1< Size > terminal_res       =
		full_model_parameters->conventional_to_full( option[ full_model::rna::terminal_res ].resnum_and_chain() );
	vector1< Size > block_stack_above_res  =
		full_model_parameters->conventional_to_full( option[ full_model::rna::block_stack_above_res ].resnum_and_chain() );
	vector1< Size > block_stack_below_res  =
		full_model_parameters->conventional_to_full( option[ full_model::rna::block_stack_below_res ].resnum_and_chain() );
	vector1< Size > preferred_root_res =
		full_model_parameters->conventional_to_full( option[ full_model::root_res ].resnum_and_chain() );
	vector1< Size > jump_res           =
		full_model_parameters->conventional_to_full( option[ full_model::jump_res ].resnum_and_chain() );
	vector1< Size > const cutpoint_closed          =
		full_model_parameters->conventional_to_full( option[ full_model::cutpoint_closed ].resnum_and_chain() );
	if ( option[ full_model::cyclize ].user() ) utility_exit_with_message( "Cannot handle cyclize yet in stepwise." );
	vector1< Size > const bulge_res                =
		full_model_parameters->conventional_to_full( option[ full_model::rna::bulge_res ].resnum_and_chain() );
	vector1< Size > const extra_minimize_jump_res  =
		full_model_parameters->conventional_to_full( option[ full_model::extra_min_jump_res ].resnum_and_chain() );
	vector1< Size > const virtual_sugar_res        =
		full_model_parameters->conventional_to_full( option[ full_model::virtual_sugar_res ].resnum_and_chain() );
	update_jump_res( jump_res, extra_minimize_jump_res );

	// calebgeniesse: override preferred_root_res if mapfile provided via cmd-line (hacky)
	if ( option[ edensity::mapfile ].user() ) {
		preferred_root_res = full_model_parameters->conventional_to_full(
			std::make_pair(
			utility::tools::make_vector1< int >(1),
			utility::tools::make_vector1< char >('z')
			)
		);
	}

	// Figure out res_list and input_domain_map.
	vector1< vector1< Size > > pose_res_lists;
	std::string const clean_desired_seq = core::pose::rna::remove_bracketed( desired_sequence );
	Size const desired_nres = clean_desired_seq.size();
	vector1< Size > input_domain_map( desired_nres, 0 );
	for ( Size n = 1; n <= pose_pointers.size(); n++ ) {
		Pose & pose = *pose_pointers[n];
		vector1< Size > res_list = full_model_parameters->conventional_to_full( make_pair( get_res_num_from_pdb_info( pose ), get_chains_from_pdb_info( pose ) ) );
		reorder_pose( pose, res_list );
		pose_res_lists.push_back( res_list );
		for ( Size k = 1; k <= pose.size(); k++ ) {
			if ( !sample_res.has_value( res_list[ k ] ) ) input_domain_map[ res_list[ k ] ] = n;
		}
	}

	// calebgeniesse: override input_domain_map if mapfile provided via cmd-line (hacky)
	if ( desired_sequence[desired_sequence.size()-1] == 'X' ) {
		input_domain_map[desired_sequence.size()] = 1;// input_domain_map.size() ? max(input_domain_map) + 1 : 1; // or just 1?
	}

	if ( option[ full_model::motif_mode ]() ) figure_out_motif_mode( extra_minimize_res, terminal_res, working_res, input_domain_map, cutpoint_open_in_full_model );
	add_block_stack_variants( pose_pointers, pose_res_lists, block_stack_above_res, block_stack_below_res );

	// calebgeniesse: hacky way of making sure virt root is not included in terminal_res/extra_minimize_res
	if ( desired_sequence[desired_sequence.size()-1] == 'X' ) {
		if ( terminal_res.has_value( desired_sequence.size() ) ) {
			terminal_res.erase( std::remove(terminal_res.begin(), terminal_res.end(), desired_sequence.size()), terminal_res.end());
		}
		if ( extra_minimize_res.has_value( desired_sequence.size() ) ) {
			extra_minimize_res.erase( std::remove(extra_minimize_res.begin(), extra_minimize_res.end(), desired_sequence.size()), extra_minimize_res.end());
		}
	}

	// everything that is not fixed is sampleable (unless -sample_res explicitly specified).
	if ( sample_res.size() == 0 )  sample_res  = figure_out_sample_res( input_domain_map, working_res );
	if ( working_res.size() == 0 ) working_res = figure_out_working_res( input_domain_map, sample_res );
	vector1< Size > fixed_domain_map = figure_out_fixed_domain_map( input_domain_map, extra_minimize_res ); //remove extra_minimize_res.

	setup_fold_trees( pose_pointers, cutpoint_open_in_full_model /* can update */, fixed_domain_map /* can update */,
		cutpoint_closed, extra_minimize_res, extra_minimize_jump_res,
		sample_res, working_res, jump_res,
		preferred_root_res, virtual_sugar_res,
		*full_model_parameters, pose_res_lists );

	vector1< Size > const dock_domain_map = figure_out_dock_domain_map( cutpoint_open_in_full_model /* can be updated */,
		pose_res_lists, working_res, sample_res, desired_nres );

	// some checks
	check_working_res( working_res, input_domain_map, sample_res );
	check_extra_minimize_res_are_input( extra_minimize_res,      input_domain_map );
	check_extra_minimize_res_are_input( extra_minimize_jump_res, input_domain_map );

	full_model_parameters->set_parameter( INPUT_DOMAIN,  input_domain_map );
	full_model_parameters->set_parameter( FIXED_DOMAIN,  fixed_domain_map );
	full_model_parameters->set_parameter( DOCK_DOMAIN,   dock_domain_map  );
	full_model_parameters->set_parameter_as_res_list( CUTPOINT_OPEN, cutpoint_open_in_full_model );
	full_model_parameters->set_parameter_as_res_list( EXTRA_MINIMIZE, extra_minimize_res );
	full_model_parameters->set_parameter_as_res_list( SAMPLE, sample_res );
	full_model_parameters->set_parameter_as_res_list( WORKING, working_res );
	full_model_parameters->set_parameter_as_res_list( CALC_RMS,
		full_model_parameters->conventional_to_full( option[ full_model::calc_rms_res ].resnum_and_chain() ) );
	full_model_parameters->set_parameter_as_res_list( PREFERRED_ROOT, preferred_root_res );
	full_model_parameters->set_parameter_as_res_list( RNA_SYN_CHI,
		full_model_parameters->conventional_to_full( option[ full_model::rna::force_syn_chi_res_list ].resnum_and_chain() ) );
	full_model_parameters->set_parameter_as_res_list( RNA_ANTI_CHI,
		full_model_parameters->conventional_to_full( option[ full_model::rna::force_anti_chi_res_list ].resnum_and_chain() ) );
	full_model_parameters->set_parameter_as_res_list( RNA_NORTH_SUGAR,
		full_model_parameters->conventional_to_full( option[ full_model::rna::force_north_sugar_list ].resnum_and_chain() ) );
	full_model_parameters->set_parameter_as_res_list( RNA_SOUTH_SUGAR,
		full_model_parameters->conventional_to_full( option[ full_model::rna::force_south_sugar_list ].resnum_and_chain() ) );
	full_model_parameters->set_parameter_as_res_list( RNA_TERMINAL, terminal_res );
	full_model_parameters->set_parameter_as_res_list( RNA_BLOCK_STACK_ABOVE, block_stack_above_res );
	full_model_parameters->set_parameter_as_res_list( RNA_BLOCK_STACK_BELOW, block_stack_below_res );
	full_model_parameters->set_parameter_as_res_list( RNA_BULGE,  bulge_res );
	full_model_parameters->set_parameter_as_res_list( RNA_SAMPLE_SUGAR,
		full_model_parameters->conventional_to_full( option[ full_model::rna::sample_sugar_res ].resnum_and_chain() ) );
	full_model_parameters->set_parameter_as_res_list_in_pairs( full_model_info::JUMP, jump_res );
	full_model_parameters->set_parameter_as_res_list_in_pairs( EXTRA_MINIMIZE_JUMP, extra_minimize_jump_res );
	full_model_parameters->read_disulfides( option[ OptionKeys::stepwise::protein::disulfide_file]() );
	if ( option[ constraints::cst_file ].user() ) full_model_parameters->read_cst_file( option[ constraints::cst_file ]()[ 1 ] );

	// move this code block somewhere else when ready.
	if ( option[ basic::options::OptionKeys::stepwise::monte_carlo::vary_loop_length_frequency ] > 0.0 ) {
		// placeholder -- testing if loops can be 'evaporated'.
		vector1< Size > full_model_res_no_loops;
		for ( Size n = 1; n <= desired_nres; ++n ) {
			if ( clean_desired_seq[ n-1 ] != 'n' || input_domain_map[ n ] ) full_model_res_no_loops.push_back( n );
		}
		full_model_parameters = full_model_parameters->slice( full_model_res_no_loops );
		for ( Size n = 1; n <= pose_res_lists.size(); n++ ) {
			vector1< Size > const & pose_res_list = pose_res_lists[ n ];
			vector1< Size > pose_res_list_new;
			for ( Size i = 1; i <= pose_res_list.size(); i++ ) {
				if ( full_model_res_no_loops.has_value( pose_res_list[ i ] ) ) {
					pose_res_list_new.push_back( full_model_res_no_loops.index( pose_res_list[ i ] ) );
				}
			}
			pose_res_lists[ n ] = pose_res_list_new;
		}
	}

	for ( Size n = 1; n <= pose_pointers.size(); n++ ) {
		Pose & pose = *pose_pointers[n];
		FullModelInfoOP full_model_info_for_pose( new FullModelInfo( full_model_parameters ) );
		full_model_info_for_pose->set_res_list( pose_res_lists[ n ] );
		full_model_info_for_pose->update_submotif_info_list();
		set_full_model_info( pose, full_model_info_for_pose );
		update_pose_objects_from_full_model_info( pose ); // for output pdb or silent file (residue numbering), constraints, disulfides
		modeler::fix_up_residue_type_variants( pose ); // for sample sugars...
	}
}


////////////////////////////////////////////////////////////////////////////////////
// In principle, all the parameters being sent in to this function could be
// encapsulated into a FullModelParameters object. However, its important to be conscious
// of what is const and what can be updated at the level of the individual variables.
void
setup_fold_trees( vector1< Pose * > & pose_pointers,
	vector1< Size > & cutpoint_open_in_full_model /* can be updated here*/,
	vector1< Size > & fixed_domain_map /* domain colors can be updated here */,
	vector1< Size > const & cutpoint_closed,
	vector1< Size > const & extra_minimize_res,
	vector1< Size > const & extra_minimize_jump_res,
	vector1< Size > const & sample_res,
	vector1< Size > const & working_res,
	vector1< Size > const & jump_res,
	vector1< Size > const & preferred_root_res,
	vector1< Size > const & virtual_sugar_res,
	FullModelParameters const & full_model_parameters,
	vector1< vector1< Size > > const & pose_res_lists ) {

	for ( Size n = 1; n <= pose_pointers.size(); n++ ) {
		Pose & pose = *pose_pointers[n];
		vector1< Size > const & res_list = pose_res_lists[ n ];
		for ( Size i = 1; i < pose.size(); i++ ) {
			if ( (res_list[ i+1 ] > res_list[ i ] + 1) && !pose.fold_tree().is_cutpoint(i) ) {
				put_in_cutpoint( pose, i );
			}
			if ( cutpoint_open_in_full_model.has_value( res_list[ i ]) ) continue;
			if ( (res_list[ i+1 ] == res_list[ i ] + 1) &&
					pose.fold_tree().is_cutpoint( i ) &&
					!pose.residue_type( i   ).has_variant_type( core::chemical::CUTPOINT_LOWER ) &&
					!pose.residue_type( i+1 ).has_variant_type( core::chemical::CUTPOINT_UPPER ) ) {
				TR << TR.Red << "There appears to be a strand boundary at " << res_list[ i ] << " so adding to cutpoint_in_full_model." << TR.Reset << std::endl;
				cutpoint_open_in_full_model.push_back( res_list[ i ] ); continue;
			}
			if ( pose.residue_type( i ).is_RNA() != pose.residue_type( i+1 ).is_RNA() ) {
				cutpoint_open_in_full_model.push_back( res_list[ i ] ); continue;
			}
			if ( pose.residue_type( i ).is_protein() != pose.residue_type( i+1 ).is_protein() ) {
				cutpoint_open_in_full_model.push_back( res_list[ i ] ); continue;
			}
		}

		add_cutpoint_closed( pose, res_list, cutpoint_closed );

		update_pose_fold_tree( pose,
			res_list,
			extra_minimize_res, sample_res, jump_res,
			full_model_parameters /* set preferences for chain connections*/ );
		add_virtual_sugar_res( pose, res_list,
			virtual_sugar_res ); // this was for checks -- no longer in use?

		update_fixed_domain_from_extra_minimize_jump_res( fixed_domain_map, pose, res_list, extra_minimize_jump_res );

		vector1< Size> root_partition_res; for ( Size n = 1; n <= pose.size(); n++ ) root_partition_res.push_back( n );
		modeler::reroot( pose, root_partition_res, res_list, preferred_root_res, fixed_domain_map,
			cutpoint_open_in_full_model, working_res );

	}
}

////////////////////////////////////////////////////////////////////////////////////
void
update_pose_fold_tree( pose::Pose & pose,
	vector1< Size > const & res_list,
	vector1< Size > const & extra_min_res,
	vector1< Size > const & sample_res,
	vector1< Size > const & jump_res,
	core::pose::full_model_info::FullModelParameters const & full_model_parameters ) {

	if ( pose.size() == 0 ) return;

	vector1< vector1< Size > > all_res_in_chain, all_fixed_res_in_chain;
	vector1< Size > moveable_res = sample_res;
	for ( Size n = 1; n <= extra_min_res.size(); n++ ) moveable_res.push_back( extra_min_res[n] );
	define_chains( pose, all_res_in_chain, all_fixed_res_in_chain, res_list, moveable_res );
	Size nchains = all_res_in_chain.size();

	vector1< Size > jump_partners1, jump_partners2, cuts, blank_vector;
	vector1< pair< Size, Size > > chain_connections;
	setup_user_defined_jumps( jump_res, jump_partners1, jump_partners2,
		chain_connections, res_list, all_res_in_chain );
	runtime_assert( jump_partners1.size() < nchains );

	// needed to determine preferences for chain connections
	pair< vector1< int >, vector1< char > > const resnum_and_chain_in_pose = full_model_parameters.full_to_conventional_resnum_and_chain( res_list );

	// choose fixed_res as jump partners first.
	setup_jumps( pose, jump_partners1, jump_partners2, chain_connections, all_fixed_res_in_chain, resnum_and_chain_in_pose );
	// use moving res if necessary
	setup_jumps( pose, jump_partners1, jump_partners2, chain_connections, all_res_in_chain, resnum_and_chain_in_pose );
	runtime_assert( jump_partners1.size() == (nchains - 1) );

	for ( Size n = 1; n < nchains; n++ ) cuts.push_back( all_res_in_chain[n][ all_res_in_chain[n].size() ] );
	FoldTree f = get_tree( pose, cuts, jump_partners1, jump_partners2 );
	pose.fold_tree( f );
}

///////////////////////////////////////////////////////////////////////////////////////
void
define_chains( pose::Pose const & pose,
	vector1< vector1< Size > > & all_res_in_chain,
	vector1< vector1< Size > > & all_fixed_res_in_chain,
	vector1< Size > const & res_list,
	vector1< Size > const & moveable_res ) {

	Size chain_start( 1 ), chain_end( 0 );
	for ( Size n = 1; n <= pose.size(); n++ ) {
		if ( !pose.fold_tree().is_cutpoint( n ) &&
				n < pose.size() ) continue;
		chain_end = n;
		vector1< Size > res_in_chain, fixed_res_in_chain;
		for ( Size i = chain_start; i <= chain_end; i++ ) {
			res_in_chain.push_back( i );
			if ( !moveable_res.has_value( res_list[ i ] ) ) fixed_res_in_chain.push_back( i );
		}
		all_res_in_chain.push_back( res_in_chain );
		all_fixed_res_in_chain.push_back( fixed_res_in_chain );
		chain_start = chain_end + 1;
	}
}


////////////////////////////////////////////////////////////////////////////////
void
setup_user_defined_jumps( vector1< Size > const & jump_res,
	vector1< Size > & jump_partners1,
	vector1< Size > & jump_partners2,
	vector1< pair< Size, Size > > & chain_connections,
	vector1< Size > const & res_list,
	vector1< vector1< Size > > const & all_res_in_chain ) {

	vector1< Size > connection_domains = get_connection_domains( chain_connections, all_res_in_chain.size() );
	// Figure out jump res.
	for ( Size n = 1; n <= jump_res.size()/2; n++ ) {
		if ( res_list.has_value( jump_res[ 2*n - 1 ] )  &&
				res_list.has_value( jump_res[ 2*n ] ) ) {
			Size const i = res_list.index( jump_res[ 2*n - 1 ] );
			Size const j = res_list.index( jump_res[ 2*n ] );
			jump_partners1.push_back( i );
			jump_partners2.push_back( j );
			Size const chain_i( get_chain( i, all_res_in_chain ) );
			Size const chain_j( get_chain( j, all_res_in_chain ) );
			runtime_assert( connection_domains[ chain_i ] != connection_domains[ chain_j ] );
			chain_connections.push_back( make_pair( chain_i, chain_j ) );
			connection_domains = get_connection_domains( chain_connections, all_res_in_chain.size() );
		}
	}

}

////////////////////////////////////////////////////////////////////////////////
Size
get_chain( Size const i, vector1< vector1< Size > > const & all_res_in_chain ) {
	for ( Size n = 1; n <= all_res_in_chain.size(); n++ ) {
		if ( all_res_in_chain[ n ].has_value( i ) ) return n;
	}
	return 0;
}

////////////////////////////////////////////////////////////////////////////////
void
setup_jumps( core::pose::Pose const & pose,
	vector1< Size > & jump_partners1,
	vector1< Size > & jump_partners2,
	vector1< pair< Size, Size > > & chain_connections,
	vector1< vector1< Size > > const & all_res_in_chain,
	pair< vector1< int >, vector1< char > > const & resnum_and_chain_in_pose ) {

	Size const num_chains = all_res_in_chain.size();
	if ( jump_partners1.size() == num_chains - 1 ) return;

	// order of preference for chain connections:
	//
	// 1. start with chains that have maximum number of contacts, then
	// 2. choose chains that are *closest in sequence*, as defined
	//      by end of chain i and start of chain j.
	// 3. break remaining ties based on chain order.
	//
	// For jump connection residues, prefer residues that have the most contacts.
	//
	vector1< int  > const & conventional_numbering_in_pose = resnum_and_chain_in_pose.first;
	vector1< char > const & conventional_chains_in_pose = resnum_and_chain_in_pose.second;

	// Data structure set up for sorting...
	//
	// chain_pair = ( -num_contacts, sequence_separation, ( ( chain_idx i, chain_idx j ), ( res in chain i, res in chain j ) ) )
	//
	utility::vector1< pair< int, pair< Size, pair< pair< Size, Size >, pair< Size, Size > > > > > chain_pairs;
	Size const max_seq_separation = max( conventional_numbering_in_pose ) - min( conventional_numbering_in_pose ) + 1;
	for ( Size i = 1; i < num_chains; i++ ) {
		vector1< Size > const & fixed_res_in_chain_i = all_res_in_chain[ i ];
		if ( fixed_res_in_chain_i.size() == 0 ) continue;

		for ( Size j = i+1; j <= num_chains; j++ ) {
			pair< Size, Size > jump_res_pair( 0, 0 );
			vector1< Size > const & fixed_res_in_chain_j = all_res_in_chain[ j ];
			if ( fixed_res_in_chain_j.size() == 0 ) continue;

			// num_contacts
			static Distance const CONTACT_DIST_CUTOFF( 4.0 );
			int num_contacts( 0 );
			vector1< pair< Size, pair< Size, Size > > > num_contacts_pairwise;
			for ( Size m = 1; m <= fixed_res_in_chain_i.size(); m++ ) {
				core::conformation::Residue rsd_i = pose.residue( fixed_res_in_chain_i[ m ] );
				for ( Size n = 1; n <= fixed_res_in_chain_j.size(); n++ ) {
					core::conformation::Residue rsd_j = pose.residue( fixed_res_in_chain_j[ n ] );
					Size count( 0 );
					for ( Size mm = 1; mm <= rsd_i.nheavyatoms(); mm++ ) {
						for ( Size nn = 1; nn <= rsd_j.nheavyatoms(); nn++ ) {
							if ( ( rsd_i.xyz( mm ) - rsd_j.xyz( nn ) ).length() < CONTACT_DIST_CUTOFF ) {
								count++;
							}
						}
					}
					num_contacts_pairwise.push_back( make_pair( count, make_pair( fixed_res_in_chain_i[m],
						fixed_res_in_chain_j[n] ) ) );
					num_contacts += count;
				}
			}
			// find residue pair with the most contacts.
			std::sort( num_contacts_pairwise.begin(), num_contacts_pairwise.end() );
			if ( num_contacts_pairwise.size() > 0 && num_contacts > 0 ) jump_res_pair = num_contacts_pairwise[ 1 ].second;

			// figure out sequence separation
			Size sequence_separation( max_seq_separation );
			Size const chain_i_end   = fixed_res_in_chain_i[ fixed_res_in_chain_i.size() ];
			Size const chain_j_begin = fixed_res_in_chain_j[ 1 ];
			if ( conventional_chains_in_pose[ chain_i_end ] == conventional_chains_in_pose[ chain_j_begin ] ) {
				//runtime_assert( conventional_numbering_in_pose[ chain_j_begin ] > conventional_numbering_in_pose[ chain_i_end ] );
				sequence_separation = std::abs( conventional_numbering_in_pose[ chain_j_begin ] - conventional_numbering_in_pose[ chain_i_end ] );
			}
			if ( jump_res_pair == make_pair( Size(0), Size(0) ) ) jump_res_pair = make_pair( chain_i_end, chain_j_begin);

			chain_pairs.push_back( make_pair( -num_contacts /*want to maximize contact*/,
				make_pair( sequence_separation,
				make_pair( make_pair( i, j ), jump_res_pair ) ) ) );
		}
	}

	std::sort( chain_pairs.begin(), chain_pairs.end() );

	utility::vector1< Size > connection_domains = get_connection_domains( chain_connections, all_res_in_chain.size() );
	for ( Size q = 1; q <= chain_pairs.size(); q++ ) {
		Size const & i = chain_pairs[ q ].second.second.first.first;
		Size const & j = chain_pairs[ q ].second.second.first.second;
		if ( connection_domains[ i ] == connection_domains[ j ] ) continue;
		Size const & jump_res_i = chain_pairs[ q ].second.second.second.first;
		Size const & jump_res_j = chain_pairs[ q ].second.second.second.second;
		jump_partners1.push_back( jump_res_i );
		jump_partners2.push_back( jump_res_j );
		chain_connections.push_back( make_pair( i, j ) );
		connection_domains = get_connection_domains( chain_connections, all_res_in_chain.size() );
		if ( jump_partners1.size() == num_chains - 1 ) break;
	}

}

/////////////////////////////////////////////////////////////////////////////////
FoldTree
get_tree( pose::Pose const & pose,
	vector1< Size > const & cuts,
	vector1< Size > const & jump_partners1,
	vector1< Size > const & jump_partners2 ) {

	vector1< std::string > jump_atoms1, jump_atoms2;
	for ( Size n = 1; n <= jump_partners1.size(); n++ ) {
		jump_atoms1.push_back( chemical::rna::default_jump_atom( pose.residue_type( jump_partners1[n] ) ) );
		jump_atoms2.push_back( chemical::rna::default_jump_atom( pose.residue_type( jump_partners2[n] ) ) );
	}
	return get_tree( pose.size(), cuts, jump_partners1, jump_partners2, jump_atoms1, jump_atoms2 );
}


/////////////////////////////////////////////////////////////////////////////////
FoldTree
get_tree( Size const nres,
	vector1< Size > const & cuts,
	vector1< Size > const & jump_partners1,
	vector1< Size > const & jump_partners2,
	vector1< std::string > const & jump_atoms1,
	vector1< std::string > const & jump_atoms2 ) {

	Size const num_cuts = cuts.size();

	FoldTree f;
	ObjexxFCL::FArray2D< int > jump_point_( 2, num_cuts, 0 );
	ObjexxFCL::FArray1D< int > cuts_( num_cuts, 0 );
	for ( Size i = 1; i <= num_cuts; i++ ) {
		jump_point_( 1, i ) = std::min( jump_partners1[ i ], jump_partners2[ i ] );
		jump_point_( 2, i ) = std::max( jump_partners1[ i ], jump_partners2[ i ] );
		cuts_( i ) = cuts[ i ];
	}
	f.tree_from_jumps_and_cuts( nres, num_cuts, jump_point_, cuts_ );

	bool const KeepStubInResidue( true );
	for ( Size i = 1; i <= num_cuts; i++ ) {
		Size const n = f.jump_nr( jump_partners1[ i ], jump_partners2[ i ] );
		f.set_jump_atoms( n,
			jump_partners1[ i ], jump_atoms1[ i ],
			jump_partners2[ i ], jump_atoms2[ i ], KeepStubInResidue );
	}
	f.reassign_atoms_for_intra_residue_stubs(); // it seems silly that we need to do this separately.

	return f;
}

///////////////////////////////////////////////////////////////////////////////////////
void
update_fixed_domain_from_extra_minimize_jump_pairs( utility::vector1< Size > & fixed_domain,
	pose::Pose const & pose,
	utility::vector1< Size > const & res_list,
	utility::vector1< pair< Size, Size > > const & extra_minimize_jump_pairs ) {


	for ( Size i = 1; i <= extra_minimize_jump_pairs.size(); i++ ) {
		Size const & res1_full = extra_minimize_jump_pairs[ i ].first;
		Size const & res2_full = extra_minimize_jump_pairs[ i ].second;

		if ( !res_list.has_value( res1_full ) ) continue;
		if ( !res_list.has_value( res2_full ) ) continue;

		runtime_assert( fixed_domain[ res1_full ] > 0 );
		runtime_assert( fixed_domain[ res2_full ] > 0 );
		if ( fixed_domain[ res1_full ] != fixed_domain[ res2_full ] ) continue; // don't need to do anything.

		Size const res1 = res_list.index( res1_full );
		Size const res2 = res_list.index( res2_full );
		runtime_assert( pose.fold_tree().jump_exists( res1, res2 ) );
		Size jump_nr = pose.fold_tree().jump_nr( res1, res2 ); // need to check both orders!
		if ( jump_nr == 0 ) jump_nr = pose.fold_tree().jump_nr( res2, res1 );
		runtime_assert( jump_nr > 0 );

		utility::vector1< bool > partition_definition = pose.fold_tree().partition_by_jump( jump_nr );
		utility::vector1< Size > residues_to_update;
		Size const downstream_res = pose.fold_tree().downstream_jump_residue( jump_nr );
		for ( Size n = 1; n <= pose.size(); n++ ) {
			if ( partition_definition[ n ] == partition_definition[ downstream_res ] &&
					fixed_domain[ res_list[ n ] ] ==         fixed_domain[ res_list[ downstream_res ] ] ) {
				residues_to_update.push_back( n );
			}
		}
		runtime_assert( residues_to_update.size() > 0 );
		Size const new_domain = max( fixed_domain ) + 1;
		for ( Size n = 1; n <= residues_to_update.size(); n++ ) {
			fixed_domain[ res_list[ residues_to_update[ n ] ] ] = new_domain;
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////
void
update_fixed_domain_from_extra_minimize_jump_res( vector1< Size > & fixed_domain,
	pose::Pose const & pose,
	vector1< Size > const & res_list,
	vector1< Size > const & extra_minimize_jump_res ) {
	utility::vector1< pair< Size, Size > > extra_minimize_jump_pairs;
	for ( Size n = 1; n <= extra_minimize_jump_res.size()/2; n++ ) {
		extra_minimize_jump_pairs.push_back( make_pair( extra_minimize_jump_res[ 2*n - 1 ],
			extra_minimize_jump_res[ 2*n     ] ) );
	}
	update_fixed_domain_from_extra_minimize_jump_pairs( fixed_domain, pose, res_list, extra_minimize_jump_pairs );
}

/////////////////////////////////////////////////////////////////////////////////////
void
add_cutpoint_closed( pose::Pose & pose,
	vector1< Size > const & res_list,
	vector1< Size > const & cutpoint_closed ) {
	using namespace core::pose::rna;
	for ( Size n = 1; n <= cutpoint_closed.size(); n++ ) {
		if ( !res_list.has_value( cutpoint_closed[ n ] ) ) continue;
		Size const i = res_list.index( cutpoint_closed[ n ] );
		// could be useful in general -- share this with TransientCutpointHandler?
		vector1< pair< TorsionID, Real > > const suite_torsion_info = get_suite_torsion_info( pose, i );
		put_in_cutpoint( pose, i );
		correctly_add_cutpoint_variants( pose, i );
		apply_suite_torsion_info( pose, suite_torsion_info );
	}
}


/////////////////////////////////////////////////////////////////////////////////////
void
put_in_cutpoint( pose::Pose & pose, Size const i ) {
	core::kinematics::FoldTree f = pose.fold_tree();
	Size const new_jump = f.new_jump( i, i+1, i );
	if ( pose.residue_type( i ).is_RNA() && pose.residue_type( i+1 ).is_RNA() ) {
		f.set_jump_atoms( new_jump,
			chemical::rna::default_jump_atom( pose.residue_type( f.upstream_jump_residue( new_jump ) ) ),
			chemical::rna::default_jump_atom( pose.residue_type( f.downstream_jump_residue( new_jump ) ) ) );
	}
	pose.fold_tree( f );
}

/////////////////////////////////////////////////////////////////////////////////////
void
add_virtual_sugar_res( pose::Pose & pose,
	vector1< Size > const & res_list,
	vector1< Size > const & virtual_sugar_res )
{
	for ( Size n = 1; n <= virtual_sugar_res.size(); n++ ) {
		if ( !res_list.has_value( virtual_sugar_res[ n ] ) ) continue;
		Size const i = res_list.index( virtual_sugar_res[ n ] );
		runtime_assert( i == 1 || pose.fold_tree().is_cutpoint( i - 1 ) );
		runtime_assert( i == pose.size() || pose.fold_tree().is_cutpoint( i ) );
		add_variant_type_to_pose_residue( pose, core::chemical::VIRTUAL_RIBOSE, i );
	}
}

/////////////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
figure_out_working_res( utility::vector1< Size > const & input_domain_map,
	utility::vector1< Size > const & sample_res ) {
	vector1< Size > working_res;
	for ( Size n = 1; n <= input_domain_map.size(); n++ ) {
		if ( input_domain_map[ n ] == 0  && !sample_res.has_value( n ) ) continue;
		working_res.push_back( n );
	}
	return working_res;
}

/////////////////////////////////////////////////////////////////////////////////////
// 'default': sample any residues that are not inputted as fixed PDBs.
utility::vector1< Size >
figure_out_sample_res( utility::vector1< Size > const & input_domain_map,
	utility::vector1< Size > const & working_res ) {
	vector1< Size > sample_res;
	for ( Size n = 1; n <= input_domain_map.size(); n++ ) {
		if ( working_res.size() > 0 && !working_res.has_value( n ) ) continue;
		if ( input_domain_map[ n ] == 0  ) sample_res.push_back( n );
	}
	return sample_res;
}

/////////////////////////////////////////////////////////////////////////////////////
void
check_working_res( utility::vector1< Size > const & working_res,
	utility::vector1< Size > const & input_domain_map,
	utility::vector1< Size > const & sample_res ) {
	for ( Size i = 1; i <= sample_res.size(); i++ ) {
		runtime_assert( working_res.has_value( sample_res[ i ] ) );
	}
	for ( Size n = 1; n <= input_domain_map.size(); n++ ) {
		if ( input_domain_map[ n ] > 0 ) {
			if ( !working_res.has_value( n ) /* && !reference_pose */ ) {
				utility_exit_with_message( "Working res does not have input_domain_map residue "+ObjexxFCL::string_of(n) );
			}
			if ( sample_res.has_value( n ) ) utility_exit_with_message( "Sample res should not have "+ObjexxFCL::string_of(n) );
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////
// motif_mode:
//
//  extra_min_res -- minimize residues in input PDBs immediately neighboring loops to be built.
//  terminal_res  -- don't stack/pair loops on any cutpoint residues in input PDBs (except right next to loops)
//
// originally worked out in python (setup_stepwise_benchmark.py) by rhiju, 2014.
//
void
figure_out_motif_mode( utility::vector1< Size > & extra_min_res,
	utility::vector1< Size > & terminal_res,
	utility::vector1< Size > const & working_res,
	utility::vector1< Size > const & input_domain_map,
	utility::vector1< Size > const & cutpoint_open_in_full_model ) {
	// for double-checking work against any user-inputted extra_min_res & terminal_res
	utility::vector1< Size > const extra_min_res_save = extra_min_res;
	utility::vector1< Size > const terminal_res_save = terminal_res;
	extra_min_res.clear();
	terminal_res.clear();

	for ( Size m = 1; m <= input_domain_map.size(); m++ ) {
		if ( input_domain_map[ m ] == 0 ) continue;
		bool const right_before_chainbreak = ( m == input_domain_map.size() || cutpoint_open_in_full_model.has_value( m ) ||
			( working_res.size() > 0 && !working_res.has_value( m + 1 ) ) );
		bool const right_after_chainbreak  = ( m == 1                 || cutpoint_open_in_full_model.has_value( m - 1 ) ||
			( working_res.size() > 0 && !working_res.has_value( m - 1 )  ) );
		bool const prev_moving = !right_after_chainbreak && m > 1 &&
			( ( input_domain_map[ m - 1 ] == 0  && ( working_res.size() == 0 || working_res.has_value( m - 1 ) ) ) ||
			( input_domain_map[ m - 1 ] != input_domain_map[ m ] ) ) ;
		bool const next_moving = m < input_domain_map.size() && !right_before_chainbreak &&
			( ( input_domain_map[ m + 1 ] == 0 &&  ( working_res.size() == 0 || working_res.has_value( m + 1 ) ) ) ||
			( input_domain_map[ m + 1 ] != input_domain_map[ m ] ) );
		if ( ( right_after_chainbreak  && !right_before_chainbreak && !next_moving ) ||
				( right_before_chainbreak && !right_after_chainbreak  && !prev_moving ) ) {
			terminal_res.push_back( m );
		}
		if ( ( prev_moving && !next_moving && !right_before_chainbreak ) ||
				( next_moving && !prev_moving && !right_after_chainbreak ) ) {
			extra_min_res.push_back( m );
		}
	}

	//  if ( looks_like_reference_pose( input_domain_map ) ) return;
	if ( extra_min_res_save.size() > 0 ) {
		if ( extra_min_res != extra_min_res_save ) {
			TR << TR.Red << "Auto-computed EXTRA_MIN_RES " << extra_min_res <<      TR.Reset << std::endl;
			TR << TR.Red << "User-input    EXTRA_MIN_RES " << extra_min_res_save << TR.Reset << std::endl;
		}
		runtime_assert( extra_min_res_save.size() == 0 || extra_min_res == extra_min_res_save );
	}
	if ( terminal_res_save.size() > 0 ) {
		if ( terminal_res != terminal_res_save ) {
			TR << TR.Red << "Auto-computed TERMINAL_RES " << terminal_res <<      TR.Reset << std::endl;
			TR << TR.Red << "User-input    TERMINAL_RES " << terminal_res_save << TR.Reset << std::endl;
		}
	}
	runtime_assert( terminal_res_save.size() == 0 || terminal_res == terminal_res_save );

}

//////////////////////////////////////////////////////////////////////////////////////////////////
// An 'advanced' version of terminal_res -- places 'repulsion atoms' above or below base faces to
// prevent stacking during both sampling and minimization.
void
add_block_stack_variants( vector1< pose::Pose * > const & pose_pointers,
	vector1< vector1< Size > > const & pose_res_lists,
	vector1< Size > const & block_stack_above_res,
	vector1< Size > const & block_stack_below_res ) {
	using namespace core::chemical;
	for ( Size n = 1; n <= pose_pointers.size(); n++ ) {
		pose::Pose & pose = *( pose_pointers[ n ] ) ;
		vector1 < Size > const & pose_res_list = pose_res_lists[ n ] ;
		for ( Size i = 1; i <= pose_res_list.size(); i++ ) {
			if ( block_stack_above_res.has_value( pose_res_list[ i ] ) ) {
				runtime_assert( pose.residue_type( i ).is_RNA() );
				add_variant_type_to_pose_residue( pose, BLOCK_STACK_ABOVE, i );
			}
			if ( block_stack_below_res.has_value( pose_res_list[ i ] ) ) {
				runtime_assert( pose.residue_type( i ).is_RNA() );
				add_variant_type_to_pose_residue( pose, BLOCK_STACK_BELOW, i );
			}
		}
	}
}

/////////////////////////////////////////////////////////////////////////////
void
update_jump_res( utility::vector1< Size > & jump_res,
	utility::vector1< Size > const & extra_minimize_jump_res ) {

	runtime_assert( jump_res.size() % 2 == 0 );
	runtime_assert( extra_minimize_jump_res.size() % 2 == 0 );

	for ( Size n = 1; n <= extra_minimize_jump_res.size()/2; n++ ) {
		Size const & res1 = extra_minimize_jump_res[ 2*n - 1 ];
		Size const & res2 = extra_minimize_jump_res[ 2*n     ];
		bool matches_existing( false );
		for ( Size q = 1; q <= jump_res.size()/2; q++ ) {
			Size const & existing_res1 = jump_res[ 2*n - 1 ];
			Size const & existing_res2 = jump_res[ 2*n     ];
			if ( ( res1 == existing_res1 && res2 == existing_res2 ) ||
					( res1 == existing_res2 && res2 == existing_res1 ) ) {
				matches_existing = true;
				break;
			}
		}
		if ( matches_existing ) continue;
		jump_res.push_back( res1 );
		jump_res.push_back( res2 );
	}
}

/////////////////////////////////////////////////////////////////////////////
void
check_extra_minimize_res_are_input( utility::vector1< core::Size > const & extra_minimize_res,
	utility::vector1< core::Size > const & input_domain_map ) {
	for ( Size n = 1; n <= extra_minimize_res.size(); n++ ) {
		runtime_assert( input_domain_map[ extra_minimize_res[ n ] ] > 0 );
	}
}

/////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
figure_out_fixed_domain_map( utility::vector1< Size > const & input_domain_map,
	utility::vector1< Size > const & extra_minimize_res ) {
	vector1< Size > fixed_domain_map = input_domain_map;
	for ( Size n = 1; n <= extra_minimize_res.size(); n++ ) {
		fixed_domain_map[ extra_minimize_res[ n ] ] = 0;
	}
	return fixed_domain_map;
}

/////////////////////////////////////////////////////////////////////////////
utility::vector1< Size >
figure_out_dock_domain_map( utility::vector1< Size > & cutpoint_open_in_full_model,
	utility::vector1< utility::vector1< Size > > const & pose_res_lists,
	utility::vector1< Size > const & working_res,
	utility::vector1< Size > const & sample_res,
	Size const nres ) {

	// supplement cutpoint_open with breaks in working_res.
	vector1< Size > chains( nres, 0 ); // 0 outside working_res
	Size chain_number( 0 );
	for ( Size n = 1; n <= nres; n++ ) {
		if ( !working_res.has_value( n ) ) continue;
		bool new_chain( false );
		if ( ( n == 1 || !working_res.has_value( n - 1 )) && working_res.has_value( n ) ) new_chain = true;
		if ( n > 1 && cutpoint_open_in_full_model.has_value( n - 1 ) ) new_chain = true;
		if ( new_chain ) chain_number++;
		chains[ n ]= chain_number; // every working_res gets a chain number.
	}

	// which segments are connected by input poses?
	vector1< pair< Size, Size > > chain_connections;
	for ( Size n = 1; n <= pose_res_lists.size(); n++ ) {
		vector1< Size > const & res_list = pose_res_lists[ n ];
		std::set< Size > chains_in_pose;
		for ( Size k = 1; k <= res_list.size(); k++ ) {
			if ( sample_res.has_value( res_list[ k ] ) ) continue;
			if ( !working_res.has_value( res_list[ k ] ) ) continue;
			chains_in_pose.insert( chains[ res_list[ k ] ] );
		}
		for ( std::set< Size >::const_iterator it1 = chains_in_pose.begin(), end = chains_in_pose.end(); it1 != end; ++it1 ) {
			for ( std::set< Size >::const_iterator it2 = it1; it2 != end; ++it2 ) {
				if ( it1 != it2 ) chain_connections.push_back( make_pair( *it1, *it2 ) );
			}
		}
	}
	// find connected clusters
	utility::vector1< Size > connection_domains = get_connection_domains( chain_connections, max( chains ) );

	vector1< Size > dock_domain_map( nres, 0 );
	for ( Size n = 1; n <= nres; n++ ) {
		// to get from 1,2,... (connection_domains) to 0,1,... (convention in FullModelParameters)
		if ( !working_res.has_value( n ) ) continue;
		dock_domain_map[ n ] = connection_domains[ chains[ n ] ];
	}

	for ( Size k = 1; k < working_res.size(); k++ ) {
		Size const n = working_res[ k ];
		Size const n_next = working_res[ k+1 ];
		if (  dock_domain_map[ n ] != dock_domain_map[ n_next ] &&
				!cutpoint_open_in_full_model.has_value( n ) ) {
			TR << TR.Red << "There appears to be a dock boundary at " << n << " so adding to cutpoint_in_full_model." << TR.Reset << std::endl;
			cutpoint_open_in_full_model.push_back( n );
		}
	}

	return dock_domain_map;
}

/////////////////////////////////////////////////////////////////////////////
void
reorder_pose( pose::Pose & pose, utility::vector1< Size > & res_list ) {
	utility::vector1< Size > res_list_ordered = res_list;
	std::sort( res_list_ordered.begin(), res_list_ordered.end() );
	if ( res_list == res_list_ordered ) return;

	vector1< Size > slice_res;
	for ( Size n = 1; n <= res_list_ordered.size(); n++ ) slice_res.push_back( res_list.index( res_list_ordered[ n ] ) );
	pose::pdbslice( pose, slice_res );  // will do the rearrangement.
	res_list = res_list_ordered;
}

/////////////////////////////////////////////////////////////////////////////
// helper function to see if this is just an RNA modeling problem, as defined
// in fasta file.
bool
just_modeling_RNA( utility::vector1< std::string > const & fasta_files ) {
	std::string sequence = core::sequence::read_fasta_file_return_str( fasta_files[1] );
	core::sequence::parse_out_non_standard_residues( sequence );
	return ( modeler::rna::just_modeling_RNA( sequence ) );
}

/////////////////////////////////////////////////////////////////////////////
//calebgeniesse: setup for density scoring
void
setup_for_density_scoring( core::pose::Pose & pose ) {
	// add virtual root for density scoring
	TR << "Adding virtual residue as root" << std::endl;
	pose::addVirtualResAsRoot( pose );
	// Note from rhiju -- should we update full_model_info & full_model_parameters?
	// I put a new function in core/pose/full_model_info/util.hh called append_virtual_residue_to_full_model_info().
}

} //setup
} //stepwise
} //protocols
