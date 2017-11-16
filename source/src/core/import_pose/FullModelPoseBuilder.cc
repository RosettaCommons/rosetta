// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/import_pose/FullModelPoseBuilder
/// @brief  create a pose from PDB in a way that sets up FullModelInfo
/// @author Andrew Watkins

// Unit headers
#include <core/import_pose/FullModelPoseBuilder.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/import_pose_options.hh>

// Project headers
#include <core/types.hh>

#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/AA.hh>

#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/PDBInfo.hh>

#include <core/io/pdb/pdb_reader.hh>
#include <core/io/pdb/RecordType.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/StructFileReaderOptions.hh>  // TODO: Rename after refactor.
#include <core/io/pose_from_sfr/PoseFromSFRBuilder.hh>
#include <core/io/mmcif/cif_reader.hh>
#include <core/io/StructFileRep.hh>
#include <core/io/StructFileRepOptions.hh>

#include <core/io/pdb/pdb_reader.hh>

#include <core/id/TorsionID.hh>

// Basic headers
#include <basic/Tracer.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/full_model.OptionKeys.gen.hh>
#include <basic/options/keys/stepwise.OptionKeys.gen.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/keys/magnesium.OptionKeys.gen.hh>
#include <basic/options/option.hh>

// Utility headers
#include <utility/exit.hh>
#include <utility/string_util.hh>
#include <utility/io/izstream.hh>
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>

// External headers
#include <ObjexxFCL/string.functions.hh>
#include <cifparse/CifFile.h>
#include <cifparse/CifParserBase.h>

static basic::Tracer TR( "core.import_pose.util" );

namespace core {
namespace import_pose {

using namespace kinematics;

using core::Size;
using core::SSize;

using basic::Error;
using basic::Warning;
using namespace core::io;
using namespace ObjexxFCL;
using namespace core::pose;
using namespace core::pose::full_model_info;

static basic::Tracer TR( "core.import_pose.import_pose" );

using utility::vector1;

FullModelPoseBuilder::FullModelPoseBuilder():
	options_( basic::options::option )
{}

void
FullModelPoseBuilder::initialize_input_poses_from_options( core::chemical::ResidueTypeSetCAP rsd_set ) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::io::silent;
	using namespace core::pose;
	using namespace utility;

	utility::vector1< std::string > const & input_pdb_files    = options_[ in::file::s ]();
	utility::vector1< std::string > const & input_silent_files = options_[ in::file::silent ]();

	for ( Size n = 1; n <= input_pdb_files.size(); n++ ) {
		input_poses_.push_back( get_pdb_and_cleanup( input_pdb_files[ n ], rsd_set ) );
	}
	for ( Size n = 1; n <= input_silent_files.size(); n++ ) {
		PoseOP pose( new Pose );
		SilentFileOptions opts;
		SilentFileData silent_file_data(opts);
		core::chemical::ResidueTypeSetCOP rsd_set_op( rsd_set );
		silent_file_data.read_file( input_silent_files[n] );
		silent_file_data.begin()->fill_pose( *pose, *rsd_set_op );
		input_poses_.push_back( pose );
	}

	if ( options_[ full_model::other_poses ].user() ) {
		get_other_poses( input_poses_, options_[ full_model::other_poses ](), rsd_set );
	}
}

// AMW: todo, if you ever change input_poses you've got to re-think-through the initialize_further
// bit too.

void
FullModelPoseBuilder::initialize_full_model_parameters() {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	//
	if ( fasta_file_ != "" ) {
		full_model_parameters_ = get_sequence_information( fasta_file_, cutpoint_open_in_full_model_ );
		// Could also have set fasta_file_ earlier.
	} else if ( options_[ in::file::fasta ].user() ) {
		fasta_file_ = options_[ in::file::fasta ]()[1];
		full_model_parameters_ = get_sequence_information( fasta_file_, cutpoint_open_in_full_model_ );
	} else {
		// guess sequence, chain, resnum from pose; not specified in fasta.
		//runtime_assert( pose_pointers.size() > 0 );
		runtime_assert( input_poses_.size() > 0 );
		vector1< Size > dummy;
		full_model_parameters_ = FullModelParametersOP( new FullModelParameters( *input_poses_[1], dummy ) );
		if ( input_poses_.size() > 1 ) {
			utility_exit_with_message( "Currently need to specify -fasta if dealing with multiple poses. Would not be hard to fix this, by merging resnum/chain across poses! Ask rhiju." );
		}
	}
}

void
FullModelPoseBuilder::initialize_further_from_options() {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::io::silent;
	using namespace core::pose;
	using namespace utility;

	std::tuple< vector1< Size >, vector1< char >, vector1< std::string > > const & input_resnum_and_chain_and_segid = options_[ in::file::input_res ].resnum_and_chain();

	set_input_resnum_and_chain_and_segid( input_resnum_and_chain_and_segid );
	set_cutpoint_open_in_full_model( options_[ full_model::cutpoint_open ]() );


	initialize_full_model_parameters();


	set_extra_minimize_res( full_model_parameters_->conventional_to_full( options_[ full_model::extra_min_res ].resnum_and_chain() ) );
	set_sample_res( full_model_parameters_->conventional_to_full( options_[ full_model::sample_res ].resnum_and_chain() ) ); //stuff that can be resampled.
	set_working_res( full_model_parameters_->conventional_to_full( options_[ full_model::working_res ].resnum_and_chain() ) ); //all working stuff
	set_terminal_res( full_model_parameters_->conventional_to_full( options_[ full_model::rna::terminal_res ].resnum_and_chain() ) );
	set_block_stack_above_res( full_model_parameters_->conventional_to_full( options_[ full_model::rna::block_stack_above_res ].resnum_and_chain() ) );
	set_block_stack_below_res( full_model_parameters_->conventional_to_full( options_[ full_model::rna::block_stack_below_res ].resnum_and_chain() ) );
	set_preferred_root_res( full_model_parameters_->conventional_to_full( options_[ full_model::root_res ].resnum_and_chain() ) );
	set_jump_res( full_model_parameters_->conventional_to_full( options_[ full_model::jump_res ].resnum_and_chain() ) );
	set_cutpoint_closed( full_model_parameters_->conventional_to_full( options_[ full_model::cutpoint_closed ].resnum_and_chain() ) );
	if ( options_[ full_model::cyclize ].user() ) utility_exit_with_message( "Cannot handle cyclize yet in stepwise." );
	set_fiveprime_res( full_model_parameters_->conventional_to_full( options_[ full_model::fiveprime_cap ].resnum_and_chain() ) );
	set_bulge_res( full_model_parameters_->conventional_to_full( options_[ full_model::rna::bulge_res ].resnum_and_chain() ) );
	set_extra_minimize_jump_res( full_model_parameters_->conventional_to_full( options_[ full_model::extra_min_jump_res ].resnum_and_chain() ) );
	set_virtual_sugar_res( full_model_parameters_->conventional_to_full( options_[ full_model::virtual_sugar_res ].resnum_and_chain() ) );
	set_alignment_anchor_res( full_model_parameters_->conventional_to_full( options_[ OptionKeys::stepwise::alignment_anchor_res ].resnum_and_chain() ) );

	set_calc_rms_res( full_model_parameters_->conventional_to_full( options_[ full_model::calc_rms_res ].resnum_and_chain() ) );
	set_rna_syn_chi( full_model_parameters_->conventional_to_full( options_[ full_model::rna::force_syn_chi_res_list ].resnum_and_chain() ) );
	set_rna_anti_chi( full_model_parameters_->conventional_to_full( options_[ full_model::rna::force_anti_chi_res_list ].resnum_and_chain() ) );
	set_rna_north_sugar( full_model_parameters_->conventional_to_full( options_[ full_model::rna::force_north_sugar_list ].resnum_and_chain() ) );
	set_rna_south_sugar( full_model_parameters_->conventional_to_full( options_[ full_model::rna::force_south_sugar_list ].resnum_and_chain() ) );
	set_rna_sample_sugar( full_model_parameters_->conventional_to_full( options_[ full_model::rna::sample_sugar_res ].resnum_and_chain() ) );

	// earlier had set this based on
	if ( options_[ full_model::cutpoint_open ].user() ) {
		set_cutpoint_open_in_full_model( full_model_parameters_->conventional_to_full( options_[ full_model::cutpoint_open ].resnum_and_chain() ) );
	}

	//set_full_model_parameters( full_model_parameters_ );

	set_disulfide_file( options_[ OptionKeys::stepwise::protein::disulfide_file ]() );
	if ( options_[ OptionKeys::constraints::cst_file ].user() ) set_constraint_file( options_[ OptionKeys::constraints::cst_file ]()[ 1 ] );
}

PoseOP FullModelPoseBuilder::build() { // can't be const because it can change input_poses_

	vector1< Size > const & input_res_list = std::get< 0 >( input_resnum_and_chain_and_segid_ );
	if ( input_res_list.size() ) {
		vector1< char > const & input_chain_list = std::get< 1 >( input_resnum_and_chain_and_segid_ );
		vector1< std::string > const & input_segid_list = std::get< 2 >( input_resnum_and_chain_and_segid_ );
		Size input_res_count = 0;
		for ( Size n = 1; n <= input_poses_.size(); n++ ) {
			Pose & pose = *input_poses_[ n ];
			PDBInfoOP pdb_info( new PDBInfo( pose ) );
			vector1< Size > input_res_for_pose;
			vector1< char > input_chain_for_pose;
			vector1< std::string> input_segid_for_pose;
			for ( Size k = 1; k <= pose.size(); k++ ) {
				input_res_count++;
				runtime_assert( input_res_count <= input_res_list.size() );
				Size const & number_in_full_model = input_res_list[ input_res_count ];
				input_res_for_pose.push_back( number_in_full_model );
				input_chain_for_pose.push_back( input_chain_list[ input_res_count ] );
				input_segid_for_pose.push_back( input_segid_list[ input_res_count ] );
			}
			pdb_info->set_numbering( input_res_for_pose );
			pdb_info->set_chains(    input_chain_for_pose );
			pdb_info->set_segment_ids(    input_segid_for_pose );
			pose.pdb_info( pdb_info );
		}
		runtime_assert( input_res_count == input_res_list.size() );
	}

	if ( input_poses_.empty() ) input_poses_.push_back( core::pose::PoseOP( new Pose ) ); // just a blank pose for now.



	// AAA

	fill_full_model_info( input_poses_ );  //FullModelInfo (minimal object needed for add/delete)



	return input_poses_[1];
}


///////////////////////////////////////////////////////////////////////////////////////
void
FullModelPoseBuilder::fill_full_model_info( pose::Pose & pose ) {
	utility::vector1< PoseOP > other_pose_ops; // dummy
	fill_full_model_info( pose, other_pose_ops );
}

///////////////////////////////////////////////////////////////////////////////////////
void
FullModelPoseBuilder::fill_full_model_info( vector1< PoseOP > & pose_ops ) {

	runtime_assert( pose_ops.size() > 0 );

	utility::vector1< PoseOP > other_pose_ops;
	for ( Size n = 2; n <= pose_ops.size(); n++ ) other_pose_ops.push_back( pose_ops[n] );

	fill_full_model_info( *(pose_ops[1]), other_pose_ops );
}


///////////////////////////////////////////////////////////////////////////////////////
// this is needed so that first pose holds pointers to other poses.
void
FullModelPoseBuilder::fill_full_model_info( pose::Pose & pose, vector1< PoseOP > & other_pose_ops ) {

	utility::vector1< Pose * > pose_pointers;
	pose_pointers.push_back( & pose );
	for ( Size n = 1; n <= other_pose_ops.size(); n++ ) pose_pointers.push_back( & (*other_pose_ops[ n ]) );

	fill_full_model_info( pose_pointers );

	core::pose::full_model_info::nonconst_full_model_info( pose ).set_other_pose_list( other_pose_ops );
}

void
FullModelPoseBuilder::fill_full_model_info( vector1< Pose * > & pose_pointers ) {

	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::pose::full_model_info;


	std::string const & desired_sequence = full_model_parameters_->full_sequence();

	// calebgeniesse: setup for edensity scoring, if map is provided via cmd-line
	if ( options_[ edensity::mapfile ].user() ) {
		// update pose
		if ( pose_pointers[1]->total_residue() > 0 ) {
			setup_for_density_scoring( *pose_pointers[1] );
		}
	}


	// calebgeniesse: override preferred_root_res if mapfile provided via cmd-line (hacky)
	if ( options_[ edensity::mapfile ].user() ) {
		preferred_root_res_ = full_model_parameters_->conventional_to_full(
			std::make_tuple(
			utility::tools::make_vector1< int >(1),
			utility::tools::make_vector1< char >('z'),
			utility::tools::make_vector1< std::string >("    ")
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
		vector1< Size > res_list = full_model_parameters_->conventional_to_full( make_tuple( get_res_num_from_pdb_info( pose ), get_chains_from_pdb_info( pose ), get_segids_from_pdb_info( pose ) ) );
		reorder_pose( pose, res_list );
		pose_res_lists.push_back( res_list );
		for ( Size k = 1; k <= pose.size(); k++ ) {
			if ( !sample_res_.has_value( res_list[ k ] ) ) input_domain_map[ res_list[ k ] ] = n;
		}
	}

	// calebgeniesse: override input_domain_map if mapfile provided via cmd-line (hacky)
	if ( desired_sequence[desired_sequence.size()-1] == 'X' ) {
		input_domain_map[desired_sequence.size()] = 1;// input_domain_map.size() ? max(input_domain_map) + 1 : 1; // or just 1?
	}


	// PLACEMENT?
	// Should figure_out_motif_mode be a method?
	if ( options_[ full_model::motif_mode ]() ) figure_out_motif_mode( extra_minimize_res_, terminal_res_, working_res_, input_domain_map, cutpoint_open_in_full_model_ );

	// Assumes motif mode has already been figured out
	add_block_stack_variants( pose_pointers, pose_res_lists, block_stack_above_res_, block_stack_below_res_ );

	// calebgeniesse: hacky way of making sure virt root is not included in terminal_res/extra_minimize_res
	if ( desired_sequence[desired_sequence.size()-1] == 'X' ) {
		if ( terminal_res_.has_value( desired_sequence.size() ) ) {
			terminal_res_.erase( std::remove(terminal_res_.begin(), terminal_res_.end(), desired_sequence.size()), terminal_res_.end());
		}
		if ( extra_minimize_res_.has_value( desired_sequence.size() ) ) {
			extra_minimize_res_.erase( std::remove(extra_minimize_res_.begin(), extra_minimize_res_.end(), desired_sequence.size()), extra_minimize_res_.end());
		}
	}

	// everything that is not fixed is sampleable (unless -sample_res explicitly specified).
	if ( sample_res_.empty() )  sample_res_  = figure_out_sample_res( input_domain_map, working_res_ );
	if ( working_res_.empty() ) working_res_ = figure_out_working_res( input_domain_map, sample_res_ );
	vector1< Size > fixed_domain_map = figure_out_fixed_domain_map( input_domain_map, extra_minimize_res_ ); //remove extra_minimize_res.

	setup_fold_trees( pose_pointers, cutpoint_open_in_full_model_ /* can update */, fixed_domain_map /* can update */,
		cutpoint_closed_, extra_minimize_res_, extra_minimize_jump_res_,
		sample_res_, working_res_, jump_res_,
		preferred_root_res_, virtual_sugar_res_,
		*full_model_parameters_, pose_res_lists );

	vector1< Size > const dock_domain_map = figure_out_dock_domain_map( cutpoint_open_in_full_model_ /* can be updated */,
		pose_res_lists, working_res_, sample_res_, desired_nres );

	// some checks
	check_working_res( working_res_, input_domain_map, sample_res_ );
	check_extra_minimize_res_are_input( extra_minimize_res_,      input_domain_map );
	check_extra_minimize_res_are_input( extra_minimize_jump_res_, input_domain_map );

	full_model_parameters_->set_parameter( INPUT_DOMAIN,  input_domain_map );
	full_model_parameters_->set_parameter( FIXED_DOMAIN,  fixed_domain_map );
	full_model_parameters_->set_parameter( DOCK_DOMAIN,   dock_domain_map  );
	full_model_parameters_->set_parameter_as_res_list( CUTPOINT_OPEN, cutpoint_open_in_full_model_ );
	full_model_parameters_->set_parameter_as_res_list( EXTRA_MINIMIZE, extra_minimize_res_ );
	full_model_parameters_->set_parameter_as_res_list( SAMPLE, sample_res_ );
	full_model_parameters_->set_parameter_as_res_list( WORKING, working_res_ );
	full_model_parameters_->set_parameter_as_res_list( CALC_RMS, calc_rms_res_ );
	full_model_parameters_->set_parameter_as_res_list( PREFERRED_ROOT, preferred_root_res_ );
	full_model_parameters_->set_parameter_as_res_list( RNA_SYN_CHI, rna_syn_chi_ );
	full_model_parameters_->set_parameter_as_res_list( RNA_ANTI_CHI, rna_anti_chi_ );
	full_model_parameters_->set_parameter_as_res_list( RNA_NORTH_SUGAR, rna_north_sugar_ );
	full_model_parameters_->set_parameter_as_res_list( RNA_SOUTH_SUGAR, rna_south_sugar_ );
	full_model_parameters_->set_parameter_as_res_list( RNA_TERMINAL, terminal_res_ );
	full_model_parameters_->set_parameter_as_res_list( RNA_BLOCK_STACK_ABOVE, block_stack_above_res_ );
	full_model_parameters_->set_parameter_as_res_list( RNA_BLOCK_STACK_BELOW, block_stack_below_res_ );
	full_model_parameters_->set_parameter_as_res_list( RNA_BULGE,  bulge_res_ );
	full_model_parameters_->set_parameter_as_res_list( ALIGNMENT_ANCHOR_RES,  alignment_anchor_res_ );
	full_model_parameters_->set_parameter_as_res_list( RNA_SAMPLE_SUGAR, rna_sample_sugar_ );

	// Temporary. 'res_list_in_pairs' assumes that different pairs involve different residues, but that's
	//  not the case, actually. Also, jump_res is not apparently used elsewhere in the code except to
	//  as a sanity check in ResampleMover.
	// full_model_parameters->set_parameter_as_res_list_in_pairs( full_model_info::JUMP, jump_res );

	full_model_parameters_->set_parameter_as_res_list_in_pairs( EXTRA_MINIMIZE_JUMP, extra_minimize_jump_res_ );
	full_model_parameters_->set_parameter_as_res_list_in_pairs( FIVEPRIME_CAP,  fiveprime_res_ );
	full_model_parameters_->read_disulfides( disulfide_file_ );
	if ( constraint_file_ != "" ) full_model_parameters_->read_cst_file( constraint_file_ );




	// AMW: can't figure out how to move this yet.
	// move this code block somewhere else when ready.
	if ( options_[ basic::options::OptionKeys::stepwise::monte_carlo::vary_loop_length_frequency ] > 0.0 ) {
		// placeholder -- testing if loops can be 'evaporated'.
		vector1< Size > full_model_res_no_loops;
		for ( Size n = 1; n <= desired_nres; ++n ) {
			if ( clean_desired_seq[ n-1 ] != 'n' || input_domain_map[ n ] ) full_model_res_no_loops.push_back( n );
		}
		full_model_parameters_ = full_model_parameters_->slice( full_model_res_no_loops );
		for ( auto & pose_res_list : pose_res_lists ) {
			vector1< Size > pose_res_list_new;
			for ( Size const res : pose_res_list ) {
				if ( full_model_res_no_loops.has_value( res ) ) {
					pose_res_list_new.push_back( full_model_res_no_loops.index( res ) );
				}
			}
			pose_res_list = pose_res_list_new;
		}
	}

	// AMW: FINAL PHASE.
	// Can't go back now.
	for ( Size n = 1; n <= pose_pointers.size(); n++ ) {
		Pose & pose = *pose_pointers[n];
		FullModelInfoOP full_model_info_for_pose( new FullModelInfo( full_model_parameters_ ) );
		full_model_info_for_pose->set_res_list( pose_res_lists[ n ] );
		full_model_info_for_pose->update_submotif_info_list();
		set_full_model_info( pose, full_model_info_for_pose );
		update_pose_objects_from_full_model_info( pose ); // for output pdb or silent file (residue numbering), constraints, disulfides
		pose::fix_up_residue_type_variants( pose ); // for sample sugars...
	}
}

/////////////////////////////////////////////////////////////////////////////
void
FullModelPoseBuilder::set_extra_minimize_jump_res( utility::vector1< Size > const & extra_minimize_jump_res ) {

	runtime_assert( jump_res_.size() % 2 == 0 );
	runtime_assert( extra_minimize_jump_res.size() % 2 == 0 );

	extra_minimize_jump_res_ = extra_minimize_jump_res;

	for ( Size n = 1; n <= extra_minimize_jump_res.size()/2; n++ ) {
		Size const & res1 = extra_minimize_jump_res[ 2*n - 1 ];
		Size const & res2 = extra_minimize_jump_res[ 2*n     ];
		bool matches_existing( false );
		for ( Size q = 1; q <= jump_res_.size()/2; q++ ) {
			Size const existing_res1 = jump_res_[ 2*n - 1 ];
			Size const existing_res2 = jump_res_[ 2*n     ];
			if ( ( res1 == existing_res1 && res2 == existing_res2 ) ||
					( res1 == existing_res2 && res2 == existing_res1 ) ) {
				matches_existing = true;
				break;
			}
		}
		if ( matches_existing ) continue;
		jump_res_.push_back( res1 );
		jump_res_.push_back( res2 );
	}
}

} // namespace import_pose
} // namespace core
