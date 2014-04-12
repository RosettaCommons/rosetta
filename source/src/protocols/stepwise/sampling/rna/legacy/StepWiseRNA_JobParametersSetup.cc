// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_JobParametersSetup
/// @brief Sets up pose and job parameters for RNA stepwise building.
/// @detailed
/// @author Rhiju Das, Parin Sripakdeevong


//////////////////////////////////
#include <protocols/stepwise/sampling/rna/legacy/StepWiseRNA_JobParametersSetup.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Classes.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Util.hh>

#include <protocols/farna/RNA_ProtocolUtil.hh>
//////////////////////////////////
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/scoring/func/Func.fwd.hh>
#include <core/scoring/func/FadeFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>
#include <core/io/pdb/pose_io.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/scoring/rms_util.tmpl.hh>
#include <ObjexxFCL/string.functions.hh>
#include <ObjexxFCL/format.hh>

#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>


#include <utility/exit.hh>
#include <time.h>

#include <string>

using namespace core;
using core::Real;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Core routine for stepwise sampling of proteins (and probably other
// biopolymers soon). Take a starting pose and a list of residues to sample,
//  and comprehensively sample all backbone torsion angles by recursion.
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.stepwise.rna.StepWiseRNA_JobParametersSetup" );

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {
namespace legacy {

  //////////////////////////////////////////////////////////////////////////
  //constructor!
	StepWiseRNA_JobParametersSetup::StepWiseRNA_JobParametersSetup( utility::vector1< core::Size > const & moving_res_list,
																																		std::string const & full_sequence,
																																		utility::vector1< core::Size > const & input_res,
																																		utility::vector1< core::Size > const & input_res2,
																																		utility::vector1< core::Size > const & cutpoint_open,
																																		Size const & user_specified_cutpoint_closed ):
		moving_res_( moving_res_list[1] ),
		moving_res_list_( moving_res_list ),
		cutpoint_open_( cutpoint_open ),
		is_cutpoint_( full_sequence.size(), false ),
		job_parameters_( new StepWiseRNA_JobParameters ),
		filter_user_alignment_res_( true ), //Generally want to keep this true!
		allow_chain_boundary_jump_partner_right_at_fixed_BP_( false ), //hacky, just to get the Square RNA working..Nov 6, 2010
		allow_fixed_res_at_moving_res_( false ), //hacky, just to get the Hermann Duplex RNA working..Nov 15, 2010
		simple_append_map_( false ),
		skip_complicated_stuff_( false ),
		force_fold_tree_( false ),
		force_user_defined_jumps_( false ),
		force_internal_( false ),
		assert_jump_point_in_fixed_res_( true )
  {
		output_title_text( "Enter StepWiseRNA_JobParametersSetup::constructor", TR.Debug );

		for ( Size n = 1; n <= cutpoint_open_.size();   n++ ) {
			is_cutpoint_( cutpoint_open_[ n ] ) = true;
			if ( cutpoint_open_[ n ] == user_specified_cutpoint_closed ) utility_exit_with_message( "Position cannot be both cutpoint_open and user_specified_cutpoint_closed" );
		}
		job_parameters_->set_cutpoint_open_list( cutpoint_open_ );

		utility::vector1< Size > cutpoint_closed_list;
		if ( user_specified_cutpoint_closed > 0 ) {
			is_cutpoint_( user_specified_cutpoint_closed ) = true;
			cutpoint_closed_list.push_back( user_specified_cutpoint_closed );
		}
		job_parameters_->set_cutpoint_closed_list( cutpoint_closed_list );

		///////////////////////////////////////////////////////
		output_seq_num_list( "input_res = ", input_res, TR.Debug, 30 );
		output_seq_num_list( "input_res2 = ", input_res2, TR.Debug, 30 );
		output_seq_num_list( "moving_res_list = ", moving_res_list_, TR.Debug, 30 );
		TR.Debug << "moving_res = " << moving_res_ << std::endl;

		utility::vector1< utility::vector1< Size > > input_res_vectors;
		input_res_vectors.push_back( input_res );
		input_res_vectors.push_back( input_res2 );
		job_parameters_->set_input_res_vectors( input_res_vectors );
		job_parameters_->set_full_sequence( full_sequence );
		job_parameters_->set_moving_res( moving_res_ );
		figure_out_working_sequence_and_mapping(); //Initialize this here since full_to_sub and is_working_res is needed by many setting functions

		output_title_text( "Exit StepWiseRNA_JobParametersSetup::constructor", TR.Debug );
	}

  //////////////////////////////////////////////////////////////////////////
  //destructor
  StepWiseRNA_JobParametersSetup::~StepWiseRNA_JobParametersSetup()
  {}

	//////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParametersSetup::apply() {

		using namespace core::pose;
		using namespace core::kinematics;
		using namespace core::chemical;

		output_title_text( "Enter StepWiseRNA_JobParameter_Setup::apply", TR.Debug );

		setup_additional_cutpoint_closed();
		setup_floating_base_jump_to_anchor(); // special case -- connect moving res & anchor res by jump.
		figure_out_chain_boundaries();
		if ( !force_fold_tree_ ){
			figure_out_cuts();
			figure_out_jump_partners();
			setup_fold_tree();
		}

		////////////////Change the order that these functions are called on May 3, 2010 Parin S./////////////////////////////////////
		Size root_res( 1 );
		if ( !skip_complicated_stuff_ ){
			// Following determines which residues in the pose will keep fixed coordinates. Its complicated
			// because of the many use cases...
			InternalWorkingResidueParameter const internal_params = figure_out_partition_definition();
			root_res = reroot_fold_tree( internal_params.fake_working_moving_suite );
			//need the final rerooted fold_tree, WARNING: this function resets the working_moving_res_list annd working_moving_res for the internal case...
			//Warning this leaves is_working_res NOT updated...
			figure_out_prepend_internal( root_res, internal_params );
			figure_out_gap_size_and_five_prime_chain_break_res(); //Need partition definition to be initialized...
			figure_out_is_prepend_map(); //Need fold_tree and fixed_res to be initialized
		} else {
			root_res = reroot_fold_tree_simple();
			filter_user_alignment_res_ = false;
		}

		figure_out_working_alignment( root_res );

		output_title_text( "Exit StepWiseRNA_JobParameter_Setup::apply", TR.Debug );
	}

///////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< core::Size >
	StepWiseRNA_JobParametersSetup::get_user_input_alignment_res_list( core::Size const root_res ){

		//////////////////////////////////////////////////////////////////////////
		// 1. Why is alignment_res a string vector? I see that there is a tokenize to
		//   parse dashes, but why? -- rhiju
		//
		//    Alignment_res is a string vector to allow user to specify multiple possible
		//    alignments. Allow us to distinguish the follow cases:
		//
		//        i. 1-16 8-9 (align poses using res 1 and 16 OR res 8 and 9)
		//        ii. 1-6 7-12 (align poses using res 1 and 6 OR res 7 and 12)
		//        iii. 1-6-7-12 (align poses using res 1, 8, 9 and 16)
		//
		//                C01-C02-A03-A04-A05-A06-G07-G08
		//                 |   |   .   .   .   .   |   |
		//                C16-C16-G14-G13-G12-G11-C10-C09
		//
		//    If multiple possible alignments are specified, the choice of the best
		//    alignment is determined by the get_user_input_alignment_res_list() function
		//    inside StepWiseRNA_JobParametersSetup class.  -- parin (2013)

		output_title_text( "Enter get_user_input_alignment_res_list", TR.Debug );

		utility::vector1< core::Size > working_best_alignment;

		if ( working_best_alignment.size() != 0 ){ //make sure the initial size of the vector is zero..this is kinda dumb...but safe.. May 4, 2010 Parin S.
			utility_exit_with_message( " empty_working_best_alignment.size() != 0" );
		}

		utility::vector1< std::string > const & alignment_res_string_list = alignment_res_string_list_;

		if ( alignment_res_string_list.size() == 0 ){
			TR.Debug << "WARNING: alignment_res_string_list.size() == 0. EARLY RETURN AN EMPTY LIST" << std::endl;
			return working_best_alignment;
		}

		if ( fixed_res_.size() == 0 ){
			utility_exit_with_message( "need to called set_fixed_res before calling set_alignment_res" );
		}

		output_boolean( "filter_user_alignment_res_ = ", filter_user_alignment_res_, TR.Debug ); TR.Debug << std::endl;
		output_seq_num_list( "fixed_res = ", fixed_res_, TR.Debug, 30 );

		for ( Size n = 1; n <= alignment_res_string_list.size(); n++ ){

			utility::vector1< std::string > alignments_res_string = tokenize( alignment_res_string_list[n], "-" );
			utility::vector1< core::Size > alignment_res;
			utility::vector1< core::Size > working_alignment;

			for ( Size ii = 1; ii <= alignments_res_string.size(); ii++ ){
				alignment_res.push_back( string_to_int( alignments_res_string[ii] ) );
			}

			if ( filter_user_alignment_res_ ){

				utility::vector1< core::Size > actual_alignment_res;

				for ( Size ii = 1; ii <= alignment_res.size(); ii++ ){
					Size seq_num = alignment_res[ii];
					if ( fixed_res_.has_value( seq_num ) ) actual_alignment_res.push_back( seq_num );
				}

				working_alignment = apply_full_to_sub_mapping( actual_alignment_res, job_parameters_ );

				ObjexxFCL::FArray1D < bool > const & partition_definition = job_parameters_->partition_definition();
				bool contain_non_root_partition_seq_num = false;
				for ( Size ii = 1; ii < working_alignment.size(); ii++ ){
					if ( partition_definition( working_alignment[ii] ) != partition_definition( root_res ) ) contain_non_root_partition_seq_num = true;
				}

				if ( contain_non_root_partition_seq_num == true ) continue;

			} else{
				working_alignment = apply_full_to_sub_mapping( alignment_res, job_parameters_ );
			}

			//utility::vector1< core::Size > working_alignment=apply_full_to_sub_mapping(alignment_res, job_parameters_);

			if ( working_alignment.size() > working_best_alignment.size() ) working_best_alignment = working_alignment;
		}


		output_seq_num_list( "best_working_align = ", apply_sub_to_full_mapping( working_best_alignment, job_parameters_ ), TR.Debug, 30 );

//		if(alignment_res_string_list.size()>0 && working_best_alignment.size()==0){ Not compatible with build from scratch mode! Sept 08, 2010
//			utility_exit_with_message( "User supplied alignment_res_string_list but working_best_alignment.size()==0!!" );
//		}


		output_title_text( "", TR.Debug );
		return working_best_alignment;

	}

	//////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParametersSetup::figure_out_working_alignment( Size const & root_res ){
		utility::vector1< core::Size > const working_best_alignment = get_user_input_alignment_res_list( root_res );
		if ( working_best_alignment.size() > 0 ){
			job_parameters_->set_working_best_alignment( working_best_alignment );
		} else{
			figure_out_best_working_alignment();
		}
	}

	//////////////////////////////////////////////////////////////////////////
	//no user specified alignment, have to figure this out internally
	void
	StepWiseRNA_JobParametersSetup::figure_out_best_working_alignment(){

		output_title_text( "Enter StepWiseRNA_JobParametersSetup::figure_out_best_working_alignment. no user specified alignment, have to figure this out internally", TR.Debug );

		core::kinematics::FoldTree const & fold_tree = job_parameters_->fold_tree();
		utility::vector1< core::Size > working_fixed_res( job_parameters_->working_fixed_res() );
		ObjexxFCL::FArray1D < bool > const & partition_definition = job_parameters_->partition_definition();
		std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );
		utility::vector1< core::Size > const & is_working_res( job_parameters_->is_working_res() );

		bool const root_partition = partition_definition( fold_tree.root() );

		utility::vector1< core::Size > working_alignment;


		for ( Size n = 1; n <= working_fixed_res.size(); n++ ){
			Size seq_num = working_fixed_res[n];
			if ( partition_definition( seq_num ) == root_partition ) {
				working_alignment.push_back( seq_num );
			}
		}



 		//Special case for building loop outward
		if ( working_alignment.size() == 0 ){

			//March 17, 2012: RE-comment out since that are possible cases where all the working_fixed_res resides in the 'non-root' partition.
			//For example at the last building-step of of build-outward mode + idealized helix.
			//if( working_fixed_res.size()!=0) utility_exit_with_message( "working_alignment.size()==0 but fixed_res_.size()!=0 !!" );

			TR.Debug << " special case of building loop outward...no fixed element. " << std::endl;

			//Basically include every working_res as an alignment_res.....INCLUDING THE MOVING RES!
			//Don't worry about virtual res, this is taken care of by the function create_alignment_id_map(). Parin Apr 23, 2010
			for ( Size full_seq_num = 1; full_seq_num <= is_working_res.size(); full_seq_num++ ){
				if ( is_working_res[full_seq_num] == false ) continue;

				Size const seq_num = full_to_sub[full_seq_num];
				working_alignment.push_back( seq_num );
			}
		}

		output_seq_num_list( "working_fixed_res = ", working_fixed_res, TR.Debug, 30 );
		output_seq_num_list( "working_alignment = ", working_alignment, TR.Debug, 30 );
		output_seq_num_list( "working_alignment = ", working_alignment, TR.Debug, 30 );

		job_parameters_->set_working_best_alignment( working_alignment );

		if ( job_parameters_->working_best_alignment().size() == 0 ){
			utility_exit_with_message( "job_parameters_->working_best_alignment().size() == 0" );
		}

		output_title_text( "", TR.Debug );

	}
	/////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParametersSetup::figure_out_working_sequence_and_mapping(){ //working_sequence is now updated internally by JP.

		output_title_text( "Enter StepWiseRNA_JobParametersSetup::figure_out_working_sequence_and_mapping", TR.Debug );

		// There are up to two input poses. Need to merge their input residues and figure out chain boundaries.
		std::string const & full_sequence = job_parameters_->full_sequence();

		Size const nres( full_sequence.size() );
		utility::vector1< utility::vector1< Size > > const & input_res_vectors = job_parameters_->input_res_vectors();

		///////////////////////////////////////////////////////////////////////////////////////////////
		//A vector to indicate if res exist in the working pose and if it belong to input_pose 1 or input_pose 2 or working_res_list.
		utility::vector1< core::Size > is_working_res( nres, 0 );

		for ( Size i = 1; i <= input_res_vectors.size(); i++ ){
			for ( Size n = 1; n <= input_res_vectors[i].size(); n++ ){
				is_working_res[ input_res_vectors[i][ n ] ] = i; //WARNING, might not work with case of overlapping residue!
			}
		}

		for ( Size i = 1; i <= moving_res_list_.size(); i++ ){
			is_working_res[ moving_res_list_[i] ] = 999;
		}

		job_parameters_->set_is_working_res( is_working_res );
		///////////////////////////////////////////////////////////////////////////////////////////////
		std::map< core::Size, core::Size > full_to_sub;

		full_to_sub.clear();
		Size count( 0 );
		for ( Size i = 1; i <= full_sequence.size(); i++ ) {
			if ( is_working_res[ i ] ) {
				count++;
				full_to_sub[ i ] = count;
			}
		}

		job_parameters_->set_full_to_sub( full_to_sub );  //res_map
		/////////////////////////////////////////////////////////////////////////////////////////////////////


		//Size const working_moving_res = full_to_sub[ moving_res_ ];
		//job_parameters_->set_working_moving_res( working_moving_res );

		utility::vector1< core::Size > working_moving_res_list;
		for ( Size i = 1; i <= moving_res_list_.size(); i++ ){
			working_moving_res_list.push_back( full_to_sub[ moving_res_list_[i] ] );
		}

		job_parameters_->set_working_moving_res_list(  working_moving_res_list ); //This automatically set working_moving_res as well (Jan 15, 2011)
		////////////////////////////////////////////////////////////////////////////////////////////////////

		output_title_text( "", TR.Debug );

	}

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// prepend map is used solely for RMSD calculations, I think -- rhiju.
	bool
	StepWiseRNA_JobParametersSetup::figure_out_is_residue_prepend( Size const seq_num ) const {

		output_title_text("Enter StepWiseRNA_JobParametersSetup::figure_out_is_residue_prepend", TR.Debug );

		using namespace ObjexxFCL;

		utility::vector1< core::Size > const & is_working_res( job_parameters_->is_working_res() );
		Size const total_residues( job_parameters_->full_sequence().size() );
		core::kinematics::FoldTree const & fold_tree = job_parameters_->fold_tree();
		std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );
		utility::vector1< core::Size > working_fixed_res( job_parameters_->working_fixed_res() );

		if ( !is_working_res[seq_num] ){
			output_bool_list( "is_working_res = ", is_working_res, TR.Debug, 30 );
			utility_exit_with_message( "In is_residue_prependable function. Input seq_num: " + string_of( seq_num ) + " is not a element of moving_res_list_" );
		}

		// special case where floating base is specified to have a jump to an anchor, with no intervening bulge.
		if ( job_parameters_->floating_base_anchor_res() > 0 && seq_num == job_parameters_->moving_res() ) {
			return ( job_parameters_->floating_base_anchor_res() > job_parameters_->moving_res() );
		}

		//if can prepend, then must be able to find a fix res as decrement without encounter chainbreak
		Size cur_seq_num = seq_num;

		//Check if can prepend
		while ( true ){
			if ( cur_seq_num < 1 ) break;
			if ( seq_num == total_residues ) break; //Cannot prepend last residue
			if ( !is_working_res[cur_seq_num] ) break;
			if ( cur_seq_num != seq_num && is_cutpoint_( cur_seq_num - 1 ) ) break;
			if ( allow_chain_boundary_jump_partner_right_at_fixed_BP_ ) {
				for ( Size i = 1; i <= jump_point_pair_list_.size(); i++ ){
					if ( cur_seq_num == jump_point_pair_list_[i].first ) return true; //is_prepend
					if ( cur_seq_num == jump_point_pair_list_[i].second ) return true; //is_prepend
				}
			}
			if ( fixed_res_.has_value( cur_seq_num ) ) return true; //is_prepend

			//For build loop outward case, in this case make every residue append except for the 1st residue in the chain..
			//check that the fold_tree is simple..could also use the is_simple_tree() function that don't really understand how this function works
			if ( fold_tree.num_cutpoint() == 0 ){

				//This line was commented out since it conflicted with Fang's electron density code.
				//if( working_fixed_res.size()!=0) utility_exit_with_message( "fold_tree.num_cutpoint()==0 but fixed_res_.size()!=0 !!" );

				if ( full_to_sub[cur_seq_num] == 1 ) return true; //prepend

			}

			cur_seq_num++;
		}

		//Check if can append
		cur_seq_num = seq_num;

		while ( true ){
			if ( cur_seq_num > total_residues ) break;
			if ( seq_num == 1 ) break; //Cannot append first residue
			if ( !is_working_res[cur_seq_num] ) break;
			if ( is_cutpoint_( cur_seq_num ) && cur_seq_num != seq_num ) break;
			if ( allow_chain_boundary_jump_partner_right_at_fixed_BP_ ){ //Hacky Nov 12, 2010
				for ( Size i = 1; i <= jump_point_pair_list_.size(); i++ ){
					if ( cur_seq_num == jump_point_pair_list_[i].first ) return false; //is_append
					if ( cur_seq_num == jump_point_pair_list_[i].second ) return false; //is_append
				}
			}
			if ( fixed_res_.has_value( cur_seq_num ) ) return false; //is_append

			//For build loop outward case, in this case make every residue append except for the 1st residue in the chain..
			//check that the fold_tree is simple..could also use the is_simple_tree() function that don't really understand how this function works
			if ( fold_tree.num_cutpoint() == 0 ){

				// commented out. It should be OK to append residue to a pose with a simple fold tree! rhiju, july 2013.
				//				if( working_fixed_res.size()!=0) utility_exit_with_message( "fold_tree.num_cutpoint()==0 but fixed_res_.size()!=0 !!" );

				if ( full_to_sub[cur_seq_num] != 1 ) return false; //append

			}
			cur_seq_num--;
		}

		//Error, if reach this point of the function
		TR.Error << "Error: figure_out_is_residue_prepend, residue seq_num: " << seq_num << std::endl;
		TR.Error << "Cannot attach residue by either prepending and appending!" << std::endl;
		exit ( 0 );

//		output_title_text("", TR.Debug );
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Choose to use map data structure due to its ability to access an arbitrary element of a sequence in equal time (random access)
	void
	StepWiseRNA_JobParametersSetup::figure_out_is_prepend_map(){
		//
		//6. What is this is_prepend_map used for?
		//   Why is its setup connected to calc_rms_res?
		//                 -- rhiju
		//
		//  void
		//  StepWiseRNA_JobParametersSetup::figure_out_is_prepend_map():
		//
		//   is_prepend_map is used during RMSD calculation [see rmsd_over_copy_dofs()
		//   function in StepWiseRNA_Util.cc]
		//
		//   It determine whether the atoms in the 5' or 3' phosphate should be
		//   included in each residue's RMSD calculation.
		//
		//   For example if prepend_map[2] is true, then it specifies that res 2 was
		//   built by prepend and its 3' phosphate should be included in the RMSD calc.
		//   Likewise if prepend_map[5] is false, then it specifies that res 5 was
		//   built by append and its 5' phosphate should be included in the RMSD calc.
		//
		//   As you might guess, this is important only for edge residues (the one right
		//   at the 5' and 3' end of the chains).
		//
		//   There is probably better ways to handle this problem: for example include
		//   both 3' and 5' phosphate atoms of residues in RMSD calc. by default and
		//   then skip atoms that are virtualized.
		//                     -- parin (2013)

		output_title_text( "Enter StepWiseRNA_JobParametersSetup::figure_out_is_residue_prepend_map", TR.Debug );

		utility::vector1< core::Size > const & calc_rms_res = job_parameters_->calc_rms_res();
		output_seq_num_list( "calc_rms_res = ", calc_rms_res, TR.Debug, 30 );
		std::map< core::Size, bool > is_prepend_map;

		if ( simple_append_map_ ){

			TR.Debug << "WARNING using simple_append_map!" << std::endl;
			//Only initial is_residue_prepend for rebuild residues;
			for ( Size n = 1; n <= calc_rms_res.size(); n++ ){
				Size const seq_num = calc_rms_res[n];
				is_prepend_map[seq_num] = false;
			}

		} else{

			//Only initial is_residue_prepend for rebuild residues;
			for ( Size n = 1; n <= calc_rms_res.size(); n++ ){
				Size const seq_num = calc_rms_res[n];
				is_prepend_map[seq_num] = figure_out_is_residue_prepend( seq_num );
			}
		}

		job_parameters_->set_is_prepend_map( is_prepend_map );

		output_title_text( "", TR.Debug );
	}
	//////////////////////////////////////////////////////////////////////////

	utility::vector1< core::Size >
	StepWiseRNA_JobParametersSetup::get_previously_closed_cutpoint_from_imported_silent_file() const{

		using namespace ObjexxFCL;

		utility::vector1< core::Size > previously_closed_cutpoint;

		utility::vector1< utility::vector1< Size > > const & input_res_vectors = job_parameters_->input_res_vectors();

		if ( silent_files_in_.size() != input_tags_.size() ) utility_exit_with_message( "silent_files_in.size() != input_tags.size()" );

		for ( Size silent_file_num = 1; silent_file_num <= silent_files_in_.size(); silent_file_num++ ){
			pose::Pose import_pose;
			import_pose_from_silent_file( import_pose, silent_files_in_[ silent_file_num ], input_tags_[silent_file_num] );

			 utility::vector1< Size > const & input_to_full_res_map = input_res_vectors[silent_file_num]; //input_res is map from input_seq_num to full_pose_seq_num

			for ( Size seq_num = 1; seq_num <= import_pose.total_residue(); seq_num++ ){
				if ( import_pose.residue( seq_num ).has_variant_type( chemical::CUTPOINT_LOWER ) ){
					if ( import_pose.residue( seq_num + 1 ).has_variant_type( chemical::CUTPOINT_UPPER ) == false ){
						utility_exit_with_message( "seq_num " + string_of( seq_num ) + " is a CUTPOINT_LOWER but seq_num " + string_of( seq_num + 1 ) + " is not a cutpoint CUTPOINT_UPPER??" );
					}

					previously_closed_cutpoint.push_back( input_to_full_res_map[seq_num] );

				}
			}
		}

		output_seq_num_list( "previously_closed_cutpoint_ = ", previously_closed_cutpoint, TR.Debug, 30 );

		return previously_closed_cutpoint;
	}

	//////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParametersSetup::setup_additional_cutpoint_closed(){
		//The setup_additional_cutpoint_closed() checks and adds any
		//       additional cutpoint_closed.
		//
		//       It is necessary to add an additional cutpoint_closed in the
		//       following situation:
		//
		//                 C1-N2-N3-N4-N5-N6-N7-G8
		//
		//       Suppose C1 and G8 is a jump-point pair. We want to keep the
		//       relative coordinates of C1 and G8 fixed.
		//
		//       The function setup_additional_cutpoint_closed() will check to see if
		//       any of the residues N2 to N7 are non_fixed_res (movable during
		//       minimization). If there is a non_fixed_res, then there must be a
		//       cutpoint somewhere in between C1 and G8 in order to keep the
		//       relative coordinates of C1 and C8 fixed.
		//
		//       The function setup_additional_cutpoint_closed() check if any such
		//       cutpoint was already specified by the user. If no such point is
		//       specified, then it adds an any additional cutpoint_closed.
		//
		//       A real use case is finding the cutpoint-closed position of the
		//       tetraloop, described in iii. in setup_rna_job_parameters in stepwise.rna_main.cc

		output_title_text( "Enter StepWiseRNA_JobParametersSetup::setup_additional_cutpoint_closed", TR.Debug );
		utility::vector1< core::Size > const & is_working_res( job_parameters_->is_working_res() );
		utility::vector1< Size > const & cutpoint_closed_list = job_parameters_->cutpoint_closed_list();
		Size const nres( job_parameters_->full_sequence().size() );

		utility::vector1< Size > new_cutpoint_closed_list = cutpoint_closed_list;
		utility::vector1< core::Size > non_fixed_res;
		added_cutpoint_closed_.clear();

		for ( Size seq_num = 1; seq_num <= nres; seq_num++ ){
			if ( fixed_res_.has_value( seq_num ) == false ){
				non_fixed_res.push_back( seq_num );
			}
		}

		//0==false, 1==true
		output_size_list( "is_working_res = ", is_working_res, TR.Debug, 30 );
		output_bool_list( "is_working_res = ", is_working_res, TR.Debug, 30 );
		output_seq_num_list( "fixed res = ", fixed_res_, TR.Debug, 30 );
		output_seq_num_list( "non_fixed res = ", non_fixed_res, TR.Debug, 30 );
		if ( new_cutpoint_closed_list.size() > 0 ){
			TR.Debug << "user_specified_cutpoint_closed = " << new_cutpoint_closed_list[1] << std::endl;
		}


		utility::vector1< Size > jump_res_list;
		for ( Size i = 1; i <= jump_point_pair_list_.size(); i++ ) {
			jump_res_list.push_back( jump_point_pair_list_[i].first );
			jump_res_list.push_back( jump_point_pair_list_[i].second );
		}

		for ( Size n = 1; n <= non_fixed_res.size(); n++ ){

			if ( !is_working_res[ non_fixed_res[n] ] ) continue;

			/////////////////////////////////////////////////////////////////////////////////////
			bool free_boundary = false;

			Size five_prime_boundary = 0;
			for ( int seq_num = non_fixed_res[n]; seq_num >= 0; seq_num-- ){ //find the 5' base_pair_list res boundary
				if ( seq_num == 0 ){ //no 5' fixed res boundary
					free_boundary = true;
					break;
				}

				if ( jump_res_list.has_value( seq_num ) && is_working_res[seq_num] ){
					five_prime_boundary = seq_num;
					break;
				}
			}

			Size three_prime_boundary = 0;
			for ( Size seq_num = non_fixed_res[n]; seq_num <= nres + 1; seq_num++ ){ //find the 3' base_pair_list res boundary

				if ( seq_num == nres + 1 ){ //no 3' fixed res boundary
					free_boundary = true;
					break;
				}

				if ( jump_res_list.has_value( seq_num ) && is_working_res[seq_num] ){
					three_prime_boundary = seq_num;
					break;
				}
			}

			if ( free_boundary == true ) continue;
			/////////////////////////////////////////////////////////////////////////////////////
			bool boundary_is_a_pair = false;
			for ( Size ii = 1; ii <= jump_point_pair_list_.size(); ii++ ){
				std::pair < core::Size, core::Size > jump_point_pair = jump_point_pair_list_[ii];
				if ( five_prime_boundary == jump_point_pair.first && three_prime_boundary == jump_point_pair.second ){
					boundary_is_a_pair = true;
					break;
				}
			}
			if ( boundary_is_a_pair == false ) continue;
			/////////////////////////////////////////////////////////////////////////////////////

			bool found_cutpoint_or_moving_res_or_nonworking_res = false;
			for ( Size seq_num = five_prime_boundary; seq_num <= three_prime_boundary - 1; seq_num++ ){
				if ( is_cutpoint_( seq_num ) || seq_num == moving_res_ || ( !is_working_res[ seq_num ] ) ){
					found_cutpoint_or_moving_res_or_nonworking_res = true;
					break;
				}
			}

			//Aug 22, 2010...include the nonworking_res check
			//This fix is to get the FOUR edge mode working..basically if there a non_working_res..the cutpoint will be automatically set up by itself..no need for additional cutpoint close


			if ( found_cutpoint_or_moving_res_or_nonworking_res == true ) continue;  //Actually check for moving res should be enough...
			///////////////////////////////////////////////////////////////////////////////////////

			//Really need to add a cutpoint to keep the fixed_base_pair fixed wrt to each other!!!

			//Preferable to put cutpoint_closed at a place where the chain is closed at a previous folding step. Sept 1, 2010

			Size cutpoint = 0;

			utility::vector1< core::Size > const previously_closed_cutpoint = get_previously_closed_cutpoint_from_imported_silent_file();

			for ( Size seq_num = five_prime_boundary; seq_num <= three_prime_boundary - 1; seq_num++ ){

				if ( previously_closed_cutpoint.has_value( seq_num ) && is_working_res[seq_num] ){
					cutpoint = seq_num;
					break;
				}
			}

			if ( cutpoint == 0 ){ //OK couldn't find the a valid previously closed cutpoint..so will put arbitrary put it 2 residues from three_prime_fixed_res

				Size three_prime_fixed_res = 0;
				for ( Size seq_num = non_fixed_res[n]; seq_num <= nres + 1; seq_num++ ){ //find the 3' base_pair_list res boundary

					if ( seq_num > three_prime_boundary ){ //no 3' fixed res boundary
						utility_exit_with_message( "three_prime_fixed_res > three_prime_boundary of fixed base pair ?? " );
					}

					if ( fixed_res_.has_value( seq_num ) && is_working_res[seq_num] ){
						three_prime_fixed_res = seq_num;
						break;
					}
				}

				cutpoint = three_prime_fixed_res - 2; //-1 is 3' most non-fixed res, use -2 so that the chain_break torsion will be minimize
			}

			if ( is_working_res[cutpoint] == false ) utility_exit_with_message( "cutpoint is not a working_res!" );

			new_cutpoint_closed_list.push_back( cutpoint );
			added_cutpoint_closed_.push_back( cutpoint );
			is_cutpoint_( cutpoint ) = true;

		}

		TR.Debug << "added_cutpoint_closed_ = ";
		for ( Size n = 1; n <= added_cutpoint_closed_.size(); n++ ){
			TR.Debug << " " << added_cutpoint_closed_[n];
		}
		TR.Debug << std::endl;

		output_seq_num_list( "cutpoint_closed_ = ", new_cutpoint_closed_list, TR.Debug, 30 );
		output_title_text( "", TR.Debug );

		job_parameters_->set_cutpoint_closed_list( new_cutpoint_closed_list );

	}

	///////////////////////////////////////////////////////////////////////////////////

	void
	StepWiseRNA_JobParametersSetup::figure_out_chain_boundaries(){

		output_title_text( "Enter StepWiseRNA_JobParametersSetup::figure_out_chain_boundaries", TR.Debug );

		Size const nres( job_parameters_->full_sequence().size() );

		utility::vector1< core::Size > const & is_working_res( job_parameters_->is_working_res() );
		utility::vector1< std::pair < core::Size, core::Size > > chain_boundaries;

		Size start_chain( 0 );
		Size end_chain( 0 );
		Size n( 0 );
		for ( Size pos = 1; pos <= nres; pos++ ) {
			if ( !is_working_res[ pos ] ) continue;
			n++;
			if ( n == 1 ) start_chain = pos;
			if ( n > 1 && ( pos > end_chain + 1 || is_cutpoint_( end_chain ) ) ) { //pos > end_chain + 1 happen if !is_working_res[ pos ]==true. i.e a gap in the chain
				TR.Debug << "start_chain = " << start_chain << " end_chain = " << end_chain << std::endl;
				chain_boundaries.push_back( std::make_pair( start_chain, end_chain ) ); //The last chain...
				start_chain = pos;
			}
			end_chain = pos;
		}

		// For now, need to have at least one chain defined in the input!
		if ( ( start_chain > 0 ) == false ) utility_exit_with_message( "start_chain > 0" );

		chain_boundaries.push_back( std::make_pair( start_chain, end_chain ) );
		TR.Debug << "start_chain = " << start_chain << " end_chain = " << end_chain << std::endl;

		job_parameters_->set_chain_boundaries( chain_boundaries );

		output_title_text( "", TR.Debug );
	}

	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParametersSetup::figure_out_jump_partners() {

		output_title_text( "StepWiseRNA_JobParametersSetup::figure_out_jump_partners", TR.Debug );
		output_boolean( "allow_chain_boundary_jump_partner_right_at_fixed_BP_ = ", allow_chain_boundary_jump_partner_right_at_fixed_BP_, TR.Debug ); TR.Debug << std::endl;

		jump_partners_.clear();

		utility::vector1< std::pair < core::Size, core::Size > > const & chain_boundaries(	 job_parameters_->chain_boundaries() );
		std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );
		utility::vector1< core::Size > const & is_working_res = job_parameters_->is_working_res();
		utility::vector1< Size > const & cutpoint_closed_list = job_parameters_->cutpoint_closed_list();

		Size const num_chains( chain_boundaries.size() );

		//  rhiju put this in out of desperation (aug. 2013) -- need to fix this (perhaps use what is in stepwise pose setup!?)
		if ( force_user_defined_jumps_ ){
			for ( Size n = 1; n <= jump_point_pair_list_.size(); n++ ){
				Size const jump_partner1 = jump_point_pair_list_[n].first;
				if ( !is_working_res[ jump_partner1 ] ) continue;
				Size const jump_partner2 = jump_point_pair_list_[n].second;
				if ( !is_working_res[ jump_partner2 ] ) continue;
				jump_partners_.push_back( std::make_pair( full_to_sub[ jump_partner1 ], full_to_sub[ jump_point_pair_list_[n].second ] ) );
			}
			TR.Debug << "Number of jumps: "  << jump_partners_.size() << std::endl;
			TR.Debug << "Number of chains: " << num_chains << std::endl;
			runtime_assert( jump_partners_.size() == (num_chains - 1) );
			return;
		}

		TR.Debug << "num_chains = " << num_chains << std::endl;
		for ( Size n = 1 ; n < num_chains; n++ ) {

			Size const jump_partner1 = chain_boundaries[n].second;
			Size const jump_partner2 = chain_boundaries[n + 1].first;

			////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			bool pass_consecutive_res_jump_partner_test = true;

			// this condition should never hold -- indeed, there's a utility_exit below.
 			if ( !is_working_res[ jump_partner1 ] || !is_working_res[ jump_partner2 ] ){
				 pass_consecutive_res_jump_partner_test = false;
			}
			if ( moving_res_ == jump_partner2 || moving_res_ == jump_partner1 ){
				 pass_consecutive_res_jump_partner_test = false;
			}

			if ( jump_partner1 + 1 == jump_partner2 &&
					 ! cutpoint_closed_list.has_value( jump_partner1 ) && //not a cutpoint closed
					 ( ( fixed_res_.has_value( jump_partner1 ) && fixed_res_.has_value( jump_partner2 ) )	 ||
						 ( allow_chain_boundary_jump_partner_right_at_fixed_BP_ && pass_consecutive_res_jump_partner_test ) ) //Nov 2010, get Square RNA to work
					 ){

 				if ( !is_working_res[ jump_partner1 ] || !is_working_res[ jump_partner2 ] ) utility_exit_with_message( "jump_partner should be working res!" );
				if ( moving_res_ == jump_partner2 || moving_res_ == jump_partner1 ) utility_exit_with_message( "jump_partner should not be moving_res_!" );
				TR.Debug << std::setw( 80 ) << "jump_partner1 + 1 = jump_partner2 case: jump_partner1 = " << jump_partner1 << "  jump_partner2 = " << jump_partner2 << std::endl;
				jump_partners_.push_back( std::make_pair( full_to_sub[ jump_partner1 ], full_to_sub[ jump_partner2 ] ) );

			} else{

				std::pair < core::Size, core::Size > fixed_base_pair;

				bool found_jump_point_pair = false;
				for ( Size i = 1; i <= jump_point_pair_list_.size(); i++ ){ //Try an exterior pair
					if ( jump_point_pair_list_[i].first < jump_partner1 && jump_partner2 < jump_point_pair_list_[i].second ){
						if ( is_working_res[jump_point_pair_list_[i].first] && is_working_res[jump_point_pair_list_[i].second] ){
							fixed_base_pair = jump_point_pair_list_[i];
							found_jump_point_pair = true;
							break;
						}
					}
				}

				///////////Nov 6, 2010...hacky mainly to get the square RNA working.../////////////
				if ( !found_jump_point_pair ){
					if ( allow_chain_boundary_jump_partner_right_at_fixed_BP_ ){
						for ( Size i = 1; i <= jump_point_pair_list_.size(); i++ ){ //Try an exterior pair
							if ( jump_point_pair_list_[i].first <= jump_partner1 && jump_partner2 <= jump_point_pair_list_[i].second ){
								if ( is_working_res[jump_point_pair_list_[i].first] && is_working_res[jump_point_pair_list_[i].second] ){
									fixed_base_pair = jump_point_pair_list_[i];
									TR.Debug << "warning allow_chain_boundary_jump_partner_right_at_fixed_BP_ = true" << std::endl;
									TR.Debug << "jump_partner1( local ) = " << fixed_base_pair.first << "  jump_partner2( local ) = " << fixed_base_pair.second << " is right at a fixed_BP " << std::endl;
									found_jump_point_pair = true;
									break;
								}
							}
						}
					}
				}
				///////////////////////////////////////////////////////////////////////////////////


				runtime_assert ( found_jump_point_pair );

				TR.Debug << std::setw( 80 ) << "exterior_fixed_base_pair_case : jump_partner1 = " << fixed_base_pair.first << "  jump_partner2 = " << fixed_base_pair.second << std::endl;
				jump_partners_.push_back( std::make_pair( full_to_sub[ fixed_base_pair.first ], full_to_sub[ fixed_base_pair.second ] ) );
			}
		}

		output_title_text( "", TR.Debug );
	}

	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParametersSetup::figure_out_cuts() {

		cuts_.clear();

		utility::vector1< std::pair < core::Size, core::Size > > const & chain_boundaries(	 job_parameters_->chain_boundaries() );
		std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );
		Size const num_chains = chain_boundaries.size();

		for ( Size n = 1; n < num_chains; n++ ) {
			cuts_.push_back( full_to_sub[ chain_boundaries[n].second ] );
		}

	}

	///////////////////////////////////////////////////////////////////////////////////
	// special case for floating base -- moving base is attached via a 'jump' to an anchor residue, specified by user.
	void
	StepWiseRNA_JobParametersSetup::setup_floating_base_jump_to_anchor(){
		if ( !job_parameters_->floating_base() ) return;
		Size const moving_res = job_parameters_->moving_res();
		Size const anchor_res = job_parameters_->floating_base_anchor_res();
		if ( anchor_res == 0 ) return;

		utility::vector1< core::Size > const & is_working_res = job_parameters_->is_working_res();
		runtime_assert( is_working_res[ moving_res ] );
		runtime_assert( is_working_res[ anchor_res ] );

		std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );
		jump_point_pair_list_.push_back( std::make_pair( moving_res, anchor_res ) );
		jump_partners_.push_back( std::make_pair( full_to_sub[moving_res], full_to_sub[anchor_res] ) );
		TR.Debug << "MOVING_RES -- ANCHOR_RES in working numbering: " << full_to_sub[ moving_res ] << " -- " << full_to_sub[ anchor_res ] << std::endl;

		// Note -- following is not necessary.
		// put a cutpoint right in between these residues. Note that they should not be adjacent in the full
		// sequence, but let's assume that they are together in the working pose, leading to a cutpoint open.
		// runtime_assert(  std::abs( int(moving_res) - int(anchor_res) ) > 1 ); // full numbering.
		// runtime_assert(  std::abs( int(full_to_sub[ moving_res ]) - int(full_to_sub[ anchor_res ]) ) == 1 );
		Size cutpoint = ( moving_res < anchor_res ) ? full_to_sub[ moving_res ] : full_to_sub[ anchor_res ];
		cuts_.push_back( cutpoint );
		is_cutpoint_( cutpoint ) = true;

	}


	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParametersSetup::setup_fold_tree(){

		output_title_text( "Enter StepWiseRNA_JobParametersSetup::setup_fold_tree", TR.Debug );

		std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );
		runtime_assert( cuts_.size() == jump_partners_.size() );

		std::string const working_sequence = job_parameters_->working_sequence();
		Size const nres( working_sequence.size() );
		Size const num_cuts( cuts_.size() );

		ObjexxFCL::FArray2D < int > jump_point( 2, num_cuts, 0 );
		ObjexxFCL::FArray1D < int > cuts( num_cuts, 0 );

		for ( Size i = 1; i <= num_cuts; i++ ) {
			jump_point( 1, i ) = std::min( jump_partners_[i].first, jump_partners_[i].second );
			jump_point( 2, i ) = std::max( jump_partners_[i].first, jump_partners_[i].second );
			cuts( i ) = cuts_[ i ];
			TR.Debug << " JUMP POINT: " << jump_point( 1, i ) << " " << jump_point( 2, i ) << " " << cuts( i ) << std::endl;
		}

		Size const root_res = ( full_to_sub[ moving_res_ ] == 1 ) ? nres : 1;

		core::kinematics::FoldTree fold_tree( nres );
		fold_tree.tree_from_jumps_and_cuts( nres, num_cuts, jump_point, cuts, root_res ); //order of element in jump_point and cuts does not have to match. Jan 29, 2010 Parin S.
		Size num_cutpoint = fold_tree.num_cutpoint();

		for ( Size i = 1; i <= num_cutpoint; i++ ) {
			Size const k = fold_tree.upstream_jump_residue( i );
			Size const m = fold_tree.downstream_jump_residue( i );

			char upstream_res = working_sequence[k - 1];
			char downstream_res = working_sequence[m - 1];

			//Base atoms...
			std::string upstream_jump_atom;
			std::string downstream_jump_atom;

			if ( upstream_res == 'u' || upstream_res == 'c' ){
				upstream_jump_atom = " C2 ";
			} else if ( upstream_res == 'a' || upstream_res == 'g' ){
				upstream_jump_atom = " C4 ";
			} else{
				utility_exit_with_message( "Invalid upstream_res!!" );
			}

			if ( downstream_res == 'u' || downstream_res == 'c' ){
				downstream_jump_atom = " C2 ";
			} else if ( downstream_res == 'a' || downstream_res == 'g' ){
				downstream_jump_atom = " C4 ";
			} else{
				utility_exit_with_message( "Invalid downstream_res!!" );
			}

			TR.Debug << "upstream_res = " << k << upstream_res << " upstream_jump_atom = " << upstream_jump_atom;
			TR.Debug << " downstream_res = " << k << downstream_res << " downstream_jump_atom = " << downstream_jump_atom << std::endl;
			fold_tree.set_jump_atoms( i, downstream_jump_atom, upstream_jump_atom );
		}
		/////////////////////////////////////////////////////////////////////////////////////////////////

		output_fold_tree_info( fold_tree, "StepWiseRNA_JobParametersSetup::setup_fold_tree()", TR.Debug );
		job_parameters_->set_fold_tree( fold_tree );
		output_title_text( "", TR.Debug );

	}


	///////////////////////////////////////////////////////////////////////////////////

	core::Size
	StepWiseRNA_JobParametersSetup::input_struct_definition( core::Size const working_seq_num ){

		utility::vector1< utility::vector1< Size > > const & input_res_vectors = job_parameters_->input_res_vectors();

		if ( input_res_vectors.size() != 2 ) utility_exit_with_message( "is_internal case but input_res_vectors.size() != 2" );
		std::map< core::Size, core::Size > & sub_to_full( job_parameters_->sub_to_full() );
		if ( input_res_vectors[1].has_value( sub_to_full[working_seq_num] ) ) return 1;
		if ( input_res_vectors[2].has_value( sub_to_full[working_seq_num] ) ) return 2;

		TR << "working_seq_num = " << working_seq_num << " full_seq_num = " << sub_to_full[working_seq_num] << std::endl;
		utility_exit_with_message( "seq_num is not part of either input_res_vectors[1] or input_res_vectors[2]" );
		return 0; //Just to prevent compiler warning...

	}

	///////////////////////////////////////////////////////////////////////////////////
	InternalWorkingResidueParameter
	StepWiseRNA_JobParametersSetup::figure_out_partition_definition(){

		/////////////////////////////////////////////////////////////////////////////////////////////
		// trick to figure out which residues are upstream vs. downstream of the moving suite --
		// there's already a fold_tree function to do this, but it partitions based on a JUMP.
		//  So put in a fake jump between the moving_residue and the neighbor it is connected to.

		Size const nres( job_parameters_->working_sequence().size() );

		ObjexxFCL::FArray1D_bool partition_definition( nres, false );

		core::kinematics::FoldTree const & fold_tree = job_parameters_->fold_tree();
//		Size const & moving_suite( job_parameters_->working_moving_suite() );

		/////////May 3, 2010/////////////////////////////////////////////////////////////////////
		// Does the order of working_res matter? If so, which residue is 'closest' to attachment point?
		utility::vector1< core::Size > const & working_moving_res_list( job_parameters_->working_moving_res_list() );

		Size const working_moving_res = working_moving_res_list[1]; //The one furthest away from the existing structure
		Size const first_working_moving_res = working_moving_res_list[working_moving_res_list.size()]; //The one attached to the existing structure

		core::Size fake_working_moving_suite;
		InternalWorkingResidueParameter internal_working_res_params;

		if ( ( working_moving_res == 1 || fold_tree.is_cutpoint( working_moving_res - 1 ) ) && !force_internal_ ) { //prepend

			fake_working_moving_suite = first_working_moving_res ;
			if ( job_parameters_->floating_base() ) fake_working_moving_suite = working_moving_res;

		} else if ( ( fold_tree.is_cutpoint( working_moving_res ) || working_moving_res == nres )  && !force_internal_ ){

			fake_working_moving_suite = first_working_moving_res - 1;
			if ( job_parameters_->floating_base() ) fake_working_moving_suite = working_moving_res - 1;

		} else { //internal case...problematic/complicated case....

			// I don't remember what this coding stands for. When would we supply the working_moving_res 'backwards'? -- rhiju
			//			So the order of the working_moving_res should only matters if there is
			//       more than 1 moving_res (e.g. dinucleotide move)
			//
			//       One possible use case that I can think of is dinucleotide + chain-chain.
			//       (although this move is current not used in default runs)
			//
			//               N1-N2-C3-C4-N5-N6.
			//
			//       Suppose we start with :
			//
			//               N1-N2       N5-N6
			//
			//       Dinucleotide move can be used to close the chain in two ways:
			//
			//           i. Append [4, 3] : Treat res 4 as sampling base and 3 as bulge.
			//           ii. prepend [3, 4] Treat res 3 as sampling base and 4 as bulge.
			//
			//           Note that the first residue in the list is always the one that
			//           is sampled (in this case by 'floating base' sampling). The
			//           second and subsequent residues are bulge res.
			//
			//       However, it should be possible to figure out which direction to built
			//       using just the cutpoint-closed position in this case. So perhaps ordering
			//       the working_moving_res is a redundant feature that should be removed?
			//
			// I thought of another (more realistic) use case for question 5:
			//
			//
			//               N6-N7
			//               N5-N8
			//               C4
			//               C3
			//               N2 N9
			//               N1-N10
			//
			//       In this case, working_moving_res [3,4] and [4,3] gives different
			//       behavior:
			//
			//           i. Append [4, 3] : Treat res 4 as sampling base and 3 as bulge.
			//           ii. prepend [3, 4] : Treat res 3 as sampling base and 4 as bulge.
			//
			//       Although, this move is also currently not used in default SWA runs.
			//
			bool const can_append  = check_can_append( working_moving_res_list ); //[14, 13, 12]
			bool const can_prepend = check_can_prepend( working_moving_res_list ); //[12, 13, 14]

			output_boolean( "can_prepend = ", can_prepend, TR.Debug ); output_boolean( " can_append = ", can_append, TR.Debug ); TR.Debug << std::endl;

			if ( !can_prepend && !can_append ){
				output_seq_num_list( "working_moving_res_list:", working_moving_res_list, TR.Debug );
				utility_exit_with_message( "Cannot prepend or append residue in working_moving_res_list" );
			}

			//Ok first find the two possible positions to put the actual_working_res.
			Size possible_working_res_1 = 0;
			Size possible_working_res_2 = 0;
			Size found_possible_working_res = 0;

			//Check if this is a suite right between two previously built moving_elements.
			if ( working_moving_res_list.size() == 1 ){

				if ( working_moving_res < nres ){
					if ( input_struct_definition( working_moving_res ) != input_struct_definition( working_moving_res + 1 ) ) {
						possible_working_res_1 = working_moving_res;
						possible_working_res_2 = working_moving_res + 1;
						found_possible_working_res++;
					}
				}

				if ( working_moving_res > 1 ){
					if ( input_struct_definition( working_moving_res ) != input_struct_definition( working_moving_res - 1 ) ) {
						possible_working_res_1 = working_moving_res - 1;
						possible_working_res_2 = working_moving_res;
						found_possible_working_res++;
					}
				}

			} else{

				if ( can_prepend ){ //[11,12,13]
					possible_working_res_1 = working_moving_res;
					possible_working_res_2 = working_moving_res_list[working_moving_res_list.size()] + 1;
					found_possible_working_res++;
				}
				if ( can_append ){ //[13,12,11]
					possible_working_res_1 = working_moving_res_list[working_moving_res_list.size()] - 1;
					possible_working_res_2 = working_moving_res;
					found_possible_working_res++;
				}
			}

			if ( found_possible_working_res != 1 ){
				TR.Debug << "found_possible_working_res = " << found_possible_working_res << std::endl;
				utility_exit_with_message( "found_possible_working_res != 1. Cannot figure out use case!" );
			}

			// RHIJU is disabling this temporarily for stepwise.monte carlo stuff -- need to cleanup this entire
			//    script soon with Parin. 8 July, 2012
			if ( !skip_complicated_stuff_ && input_struct_definition( possible_working_res_1 ) == input_struct_definition( possible_working_res_2 ) ){
				utility_exit_with_message( "input_struct_definition( possible_working_res_1 ) == input_struct_definition( possible_working_res_2 )" );
			}

			fake_working_moving_suite = possible_working_res_1; //This is kinda adhoc...can choose any res between possible_working_res_1 and (possible_working_res_2-1)

			// rhiju -- found boundary case -- this whole setup stuff needs to be cleaned up.
			TR.Debug << "POSSIBLE_WORKING_RES1 " << possible_working_res_1 << " POSSIBLE_WORKING_RES2 " << possible_working_res_2 << std::endl;
			while ( fold_tree.is_cutpoint( fake_working_moving_suite ) && fake_working_moving_suite < possible_working_res_2 ) fake_working_moving_suite++;

			internal_working_res_params.possible_working_res_1 = possible_working_res_1;
			internal_working_res_params.possible_working_res_2 = possible_working_res_2;
		}

		internal_working_res_params.fake_working_moving_suite = fake_working_moving_suite;


		///////////////////////////////////////////////////////////////////////////////
		fold_tree.partition_by_residue( fake_working_moving_suite, partition_definition );
		job_parameters_->set_partition_definition( partition_definition ); //this is a useful decomposition.

		return internal_working_res_params;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Figure out a good root residue -- ideally in a fixed segment!
	core::Size
	StepWiseRNA_JobParametersSetup::reroot_fold_tree_simple(){

		core::kinematics::FoldTree const & fold_tree = job_parameters_->fold_tree();

		utility::vector1< core::Size > working_fixed_res( job_parameters_->working_fixed_res() );

		for ( Size i = 1; i <= working_fixed_res.size(); i++ ){
			Size const n = working_fixed_res[i];
			if ( fold_tree.possible_root( n )  ){
				core::kinematics::FoldTree rerooted_fold_tree = fold_tree;
				bool reorder_went_OK = rerooted_fold_tree.reorder( n );
				if ( reorder_went_OK ) {
					job_parameters_->set_fold_tree( rerooted_fold_tree );
					return n;
				}
			}
		}

		return 1;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Figure out a good root residue -- which partition of the pose has the fewest residues to move around?
	core::Size
	StepWiseRNA_JobParametersSetup::reroot_fold_tree( core::Size const fake_working_moving_suite ){

		output_title_text( "Enter StepWiseRNA_JobParametersSetup::reroot_fold_tree", TR.Debug );

		ObjexxFCL::FArray1D < bool > const & partition_definition = job_parameters_->partition_definition();
		core::kinematics::FoldTree const & fold_tree = job_parameters_->fold_tree();
		Size const nres = job_parameters_->working_sequence().size();

		std::map< core::Size, core::Size > & sub_to_full( job_parameters_->sub_to_full() );

		Size num_partition_0( 0 ), num_partition_1( 0 );
		Size possible_new_root_residue_in_partition_0( 0 ), possible_new_root_residue_in_partition_1( 0 ), root_res( 0 );

		for ( Size n = 1; n <= nres; n++ ) {
			if ( partition_definition( n ) ) {
				num_partition_1 += 1;
				if ( fold_tree.possible_root( n ) ) possible_new_root_residue_in_partition_1 = n;
			} else {
				num_partition_0 += 1;
				if ( fold_tree.possible_root( n ) ) possible_new_root_residue_in_partition_0 = n;
			}
		}

		///////////////////////////////////////////////////////////////////////////////////////////////////////
		//If moving_res = 1 or nres, try putting root at fix res (second best choice)
		utility::vector1< core::Size > working_fixed_res( job_parameters_->working_fixed_res() );
		for ( Size i = 1; i <= working_fixed_res.size(); i++ ){
			Size const seq_num = working_fixed_res[i];
			if ( partition_definition( seq_num ) && fold_tree.possible_root( seq_num ) ){
				possible_new_root_residue_in_partition_1 = seq_num;
				break;
			}
		}

		for ( Size i = 1; i <= working_fixed_res.size(); i++ ){
			Size const seq_num = working_fixed_res[i];
			if ( !partition_definition( seq_num ) &&  fold_tree.possible_root( seq_num ) ){
				possible_new_root_residue_in_partition_0 = seq_num;
				break;
			}
		}

		//If moving_res=1 or nres, try putting root at terminus res (best choice)
		utility::vector1< core::Size > working_terminal_res( job_parameters_->working_terminal_res() );

		for ( Size i = 1; i <= working_terminal_res.size(); i++ ){
			Size const seq_num = working_terminal_res[i];
				if ( partition_definition( seq_num ) &&  fold_tree.possible_root( seq_num ) ){
				possible_new_root_residue_in_partition_1 = seq_num;
				break;
			}
		}

		for ( Size i = 1; i <= working_terminal_res.size(); i++ ){
			Size const seq_num = working_terminal_res[i];
			if ( !partition_definition( seq_num ) &&  fold_tree.possible_root( seq_num ) ){
				possible_new_root_residue_in_partition_0 = seq_num;
				break;
			}
		}
		///////////////////////////////////////////////////////////////////////////////////////////////////////
		runtime_assert( num_partition_0 > 0 );
		runtime_assert( num_partition_1 > 0 );

		Size const moving_res( job_parameters_->working_moving_res() );
		if ( num_partition_1 == num_partition_0 ){
			if ( partition_definition( moving_res ) == 0 ){ //put root_res in partition 1 away from the moving_res.
				if ( partition_definition( 1 ) && partition_definition( nres ) ) {
					root_res = 1;
				} else {
					root_res = possible_new_root_residue_in_partition_1;
				}
			} else{ //put root_res in partition 0 away from the moving_res.
				if ( !partition_definition( 1 ) && !partition_definition( nres ) ) {
					root_res = 1;
				} else {
					root_res = possible_new_root_residue_in_partition_0;
				}
			}
		} else if ( num_partition_1 > num_partition_0 ){
			// best to put the root in partition 1 -- it is bigger, and will stay anchored.
			if ( partition_definition( 1 ) && partition_definition( nres ) ) {
				root_res = 1;
			} else {
				root_res = possible_new_root_residue_in_partition_1;
			}
		} else {
			if ( !partition_definition( 1 ) && !partition_definition( nres ) ) {
				root_res = 1;
			} else {
				root_res = possible_new_root_residue_in_partition_0;
			}
		}

		TR.Debug << "Num. res in partition 0:  " << num_partition_0 << ".   Num. res in partition 1: " << num_partition_1 << std::endl;
		TR.Debug << "Moving_res full seq_num   = "   << sub_to_full[moving_res] << ", partition "   << partition_definition( moving_res ) << std::endl;
		TR.Debug << "fake_working_moving_suite = "   << sub_to_full[fake_working_moving_suite ] << ", partition "   << partition_definition( fake_working_moving_suite ) << std::endl;
		TR.Debug << "New root res full seq_num = " << sub_to_full[root_res]   << ", partition " << partition_definition( root_res ) << std::endl;
		runtime_assert ( root_res > 0 );

		TR.Debug << "Before reroot: " << fold_tree << std::endl;
		core::kinematics::FoldTree rerooted_fold_tree = fold_tree;
		bool reorder_went_OK = rerooted_fold_tree.reorder( root_res );
		if ( !reorder_went_OK ) utility_exit_with_message( "!reorder_went_OK" );

		if ( force_fold_tree_ ) rerooted_fold_tree = job_parameters_->fold_tree();;

		// moving positions
		utility::vector1< Size > working_moving_partition_pos;
		bool const root_partition = partition_definition( rerooted_fold_tree.root() );

		for ( Size seq_num = 1; seq_num <= nres; seq_num++ ){
			if ( partition_definition( seq_num ) != root_partition ) 	working_moving_partition_pos.push_back( seq_num );
		}

		job_parameters_->set_working_moving_partition_pos( working_moving_partition_pos );
		job_parameters_->set_fold_tree( rerooted_fold_tree );

		output_fold_tree_info( fold_tree, "StepWiseRNA_JobParametersSetup::reroot_fold_tree()", TR.Debug );

		output_title_text( "", TR.Debug );

		return root_res;
	}

	///////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParametersSetup::figure_out_gap_size_and_five_prime_chain_break_res(){

		utility::vector1< std::pair < core::Size, core::Size > > const & chain_boundaries(	 job_parameters_->chain_boundaries() );

 		Size gap_size = GAP_SIZE_DUMMY; // junk value... totally "free" end.
 		job_parameters_->set_gap_size( gap_size /*DUMMY*/ );
 		job_parameters_->set_five_prime_chain_break_res( 0 );

		/////////////////////////////////////////////////////////////////////////////////////////////
		// Need to look for a chainbreak whose ends actually will move relative to each other if
		// we change degrees of freedom in the "moving residues".
		ObjexxFCL::FArray1D < bool > const & partition_definition = job_parameters_->partition_definition();
		Size const moving_res( job_parameters_->moving_res() );
		Size const anchor_res( job_parameters_->floating_base_anchor_res() ); // 0 if not floating_base by jump.

		/////////////////////////////////////////////////////////////////////////////////////////////
 		std::map< core::Size, core::Size > & full_to_sub( job_parameters_->full_to_sub() );
		Size found_moving_gap( 0 );
		Size const num_chains = chain_boundaries.size();
		for ( Size n = 1; n < num_chains; n++ ) {
			Size chain_end = chain_boundaries[ n ].second;
			Size next_chain_start = chain_boundaries[ n + 1 ].first;
			if ( chain_end == anchor_res && next_chain_start == moving_res ) continue; // attachment by jump
			if ( chain_end == moving_res && next_chain_start == anchor_res ) continue; // attachment by jump
			if ( partition_definition( full_to_sub[ chain_end ] ) !=
					 partition_definition( full_to_sub[ next_chain_start ] ) ){

				bool found_cutpoint_open( false );
				for ( Size i = 1; i <= cutpoint_open_.size(); i++ ){
					if ( cutpoint_open_[i] >= chain_end && cutpoint_open_[i] < next_chain_start ) {
						found_cutpoint_open = true;
						break;
					}
				}
				TR.Debug << "CHECK CHAIN FOR GAP SIZE " << chain_end << " " << next_chain_start << " " << found_cutpoint_open << std::endl;
				if ( found_cutpoint_open ) continue;

				bool found_added_cutpoint_closed_( false );
				for ( Size i = 1; i <= added_cutpoint_closed_.size(); i++ ){
					if ( added_cutpoint_closed_[i] >= chain_end && added_cutpoint_closed_[i] < next_chain_start ) {
						found_added_cutpoint_closed_ = true;
						break;
					}
				}
				if ( found_added_cutpoint_closed_ ) continue;

				job_parameters_->set_gap_size(  next_chain_start - chain_end - 1 );
				job_parameters_->set_five_prime_chain_break_res( full_to_sub[ chain_end ] );
				found_moving_gap++;

			}
		}

		if ( found_moving_gap > 1 ){
			utility_exit_with_message( "Had trouble figure out which gap corresponds to the user specified cutpoint_closed! Try to renumber input poses sequentially." );
		}

	}
	///////////////////////////////////////////////////////////////////////////////////////////////////
	StepWiseRNA_JobParametersOP & StepWiseRNA_JobParametersSetup::job_parameters(){
		return job_parameters_;
	}

	/////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParametersSetup::figure_out_prepend_internal( core::Size const root_res, InternalWorkingResidueParameter const & internal_params ){

		output_title_text( "Enter StepWiseRNA_JobParametersSetup::figure_out_prepend_internal", TR.Debug );
		std::map< core::Size, core::Size > & sub_to_full( job_parameters_->sub_to_full() );

		ObjexxFCL::FArray1D < bool > const & partition_definition = job_parameters_->partition_definition();
		Size const working_moving_res( job_parameters_->working_moving_res() ); //hard copy, since this is set later in the function
		utility::vector1< core::Size > const working_moving_res_list( job_parameters_->working_moving_res_list() ); //hard copy, since this is set later in the function
		core::kinematics::FoldTree const & fold_tree = job_parameters_->fold_tree();
		Size const nres = job_parameters_->working_sequence().size();
		Size const fake_working_moving_suite = internal_params.fake_working_moving_suite;

		bool is_prepend( true ), is_internal( false );

		// some defaults
		job_parameters_->set_is_prepend( is_prepend );
		job_parameters_->set_is_internal( is_internal );
		if ( skip_complicated_stuff_ ) return;

		if ( job_parameters_->floating_base_anchor_res() > 0 ) { // residue is connected by jump ){
			is_prepend = ( job_parameters_->moving_res() < job_parameters_->floating_base_anchor_res() );
			is_internal = false;
		}	else if ( ( working_moving_res == 1 ||
									fold_tree.is_cutpoint( working_moving_res - 1 ) ) && !force_internal_ ){
			is_prepend = true;
			is_internal = false;
		} else if ( ( fold_tree.is_cutpoint( working_moving_res ) ||
									working_moving_res == nres ) && !force_internal_ ){
			//Are we sure that the 3' res of the chain will always be a cutpoint? No, that's why I
			// put in the force_internal_ override. There are cases where moving residue has
			// a covalent connection to next residue and then a jump to previous (which could come
			// from a skip bulge move)
			is_prepend = false;
			is_internal = false;
		} else { //The problematic internal case....this case is quite complicated....

			TR.Debug << "is_internal case " << std::endl;
			is_internal = true;

			bool const can_append = check_can_append( working_moving_res_list ); //[14, 13, 12]
			bool const can_prepend = check_can_prepend( working_moving_res_list ); //[12, 13, 14]

			Size const possible_working_res_1 = internal_params.possible_working_res_1; //lower
			Size const possible_working_res_2 = internal_params.possible_working_res_2; //upper

			TR.Debug << can_append << " " << can_prepend << " " << possible_working_res_1 << " " << possible_working_res_2 << std::endl;

			Size found_actual_working_res = 0;
			utility::vector1< core::Size > actual_working_moving_res_list;
			core::Size actual_working_moving_res;

			//OK have to put moving_res AWAY from the root res...
			if ( partition_definition( possible_working_res_1 ) != partition_definition( root_res ) ){
				// prepend case -- the parts of the pose 3' to the moving suites are fixed (containing root_res).
				found_actual_working_res++;
				actual_working_moving_res = possible_working_res_1;
				is_prepend = true;

				if ( working_moving_res_list.size() == 1 ){
					actual_working_moving_res_list.push_back( actual_working_moving_res );
				} else{
					if ( can_prepend ){
						actual_working_moving_res_list = working_moving_res_list;
					} else{
						//Ok should prepend base on partition definition, but user input an append working_moving_res_list...need to convert to prepend
						for ( Size n = working_moving_res_list.size(); n >= 1; n-- ){ //Convert [14,13,12] to [11,12,13]
							actual_working_moving_res_list.push_back( working_moving_res_list[n] - 1 );
						}

						if ( check_can_prepend( actual_working_moving_res_list ) == false ){
							output_seq_num_list( "actual_working_moving_res_list = ", actual_working_moving_res_list, TR.Debug );
							utility_exit_with_message( "actual_working_moving_res_list fails can_prepend assertion" );
						}
					}
				}
			}

			//OK have to put moving_res AWAY from the root res...
			if ( partition_definition( possible_working_res_2 ) != partition_definition( root_res ) ){
				found_actual_working_res++;
				actual_working_moving_res = possible_working_res_2;
				is_prepend = false;

				if ( working_moving_res_list.size() == 1 ){
					actual_working_moving_res_list.push_back( actual_working_moving_res );
				} else{
					if ( can_append ){
						actual_working_moving_res_list = working_moving_res_list;
					} else{
					//Ok should append base on partition definition, but user input an prepend working_moving_res_list...need to convert to append
						for ( Size n = working_moving_res_list.size(); n >= 1; n-- ){ //Convert [11,12,13] to [14,13,12]
							actual_working_moving_res_list.push_back( working_moving_res_list[n] + 1 );
						}

						if ( check_can_append( actual_working_moving_res_list ) == false ){
							output_seq_num_list( "actual_working_moving_res_list = ", actual_working_moving_res_list, TR.Debug );
							utility_exit_with_message( "actual_working_moving_res_list fails can_append assertion" );
						}
					}
				}
			}
			if ( actual_working_moving_res_list[1] != actual_working_moving_res ){
				TR.Debug << "actual_working_moving_res = " << actual_working_moving_res << std::endl;
				output_seq_num_list( "actual_working_moving_res_list = ", actual_working_moving_res_list, TR.Debug );
				utility_exit_with_message( "actual_working_moving_res_list[1] != actual_working_moving_res" );
			}

			if ( found_actual_working_res != 1 ){
				TR.Debug << "found_actual_working_res = " << found_actual_working_res << std::endl;
				utility_exit_with_message( "found_actual_working_res != 1" );
			}
//////////////////////////////////////////////////////////////

			output_seq_num_list( "actual_working_moving_res_list = ", actual_working_moving_res_list, TR.Debug );
			output_seq_num_list( "user_input_working_moving_res_list = ", working_moving_res_list, TR.Debug );

			TR.Debug << "User input moving_res full seq_num = " << sub_to_full[working_moving_res] << ", partition "   << partition_definition( working_moving_res ) << std::endl;
			TR.Debug << "Actual moving_res full seq_num = " << sub_to_full[actual_working_moving_res_list[1]] << ", partition "   << partition_definition( actual_working_moving_res_list[1] ) << std::endl;
			TR.Debug << "New root res full seq_num =      " << sub_to_full[root_res]   << ", partition " << partition_definition( root_res ) << std::endl;

			//reset working_moving_res and working_moving_res_list
			//job_parameters_->set_working_moving_res( actual_working_moving_res);
			job_parameters_->set_working_moving_res_list( actual_working_moving_res_list ); //This update working_move_res as well

			//final check!
			if ( partition_definition( root_res ) == partition_definition( actual_working_moving_res_list[1] ) ){
				utility_exit_with_message( "partition_definition( root_res ) == partition_definition( actual_working_moving_res_list[1]!!" );
			}

		}

		job_parameters_->set_is_prepend( is_prepend );
		job_parameters_->set_is_internal( is_internal );
		//Check...
		if ( !job_parameters_->is_internal() ){
			bool const should_be_prepend = ( partition_definition( fake_working_moving_suite ) != partition_definition( root_res ) );
			if ( should_be_prepend != job_parameters_->is_prepend() ) {
				TR.Debug << "sub_to_full[fake_working_moving_suite] = " << sub_to_full[fake_working_moving_suite];
				TR.Debug << " partition_definition( fake_working_moving_suite ) = " << partition_definition( fake_working_moving_suite );
				TR.Debug << " sub_to_full[root_res] = " << sub_to_full[root_res];
				TR.Debug << " partition_definition( root_res ) = " << partition_definition( root_res ) << std::endl;
				output_boolean( "should_be_prepend = ", should_be_prepend, TR.Debug );
				output_boolean( " job_parameters_->is_prepend() = ", job_parameters_->is_prepend(), TR.Debug ); TR.Debug << std::endl;
				//				utility_exit_with_message( "Possible problem with prepend/append assignment!!" );
			}
		}
	}



  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParametersSetup::set_fixed_res( utility::vector1 < core::Size > const & fixed_res ){

		output_title_text( "Enter StepWiseRNA_JobParametersSetup::set_fixed_res", TR.Debug );

		fixed_res_ = fixed_res;

		using namespace ObjexxFCL;

		output_seq_num_list( "moving_res_list = ", moving_res_list_, TR.Debug, 30 );
		output_boolean( "allow_fixed_res_at_moving_res_ = ", allow_fixed_res_at_moving_res_, TR.Debug ); TR.Debug << std::endl;

		if ( !allow_fixed_res_at_moving_res_ ){ //Nov 15, 2010 (the outer if statement)
			for ( Size n = 1; n <= moving_res_list_.size(); n++ ){
				if ( fixed_res_.has_value( moving_res_list_[n] ) ){
					TR.Debug << "moving_res " + string_of( moving_res_list_[n] ) + " should not be in fixed_res_list!" << std::endl;
					//utility_exit_with_message( "moving_res " + string_of( moving_res_list_[n] ) + " should not be in fixed_res_list!" );
				}
			}
		}

		job_parameters_->set_working_fixed_res(  apply_full_to_sub_mapping( fixed_res_, StepWiseRNA_JobParametersCOP( job_parameters_ ) ) );

		output_title_text( "", TR.Debug );

	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParametersSetup::set_terminal_res( utility::vector1 < core::Size > const & terminal_res ){

		job_parameters_->set_terminal_res( terminal_res );
		job_parameters_->set_working_terminal_res(  apply_full_to_sub_mapping( terminal_res, StepWiseRNA_JobParametersCOP( job_parameters_ ) ) );


//		output_seq_num_list("working_terminal_res= ", job_parameters_->working_terminal_res(), TR.Debug, 30 );

	}

  //////////////////////////////////////////////////////////////////////////

	void
	StepWiseRNA_JobParametersSetup::set_calc_rms_res( utility::vector1< core::Size > const & input_calc_rms_res ){

		utility::vector1< core::Size > const & is_working_res( job_parameters_->is_working_res() );
		utility::vector1< core::Size > actual_calc_rms_res;

		for ( Size n = 1; n <= input_calc_rms_res.size(); n++ ){
			Size seq_num = input_calc_rms_res[n];
			if ( !is_working_res[seq_num] ) continue;
			actual_calc_rms_res.push_back( seq_num );
		}

		job_parameters_->set_calc_rms_res( actual_calc_rms_res );

	}

  //////////////////////////////////////////////////////////////////////////
	//example input: "1-2 19-20"
	void
	StepWiseRNA_JobParametersSetup::set_jump_point_pair_list( utility::vector1< std::string > const & jump_point_pairs_string ){
		// 2. Why is jump_point_pairs a string vector? OK to allow pairs readin "1,2 3,4"
		//  (but also make backwards compatible with "1-2 3-4"?  -- rhiju
		//
		//    Jump_point_pairs is a string vector to specify integer pairs. Allow us
		//    to distinguish:
		//
		//        i.  1-16 8-9  (jump points between res 1 and 16 and between res 8 and 9)
		//        ii. 1-8  9-16 (jump points between res 1 and 8 and between res 9 and 16)
		//
		//   OK to a new PairIntegerOption type?
		//
		//    Yes, the "PairIntegerOption" seems like a good idea.

		jump_point_pair_list_.clear();

		if ( jump_point_pairs_string.size() == 0 ) return;

		if ( fixed_res_.size() == 0 ) utility_exit_with_message( "need to called set_fixed_res before calling set_jump_point_pair_list" );

		bool assert_jump_point_in_fixed_res = assert_jump_point_in_fixed_res_;
		for ( Size n = 1; n <= jump_point_pairs_string.size(); n++ ){

			if ( ( n == 1 ) && ( jump_point_pairs_string[n] == "NOT_ASSERT_IN_FIXED_RES" ) ){ //Feb 09, 2012: FIXED BUG. Used to be "and" instead of "&&"
				assert_jump_point_in_fixed_res = false;
				continue;
			}

			utility::vector1< std::string > JP_pair_string = tokenize( jump_point_pairs_string[n], "-" );
			runtime_assert( JP_pair_string.size() == 2 /* Each jump_point_pair needs to have two elements! ( e.g. 1 - 10 ) */);
			std::pair < core::Size, core::Size > jump_point_pair = std::make_pair( string_to_int( JP_pair_string[1] ), string_to_int( JP_pair_string[2] ) );
			runtime_assert ( jump_point_pair.first != jump_point_pair.second );
			if ( jump_point_pair.first > jump_point_pair.second )	jump_point_pair = std::make_pair( string_to_int( JP_pair_string[2] ), string_to_int( JP_pair_string[1] ) );

			jump_point_pair_list_.push_back( jump_point_pair );
		}

		sort_pair_list( jump_point_pair_list_ ); //sort the BP list by the 1st element

		output_pair_size( jump_point_pair_list_, "jump_point_pair_list: ", TR.Debug );

		output_boolean( "assert_jump_point_in_fixed_res = ", assert_jump_point_in_fixed_res, TR.Debug );

		//Ensure that every seq_num in jump_point_pair_list_ is a fixed_res
		if ( assert_jump_point_in_fixed_res ){
			for ( Size n = 1; n <= jump_point_pair_list_.size(); n++ ){
				runtime_assert ( fixed_res_.has_value( jump_point_pair_list_[n].first ) );
				runtime_assert ( fixed_res_.has_value( jump_point_pair_list_[n].second ) );
			}
		}
	}

  //////////////////////////////////////////////////////////////////////////

	void
	StepWiseRNA_JobParametersSetup::set_native_alignment_res( utility::vector1< Size > const & setting ){

		utility::vector1< Size > const native_alignment = setting;

		job_parameters_->set_native_alignment(  native_alignment );
		job_parameters_->set_working_native_alignment(  apply_full_to_sub_mapping( native_alignment, StepWiseRNA_JobParametersCOP( job_parameters_ ) ) );

		output_title_text( "StepWiseRNA_JobParametersSetup::set_native_alignment_res", TR.Debug );
		output_seq_num_list( "native_alignment = ", job_parameters_->native_alignment(), TR.Debug );
		output_seq_num_list( "working_native_alignment = ", job_parameters_->working_native_alignment(), TR.Debug );
		output_title_text( "", TR.Debug );

	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParametersSetup::set_global_sample_res_list( utility::vector1 < core::Size > const & setting ){

		utility::vector1< Size > const global_sample_res_list = setting;
		job_parameters_->set_global_sample_res_list(  global_sample_res_list );

	}

  //////////////////////////////////////////////////////////////////////////

	void
	StepWiseRNA_JobParametersSetup::set_force_syn_chi_res_list( utility::vector1 < core::Size > const & setting ){

		utility::vector1< Size > const force_syn_chi_res_list = setting;
		job_parameters_->set_force_syn_chi_res_list(  force_syn_chi_res_list );

	}

  //////////////////////////////////////////////////////////////////////////

	void
	StepWiseRNA_JobParametersSetup::set_force_north_sugar_list( utility::vector1 < core::Size > const & setting ){

		utility::vector1 < core::Size > const force_north_sugar_list = setting;
		job_parameters_->set_force_north_sugar_list(  force_north_sugar_list );

	}

  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParametersSetup::set_force_south_sugar_list( utility::vector1 < core::Size > const & setting ){

		utility::vector1 < core::Size > const force_south_sugar_list = setting;
		job_parameters_->set_force_south_sugar_list(  force_south_sugar_list );

	}

  //////////////////////////////////////////////////////////////////////////

	void
	StepWiseRNA_JobParametersSetup::set_protonated_H1_adenosine_list( utility::vector1 < core::Size > const & setting ){

		utility::vector1 < core::Size > const protonated_H1_adenosine_list = setting;
		job_parameters_->set_protonated_H1_adenosine_list(  protonated_H1_adenosine_list );

	}


  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParametersSetup::set_output_extra_RMSDs( bool const setting ){

		job_parameters_->set_output_extra_RMSDs( setting );

	}
  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParametersSetup::set_add_virt_res_as_root( bool const setting ){

		job_parameters_->set_add_virt_res_as_root( setting );

	}
  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParametersSetup::set_floating_base( bool const setting ){

		job_parameters_->set_floating_base( setting );

	}
  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParametersSetup::set_floating_base_anchor_res( Size const setting ){

		job_parameters_->set_floating_base_anchor_res( setting );

	}
  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParametersSetup::set_rebuild_bulge_mode( bool const setting ){

		job_parameters_->set_rebuild_bulge_mode( setting );

	}
  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParametersSetup::set_sample_both_sugar_base_rotamer( bool const setting ){

		job_parameters_->set_sample_both_sugar_base_rotamer( setting );

	}
  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParametersSetup::force_fold_tree( core::kinematics::FoldTree const & fold_tree ){

		job_parameters_->set_fold_tree( fold_tree );
		force_fold_tree_ = true;

	}
  //////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_JobParametersSetup::set_cutpoint_closed_list( utility::vector1 < core::Size > const & setting ){
		job_parameters_->set_cutpoint_closed_list( setting );
	}
  //////////////////////////////////////////////////////////////////////////

} //legacy
} //rna
} //sampling
} //stepwise
} //protocols
