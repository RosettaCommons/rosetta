// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_CombineLongLoopFilterer
/// @brief Sets up pose and job parameters for RNA stepwise building.
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_CombineLongLoopFilterer.hh>

#include <protocols/swa/rna/StepWiseRNA_Classes.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>

#include <protocols/rna/RNA_ProtocolUtil.hh>
//////////////////////////////////

#include <core/types.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/util.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/VariantType.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/constraints/Func.fwd.hh>
#include <core/scoring/constraints/FadeFunc.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.fwd.hh>
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/pose/util.tmpl.hh>
#include <basic/Tracer.hh>

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
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace core;
using core::Real;

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Add class description here!
// 
// 
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

static basic::Tracer TR( "protocols.swa.rna_stepwise_rna_CombineLongLoopFilterer" ) ;

namespace protocols {
namespace swa {
namespace rna {

  //////////////////////////////////////////////////////////////////////////
  //constructor!
	StepWiseRNA_CombineLongLoopFilterer::StepWiseRNA_CombineLongLoopFilterer( StepWiseRNA_JobParametersCOP const & job_parameters, bool const combine_helical_silent_file ):
		rsd_set_( core::chemical::ChemicalManager::get_instance()->residue_type_set( core::chemical::RNA ) ),
		job_parameters_( job_parameters ),
		verbose_( true ), //Parin Mar 22, 2010		
		filter_for_previous_contact_( false ),
		filter_for_previous_clash_( false ),
		undercount_ribose_rotamers_( false ), //July 29, 2011
		filter_for_chain_closable_( true ),
		filter_for_moving_res_contact_( true ),
		moving_res_to_base_contact_only_( true ), //used if filter_for_movign_res_contact_ is true, buggy right now since ignoring hydrogen (non-heavy) atom in the base..
		total_input_struct_pair_( 0 ),
		pass_screen_struct_pair_( 0 ),
		output_filename_( "filter_struct.txt" ), //will have to change later...perhap pass in this from the python script?
		//best_combine_score_(99999999.99),   //Feb 02, 2012; This might lead to server-test error at R47198
		//worst_combine_score_(-99999999.99), //Feb 02, 2012; This might lead to server-test error at R47198
		best_combine_score_( 999999.99 ), //Feb 02, 2012;
		worst_combine_score_(  - 999999.99 ), //Feb 02, 2012;
		//score_diff_cut_(1000000.0), //Remove all score filtering on Jan 12, 2012
		contact_dist_cutoff_(  - 1.0 ), //two atoms are considered in contact if their VDW radius edge is within 1.0 Angstrom of each other
		clash_dist_cutoff_( 0.8 ),    //two atoms are considered clash if their VDW radius edge overlap by 0.8 Angstrom. 
															 //0.8 is appropriate for VDW clash screen, although value about 1.2 would be more appropriate if we consider minimum H-bond distance
															 //See StepWiseRNA_VDW_BinScreener.cc for details.
//		num_contact_cutoff_(9), //num of contact between the two sides before discarding pose (Used to be 1 before Nov 18, 2010)
		num_contact_cutoff_( 1 ), //num of contact between the two sides before discarding pose (Used to be 1 before Nov 18, 2010)
		num_clash_cutoff_( 1 ), // num of clash between the two sides before discarding pose.
		max_pose_data_list_size_( 200 ),
		side_ONE_NUM_pose_list_( 0 ),
		side_TWO_NUM_pose_list_( 0 ),
		side_ONE_pose_list_id_( 1 ),
		side_TWO_pose_list_id_( 1 ),
		//max_decoys_(9999999999), //Feb 02, 2012; This might lead to server-test error at R47198
		max_decoys_( 999999 ), //Feb 02, 2012;
		combine_helical_silent_file_( combine_helical_silent_file ) //Hacky mode to build VC_two 3 way junction
  {

		if ( combine_helical_silent_file_ == false ){

			utility::vector1< utility::vector1< Size > > const & input_res_vectors = job_parameters_->input_res_vectors();

			if ( input_res_vectors.size() != 2 ) utility_exit_with_message( "input_res_vectors.size() != 2" );

			full_to_input_res_map_ONE_ = create_full_to_input_res_map( input_res_vectors[1] ); //right now, only use this in figure_out_appended_and_prepended_res_list function
			full_to_input_res_map_TWO_ = create_full_to_input_res_map( input_res_vectors[2] ); //right now, only use this in figure_out_appended_and_prepended_res_list function

			figure_out_appended_and_prepended_res_list();
			figure_out_last_appended_and_last_prepended_res();

			//O3I_C5I_PLUS_ONE_MAX_DIST=3.968000
			//max_centroid_to_atom_distance for atom:  C5' base RAD: 6.18959
			//max_centroid_to_atom_distance for atom:  O3' base RAD: 6.65341
			//For is for a A-nucleotides, more chi rotamers
			//Suppose that moving_res is making base-stack contact to the last SWA-built residue from the another side:
			//	Then atom-atom distance must be lesser than minus_contact_dist_cutoff_(1)+atom vanderWaal radius -> for Carbon: 	1.70

			//so 3.968000 + (6.18959 or 6.65341) + (1) + (3.4)=
			moving_res_contact_dist_cutoff_ = 3.968000 + 6.65341 + 1 + 3.4;
			TR << "moving_res_contact_dist_cutoff_ = " << moving_res_contact_dist_cutoff_ << std::endl;
			TR << "contact_dist_cutoff_ = " << contact_dist_cutoff_ << std::endl;
			TR << "clash_dist_cutoff_ = " << clash_dist_cutoff_ << std::endl;
			TR << "num_contact_cutoff_ = " << num_contact_cutoff_ << std::endl;
			TR << "num_clash_cutoff_ = " << num_clash_cutoff_ << std::endl;

		}

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  //destructor
	StepWiseRNA_CombineLongLoopFilterer::~StepWiseRNA_CombineLongLoopFilterer()
  {}
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	void
	StepWiseRNA_CombineLongLoopFilterer::figure_out_appended_and_prepended_res_list(){

		//OK first find that residues that are common between two pose...there are user inputted residues

		utility::vector1< utility::vector1< Size > > const & input_res_vectors = job_parameters_->input_res_vectors();


		utility::vector1< core::Size > common_res_list;

		for ( Size n = 1; n <= input_res_vectors[1].size(); n++ ){

			Size const seq_num = input_res_vectors[1][n];

			if ( input_res_vectors[2].has_value( seq_num ) ){
				common_res_list.push_back( seq_num );
			}

		}		

		utility::vector1< core::Size > full_pose_appended_res_list;
		input_pose_ONE_appended_res_list_.clear();

		for ( Size n = 1; n <= input_res_vectors[1].size(); n++ ){
			Size const seq_num = input_res_vectors[1][n];
			if ( common_res_list.has_value( seq_num ) == false ){
			 
				full_pose_appended_res_list.push_back( seq_num );
				input_pose_ONE_appended_res_list_.push_back( full_to_input_res_map_ONE_.find( seq_num )->second );

			}
		}

		utility::vector1< core::Size > full_pose_prepended_res_list;
		input_pose_TWO_prepended_res_list_.clear();


		for ( Size n = 1; n <= input_res_vectors[2].size(); n++ ){
			Size const seq_num = input_res_vectors[2][n];
			if ( common_res_list.has_value( seq_num ) == false ){

				full_pose_prepended_res_list.push_back( seq_num );
				input_pose_TWO_prepended_res_list_.push_back( full_to_input_res_map_TWO_.find( seq_num )->second );

			}
		}

		Output_seq_num_list( "full_pose_appended_res:", full_pose_appended_res_list, TR ); 
		Output_seq_num_list( "full_pose_prepended_res:", full_pose_prepended_res_list, TR ); 

		Output_seq_num_list( "input_ONE_appended_res:", input_pose_ONE_appended_res_list_, TR ); 
		Output_seq_num_list( "input_TWO_prepended_res:", input_pose_TWO_prepended_res_list_, TR ); 

		
	}


  /////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_CombineLongLoopFilterer::figure_out_last_appended_and_last_prepended_res(){

		utility::vector1< utility::vector1< Size > > const & input_res_vectors = job_parameters_->input_res_vectors();

		//////////////////////////////////////////////////////////////////
		
		input_pose_ONE_last_appended_res_ = 0;
	
		if ( input_pose_ONE_appended_res_list_.size() == 0 ) utility_exit_with_message( "input_pose_ONE_appended_res_list_.size() == 0" );

		for ( Size n = 1; n <= input_pose_ONE_appended_res_list_.size(); n++ ){
			if ( input_pose_ONE_last_appended_res_ <= input_pose_ONE_appended_res_list_[n] ){

				input_pose_ONE_last_appended_res_ = input_pose_ONE_appended_res_list_[n];

			}
		}	

		//////////////////////////////////////////////////////////////////

		//input_pose_TWO_last_prepended_res_=9999999; //Feb 02, 2012: This might lead to server-test error at R47200 
		input_pose_TWO_last_prepended_res_ = 999999;  //Feb 02, 2012

		if ( input_pose_TWO_prepended_res_list_.size() == 0 ) utility_exit_with_message( "input_pose_TWO_prepended_res_list_.size() == 0" );

		for ( Size n = 1; n <= input_pose_TWO_prepended_res_list_.size(); n++ ){
			if ( input_pose_TWO_last_prepended_res_ >= input_pose_TWO_prepended_res_list_[n] ){

				input_pose_TWO_last_prepended_res_ = input_pose_TWO_prepended_res_list_[n];

			}
		}			

		//////////////////////////////////////////////////////////////////

		Size const full_last_appended_res = input_res_vectors[1][input_pose_ONE_last_appended_res_];
		Size const full_last_prepended_res = input_res_vectors[2][input_pose_TWO_last_prepended_res_];


		TR << "full_pose_last_appended_res_ = " << full_last_appended_res << " input_pose_ONE_last_appended_res_ = " << input_pose_ONE_last_appended_res_;
		TR << " full_pose_last_prepended_res_ = " << full_last_prepended_res << " input_pose_TWO_last_prepended_res_ = " << input_pose_TWO_last_prepended_res_ << std::endl;
	

/*
		using namespace ObjexxFCL;

		bool const Is_prepend(  job_parameters_->Is_prepend() ); 

		std::map< core::Size, core::Size > const & sub_to_full = job_parameters_->const_sub_to_full(); //make these const
	
		Size const moving_res(  job_parameters_->working_moving_res() ); // Might not corresponds to user input.
		Size const reference_res( job_parameters_->working_reference_res() ); //the last static_residues that this attach to the moving residues
	


		Size working_last_prepended_res, working_last_appended_res;

		if ( Is_prepend ){
			working_last_prepended_res = reference_res;
			working_last_appended_res = moving_res - 1; 
		} else{
			working_last_prepended_res = reference_res;
			working_last_appended_res = moving_res + 1;
		}
	

		if ( sub_to_full.count( working_last_prepended_res ) == 0 ) utility_exit_with_message( "sub_to_full.count( working_last_prepended_res ) == 0, working_last_prepended_res = " + string_of( working_last_prepended_res ) );
		if ( sub_to_full.count( working_last_appended_res ) == 0 ) utility_exit_with_message( "sub_to_full.count( working_last_appended_res ) == 0, working_last_appended_res = " + string_of( working_last_appended_res ) );


		Size const full_last_appended_res = sub_to_full.find( working_last_prepended_res )->second;
		Size const full_last_prepended_res = sub_to_full.find( working_last_appended_res )->second;

		if ( full_to_input_res_map_ONE_.count( full_last_appended_res ) == 0 ) utility_exit_with_message( "full_to_input_res_map_ONE_[full_last_appended_res].count == 0" );
		if ( full_to_input_res_map_TWO_.count( full_last_prepended_res ) == 0 ) utility_exit_with_message( "full_to_input_res_map_TWO_[full_last_prepended_res].count == 0" );

		input_pose_ONE_last_appended_res_ = full_to_input_res_map_ONE_.find( full_last_appended_res )->second;
		input_pose_TWO_last_prepended_res_ = full_to_input_res_map_TWO_.find( full_last_prepended_res )->second;
*/

	
	}

  /////////////////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< pose_data_struct2 >
	StepWiseRNA_CombineLongLoopFilterer::convert_silent_file_to_pose_data_list( core::import_pose::pose_stream::SilentFilePoseInputStreamOP & silent_file_stream, Size const pose_list_id ){

		using namespace core::pose;
		using namespace ObjexxFCL;

//		Output_title_text("Enter StepWiseRNA_CombineLongLoopFilterer::convert_silent_file_to_pose_data_list (" + path_basename(silent_file) +  ")", TR );


		utility::vector1< pose_data_struct2 > pose_data_list;



/*
		core::io::silent::SilentFileData silent_file_data;
		silent_file_data.read_file( silent_file );
		for ( core::io::silent::SilentFileData::iterator iter = silent_file_data.begin(), end = silent_file_data.end(); iter != end; ++iter ) {

		//iter->fill_pose(  *pose_op, *(rsd_set_) ); 

*/

		Size min_pose_ID = ( ( pose_list_id - 1 )*( max_pose_data_list_size_ ) ) + 1;
		Size max_pose_ID = ( ( pose_list_id )*( max_pose_data_list_size_ ) );

		Size pose_ID = 0;
		Size num_struct_in_range = 0;

		silent_file_stream->reset();

		while ( silent_file_stream->has_another_pose() ) {

			pose_ID++;

			core::io::silent::SilentStructOP const silent_struct( silent_file_stream->next_struct() );

			if ( ( pose_ID < min_pose_ID ) || ( pose_ID > max_pose_ID ) ) continue;

			num_struct_in_range++;
	
			PoseOP pose_op( new Pose );

			silent_struct->fill_pose( *pose_op, *( rsd_set_ )  ); 

			Real score( 0.0 );
			getPoseExtraScores( *pose_op, "score", score );

			std::string const & tag( silent_struct->decoy_tag() );

			if ( protocols::swa::rna::check_for_messed_up_structure( (*pose_op), tag ) ) {
				utility_exit_with_message( "tag = " + tag  + " is messed up!" );
			}

			pose_data_struct2 pose_data;
			pose_data.pose_OP = pose_op;
			pose_data.score = score;
			pose_data.tag = tag;

			pose_data_list.push_back( pose_data );


		}

		if ( num_struct_in_range == 0 ) utility_exit_with_message( "num_struct_in_range == 0! for pose_list_id = " +  string_of( pose_list_id ) + " min_pose_ID = " + string_of( min_pose_ID ) + " max_pose_ID = "  + string_of( max_pose_ID ) );

//		if(pose_data_list.size()==0) utility_exit_with_message("pose_data_list.size()==0" );

//		Output_title_text("Exit StepWiseRNA_CombineLongLoopFilterer::convert_silent_file_to_pose_data_list (" + path_basename(silent_file) +  ")", TR );


		return pose_data_list;

	}

  /////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_CombineLongLoopFilterer::align_all_pose(	utility::vector1< pose_data_struct2 > const & side_ONE_pose_data_list, 
																									  utility::vector1< pose_data_struct2 > const & side_TWO_pose_data_list ){
 		using namespace chemical;

		Output_title_text( "Enter StepWiseRNA_CombineLongLoopFilterer::align_all_pose ", TR );

		if ( side_ONE_pose_data_list.size() == 0 ) utility_exit_with_message( "side_ONE_pose_data_list.size() == 0" );
		if ( side_TWO_pose_data_list.size() == 0 ) utility_exit_with_message( "side_TWO_pose_data_list.size() == 0" );


		core::pose::Pose const alignment_pose = ( *side_ONE_pose_data_list[1].pose_OP );


		for ( Size n = 1; n <= side_ONE_pose_data_list.size(); n++ ){

			core::pose::Pose & side_ONE_pose = ( *side_ONE_pose_data_list[n].pose_OP );

			id::AtomID_Map < id::AtomID > atom_ID_map; //Align the first and last residues of the two pose (which should be the same residue)
			pose::initialize_atomid_map( atom_ID_map, side_ONE_pose, id::BOGUS_ATOM_ID ); 

			setup_suite_atom_id_map( side_ONE_pose.residue( 1 ),                             alignment_pose.residue( 1 ),                               atom_ID_map );
			setup_suite_atom_id_map( side_ONE_pose.residue( side_ONE_pose.total_residue() ), alignment_pose.residue( alignment_pose.total_residue() ),  atom_ID_map );
	
			core::scoring::superimpose_pose( side_ONE_pose, alignment_pose, atom_ID_map );

			//if(n==1) side_ONE_pose.dump_pdb( "Filterer_side_ONE_pose_"+ side_ONE_pose_data_list[n].tag +".pdb" );
			
		}

		
		for ( Size n = 1; n <= side_TWO_pose_data_list.size(); n++ ){

			core::pose::Pose & side_TWO_pose = ( *side_TWO_pose_data_list[n].pose_OP );

			id::AtomID_Map < id::AtomID > atom_ID_map; //Align the first and last residues of the two pose (which should be the same residue)
			pose::initialize_atomid_map( atom_ID_map, side_TWO_pose, id::BOGUS_ATOM_ID );
	
			setup_suite_atom_id_map( side_TWO_pose.residue( 1 ), 														alignment_pose.residue( 1 ), 															atom_ID_map );
			setup_suite_atom_id_map( side_TWO_pose.residue( side_TWO_pose.total_residue() ), alignment_pose.residue( alignment_pose.total_residue() ),  atom_ID_map );

			core::scoring::superimpose_pose( side_TWO_pose, alignment_pose, atom_ID_map );		

			//if(n==1) side_TWO_pose.dump_pdb( "Filterer_side_TWO_pose_" + side_TWO_pose_data_list[n].tag + ".pdb" );
			

		}



		Output_title_text( "Exit StepWiseRNA_CombineLongLoopFilterer::align_all_pose ", TR );


	}

  /////////////////////////////////////////////////////////////////////////////////////////////////



	bool
	StepWiseRNA_CombineLongLoopFilterer::previously_builded_res_VDW_filter( pose_data_struct2 const & side_ONE_pose_data, 
																																		 pose_data_struct2 const & side_TWO_pose_data, 
																																		 core::Real const overlap_dist_cutoff, 
																																		 core::Size const num_atom_contacts_cutoff ){

		core::pose::Pose const & side_ONE_pose = ( *side_ONE_pose_data.pose_OP );
		core::pose::Pose const & side_TWO_pose = ( *side_TWO_pose_data.pose_OP );

		for ( Size i = 1; i <= input_pose_ONE_appended_res_list_.size(); i++ ){
			for ( Size j = 1; j <= input_pose_TWO_prepended_res_list_.size(); j++ ){
	
				Size const input_pose_ONE_appended_res = input_pose_ONE_appended_res_list_[i];
				Size const input_pose_TWO_prepended_res = input_pose_TWO_prepended_res_list_[j];

				bool const residues_in_contact = Is_residues_in_contact( input_pose_ONE_appended_res, side_ONE_pose, input_pose_TWO_prepended_res, side_TWO_pose, overlap_dist_cutoff, num_atom_contacts_cutoff );

				if ( residues_in_contact ) return false;

			}
		}
	
		return true;

	}

  /////////////////////////////////////////////////////////////////////////////////////////////////


/*
		core::pose::Pose const & side_ONE_pose = ( *side_ONE_pose_data.pose_OP );
		core::pose::Pose const & side_TWO_pose = ( *side_TWO_pose_data.pose_OP );

		for ( Size i = 1; i <= input_pose_ONE_appended_res_list_.size(); i++ ){
			for ( Size j = 1; j <= input_pose_TWO_prepended_res_list_.size(); j++ ){
	
				Size const input_pose_ONE_appended_res = input_pose_ONE_appended_res_list_[i];
				Size const input_pose_TWO_prepended_res = input_pose_TWO_prepended_res_list_[j];

				core::conformation::Residue const & rsd_ONE = side_ONE_pose.residue( input_pose_ONE_appended_res );
				core::conformation::Residue const & rsd_TWO = side_TWO_pose.residue( input_pose_TWO_prepended_res );


				for ( Size at_ONE = 1; at_ONE <= rsd_ONE.natoms(); at_ONE++ ){
					for ( Size at_TWO = 1; at_TWO <= rsd_TWO.natoms(); at_TWO++ ){

						if ( rsd_ONE.atom_type( at_ONE ).name() == "VIRT" ) continue;
						if ( rsd_TWO.atom_type( at_TWO ).name() == "VIRT" ) continue;

						Real const VDW_radius_ONE = rsd_ONE.atom_type( at_ONE ).lj_radius();
						Real const VDW_radius_TWO = rsd_TWO.atom_type( at_TWO ).lj_radius();

						Real const cutoff_sum_VDW_radius = VDW_radius_ONE + VDW_radius_TWO - overlap_dist_cutoff;

						if ( cutoff_sum_VDW_radius < 0 ) utility_exit_with_message( "( VDW_radius_ONE + VDW_radius_TWO - overlap_dist_cutoff ) < 0!!" );


						Real const atom_atom_dist_squared = ( rsd_ONE.xyz( at_ONE ) - rsd_TWO.xyz( at_TWO ) ).length_squared();

						if ( atom_atom_dist_squared < ( cutoff_sum_VDW_radius*cutoff_sum_VDW_radius ) ) return false;


					}
				}

			}
		}
*/


  /////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_CombineLongLoopFilterer::previously_builded_res_contact_filter( pose_data_struct2 const & side_ONE_pose_data, pose_data_struct2 const & side_TWO_pose_data ){

		return ( previously_builded_res_VDW_filter( side_ONE_pose_data, side_TWO_pose_data, contact_dist_cutoff_, num_contact_cutoff_ ) );

	}

 	/////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_CombineLongLoopFilterer::previously_builded_res_clash_filter( pose_data_struct2 const & side_ONE_pose_data, pose_data_struct2 const & side_TWO_pose_data ){

		return ( previously_builded_res_VDW_filter( side_ONE_pose_data, side_TWO_pose_data, clash_dist_cutoff_, num_clash_cutoff_ ) );

	}
  /////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_CombineLongLoopFilterer::moving_res_contact_filter( pose_data_struct2 const & side_ONE_pose_data, pose_data_struct2 const & side_TWO_pose_data ){


		bool const Is_prepend(  job_parameters_->Is_prepend() ); 


		core::pose::Pose const & side_ONE_pose = ( *side_ONE_pose_data.pose_OP );
		core::pose::Pose const & side_TWO_pose = ( *side_TWO_pose_data.pose_OP );


		core::conformation::Residue const & last_append_rsd = side_ONE_pose.residue( input_pose_ONE_last_appended_res_ );
		core::conformation::Residue const & last_prepend_rsd = side_TWO_pose.residue( input_pose_TWO_last_prepended_res_ );

		numeric::xyzVector< core::Real > const anchor_atom_xyz = ( Is_prepend ) ? last_prepend_rsd.xyz( "C5'" ) : last_append_rsd.xyz( "O3'" );

		core::conformation::Residue const & enforce_contact_rsd = ( Is_prepend ) ? last_append_rsd : last_prepend_rsd;

		if ( enforce_contact_rsd.type().atom_name( enforce_contact_rsd.first_sidechain_atom() ) != " O2'" ){ 
			utility_exit_with_message( "enforce_contact_rsd.type().atom_name( enforce_contact_rsd.first_sidechain_atom() ) != \" O2'\" " );
		}

		Size const first_at = ( moving_res_to_base_contact_only_ ) ? ( enforce_contact_rsd.first_sidechain_atom() + 1 ) : 1;
		Size const last_at = ( moving_res_to_base_contact_only_ ) ? ( enforce_contact_rsd.nheavyatoms() ) : enforce_contact_rsd.natoms();
		//crap this ignore hydrogen atoms in the base...rsd.nheavyatoms()

		for ( Size at = first_at; at <= last_at; at++ ){

			Real const atom_atom_dist_squared = ( enforce_contact_rsd.xyz( at ) - anchor_atom_xyz ).length_squared();

			if ( atom_atom_dist_squared < ( moving_res_contact_dist_cutoff_* moving_res_contact_dist_cutoff_ ) ) return true;

		}

		return false;

	}

  /////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_CombineLongLoopFilterer::pass_all_filters( pose_data_struct2 const & side_ONE_pose_data, pose_data_struct2 const & side_TWO_pose_data ){

		using namespace ObjexxFCL;

		//if(verbose_) Output_title_text("Enter StepWiseRNA_CombineLongLoopFilterer::do_some_filtering(" + side_ONE_pose_data.tag + "," +  side_TWO_pose_data.tag  +")", TR );


		//Would this slow down the code?

		Size const num_nucleotides(  job_parameters_->working_moving_res_list().size() );
		Size const previous_step_gap_size = ( job_parameters_->gap_size() ) + num_nucleotides;

		if ( previous_step_gap_size == 0 ) utility_exit_with_message( "previous_step_gap_size == 0!!" );
		///////////////////////////////////////////////////////////////////////////


		core::pose::Pose const & side_ONE_pose = ( *side_ONE_pose_data.pose_OP );
		core::pose::Pose const & side_TWO_pose = ( *side_TWO_pose_data.pose_OP );

		Real const curr_combine_score = side_ONE_pose_data.score + side_TWO_pose_data.score;

		filterer_count_.total_count++;

		//if(curr_combine_score > (best_combine_score_ + score_diff_cut_) ) return false;
 
		filterer_count_.score_cut_count++;


		if ( ( filter_for_chain_closable_ == true ) ){ 

			//OK STILL HAVE TO WRITE CODE FOR THE DINUCLEOTIDE case....
			//July 19th, 2011..This assumes that the 5' and 3' ribose is already built (not virtual!) but this might not be the case!

			if ( Check_chain_closable( side_TWO_pose.residue( input_pose_TWO_last_prepended_res_ ).xyz( "C5'" ), side_ONE_pose.residue( input_pose_ONE_last_appended_res_ ).xyz( "O3'" ), previous_step_gap_size ) == false ) return false;
			
			filterer_count_.chain_closable_screen++;

		}

		if ( filter_for_previous_contact_ == true ){
			if ( previously_builded_res_contact_filter( side_ONE_pose_data, side_TWO_pose_data ) == false ) return false;

			filterer_count_.filter_for_previous_contact++;
		}

		if ( filter_for_previous_clash_ == true ){
			if ( previously_builded_res_clash_filter( side_ONE_pose_data, side_TWO_pose_data ) == false ) return false;

			filterer_count_.filter_for_previous_clash++;
		}


		if ( best_combine_score_ > curr_combine_score ) best_combine_score_ = curr_combine_score; //best_combine_score that still PASS the screen..


		//DON'T include this screen to determine best_combine_score since it can lead to an artificially bad best_combine_score (still need to verify this!)
		if ( ( previous_step_gap_size != 1 ) && ( filter_for_moving_res_contact_ == true ) ){ 

			//INCLUDE the dinucleotide Chainbreak code since if dinucleotide...moving_res (the floating base) cannot be bulge and need to make contact.
			if ( moving_res_contact_filter( side_ONE_pose_data, side_TWO_pose_data ) == false ) return false;

			filterer_count_.filter_for_moving_res_contact++;
		}

		if ( worst_combine_score_ < curr_combine_score ) worst_combine_score_ = curr_combine_score; //worst_combine_score that still PASS the screen


		return true;

	}


  /////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_CombineLongLoopFilterer::do_some_filtering() {


		using namespace ObjexxFCL;

		Output_title_text( "Enter StepWiseRNA_CombineLongLoopFilterer::do_some_filtering for side_ONE_pose_list_id_ = " + string_of( side_ONE_pose_list_id_ ) + "/" + string_of( side_ONE_NUM_pose_list_ ) + " side_ONE_pose_list_id_ = " + string_of( side_TWO_pose_list_id_ ) + "/" + string_of( side_TWO_NUM_pose_list_ ), TR );


//		std::ofstream outfile; 
//		outfile.open(output_filename_.c_str(), std::ios_base::out | std::ios_base::app); //This does not delete existing content
		
		utility::vector1< pose_data_struct2 > const side_ONE_pose_data_list = convert_silent_file_to_pose_data_list( silent_file_stream_ONE_, side_ONE_pose_list_id_ ); //assume this is the one build loop by append from 5' side
		utility::vector1< pose_data_struct2 > const side_TWO_pose_data_list = convert_silent_file_to_pose_data_list( silent_file_stream_TWO_, side_TWO_pose_list_id_ ); //assume this is the one build loop by prepend from 3' side

		if ( side_ONE_pose_data_list.size() == 0 ) return;
		if ( side_TWO_pose_data_list.size() == 0 ) return;


		if ( combine_helical_silent_file_ == false ){ //Nov 27 2010
			align_all_pose( side_ONE_pose_data_list, side_TWO_pose_data_list );
		}

		for ( Size i = 1; i <= side_ONE_pose_data_list.size(); i++ ){
			for ( Size j = 1; j <= side_TWO_pose_data_list.size(); j++ ){

				total_input_struct_pair_++;

				pose_data_struct2 const & side_ONE_pose_data = side_ONE_pose_data_list[i];
				pose_data_struct2 const & side_TWO_pose_data = side_TWO_pose_data_list[j];

				if ( combine_helical_silent_file_ == false ){ 
					if ( pass_all_filters( side_ONE_pose_data, side_TWO_pose_data ) == false ) continue;
				} else{ //Nov 27 2010

					Real const curr_combine_score = side_ONE_pose_data.score + side_TWO_pose_data.score;

					filterer_count_.total_count++;

					//if(curr_combine_score > (best_combine_score_ + score_diff_cut_) ) continue;
			 
					filterer_count_.score_cut_count++;

					if ( best_combine_score_ > curr_combine_score ) best_combine_score_ = curr_combine_score; //best_combine_score that still PASS the screen..

				}
	
				pass_screen_struct_pair_++;

				//ok will print to file after finish debugging 
				TR << "struct pair: ( " <<  side_ONE_pose_data.tag  << ", " <<  side_TWO_pose_data.tag << " ) pass screening test. ";
				TR << pass_screen_struct_pair_ << " out of " << total_input_struct_pair_ << " passed screen so far" << std::endl;


				Combine_Tags_Info combine_tag_info;
				combine_tag_info.side_one_tag = side_ONE_pose_data.tag;
				combine_tag_info.side_two_tag = side_TWO_pose_data.tag;
				combine_tag_info.combine_score = side_ONE_pose_data.score + side_TWO_pose_data.score;

				filterered_combine_tag_info_list_.push_back( combine_tag_info );
			
			}
		}

//		outfile.flush();
//		outfile.close();

		Output_title_text( "Exit StepWiseRNA_CombineLongLoopFilterer::do_some_filtering", TR );

		 	

	}

  /////////////////////////////////////////////////////////////////////////////////////////////////

	void
	StepWiseRNA_CombineLongLoopFilterer::figure_out_NUM_pose_list(){

		silent_file_stream_ONE_->reset();
		silent_file_stream_TWO_->reset();


		Size total_pose_side_ONE = 0;

		while ( silent_file_stream_ONE_->has_another_pose() ) {
			total_pose_side_ONE++;
			core::io::silent::SilentStructOP const silent_struct( silent_file_stream_ONE_->next_struct() );
		}

		Size total_pose_side_TWO = 0;

		while ( silent_file_stream_TWO_->has_another_pose() ) {
			total_pose_side_TWO++;
			core::io::silent::SilentStructOP const silent_struct( silent_file_stream_TWO_->next_struct() );
		}


		silent_file_stream_ONE_->reset();
		silent_file_stream_TWO_->reset();


		side_ONE_NUM_pose_list_ = 0;

		while ( true ){

			if ( total_pose_side_ONE <= ( side_ONE_NUM_pose_list_*max_pose_data_list_size_ ) ) break;

			side_ONE_NUM_pose_list_++;

		}


		side_TWO_NUM_pose_list_ = 0;

		while ( true ){

			if ( total_pose_side_TWO <= ( side_TWO_NUM_pose_list_*max_pose_data_list_size_ ) ) break;

			side_TWO_NUM_pose_list_++;

		}


		TR << "max_pose_data_list_size_ = " <<  max_pose_data_list_size_ << std::endl;
		TR << "total_pose_side_ONE = " << total_pose_side_ONE << " side_ONE_NUM_pose_list_ = " << side_ONE_NUM_pose_list_ << std::endl;
		TR << "total_pose_side_TWO = " << total_pose_side_TWO << " side_TWO_NUM_pose_list_ = " << side_TWO_NUM_pose_list_ << std::endl;


	}

  /////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_CombineLongLoopFilterer::setup_tag_to_source_map(){

		silent_file_stream_ONE_->reset();
		silent_file_stream_TWO_->reset();

		tag_to_source_map_ONE_.clear();
		tag_to_source_map_TWO_.clear();


		TR << "silent_file_stream_ONE_ parent_remarks" << std::endl;
		while ( silent_file_stream_ONE_->has_another_pose() ) {
			core::io::silent::SilentStructOP const silent_struct( silent_file_stream_ONE_->next_struct() );
			//silent_struct->print_parent_remarks(TR);

			std::string const & tag( silent_struct->decoy_tag() );

			if ( tag_to_source_map_TWO_.count( tag ) != 0 ) utility_exit_with_message( "tag " + tag + " already exist in tag_to_source_map_TWO_" );

			if ( silent_struct->has_parent_remark( "SOURCE" ) == false ){
				//utility_exit_with_message("silent_struct->has_parent_remark(\"SOURCE\")==false for tag: " + tag);
				tag_to_source_map_ONE_[tag] = "MISSING";	

			} else{

				std::string const source_file = silent_struct->get_parent_remark( "SOURCE" );
				tag_to_source_map_ONE_[tag] = source_file;	

			}

		}
		TR << "--------------------------------" << std::endl;

		TR << "silent_file_stream_TWO_ parent_remarks" << std::endl;
		while ( silent_file_stream_TWO_->has_another_pose() ) {
			core::io::silent::SilentStructOP const silent_struct( silent_file_stream_TWO_->next_struct() );
			//silent_struct->print_parent_remarks(TR);

			std::string const & tag( silent_struct->decoy_tag() );

			if ( tag_to_source_map_TWO_.count( tag ) != 0 ) utility_exit_with_message( "tag " + tag + " already exist in tag_to_source_map_TWO_" );

			if ( silent_struct->has_parent_remark( "SOURCE" ) == false ){
				//utility_exit_with_message("silent_struct->has_parent_remark(\"SOURCE\")==false for tag: " + tag);
				tag_to_source_map_TWO_[tag] = "MISSING";	

			} else{

				std::string const source_file = silent_struct->get_parent_remark( "SOURCE" );
				tag_to_source_map_TWO_[tag] = source_file;	

			}



		}
		TR << "--------------------------------" << std::endl;


		silent_file_stream_ONE_->reset();
		silent_file_stream_TWO_->reset();



	}


  /////////////////////////////////////////////////////////////////////////////////////////////////

	void
	StepWiseRNA_CombineLongLoopFilterer::setup_silent_file_stream(){

		using namespace ObjexxFCL;

		if ( silent_files_in_.size() != 2 ) utility_exit_with_message( "silent_files_in_.size() != 2!, silent_files_in_.size() = " + string_of( silent_files_in_.size() ) );
	
		silent_file_stream_ONE_ = new core::import_pose::pose_stream::SilentFilePoseInputStream();
		silent_file_stream_ONE_->set_order_by_energy( true );
		silent_file_stream_ONE_->set_record_source( false ); //change to false on July 29, 2011


		utility::vector1< std::string > singleton_list_ONE;
		singleton_list_ONE.push_back( silent_files_in_[1] );
		silent_file_stream_ONE_->filenames( singleton_list_ONE ); //triggers read in of files, too.

		silent_file_stream_TWO_ = new core::import_pose::pose_stream::SilentFilePoseInputStream();
		silent_file_stream_TWO_->set_order_by_energy( true );
		silent_file_stream_TWO_->set_record_source( false ); //change to false on July 29, 2011 

		utility::vector1< std::string > singleton_list_TWO;
		singleton_list_TWO.push_back( silent_files_in_[2] );
		silent_file_stream_TWO_->filenames( singleton_list_TWO ); //triggers read in of files, too.


	}
  /////////////////////////////////////////////////////////////////////////////////////////////////

//	bool
//	StepWiseRNA_CombineLongLoopFilterer::
	bool
	score_sort_citeria( Combine_Tags_Info tag_info_1, Combine_Tags_Info tag_info_2 ){
		return ( tag_info_1.combine_score < tag_info_2.combine_score );
	}

	void
	StepWiseRNA_CombineLongLoopFilterer::sort_Combine_Tags_Info( utility::vector1< Combine_Tags_Info > & combine_tags_info_list ) { 	//Lowest combine score on the top of the list
		sort( combine_tags_info_list.begin(), combine_tags_info_list.end(), score_sort_citeria );
	}

  /////////////////////////////////////////////////////////////////////////////////////////////////

	std::string
	StepWiseRNA_CombineLongLoopFilterer::get_parent_tag( utility::vector1< std::string > const & tag_token ) const {

		if ( tag_token.size() < 2 ) utility_exit_with_message( "tag_token.size() < 2" );

		return ( tag_token[1] + "_" + tag_token[2] ); //Misnomer!

	}

	bool
	StepWiseRNA_CombineLongLoopFilterer::Is_virt_sample_ribose_tag( std::string const & tag, utility::vector1< std::string > const & tag_token ) const {

		
		if ( tag_token.size() >= 5 ){ //VIRT_RIBOSE_SAMPLED tag

			//consistency check!
			if ( tag_token[3] != "sample" ) utility_exit_with_message( "tag_token[3] != \"sample\" for tag = " + tag );

			if ( tag_token[4] != "ribose" ) utility_exit_with_message( "tag_token[4] != \"ribose\" for tag = " + tag );

			return true;

		} else if ( tag_token.size() == 2 ){ //standard tag
			return false;
		}

		utility_exit_with_message( "tag_token.size() != 5 and tag_token.size() != 2 for tag = " + tag );
		return false; //PREVENT COMPILER COMPLAINT!	

	}


	bool
	StepWiseRNA_CombineLongLoopFilterer::Is_sibling_ribose_rotamer_pose( std::string const & curr_tag, 
																																	 std::string const & prev_tag, 
																																	 std::map< std::string, std::string > const & tag_to_source_map ) const{


		utility::vector1< std::string > const curr_tag_token = Tokenize( curr_tag, "_" );
		std::string const curr_parent_tag = get_parent_tag( curr_tag_token );

		utility::vector1< std::string > const prev_tag_token = Tokenize( prev_tag, "_" );
		std::string const prev_parent_tag = get_parent_tag( prev_tag_token );

		bool const curr_tag_is_virt_sample_ribose = Is_virt_sample_ribose_tag( curr_parent_tag, curr_tag_token );

		bool const prev_tag_is_virt_sample_ribose = Is_virt_sample_ribose_tag( prev_parent_tag, prev_tag_token );

		//The parent_tag (e.g S_0) should be the same only if pose originate from the same struct before the VIRT_RIBOSE_SAMPLING!
		if ( prev_parent_tag != curr_parent_tag ) return false;

		/////More Consistency check/////
		//If same src struct, then if one pose have the virt_ribose sampled then so must the other.
		if ( curr_tag_is_virt_sample_ribose != prev_tag_is_virt_sample_ribose ){
			utility_exit_with_message( "curr_tag_is_virt_sample_ribose != prev_tag_is_virt_sample_ribose" );
		}

		if ( tag_to_source_map.count( curr_tag ) == 0 ){
			 utility_exit_with_message( "tag " + curr_tag + " doesn't exist in tag_to_source_map" );
		}

		if ( tag_to_source_map.count( prev_tag ) == 0 ){
			 utility_exit_with_message( "tag " + prev_tag + " doesn't exist in tag_to_source_map" );
		}

		std::string const curr_source_file = tag_to_source_map.find( curr_tag )->second;
		std::string const prev_source_file = tag_to_source_map.find( prev_tag )->second;

		//if(source_file_ONE!="MISSING" && prev_source_file_ONE!="MISSING"){
		if ( curr_source_file != prev_source_file ) utility_exit_with_message( "curr_source_file != prev_source_file" );
		//}
		////////////////////////////

		return true;


	}


  /////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_CombineLongLoopFilterer::filter() {

		using namespace ObjexxFCL;

		clock_t const time_start( clock() ); 

		Output_title_text( "Enter StepWiseRNA_CombineLongLoopFilterer:apply", TR );


		Output_boolean( "parin_favorite_ouput = ", parin_favorite_output_, TR ); TR << std::endl;
		Output_boolean( " combine_helical_silent_file_ = ", combine_helical_silent_file_, TR ); TR << std::endl;
		Output_boolean( " filter_for_previous_contact_ = ", filter_for_previous_contact_, TR ); TR << std::endl;
		Output_boolean( " filter_for_previous_clash_ = ", filter_for_previous_clash_, TR ); TR << std::endl;
		Output_boolean( " filter_for_chain_closable_ = ", filter_for_chain_closable_, TR ); TR << std::endl;
		Output_boolean( " filter_for_moving_res_contact_ = ", filter_for_moving_res_contact_, TR ); TR << std::endl;
		Output_boolean( " moving_res_to_base_contact_only_ = ", moving_res_to_base_contact_only_, TR ); TR << std::endl;
		TR << "max_decoys_( nstruct ) = " << max_decoys_ << std::endl;

		setup_silent_file_stream();
		figure_out_NUM_pose_list();

		if ( undercount_ribose_rotamers_ ) setup_tag_to_source_map();


		side_ONE_pose_list_id_ = 1;
		side_TWO_pose_list_id_ = 1;

		filterered_combine_tag_info_list_.clear();

		while ( side_ONE_pose_list_id_ <= side_ONE_NUM_pose_list_ ){
	
			do_some_filtering();

			side_TWO_pose_list_id_++;

			if ( side_TWO_pose_list_id_ > side_TWO_NUM_pose_list_ ){
				side_TWO_pose_list_id_ = 1;
				side_ONE_pose_list_id_++;
			}
		}


		TR << "CombineLongLoopFilterer COUNTS ( BEFORE FINAL SCORE SCREENING )" << std::endl;
		TR << pass_screen_struct_pair_ << " out of " << total_input_struct_pair_ << " passed screen" << std::endl;


		TR << "total_count = " << filterer_count_.total_count;
		TR << " score_cut_count = " << filterer_count_.score_cut_count;
		TR << " chain_closable = " << filterer_count_.chain_closable_screen;
		TR << " filter_for_previous_contact = " << filterer_count_.filter_for_previous_contact;
		TR << " filter_for_previous_clash = " << filterer_count_.filter_for_previous_clash;
		TR << " filter_for_moving_res_contact = " << filterer_count_.filter_for_moving_res_contact;
		TR << std::endl;

			

		TR << "StepWiseRNA_CombineLongLoopFilterer::apply: " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;

		if ( ( filterered_combine_tag_info_list_.size() ) != pass_screen_struct_pair_ ){
			utility_exit_with_message( " ( filterered_combine_tag_info_list_.size() ) != pass_screen_struct_pair_" );
		}

		std::ofstream outfile;
		outfile.open( output_filename_.c_str() ); //Opening the file with this command removes all prior content..

		Size pass_screen_struct_pair_ACT = 0;
		Size pass_screen_struct_pair_undercount_ribose_rotamers = 0;


		//TR << "best_combine_score_= " << best_combine_score_ << " score_diff_cut_= " << score_diff_cut_ << std::endl;

		//OK have to sort the filtererer_combine_tag_info_list.
		sort_Combine_Tags_Info( filterered_combine_tag_info_list_ ); //Oct 19,2010

		clock_t const time_start_FINAL_output( clock() ); 

		for ( Size n = 1; n <= filterered_combine_tag_info_list_.size(); n++ ){
		
			Combine_Tags_Info const combine_tag_info = filterered_combine_tag_info_list_[n];

			//if( (combine_tag_info.combine_score) > (best_combine_score_ + score_diff_cut_) ) continue;

			if ( undercount_ribose_rotamers_ ){

				bool match_existing_pair = false;

				for ( Size ii = ( n - 1 ); ii >= 1; ii-- ){

					Combine_Tags_Info const prev_combine_tag_info = filterered_combine_tag_info_list_[ii];

					if ( Is_sibling_ribose_rotamer_pose( combine_tag_info.side_one_tag, prev_combine_tag_info.side_one_tag, tag_to_source_map_ONE_ ) == false ) continue;

					if ( Is_sibling_ribose_rotamer_pose( combine_tag_info.side_two_tag, prev_combine_tag_info.side_two_tag, tag_to_source_map_TWO_ ) == false ) continue;

					TR << "tag_pair: ";	
					TR << "( " << std::setw( 28 ) << std::left << combine_tag_info.side_one_tag       << ", ";
					TR <<         std::setw( 28 ) << std::left << combine_tag_info.side_two_tag       << " )" ;
					TR << " is a sibling of prev_tag_pair: ";
					TR << "( " << std::setw( 28 ) << std::left << prev_combine_tag_info.side_one_tag  << ", ";
					TR <<         std::setw( 28 ) << std::left << prev_combine_tag_info.side_two_tag  << " )" ;
					TR << std::endl;

					match_existing_pair = true;

					break;					

				} 


				if ( match_existing_pair == false ) pass_screen_struct_pair_undercount_ribose_rotamers++;
				
			}

			pass_screen_struct_pair_ACT++;

			//new formatting, Oct 19, 2010
			outfile << std::setw( 40 ) << std::left << combine_tag_info.side_one_tag; //40 spacing just to be safe!
			outfile << std::setw( 40 ) << std::left << combine_tag_info.side_two_tag; //40 spacing just to be safe!
			outfile << std::setw( 15 ) << std::left << combine_tag_info.combine_score;
			outfile << std::setw( 15 ) << std::left << pass_screen_struct_pair_ACT;
			outfile << "\n"; 

			bool max_decoys_reached = false;

			if ( undercount_ribose_rotamers_ ){
				if ( pass_screen_struct_pair_undercount_ribose_rotamers >= max_decoys_ ) max_decoys_reached = true;
			} else{
				if ( pass_screen_struct_pair_ACT >= max_decoys_ ) max_decoys_reached = true;
			}

			if ( max_decoys_reached ){
				TR << " max_decoys_ ( " << max_decoys_ << " ), early break! " << std::endl;
				break;
			}

		}

		if ( pass_screen_struct_pair_ACT == 0 ){
			outfile << "Empty filterer_outfile. No struct_pair passed screen.\n";
		}

		outfile.flush();
		outfile.close();

		TR << "CombineLongLoopFilterer COUNTS ( AFTER FINAL SCORE SCREENING )" << std::endl;
		TR << "pass_screen_struct_pair_ACT = " << pass_screen_struct_pair_ACT << std::endl;
		TR << "pass_screen_struct_pair_undercount_ribose_rotamers = " << pass_screen_struct_pair_undercount_ribose_rotamers << std::endl;
		TR << "StepWiseRNA_CombineLongLoopFilterer::final_output: " << static_cast< Real > ( clock() - time_start_FINAL_output ) / CLOCKS_PER_SEC << std::endl;
		TR << "StepWiseRNA_CombineLongLoopFilterer::apply: " << static_cast< Real > ( clock() - time_start ) / CLOCKS_PER_SEC << std::endl;


		Output_title_text( "Exit StepWiseRNA_CombineLongLoopFilterer:apply", TR );

	}
}
}
}
