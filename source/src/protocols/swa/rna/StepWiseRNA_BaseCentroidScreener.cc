// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_BaseCentroidScreener
/// @brief Parameters to be passed between different modules of stepwise RNA building.
/// @detailed
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/swa/rna/StepWiseRNA_Classes.hh>
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rna/RNA_CentroidInfo.hh>
#include <basic/Tracer.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <string>

using namespace core;
using basic::T;
using core::Real;
using ObjexxFCL::fmt::F;

static basic::Tracer TR( "protocols.swa.rna.base_centroid_screener" );

namespace protocols {
namespace swa {
namespace rna {

	//////////////////////////////////////////////////////////////////////////////////////////
	// Constructor
	StepWiseRNA_BaseCentroidScreener::StepWiseRNA_BaseCentroidScreener( core::pose::Pose const & pose, StepWiseRNA_JobParametersCOP & job_parameters ):
		job_parameters_( job_parameters ),
		rna_centroid_info_( new core::scoring::rna::RNA_CentroidInfo ),
		base_stack_dist_cutoff_( 6.364 ),
		base_stack_z_offset_max_( 4.5 ),
		base_stack_z_offset_min_( 2.5 ),
		base_stack_axis_cutoff_( 0.707 /*Rhiju value is 0.650*/ ),
		base_stack_planarity_cutoff_( 0.707 /*Rhiju value is 0.5*/ ),
		base_pair_dist_min_( 5.0 ),
		base_pair_dist_max_( 12.0 ),
		base_pair_z_offset_cutoff_( 3.0 ),
		base_pair_axis_cutoff_( 0.5 ),
		base_pair_planarity_cutoff_( 0.866 ),
		base_pair_rho_min_( 5 ),
		base_pair_rho_max_( 10 )
	{
		Initialize_is_virtual_base( pose,  true /*verbose*/ );
		Initialize_base_stub_list( pose, true /*verbose*/ );
		Initialize_terminal_res( pose );
	}
	////////////////////////////////////////////////////////////////////////
	StepWiseRNA_BaseCentroidScreener::~StepWiseRNA_BaseCentroidScreener(){}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_BaseCentroidScreener::Initialize_is_virtual_base( pose::Pose const & pose, bool const ){

		Size const & nres = pose.total_residue();
		is_virtual_base_.dimension( nres, false );

		for ( Size seq_num = 1; seq_num <= nres; seq_num++ ){

			conformation::Residue const & residue_object = pose.residue( seq_num );

			if ( residue_object.has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ){
				TR << "Residue " << seq_num << " is a VIRTUAL_RNA_RESIDUE!" << std::endl;
				is_virtual_base_( seq_num ) = true;
			}

			if ( residue_object.has_variant_type( "BULGE" ) ){
				TR << "Residue " << seq_num << " is a BULGE!" << std::endl;
				is_virtual_base_( seq_num ) = true;
			}
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_BaseCentroidScreener::Initialize_base_stub_list( pose::Pose const & pose, bool const verbose ){

		ObjexxFCL::FArray1D < bool > const & partition_definition = job_parameters_->partition_definition();

		//Used to be this before March 19, 2012. Fang switched to the version below to fix a memory leak problem.
		//bool const root_partition = partition_definition( pose.fold_tree().root() );

		//To Fang: This version is buggy since it assumes full-length pose which is not always the case!
		//Fang modified to this version On March 18, 2012
		//Size const moving_res = job_parameters_ -> moving_res();
		//bool const moving_partition = partition_definition( moving_res );

		Size const working_moving_res = job_parameters_->working_moving_res();
		utility::vector1 < core::Size > const & working_moving_partition_pos = job_parameters_->working_moving_partition_pos();

		if ( working_moving_partition_pos.size() == 0 ) utility_exit_with_message( "working_moving_partition_pos.size() == 0!" );

		bool const moving_partition = partition_definition( working_moving_res );
		bool const moving_partition_check = partition_definition( working_moving_partition_pos[1] );

		if ( moving_partition != moving_partition_check ){
			TR << "working_moving_res = " << working_moving_res << std::endl;
			Output_seq_num_list( "working_moving_partition_pos = ", job_parameters_->working_moving_partition_pos(), TR );
			TR << "moving_partition = " << moving_partition << std::endl;
			TR << "moving_partition_check = " << moving_partition_check << std::endl;
			utility_exit_with_message( "moving_partition != moving_partition_check!" );
		}

		moving_residues_.clear();
		fixed_residues_.clear();
		base_stub_list_.clear();


		Size const & nres = pose.total_residue();
		is_moving_res_.dimension( nres, false );
		is_fixed_res_.dimension( nres, false );

		for ( Size seq_num = 1; seq_num <= nres; seq_num++ ){
			conformation::Residue const & residue_object = pose.residue( seq_num );
			if ( residue_object.aa() == core::chemical::aa_vrt ) continue;
			core::kinematics::Stub base_stub;

			if ( is_virtual_base_( seq_num ) == true ){
				base_stub = core::kinematics::Stub();	//"default" tub, this will never be called
			} else{
				Vector const centroid = rna_centroid_info_->get_base_centroid( residue_object );
				base_stub = rna_centroid_info_->get_base_coordinate_system( residue_object, centroid );
			}

			base_stub_list_.push_back( base_stub );

			if ( is_virtual_base_( seq_num ) == true ) continue;

			if ( partition_definition( seq_num ) != moving_partition ) {
				// This is a "fixed" residue -- on the same side of the moving suite as the root.
				fixed_residues_.push_back( seq_num );
				//if ( verbose ) TR << " FIXED POSITION  --> " << seq_num << std::endl;
				is_fixed_res_( seq_num ) = true;
			} else {
				moving_residues_.push_back( seq_num );
				//				if ( verbose ) TR << " MOVING POSITION --> " << seq_num << std::endl;
				is_moving_res_( seq_num ) = true;
			}
		}


//		moving_residues_ does not necessarily equal job_parameters_->working_moving_partition_pos since job_parameters_->working_moving_partition_pos include virtual residues. May 25, 2010
//		if(Is_equivalent_vector(moving_residues_,job_parameters_->working_moving_partition_pos())==false){
//			Output_seq_num_list("moving_residues_= ", moving_residues_, 50, TR );
//			Output_seq_num_list("job_parameters_->working_moving_partition_pos()= ", job_parameters_->working_moving_partition_pos(), 50, TR );
//			utility_exit_with_message( "moving_residues_,job_parameters_->working_moving_partition_pos()) ==false ");
//		}

	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_BaseCentroidScreener::Initialize_terminal_res( pose::Pose const & pose ){

		using namespace ObjexxFCL;


		terminal_res_ = job_parameters_->working_terminal_res();

		Size const & nres = pose.total_residue();
		is_terminal_res_.dimension( nres, false );
		stacked_on_terminal_res_in_original_pose_.dimension( nres, nres, false );

		for ( Size n = 1; n <= terminal_res_.size(); n++ ) {

			//			TR << "NRES " << pose.total_residue() << " " << is_terminal_res_.size() << "      " << terminal_res_[ n ] << std::endl;
			Size const terminal_res = terminal_res_[ n ];

			if ( is_virtual_base_( terminal_res ) == true ){
				utility_exit_with_message( "working_res: " + string_of( terminal_res ) + " is a terminal res but has a virtual! " );
			}

			is_terminal_res_( terminal_res ) = true;

			for ( Size m = 1; m <= nres; m++ ) {

				//				TR << " about to check stack: --- " << std::endl;
				//				TR << " TERMINAL_RES " << terminal_res << " " << is_moving_res_( terminal_res ) <<  " " << is_fixed_res_( terminal_res ) << std::endl;
				//				TR << " M            " << m << " " << is_moving_res_( m ) <<  " " << is_fixed_res_( m ) << std::endl;

				if ( ( is_moving_res_( terminal_res )  && is_moving_res_( m ) ) ||
						 ( is_fixed_res_(  terminal_res )  && is_fixed_res_(  m ) ) ){

					stacked_on_terminal_res_in_original_pose_( terminal_res, m )  = check_stack_base( terminal_res, m );
					//if ( stacked_on_terminal_res_in_original_pose_( terminal_res, m ) ) TR << "ALREADY STACKED: " << terminal_res << " " << m << std::endl;

				}
			}
		}

	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseCentroidScreener::check_stack_base( core::kinematics::Stub const & rebuild_residue_base_stub, core::kinematics::Stub const & base_stub, bool const verbose  ) const{

		numeric::xyzVector< Real > const other_z_vector = base_stub.M.col_z();
		numeric::xyzVector< Real > rebuild_z_vector = rebuild_residue_base_stub.M.col_z();

		numeric::xyzVector< Real > centroid_diff;
		subtract( rebuild_residue_base_stub.v, base_stub.v, centroid_diff );
		Real centroid_distance = centroid_diff.length();
		if ( verbose ) TR << "Centroid Distance: " << centroid_distance << std::endl;
		if ( centroid_distance > base_stack_dist_cutoff_ ) return false;

		Real base_z_offset_one = std::abs( dot( centroid_diff, other_z_vector ) );
		Real base_z_offset_two = std::abs( dot( centroid_diff, rebuild_z_vector ) );

		if ( verbose ) TR << "Base Z offset 1: " << base_z_offset_one << std::endl;
		if ( verbose ) TR << "Base Z offset 2: " << base_z_offset_two << std::endl;

		if ( ( base_z_offset_one > base_stack_z_offset_max_ || base_z_offset_one <  base_stack_z_offset_min_ ) &&
				 ( base_z_offset_two > base_stack_z_offset_max_ || base_z_offset_two <  base_stack_z_offset_min_ ) ) return false;

		Real base_axis_one = base_z_offset_one/centroid_distance;
		Real base_axis_two = base_z_offset_two/centroid_distance;

		if ( verbose ) TR << "Base Axis 1: " << base_axis_one << std::endl;
		if ( verbose ) TR << "Base Axis 2: " << base_axis_two << std::endl;

		if ( base_axis_one < base_stack_axis_cutoff_ && base_axis_two < base_stack_axis_cutoff_ ) return false;

		Real base_planarity = std::abs( dot( other_z_vector, rebuild_z_vector ) );

		if ( verbose ) TR << "Base planarity: " << base_planarity << std::endl;

		if ( base_planarity < base_stack_planarity_cutoff_ ) return false;

		return true;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseCentroidScreener::check_base_pairing( core::kinematics::Stub const & rebuild_residue_base_stub, core::kinematics::Stub const & base_stub   ) const{

		numeric::xyzVector< Real > const other_z_vector = base_stub.M.col_z();
		numeric::xyzVector< Real > rebuild_z_vector = rebuild_residue_base_stub.M.col_z();

		numeric::xyzVector< Real > centroid_diff;
		subtract( rebuild_residue_base_stub.v, base_stub.v, centroid_diff );

		Real centroid_distance = centroid_diff.length();
		if ( centroid_distance < base_pair_dist_min_ || centroid_distance > base_pair_dist_max_ ) return false;

		Real base_z_offset_one = std::abs( dot( centroid_diff, other_z_vector ) );
		Real base_z_offset_two = std::abs( dot( centroid_diff, rebuild_z_vector ) );

		if ( base_z_offset_one > base_pair_z_offset_cutoff_  && base_z_offset_two > base_pair_z_offset_cutoff_ ) return false;

		Real base_axis_one = base_z_offset_one/centroid_distance;
		Real base_axis_two = base_z_offset_two/centroid_distance;
		if ( base_axis_one > base_pair_axis_cutoff_ && base_axis_two > base_pair_axis_cutoff_ ) return false; //This is a stronger condition compare to baze_z_offset check

		Real base_planarity = std::abs( dot( rebuild_z_vector, other_z_vector ) );
		if ( base_planarity < base_pair_planarity_cutoff_ ) return false;

		numeric::xyzVector< Real > centroid_diff_parallel_one = dot( centroid_diff, other_z_vector )*other_z_vector;
		numeric::xyzVector< Real > centroid_diff_perpendicular_one = centroid_diff - centroid_diff_parallel_one;
		Real rho_one = centroid_diff_perpendicular_one.length();
		numeric::xyzVector< Real > centroid_diff_parallel_two = dot( centroid_diff, rebuild_z_vector )*rebuild_z_vector;
		numeric::xyzVector< Real > centroid_diff_perpendicular_two = centroid_diff - centroid_diff_parallel_two;
		Real rho_two = centroid_diff_perpendicular_two.length();

		if ( ( rho_one < base_pair_rho_min_ || rho_one > base_pair_rho_max_ ) &&
				 ( rho_two < base_pair_rho_min_ || rho_two > base_pair_rho_max_ ) ) return false;

		//If reach this point means success!
		return true;
	}


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseCentroidScreener::Update_base_stub_list_and_Check_centroid_interaction( core::pose::Pose const & pose, SillyCountStruct & count_data ){
		Update_base_stub_list( pose );
		return Check_centroid_interaction( count_data );
	}


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_BaseCentroidScreener::Update_base_stub_list( core::pose::Pose const & pose ){

		for ( Size m = 1; m <= moving_residues_.size(); m++ ) {

			Size const moving_res( moving_residues_[ m ] );

			core::conformation::Residue const & residue_object( pose.residue( moving_res ) );

			Vector const centroid = rna_centroid_info_->get_base_centroid( residue_object );
			core::kinematics::Stub base_stub = rna_centroid_info_->get_base_coordinate_system( residue_object, centroid );
			base_stub_list_[ moving_res ] = base_stub;

		}

	}


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseCentroidScreener::Check_centroid_interaction( SillyCountStruct & count_data ) const{

		bool stack_base( false ), base_pairing( false );

		for ( Size m = 1; m <= moving_residues_.size(); m++ ) {

			core::kinematics::Stub const & rebuild_residue_base_stub = base_stub_list_[ moving_residues_[ m ] ];
			stack_base = false;

			for ( Size i = 1; i <= fixed_residues_.size(); i++ ){
				stack_base = check_stack_base( rebuild_residue_base_stub, base_stub_list_[ fixed_residues_[ i ] ] );
				if ( stack_base ) break;
			}

			base_pairing = false;

			for ( Size i = 1; i <= fixed_residues_.size(); i++ ) {
				base_pairing = check_base_pairing( rebuild_residue_base_stub, base_stub_list_[ fixed_residues_[ i ] ] );
				if ( base_pairing ) break;
			}

			if ( base_pairing || stack_base ) break; // found an interaction!

		}

		if ( base_pairing ) count_data.base_pairing_count++;
		if ( stack_base ) count_data.base_stack_count++;
		if ( base_pairing || stack_base ) count_data.pass_base_centroid_screen++;

		//		if ( stack_base ) count_data_.base_stack_count++;
		//		if ( base_pairing ) count_data_.base_pairing_count++;

		if ( !base_pairing && !stack_base ) return false;


		//		TR << " BASE_PAIRING " << base_pairing << "  BASE_STACKING " << stack_base << std::endl;
		return true;

	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseCentroidScreener::Update_base_stub_list_and_Check_that_terminal_res_are_unstacked( core::pose::Pose const & pose, bool const reinitialize /* = false */ ){
		if ( reinitialize ) {
			Initialize_base_stub_list( pose );
		} else {
			Update_base_stub_list( pose );
		}
		bool const passed = Check_that_terminal_res_are_unstacked( false /*verbose*/ );
		return passed;
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseCentroidScreener::check_stack_base( Size const & pos1, Size const & pos2, bool const verbose /* = false */  ) {

		if ( is_virtual_base_( pos1 ) == true || is_virtual_base_( pos2 ) == true ){
			utility_exit_with_message( "is_virtual_base_( pos1 ) == true || is_virtual_base_( pos2 ) == true !" );
		}

		if ( pos1 == pos2 ) return true;
		return check_stack_base(  base_stub_list_[ pos1 ], base_stub_list_[ pos2 ], verbose );

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseCentroidScreener::Check_that_terminal_res_are_unstacked( bool const verbose ){

		//		bool stack_base( false );//, base_pairing( false );

		// Look through all terminal_res
		for ( Size i = 1; i <= terminal_res_.size(); i++ ) {
			Size const & terminal_res = terminal_res_[ i ];

			for ( Size m = 1; m <= moving_residues_.size(); m++ ) {
				Size const & moving_res = moving_residues_[ m ];
				if ( verbose ) TR << "about to check stack: " << terminal_res << " " << moving_res << " " << stacked_on_terminal_res_in_original_pose_( terminal_res, moving_res ) << std::endl;
				if ( !stacked_on_terminal_res_in_original_pose_( terminal_res, moving_res ) &&
						 check_stack_base( terminal_res, moving_res, verbose  ) ) return false;
			}

			for ( Size m = 1; m <= fixed_residues_.size(); m++ ) {
				Size const & fixed_res = fixed_residues_[ m ];
				if ( verbose ) TR << "about to check stack: " << terminal_res << " " << fixed_res << " " << stacked_on_terminal_res_in_original_pose_( terminal_res, fixed_res ) << std::endl;
				if ( !stacked_on_terminal_res_in_original_pose_( terminal_res, fixed_res ) &&
						 check_stack_base( terminal_res, fixed_res, verbose  ) ) return false;
			}

		}

		return true;

	}


}
}
}
