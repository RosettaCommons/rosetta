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
#include <protocols/swa/rna/screener/StepWiseRNA_BaseCentroidScreener.hh>
#include <protocols/swa/rna/StepWiseRNA_Util.hh>
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
using ObjexxFCL::format::F;

static basic::Tracer TR( "protocols.swa.rna.screener.StepWiseRNA_BaseCentroidScreener" );

namespace protocols {
namespace swa {
namespace rna {
namespace screener {

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
		base_pair_rho_max_( 10 ),
		allow_base_pair_only_screen_( false )
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
				TR.Debug << "Residue " << seq_num << " is a VIRTUAL_RNA_RESIDUE!" << std::endl;
				is_virtual_base_( seq_num ) = true;
			}

			if ( residue_object.has_variant_type( "BULGE" ) ){
				TR.Debug << "Residue " << seq_num << " is a BULGE!" << std::endl;
				is_virtual_base_( seq_num ) = true;
			}
		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_BaseCentroidScreener::Initialize_base_stub_list( pose::Pose const & pose, bool const verbose ){

		ObjexxFCL::FArray1D < bool > const & partition_definition = job_parameters_->partition_definition();

		Size const working_moving_res = job_parameters_->working_moving_res();
		utility::vector1 < core::Size > const & working_moving_partition_pos = job_parameters_->working_moving_partition_pos();
		runtime_assert( working_moving_partition_pos.size() > 0 );

		bool const moving_partition_based_on_res = partition_definition( working_moving_res );
		bool const moving_partition = partition_definition( working_moving_partition_pos[1] );
		// this had to be disabled for stepwise monte carlo stuff -- root of fold tree sometimes is on same partition
		// as moving residue.
		//		runtime_assert( moving_partition_based_on_res == moving_partition );

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

			if ( is_virtual_base_( seq_num ) ){
				base_stub = core::kinematics::Stub();	//"default" stub, this will never be called
			} else{
				Vector const centroid = rna_centroid_info_->get_base_centroid( residue_object );
				base_stub = rna_centroid_info_->get_base_coordinate_system( residue_object, centroid );
			}
			base_stub_list_.push_back( base_stub );

			if ( is_virtual_base_( seq_num ) == true ) continue;

			if ( partition_definition( seq_num ) != moving_partition ) {
				// This is a "fixed" residue -- on the same side of the moving suite as the root.
				fixed_residues_.push_back( seq_num );
				is_fixed_res_( seq_num ) = true;
			} else {
				moving_residues_.push_back( seq_num );
				is_moving_res_( seq_num ) = true;
			}
		}
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

			Size const terminal_res = terminal_res_[ n ];
			if ( is_virtual_base_( terminal_res ) == true ){
				utility_exit_with_message( "working_res: " + string_of( terminal_res ) + " is a terminal res but has a virtual! " );
			}

			is_terminal_res_( terminal_res ) = true;

			for ( Size m = 1; m <= nres; m++ ) {

				if ( ( is_moving_res_( terminal_res )  && is_moving_res_( m ) ) ||
						 ( is_fixed_res_(  terminal_res )  && is_fixed_res_(  m ) ) ){

					stacked_on_terminal_res_in_original_pose_( terminal_res, m )  = check_base_stack( terminal_res, m );
				}
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseCentroidScreener::check_base_stack( core::kinematics::Stub const & moving_residue_base_stub,
																											core::kinematics::Stub const & other_base_stub,
																											core::Real const base_axis_CUTOFF,
																											core::Real const base_planarity_CUTOFF,
																											bool const verbose  /* = false */ ) const{

		numeric::xyzVector< Real > const other_z_vector = other_base_stub.M.col_z();
		numeric::xyzVector< Real > rebuild_z_vector = moving_residue_base_stub.M.col_z();

		numeric::xyzVector< Real > centroid_diff;
		subtract( moving_residue_base_stub.v, other_base_stub.v, centroid_diff );
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

		if ( base_axis_one < base_axis_CUTOFF && base_axis_two < base_axis_CUTOFF ) return false;

		Real base_planarity = std::abs( dot( other_z_vector, rebuild_z_vector ) );

		if ( verbose ) TR << "Base planarity: " << base_planarity << std::endl;

		if ( base_planarity < base_planarity_CUTOFF ) return false;

		return true;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseCentroidScreener::check_base_stack( core::kinematics::Stub const & moving_residue_base_stub,
																											core::kinematics::Stub const & other_base_stub,
																											bool const verbose  /* = false */ ) const{
		return check_base_stack( moving_residue_base_stub, other_base_stub, base_stack_axis_cutoff_, base_stack_planarity_cutoff_ );
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseCentroidScreener::check_base_stack( core::kinematics::Stub const & moving_res_base_stub,
																											utility::vector1 < core::kinematics::Stub > const & other_residues_base_list,
																											core::Real const base_axis_CUTOFF,
																											core::Real const base_planarity_CUTOFF ) const {

		for ( Size i = 1; i <= other_residues_base_list.size(); i++ ){
			core::kinematics::Stub const & other_base_stub = other_residues_base_list[i];
			if ( check_base_stack( moving_res_base_stub, other_base_stub, base_axis_CUTOFF, base_planarity_CUTOFF ) ) return true;
		}
		return false;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseCentroidScreener::check_base_pair( core::kinematics::Stub const & moving_residue_base_stub,
																										 core::kinematics::Stub const & other_base_stub,
																										 core::Real const base_axis_CUTOFF,
																										 core::Real const base_planarity_CUTOFF   ) const{

		numeric::xyzVector< Real > const other_z_vector = other_base_stub.M.col_z();
		numeric::xyzVector< Real > rebuild_z_vector = moving_residue_base_stub.M.col_z();

		numeric::xyzVector< Real > centroid_diff;
		subtract( moving_residue_base_stub.v, other_base_stub.v, centroid_diff );

		Real centroid_distance = centroid_diff.length();
		if ( centroid_distance < base_pair_dist_min_ || centroid_distance > base_pair_dist_max_ ) return false;

		Real base_z_offset_one = std::abs( dot( centroid_diff, other_z_vector ) );
		Real base_z_offset_two = std::abs( dot( centroid_diff, rebuild_z_vector ) );

		if ( base_z_offset_one > base_pair_z_offset_cutoff_  && base_z_offset_two > base_pair_z_offset_cutoff_ ) return false;

		Real base_axis_one = base_z_offset_one/centroid_distance;
		Real base_axis_two = base_z_offset_two/centroid_distance;
		if ( base_axis_one > base_axis_CUTOFF && base_axis_two > base_axis_CUTOFF ) return false; //This is a stronger condition compare to baze_z_offset check

		Real base_planarity = std::abs( dot( rebuild_z_vector, other_z_vector ) );
		if ( base_planarity < base_planarity_CUTOFF ) return false;

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

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseCentroidScreener::check_base_pair( core::kinematics::Stub const & moving_res_base_stub,
																											utility::vector1 < core::kinematics::Stub > const & other_residues_base_list,
																											core::Real const base_axis_CUTOFF,
																											core::Real const base_planarity_CUTOFF ) const {

		for ( Size i = 1; i <= other_residues_base_list.size(); i++ ){
			core::kinematics::Stub const & other_base_stub = other_residues_base_list[i];
			if ( check_base_pair( moving_res_base_stub, other_base_stub, base_axis_CUTOFF, base_planarity_CUTOFF ) ) return true;
		}
		return false;
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseCentroidScreener::is_strong_base_stack( core::kinematics::Stub const & moving_res_base,
																													utility::vector1 < core::kinematics::Stub > const & other_residues_base_list ) const{

		Real const base_axis_CUTOFF = 0.9000;
		Real const base_planarity_CUTOFF = 0.9000;

		return check_base_stack( moving_res_base, other_residues_base_list, base_axis_CUTOFF, base_planarity_CUTOFF );

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseCentroidScreener::is_medium_base_stack_and_medium_base_pair( core::kinematics::Stub const & moving_res_base,
																																							 utility::vector1 < core::kinematics::Stub > const & other_residues_base_list ) const{

		bool base_stack = check_base_stack( moving_res_base, other_residues_base_list, 0.7070 /*base_axis_CUTOFF*/, 0.7070 /*base_planarity_CUTOFF*/ );

		bool base_pair = check_base_pair( moving_res_base, other_residues_base_list, 0.5000 /*base_axis_CUTOFF*/, 0.7070 /*base_planarity_CUTOFF*/ );
		//value in Base_screener_class is 0.866 Sept 16 2010, Parin S.

		return ( base_stack && base_pair );
	}


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseCentroidScreener::update_base_stub_list_and_check_centroid_interaction( core::pose::Pose const & pose,
																																													StepWiseRNA_CountStruct & count_data ){
		update_base_stub_list( pose );
		return check_centroid_interaction( count_data );
	}


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	void
	StepWiseRNA_BaseCentroidScreener::update_base_stub_list( core::pose::Pose const & pose ){

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
	StepWiseRNA_BaseCentroidScreener::check_centroid_interaction_floating_base( core::kinematics::Stub const &  moving_res_base_stub,
																																							StepWiseRNA_CountStruct & count_data ) const{

		utility::vector1< core::kinematics::Stub > other_residues_base_list;
		for ( Size i = 1; i <= fixed_residues_.size(); i++ ) other_residues_base_list.push_back( base_stub_list_[ fixed_residues_[ i ] ] );

		bool const strong_stack_base = is_strong_base_stack( moving_res_base_stub, other_residues_base_list );
		if ( strong_stack_base ) count_data.base_stack_count++;

		bool const medium_base_stack_and_medium_base_pair = is_medium_base_stack_and_medium_base_pair( moving_res_base_stub, other_residues_base_list );
		if ( medium_base_stack_and_medium_base_pair ) count_data.base_pairing_count++;

		bool strict_base_pair = false;
		if ( allow_base_pair_only_screen_ ){
			strict_base_pair = check_base_pair( moving_res_base_stub, other_residues_base_list, 0.2588 /*base_axis_CUTOFF*/, 0.8660 /*base_planarity_CUTOFF*/ );
			if ( strict_base_pair ) count_data.strict_base_pairing_count++;
		}

		if ( strong_stack_base || medium_base_stack_and_medium_base_pair || ( allow_base_pair_only_screen_ && strict_base_pair ) ){
			count_data.pass_base_centroid_screen++;
			return true;
		}

		return false;
	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseCentroidScreener::check_centroid_interaction( core::kinematics::Stub const &  moving_res_base_stub,
																																StepWiseRNA_CountStruct & count_data ) {

		if ( floating_base_ )	return check_centroid_interaction_floating_base( moving_res_base_stub, count_data );

		runtime_assert( moving_residues_.size() == 1 );
		base_stub_list_[ moving_residues_[ 1 ] ] = moving_res_base_stub;
		return check_centroid_interaction( count_data );

	}

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseCentroidScreener::check_centroid_interaction( StepWiseRNA_CountStruct & count_data ) {

		bool stack_base( false ), base_pair( false );

		for ( Size m = 1; m <= moving_residues_.size(); m++ ) {

			core::kinematics::Stub const & moving_residue_base_stub = base_stub_list_[ moving_residues_[ m ] ];
			stack_base = false;

			for ( Size i = 1; i <= fixed_residues_.size(); i++ ){
				stack_base = check_base_stack( moving_residue_base_stub, base_stub_list_[ fixed_residues_[ i ] ] );
				if ( stack_base ) break;
			}

			base_pair = false;

			for ( Size i = 1; i <= fixed_residues_.size(); i++ ) {
				base_pair = check_base_pair( moving_residue_base_stub, base_stub_list_[ fixed_residues_[ i ] ], base_pair_axis_cutoff_, base_pair_planarity_cutoff_ );
				if ( base_pair ) break;
			}

			if ( base_pair || stack_base ) break; // found an interaction!

		}

		if ( base_pair ) count_data.base_pairing_count++;
		if ( stack_base ) count_data.base_stack_count++;
		if ( base_pair || stack_base ) count_data.pass_base_centroid_screen++;

		//		if ( stack_base ) count_data_.base_stack_count++;
		//		if ( base_pair ) count_data_.base_pairing_count++;

		if ( !base_pair && !stack_base ) return false;


		//		TR << " BASE_PAIR " << base_pair << "  BASE_STACKING " << stack_base << std::endl;
		return true;

	}


	////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseCentroidScreener::update_base_stub_list_and_check_that_terminal_res_are_unstacked( core::pose::Pose const & pose, bool const reinitialize /* = false */ ){
		if ( reinitialize ) {
			Initialize_base_stub_list( pose );
		} else {
			update_base_stub_list( pose );
		}
		bool const passed = check_that_terminal_res_are_unstacked( false /*verbose*/ );
		return passed;
	}


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseCentroidScreener::check_base_stack( Size const & pos1, Size const & pos2, bool const verbose /* = false */  ) {

		if ( is_virtual_base_( pos1 ) == true || is_virtual_base_( pos2 ) == true ){
			utility_exit_with_message( "is_virtual_base_( pos1 ) == true || is_virtual_base_( pos2 ) == true !" );
		}

		if ( pos1 == pos2 ) return true;
		return check_base_stack(  base_stub_list_[ pos1 ], base_stub_list_[ pos2 ], verbose );

	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	bool
	StepWiseRNA_BaseCentroidScreener::check_that_terminal_res_are_unstacked( bool const verbose ){

		// Look through all terminal_res
		for ( Size i = 1; i <= terminal_res_.size(); i++ ) {
			Size const & terminal_res = terminal_res_[ i ];

			for ( Size m = 1; m <= moving_residues_.size(); m++ ) {
				Size const & moving_res = moving_residues_[ m ];
				if ( verbose ) TR << "about to check stack: " << terminal_res << " " << moving_res << " " << stacked_on_terminal_res_in_original_pose_( terminal_res, moving_res ) << std::endl;
				if ( !stacked_on_terminal_res_in_original_pose_( terminal_res, moving_res ) &&
						 check_base_stack( terminal_res, moving_res, verbose  ) ) return false;
			}

			for ( Size m = 1; m <= fixed_residues_.size(); m++ ) {
				Size const & fixed_res = fixed_residues_[ m ];
				if ( verbose ) TR << "about to check stack: " << terminal_res << " " << fixed_res << " " << stacked_on_terminal_res_in_original_pose_( terminal_res, fixed_res ) << std::endl;
				if ( !stacked_on_terminal_res_in_original_pose_( terminal_res, fixed_res ) &&
						 check_base_stack( terminal_res, fixed_res, verbose  ) ) return false;
			}

		}

		return true;

	}

}
}
}
}
