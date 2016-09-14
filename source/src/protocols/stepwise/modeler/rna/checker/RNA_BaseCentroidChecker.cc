// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file RNA_BaseCentroidChecker
/// @brief Parameters to be passed between different modules of stepwise RNA building.
/// @details
/// @author Rhiju Das


//////////////////////////////////
#include <protocols/stepwise/modeler/rna/checker/RNA_BaseCentroidChecker.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/rna/StepWiseRNA_Classes.hh>
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/Stub.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/rna/RNA_CentroidInfo.hh>
#include <basic/Tracer.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.hh>
#include <ObjexxFCL/FArray1D.hh>
#include <ObjexxFCL/FArray2D.hh>
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <string>

using namespace core;
using basic::T;
using core::Real;
using ObjexxFCL::format::F;
using core::chemical::rna::BaseStackWhichSide;

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.modeler.rna.checker.RNA_BaseCentroidChecker" );

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace checker {


//////////////////////////////////////////////////////////////////////////////////////////
// Constructor
RNA_BaseCentroidChecker::RNA_BaseCentroidChecker( core::pose::Pose const & pose,
	working_parameters::StepWiseWorkingParametersCOP & working_parameters,
	bool const tether_jump /* = false */ ):
	working_parameters_( working_parameters ),
	rna_centroid_info_( core::scoring::rna::RNA_CentroidInfoOP( new core::scoring::rna::RNA_CentroidInfo ) ),
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
	allow_base_pair_only_screen_( false ),
	floating_base_( false ),
	found_centroid_interaction_( false ),
	tether_jump_( tether_jump )
{
	Initialize_is_virtual_base( pose,  true /*verbose*/ );
	Initialize_base_stub_list( pose, true /*verbose*/ );
	Initialize_terminal_res( pose );
}
////////////////////////////////////////////////////////////////////////
RNA_BaseCentroidChecker::~RNA_BaseCentroidChecker(){}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_BaseCentroidChecker::Initialize_is_virtual_base( pose::Pose const & pose, bool const ){

	Size const & nres = pose.size();
	is_virtual_base_.dimension( nres, false );

	for ( Size seq_num = 1; seq_num <= nres; seq_num++ ) {

		conformation::Residue const & residue_object = pose.residue( seq_num );

		if ( residue_object.has_variant_type( core::chemical::VIRTUAL_RNA_RESIDUE ) ) {
			TR.Debug << "Residue " << seq_num << " is a VIRTUAL_RNA_RESIDUE!" << std::endl;
			is_virtual_base_( seq_num ) = true;
		}

		if ( residue_object.has_variant_type( core::chemical::BULGE ) ) {
			TR.Debug << "Residue " << seq_num << " is a BULGE!" << std::endl;
			is_virtual_base_( seq_num ) = true;
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_BaseCentroidChecker::Initialize_base_stub_list( pose::Pose const & pose, bool const  ){

	ObjexxFCL::FArray1D < bool > const & partition_definition = working_parameters_->partition_definition();
	utility::vector1 < core::Size > const & working_moving_partition_res = working_parameters_->working_moving_partition_res();

	//  runtime_assert( working_moving_partition_res.size() > 0 );
	if ( working_moving_partition_res.size() == 0 ) return;
	bool const moving_partition = partition_definition( working_moving_partition_res[1] );
	// this had to be disabled for stepwise monte carlo stuff -- root of fold tree sometimes is on same partition
	// as moving residue.
	// bool const moving_partition_based_on_res = partition_definition( working_moving_res );
	//  runtime_assert( moving_partition_based_on_res == moving_partition );

	moving_residues_.clear();
	fixed_residues_.clear();
	base_stub_list_.clear();

	Size const & nres = pose.size();
	is_moving_res_.dimension( nres, false );
	is_fixed_res_.dimension( nres, false );

	for ( Size seq_num = 1; seq_num <= nres; seq_num++ ) {
		conformation::Residue const & residue_object = pose.residue( seq_num );
		if ( residue_object.aa() == core::chemical::aa_vrt ) continue;
		core::kinematics::Stub base_stub;

		base_stub = core::kinematics::Stub(); //"default" stub, this will never be called
		if ( !is_virtual_base_( seq_num ) &&
				residue_object.is_RNA() ) {
			Vector const centroid = rna_centroid_info_->get_base_centroid( residue_object );
			base_stub = rna_centroid_info_->get_base_coordinate_system( residue_object, centroid );
		}
		base_stub_list_.push_back( base_stub );

		if ( is_virtual_base_( seq_num ) ) continue;

		if ( partition_definition( seq_num ) != moving_partition ) {
			if ( tether_jump_ && seq_num != working_parameters_->working_reference_res() ) continue;
			// This is a "fixed" residue -- on the same side of the moving suite as the root.
			fixed_residues_.push_back( seq_num );
			is_fixed_res_( seq_num ) = true;
			//    TR << "Fixed partition " << seq_num << std::endl;
		} else {
			if ( tether_jump_ && seq_num != working_parameters_->working_moving_res() ) continue;
			moving_residues_.push_back( seq_num );
			is_moving_res_( seq_num ) = true;
			//    TR << "Moving partition " << seq_num << std::endl;
		}
	}
}


//////////////////////////////////////////////////////////////////////////////////////////
//
// Note that terminal_res setup is kind of complicated below --
//  Some terminal res are allowed to be stacked on some other residues if they
//  are in the same partition.
// This should soon be deprecated by block_stack_above, block_stack_below...
//  Sorry, rhiju.
//
//////////////////////////////////////////////////////////////////////////////////////////
void
RNA_BaseCentroidChecker::Initialize_terminal_res( pose::Pose const & pose ){

	using namespace ObjexxFCL;

	terminal_res_ = working_parameters_->working_terminal_res();

	Size const & nres = pose.size();
	is_terminal_res_.dimension( nres, false );
	stacked_on_terminal_res_in_original_pose_.dimension( nres, nres, false );

	for ( Size const terminal_res : terminal_res_ ) {

		if ( is_virtual_base_( terminal_res ) ) {
			TR << TR.Red << pose.fold_tree() << TR.Reset << std::endl;
			TR << TR.Red << pose.annotated_sequence() << TR.Reset << std::endl;
			utility_exit_with_message( "working_res: " + string_of( terminal_res ) + " is a terminal res but has a virtual! " );
		}

		is_terminal_res_( terminal_res ) = true;

		for ( Size m = 1; m <= nres; m++ ) {
			if ( ( is_moving_res_( terminal_res )  && is_moving_res_( m ) ) ||
					( is_fixed_res_(  terminal_res )  && is_fixed_res_( m ) ) ) {
				stacked_on_terminal_res_in_original_pose_( terminal_res, m )  = check_base_stack( terminal_res, m );
			}
		}
	}

	block_stack_above_res_ = working_parameters_->working_block_stack_above_res();
	block_stack_below_res_ = working_parameters_->working_block_stack_below_res();
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_BaseCentroidChecker::check_base_stack(
	core::kinematics::Stub const & moving_residue_base_stub,
	core::kinematics::Stub const & other_base_stub,
	core::Real const base_axis_CUTOFF,
	core::Real const base_planarity_CUTOFF,
	BaseStackWhichSide & base_stack_side,
	bool const verbose  /* = false */ ) const
{
	numeric::xyzVector< Real > const other_z_vector = other_base_stub.M.col_z();
	numeric::xyzVector< Real > rebuild_z_vector = moving_residue_base_stub.M.col_z();

	numeric::xyzVector< Real > centroid_diff;
	subtract( moving_residue_base_stub.v, other_base_stub.v, centroid_diff );
	Real centroid_distance = centroid_diff.length();
	if ( verbose ) TR << "Centroid Distance: " << centroid_distance << std::endl;
	if ( centroid_distance > base_stack_dist_cutoff_ ) return false;

	Real base_z_offset_one = std::abs( dot( centroid_diff, other_z_vector ) );
	Real base_z_offset_two = std::abs( dot( centroid_diff, rebuild_z_vector ) );

	using namespace core::chemical::rna;
	base_stack_side = ( dot( centroid_diff, rebuild_z_vector ) < 0 ) ? ABOVE : BELOW;

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

//////////////////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_BaseCentroidChecker::check_base_stack( core::kinematics::Stub const & moving_residue_base_stub,
	core::kinematics::Stub const & other_base_stub,
	core::Real const base_axis_CUTOFF,
	core::Real const base_planarity_CUTOFF,
	bool const verbose  /* = false */ ) const
{
	BaseStackWhichSide base_stack_side( core::chemical::rna::ANY_BASE_STACK_SIDE );
	return check_base_stack( moving_residue_base_stub, other_base_stub, base_axis_CUTOFF, base_planarity_CUTOFF, base_stack_side, verbose );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_BaseCentroidChecker::check_base_stack(
	core::kinematics::Stub const & moving_residue_base_stub,
	core::kinematics::Stub const & other_base_stub,
	BaseStackWhichSide & base_stack_side )
{
	return check_base_stack( moving_residue_base_stub, other_base_stub, base_stack_axis_cutoff_, base_stack_planarity_cutoff_, base_stack_side, false );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_BaseCentroidChecker::check_base_stack( core::kinematics::Stub const & moving_residue_base_stub,
	core::kinematics::Stub const & other_base_stub,
	bool const /* verbose = false */ ) const{
	return check_base_stack( moving_residue_base_stub, other_base_stub, base_stack_axis_cutoff_, base_stack_planarity_cutoff_ );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_BaseCentroidChecker::check_base_stack( core::kinematics::Stub const & moving_res_base_stub,
	core::Real const base_axis_CUTOFF,
	core::Real const base_planarity_CUTOFF ) const {

	for ( Size const fixed_res : fixed_residues_ ) {
		core::kinematics::Stub const & other_base_stub = base_stub_list_[ fixed_res ];
		if ( check_base_stack( moving_res_base_stub, other_base_stub, base_axis_CUTOFF, base_planarity_CUTOFF ) ) return true;
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_BaseCentroidChecker::check_base_pair( core::kinematics::Stub const & moving_residue_base_stub,
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
RNA_BaseCentroidChecker::check_base_pair( core::kinematics::Stub const & moving_res_base_stub,
	core::Real const base_axis_CUTOFF,
	core::Real const base_planarity_CUTOFF ) const {

	for ( Size const fixed_residue : fixed_residues_ ) {
		core::kinematics::Stub const & other_base_stub = base_stub_list_[ fixed_residue ];
		if ( check_base_pair( moving_res_base_stub, other_base_stub, base_axis_CUTOFF, base_planarity_CUTOFF ) ) return true;
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_BaseCentroidChecker::is_strong_base_stack( core::kinematics::Stub const & moving_res_base ) const{

	Real const base_axis_CUTOFF = 0.9000;
	Real const base_planarity_CUTOFF = 0.9000;
	return check_base_stack( moving_res_base, base_axis_CUTOFF, base_planarity_CUTOFF );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_BaseCentroidChecker::is_medium_base_stack_and_medium_base_pair( core::kinematics::Stub const & moving_res_base ) const{

	bool base_stack = check_base_stack( moving_res_base, 0.7070 /*base_axis_CUTOFF*/, 0.7070 /*base_planarity_CUTOFF*/ );
	bool base_pair = check_base_pair( moving_res_base, 0.5000 /*base_axis_CUTOFF*/, 0.7070 /*base_planarity_CUTOFF*/ );
	//value in Base_checker_class is 0.866 Sept 16 2010, Parin S.

	return ( base_stack && base_pair );
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_BaseCentroidChecker::update_base_stub_list_and_check_centroid_interaction( core::pose::Pose const & pose,
	StepWiseRNA_CountStruct & count_data ){
	update_base_stub_list( pose );
	return check_centroid_interaction( count_data );
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void
RNA_BaseCentroidChecker::update_base_stub_list( core::pose::Pose const & pose ){

	for ( Size const moving_res : moving_residues_ ) {
		if ( !pose.residue_type( moving_res ).is_RNA() ) continue;
		core::conformation::Residue const & residue_object( pose.residue( moving_res ) );
		Vector const centroid = rna_centroid_info_->get_base_centroid( residue_object );
		core::kinematics::Stub base_stub = rna_centroid_info_->get_base_coordinate_system( residue_object, centroid );
		base_stub_list_[ moving_res ] = base_stub;
	}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_BaseCentroidChecker::check_centroid_interaction_floating_base( core::kinematics::Stub const & moving_res_base_stub,
	StepWiseRNA_CountStruct & count_data ) const{

	bool const strong_stack_base = is_strong_base_stack( moving_res_base_stub );
	if ( strong_stack_base ) count_data.base_stack_count++;

	bool const medium_base_stack_and_medium_base_pair = is_medium_base_stack_and_medium_base_pair( moving_res_base_stub );
	if ( medium_base_stack_and_medium_base_pair ) count_data.base_pairing_count++;

	bool strict_base_pair = false;
	if ( allow_base_pair_only_screen_ ) {
		strict_base_pair = check_base_pair( moving_res_base_stub, 0.2588 /*base_axis_CUTOFF*/, 0.8660 /*base_planarity_CUTOFF*/ );
		if ( strict_base_pair ) count_data.strict_base_pairing_count++;
	}

	if ( strong_stack_base || medium_base_stack_and_medium_base_pair || ( allow_base_pair_only_screen_ && strict_base_pair ) ) {
		count_data.pass_base_centroid_screen++;
		return true;
	}

	return false;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_BaseCentroidChecker::check_centroid_interaction( core::kinematics::Stub const & moving_res_base_stub,
	StepWiseRNA_CountStruct & count_data ) {

	if ( floating_base_ ) {
		found_centroid_interaction_ = check_centroid_interaction_floating_base( moving_res_base_stub, count_data );
		return found_centroid_interaction_;
	}

	runtime_assert( moving_residues_.size() == 1 );
	base_stub_list_[ moving_residues_[ 1 ] ] = moving_res_base_stub;
	found_centroid_interaction_ = check_centroid_interaction( count_data );
	return found_centroid_interaction_;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_BaseCentroidChecker::check_centroid_interaction( StepWiseRNA_CountStruct & count_data ) {

	bool stack_base( false ), base_pair( false );

	for ( Size const res : moving_residues_ ) {

		core::kinematics::Stub const & moving_residue_base_stub = base_stub_list_[ res ];
		stack_base = false;

		for ( Size i = 1; i <= fixed_residues_.size(); i++ ) {
			stack_base = check_base_stack( moving_residue_base_stub, base_stub_list_[ fixed_residues_[ i ] ] );
			if ( stack_base ) break;
		}

		base_pair = false;

		for ( Size const fixed_res : fixed_residues_ ) {
			base_pair = check_base_pair( moving_residue_base_stub, base_stub_list_[ fixed_res ], base_pair_axis_cutoff_, base_pair_planarity_cutoff_ );
			if ( base_pair ) break;
		}

		if ( base_pair || stack_base ) break; // found an interaction!
	}

	if ( base_pair ) count_data.base_pairing_count++;
	if ( stack_base ) count_data.base_stack_count++;
	if ( base_pair || stack_base ) count_data.pass_base_centroid_screen++;

	if ( !base_pair && !stack_base ) return false;
	return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_BaseCentroidChecker::update_base_stub_list_and_check_that_terminal_res_are_unstacked( core::pose::Pose const & pose, bool const reinitialize /* = false */ ){
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
RNA_BaseCentroidChecker::check_base_stack( Size const & pos1, Size const & pos2, bool const verbose /* = false */  ) {

	if ( is_virtual_base_( pos1 ) == true || is_virtual_base_( pos2 ) == true ) {
		utility_exit_with_message( "is_virtual_base_( pos1 ) == true || is_virtual_base_( pos2 ) == true !" );
	}

	if ( pos1 == pos2 ) return true;
	return check_base_stack(  base_stub_list_[ pos1 ], base_stub_list_[ pos2 ], verbose );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
bool
RNA_BaseCentroidChecker::check_that_terminal_res_are_unstacked( bool const verbose ){

	// Look through all terminal_res
	for ( Size i = 1; i <= terminal_res_.size(); i++ ) {
		Size const & terminal_res = terminal_res_[ i ];

		for ( Size const moving_res : moving_residues_ ) {
			if ( verbose ) TR << "about to check stack: " << terminal_res << " " << moving_res << " " << stacked_on_terminal_res_in_original_pose_( terminal_res, moving_res ) << " " <<  check_base_stack( terminal_res, moving_res, false ) << std::endl;
			if ( !stacked_on_terminal_res_in_original_pose_( terminal_res, moving_res ) &&
					check_base_stack( terminal_res, moving_res, verbose  ) ) return false;
		}

		// what is this? seems gratuitous -- rhiju.
		for ( Size const fixed_res : fixed_residues_ ) {
			if ( !is_fixed_res_( fixed_res ) ) continue; // in -tether_jump condition, is_fixed_res may be 0 at fixed_res. Confusing.
			if ( verbose ) TR << "about to check stack: " << terminal_res << " " << fixed_res << " " << stacked_on_terminal_res_in_original_pose_( terminal_res, fixed_res ) << " " << check_base_stack( terminal_res, fixed_res, verbose  ) << std::endl;
			if ( !stacked_on_terminal_res_in_original_pose_( terminal_res, fixed_res ) &&
					check_base_stack( terminal_res, fixed_res, verbose  ) ) return false;
		}

	}

	if ( check_block_stack_res( block_stack_above_res_, core::chemical::rna::ABOVE ) ) return false;
	if ( check_block_stack_res( block_stack_below_res_, core::chemical::rna::BELOW ) ) return false;

	return true;
}

/////////////////////////////////////////////////////////////////////////////////
bool
RNA_BaseCentroidChecker::check_block_stack_res(
	utility::vector1< Size > const & block_stack_res,
	BaseStackWhichSide const & block_stack_side ) const
{
	for ( Size const res : block_stack_res ) {
		// if not in moving_res, look for stack with moving partitions
		utility::vector1< Size > other_partition = ( is_moving_res_( res ) ? fixed_residues_ : moving_residues_ );
		if ( check_base_stack_in_partition( res, other_partition, block_stack_side ) ) return true;
	}
	return false;
}

/////////////////////////////////////////////////////////////////////////////////
bool
RNA_BaseCentroidChecker::check_base_stack_in_partition(
	Size const & block_stack_res,
	utility::vector1< Size > const & other_res,
	BaseStackWhichSide const & block_stack_side ) const
{
	BaseStackWhichSide base_stack_side( core::chemical::rna::ANY_BASE_STACK_SIDE );
	for ( Size const res : other_res ) {
		if ( check_base_stack( base_stub_list_[ block_stack_res ],
				base_stub_list_[ res ],
				base_stack_side ) &&
				base_stack_side == block_stack_side ) {
			return true;
		}
	}
	return false;
}

} //checker
} //rna
} //modeler
} //stepwise
} //protocols
