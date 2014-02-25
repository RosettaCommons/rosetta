// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/rna/phosphate/MultiPhosphateSampler.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/sampling/rna/phosphate/MultiPhosphateSampler.hh>
#include <protocols/stepwise/sampling/rna/phosphate/PhosphateMover.hh>
#include <protocols/stepwise/sampling/rna/phosphate/PhosphateUtil.hh>
#include <protocols/stepwise/StepWiseUtil.hh>
#include <core/chemical/rna/RNA_Util.hh>
#include <core/conformation/Residue.hh>
#include <core/id/TorsionID.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/pose/Pose.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/scoring/ScoreFunction.hh>
#include <utility/vector1.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.sampling.rna.phosphate.MultiPhosphateSampler" );

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {
namespace phosphate {

//Constructor
	MultiPhosphateSampler::MultiPhosphateSampler( pose::Pose const & reference_pose ):
		pose_with_original_phosphates_( reference_pose.clone() ),
		phosphate_sample_pose_( reference_pose.clone() ),
		phosphate_takeoff_donor_distance_cutoff2_( 8.0 * 8.0 ),
		screen_all_( true )
	{
		scorefxn_ = new scoring::ScoreFunction;
		scorefxn_->set_weight( scoring::fa_atr, 0.1);
		scorefxn_->set_weight( scoring::fa_rep, 0.1);
		scorefxn_->set_weight( scoring::free_suite, 1.0 );
		scorefxn_->set_weight( scoring::rna_torsion, 0.5 ); // trying to increase accepts
		scorefxn_->set_weight( scoring::hbond_lr_bb_sc, 3.0 ); // trying to increase accepts
		scorefxn_->set_weight( scoring::hbond_sr_bb_sc, 3.0 ); // trying to increase accepts
	}

	//Destructor
	MultiPhosphateSampler::~MultiPhosphateSampler()
	{}

	MultiPhosphateSamplerOP
	MultiPhosphateSampler::clone() {
		MultiPhosphateSamplerOP multi_phosphate_sampler = new MultiPhosphateSampler( *pose_with_original_phosphates_ );
		multi_phosphate_sampler->set_screen_all( screen_all_ );
		multi_phosphate_sampler->set_phosphate_move_list( phosphate_move_list_ );
		return multi_phosphate_sampler;
	}

	////////////////////////////////////////////////////////////////////
	void
	MultiPhosphateSampler::sample_phosphates(){

		phosphate_move_list_ = initialize_phosphate_move_list( *phosphate_sample_pose_ );
		copy_over_phosphate_variants( *phosphate_sample_pose_, *pose_with_original_phosphates_, phosphate_move_list_ );

		actual_phosphate_move_list_ = check_moved( phosphate_move_list_, *phosphate_sample_pose_ );

		// note that the mover will instantiate or virtualize phosphates as it sees fit.
		instantiated_some_phosphate_ = false;
		for ( Size i = 1; i <= actual_phosphate_move_list_.size(); i++ ){
			PhosphateMover phosphate_mover( actual_phosphate_move_list_[i], scorefxn_ );
			phosphate_mover.apply( *phosphate_sample_pose_ );
			if ( phosphate_mover.instantiated_phosphate() ) instantiated_some_phosphate_ = true;
		}

	}

	////////////////////////////////////////////////////////////////////
	void
	MultiPhosphateSampler::sample_phosphates( pose::PoseOP & viewer_pose_op ){
		phosphate_sample_pose_ = viewer_pose_op;
		sample_phosphates();
	}

	////////////////////////////////////////////////////////////////////
	void
	MultiPhosphateSampler::reset_to_original_pose(){
		phosphate_sample_pose_ = pose_with_original_phosphates_->clone();
	}

	////////////////////////////////////////////////////////////////////
	void
	MultiPhosphateSampler::copy_phosphates( pose::Pose & mod_pose ) const {
		copy_over_phosphate_variants( mod_pose, *phosphate_sample_pose_, phosphate_move_list_ );
	}

	////////////////////////////////////////////////////////////////////
	utility::vector1< PhosphateMove >
	MultiPhosphateSampler::initialize_phosphate_move_list( pose::Pose & pose ){

		using namespace pose::full_model_info;
		utility::vector1< PhosphateMove > phosphate_move_list;

		if ( screen_all_ ){

			FullModelInfo const & full_model_info = const_full_model_info( pose );
			utility::vector1< Size > const & res_list = full_model_info.res_list();
			utility::vector1< Size > const & cutpoint_open_in_full_model = full_model_info.cutpoint_open_in_full_model();

			Size const nres_full_model = full_model_info.size();

			for ( Size n = 1; n <= pose.total_residue(); n++ ){
				if ( !pose.residue( n ).is_RNA() ) continue;
				Size const seqpos_in_full_model = res_list[ n ];
				if ( pose.residue_type( n ).has_variant_type( "VIRTUAL_RNA_RESIDUE" ) ) continue;
				if ( ( n == 1 || pose.fold_tree().is_cutpoint( n - 1 ) ) &&
						 seqpos_in_full_model > 1 &&
						 !cutpoint_open_in_full_model.has_value( seqpos_in_full_model - 1 ) &&
						 !pose.residue_type( n ).has_variant_type( "CUTPOINT_UPPER" ) &&
						 !pose.residue_type( n ).has_variant_type( "VIRTUAL_RNA_RESIDUE_UPPER" ) ){
					phosphate_move_list.push_back( PhosphateMove( n, FIVE_PRIME_PHOSPHATE ) );
				}
				if ( ( n == pose.total_residue() || pose.fold_tree().is_cutpoint( n ) ) &&
						 seqpos_in_full_model < nres_full_model &&
						 !cutpoint_open_in_full_model.has_value( seqpos_in_full_model ) &&
						 !pose.residue_type( n ).has_variant_type( "CUTPOINT_LOWER" ) ){
					phosphate_move_list.push_back( PhosphateMove( n, THREE_PRIME_PHOSPHATE ) );
				}
			}

		} else {
			for ( Size i = 1; i <= five_prime_phosphate_res_input_.size(); i++ ){
				phosphate_move_list.push_back(  PhosphateMove( five_prime_phosphate_res_input_[i], FIVE_PRIME_PHOSPHATE ) );
			}
			for ( Size i = 1; i <= three_prime_phosphate_res_input_.size(); i++ ){
				phosphate_move_list.push_back(  PhosphateMove( three_prime_phosphate_res_input_[i], THREE_PRIME_PHOSPHATE ) );
			}
		}

		return phosphate_move_list;
	}

	///////////////////////////////////////////////////////////////////////////////////////
	utility::vector1< PhosphateMove >
	MultiPhosphateSampler::check_moved( utility::vector1< PhosphateMove > phosphate_move_list,
																			pose::Pose const & pose ) const {

		utility::vector1< PhosphateMove > actual_phosphate_move_list;

		if ( moving_partition_res_.size() == 0 ) return phosphate_move_list;

		utility::vector1< Size > partition_res1, partition_res2;
		for ( Size n = 1; n <= pose.total_residue(); n++ ){
			if ( !pose.residue( n ).is_RNA() ) continue;
			if ( moving_partition_res_.has_value( n ) ){
				partition_res1.push_back( n );
			} else {
				partition_res2.push_back( n );
			}
		}

		// check phosphates that make cross-partition contacts.
		find_phosphate_contacts_other_partition( partition_res1, partition_res2, pose,
																						phosphate_move_list, actual_phosphate_move_list );
		find_phosphate_contacts_other_partition( partition_res2, partition_res1, pose,
																						phosphate_move_list, actual_phosphate_move_list );

		// need to also add in any phosphates that made cross-contact partitions in *original* pose.
		// those need to be re-evaluated.
		find_phosphate_contacts_other_partition( partition_res1, partition_res2, *pose_with_original_phosphates_,
																						phosphate_move_list, actual_phosphate_move_list );
		find_phosphate_contacts_other_partition( partition_res2, partition_res1, *pose_with_original_phosphates_,
																						phosphate_move_list, actual_phosphate_move_list );

		return actual_phosphate_move_list;
	}

	///////////////////////////////////////////////////////////////////////////////////////
	void
	MultiPhosphateSampler::find_phosphate_contacts_other_partition( utility::vector1< Size > const & partition_res1,
																																	utility::vector1< Size > const & partition_res2,
																																	pose::Pose const & pose,
																																	utility::vector1< PhosphateMove > const & phosphate_move_list,
																																	utility::vector1< PhosphateMove > & actual_phosphate_move_list ) const {

		for ( Size i = 1; i <= phosphate_move_list.size(); i++ ){
			PhosphateMove const & phosphate_move = phosphate_move_list[ i ];
			if ( actual_phosphate_move_list.has_value( phosphate_move ) ) continue;
			Size const n = phosphate_move.rsd();
			if ( !partition_res1.has_value( n ) ) continue;
			Vector const phosphate_takeoff_xyz = (phosphate_move.terminus() == FIVE_PRIME_PHOSPHATE) ?
				pose.residue( n ).xyz( " C5'" ) :
				pose.residue( n ).xyz( " O3'" );
			if ( !check_other_partition_for_contact( pose, partition_res2, phosphate_takeoff_xyz ) ) continue;
			actual_phosphate_move_list.push_back( phosphate_move );
		}
	}

	////////////////////////////////////////////////////////////////////
	bool
	MultiPhosphateSampler::check_other_partition_for_contact( pose::Pose const & pose,
																														utility::vector1< Size > const & other_partition_res,
																														Vector const & takeoff_xyz ) const {

		for ( Size i = 1; i <= other_partition_res.size(); i++ ){

			Size const & m = other_partition_res[ i ];
			core::chemical::AtomIndices Hpos_polar = pose.residue_type( m ).Hpos_polar();

			for ( Size ii = 1; ii <= Hpos_polar.size(); ii++ ){
				Size const & q = Hpos_polar[ ii ];
				if ( ( pose.residue( m ).xyz( q ) - takeoff_xyz ).length_squared() < phosphate_takeoff_donor_distance_cutoff2_ ){
					return true;
				}
			}
		}

		return false;
	}

	////////////////////////////////////////////////////////////////////////
	pose::Pose &
	MultiPhosphateSampler::pose(){ return *phosphate_sample_pose_; }

	////////////////////////////////////////////////////////////////////////
	void
	MultiPhosphateSampler::set_scorefxn( scoring::ScoreFunctionOP scorefxn ){
		scorefxn_ = scorefxn;
	}


} //phosphate
} //rna
} //sampling
} //stepwise
} //protocols
