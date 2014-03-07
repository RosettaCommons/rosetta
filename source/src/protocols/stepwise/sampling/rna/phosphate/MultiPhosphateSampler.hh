// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/rna/phosphate/MultiPhosphateSampler.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_sampling_rna_phosphate_MultiPhosphateSampler_HH
#define INCLUDED_protocols_stepwise_sampling_rna_phosphate_MultiPhosphateSampler_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/sampling/rna/phosphate/MultiPhosphateSampler.fwd.hh>
#include <protocols/stepwise/sampling/rna/phosphate/PhosphateMove.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {
namespace phosphate {

	class MultiPhosphateSampler: public protocols::moves::Mover {

	public:

		//constructor -- prepacks pose (after temporary splitting into partitions)
		MultiPhosphateSampler( pose::Pose & pose_to_prepack,
													 Size const moving_res /*sets partition*/ );

		//constructor
		MultiPhosphateSampler( pose::Pose const & reference_pose );


		//destructor
		~MultiPhosphateSampler();

	public:

		MultiPhosphateSamplerOP clone();

		void
		apply( core::pose::Pose & pose ){ copy_phosphates( pose ); }

		std::string get_name() const{ return "MultiPhosphateSampler"; }

		void
		sample_phosphates();

		void
		sample_phosphates( pose::PoseOP & viewer_pose_op );

		void
		copy_phosphates( pose::Pose & mod_pose ) const;

		void reset_to_original_pose();

		void set_screen_all( bool const & setting ){ screen_all_ = setting; }
		bool screen_all() const{ return screen_all_; }

		bool instantiated_some_phosphate() const{ return instantiated_some_phosphate_; }

		void set_moving_partition_res( utility::vector1< Size > const & setting ){ moving_partition_res_  = setting; }
		void set_five_prime_phosphate_res( utility::vector1< Size > const & setting ){ five_prime_phosphate_res_input_ = setting; }
		void set_three_prime_phosphate_res( utility::vector1< Size > const & setting ){ three_prime_phosphate_res_input_ = setting; }

		pose::Pose & pose();

		void set_scorefxn( scoring::ScoreFunctionOP scorefxn );

		utility::vector1< PhosphateMove> phosphate_move_list() const { return phosphate_move_list_; }
		void set_phosphate_move_list( utility::vector1< PhosphateMove> const & setting ) { phosphate_move_list_ = setting; }

		utility::vector1< Size > const & moving_partition_res(){ return moving_partition_res_; }

	private:

		void initialize_parameters();

		void
		initialize_by_prepack( pose::Pose & pose, Size const moving_res );

		utility::vector1< PhosphateMove >
		initialize_phosphate_move_list( pose::Pose & pose );

		utility::vector1< PhosphateMove >
		check_moved( utility::vector1< PhosphateMove > phosphate_move_list,
								 pose::Pose const & pose ) const;

		void
		find_phosphate_contacts_other_partition( utility::vector1< Size > const & partition_res1,
																						 utility::vector1< Size > const & partition_res2,
																						 pose::Pose const & pose,
																						 utility::vector1< PhosphateMove > const & phosphate_move_list,
																						 utility::vector1< PhosphateMove > & actual_phosphate_move_list ) const;

		bool
		check_other_partition_for_contact( pose::Pose const & pose,
																			 utility::vector1< Size > const & other_partition_res,
																			 Vector const & takeoff_xyz ) const;


	private:

		pose::PoseCOP pose_with_original_phosphates_;
		pose::PoseOP phosphate_sample_pose_;
		scoring::ScoreFunctionOP scorefxn_;
		Real phosphate_takeoff_donor_distance_cutoff2_;
		bool screen_all_;

		utility::vector1< Size > moving_partition_res_;
		utility::vector1 < Size > five_prime_phosphate_res_input_;
		utility::vector1 < Size > three_prime_phosphate_res_input_;

		utility::vector1< PhosphateMove> phosphate_move_list_;
		utility::vector1< PhosphateMove> actual_phosphate_move_list_;

		bool instantiated_some_phosphate_;
		bool prepacked_;
	};

} //phosphate
} //rna
} //sampling
} //stepwise
} //protocols

#endif
