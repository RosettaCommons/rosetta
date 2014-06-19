// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/align/StepWiseClusterer.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_sampling_align_StepWiseClusterer_HH
#define INCLUDED_protocols_stepwise_sampling_align_StepWiseClusterer_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/sampling/modeler_options/StepWiseModelerOptions.fwd.hh>
#include <protocols/stepwise/sampling/align/StepWiseClusterer.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace stepwise {
namespace sampling {
namespace align {

	class StepWiseClusterer: public utility::pointer::ReferenceCount {

	public:

		//constructor
		StepWiseClusterer();

		//constructor
		StepWiseClusterer( modeler_options::StepWiseModelerOptionsCOP options );

		//destructor
		~StepWiseClusterer();

	public:

		void cluster();

		Size size() const { return pose_list_.size(); }

		void set_max_decoys( core::Size const & setting ){ max_decoys_ = setting; }

		void set_pose_list( utility::vector1< core::pose::PoseOP > const & setting ){ pose_list_ = setting; }
		utility::vector1< core::pose::PoseOP > pose_list() { return pose_list_; }

		void set_calc_rms_res( utility::vector1< core::Size > const & setting ){ calc_rms_res_ = setting; }

		void
		apply( core::pose::Pose const & pose );

	private:

		void
		initialize_parameters( core::pose::Pose const & pose );

		void
		sort_pose_list();

		bool
		check_screen_and_kick_out_displaced_model( core::pose::Pose const & pose );

		void
		kick_out_pose_at_idx( Size const n );

		bool
		check_for_closeness( core::pose::Pose const & pose1, core::pose::Pose const & pose2 );

	private:

		utility::vector1< core::pose::PoseOP >  pose_list_;
		utility::vector1< core::Size > calc_rms_res_;
		core::Size max_decoys_;
		core::Real rmsd_;
		core::Real cluster_rmsd_;
		core::Real score_diff_cut_;
		bool initialized_;

	};

} //align
} //sampling
} //stepwise
} //protocols

#endif
