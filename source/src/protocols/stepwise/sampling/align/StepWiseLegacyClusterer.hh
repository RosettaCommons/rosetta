// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseLegacyClusterer.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_StepWiseLegacyClusterer_HH
#define INCLUDED_protocols_stepwise_StepWiseLegacyClusterer_HH

#include <protocols/stepwise/sampling/align/StepWiseLegacyClusterer.fwd.hh>
#include <protocols/stepwise/sampling/util.hh> // For PoseList. This is a bit silly.
#include <protocols/stepwise/sampling/modeler_options/StepWiseModelerOptions.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>

#include <map>
#include <string>

namespace protocols {
namespace stepwise {
namespace sampling {
namespace align {

  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  class StepWiseLegacyClusterer: public utility::pointer::ReferenceCount {
  public:

    //constructor!
		StepWiseLegacyClusterer( utility::vector1< PoseOP > const & pose_list );

		StepWiseLegacyClusterer( utility::vector1< PoseOP > const & pose_list,
											 utility::vector1< Size > const & moving_res_list,
											 StepWiseModelerOptionsCOP options,
											 bool const force_align );

    //destructor -- necessary?
    virtual ~StepWiseLegacyClusterer();

    /// @brief Filter a list of poses by score.

		void set_max_decoys( core::Size const & setting ){ max_decoys_ = setting; }

		void set_cluster_radius( core::Real const & setting ){ cluster_radius_ = setting; }

		void set_cluster_by_all_atom_rmsd( core::Real const & setting ){ cluster_by_all_atom_rmsd_ = setting; }

		void set_rename_tags( core::Real const & setting ){ rename_tags_ = setting; }

		void set_score_diff_cut( core::Real const & setting ){ score_diff_cut_ = setting; }

		void set_auto_tune( bool const auto_tune ){ auto_tune_ = auto_tune; }

		void set_force_align( bool const force_align ){ force_align_ = force_align; }

    void set_calc_rms_res( utility::vector1< core::Size > const & calc_rms_res );

		void cluster();

		utility::vector1< PoseOP > rms_pose_list() { return rms_pose_list_; }

		utility::vector1< PoseOP > get_pose_list() { return output_pose_list_; }

  private:

		void
		initialize_parameters_and_input();

		void
		initialize_cluster_list();

		void
		initialize_corresponding_atom_id_map( core::pose::Pose const & pose );

		void
		do_some_clustering();

		Size
		check_for_closeness( core::pose::PoseOP rms_pose );

		void
		cluster_with_auto_tune();

		void
		recluster_current_pose_list();

		void
		initialize_auto_tune_cluster_rmsds();

		core::io::silent::SilentStructOP
		setup_silent_struct( Size const n );

	private:

		Size max_decoys_;
		core::Real cluster_radius_;
		bool cluster_by_all_atom_rmsd_;
		core::Real score_diff_cut_;
		bool auto_tune_;
		bool rename_tags_;
		bool force_align_;

		core::Real score_min_;

		utility::vector1< core::pose::PoseOP > input_pose_list_;

		utility::vector1< core::pose::PoseOP > rms_pose_list_;
		utility::vector1< core::pose::PoseOP > output_pose_list_;
		utility::vector1< core::Size > num_pose_in_cluster_;

		utility::vector1< core::Real > cluster_rmsds_to_try_with_auto_tune_;

		Size input_pose_counter_;
		bool hit_score_cutoff_;
		bool initialized_atom_id_map_for_rmsd_;
    utility::vector1< core::Size > calc_rms_res_;

		std::map< core::id::AtomID, core::id::AtomID > corresponding_atom_id_map_;

  };

} //align
} //sampling
} //stepwise
} //protocols

#endif
