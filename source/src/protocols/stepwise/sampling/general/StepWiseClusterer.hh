// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseClusterer.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_StepWiseClusterer_HH
#define INCLUDED_protocols_stepwise_StepWiseClusterer_HH

#include <core/pose/Pose.fwd.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <core/types.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <protocols/stepwise/StepWiseUtil.hh> // For PoseList. This is a bit silly.
#include <protocols/stepwise/sampling/general/StepWiseClusterer.fwd.hh>

#include <map>
#include <string>

namespace protocols {
namespace stepwise {
namespace sampling {
namespace general {

  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  class StepWiseClusterer: public utility::pointer::ReferenceCount {
  public:

    //constructor!
		StepWiseClusterer( utility::vector1< std::string > const & silent_files_in );

		StepWiseClusterer( std::string const & silent_file_in );

		StepWiseClusterer(  core::io::silent::SilentFileDataOP & sfd );

    //destructor -- necessary?
    virtual ~StepWiseClusterer();

    /// @brief Filter a list of poses by score.

		void set_max_decoys( core::Size const & setting ){ max_decoys_ = setting; }

		void set_cluster_radius( core::Real const & setting ){ cluster_radius_ = setting; }

		void set_cluster_by_all_atom_rmsd( core::Real const & setting ){ cluster_by_all_atom_rmsd_ = setting; }

		void set_rename_tags( core::Real const & setting ){ rename_tags_ = setting; }

		void set_score_diff_cut( core::Real const & setting ){ score_diff_cut_ = setting; }

		void set_auto_tune( bool const auto_tune ){ auto_tune_ = auto_tune; }

		void set_force_align( bool const force_align ){ force_align_ = force_align; }

    void set_calc_rms_res( utility::vector1< core::Size > const & calc_rms_res );

		void set_silent_file_data( core::io::silent::SilentFileDataOP & sfd );

		void set_rsd_type_set( std::string const & setting ){ rsd_type_set_ = setting; }

		void cluster();

		core::io::silent::SilentFileDataOP silent_file_data();

		void
		output_silent_file( std::string const & silent_file );

		PoseList
		clustered_pose_list();

		utility::vector1< core::io::silent::SilentStructOP > &
		silent_struct_output_list(){ return silent_struct_output_list_; };

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
		check_for_closeness( core::pose::PoseOP const & pose_op );

		void
		cluster_with_auto_tune();

		void
		recluster_current_pose_list();

		void
		initialize_auto_tune_cluster_rmsds();

		utility::vector1< std::string > silent_files_;
		Size max_decoys_;
		core::Real cluster_radius_;
		bool cluster_by_all_atom_rmsd_;
		core::Real score_diff_cut_;
		bool auto_tune_;
		bool rename_tags_;
		bool force_align_;
		std::string rsd_type_set_;

		core::Real score_min_;
		bool score_min_defined_;
		//		core::scoring::ScoreFunctionOP scorefxn_;

		utility::vector1< core::pose::PoseOP > pose_output_list_;
		utility::vector1< std::string > tag_output_list_;
		utility::vector1< core::io::silent::SilentStructOP > silent_struct_output_list_;
		utility::vector1< core::Size > num_pose_in_cluster_;

		core::import_pose::pose_stream::SilentFilePoseInputStreamOP input_;

		utility::vector1< core::Real > cluster_rmsds_to_try_with_auto_tune_;

		bool hit_score_cutoff_;
		bool initialized_atom_id_map_for_rmsd_;
    utility::vector1< core::Size > calc_rms_res_;

		std::map< core::id::AtomID, core::id::AtomID > corresponding_atom_id_map_;

		core::chemical::ResidueTypeSetCAP rsd_set_;

  };

} //general
} //sampling
} //stepwise
} //protocols

#endif
