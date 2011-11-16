// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_Clusterer.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_swa_StepWiseClusterer_hh
#define INCLUDED_protocols_swa_StepWiseClusterer_hh

#include <core/pose/Pose.fwd.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <utility/vector1.hh>
#include <core/types.hh>
#include <core/io/silent/SilentStruct.fwd.hh>

#include <map>
#include <string>

namespace protocols {
namespace swa {

	typedef std::map< std::string, core::pose::PoseOP > PoseList;

  /////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////
  class StepWiseClusterer: public utility::pointer::ReferenceCount {
  public:

    //constructor!
		StepWiseClusterer( utility::vector1< std::string > const & silent_files_in );

		StepWiseClusterer( std::string const & silent_file_in );

		StepWiseClusterer(  core::io::silent::SilentFileDataOP & sfd );

    //destructor -- necessary?
    ~StepWiseClusterer();

    /// @brief Filter a list of poses by score.

		void set_max_decoys( core::Size const & setting ){ max_decoys_ = setting; }

		void set_cluster_radius( core::Real const & setting ){ cluster_radius_ = setting; }

		void set_cluster_by_all_atom_rmsd( core::Real const & setting ){ cluster_by_all_atom_rmsd_ = setting; }

		void set_rename_tags( core::Real const & setting ){ rename_tags_ = setting; }

		void set_score_diff_cut( core::Real const & setting ){ score_diff_cut_ = setting; }

		void set_rsd_type_set( std::string const & setting ){ rsd_type_set_ = setting; }

		void cluster();

		void
		output_silent_file( std::string const & silent_file );

		PoseList
		clustered_pose_list();

  private:

		bool
		check_for_closeness( core::pose::PoseOP const & pose_op );


		utility::vector1< std::string > silent_files_;
		Size max_decoys_;
		core::Real cluster_radius_;
		bool cluster_by_all_atom_rmsd_;
		core::Real score_diff_cut_;
		bool rename_tags_;
		std::string rsd_type_set_;

		//		core::scoring::ScoreFunctionOP scorefxn_;

		utility::vector1< core::pose::PoseOP > pose_output_list_;
		utility::vector1< std::string > tag_output_list_;
		utility::vector1< core::io::silent::SilentStructOP > silent_struct_output_list_;

		core::import_pose::pose_stream::SilentFilePoseInputStreamOP input_;

  };

} //swa
} // protocols

#endif
