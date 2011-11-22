// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_ResidueSampler.hh
/// @brief
/// @detailed
///
/// @author Parin Sripakdeevong


#ifndef INCLUDED_protocols_swa_SWA_RNA_Clusterer_HH
#define INCLUDED_protocols_swa_SWA_RNA_Clusterer_HH

#include <string>
#include <map>

#include <protocols/swa/rna/StepWiseRNA_Util.hh>
#include <protocols/swa/rna/StepWiseRNA_ResidueInfo.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
// AUTO-REMOVED #include <core/io/silent/SilentFileData.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace swa {
namespace rna {

	class StepWiseRNA_Clusterer: public utility::pointer::ReferenceCount {
  	public:

		//constructor!
    StepWiseRNA_Clusterer(utility::vector1<std::string> const & input_silent_file_list);

	  //destructor -- necessary?
    ~StepWiseRNA_Clusterer();

		utility::vector1< pose_data_struct2 >
		cluster();

	  void
  	set_output_silent_file( std::string const & output_silent_file);

		void
		set_cluster_mode( std::string const & cluster_mode);

		// Undefinded, commenting out to fix PyRosetta build  void set_input_pose_data_list(utility::vector1 <pose_data_struct2> const & input_pose_data_list);

		void
		create_cluster_residue_list(utility::vector1< core::Size > const & cluster_res_seq_num_list,
																std::string const & full_fasta_sequence);

		void
		set_is_prepend_map(std::map< core::Size, bool > const & Is_prepend_map);

		void
		set_res_map(std::map< core::Size, core::Size > const & res_map);

		private:

		bool
		Individual_residue_clustering(core::pose::Pose const & current_pose, core::pose::Pose const & cluster_center_pose) const;

		bool
		Is_new_cluster_center(pose_data_struct2 const & current_pose_data, utility::vector1< pose_data_struct2> const & cluster_pose_data_list) const;

		void
		Check_residue_list_parameters() const;


  	private:

		utility::vector1 <std::string> const input_silent_file_list_;


		core::Size const num_pose_kept_;
		core::Real const cluster_rmsd_;
		std::string cluster_mode_;
		std::string output_silent_file_;



		std::map< core::Size, bool > Is_prepend_map_; //full_res_to_is_prepend_map
		std::map< core::Size, core::Size > res_map_; //full_res_to_partial_res_map_
		utility::vector1< Residue_info > cluster_residue_list_;

//		utility::vector1< Residue_info > full_residue_list_;

		bool const verbose_;


  };

}
} //swa
} // protocols

#endif

