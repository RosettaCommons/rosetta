// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/swa/rna/StepWiseRNA_PoseSelection.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_swa_rna_StepWiseRNA_PoseSelection_HH
#define INCLUDED_protocols_swa_rna_StepWiseRNA_PoseSelection_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/swa/rna/StepWiseRNA_PoseSelection.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/swa/rna/StepWiseRNA_Classes.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

using namespace core;

namespace protocols {
namespace swa {
namespace rna {

	class StepWiseRNA_PoseSelection: public utility::pointer::ReferenceCount {

	public:

	//constructor
		StepWiseRNA_PoseSelection( StepWiseRNA_JobParametersCOP & job_parameters,
															 core::scoring::ScoreFunctionCOP scorefxn );

		//destructor
		~StepWiseRNA_PoseSelection();

	public:

		void
		initialize_sampling_scorefxn( core::scoring::ScoreFunctionCOP & scorefxn );

		Real
		pose_selection_by_full_score( pose::Pose & current_pose, std::string const & tag );

		void
		cluster_pose_list();

		void set_num_pose_kept( core::Size const & num_pose_kept );

		void
		set_cluster_rmsd( core::Real const & setting );

		void
		set_distinguish_pucker( bool const & setting ){ distinguish_pucker_ = setting ; }

		void
		set_PBP_clustering_at_chain_closure( bool const & setting ){ PBP_clustering_at_chain_closure_ = setting; }

		void
		update_pose_list(
													std::string const & tag,
													pose::Pose const & current_pose,
													Real const & current_score );

		utility::vector1< pose::PoseOP > pose_list();

		void set_pose_list( utility::vector1< pose::PoseOP > &	pose_list );

		void
		finalize( bool const do_clustering = true );

	private:

		StepWiseRNA_JobParametersCOP job_parameters_;
		core::scoring::ScoreFunctionOP sampling_scorefxn_;

		Size num_pose_kept_;
		Size const multiplier_;
		core::Real cluster_rmsd_;
		bool PBP_clustering_at_chain_closure_;
		bool distinguish_pucker_;
		core::Real current_score_cutoff_;

		bool verbose_;

		StepWiseRNA_CountStruct count_data_;

		utility::vector1< pose::PoseOP > pose_list_;

	};

} //rna
} //swa
} //protocols

#endif
