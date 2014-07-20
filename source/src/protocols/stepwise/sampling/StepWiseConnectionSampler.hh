// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/StepWiseConnectionSampler.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_sampling_StepWiseConnectionSampler_HH
#define INCLUDED_protocols_stepwise_sampling_StepWiseConnectionSampler_HH

#include <protocols/moves/MoverForPoseList.hh>
#include <protocols/stepwise/sampling/align/StepWiseClusterer.hh>
#include <protocols/stepwise/sampling/modeler_options/StepWiseModelerOptions.fwd.hh>
#include <protocols/stepwise/sampling/protein/loop_close/StepWiseProteinCCD_Closer.fwd.hh>
#include <protocols/stepwise/sampling/protein/InputStreamWithResidueInfo.fwd.hh>
#include <protocols/stepwise/sampling/packer/StepWiseMasterPacker.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_Classes.hh>
#include <protocols/stepwise/sampling/rna/rigid_body/FloatingBaseClasses.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_BaseCentroidChecker.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_VDW_BinChecker.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_AtrRepChecker.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_ChainClosureChecker.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_ChainClosableGeometryChecker.fwd.hh>
#include <protocols/stepwise/sampling/rna/o2prime/O2PrimePacker.fwd.hh>
#include <protocols/stepwise/sampling/rna/phosphate/MultiPhosphateSampler.fwd.hh>
#include <protocols/stepwise/sampling/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/screener/StepWiseScreener.fwd.hh>
#include <protocols/stepwise/screener/TagDefinition.fwd.hh>
#include <protocols/rotamer_sampler/copy_dofs/ResidueAlternativeRotamerSamplerComb.fwd.hh>
#include <protocols/rotamer_sampler/copy_dofs/ResidueAlternativeSet.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamerSampler.fwd.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamerSamplerWithResidueList.fwd.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamerSamplerWithResidueAlternatives.fwd.hh>
#include <protocols/rotamer_sampler/RotamerSamplerSized.fwd.hh>
#include <protocols/rotamer_sampler/RotamerSamplerBase.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <protocols/stepwise/sampling/StepWiseConnectionSampler.fwd.hh>

#ifdef PYROSETTA
	#include <core/scoring/ScoreFunction.hh>
	#include <protocols/stepwise/sampling/packer/StepWiseMasterPacker.hh>
	#include <protocols/stepwise/sampling/protein/InputStreamWithResidueInfo.hh>
#endif

namespace protocols {
namespace stepwise {
namespace sampling {

	class StepWiseConnectionSampler: public protocols::moves::MoverForPoseList {

	private:

		StepWiseConnectionSampler(); // Can't use a default constructor

	public:

		//constructor
		StepWiseConnectionSampler( working_parameters::StepWiseWorkingParametersCOP & working_parameters_ );

		//destructor
		~StepWiseConnectionSampler();

		virtual void apply( core::pose::Pose & pose_to_visualize );

		virtual std::string get_name() const;

		using MoverForPoseList::apply;

	public:

		void set_scorefxn( core::scoring::ScoreFunctionCOP const & scorefxn ){ scorefxn_ = scorefxn; }

		void set_pose_list( utility::vector1< core::pose::PoseOP > &	pose_list );

		utility::vector1< core::pose::PoseOP > & get_pose_list();

		void
		set_user_input_VDW_bin_checker( rna::checker::RNA_VDW_BinCheckerOP const & user_input_VDW_bin_checker );

		void
		set_master_packer( packer::StepWiseMasterPackerOP master_packer ){ master_packer_ = master_packer; }

		void
		add_residue_alternative_set( rotamer_sampler::copy_dofs::ResidueAlternativeSet const & residue_alternative_set );

		void
		set_options( modeler_options::StepWiseModelerOptionsCOP options );

		void
		set_input_streams( utility::vector1< protein::InputStreamWithResidueInfoOP > const & input_streams ){ input_streams_ = input_streams; }

	private:

		void
		figure_out_reference_res( core::pose::Pose const & pose );

		void
		figure_out_reference_res_with_rigid_body_rotamer( core::pose::Pose const & pose );

		void
		initialize_useful_info( core::pose::Pose const & pose );

		bool
		initialize_pose( core::pose::Pose & pose );

		void
		initialize_checkers( core::pose::Pose const & pose );

		core::Size
		get_max_ntries();

		core::Size
		get_num_pose_kept();

		void
		initialize_euler_angle_grid_parameters();

		void
		initialize_xyz_grid_parameters();

		bool
		initialize_sampler( core::pose::Pose const & pose );

		rotamer_sampler::RotamerSamplerSizedOP
		initialize_protein_bond_sampler( core::pose::Pose const & pose );

		rotamer_sampler::RotamerSamplerBaseOP
		initialize_rna_bond_sampler( core::pose::Pose const & pose );

		void
		initialize_full_rigid_body_sampler();

		rotamer_sampler::RotamerSamplerBaseOP
		get_full_bond_sampler();

		rotamer_sampler::copy_dofs::ResidueAlternativeRotamerSamplerCombOP
		get_rsd_alternatives_rotamer();

		void
		initialize_rigid_body_rotamer();

		void
		initialize_screeners( core::pose::Pose & pose );

		void
		initialize_residue_level_screeners( core::pose::Pose & pose );

		void
		initialize_pose_level_screeners( core::pose::Pose & pose );

		void
		initialize_moving_residue_pose_list( core::pose::Pose const & pose );

		void
		initialize_protein_packer( core::pose::Pose & pose );

		core::Size
		which_residue_alternative_set_is_moving_residue() const;

		utility::vector1< core::conformation::ResidueOP >	get_moving_rsd_list() const;

		core::Size truly_floating_base();

		void check_working_parameters( core::pose::Pose const & pose );

		bool
		presample_virtual_sugars( core::pose::Pose & pose );

	private:

		working_parameters::StepWiseWorkingParametersCOP working_parameters_;
		utility::vector1< core::Size > const moving_res_list_;
		core::Size const moving_res_;
		utility::vector1< core::Size > const moving_partition_res_;
		core::scoring::ScoreFunctionCOP scorefxn_;
		modeler_options::StepWiseModelerOptionsCOP options_;

		utility::vector1< core::pose::PoseOP > pose_list_;

		bool rigid_body_sampling_;
		core::Size reference_res_;
		bool kic_sampling_;
		bool protein_connection_; // should be able to deprecate soon

		rotamer_sampler::rigid_body::RigidBodyRotamerSamplerOP rigid_body_rotamer_;
		protocols::rotamer_sampler::RotamerSamplerBaseOP sampler_;
		utility::vector1< screener::StepWiseScreenerOP > screeners_;
		screener::TagDefinitionOP tag_definition_;

		utility::vector1< rotamer_sampler::copy_dofs::ResidueAlternativeSet > residue_alternative_sets_;

		packer::StepWiseMasterPackerOP master_packer_;

		// atr/rep checks
		core::pose::PoseOP protein_atr_rep_screening_pose_;
		rna::checker::RNA_AtrRepCheckerOP atr_rep_checker_;
		rna::checker::RNA_AtrRepCheckerOP virt_sugar_atr_rep_checker_;
		rna::checker::RNA_VDW_BinCheckerOP VDW_bin_checker_;
		rna::checker::RNA_VDW_BinCheckerOP user_input_VDW_bin_checker_;

		// chain closure.
		utility::vector1< core::Size > protein_cutpoints_closed_;
		utility::vector1< core::pose::PoseOP > protein_ccd_poses_;
		utility::vector1< protein::loop_close::StepWiseProteinCCD_CloserOP > protein_ccd_closers_;
		utility::vector1< core::Size > rna_five_prime_chain_breaks_;
		utility::vector1< core::Size > rna_three_prime_chain_breaks_;
		utility::vector1< core::Size > rna_chain_break_gap_sizes_;
		utility::vector1< rna::checker::RNA_ChainClosableGeometryCheckerOP > rna_chain_closable_geometry_checkers_;
		utility::vector1< core::Size > rna_cutpoints_closed_;
		utility::vector1< rna::checker::RNA_ChainClosureCheckerOP > rna_chain_closure_checkers_;

		rna::checker::RNA_BaseCentroidCheckerOP base_centroid_checker_;
		core::pose::PoseOP screening_pose_, virt_sugar_screening_pose_;

		align:: StepWiseClustererOP clusterer_;

		core::kinematics::Stub moving_res_base_stub_;

		core::Real max_distance_squared_;

		utility::vector1< protein::InputStreamWithResidueInfoOP > input_streams_; // this is now awkward.
		bool const virt_sugar_atr_rep_screen_; // deprecate?
		rna::StepWiseRNA_CountStruct count_data_;

	};

} //sampling
} //stepwise
} //protocols

#endif
