// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/sampler/rna/MC_RNA_OneJump.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_recces_sampler_rna_MC_RNA_OneJump_HH
#define INCLUDED_protocols_recces_sampler_rna_MC_RNA_OneJump_HH

// Unit headers
#include <protocols/recces/sampler/rna/MC_RNA_OneJump.fwd.hh>

// Package headers
#include <protocols/recces/sampler/MC_Sampler.hh>

// Project headers
#include <core/kinematics/Jump.hh>


namespace protocols {
namespace recces {
namespace sampler {
namespace rna {

	class MC_RNA_OneJump: public MC_Sampler {

	public:

		//constructor
		MC_RNA_OneJump( core::pose::Pose const & pose,
																		core::Size const & jump_num );

		//destructor
		~MC_RNA_OneJump();

	/// @brief Initialization
		void init() {
			set_init( true );
			reset();
		}

		/// @brief Reset to current angle
		void reset() { active_jump_ = stored_jump_; }

		/// @brief Update the active jump based on stored
		/// (do not update stored_jump_)
		void operator++();

		/// @brief Update the stored jump to match active jump
		void update();

		/// @brief Apply the active jump to pose
		void apply( core::pose::Pose & pose );

		/// @brief Get the stored jump
		core::kinematics::Jump stored_jump() const { return stored_jump_; }

		/// @brief Get the active jump
		core::kinematics::Jump active_jump() const { return active_jump_; }

		/// @brief Set the stored jump
		void set_jump( core::kinematics::Jump const & setting ) { stored_jump_ = setting; }

		/// @brief Set the jump range (defined by RMSD of base from reference location )
		void set_rmsd_cutoff( core::Real const & setting ) { rmsd_cutoff_ = setting; }

		/// @brief Set the standard deviation of translation -- LATER: split out Gaussian vs. uniform
		void set_translation_mag( core::Real const & setting ) {
			translation_mag_ = setting;
		}

		/// @brief Set the standard deviation of rotation -- LATER: split out Gaussian vs. uniform
		void set_rotation_mag( core::Real const & setting ) {
			rotation_mag_ = setting;
		}

		/// @brief Name of the class
		std::string get_name() const { return "MC_RNA_OneJump"; }

		/// @brief Type of class (see enum in toolbox::SamplerPlusPlusTypes.hh)
		virtual toolbox::SamplerPlusPlusType type() const { return toolbox::MC_RNA_ONE_JUMP; }

		/// @brief output summary of class
		virtual
		void show( std::ostream & out, Size const indent) const;

		/// @brief return OP to the subsampler that controls exactly this torsion_id (assume only one).
		virtual
		MC_SamplerOP
		find( core::id::TorsionID const & torsion_id );

	private:

	  bool check_jump_in_range();

		Size const jump_num_;

		core::kinematics::Jump stored_jump_, active_jump_;

		// user defined parameters
		core::Real translation_mag_, rotation_mag_, rmsd_cutoff_;

		// scratchwork to check if Jump's are within range
		core::pose::PoseOP scratch_pose_;
		core::kinematics::Stub stored_upstream_stub_;
		core::Vector stored_base_centroid_;
		// original pose that serves to define reference geometry.
		core::pose::PoseCOP ref_scratch_pose_;


	};

} //rna
} //sampler
} //recces
} //protocols

#endif
