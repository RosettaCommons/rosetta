// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_LoopBridger.hh
/// @brief
/// @details
///
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_protein_StepWiseProteinKIC_LoopBridger_HH
#define INCLUDED_protocols_stepwise_protein_StepWiseProteinKIC_LoopBridger_HH

#include <core/pose/Pose.fwd.hh>
#include <core/scoring/Ramachandran.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/sampler/StepWiseSamplerSized.fwd.hh>
#include <protocols/moves/Mover.hh>
#include <protocols/loops/Loop.hh>
#include <core/id/TorsionID.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace stepwise {
namespace modeler {
namespace protein {
namespace loop_close {

  /////////////////////////////////////////////////////////////////////////////////////////////////
  class StepWiseProteinKIC_LoopBridger: public protocols::moves::Mover {
  public:

    //constructor!
		StepWiseProteinKIC_LoopBridger( sampler::StepWiseSamplerSizedOP sampler,
																protocols::stepwise::modeler::working_parameters::StepWiseWorkingParametersCOP working_parameters );

    //destructor -- necessary?
    ~StepWiseProteinKIC_LoopBridger();

    /// @brief Apply the minimizer to one pose
    virtual void apply( core::pose::Pose & pose_to_visualize );

	virtual std::string get_name() const;


		utility::vector1< utility::vector1< core::Real > > const & main_chain_torsion_set_lists() const;
		utility::vector1< core::id::TorsionID > const & which_torsions() const;

	private:

		void
		setup_torsions( core::pose::Pose const & pose );

		void
		figure_out_loop( core::pose::Pose const & pose );

		void
		output_chainTORS( utility::vector1< core::Real > const & dt_ang,
											utility::vector1< core::Real > const & db_ang,
											utility::vector1< core::Real > const & db_len ) const;

		void
		fill_chainTORS_info( core::pose::Pose const & pose,
												 utility::vector1<utility::vector1<core::Real> > & atoms,
												 utility::vector1<core::Real> & dt_ang,
												 utility::vector1<core::Real> & db_ang,
												 utility::vector1<core::Real> & db_len,
												 core::Size const & start_res_ ,
												 core::Size const & end_res_ ) const;

		void
		KIC_loop_close_with_perturbations( core::pose::Pose & pose );

		void
		KIC_loop_close( core::pose::Pose & pose );

		void
		grab_main_chain_torsion_set_list( core::pose::Pose const & pose );

		void save_phi_psi_offsets( core::pose::Pose const & pose );
		void fix_phi_psi_offsets( core::pose::Pose & pose ) const;

		void
		sample_omega_recursively( core::pose::Pose & pose, int const offset, utility::vector1<utility::vector1< core::Real> > & atoms, utility::vector1< core::Real> & dt_ang, utility::vector1< core::Real> & db_ang, utility::vector1< core::Real> & db_len, utility::vector1< core::Size > const & pivots, utility::vector1< core::Size > const & order );

		void
		initialize_is_fixed_res( utility::vector1< core::Size > const & fixed_res, std::string const & working_sequence );

	private:

		sampler::StepWiseSamplerSizedOP sampler_;

		utility::vector1< Size > working_bridge_res_;
		utility::vector1< bool > is_pre_proline_;
		utility::vector1< bool > is_fixed_res_;

		loops::Loop loop_;

		utility::vector1< core::id::TorsionID >  which_torsions_;

		utility::vector1< utility::vector1< core::Real > > main_chain_torsion_sets_for_moving_residues_;

		core::scoring::Ramachandran ramachandran_;

		int const num_perturb_steps_;
		core::Real const perturb_torsion_;

		core::Real const idl_CA_C_N_, idl_C_N_CA_, idl_C_N_, OMEGA_MEAN_;
		bool const use_icoor_geometry_;

		bool const verbose_;

		utility::vector1< core::Real > phi_offsets_;
		utility::vector1< core::Real > psi_offsets_;

		Size start_res_, middle_res_, end_res_, middle_offset_, seg_len_;
		Size pose_count_;

  };

} //loop_close
} //protein
} //modeler
} //stepwise
} //protocols

#endif
