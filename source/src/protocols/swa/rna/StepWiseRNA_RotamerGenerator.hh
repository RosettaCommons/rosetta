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
/// @author Rhiju Das
/// @author Parin Sripakdeevong


#ifndef INCLUDED_protocols_swa_rna_StepWiseRNA_RotamerGenerator_hh
#define INCLUDED_protocols_swa_rna_StepWiseRNA_RotamerGenerator_hh

#include <protocols/swa/rna/StepWiseRNA_Classes.hh>

#include <core/types.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <string>
#include <map>



namespace protocols {
namespace swa {
namespace rna {

	class StepWiseRNA_RotamerGenerator: public utility::pointer::ReferenceCount {
	public:

    StepWiseRNA_RotamerGenerator( Size const moving_suite,
																bool const sample_lower_sugar_and_base,
																bool const sample_upper_sugar_and_base,
																PuckerState const pucker1,
																PuckerState const pucker2 );


    virtual ~StepWiseRNA_RotamerGenerator();

		utility::vector1< core::id::TorsionID > const & torsion_ids() const;

		bool has_another_rotamer() const;

		void update_to_next_rotamer();


		utility::vector1< Torsion_Info > const & get_current_rotamer(){ return rotamer_list_;	}

		void reset();

//		void set_Is_bulge( bool const & setting ){ Is_bulge_ = setting; }

		void set_sample_extra_rotamers( bool const setting ){ sample_extra_rotamers_ = setting; }

		void set_fast( bool const setting ){ fast_ = setting; }

		void set_sample_chi_torsion( bool const setting ) { sample_chi_torsion_ = setting; }

		void set_include_syn_chi( bool const setting ){ include_syn_chi_ = setting; }

		void set_force_syn_chi_res_list( utility::vector1< core::Size > const & setting ){ force_syn_chi_res_list_ = setting; } //April 29, 2011

		void set_bin_size( core::Size const setting ){ bin_size_ = setting; }

		void set_extra_epsilon( bool const setting ){ extra_epsilon_ = setting; }

		void set_extra_beta( bool const setting ){ extra_beta_ = setting; }

		void set_extra_anti_chi( bool const setting ){ extra_anti_chi_ = setting; }

		void set_extra_syn_chi( bool const setting ){ extra_syn_chi_ = setting; }

		void set_exclude_alpha_beta_gamma_sampling(  bool const setting ){ exclude_alpha_beta_gamma_sampling_ = setting; }

		void set_allow_syn_pyrimidine(  bool const setting ){ allow_syn_pyrimidine_ = setting; }

		void initialize_rotamer_generator( core::pose::Pose const & pose );

		core::Size num_rotamer_centers();
		core::Size const & group_rotamer();
		core::Size const & subgroup_rotamer();
		core::Size const & moving_suite();

		PuckerState	pucker_state( std::string const which_sugar );

		void set_choose_random( bool const setting ){ choose_random_ = setting; }

	private:
/*
		void
		initialize_puckers(
											 core::pose::Pose const & pose,
											 bool const & sample_sugar_and_base1,
											 bool const & sample_sugar_and_base2 );
*/
		void
		initialize_sample_base_states( core::pose::Pose const & pose );

		void
		initialize_rotamers();

		void
		initialize_extra_rotamer_perturbations();

		void
		add_torsion_id( core::id::TorsionID const torsion_id );

		void
		update_to_random_rotamer();

	private:

		core::Size const moving_suite_;

		bool const sample_lower_sugar_and_base_;
		bool const sample_upper_sugar_and_base_;

		PuckerState const pucker1_specified_;
		PuckerState const pucker2_specified_;

		bool sample_extra_rotamers_;
		bool fast_;

		core::Size bin_size_;
		int bins1_, bins2_, bins3_, eps_bins_, beta_bins_; //int because they are compared to ints in torsion definition loops.
		bool verbose_;

		utility::vector1< core::id::TorsionID > torsion_ids_;
		utility::vector1< utility::vector1< core::Real > > rotamer_centers_;

		utility::vector1< core::id::TorsionID > perturb_torsion_ids_;
		utility::vector1< utility::vector1< core::Real > > extra_rotamer_perturbations_;

		core::Size group_rotamer_;
		core::Size subgroup_rotamer_;

		utility::vector1< Torsion_Info > rotamer_list_;
		bool sample_chi_torsion_;
		bool include_syn_chi_;
		utility::vector1< core::Size > force_syn_chi_res_list_; //April 29, 2011
		bool extra_epsilon_;
		bool extra_beta_;
		bool extra_anti_chi_;
		bool extra_syn_chi_;
		bool exclude_alpha_beta_gamma_sampling_;
		bool allow_syn_pyrimidine_;
		BaseState lower_base_state_;
		BaseState upper_base_state_;

		bool choose_random_;

  };

} //rna
} //swa
} // protocols

#endif

