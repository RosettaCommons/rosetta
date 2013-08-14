// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_RotamerGeneratorWrapper.hh
/// @brief
/// @detailed
///
/// @author Rhiju Das
/// @author Parin Sripakdeevong


#ifndef INCLUDED_protocols_swa_rna_StepWiseRNA_RotamerGeneratorWrapper_HH
#define INCLUDED_protocols_swa_rna_StepWiseRNA_RotamerGeneratorWrapper_HH

#include <protocols/swa/rna/StepWiseRNA_Classes.hh>

#include <protocols/swa/rna/StepWiseRNA_RotamerGenerator.hh>  //Feb 09, 2012: Uncomment this line. Necessary for BOINC build? (R47296 by cmiles)

#include <protocols/swa/rna/StepWiseRNA_RotamerGenerator.fwd.hh>

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

	class StepWiseRNA_RotamerGeneratorWrapper: public utility::pointer::ReferenceCount {

	public:
		StepWiseRNA_RotamerGeneratorWrapper(
								core::pose::Pose const & pose,
								utility::vector1< core::Size > const & moving_suite_list,
								bool const & sample_sugar_and_base1,
								bool const & sample_sugar_and_base2 );

	  virtual ~StepWiseRNA_RotamerGeneratorWrapper();

		bool has_another_rotamer() const;

		utility::vector1< Torsion_Info > get_next_rotamer();

		void set_sample_extra_rotamers( bool const & setting ){ sample_extra_rotamers_ = setting; }

		void set_fast( bool const & setting );

		void set_sample_chi_torsion( bool const & setting ){ sample_chi_torsion_ = setting; }

		void set_include_syn_chi( bool const include_syn_chi ) { include_syn_chi_ = include_syn_chi; }

		void set_force_syn_chi_res_list( utility::vector1< core::Size > const & setting ){ force_syn_chi_res_list_ = setting; }

		void set_force_north_ribose_list( utility::vector1< core::Size > const & setting ){ force_north_ribose_list_ = setting; }

		void set_force_south_ribose_list( utility::vector1< core::Size > const & setting ){ force_south_ribose_list_ = setting; }

		void set_bin_size( core::Size const setting ){ bin_size_ = setting; }

		void set_extra_epsilon( bool const setting ){ extra_epsilon_ = setting; }

		void set_extra_beta( bool const setting ){ extra_beta_ = setting; }

		void set_extra_anti_chi( bool const setting ){ extra_anti_chi_ = setting; }

		void set_extra_syn_chi( bool const setting ){ extra_syn_chi_ = setting; }

		void set_exclude_alpha_beta_gamma_sampling(  bool const setting ){ exclude_alpha_beta_gamma_sampling_ = setting; }

		void set_allow_syn_pyrimidine(  bool const setting ){ allow_syn_pyrimidine_ = setting; }


		core::Size group_rotamer( core::Size const list_position );
		core::Size subgroup_rotamer( core::Size const list_position );

//		core::Size group_rotamer(core::Size const list_position) { return 1; }
//		core::Size subgroup_rotamer(core::Size const list_position) { return 1; }


		core::Size rotamer_generator_list_size(){ return rotamer_generator_list_.size(); }

		void
		initialize_rotamer_generator_list();

		void set_choose_random( bool const setting );

	private:

		StepWiseRNA_RotamerGeneratorOP const
		setup_rotamer_generator( core::Size const list_position );

		core::Size
		Get_residue_pucker_state_internal( core::pose::Pose const & pose, Size list_position, std::string which_sugar, bool sample_sugar_pucker ) const;


	private:
		core::pose::Pose const & pose_;
		utility::vector1< core::Size > const moving_suite_list_;
		bool const sample_sugar_and_base1_;
		bool const sample_sugar_and_base2_;
		bool sample_extra_rotamers_;
		bool fast_;
		bool verbose_;
		bool sample_chi_torsion_;
		bool include_syn_chi_;
		bool Is_prepend_;

		utility::vector1< core::Size > force_syn_chi_res_list_;
		utility::vector1< core::Size > force_north_ribose_list_;
		utility::vector1< core::Size > force_south_ribose_list_;

		core::Size bin_size_;
		bool extra_epsilon_;
		bool extra_beta_;
		bool extra_anti_chi_;
		bool extra_syn_chi_;
		bool exclude_alpha_beta_gamma_sampling_;
		bool allow_syn_pyrimidine_;
		utility::vector1< StepWiseRNA_RotamerGeneratorOP > rotamer_generator_list_;

		bool choose_random_;

	};

}
} //swa
} // protocols

#endif

