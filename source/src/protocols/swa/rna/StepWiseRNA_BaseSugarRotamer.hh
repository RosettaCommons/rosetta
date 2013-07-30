// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file StepWiseRNA_BaseSugarRotamer
/// @brief
/// @detailed
/// @author Parin Sripakdeevong


#ifndef INCLUDED_protocols_swa_rna_StepWiseRNA_BaseSugarRotamer_HH
#define INCLUDED_protocols_swa_rna_StepWiseRNA_BaseSugarRotamer_HH

#include <protocols/swa/rna/StepWiseRNA_Classes.hh>
#include <protocols/swa/rna/StepWiseRNA_RotamerGenerator.fwd.hh>
//#include <protocols/swa/rna/StepWiseRNA_RotamerGenerator.hh>

#include <core/types.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <core/scoring/rna/RNA_FittedTorsionInfo.hh>

#include <string>
#include <map>

namespace protocols {
namespace swa {
namespace rna {

	class StepWiseRNA_BaseSugarRotamer: public utility::pointer::ReferenceCount {
	public:

		//constructor!
		StepWiseRNA_BaseSugarRotamer(
											BaseState const & base_state,
											PuckerState const & pucker_state,
											core::scoring::rna::RNA_FittedTorsionInfo const & rna_fitted_torsion_info,
											core::Size const bin_size=20);

    virtual ~StepWiseRNA_BaseSugarRotamer();

		void reset();

		void initialize_master_rotamer_list();

		bool get_next_rotamer();

		bool get_next_rotamer_original();

		PuckerState const & current_pucker_state() const;
		std::string const current_base_state() const;
		std::string const current_tag() const;

		core::Real const & chi()   const {return chi_;}
		core::Real const & delta() const {return delta_;}
		core::Real const & nu2() 	 const {return nu2_;}
		core::Real const & nu1() 	 const {return nu1_;}

		void set_extra_syn_chi( bool const setting){ extra_syn_chi_ =setting; }
		void set_extra_anti_chi( bool const setting){ extra_anti_chi_ =setting; }
		void set_choose_random( bool const setting){ choose_random_ =setting; }

	private:

	private:


		BaseState const base_state_;
		PuckerState const pucker_state_;
		core::scoring::rna::RNA_FittedTorsionInfo const rna_fitted_torsion_info_;
		core::Size const inputted_bin_size_; // must be 20, 10, or 5
		core::Size bin_size_;
		core::Size num_base_std_ID_;

		core::Size num_base_ID_;  //Should make this a const
		utility::vector1 < PuckerState > pucker_state_list_; //Should make this a const
		utility::vector1 < BaseState > base_state_list_; //April 30, 2011

		core::Size pucker_ID_;
		core::Size base_ID_;
		core::Size base_std_ID_;

		core::Size pucker_ID_old_;
		core::Size base_ID_old_;
		core::Size base_std_ID_old_;

		core::Real chi_;
		core::Real delta_;
		core::Real nu2_;
		core::Real nu1_;

		core::Real total_variation_;
		bool extra_anti_chi_;
		bool extra_syn_chi_;

		bool choose_random_;
		core::Size rotamer_count_;
		utility::vector1< utility::vector1< core::Size > > master_rotamer_list_;

	};
}
} //swa
} // protocols

#endif
