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


#ifndef INCLUDED_protocols_swa_SWA_RNA_RotamerGenerator_HH
#define INCLUDED_protocols_swa_SWA_RNA_RotamerGenerator_HH

#include <protocols/swa/rna/StepWiseRNA_Classes.hh>

#include <core/types.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <string>
#include <map>



#ifdef WIN32
	#include <core/id/TorsionID.hh>
#endif


namespace protocols {
namespace swa {
namespace rna {

	class StepWiseRNA_RotamerGenerator: public utility::pointer::ReferenceCount {
	public:

    StepWiseRNA_RotamerGenerator(
																	core::pose::Pose const & pose,
																	Size const moving_suite,
																	PuckerState const & pucker1,
																	PuckerState const & pucker2,
																	bool const & Is_bulge=false,
																	core::Real const binsize = 20.0
																	);


    ~StepWiseRNA_RotamerGenerator();

		utility::vector1< core::id::TorsionID > const & torsion_ids() const;

		bool has_another_rotamer() const;

		void update_to_next_rotamer();

		utility::vector1< Torsion_Info > const & get_current_rotamer(){ return rotamer_list_;	}

		void reset();

//		void set_Is_bulge( bool const & setting ){ Is_bulge_ = setting; }

		void set_sample_extra_rotamers( bool const & setting ){ sample_extra_rotamers_ = setting; }
		void set_fast( bool const & setting ){ fast_ = setting; }

		core::Size num_rotamer_centers();
		core::Size const & group_rotamer();
		core::Size const & subgroup_rotamer();
		core::Size const & moving_suite();

		PuckerState	pucker_state(std::string const which_sugar);


	private:
/*
		void
		initialize_puckers(
											 core::pose::Pose const & pose,
											 bool const & sample_sugar_and_base1,
											 bool const & sample_sugar_and_base2 );
*/
		void
		initialize_syn_chi( core::pose::Pose const & pose );

		void
		initialize_rotamers();

		void
		initialize_extra_rotamer_perturbations();

		void
		add_torsion_id(core::id::TorsionID const torsion_id);

		void
		print_pucker_state(PuckerState const & pucker_state, std::string const tag) const;

		void
		print_base_state(BaseState const & base_state, std::string const tag) const;

	private:

		core::Size const moving_suite_;
		bool const Is_bulge_;

		PuckerState pucker1_specified_;
		PuckerState pucker2_specified_;
		bool sample_syn_chi1_;
		bool sample_syn_chi2_;

		bool sample_extra_rotamers_;
		bool fast_;

		core::Size bin_size_;
		int bins1_, bins2_, bins3_, bins4_; //int because they are compared to ints in torsion definition loops.
		bool verbose_;

		utility::vector1< core::id::TorsionID > torsion_ids_;
		utility::vector1< utility::vector1< core::Real > > rotamer_centers_;

		utility::vector1< core::id::TorsionID > perturb_torsion_ids_;
		utility::vector1< utility::vector1< core::Real > > extra_rotamer_perturbations_;

		core::Size group_rotamer_;
		core::Size subgroup_rotamer_;

		utility::vector1< Torsion_Info > rotamer_list_;
  };

}
} //swa
} // protocols

#endif

