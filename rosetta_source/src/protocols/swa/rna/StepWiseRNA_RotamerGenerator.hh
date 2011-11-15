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

#include <core/types.hh>
#include <core/id/TorsionID.fwd.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>


// AUTO-REMOVED #include <string>
// AUTO-REMOVED #include <map>

#ifdef WIN32
	#include <core/id/TorsionID.hh>
#endif


namespace protocols {
namespace swa {
namespace rna {

	enum PuckerState{ ALL, NORTH, SOUTH };

	class StepWiseRNA_RotamerGenerator: public utility::pointer::ReferenceCount {
	public:

		//constructor!
    StepWiseRNA_RotamerGenerator( Size const moving_suite, PuckerState const & pucker1, PuckerState const & pucker2 );

    ~StepWiseRNA_RotamerGenerator();

		utility::vector1< core::id::TorsionID > const & torsion_ids() const;

		bool has_another_rotamer() const;

		utility::vector1< core::Real > const & get_next_rotamer();

		void reset();

		void set_sample_extra_rotamers( bool const & setting ){ sample_extra_rotamers_ = setting; }
		void set_fast( bool const & setting ){ fast_ = setting; }

		core::Size const & group_rotamer();
		core::Size const & subgroup_rotamer();

	private:

		void
		initialize_rotamers();

		void
		initialize_extra_rotamer_perturbations();


		core::Size const moving_suite_;
		PuckerState const pucker1_specified_;
		PuckerState const pucker2_specified_;

		bool sample_extra_rotamers_;
		bool fast_;

		Size bin_size_;
		int bins1_, bins2_, bins3_, bins4_; //int because they are compared to ints in torsion definition loops.

		utility::vector1< core::id::TorsionID > torsion_ids_;
		utility::vector1< utility::vector1< core::Real > > rotamer_centers_;

		utility::vector1< core::id::TorsionID > perturb_torsion_ids_;
		utility::vector1< utility::vector1< core::Real > > extra_rotamer_perturbations_;

		Size group_rotamer_;
		Size subgroup_rotamer_;

		utility::vector1< core::Real > rotamer_values_;
  };

}
} //swa
} // protocols

#endif

