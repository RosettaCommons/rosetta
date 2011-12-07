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


#ifndef INCLUDED_protocols_swa_SWA_RNA_RotamerGenerator_Wrapper_HH
#define INCLUDED_protocols_swa_SWA_RNA_RotamerGenerator_Wrapper_HH

#include <protocols/swa/rna/StepWiseRNA_Classes.hh>

#ifdef WIN32
#include <protocols/swa/rna/StepWiseRNA_RotamerGenerator.hh>
#endif

#include <protocols/swa/rna/StepWiseRNA_RotamerGenerator.fwd.hh>

#include <core/types.hh>
#include <core/id/TorsionID.fwd.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <string>
// AUTO-REMOVED #include <map>



namespace protocols {
namespace swa {
namespace rna {

	class StepWiseRNA_RotamerGenerator_Wrapper: public utility::pointer::ReferenceCount {

	public:
		StepWiseRNA_RotamerGenerator_Wrapper(
								core::pose::Pose const & pose,
								utility::vector1< core::Size > const & moving_suite_list,
								bool const & sample_sugar_and_base1,
								bool const & sample_sugar_and_base2,
								core::Real const bin_size_ = 20.0 );

	  ~StepWiseRNA_RotamerGenerator_Wrapper();

		bool has_another_rotamer() const;

		utility::vector1< Torsion_Info > get_next_rotamer();

		void set_sample_extra_rotamers( bool const & setting ){ sample_extra_rotamers_ = setting; }

		void set_fast( bool const & setting );

		core::Size group_rotamer(core::Size const list_position);
		core::Size subgroup_rotamer(core::Size const list_position);

//		core::Size group_rotamer(core::Size const list_position) { return 1; }
//		core::Size subgroup_rotamer(core::Size const list_position) { return 1; }


		core::Size rotamer_generator_list_size(){ return rotamer_generator_list_.size(); }

	private:

		StepWiseRNA_RotamerGeneratorOP const
		setup_rotamer_generator(core::Size const list_position);

		void
		initialize_rotamer_generator_list();

		PuckerState
		Get_residue_pucker_state_internal( core::pose::Pose const & pose, Size list_position, std::string which_sugar) const;


	private:
		core::pose::Pose const & pose_;
		utility::vector1< core::Size > const moving_suite_list_;
		bool const sample_sugar_and_base1_;
		bool const sample_sugar_and_base2_;

		bool sample_extra_rotamers_;
		bool fast_;
		bool verbose_;

		utility::vector1< StepWiseRNA_RotamerGeneratorOP > rotamer_generator_list_;

		core::Real const bin_size_;
	};

}
} //swa
} // protocols

#endif

