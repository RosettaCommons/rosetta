// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampler/copy_dofs/ResidueAlternativeStepWiseSamplerComb.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_sampler_copy_dofs_ResidueAlternativeStepWiseSamplerComb_HH
#define INCLUDED_protocols_sampler_copy_dofs_ResidueAlternativeStepWiseSamplerComb_HH

#include <protocols/stepwise/sampler/StepWiseSamplerSizedComb.hh>
#include <protocols/stepwise/sampler/copy_dofs/ResidueAlternativeStepWiseSamplerComb.fwd.hh>
#include <protocols/stepwise/sampler/copy_dofs/ResidueAlternativeStepWiseSampler.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

#ifdef WIN32
	#include <protocols/stepwise/sampler/copy_dofs/ResidueAlternativeStepWiseSampler.hh>
#endif


/*
using namespace core::conformation;

Commented out because “using namespace X” in header files outside of class declaration is explicitly forbidden
by our coding convention due to problems it create on modern compilers and because of the name clashing.
For more information please see: https://wiki.rosettacommons.org/index.php/Coding_conventions#Using
*/

namespace protocols {
namespace stepwise {
namespace sampler {
namespace copy_dofs {

	class ResidueAlternativeStepWiseSamplerComb: public StepWiseSamplerSizedComb {

	public:

		//constructor
    ResidueAlternativeStepWiseSamplerComb();

		//destructor
		~ResidueAlternativeStepWiseSamplerComb();

	public:

		core::conformation::Residue const &
		get_residue_at_origin( Size const seqpos );

		core::conformation::Residue const &
		get_residue_at_origin_with_matching_type( Size const seqpos, core::conformation::Residue const & rsd_in );

		void
		add_residue_alternative_rotamer( ResidueAlternativeStepWiseSamplerOP const & rotamer );

		/// @brief Name of the class
		virtual std::string get_name() const { return "ResidueAlternativeStepWiseSamplerComb"; }

		/// @brief Type of class (see enum in StepWiseSamplerTypes.hh)
		virtual StepWiseSamplerType type() const { return RESIDUE_ALTERNATIVE_COMB; }


		bool
		has_resnum( Size const seqpos );

		Size
		find_resnum( Size const seqpos );

		Size
		id_for_resnum( Size const seqpos );

		void
		fast_forward_to_next_residue_pair( Size const i, Size const j);

		void
		fast_forward_to_next_residue( Size const i );

	private:

    using StepWiseSamplerSizedComb::add_external_loop_rotamer; // make it private.

		std::map< Size, ResidueAlternativeStepWiseSamplerOP > residue_alternative_rotamer_map_;

		// following is redundant with rotamer_list_ in parent class, but is useful since it retains ResidueAlternativeStepWiseSampler type.
		utility::vector1< ResidueAlternativeStepWiseSamplerOP > residue_alternative_rotamer_list_;

	};

} //copy_dofs
} //sampler
} //stepwise
} //protocols

#endif
