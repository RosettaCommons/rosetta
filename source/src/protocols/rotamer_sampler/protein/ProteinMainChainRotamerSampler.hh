// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/protein/ProteinMainChainRotamerSampler.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rotamer_sampler_protein_ProteinMainChainRotamerSampler_HH
#define INCLUDED_protocols_rotamer_sampler_protein_ProteinMainChainRotamerSampler_HH

#include <protocols/rotamer_sampler/RotamerSamplerSized.hh>
#include <protocols/rotamer_sampler/protein/ProteinMainChainRotamerSampler.fwd.hh>
#include <core/types.hh>
#include <core/id/TorsionID.fwd.hh>
#include <utility/vector1.hh>

namespace protocols {
namespace rotamer_sampler {
namespace protein {

	class ProteinMainChainRotamerSampler: public protocols::rotamer_sampler::RotamerSamplerSized {

	public:

		//constructor
		ProteinMainChainRotamerSampler( utility::vector1< core::id::TorsionID > const & which_torsions,
														 utility::vector1< utility::vector1< core::Real > > const & main_chain_torsion_set_lists,
														 bool const choose_random = false );

		ProteinMainChainRotamerSampler();

		//destructor
		~ProteinMainChainRotamerSampler();

		/// @brief Get the total number of rotamers in sampler
		virtual core::Size size() const { return main_chain_torsion_set_lists_.size(); }

		/// @brief Apply the i-th rotamer to pose
		virtual void apply( core::pose::Pose&, core::Size const );

		/// @brief Name of the class
		virtual std::string get_name() const { return "ProteinMainChainRotamerSampler"; }

		/// @brief Type of class (see enum in RotamerSamplerTypes.hh)
		virtual RotamerSamplerType type() const { return PROTEIN_MAIN_CHAIN; }

	private:

		utility::vector1< core::id::TorsionID > const which_torsions_;
		utility::vector1< utility::vector1< core::Real > > const main_chain_torsion_set_lists_;

	};

} //protein
} //rotamer_sampler
} //protocols

#endif
