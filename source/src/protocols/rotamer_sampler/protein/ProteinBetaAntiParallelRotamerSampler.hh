// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/protein/ProteinBetaAntiParallelRotamerSampler.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_rotamer_sampler_protein_ProteinBetaAntiParallelRotamerSampler_HH
#define INCLUDED_protocols_rotamer_sampler_protein_ProteinBetaAntiParallelRotamerSampler_HH

#include <protocols/rotamer_sampler/JumpRotamerSampler.hh>
#include <protocols/rotamer_sampler/protein/ProteinBetaAntiParallelRotamerSampler.fwd.hh>

using namespace core;

namespace protocols {
namespace rotamer_sampler {
namespace protein {

	class ProteinBetaAntiParallelRotamerSampler: public JumpRotamerSampler {

	public:

		//constructor
		ProteinBetaAntiParallelRotamerSampler( core::pose::Pose const & pose,
																		Size const moving_residue );

		//destructor
		~ProteinBetaAntiParallelRotamerSampler();

	public:

		/// @brief Name of the class
		virtual std::string get_name() const { return "ProteinBetaAntiParallelRotamerSampler"; }

		/// @brief Type of class (see enum in RotamerSamplerTypes.hh)
		virtual RotamerSamplerType type() const { return PROTEIN_BETA_ANTIPARALLEL; }

	private:

		Size
		get_antiparallel_beta_jumps( pose::Pose const & pose, int const sample_res );

	};

} //protein
} //rotamer_sampler
} //protocols

#endif
