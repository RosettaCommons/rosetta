// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/recces/sampler/MC_Loop.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_recces_sampler_MC_Loop_HH
#define INCLUDED_protocols_recces_sampler_MC_Loop_HH

#include <protocols/recces/sampler/MC_Loop.fwd.hh>
#include <protocols/recces/sampler/MC_Any.hh>

namespace protocols {
namespace recces {
namespace sampler {

	class MC_Loop: public MC_Any {

	public:

		//constructor
		MC_Loop();

		//destructor
		~MC_Loop();

	public:

		/// @brief Move to next rotamer
		virtual void operator++();

		/// @brief Type of class (see enum in toolbox::SamplerPlusPlusTypes.hh)
		virtual toolbox::SamplerPlusPlusType type() const { return toolbox::MC_LOOP; }

	};

} //sampler
} //recces
} //protocols

#endif
