// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/PhosphateScreener.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_PhosphateScreener_HH
#define INCLUDED_protocols_stepwise_screener_PhosphateScreener_HH

#include <protocols/stepwise/screener/SampleApplier.hh>
#include <protocols/stepwise/screener/PhosphateScreener.fwd.hh>
#include <protocols/stepwise/sampling/rna/phosphate/MultiPhosphateSampler.fwd.hh>

namespace protocols {
namespace stepwise {
namespace screener {

	class PhosphateScreener: public SampleApplier {

	public:

		//constructor
		PhosphateScreener( sampling::rna::phosphate::MultiPhosphateSamplerOP phosphate_sampler );

		//destructor
		~PhosphateScreener();

	public:

		std::string
		name() const { return "PhosphateScreener"; }

		StepWiseScreenerType
		type() const { return PHOSPHATE_PACK; }

		bool
		check_screen();

		void
		add_mover( moves::CompositionMoverOP update_mover, moves::CompositionMoverOP restore_mover );

	private:

		sampling::rna::phosphate::MultiPhosphateSamplerOP phosphate_sampler_;
		sampling::rna::phosphate::MultiPhosphateSamplerOP phosphate_sampler_for_restoration_;

	};

} //screener
} //stepwise
} //protocols

#endif
