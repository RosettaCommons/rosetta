// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/O2PrimeScreener.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_O2PrimeScreener_HH
#define INCLUDED_protocols_stepwise_screener_O2PrimeScreener_HH

#include <protocols/stepwise/screener/SampleApplier.hh>
#include <protocols/stepwise/screener/O2PrimeScreener.fwd.hh>
#include <protocols/stepwise/sampling/rna/o2prime/O2PrimePacker.fwd.hh>

namespace protocols {
namespace stepwise {
namespace screener {

	class O2PrimeScreener: public SampleApplier {

	public:

		//constructor
		O2PrimeScreener( sampling::rna::o2prime::O2PrimePackerOP o2prime_packer );

		//destructor
		~O2PrimeScreener();

	public:

		bool
		check_screen();

		void
		add_mover( moves::CompositionMoverOP update_mover, moves::CompositionMoverOP restore_mover );

		std::string
		name() const { return "O2PrimeScreener"; }

		StepWiseScreenerType
		type() const { return O2PRIME_PACK; }

	private:

		sampling::rna::o2prime::O2PrimePackerOP o2prime_packer_;

	};

} //screener
} //stepwise
} //protocols

#endif
