// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/ChainClosableGeometryResidueBasedScreener.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_ChainClosableGeometryResidueBasedScreener_HH
#define INCLUDED_protocols_stepwise_screener_ChainClosableGeometryResidueBasedScreener_HH

#include <protocols/stepwise/screener/StepWiseResiduePairScreener.hh>
#include <protocols/stepwise/screener/ChainClosableGeometryResidueBasedScreener.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/ChainClosableGeometryChecker.fwd.hh>
#include <core/conformation/Residue.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace screener {

	class ChainClosableGeometryResidueBasedScreener: public StepWiseResiduePairScreener {

	public:

		//constructor
		ChainClosableGeometryResidueBasedScreener( sampling::rna::checker::ChainClosableGeometryCheckerOP chain_closable_geometry_checker );

		//destructor
		~ChainClosableGeometryResidueBasedScreener();

	public:

		void
		get_update( rotamer_sampler::RotamerBaseOP sampler );

		bool
		check_screen();

		std::string
		name() const { return "ChainClosableGeometryResidueBasedScreener"; }

		StepWiseScreenerType
		type() const { return CHAIN_CLOSABLE_GEOMETRY_RESIDUE_BASED; }

	private:

		sampling::rna::checker::ChainClosableGeometryCheckerOP chain_closable_geometry_checker_;
//		core::conformation::ResidueCOP rsd1_, rsd2_;
Vector five_prime_xyz_, three_prime_xyz_;
	};

} //screener
} //stepwise
} //protocols

#endif
