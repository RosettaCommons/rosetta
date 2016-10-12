// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/screener/RNA_ChainClosableGeometryResidueBasedScreener.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_RNA_ChainClosableGeometryResidueBasedScreener_HH
#define INCLUDED_protocols_stepwise_screener_RNA_ChainClosableGeometryResidueBasedScreener_HH

#include <protocols/stepwise/screener/StepWiseResiduePairScreener.hh>
#include <protocols/stepwise/screener/RNA_ChainClosableGeometryResidueBasedScreener.fwd.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_ChainClosableGeometryChecker.fwd.hh>
#include <core/conformation/Residue.hh>

namespace protocols {
namespace stepwise {
namespace screener {

class RNA_ChainClosableGeometryResidueBasedScreener: public StepWiseResiduePairScreener {

public:

	//constructor
	RNA_ChainClosableGeometryResidueBasedScreener( modeler::rna::checker::RNA_ChainClosableGeometryCheckerOP chain_closable_geometry_checker );

	//destructor
	~RNA_ChainClosableGeometryResidueBasedScreener();

public:

	void
	get_update( sampler::StepWiseSamplerBaseOP sampler );

	bool
	check_screen();

	std::string
	name() const { return "RNA_ChainClosableGeometryResidueBasedScreener"; }

	StepWiseScreenerType
	type() const { return RNA_CHAIN_CLOSABLE_GEOMETRY_RESIDUE_BASED; }

private:

	modeler::rna::checker::RNA_ChainClosableGeometryCheckerOP chain_closable_geometry_checker_;
	//  core::conformation::ResidueCOP rsd1_, rsd2_;
	core::Vector five_prime_xyz_, three_prime_xyz_;
};

} //screener
} //stepwise
} //protocols

#endif
