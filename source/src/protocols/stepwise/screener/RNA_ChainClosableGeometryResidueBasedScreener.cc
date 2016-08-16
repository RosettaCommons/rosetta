// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/screener/RNA_ChainClosableGeometryResidueBasedScreener.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/RNA_ChainClosableGeometryResidueBasedScreener.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_ChainClosableGeometryChecker.hh>
#include <protocols/stepwise/sampler/rigid_body/RigidBodyStepWiseSamplerWithResidueAlternatives.hh>

#include <basic/Tracer.hh>

static THREAD_LOCAL basic::Tracer TR( "protocols.stepwise.screener.RNA_ChainClosableGeometryResidueBasedScreener" );

namespace protocols {
namespace stepwise {
namespace screener {

//Constructor
RNA_ChainClosableGeometryResidueBasedScreener::RNA_ChainClosableGeometryResidueBasedScreener( modeler::rna::checker::RNA_ChainClosableGeometryCheckerOP chain_closable_geometry_checker ):
	StepWiseResiduePairScreener( chain_closable_geometry_checker->five_prime_chain_break_res(),
	chain_closable_geometry_checker->three_prime_chain_break_res() ),
	chain_closable_geometry_checker_( chain_closable_geometry_checker )
{}

//Destructor
RNA_ChainClosableGeometryResidueBasedScreener::~RNA_ChainClosableGeometryResidueBasedScreener()
{}

///////////////////////////////////////////////////////////////////
void
RNA_ChainClosableGeometryResidueBasedScreener::get_update( sampler::StepWiseSamplerBaseOP sampler ){

	using namespace sampler;
	using namespace sampler::rigid_body;

	runtime_assert( sampler->type() == RIGID_BODY_WITH_RESIDUE_ALTERNATIVES );
	RigidBodyStepWiseSamplerWithResidueAlternatives & rigid_body_rotamer_with_residue_alternatives = *( static_cast< RigidBodyStepWiseSamplerWithResidueAlternatives * >( sampler.get() ) );

	five_prime_xyz_ = rigid_body_rotamer_with_residue_alternatives.get_xyz( res1_, " O3'" ); // 5' residue
	three_prime_xyz_ = rigid_body_rotamer_with_residue_alternatives.get_xyz( res2_, " C5'" ); // 3' residue
	//  TR << "RES: " << res1_ << " " << res2_ <<  "  FIVE_PRIME X: " << five_prime_xyz_.x() << " -- " << "THREE_PRIME X: " << three_prime_xyz_.x() << std::endl;
	return;
}


///////////////////////////////////////////////////////////////////
bool
RNA_ChainClosableGeometryResidueBasedScreener::check_screen() {
	bool const ok = chain_closable_geometry_checker_->check_chain_closable_geometry( five_prime_xyz_, three_prime_xyz_ );
	//  TR << chain_closable_geometry_checker_->dist_squared() << " < " << chain_closable_geometry_checker_->max_dist_squared() << " " << ok << std::endl;
	return ok;
}

} //screener
} //stepwise
} //protocols
