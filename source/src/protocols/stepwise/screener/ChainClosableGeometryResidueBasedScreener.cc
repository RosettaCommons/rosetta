// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/ChainClosableGeometryResidueBasedScreener.cc
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/ChainClosableGeometryResidueBasedScreener.hh>
#include <protocols/stepwise/sampling/rna/checker/ChainClosableGeometryChecker.hh>
#include <protocols/rotamer_sampler/rigid_body/RigidBodyRotamerWithResidueAlternatives.hh>

#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.stepwise.screener.ChainClosableGeometryResidueBasedScreener" );

namespace protocols {
namespace stepwise {
namespace screener {

	//Constructor
	ChainClosableGeometryResidueBasedScreener::ChainClosableGeometryResidueBasedScreener( sampling::rna::checker::ChainClosableGeometryCheckerOP chain_closable_geometry_checker ):
		StepWiseResiduePairScreener( chain_closable_geometry_checker->five_prime_chain_break_res(),
																 chain_closable_geometry_checker->three_prime_chain_break_res() ),
		chain_closable_geometry_checker_( chain_closable_geometry_checker )
	{}

	//Destructor
	ChainClosableGeometryResidueBasedScreener::~ChainClosableGeometryResidueBasedScreener()
	{}

	///////////////////////////////////////////////////////////////////
	void
	ChainClosableGeometryResidueBasedScreener::get_update( rotamer_sampler::RotamerBaseOP sampler ){

		using namespace rotamer_sampler;
		using namespace rotamer_sampler::rigid_body;

		runtime_assert( sampler->type() == RIGID_BODY_WITH_RESIDUE_ALTERNATIVES );
		RigidBodyRotamerWithResidueAlternatives & rigid_body_rotamer_with_residue_alternatives = *( static_cast< RigidBodyRotamerWithResidueAlternatives * >( sampler.get() ) );

		five_prime_xyz_ = rigid_body_rotamer_with_residue_alternatives.get_xyz( res1_, " O3'" ); // 5' residue
		three_prime_xyz_ = rigid_body_rotamer_with_residue_alternatives.get_xyz( res2_, " C5'" ); // 3' residue
		return;
	}


	///////////////////////////////////////////////////////////////////
	bool
	ChainClosableGeometryResidueBasedScreener::check_screen() {
		bool const ok = chain_closable_geometry_checker_->check_chain_closable_geometry( five_prime_xyz_, three_prime_xyz_ );
		//		TR << chain_closable_geometry_checker_->dist_squared() << " " << ok << std::endl;
		return ok;
	}

} //screener
} //stepwise
} //protocols
