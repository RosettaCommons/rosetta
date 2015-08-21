// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/BulgeApplier.cc
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#include <protocols/stepwise/screener/BulgeApplier.hh>
#include <protocols/stepwise/modeler/rna/bulge/BulgeApplyMover.hh>
#include <protocols/stepwise/modeler/rna/bulge/BulgeUnApplyMover.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_AtrRepChecker.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_BaseCentroidChecker.hh>
#include <protocols/moves/CompositionMover.hh>
#include <basic/Tracer.hh>

static thread_local basic::Tracer TR( "protocols.stepwise.screener.BulgeApplier" );

using namespace protocols::stepwise::modeler::rna::checker;

namespace protocols {
namespace stepwise {
namespace screener {

//Constructor
BulgeApplier::BulgeApplier( RNA_AtrRepCheckerOP atr_rep_checker,
	RNA_BaseCentroidCheckerOP base_centroid_checker,
	Size const moving_res ):
	atr_rep_checker_( atr_rep_checker ),
	base_centroid_checker_( base_centroid_checker ),
	bulge_apply_mover_( modeler::rna::bulge::BulgeApplyMoverOP( new modeler::rna::bulge::BulgeApplyMover( moving_res ) ) ),
	bulge_unapply_mover_( modeler::rna::bulge::BulgeUnApplyMoverOP( new modeler::rna::bulge::BulgeUnApplyMover( moving_res ) ) )
{}

//Destructor
BulgeApplier::~BulgeApplier()
{}

////////////////////////////////////////////////////////////////////////////////////////
void
BulgeApplier::add_mover( moves::CompositionMoverOP update_mover, moves::CompositionMoverOP restore_mover ){
	using namespace protocols::stepwise::modeler::rna::bulge;
	if ( bulge_variant_decision() ) {
		update_mover->add_mover( bulge_apply_mover_ );
		restore_mover->add_mover( bulge_unapply_mover_ );
	} else {
		update_mover->add_mover( 0 );
		restore_mover->add_mover( 0 );
	}
}

////////////////////////////////////////////////////////////////////////////////////////
bool
BulgeApplier::bulge_variant_decision(){

	if ( base_centroid_checker_->found_centroid_interaction() ) return false;

	static Real const atr_cutoff_for_bulge( -999999.0 ); //Feb 02, 2012
	Real const delta_atr_score = atr_rep_checker_->delta_atr_score();
	runtime_assert ( delta_atr_score <= (  + 0.01 ) );

	return ( delta_atr_score >= atr_cutoff_for_bulge );
}


} //screener
} //stepwise
} //protocols
