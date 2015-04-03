// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/BulgeApplier.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_BulgeApplier_HH
#define INCLUDED_protocols_stepwise_screener_BulgeApplier_HH

#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_AtrRepChecker.fwd.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_BaseCentroidChecker.fwd.hh>
#include <protocols/stepwise/screener/BulgeApplier.fwd.hh>
#include <protocols/stepwise/modeler/rna/bulge/BulgeApplyMover.fwd.hh>
#include <protocols/stepwise/modeler/rna/bulge/BulgeUnApplyMover.fwd.hh>

using namespace protocols::stepwise::modeler::rna::checker;

namespace protocols {
namespace stepwise {
namespace screener {

	class BulgeApplier: public StepWiseScreener {

	public:

		//constructor
		BulgeApplier( RNA_AtrRepCheckerOP atr_rep_checker, RNA_BaseCentroidCheckerOP base_centroid_checker,
									Size const moving_res );

		//destructor
		~BulgeApplier();

	public:

		std::string
		name() const { return "BulgeApplier"; }

		StepWiseScreenerType
		type() const { return BULGE_APPLIER; }

		virtual
		void
		add_mover( moves::CompositionMoverOP update_mover, moves::CompositionMoverOP restore_mover );

	private:

		bool bulge_variant_decision();

	private:

		RNA_AtrRepCheckerOP atr_rep_checker_;
		RNA_BaseCentroidCheckerOP base_centroid_checker_;
		modeler::rna::bulge::BulgeApplyMoverOP   bulge_apply_mover_;
		modeler::rna::bulge::BulgeUnApplyMoverOP bulge_unapply_mover_;

	};

} //screener
} //stepwise
} //protocols

#endif
