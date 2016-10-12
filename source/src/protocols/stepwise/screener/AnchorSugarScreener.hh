// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/screener/AnchorSugarScreener.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_AnchorSugarScreener_HH
#define INCLUDED_protocols_stepwise_screener_AnchorSugarScreener_HH

#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/screener/AnchorSugarScreener.fwd.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_ChainClosableGeometryChecker.fwd.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_AtrRepChecker.fwd.hh>
#include <protocols/stepwise/screener/TagDefinition.fwd.hh>
#include <protocols/stepwise/modeler/rna/sugar/SugarModeling.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

#ifdef WIN32
#include <protocols/stepwise/modeler/rna/checker/RNA_AtrRepChecker.hh>
#include <protocols/stepwise/screener/AnchorSugarScreener.hh>
#include <protocols/stepwise/modeler/rna/sugar/SugarModeling.hh>
#endif

// functionality in original floating_base modeler, setup for sample-and-screen framework.
// This will likely be deprecated in favor of more general ResidueAlternative modeler

namespace protocols {
namespace stepwise {
namespace screener {

class AnchorSugarScreener: public StepWiseScreener {

public:

	//constructor
	AnchorSugarScreener(modeler::rna::sugar::SugarModeling const & anchor_sugar_modeling,
		modeler::rna::checker::RNA_ChainClosableGeometryCheckerOP chain_closable_geometry_to_anchor_checker,
		core::pose::Pose & sugar_screening_pose,
		bool const is_prepend,
		modeler::rna::checker::RNA_AtrRepCheckerOP atr_rep_checker_with_instantiated_sugar,
		utility::vector1< modeler::rna::checker::RNA_AtrRepCheckerOP > atr_rep_checkers_for_anchor_sugar_models,
		TagDefinitionOP tag_definition );

	//destructor
	~AnchorSugarScreener();

public:

	bool
	check_screen();

	std::string
	name() const { return "AnchorSugarScreener"; }

	StepWiseScreenerType
	type() const { return ANCHOR_SUGAR; }

	void
	add_mover( moves::CompositionMoverOP update_mover, moves::CompositionMoverOP restore_mover );

	Size const &
	anchor_sugar_solution_number() const { return anchor_sugar_solution_number_; }

private:

	modeler::rna::sugar::SugarModeling const & anchor_sugar_modeling_;
	modeler::rna::checker::RNA_ChainClosableGeometryCheckerOP chain_closable_geometry_to_anchor_checker_;
	core::pose::Pose & sugar_screening_pose_;
	modeler::rna::checker::RNA_AtrRepCheckerOP atr_rep_checker_with_instantiated_sugar_;
	utility::vector1< modeler::rna::checker::RNA_AtrRepCheckerOP > atr_rep_checkers_for_anchor_sugar_models_;
	TagDefinitionOP tag_definition_;

	bool const is_prepend_;
	std::string const moving_atom_name_;
	std::string const reference_atom_name_;

	Size anchor_sugar_solution_number_;
};

} //screener
} //stepwise
} //protocols

#endif
