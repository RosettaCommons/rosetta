// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/RNA_ChainClosureScreener.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_RNA_ChainClosureScreener_HH
#define INCLUDED_protocols_stepwise_screener_RNA_ChainClosureScreener_HH

#include <protocols/stepwise/screener/SampleApplier.hh>
#include <protocols/stepwise/screener/RNA_ChainClosureScreener.fwd.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_ChainClosureChecker.fwd.hh>
#include <core/pose/Pose.fwd.hh>

namespace protocols {
namespace stepwise {
namespace screener {

	class RNA_ChainClosureScreener: public SampleApplier {

	public:

		//constructor
		RNA_ChainClosureScreener( modeler::rna::checker::RNA_ChainClosureCheckerOP chain_closure_checker,
													core::pose::Pose & screening_pose,
													bool const just_do_closure_check = false );

		RNA_ChainClosureScreener( modeler::rna::checker::RNA_ChainClosureCheckerOP chain_closure_checker );

		//destructor
		~RNA_ChainClosureScreener();

	public:

		std::string
		name() const { return "RNA_ChainClosureScreener"; }

		StepWiseScreenerType
		type() const { return RNA_CHAIN_CLOSURE; }

		bool
		check_screen();

		void
		add_mover( moves::CompositionMoverOP update_mover, moves::CompositionMoverOP restore_mover );

		void
		fast_forward( sampler::StepWiseSamplerBaseOP sampler );

	private:

		modeler::rna::checker::RNA_ChainClosureCheckerOP chain_closure_checker_;
		core::pose::Pose & screening_pose_;
		bool const just_do_closure_check_;

	};

} //screener
} //stepwise
} //protocols

#endif
