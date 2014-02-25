// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/ChainClosableGeometryScreener.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_ChainClosableGeometryScreener_HH
#define INCLUDED_protocols_stepwise_screener_ChainClosableGeometryScreener_HH

#include <protocols/stepwise/screener/StepWiseResiduePairScreener.hh>
#include <protocols/stepwise/screener/ChainClosableGeometryScreener.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/ChainClosableGeometryChecker.fwd.hh>
#include <core/pose/Pose.fwd.hh>


namespace protocols {
namespace stepwise {
namespace screener {

	class ChainClosableGeometryScreener: public StepWiseResiduePairScreener {

	public:

		//constructor
		ChainClosableGeometryScreener( sampling::rna::checker::ChainClosableGeometryCheckerOP chain_closable_geometry_checker,
													 core::pose::PoseOP screening_pose,
													 bool const finer_sampling_at_chain_closure = false );

		//destructor
		~ChainClosableGeometryScreener();

	public:

		bool
		check_screen();

		std::string
		name() const { return "ChainClosableGeometryScreener"; }

		StepWiseScreenerType
		type() const { return CHAIN_CLOSABLE_GEOMETRY; }

	private:

		sampling::rna::checker::ChainClosableGeometryCheckerOP chain_closable_geometry_checker_;
		core::pose::PoseOP screening_pose_;
		bool const finer_sampling_at_chain_closure_;

	};

} //screener
} //stepwise
} //protocols

#endif
