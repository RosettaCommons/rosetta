// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/legacy/ProteinAtrRepScreener.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_ProteinAtrRepScreener_HH
#define INCLUDED_protocols_stepwise_screener_ProteinAtrRepScreener_HH

#include <protocols/stepwise/screener/SampleApplier.hh>
#include <protocols/stepwise/screener/legacy/ProteinAtrRepScreener.fwd.hh>
#include <protocols/stepwise/sampling/protein/checker/ProteinAtrRepChecker.fwd.hh>
#include <core/pose/Pose.fwd.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace screener {

	class ProteinAtrRepScreener: public SampleApplier {

	public:
		//constructor
		ProteinAtrRepScreener( pose::Pose & pose_atr_rep_screen,
													 sampling::protein::checker::ProteinAtrRepCheckerOP atr_rep_checker );

		//destructor
		~ProteinAtrRepScreener();

	public:

		std::string
		name() const { return "ProteinAtrRepScreener"; }

		StepWiseScreenerType
		type() const { return PROTEIN_ATR_REP; }

		bool
		check_screen();

	private:

		sampling::protein::checker::ProteinAtrRepCheckerOP atr_rep_checker_;

	};

} //screener
} //stepwise
} //protocols

#endif
