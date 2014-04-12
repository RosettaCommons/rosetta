// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/ProteinPackScreener.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_ProteinPackScreener_HH
#define INCLUDED_protocols_stepwise_screener_ProteinPackScreener_HH

#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/screener/ProteinPackScreener.fwd.hh>
#include <protocols/stepwise/sampling/protein/StepWiseProteinPacker.fwd.hh>
#include <core/pose/Pose.fwd.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace screener {

	class ProteinPackScreener:  public StepWiseScreener {

	public:

		//constructor
		ProteinPackScreener( pose::Pose & pose,
												 sampling::protein::StepWiseProteinPackerOP stepwise_packer );

		//destructor
		~ProteinPackScreener();

	public:

		bool
		check_screen();

		std::string
		name() const { return "ProteinPackScreener"; }

		StepWiseScreenerType
		type() const { return PROTEIN_PACK; }

	private:

		pose::Pose & pose_;
		sampling::protein::StepWiseProteinPackerOP stepwise_packer_;

	};

} //screener
} //stepwise
} //protocols

#endif
