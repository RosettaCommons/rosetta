// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/legacy/screener/RNA_AtrRepScreener.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_RNA_AtrRepScreener_HH
#define INCLUDED_protocols_stepwise_screener_RNA_AtrRepScreener_HH

#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/legacy/screener/RNA_AtrRepScreener.fwd.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_AtrRepChecker.fwd.hh>
#include <core/pose/Pose.fwd.hh>


#ifdef WIN32
	#include <protocols/stepwise/modeler/rna/checker/RNA_AtrRepChecker.hh>
	#include <protocols/stepwise/legacy/screener/RNA_AtrRepScreener.hh>
#endif


namespace protocols {
namespace stepwise {
namespace legacy {
namespace screener {

	class RNA_AtrRepScreener: public stepwise::screener::StepWiseScreener {

	public:

		//constructor
		RNA_AtrRepScreener( protocols::stepwise::modeler::rna::checker::RNA_AtrRepCheckerOP atr_rep_checker,
										core::pose::Pose & screening_pose );

		//destructor
		~RNA_AtrRepScreener();

	public:

		std::string
		name() const { return "RNA_AtrRepScreener"; }

		stepwise::screener::StepWiseScreenerType
		type() const { return stepwise::screener::ATR_REP; }

		bool
		check_screen();

		void
		set_exit_on_fail( bool const setting ){ exit_on_fail_ = setting; }

	private:

		protocols::stepwise::modeler::rna::checker::RNA_AtrRepCheckerOP atr_rep_checker_;
		core::pose::Pose & screening_pose_;
		bool exit_on_fail_;

	};

} //screener
} //legacy
} //stepwise
} //protocols

#endif
