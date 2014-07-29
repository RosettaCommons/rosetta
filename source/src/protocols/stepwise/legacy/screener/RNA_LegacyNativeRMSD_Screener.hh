// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/RNA_LegacyNativeRMSD_Screener.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_RNA_LegacyNativeRMSD_Screener_HH
#define INCLUDED_protocols_stepwise_screener_RNA_LegacyNativeRMSD_Screener_HH

#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/screener/RNA_LegacyNativeRMSD_Screener.fwd.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>

using namespace core;
using namespace protocols::stepwise::modeler::working_parameters;

namespace protocols {
namespace stepwise {
namespace legacy {
namespace screener {

	class RNA_LegacyNativeRMSD_Screener: public stepwise::screener::StepWiseScreener {

	public:

		//constructor
		RNA_LegacyNativeRMSD_Screener( pose::Pose const & native_pose,
												 pose::Pose & screening_pose,
												 StepWiseWorkingParametersCOP working_parameters,
												 Real const native_screen_rmsd_cutoff,
												 bool const do_screen = true );

		//destructor
		~RNA_LegacyNativeRMSD_Screener();

	public:

		bool
		check_screen();

		void set_do_screen( bool const & setting ) { do_screen_ = setting; }
		bool do_screen() const { return do_screen_; }

		std::string
		name() const { return "RNA_LegacyRNA_LegacyNativeRMSD_Screener"; }

		stepwise::screener::StepWiseScreenerType
		type() const { return RNA_LEGACY_NATIVE_RMSD; }

		Size pass_count() const { return pass_count_; }

	private:

		pose::Pose const & native_pose_;
		pose::Pose & screening_pose_;
		StepWiseWorkingParametersCOP working_parameters_;
		Real const native_screen_rmsd_cutoff_;
		bool do_screen_;
		Size pass_count_;

	};

} //screener
} //legacy
} //stepwise
} //protocols

#endif
