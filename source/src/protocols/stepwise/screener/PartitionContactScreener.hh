// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/screener/PartitionContactScreener.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_screener_PartitionContactScreener_HH
#define INCLUDED_protocols_stepwise_screener_PartitionContactScreener_HH

#include <protocols/stepwise/screener/StepWiseScreener.hh>
#include <protocols/stepwise/screener/PartitionContactScreener.fwd.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/etable/EtableEnergy.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/types.hh>

namespace protocols {
namespace stepwise {
namespace screener {

class PartitionContactScreener: public StepWiseScreener {

public:

	//constructor
	PartitionContactScreener( core::pose::Pose const & pose,
		modeler::working_parameters::StepWiseWorkingParametersCOP working_parameters,
		bool const use_loose_rep_cutoff,
		core::scoring::methods::EnergyMethodOptions const & options /* how to setup etable */ );

	//destructor
	~PartitionContactScreener();

public:

	std::string
	name() const { return "PartitionContactScreener"; }

	StepWiseScreenerType
	type() const { return PARTITION_CONTACT; }

	bool
	check_screen();

private:

	void
	initialize_actual_rep_cutoff();

	void
	initialize_evaluator( core::scoring::methods::EnergyMethodOptions const & options );

	void
	check_screen( Size const moving_res, bool & atr_ok, bool & rep_ok ) const;

private:

	core::pose::Pose const & pose_;
	utility::vector1< core::Size > const & moving_res_list_;
	core::Real const fa_atr_weight_, fa_rep_weight_, rep_cutoff_, atr_cutoff_;
	bool const use_loose_rep_cutoff_, close_chain_;
	core::Real actual_rep_cutoff_;
	core::scoring::etable::AnalyticEtableEvaluatorOP eval_;
};

} //screener
} //stepwise
} //protocols

#endif
