// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software res and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file SWA_ResidueSampler.hh
/// @brief
/// @detailed
///
/// @author Parin Sripakdeevong
/// @author Rhiju Das


#ifndef INCLUDED_protocols_stepwise_rna_StepWiseRNA_ResidueSampler_HH
#define INCLUDED_protocols_stepwise_rna_StepWiseRNA_ResidueSampler_HH

#include <protocols/stepwise/sampling/rna/StepWiseRNA_JobParameters.fwd.hh>
#include <protocols/stepwise/sampling/rna/StepWiseRNA_ModelerOptions.fwd.hh>
#include <protocols/stepwise/sampling/rna/sugar/StepWiseRNA_VirtualSugarJustInTimeInstantiator.fwd.hh>
#include <protocols/stepwise/sampling/rna/sugar/SugarModeling.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_BaseCentroidChecker.fwd.hh>
#include <protocols/stepwise/sampling/rna/checker/RNA_VDW_BinChecker.fwd.hh>
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/io/silent/SilentFileData.fwd.hh>
#include <utility/vector1.hh>
#include <protocols/moves/Mover.hh>
#include <string>
#include <map>

using namespace core;

namespace protocols {
namespace stepwise {
namespace sampling {
namespace rna {

class StepWiseRNA_ResidueSampler: public protocols::moves::Mover {

public:

	//constructor!
	StepWiseRNA_ResidueSampler( StepWiseRNA_JobParametersCOP & job_parameters_ );

	//destructor -- necessary?
	~StepWiseRNA_ResidueSampler();

	virtual void apply( core::pose::Pose & pose_to_visualize );

	virtual std::string get_name() const;

	void
	set_scorefxn( core::scoring::ScoreFunctionOP const & scorefxn );

	utility::vector1< core::pose::PoseOP > &
	get_pose_list();

	void output_pose_list( std::string const final_sampler_output_silent_file ) const;

	void
	set_base_centroid_checker( checker::RNA_BaseCentroidCheckerOP & checker );

	void
	set_user_input_VDW_bin_checker( checker::RNA_VDW_BinCheckerOP const & user_input_VDW_bin_checker );

	void
	set_options( StepWiseRNA_ModelerOptionsCOP options );

	void set_sampling_silent_file( std::string const & setting ){ sampling_silent_file_ = setting; }
	std::string sampling_silent_file() const{ return sampling_silent_file_; }

private:

	void
	initialize_scorefunctions();

	void
	check_res_not_bulged();

	void
	unified_sampling( core::pose::Pose & pose );

	// will be deprecated soon:
	void	legacy_sampling( core::pose::Pose & pose );
	void	standard_sampling( core::pose::Pose & pose );
	void	floating_base_sampling( core::pose::Pose & pose );

	bool
	instantiate_any_virtual_sugars( core::pose::Pose & pose );

	void
	output_options();

private:

	StepWiseRNA_JobParametersCOP job_parameters_;
	utility::vector1< pose::PoseOP > pose_list_;
	core::scoring::ScoreFunctionOP scorefxn_;
	StepWiseRNA_ModelerOptionsCOP options_;
	checker::RNA_BaseCentroidCheckerOP base_centroid_checker_;
	checker::RNA_VDW_BinCheckerOP user_input_VDW_bin_checker_;
	sugar::StepWiseRNA_VirtualSugarJustInTimeInstantiatorOP virtual_sugar_just_in_time_instantiator_;
	std::string sampling_silent_file_;

};

} //rna
} //sampling
} //stepwise
} //protocols

#endif
