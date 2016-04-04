// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/rna/checker/RNA_AtrRepChecker.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_rna_checker_RNA_AtrRepChecker_HH
#define INCLUDED_protocols_stepwise_rna_checker_RNA_AtrRepChecker_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_AtrRepChecker.fwd.hh>
#include <protocols/stepwise/modeler/rna/StepWiseRNA_Classes.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/methods/EnergyMethodOptions.fwd.hh>
#include <core/pose/Pose.hh>

using namespace core;

namespace protocols {
namespace stepwise {
namespace modeler {
namespace rna {
namespace checker {

class RNA_AtrRepChecker: public utility::pointer::ReferenceCount {

public:

	//Constructor
	RNA_AtrRepChecker( pose::Pose const & pose,
		working_parameters::StepWiseWorkingParametersCOP & working_parameters,
		bool const loose_rep_cutoff = false,
					   core::scoring::methods::EnergyMethodOptionsCOP energy_method_options = 0 );

	RNA_AtrRepChecker( pose::Pose const & pose,
		Size const moving_res,
		Size const reference_res,
		Size const gap_size,
		core::scoring::methods::EnergyMethodOptionsCOP energy_method_options = 0
	);

	//destructor
	~RNA_AtrRepChecker();

	Real delta_atr_score() const{ return delta_atr_score_; }
	Real delta_rep_score() const{ return delta_rep_score_; }
	Real base_atr_score() const{ return base_atr_score_; }
	Real base_rep_score() const{ return base_rep_score_; }

public:

	bool
	check_screen( pose::Pose & current_pose_screen );

	void
	set_loose_rep_cutoff( bool const & setting ) { loose_rep_cutoff_ = setting; }

	void
	set_extra_loose_rep_cutoff( bool const & setting ) { extra_loose_rep_cutoff_ = setting; }

	StepWiseRNA_CountStruct const &
	count_data() const { return count_data_; }

private:

	void
	get_base_atr_rep_score( core::pose::Pose const & pose );

	void
	initialize_scorefxn( core::scoring::methods::EnergyMethodOptionsCOP energy_method_options = 0 );

	void
	initialize_parameters();

	void
	output_rep( core::pose::Pose const & pose, std::string const tag );

private:

	Size const moving_res_;
	Size const reference_res_;
	Size const gap_size_;
	bool const is_prepend_;
	bool const is_internal_;
	bool const sample_both_sugar_base_rotamer_;
	bool const separate_moving_residue_to_estimate_baseline_;

	Real rep_cutoff_;
	Real base_atr_score_;
	Real base_rep_score_;
	Real delta_atr_score_;
	Real delta_rep_score_;
	bool output_pdb_;
	bool verbose_;
	bool loose_rep_cutoff_;
	bool extra_loose_rep_cutoff_;

	core::scoring::ScoreFunctionOP atr_rep_screening_scorefxn_;

	StepWiseRNA_CountStruct count_data_;

};

} //checker
} //rna
} //modeler
} //stepwise
} //protocols

#endif
