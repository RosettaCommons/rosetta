// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/legacy/modeler/StepWiseLegacyMinimizer.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_StepWiseLegacyMinimizer_HH
#define INCLUDED_protocols_stepwise_modeler_StepWiseLegacyMinimizer_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/legacy/modeler/StepWiseLegacyMinimizer.fwd.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.fwd.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_BaseCentroidChecker.fwd.hh>
#include <protocols/stepwise/modeler/rna/checker/RNA_VDW_BinChecker.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.hh>
#include <core/types.hh>

namespace protocols {
namespace stepwise {
namespace legacy {
namespace modeler {

	class StepWiseLegacyMinimizer: public moves::Mover {

	public:

		//constructor
		StepWiseLegacyMinimizer( utility::vector1< core::pose::PoseOP > & pose_list,
											 stepwise::modeler::working_parameters::StepWiseWorkingParametersCOP working_parameters,
											 options::StepWiseModelerOptionsCOP options,
											 core::scoring::ScoreFunctionCOP scorefxn);

		//destructor
		~StepWiseLegacyMinimizer();

	public:

		virtual void apply( core::pose::Pose & pose );

		virtual std::string get_name() const{ return "StepWiseLegacyMinimizer"; }

		void
		do_unified_minimizing( core::pose::Pose & pose );

		void
		do_protein_minimizing( core::pose::Pose & pose );

		void
		do_rna_minimizing( core::pose::Pose & pose );

    void set_working_obligate_pack_res( utility::vector1< core::Size > const & setting ){ working_obligate_pack_res_ = setting; }
		void
		set_base_centroid_checker( rna::checker::RNA_BaseCentroidCheckerOP checker ) { base_centroid_checker_ = checker; }

		void
		set_user_input_VDW_bin_checker( rna::checker::RNA_VDW_BinCheckerOP user_input_VDW_bin_checker ){ user_input_VDW_bin_checker_ = user_input_VDW_bin_checker; }

		void set_protein_full_optimize( bool const & setting ){ protein_full_optimize_ = setting; }
		bool protein_full_optimize() const{ return protein_full_optimize_; }

	private:

		utility::vector1< core::pose::PoseOP > & pose_list_; // where work will be saved.
		stepwise::modeler::working_parameters::StepWiseWorkingParametersCOP working_parameters_;
		options::StepWiseModelerOptionsCOP options_;
		core::scoring::ScoreFunctionCOP scorefxn_;

		utility::vector1< core::Size > moving_res_list_;
    utility::vector1< core::Size > working_obligate_pack_res_;
		bool protein_full_optimize_;

		rna::checker::RNA_BaseCentroidCheckerOP base_centroid_checker_;
		rna::checker::RNA_VDW_BinCheckerOP user_input_VDW_bin_checker_;

	};

} //modeler
} //legacy
} //stepwise
} //protocols

#endif
