// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/sampling/packer/StepWiseMasterPacker.hh
/// @brief
/// @detailed
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_sampling_StepWiseMasterPacker_HH
#define INCLUDED_protocols_stepwise_sampling_StepWiseMasterPacker_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/sampling/packer/StepWiseMasterPacker.fwd.hh>
#include <protocols/stepwise/sampling/packer/StepWisePacker.fwd.hh>
#include <protocols/stepwise/sampling/modeler_options/StepWiseModelerOptions.fwd.hh>
#include <protocols/stepwise/sampling/rna/phosphate/MultiPhosphateSampler.fwd.hh>
#include <protocols/stepwise/sampling/rna/o2prime/O2PrimePacker.fwd.hh>
#include <protocols/stepwise/sampling/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/screener/StepWiseScreener.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/types.hh>

#ifdef PYROSETTA
	#include <core/scoring/ScoreFunction.hh>
#endif


namespace protocols {
namespace stepwise {
namespace sampling {
namespace packer {

	class StepWiseMasterPacker: public utility::pointer::ReferenceCount {

	public:

		//constructor
		StepWiseMasterPacker( working_parameters::StepWiseWorkingParametersCOP working_parameters,
													StepWiseModelerOptionsCOP options );

		//destructor
		~StepWiseMasterPacker();

	public:

		void
		initialize( core::pose::Pose const & pose );

		void
		add_packer_screeners( utility::vector1< screener::StepWiseScreenerOP > & screeners,
													core::pose::Pose const & pose,
													core::pose::PoseOP sugar_instantiation_pose );

		void
		reset( core::pose::Pose const & pose );

		void
		do_prepack( core::pose::Pose & pose );

		void	set_working_pack_res( utility::vector1< core::Size > const & setting );
		utility::vector1< core::Size > const & working_pack_res() const { return working_pack_res_; }

		bool working_pack_res_was_inputted() const;

		void set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn ){ scorefxn_ = scorefxn; }

		core::scoring::ScoreFunctionCOP scorefxn() const { return scorefxn_; }

		packer::StepWisePackerCOP packer();

	private:

		void
		initialize_packer();

	private:

		working_parameters::StepWiseWorkingParametersCOP working_parameters_;
		StepWiseModelerOptionsCOP options_;
		core::scoring::ScoreFunctionCOP scorefxn_;
		core::scoring::ScoreFunctionCOP phosphate_scorefxn_;
		utility::vector1< core::Size > working_pack_res_;

		core::pose::PoseOP packer_pose_;
    packer::StepWisePackerOP packer_;
		rna::o2prime::O2PrimePackerOP o2prime_packer_; // deprecate after refactoring of packer
		rna::phosphate::MultiPhosphateSamplerOP phosphate_sampler_; // deprecate after refactoring of packer

	};

} //packer
} //sampling
} //stepwise
} //protocols

#endif
