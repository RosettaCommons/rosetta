// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/stepwise/modeler/packer/StepWiseMasterPacker.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_StepWiseMasterPacker_HH
#define INCLUDED_protocols_stepwise_modeler_StepWiseMasterPacker_HH

#include <utility/pointer/ReferenceCount.hh>
#include <protocols/stepwise/modeler/packer/StepWiseMasterPacker.fwd.hh>
#include <protocols/stepwise/modeler/packer/StepWisePacker.fwd.hh>
#include <protocols/stepwise/modeler/options/StepWiseModelerOptions.fwd.hh>
#include <protocols/stepwise/modeler/rna/phosphate/MultiPhosphateSampler.fwd.hh>
#include <protocols/stepwise/modeler/rna/o2prime/O2PrimePacker.fwd.hh>
#include <protocols/stepwise/modeler/working_parameters/StepWiseWorkingParameters.fwd.hh>
#include <protocols/stepwise/screener/StepWiseScreener.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/types.hh>


namespace protocols {
namespace stepwise {
namespace modeler {
namespace packer {

class StepWiseMasterPacker: public utility::pointer::ReferenceCount {

public:

	//constructor
	StepWiseMasterPacker( working_parameters::StepWiseWorkingParametersCOP working_parameters,
		options::StepWiseModelerOptionsCOP options );

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

	void set_working_pack_res( utility::vector1< core::Size > const & setting );
	utility::vector1< core::Size > const & working_pack_res() const { return working_pack_res_; }

	// Undefined, commenting out to fix PyRosetta build  bool working_pack_res_was_inputted() const;

	void set_scorefxn( core::scoring::ScoreFunctionCOP scorefxn ){ scorefxn_ = scorefxn; }

	core::scoring::ScoreFunctionCOP scorefxn() const { return scorefxn_; }

	packer::StepWisePackerCOP packer();

private:

	void
	initialize_packer();

private:

	working_parameters::StepWiseWorkingParametersCOP working_parameters_;
	options::StepWiseModelerOptionsCOP options_;
	core::scoring::ScoreFunctionCOP scorefxn_;
	core::scoring::ScoreFunctionCOP phosphate_scorefxn_;
	utility::vector1< core::Size > working_pack_res_;

	core::pose::PoseOP packer_pose_;
	packer::StepWisePackerOP packer_;
	rna::o2prime::O2PrimePackerOP o2prime_packer_; // deprecate after refactoring of packer
	rna::phosphate::MultiPhosphateSamplerOP phosphate_sampler_; // deprecate after refactoring of packer

};

} //packer
} //modeler
} //stepwise
} //protocols

#endif
