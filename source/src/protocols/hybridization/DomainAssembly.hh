// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief

#ifndef INCLUDED_protocols_hybridization_DomainAssembly_hh
#define INCLUDED_protocols_hybridization_DomainAssembly_hh

// libRosetta headers
#include <core/pose/Pose.fwd.hh>


#include <protocols/moves/Mover.hh>

// C++ headers
#include <string>

#include <core/types.hh>
#include <utility/vector1.hh>

#include <core/scoring/ScoreFunction.fwd.hh>

namespace protocols {
//namespace comparative_modeling {
namespace hybridization {

class DomainAssembly : public protocols::moves::Mover {

public:
	DomainAssembly(
		utility::vector1 < core::pose::PoseOP > & poses,
		utility::vector1 < core::Real > & domain_assembly_weights
	);

	void run();

	DomainAssembly(core::pose::PoseOP pose1,
		core::pose::PoseOP pose2,
		core::scoring::ScoreFunctionOP scorefxn)
	{
		pose1_ = pose1;
		pose2_ = pose2;
		scorefxn_ = scorefxn;
	}

	void apply(
		core::pose::Pose & pose
	) override;

	std::string
	get_name() const override {
		return "DomainAssembly";
	}

private:
	utility::vector1 < core::pose::PoseOP > poses_;

	core::pose::PoseOP pose1_;
	core::pose::PoseOP pose2_;
	core::scoring::ScoreFunctionOP scorefxn_;
}; // class DomainAssembly

}  //  namespace hybridization
//}  //  //namespace comparative_modeling
}  //  namespace protocols

#endif  // PROTOCOLS_COMPARATIVE_MODELING_HYBRIDIZE_HYBRIDIZEFOLDTREEDYNAMIC_HH_
