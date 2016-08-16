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
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/loops/loops_main.hh>
#include <core/pose/Pose.hh>
#include <core/fragment/FragSet.hh>
#include <core/fragment/FragSet.fwd.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>
#include <core/scoring/constraints/BoundConstraint.hh>

#include <protocols/relax/FastRelax.hh>
#include <protocols/viewer/viewers.hh>

#include <protocols/moves/Mover.hh>
#include <protocols/moves/CompositionMover.hh>

// C++ headers
#include <iostream>
#include <string>

#include <core/types.hh>
#include <utility/vector1.hh>

// option key includes
#include <basic/options/option.hh>
#include <basic/options/keys/cm.OptionKeys.gen.hh>
#include <basic/options/keys/loops.OptionKeys.gen.hh>

#include <core/kinematics/Jump.hh>
#include <numeric/random/random.hh>
#include <core/scoring/ScoreFunction.hh>

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
	);

	std::string
	get_name() const {
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
