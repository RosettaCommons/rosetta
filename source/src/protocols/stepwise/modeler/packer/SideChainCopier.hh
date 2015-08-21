// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/stepwise/modeler/packer/SideChainCopier.hh
/// @brief
/// @details
/// @author Rhiju Das, rhiju@stanford.edu


#ifndef INCLUDED_protocols_stepwise_modeler_packer_SideChainCopier_HH
#define INCLUDED_protocols_stepwise_modeler_packer_SideChainCopier_HH

#include <protocols/moves/Mover.hh>
#include <protocols/stepwise/modeler/packer/SideChainCopier.fwd.hh>
#include <protocols/stepwise/modeler/util.hh>

namespace protocols {
namespace stepwise {
namespace modeler {
namespace packer {

class SideChainCopier: public protocols::moves::Mover {

public:

	//constructor
	SideChainCopier( core::pose::Pose const & reference_pose,
		bool const copy_o2prime_hydrogens = false  );

	//constructor
	SideChainCopier( core::pose::Pose const & reference_pose,
		utility::vector1< Size > const & copy_res,
		bool const copy_o2prime_hydrogens = false  );

	//destructor
	~SideChainCopier();

public:

	void
	apply( core::pose::Pose & pose );

	std::string get_name() const { return "SideChainCopier"; }

private:

	core::pose::Pose const & reference_pose_;
	utility::vector1< Size > copy_res_;
	bool copy_o2prime_hydrogens_;

};

} //packer
} //modeler
} //stepwise
} //protocols

#endif
