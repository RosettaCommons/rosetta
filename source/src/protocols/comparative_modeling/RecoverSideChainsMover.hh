// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @brief simple mover that applies another Mover, and then recovers the input
/// sidechains from the input Pose.
/// @author James Thompson

#ifndef INCLUDED_protocols_comparative_modeling_RecoverSideChainsMover_hh
#define INCLUDED_protocols_comparative_modeling_RecoverSideChainsMover_hh

#include <core/pose/Pose.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace comparative_modeling {

class RecoverSideChainsMover : public protocols::moves::Mover {

public:

	RecoverSideChainsMover(
		Mover & mover
	);

	virtual void apply( core::pose::Pose & pose );
	virtual std::string get_name() const;

private:
	protocols::moves::Mover & mover_;
};

} // comparative_modeling
} // protocols

#endif
