// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loophash/LoopHashSampler.hh
/// @brief
/// @author Mike Tyka


#ifndef INCLUDED_protocols_loophash_LoopHashRelaxProtocol_hh
#define INCLUDED_protocols_loophash_LoopHashRelaxProtocol_hh

#include <protocols/loophash/LoopHashSampler.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LocalInserter.fwd.hh>
#include <protocols/moves/Mover.hh>

#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/io/silent/SilentStruct.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>
#include <string>
#include <vector>

#include <utility/vector1.hh>


namespace protocols {
namespace loophash {

class LoopHashRelaxProtocol;
typedef utility::pointer::shared_ptr< LoopHashRelaxProtocol > LoopHashRelaxProtocolOP;
typedef utility::pointer::shared_ptr< LoopHashRelaxProtocol const > LoopHashRelaxProtocolCOP;


class LoopHashRelaxProtocol: public protocols::moves::Mover {
public:

	LoopHashRelaxProtocol(
		LoopHashLibraryOP library
	);

	virtual void apply( core::pose::Pose& pose );

	void manual_call( core::pose::Pose& pose );

	virtual protocols::moves::MoverOP clone() const;

	virtual std::string get_name() const {
		return "LoopHashRelaxProtocol";
	}

	virtual protocols::moves::MoverOP fresh_instance() const;

private:
	LoopHashLibraryOP library_;

};

}
}

#endif
