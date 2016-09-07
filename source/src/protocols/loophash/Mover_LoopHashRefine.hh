// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/loophash/Mover_LoopHashRefine.hh
/// @brief
/// @author Mike Tyka


#ifndef INCLUDED_protocols_loophash_Mover_LoopHashRefine_hh
#define INCLUDED_protocols_loophash_Mover_LoopHashRefine_hh


// libRosetta headers

#include <protocols/moves/Mover.hh>
#include <protocols/loophash/LoopHashLibrary.fwd.hh>
#include <protocols/loophash/LoopHashLibrary.hh>
#include <protocols/loophash/LoopHashMap.hh>
#include <protocols/loophash/BackboneDB.hh>
#include <protocols/checkpoint/CheckPointer.hh>

#include <core/kinematics/Jump.hh>
#include <utility>
#include <utility/vector1.hh>

//Auto Headers

namespace protocols {
namespace loophash {

class Mover_LoopHashRefine;
typedef utility::pointer::shared_ptr< Mover_LoopHashRefine > Mover_LoopHashRefineOP;
typedef utility::pointer::shared_ptr< Mover_LoopHashRefine const > Mover_LoopHashRefineCOP;

class Mover_LoopHashRefine: public protocols::moves::Mover {
public:

	Mover_LoopHashRefine(
		protocols::loophash::LoopHashLibraryOP library
	):
		library_(std::move(library)),
		checkpoints_("LoopHash")
	{
	}

	void apply( core::pose::Pose& pose ) override;

	protocols::moves::MoverOP clone() const override {
		return protocols::moves::MoverOP( new Mover_LoopHashRefine( *this ) );
	}


	std::string get_name() const override {
		return "Mover_LoopHashRefine";
	}

	protocols::moves::MoverOP fresh_instance() const override {
		return protocols::moves::MoverOP( new Mover_LoopHashRefine( library_ ) );
	}

private:
	protocols::loophash::LoopHashLibraryOP library_;
	protocols::checkpoint::CheckPointer checkpoints_;

};


int loophash_main();


}
}

#endif
