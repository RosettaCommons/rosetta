// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file AbrelaxMover
/// @brief  this class will be handled to a SampleProtocol as a control instance
/// @details responsibilities:
///           know which chainbreaks to penalize and close
///           know which jumps to use during sampling, which (if any) to keep after loop-closing
///           supply a JumpMover if jumps should be moved
///           supply a MoveMap
///           supply a "StrictMoveMap": the protocol should not move anything that is dissallowed in strict_movemap(),
///                      it should try to move just stuff in movemap()
/// should this class also know how to ramp score terms ?
/// handle the titration of constraints ?
/// @author Oliver Lange


#ifndef INCLUDED_protocols_abinitio_IterativeFullatom_hh
#define INCLUDED_protocols_abinitio_IterativeFullatom_hh

// Unit Headers
//#include <protocols/abinitio/IterativeFullatom.fwd.hh>

// Package Headers
#include <protocols/abinitio/IterativeBase.hh>

// Project Headers
#include <core/types.hh>
#include <core/pose/Pose.fwd.hh>

//// C++ headers
#include <string>

#include <utility/vector1.hh>


namespace protocols {
namespace abinitio {

class IterativeFullatom : public IterativeBase {
	typedef IterativeBase Parent;
public:
	static void register_options();

	IterativeFullatom();

	/// @brief do initializing work that requires fully setup object here
	void initialize() override;

	bool ready_for_batch() const override;

	void generate_batch() override {
		Parent::generate_batch();
	}

	core::Size generate_batch( jd2::archive::Batch&, core::Size repeat_id ) override;

protected:
	void gen_resample_core( jd2::archive::Batch& batch, bool flex );

private:
	static bool options_registered_;
	core::Real perturb_start_structures_;

};


}
}

#endif
