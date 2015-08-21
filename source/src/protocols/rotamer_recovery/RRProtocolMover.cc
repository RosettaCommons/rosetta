// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/rotamer_recovery/RRProtocolMover.cc
/// @brief  Preform the rotamer recovery after applying mover
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/rotamer_recovery/RRProtocolMover.hh>

// Project Headers
#include <protocols/rotamer_recovery/RRComparer.hh>
#include <protocols/rotamer_recovery/RRReporter.hh>

// Platform Headers
#include <basic/Tracer.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/min_pack.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// C++ Headers
#include <string>

//Auto Headers
#include <protocols/moves/Mover.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>


using std::string;
using core::Size;
using core::pack::min_pack;
using core::pose::Pose;
using core::scoring::ScoreFunction;
using core::scoring::get_score_function;
using core::pack::task::PackerTask;
using protocols::moves::MoverOP;
using basic::Tracer;

namespace protocols {
namespace rotamer_recovery {

static Tracer TR("protocol.rotamer_recovery.RRProtocolMover");

RRProtocolMover::RRProtocolMover() :
	mover_(/* NULL */)
{}

RRProtocolMover::RRProtocolMover(
	MoverOP mover
) :
	mover_(mover)
{}

RRProtocolMover::RRProtocolMover( RRProtocolMover const & src) :
	RRProtocol(),
	mover_(src.mover_)
{}

RRProtocolMover::~RRProtocolMover() {}

string
RRProtocolMover::get_name() const {
	return "RRProtocolMover";
}

string
RRProtocolMover::get_parameters() const {
	return "";
}


/// @details apply Mover and measure rotamer recovery for each residue
void
RRProtocolMover::run(
	RRComparerOP comparer,
	RRReporterOP reporter,
	Pose const & pose,
	ScoreFunction const &,
	PackerTask const & packer_task
) {
	// Assume score_function.setup_for_scoring(pose) has already been called.

	Pose working_pose = pose; // deep copy

	if ( !mover_ ) {
		utility_exit_with_message("Attempting to run RotamerRecovery with the 'RRProtocolMover' protocol, but no mover was specified.");
	}

	mover_->apply(working_pose);

	for ( Size ii = 1; ii <= pose.total_residue(); ++ii ) {
		if ( !packer_task.pack_residue(ii) ) continue;
		measure_rotamer_recovery(
			comparer, reporter,
			pose, working_pose,
			pose.residue(ii), working_pose.residue(ii) );
	}
}

} // rotamer_recovery
} // protocols
