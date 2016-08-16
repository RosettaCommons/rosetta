// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rotamer_recovery/RRProtocolMinPack.cc
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/rotamer_recovery/RRProtocolMinPack.hh>

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
#include <utility/vector1.hh>
using std::string;
using core::Size;
using core::pack::min_pack;
using core::pose::Pose;
using core::scoring::ScoreFunction;
using core::scoring::get_score_function;
using core::pack::task::PackerTask;
using core::pack::task::PackerTaskOP;
using basic::Tracer;

namespace protocols {
namespace rotamer_recovery {

static Tracer TR("protocol.moves.RRProtocolMinPack");

RRProtocolMinPack::RRProtocolMinPack() {}

RRProtocolMinPack::RRProtocolMinPack( RRProtocolMinPack const & ) :
	RRProtocol()
{}

RRProtocolMinPack::~RRProtocolMinPack() {}

string
RRProtocolMinPack::get_name() const {
	return "RRProtocolMinPack";
}

string
RRProtocolMinPack::get_parameters() const {
	return "";
}


/// @details apply MinPack and measure rotamer recovery for each residue
void
RRProtocolMinPack::run(
	RRComparerOP comparer,
	RRReporterOP reporter,
	Pose const & pose,
	ScoreFunction const & score_function,
	PackerTask const & packer_task
) {
	// Assume score_function.setup_for_scoring(pose) has already been called.

	Pose working_pose = pose; // deep copy
	min_pack(working_pose, score_function, packer_task.get_self_ptr());

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
