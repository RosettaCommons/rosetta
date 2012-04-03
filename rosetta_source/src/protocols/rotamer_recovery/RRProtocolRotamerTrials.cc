// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/rotamer_recovery/RRProtocol.cc
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/rotamer_recovery/RRProtocolRotamerTrials.hh>

// Project Headers
#include <protocols/rotamer_recovery/RRComparer.hh>
#include <protocols/rotamer_recovery/RRReporter.hh>

// Platform Headers
#include <basic/Tracer.hh>
// AUTO-REMOVED #include <core/chemical/ResidueType.hh>
#include <core/chemical/AA.hh>
// AUTO-REMOVED #include <core/conformation/Residue.hh>
#include <core/pack/rotamer_trials.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>

// C++ Headers
#include <string>

//Auto Headers
#include <utility/vector1.hh>
using std::string;
using core::Size;
using core::pack::task::PackerTask;
using core::pack::task::PackerTaskOP;
using core::pose::Pose;
using core::scoring::ScoreFunction;
using core::scoring::getScoreFunction;
using basic::Tracer;

namespace protocols {
namespace rotamer_recovery {

static Tracer TR("protocol.moves.RRProtocolRotamerTrials");

RRProtocolRotamerTrials::RRProtocolRotamerTrials() {}

RRProtocolRotamerTrials::RRProtocolRotamerTrials(RRProtocolRotamerTrials const & ) {}

RRProtocolRotamerTrials::~RRProtocolRotamerTrials() {}

string
RRProtocolRotamerTrials::get_name() const {
	return "RRProtocolRotamerTrials";
}

string
RRProtocolRotamerTrials::get_parameters() const {
	return "";
}


/// @details For each residue, minimize it, and measure the rotamer
/// compared to where it started
void
RRProtocolRotamerTrials::run(
  RRComparerOP comparer,
  RRReporterOP reporter,
  Pose const & pose,
	ScoreFunction const & score_function,
  PackerTask const & packer_task
) {

	using core::chemical::aa_unk;

	// Assume score_function.setup_for_scoring(pose) has already been called.

	PackerTaskOP one_res_task( packer_task.clone() );

	// I don't know if rtmin looks at more than pack_residue(..)
	one_res_task->temporarily_fix_everything();

	// For each residue in the packer task,
	// rtmin residue -> and measure recovery
	for( Size ii = 1; ii <= pose.total_residue(); ++ii ){
		if ( !packer_task.pack_residue(ii) ) continue;
		Pose working_pose = pose;  // deep copy
		one_res_task->temporarily_set_pack_residue( ii, true );
		core::pack::rotamer_trials( working_pose, score_function, one_res_task );
		measure_rotamer_recovery(
			comparer, reporter,
			pose, working_pose,
			pose.residue(ii), working_pose.residue(ii) );
		one_res_task->temporarily_set_pack_residue( ii, false );
	}
}

} // rotamer_recovery
} // protocols
