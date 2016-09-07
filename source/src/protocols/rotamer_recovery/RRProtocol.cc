// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rotamer_recovery/RRProtocol.cc
/// @author Matthew O'Meara (mattjomeara@gmail.com)

// Unit Headers
#include <protocols/rotamer_recovery/RRProtocol.hh>
#include <protocols/rotamer_recovery/RRComparer.hh>
#include <protocols/rotamer_recovery/RRReporter.hh>

// Project Headers
#include <basic/Tracer.hh>

// C++ Headers

#include <utility/vector1.hh>


using std::string;
using std::endl;
using core::Size;
using core::Real;
using core::conformation::Residue;
using core::pose::Pose;
using basic::Tracer;

namespace protocols {
namespace rotamer_recovery {

/// @details Auto-generated virtual destructor
RRProtocol::~RRProtocol() = default;

static Tracer TR("protocol.moves.RRProtocol");

bool
RRProtocol::measure_rotamer_recovery(
	RRComparerOP comparer,
	RRReporterOP reporter,
	Pose const & pose1,
	Pose const & pose2,
	Residue const & res1,
	Residue const & res2
) {

	Real score;
	bool recovered;
	if ( comparer->measure_rotamer_recovery(
			pose1, pose2, res1, res2, score, recovered) ) {
		reporter->report_rotamer_recovery(
			pose1, pose2, res1, res2, score, recovered );
	}
	return recovered;
}

} // rotamer_recovery
} // protocols
