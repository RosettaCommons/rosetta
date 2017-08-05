// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/rotamer_recovery/RRComparer.cc
/// @author Matthew O'Meara (mattjomeara@gmail.com)
/// Adapted from:
/// protocols::optimize_weights::IterativeOptEDriver::measure_rotamer_recovery()
/// and apps::pilot::doug::rotamer_prediction_benchmark()

// Unit Headers
#include <protocols/rotamer_recovery/RRReporter.hh>

// Project Headers
#include <core/types.hh>
#include <basic/Tracer.hh>
#include <utility/vector1.hh>

//C++ Headers
#include <ostream>

namespace protocols {
namespace rotamer_recovery {

/// @details Auto-generated virtual destructor
RRReporter::~RRReporter() = default;

using std::endl;
using std::ostream;
using std::string;
using core::Real;
using core::pose::Pose;
using core::conformation::Residue;
using basic::Tracer;

static THREAD_LOCAL Tracer TR("protocol.rotamer_recovery.RRReporter");

RRReporterSimple::RRReporterSimple() :
	residues_considered_( 0 ),
	rotamers_recovered_( 0 )
{}

RRReporterSimple::RRReporterSimple( RRReporterSimple const & src ) :
	RRReporter(),
	residues_considered_(src.residues_considered_),
	rotamers_recovered_(src.rotamers_recovered_)
{}

RRReporterSimple::~RRReporterSimple() = default;

void
RRReporterSimple::reset_recovery(){
	residues_considered_=0;
	rotamers_recovered_=0;
}

void
RRReporterSimple::report_rotamer_recovery(
	Pose const & /*pose1*/,
	Pose const & /*pose2*/,
	Residue const & /*res1*/,
	Residue const & /*res2*/,
	Real const /*score*/,
	bool const recovered
) {
	residues_considered_++;
	rotamers_recovered_ += recovered;
}

void
RRReporterSimple::show( ostream & out ) const {

	out
		<< "Recovered " << rotamers_recovered_
		<< " at " << residues_considered_ << " residues considered"
		<< " for a recovery rate of " << recovery_rate() << "." << endl;
}

void
RRReporterSimple::show( ) const {

	TR << "Recovered " << rotamers_recovered_ << " rotamers"
		<< " at " << residues_considered_ << " residues considered"
		<< " for a recovery rate of " << recovery_rate() << "." << endl;
}


Real
RRReporterSimple::recovery_rate() const {
	return Real(rotamers_recovered_) / Real(residues_considered_);
}

} // rotamer_recovery
} // protocols

