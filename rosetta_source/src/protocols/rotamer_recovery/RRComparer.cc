// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/rotamer_recovery/RRComparer.cc
/// @author Matthew O'Meara (mattjomeara@gmail.com)
/// Adapted from:
/// protocols::optimize_weights::IterativeOptEDriver::measure_rotamer_recovery()
/// and apps::pilot::doug::rotamer_prediction_benchmark()


// Unit Headers
#include <protocols/rotamer_recovery/RRComparer.hh>

// Project Headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <basic/Tracer.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>


using std::string;
using std::endl;
using core::Size;
using core::Real;
using core::chemical::num_canonical_aas;
using core::conformation::Residue;
using core::pose::Pose;
using core::pack::dunbrack::RotVector;
using core::pack::dunbrack::rotamer_from_chi;
using basic::Tracer;

namespace protocols {
namespace rotamer_recovery {

static Tracer TR("protocol.moves.RRComparer");

RRComparerRotBins::RRComparerRotBins() {}

RRComparerRotBins::RRComparerRotBins( RRComparerRotBins const & ) :
	RRComparer()
{}

RRComparerRotBins::~RRComparerRotBins() {}

/// @details measure the rotamer recovery by comparing rotamer bins
/// @return  true, if the measurement was successful, false otherwise
bool
RRComparerRotBins::measure_rotamer_recovery(
	Pose const & /* pose1 */,
	Pose const & /* pose2 */,
	Residue const & res1,
	Residue const & res2,
	Real & score,
	bool & recovered
) {


	if( res1.aa() != res2.aa() ) {
		TR << "Cannot measure rotamer recovery because" << endl;
		TR << "residue 1 has type '" << res1.type().name() << "'" << endl;
		TR << "residue 2 has type '" << res2.type().name() << "'" << endl;
		TR << "Make sure the protocol to generate the conformations did not 'design' the sequence identity too." << endl;
		utility_exit();
	}

	if( res1.aa() > num_canonical_aas ){
		TR << "WARNING: trying to compare rotamer bins for non-canonical amino acid '" << res1.name() << "'" << endl;
		return false;
	}

	RotVector res1_rotbins, res2_rotbins;
	rotamer_from_chi( res1, res1_rotbins );
	rotamer_from_chi( res2, res2_rotbins );

	Size chi_match(0);
	for ( Size chi_index=1; chi_index <= res1_rotbins.size(); ++chi_index ){
		if ( res1_rotbins[ chi_index ] == res2_rotbins[ chi_index ] ) {
			chi_match++;
		}
	}
	score = res1_rotbins.size() - chi_match;
	recovered = (score == 0);
	return true;
}

string
RRComparerRotBins::get_name() const {
	return "RRComparerRotBins";
}

string
RRComparerRotBins::get_parameters() const {
	return "";
}

} // rotamer_recovery
} // protocols
