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

/// @details Auto-generated virtual destructor
RRComparer::~RRComparer() = default;

static Tracer TR("protocol.moves.RRComparer");

RRComparerRotBins::RRComparerRotBins() :
	recovery_threshold_(0)
{}

RRComparerRotBins::RRComparerRotBins( RRComparerRotBins const & src) :
	RRComparer(),
	recovery_threshold_(src.recovery_threshold_)
{}

RRComparerRotBins::~RRComparerRotBins() = default;

void
RRComparerRotBins::set_recovery_threshold(
	Real const recovery_threshold
) {
	recovery_threshold_ = recovery_threshold;
}

void
RRComparerRotBins::set_absolute_threshold(
	Real const absolute_threshold
) {
	absolute_threshold_ = absolute_threshold;
}


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


	if ( res1.aa() != res2.aa() ) {
		TR.Fatal << "Cannot measure rotamer recovery because" << endl;
		TR.Fatal << "residue 1 has type '" << res1.type().name() << "'" << endl;
		TR.Fatal << "residue 2 has type '" << res2.type().name() << "'" << endl;
		TR.Fatal << "Make sure the protocol to generate the conformations did not 'design' the sequence identity too." << endl;
		utility_exit();
	}

	if ( res1.aa() > num_canonical_aas ) {
		TR.Warning << "trying to compare rotamer bins for non-canonical amino acid '" << res1.name() << "'" << endl;
		return false;
	}

	RotVector res1_rotbins, res2_rotbins;
	rotamer_from_chi( res1, res1_rotbins );
	rotamer_from_chi( res2, res2_rotbins );

	Size chi_match(0);
	for ( Size chi_index=1; chi_index <= res1_rotbins.size(); ++chi_index ) {
		if ( res1_rotbins[ chi_index ] == res2_rotbins[ chi_index ] ) {
			chi_match++;
		}
	}
	score = res1_rotbins.size() - chi_match;
	recovered = (score <= recovery_threshold_);
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


RRComparerChiDiff::RRComparerChiDiff() : tolerance_( 20.0 ), limit_chi_angles_( false ), max_chi_considered_ ( 100 ) {}

RRComparerChiDiff::RRComparerChiDiff( RRComparerChiDiff const & src ) :
	RRComparer(),
	tolerance_( src.tolerance_ ),
	limit_chi_angles_( src.limit_chi_angles_ ),
	max_chi_considered_( src.max_chi_considered_ )
{}

RRComparerChiDiff::~RRComparerChiDiff() = default;

void
RRComparerChiDiff::set_recovery_threshold( core::Real const setting ) { tolerance_ = setting; }

void
RRComparerChiDiff::set_absolute_threshold( core::Real const setting ) { absolute_threshold_ = setting; }

// AS March 2013: enable comparison of a limited set of chi angles, no matter how many the residue actually has
void
RRComparerChiDiff::set_max_chi_considered( core::Size const max_chi ) {
	limit_chi_angles_ = true;
	max_chi_considered_ = max_chi;
}


/// @details measure the rotamer recovery by comparing rotamer bins
/// @return  true, if the measurement was successful, false otherwise
bool
RRComparerChiDiff::measure_rotamer_recovery(
	Pose const & /* pose1 */,
	Pose const & /* pose2 */,
	Residue const & res1,
	Residue const & res2,
	Real & score,
	bool & recovered
) {

	using core::pack::dunbrack::subtract_chi_angles;

	if ( res1.aa() != res2.aa() ) {
		TR.Fatal << "Cannot measure rotamer recovery because" << endl;
		TR.Fatal << "residue 1 has type '" << res1.type().name() << "'" << endl;
		TR.Fatal << "residue 2 has type '" << res2.type().name() << "'" << endl;
		TR.Fatal << "Make sure the protocol to generate the conformations did not 'design' the sequence identity too." << endl;
		utility_exit();
	}

	if ( res1.aa() > num_canonical_aas ) {
		TR.Warning << "trying to compare rotamer bins for non-canonical amino acid '" << res1.name() << "'" << endl;
		return false;
	}

	score = 0;
	recovered=true;
	Size max_chi = res1.nchi();
	if ( limit_chi_angles_ ) {
		max_chi = std::min(max_chi, max_chi_considered_);
	}

	for ( Size chi_index=1; chi_index <= max_chi; ++chi_index ) {
		if ( res1.type().chi_2_proton_chi( chi_index ) == 0 ) { // ignore proton chi (tyr,ser,thr)
			Real chidiff = std::abs( subtract_chi_angles( res1.chi(chi_index), res2.chi(chi_index), res1.aa(), chi_index ));
			if ( score < chidiff ) { score = chidiff; }
			if ( chidiff > tolerance_ ) { recovered = false; }
		}
	}

	return true;
}

string
RRComparerChiDiff::get_name() const {
	return "RRComparerChiDiff";
}

string
RRComparerChiDiff::get_parameters() const {
	return "";
}

} // rotamer_recovery
} // protocols
