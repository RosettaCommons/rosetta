// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/protocols/rotamer_recovery/RRComparerElecDensDiff.cc
/// @author Patrick Conway (ptconway@uw.edu)


// Unit Headers
#include <protocols/rotamer_recovery/RRComparer.hh>
#include <protocols/rotamer_recovery/RRComparerElecDensDiff.hh>

// Project Headers
#include <core/chemical/ResidueType.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/AA.hh>

#include <core/conformation/Residue.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/conformation/util.hh>
#include <core/scoring/rms_util.hh>
#include <basic/Tracer.hh>

#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/electron_density/ElectronDensity.hh>
#include <protocols/electron_density/util.hh>
#include <protocols/electron_density/SetupForDensityScoringMover.hh>

#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

// C++ Headers
#include <string>
#include <sstream>

#include <core/pose/Pose.hh>
#include <utility/vector1.hh>


using std::string;
using std::stringstream;
using std::endl;
using core::Size;
using core::Real;
using core::chemical::num_canonical_aas;
using core::chemical::ResidueType;
using core::chemical::ResidueTypeSet;
using core::conformation::Residue;
using core::conformation::ResidueOP;
using core::conformation::ResidueFactory;
using core::pose::Pose;
using basic::Tracer;

namespace protocols {
namespace rotamer_recovery {

static Tracer TR("protocol.moves.RRComparerElecDensDiff");

RRComparerElecDensDiff::RRComparerElecDensDiff() :
	recovery_threshold_( 0.13 )
{}

RRComparerElecDensDiff::RRComparerElecDensDiff( RRComparerElecDensDiff const & src ) :
	RRComparer(),
	recovery_threshold_( src.recovery_threshold_ )
{}

RRComparerElecDensDiff::~RRComparerElecDensDiff() {}

bool
RRComparerElecDensDiff::measure_rotamer_recovery(
	Pose const & const_pose1,
	Pose const & const_pose2,
	Residue const & res1,
	Residue const & res2,
	Real & score,
	bool & recovered
) {

	if ( res1.aa() != res2.aa()  || res1.nheavyatoms() != res2.nheavyatoms() ) {
		TR << "Cannot measure rotamer recovery of residue " << res1.seqpos() << " because" << endl;
		TR << "\nresidue 1 has type '" << res1.type().name() << "'" << endl;
		TR << "\nresidue 2 has type '" << res2.type().name() << "'" << endl;
		TR << "\nMake sure the protocol to generate the conformations did not 'design' the sequence identity too." << endl;
		score = -1; recovered = false; return false;
	}

	// TODO: Can this restriction be relaxed? What about using 'is_polymer()'?
	if ( res1.aa() > num_canonical_aas ) {
		TR << "WARNING: trying to compare rotamer bins for non-canonical amino acid '" << res1.name() << "'" << endl;
		score = -1; recovered = false; return false;
	}

	core::pose::Pose pose1 = const_pose1;
	core::pose::Pose pose2 = const_pose2;
	core::Size nres = pose1.total_residue();

	protocols::electron_density::SetupForDensityScoringMoverOP dockindens( new protocols::electron_density::SetupForDensityScoringMover );
	dockindens->apply( pose1 );
	dockindens->apply( pose2 );

	// setup density scoring
	core::scoring::electron_density::getDensityMap().set_nres( nres );
	core::scoring::electron_density::getDensityMap().setScoreWindowContext( true );
	if ( basic::options::option[ basic::options::OptionKeys::edensity::sliding_window ].user() ) {
		core::scoring::electron_density::getDensityMap().setWindow( basic::options::option[ basic::options::OptionKeys::edensity::sliding_window ] );
	}

	// score to set energy graph
	core::scoring::ScoreFunctionOP scorefxn ( new core::scoring::ScoreFunction );
	scorefxn->set_weight( core::scoring::fa_atr, 1.0 );
	(*scorefxn)(pose1);
	(*scorefxn)(pose2);

	// electron density correlation
	Real pose1_corr = core::scoring::electron_density::getDensityMap().matchRes( res1.seqpos(), res1, pose1, NULL , false);
	Real pose2_corr = core::scoring::electron_density::getDensityMap().matchRes( res2.seqpos(), res2, pose2, NULL , false);
	Real corr_diff = pose1_corr - pose2_corr;    //if Rosetta fixes an error in native density fitting, count as recovered
	//pose1 must be native

	TR << "type: " << res1.name3() << " seqpos: " << res1.seqpos() << " corr diff: " << corr_diff << std::endl;

	if ( corr_diff > recovery_threshold_ ) {
		recovered = false;
	} else {
		recovered = true;
	}

	return true;
}

string
RRComparerElecDensDiff::get_name() const {
	return "RRComparerElecDensDiff";
}

string
RRComparerElecDensDiff::get_parameters() const {
	return "";
}

void
RRComparerElecDensDiff::set_recovery_threshold(
	Real const recovery_threshold
) {
	recovery_threshold_ = recovery_threshold;
}

Real
RRComparerElecDensDiff::get_recovery_threshold() const {
	return recovery_threshold_;
}

} // rotamer_recovery
} // protocols
