// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/simple_moves/TorsionDOFMover.cc
/// @brief TorsionDOFMover methods implemented
/// @author Steven Lewis

// Unit Headers
#include <protocols/simple_moves/TorsionDOFMover.hh>

// Project Headers
#include <core/conformation/Conformation.hh>
#include <core/pose/Pose.hh>

//#include <core/scoring/methods/MMTorsionEnergy.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreType.hh>

// Numeric Headers
#include <numeric/conversions.hh> //degrees-radians
#include <numeric/random/random.hh>

// Utility Headers
#include <basic/Tracer.hh>
#include <core/types.hh>

// C++ Headers
#include <string>

#include <utility/vector1.hh>

//Auto Headers
#include <core/kinematics/AtomTree.hh>

using basic::Error;
using basic::Warning;

static THREAD_LOCAL basic::Tracer TR( "protocols.simple_moves.TorsionDOFMover" );

namespace protocols {
namespace simple_moves {

TorsionDOFMover::TorsionDOFMover() :
	protocols::moves::Mover(),
	atom1_(core::id::BOGUS_ATOM_ID),
	atom2_(core::id::BOGUS_ATOM_ID),
	atom3_(core::id::BOGUS_ATOM_ID),
	atom4_(core::id::BOGUS_ATOM_ID),
	upper_angle_(0),
	lower_angle_(0),
	check_MMT_(false),
	mmt_(/* NULL */),
	temp_(0),
	tries_(0)
{ protocols::moves::Mover::type( "TorsionDOFMover" ); }

/// @details random angle constructor.  Magic numbers 180 and -179.9999999... maintain the uniform range.  I'm sure there's a better way to get [180, -180) but I can't figure out what it is.
TorsionDOFMover::TorsionDOFMover(
	core::id::AtomID const & atom1,
	core::id::AtomID const & atom2,
	core::id::AtomID const & atom3,
	core::id::AtomID const & atom4
) :
	protocols::moves::Mover(),
	atom1_(atom1),
	atom2_(atom2),
	atom3_(atom3),
	atom4_(atom4),
	upper_angle_(180.0),
	lower_angle_(-179.9999999999999999999999999999999999999999999999),
	check_MMT_(false),
	mmt_(/* NULL */),
	temp_(0.8),
	tries_(1)
{ protocols::moves::Mover::type( "TorsionDOFMover" ); }

/// @details range of angles constructor - takes DEGREES not RADIANS.
TorsionDOFMover::TorsionDOFMover(
	core::id::AtomID const & atom1,
	core::id::AtomID const & atom2,
	core::id::AtomID const & atom3,
	core::id::AtomID const & atom4,
	core::Angle const upper,
	core::Angle const lower
) :
	protocols::moves::Mover(),
	atom1_(atom1),
	atom2_(atom2),
	atom3_(atom3),
	atom4_(atom4),
	upper_angle_(upper),
	lower_angle_(lower),
	check_MMT_(false),
	mmt_(/* NULL */),
	temp_(0.8),
	tries_(1)
{ protocols::moves::Mover::type( "TorsionDOFMover" ); }

/// @details particular angle constructor - takes DEGREES not RADIANS.
TorsionDOFMover::TorsionDOFMover( core::id::AtomID const & atom1, core::id::AtomID const & atom2, core::id::AtomID const & atom3, core::id::AtomID const & atom4, core::Angle const angle )
: protocols::moves::Mover(),
	atom1_(atom1),
	atom2_(atom2),
	atom3_(atom3),
	atom4_(atom4),
	upper_angle_(angle),
	lower_angle_(angle),
	check_MMT_(false),
	mmt_(/* NULL */),
	temp_(0.8),
	tries_(1)
{ protocols::moves::Mover::type( "TorsionDOFMover" ); }

TorsionDOFMover::~TorsionDOFMover()= default;

void TorsionDOFMover::apply( core::pose::Pose & pose ){

	if ( !(pose.atom_tree().torsion_angle_dof_id( atom1_, atom2_, atom3_, atom4_ ).valid()) ) {

		Warning() << "In TorsionDOFMover, atoms not valid against pose; atoms:"
			<< " atom1 " << atom1_
			<< " atom2 " << atom2_
			<< " atom3 " << atom3_
			<< " atom4 " << atom4_ << std::endl;
		return;
	}

	//TR << "atoms:" << " atom1 " << atom1_ << " atom2 " << atom2_ << " atom3 " << atom3_ << " atom4 " << atom4_ << std::endl;

	//if we want this score, fill the pointer!
	if ( check_MMT_ && !mmt_ ) {
		mmt_ = core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction );
		mmt_->set_weight(core::scoring::mm_twist, 1.0);
	}

	//if scoring, pre-score
	core::Energy pre_score(0), post_score(0);
	if ( check_MMT_ ) pre_score = score_torsion(pose);

	//n_tries loop: continue rotating until a good angle is found
	core::Size ntries(1);
	for ( ; ntries <= tries_; ++ntries ) {

		//make a move
		core::Angle const pre_torsion(pose.atom_tree().torsion_angle(atom1_, atom2_, atom3_, atom4_));
		pose.conformation().set_torsion_angle( atom1_, atom2_, atom3_, atom4_, pre_torsion+calc_angle() );

		//if scoring, post-score and boltzmann
		if ( check_MMT_ ) {
			post_score = score_torsion(pose);

			//if scoring, decide if try again
			if ( boltzmann( pre_score, post_score ) ) break;
			else pose.conformation().set_torsion_angle( atom1_, atom2_, atom3_, atom4_, pre_torsion );
		}

	}

	if ( ntries > tries_ ) {
		Error() << "TorsionDOFMover gave up after " << tries_ << " attempts, no move made" << std::endl;
	}

	//TR << pre_score << " " << post_score << " " <<  pose.atom_tree().torsion_angle(atom1_, atom2_, atom3_, atom4_) << std::endl;

	//TR << pose.atom_tree().torsion_angle(atom1_, atom2_, atom3_, atom4_) << std::endl;
	return;
}//apply

std::string
TorsionDOFMover::get_name() const {
	return "TorsionDOFMover";
}

/// @brief calculate angle for perturbation - call to RNG
core::Angle TorsionDOFMover::calc_angle()
{ return numeric::conversions::radians(lower_angle_ + ((upper_angle_ - lower_angle_) * numeric::random::rg().uniform())); }

/// @brief calculate mmt score for the moving bond
///This is the stupidest possible method - score the whole pose.  It would be much better if this could directly use MMTorsionEnergy to calculate about just the one bond in question, but I can't figure out how to reliably look up exactly the set of residues modified by this torsion (remember that changing this 4-body torsion affects other 4-bodies with a shared central bond).
core::Energy TorsionDOFMover::score_torsion(core::pose::Pose & pose){ return ((*mmt_)(pose)); }

/// @brief boltzmann calculation - is the new score acceptable?
bool TorsionDOFMover::boltzmann( core::Energy const pre_score, core::Energy const post_score ){

	//TR << pre_score << " " << post_score << std::endl;

	//borrowed from BackboneMover check_rama
	if ( post_score > pre_score ) {
		core::Real const boltz_factor = ((pre_score - post_score)/temp_);
		core::Real const probability = std::exp(std::max(core::Real(-40.0),boltz_factor) );
		if ( numeric::random::rg().uniform() >= probability ) return false;
	}

	//TR << "accepting" << std::endl;

	return true;
}

}//moves
}//protocols
