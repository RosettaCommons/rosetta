// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/methods/MotifDockEnergy.cc
/// @brief  Adaptation of Motif score for Docking
/// @author Nick Marze (nickmarze@gmail.com)


// Unit headers
#include <core/scoring/methods/MotifDockEnergy.hh>
#include <core/scoring/methods/MotifDockEnergyCreator.hh>

// Package headers
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/DerivVectorPair.hh>

// Project headers
#include <core/chemical/AA.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/VariantType.hh>

#include <numeric/xyzVector.hh>
#include <numeric/xyzTransform.hh>
#include <numeric/deriv/distance_deriv.hh>
#include <numeric/deriv/angle_deriv.hh>

#include <core/scoring/EnergyMap.hh>
#include <utility/vector1.hh>
#include <basic/Tracer.hh>

#include <core/pose/motif/reference_frames.hh>
#include <core/scoring/motif/motif_hash_stuff.hh>

static basic::Tracer TR( "MotifDock" );

typedef numeric::xyzTransform<core::Real> Xform;


namespace core {
namespace scoring {
namespace methods {

methods::EnergyMethodOP
MotifDockEnergyCreator::create_energy_method(
	methods::EnergyMethodOptions const &
) const {
	return methods::EnergyMethodOP( new MotifDockEnergy );
}

/// @brief Return the set of score types claimed by the EnergyMethod
/// this EnergyMethodCreator creates in its create_energy_method() function
ScoreTypes
MotifDockEnergyCreator::score_types_for_method() const {
	ScoreTypes sts;
	sts.push_back( motif_dock );
	return sts;
}

/// c-tor
MotifDockEnergy::MotifDockEnergy() :
	parent( methods::EnergyMethodCreatorOP( new MotifDockEnergyCreator ) )
{ }


/// clone
EnergyMethodOP
MotifDockEnergy::clone() const {
	return EnergyMethodOP( new MotifDockEnergy() );
}


/////////////////////////////////////////////////////////////////////////////
// score
/////////////////////////////////////////////////////////////////////////////

void
MotifDockEnergy::residue_pair_energy(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	pose::Pose const & pose,
	ScoreFunction const &,
	EnergyMap & emap
) const {

	Real score = 0;
	Size ir = rsd1.seqpos();
	Size jr = rsd2.seqpos();
	Size chain1 = rsd1.chain();
	Size chain2 = rsd2.chain();
	Xform const xb1 = core::pose::motif::get_backbone_reference_frame(pose,ir);
	Xform const xb2 = core::pose::motif::get_backbone_reference_frame(pose,jr);
	if ( xb1.t.distance_squared(xb2.t) > 100.0 ) {
		emap[ motif_dock ] = score;
	} else if ( chain1 == chain2 ) {
		emap[ motif_dock ] = score;
	} else {
		Xform const xbb = xb1.inverse() * xb2;
		char const & ss1 = pose.secstruct(ir);
		char const & ss2 = pose.secstruct(jr);
		char const & aa1 = pose.residue(ir).name1();
		char const & aa2 = pose.residue(jr).name1();
		core::scoring::motif::MotifHashManager & mman = *core::scoring::motif::MotifHashManager::get_instance();
		Real bb_motif = 0;
		core::scoring::motif::XformScoreCOP xs_bb_fxn1 = mman.get_xform_score_BB_BB(ss1,ss2,aa1,aa2);
		core::scoring::motif::XformScoreCOP xs_bb_fxn2 = mman.get_xform_score_BB_BB(ss2,ss1,aa2,aa1);
		if ( xs_bb_fxn1 ) bb_motif += xs_bb_fxn1->score_of_bin(xbb          .rt6());
		if ( xs_bb_fxn2 ) bb_motif += xs_bb_fxn2->score_of_bin(xbb.inverse().rt6());
		//TR << "bb_motif " << bb_motif << std::endl;
		score += bb_motif;
		//if( bb_motif == 0 ) score -= 0.5;

		emap[ motif_dock ] = -score;
	}
	//TR << "residue pair: A: " << ir << " B: " << jr << " score " << score << std::endl;

}


void
MotifDockEnergy::eval_intrares_energy(
	conformation::Residue const &,
	pose::Pose const &,
	ScoreFunction const &,
	EnergyMap &
) const {
	return;
}

Distance
MotifDockEnergy::atomic_interaction_cutoff() const {return 8;}

bool
MotifDockEnergy::defines_intrares_energy( EnergyMap const & ) const { return false; }

core::Size
MotifDockEnergy::version() const
{
	return 1; // Initial versioning
}

}
}
}
