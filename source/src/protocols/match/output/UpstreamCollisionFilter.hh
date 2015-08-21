// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/output/UpstreamCollisionFilter.hh
/// @brief  Declaration for class to filter matches where the upstream residues collide.
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_output_UpstreamCollisionFilter_hh
#define INCLUDED_protocols_match_output_UpstreamCollisionFilter_hh

// Unit headers
#include <protocols/match/output/UpstreamCollisionFilter.fwd.hh>

// Package headers
#include <protocols/match/BumpGrid.hh>
#include <protocols/match/output/MatchFilter.hh>
#include <protocols/match/output/UpstreamHitCacher.fwd.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/methods/ShortRangeTwoBodyEnergy.fwd.hh>


// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace output {


class MatchCollisionFilter : public MatchFilter {

public:
	typedef core::Real Real;

public:
	//MatchCollisionFilter( std::string filter_name );
	MatchCollisionFilter( std::string filter_name, UpstreamHitCacherOP coordinate_cacher );

	virtual
	~MatchCollisionFilter();

	virtual
	bool
	passes_filter(
		match const & m
	) const = 0;

	virtual
	bool
	passes_filter(
		match_dspos1 const & m
	) const = 0;

	//setters
	void set_filter_by_lj( bool setting );
	void set_lj_cutoff( Real setting );
	void set_lj_atr_weight( Real setting );
	void set_lj_rep_weight( Real setting );
	void set_lj_sol_weight( Real setting );

	virtual
	void set_tolerated_overlap( Real setting );

	//getters
	bool filter_by_lj() const;
	Real wfa_atr() const;
	Real wfa_rep() const;
	Real wfa_sol() const;
	Real lj_cutoff() const;
	Real tolerated_overlap() const;

	core::pose::PoseCOP empty_pose() const;
	core::scoring::ScoreFunctionOP empty_sfxn() const;
	core::scoring::methods::ShortRangeTwoBodyEnergyCOP etable_energy() const;
	BumpGridCOP bump_grid() const;


protected:
	UpstreamHitCacherOP cacher_;

private:
	bool filter_by_lj_;
	Real wfa_atr_;
	Real wfa_rep_;
	Real wfa_sol_;
	Real lj_cutoff_;
	Real tolerated_overlap_;

	core::pose::PoseOP empty_pose_;
	core::scoring::ScoreFunctionOP empty_sfxn_;

	core::scoring::methods::ShortRangeTwoBodyEnergyOP   etable_energy_;
	BumpGridOP bump_grid_;


};

/// @brief This class is used to detect collisions between the upstream residues
/// and filter out matches that have too much collision.  It can perform either
/// hard-sphere collision detection, or score-function (Etable) driven collision
/// detection.  Four command-line flags are read by the MatcherTask to initialize
/// this class:
/// match::filter_colliding_upstream_residues
/// match::upstream_residue_collision_tolerance
/// match::upstream_residue_collision_score_cutoff
/// match::upstream_residue_collision_Wfa_atr
/// match::upstream_residue_collision_Wfa_rep
/// match::upstream_residue_collision_Wfa_sol
class UpstreamCollisionFilter : public MatchCollisionFilter {
public:
	typedef core::Real Real;
	typedef MatchCollisionFilter parent;

public:
	//UpstreamCollisionFilter( std::string filter_name );
	UpstreamCollisionFilter( std::string filter_name, UpstreamHitCacherOP coordinate_cacher );

	virtual
	~UpstreamCollisionFilter();

	/// @brief Returns true if the given match does not contain too much residue-pair collision.
	virtual
	bool
	passes_filter(
		match const & m
	) const;

	/// @brief Returns true if the given match does not contain too much residue-pair collision.
	virtual
	bool
	passes_filter(
		match_dspos1 const & m
	) const;

	bool
	passes_hardsphere_filter(
		core::Size geomcst_a,
		core::Size geomcst_b,
		Hit const & upstream_hit_a,
		Hit const & upstream_hit_b
	) const;

	inline
	bool
	passes_hardsphere_filter(
		core::conformation::Residue const & resA,
		core::conformation::Residue const & resB
	) const
	{
		for ( Size ii = resA.first_sidechain_atom(); ii <= resA.nheavyatoms(); ++ii ) {
			ProbeRadius ii_rad = probe_radius_for_atom_type( resA.atom_type_index( ii ) );
			for ( Size jj = resB.first_sidechain_atom(); jj <= resB.nheavyatoms(); ++jj ) {
				ProbeRadius jj_rad = probe_radius_for_atom_type( resB.atom_type_index( jj ) );
				Real minsep = bump_grid()->required_separation_distance( ii_rad, jj_rad );
				if ( resA.xyz( ii ).distance_squared( resB.xyz( jj )) < minsep * minsep ) {
					return false;
				}
			}
		}
		return true;
	}

private:

};


}
}
}

#endif
