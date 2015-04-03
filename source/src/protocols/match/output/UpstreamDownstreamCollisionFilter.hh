// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/output/UpstreamDownstreamCollisionFilter.hh
/// @brief  Declaration for class to filter matches where the upstream residues collide with the downstream partner
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_output_UpstreamDownstreamCollisionFilter_hh
#define INCLUDED_protocols_match_output_UpstreamDownstreamCollisionFilter_hh

// Unit headers
#include <protocols/match/output/UpstreamDownstreamCollisionFilter.fwd.hh>

// Package headers
#include <protocols/match/BumpGrid.hh>
#include <protocols/match/output/UpstreamCollisionFilter.hh>
#include <protocols/match/output/UpstreamHitCacher.hh>
#include <protocols/match/downstream/DownstreamBuilder.fwd.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/methods/ShortRangeTwoBodyEnergy.fwd.hh>


// Utility headers
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>

#ifdef WIN32
	#include <protocols/match/downstream/DownstreamBuilder.hh>
#endif


namespace protocols {
namespace match {
namespace output {

/// @brief This class is used to detect collisions between upstream residues
/// and downstream poses.
class UpstreamDownstreamCollisionFilter : public MatchCollisionFilter {
public:
	typedef core::Real Real;
	typedef MatchCollisionFilter parent;

public:
	UpstreamDownstreamCollisionFilter(
		std::string filter_name,
		UpstreamHitCacherOP coordinate_cacher );

	virtual
	~UpstreamDownstreamCollisionFilter();

	void
	set_downstream_pose( core::pose::Pose const & downstream_pose );

	void
	set_num_geometric_constraints( Size n_geomcst );

	/// @brief Do not perform any collision detection between an upstream residue and the downstream
	/// pose if there is a chemical bond between them
	void
	set_chemical_bond_from_upstream_to_downstream( Size geomcst_id );

	void
	set_downstream_builder( Size geomcst_id, downstream::DownstreamBuilderCOP builder );

	void set_tolerated_overlap( Real setting );

	/// @brief Returns true if there is at least one placement of the
	/// downstream partner which does not exceed the collision tolerance
	/// with any of the upstream residues.
	virtual
	bool
	passes_filter(
		match const & m
	) const;

	/// @brief Returns true if the specified downstream position does not
	/// collide more than tolerated with any of the upstream residues
	virtual
	bool
	passes_filter(
		match_dspos1 const & m
	) const;

	bool
	passes_hardsphere_filter(
		core::Size geomcst_up,
    core::Size geomcst_down,
		Hit const & upstream_hit,
		Hit const & downstream_hit
	) const;

	inline
	bool
	passes_hardsphere_filter(
		core::Size geomcst,
		Hit const & upstream_hit,
		utility::vector1< core::Vector > const & coords
	) const
	{
		core::conformation::ResidueCOP iires = cacher_->upstream_conformation_for_hit( geomcst, upstream_hit );
		Size ii_first_sc = iires->first_sidechain_atom();
		for ( Size jj = 1; jj <= downstream_pose_->total_residue(); ++jj ) {
			core::conformation::Residue const & jjres = downstream_pose_->residue( jj );
			Real intxn_dis = iires->nbr_radius() + jjres.nbr_radius() + max_overlap_dis_;
			if ( iires->xyz( iires->nbr_atom() ).distance_squared(
					coords[ per_res_atom_ind_[ jj ][ jjres.nbr_atom() ]] ) >
				intxn_dis * intxn_dis ) {
				continue;
			}

			for ( Size kk = ii_first_sc; kk <= iires->nheavyatoms(); ++kk ) {
				ProbeRadius kk_rad = probe_radius_for_atom_type( iires->atom_type_index( kk ) );
				for ( Size ll = 1; ll <= jjres.nheavyatoms(); ++ll ) {
					ProbeRadius ll_rad = probe_radius_for_atom_type( jjres.atom_type_index( ll ) );
					Real minsep = bump_grid()->required_separation_distance( kk_rad, ll_rad );
					if ( iires->xyz( kk ).distance_squared( coords[ per_res_atom_ind_[ jj ][ ll ]]) < minsep * minsep ) {
						return false;
					}
				}
			}
		}//loop over residues of downstream pose
		return true;
	}

private:
	bool passes_etable_filter( match_dspos1 const & m ) const;
	bool passes_hardsphere_filter( match_dspos1 const & m ) const;

private:
	Real max_overlap_dis_;

	utility::vector1< downstream::DownstreamBuilderCOP > dsbuilders_;

	utility::vector1< bool > us_ds_chemical_bond_; // skip residues that form a chemical bond to the downstream partner

	/// NOT THREADSAFE!
	mutable core::pose::PoseOP downstream_pose_;
	mutable utility::vector1< core::Vector > coords_;
	utility::vector1< core::id::AtomID > downstream_atoms_;
	utility::vector1< utility::vector1< Size > > per_res_atom_ind_;
};


}
}
}

#endif
