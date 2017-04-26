// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/downstream/ClassicMatchAlgorithm.hh
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

#ifndef INCLUDED_protocols_match_downstream_ClassicMatchAlgorithm_hh
#define INCLUDED_protocols_match_downstream_ClassicMatchAlgorithm_hh

// Unit headers
#include <protocols/match/downstream/ClassicMatchAlgorithm.fwd.hh>

// Package headers
#include <protocols/match/BumpGrid.fwd.hh>
#include <protocols/match/downstream/ActiveSiteGrid.fwd.hh>
#include <protocols/match/downstream/DownstreamAlgorithm.hh>
#include <protocols/match/downstream/DownstreamBuilder.fwd.hh>
#include <protocols/toolbox/match_enzdes_util/ExternalGeomSampler.hh>

// Project headers
#include <core/types.hh>
#include <core/chemical/ResidueType.fwd.hh>
#include <core/conformation/Residue.fwd.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <list>
#include <string>

#include <utility/fixedsizearray1.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace downstream {

/// @brief Produce hits by hashing building the coordinates of the downstream partner
/// The downstream partner is responsible for building itself from the coordinate frame of
/// three of its atoms.  The ExternalGeomSampler describes the ways to orient the downstream
/// partner given the coordinates of the upstream partner.
class ClassicMatchAlgorithm : public DownstreamAlgorithm
{
public:
	typedef DownstreamAlgorithm parent;
	typedef core::Size Size;
	typedef core::Real Real;
	typedef numeric::HomogeneousTransform< Real > HTReal;

public:
	ClassicMatchAlgorithm( Size geom_cst_id );
	virtual ~ClassicMatchAlgorithm();

	virtual
	DownstreamAlgorithmOP
	clone() const;

	/// @brief Enable a strategy where the first round hits are discarded after they are generated
	/// and then, after the second round completes, they are regenerated but respecting the occ-space
	/// hash, thereby decreasing the memory use dramatically.  That is, the hits for the first
	/// geometric constraint (round 1 hits) are discarded, but their presence in the occupied space
	/// hash is recorded.  Then the hits for the second round are collected, but only the ones that
	/// fall into the regions of 6D where the hits from the first round fell.  At the end of round 2,
	/// the occ-space hash is updated to reflect the regions of 6D where both rounds produced hits.  Then
	/// the round 1 hits are generated a second time, and this time saved.
	void
	set_build_round1_hits_twice();

	virtual
	std::list< Hit >
	build_hits_at_all_positions(
		Matcher & matcher
	);

	std::list< Hit >
	build_and_discard_first_round_hits_at_all_positions(
		Matcher & matcher
	);

	/// @brief Reset the occupied space grid for the matcher so that only those
	/// regions which contain hits from this geometric constraint are marked as occupied.
	virtual
	void
	respond_to_primary_hitlist_change( Matcher & matcher, Size round_just_completed );

	/// @brief Delete hits for this geometric constraint if they fall into
	/// now-empty regions of 6D.  This step can be avoided if the occupied-space-grid's
	/// revision ID has not changed since the last time this function was invoked.
	virtual
	void
	respond_to_peripheral_hitlist_change( Matcher & matcher );

	/// @brief Iterate across the external geom samplers that describe the rigid body orientations
	/// of the downstream partner from the coordinates of the upstream partner.
	virtual
	std::list< Hit >
	build(
		Size const scaffold_build_point_id,
		Size const upstream_conf_id,
		core::conformation::Residue const & upstream_residue
	) const;

	/// @brief This method returns 'false' since the classic match algorithm
	/// builds coordinates of the downstream partner and its hits should
	/// be hashed in 6D to generate matches
	virtual
	bool
	upstream_only() const;

	/// @brief This method returns 'true' since the classic matcher builds the
	/// downstream coordinates from scratch.
	virtual
	bool
	generates_primary_hits() const;


	/// @brief This method should not be invoked on the ClassicMatchAlgorithm,
	/// since it returns "false" in its upstream_only method.
	virtual
	HitPtrListCOP
	hits_to_include_with_partial_match( match_dspos1 const & m ) const;

	virtual
	Size
	n_possible_hits_per_upstream_conformation() const;


	/// @brief This function completes the building of the downstream conformation
	/// once the coordinates of the upstream conformation are known (and deemed
	/// non-colliding or, generally, pass any filter the upstream builder would use).
	std::list< Hit >
	build_from_three_coords(
		Size const which_external_sampler,
		Size const scaffold_build_point_id,
		Size const upstream_conf_id,
		core::conformation::Residue const & upstream_residue
	) const;

public:

	void set_residue_type( core::chemical::ResidueTypeCOP restype );

	void add_external_geom_sampler(
		utility::vector1< toolbox::match_enzdes_util::ExternalGeomSampler > const & sampler,
		Size const exgeom_id,
		std::string const & atom1,
		std::string const & atom2,
		std::string const & atom3,
		DownstreamBuilderCOP downstream_builder
	);

	void clear_external_geom_samplers();

public:
	/// Accessors

	core::chemical::ResidueType const &
	restype() const {
		return *restype_;
	}

	utility::vector1< toolbox::match_enzdes_util::ExternalGeomSampler > const &
	external_sampler( Size external_geom_id ) const
	{
		return external_samplers_[ external_geom_id ];
	}

	Size
	launch_atom( Size external_geom_id, Size which_point ) const {
		return launch_points_[ external_geom_id ][ which_point ];
	}

	Size
	n_external_samplers() const {
		return external_samplers_.size();
	}

private:
	core::chemical::ResidueTypeCOP  restype_;

	// each elelment of the vector is a list of closely-related ExternalGeomSamplers for the same residue
	utility::vector1< utility::vector1< toolbox::match_enzdes_util::ExternalGeomSampler > > external_samplers_;
	utility::vector1< utility::fixedsizearray1< Size, 3 > > launch_points_;
	utility::vector1< DownstreamBuilderCOP > dsbuilders_;
	utility::vector1< Size > exgeom_ids_;

	Size occspace_rev_id_at_last_update_;

	/// @brief control whether this algorithm expects to rebuild round 1 hits
	/// after round 2 completes.  Only used if the geom_cst_id() for this instance
	/// is 1.
	bool build_round1_hits_twice_;
	bool completed_first_round1_hit_building_;
};

}
}
}

#endif
