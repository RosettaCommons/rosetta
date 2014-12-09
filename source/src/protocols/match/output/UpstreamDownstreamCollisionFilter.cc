// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/output/UpstreamDownstreamCollisionFilter.hh
/// @brief  Implementation for class to filter matches where the upstream residues collide.
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/output/UpstreamDownstreamCollisionFilter.hh>

// Package headers
#include <protocols/match/BumpGrid.hh>
#include <protocols/match/downstream/DownstreamBuilder.hh>
#include <protocols/match/output/UpstreamHitCacher.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include	<core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
// AUTO-REMOVED #include <core/scoring/ScoringManager.hh>
// AUTO-REMOVED #include <core/scoring/etable/EtableEnergy.hh>
// AUTO-REMOVED #include <core/scoring/methods/EnergyMethodOptions.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <iostream>

#include <core/id/AtomID.hh>
#include <core/scoring/methods/ShortRangeTwoBodyEnergy.hh>
#include <protocols/match/Hit.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace output {

UpstreamDownstreamCollisionFilter::UpstreamDownstreamCollisionFilter(
	std::string filter_name,
	UpstreamHitCacherOP coordinate_cacher
) :
	parent( filter_name, coordinate_cacher ),
	max_overlap_dis_(0.0)
{
	std::cout << "Created UpstreamDownstreamCollisionFilter" << std::endl;
}

UpstreamDownstreamCollisionFilter::~UpstreamDownstreamCollisionFilter()
{}

void
UpstreamDownstreamCollisionFilter::set_downstream_pose( core::pose::Pose const & downstream_pose )
{
	downstream_pose_ = core::pose::PoseOP( new core::pose::Pose( downstream_pose ) );
	Size count_atoms( 0 );
	per_res_atom_ind_.resize( downstream_pose.total_residue() );
	for ( Size ii = 1; ii <= downstream_pose.total_residue(); ++ii ) {
		Size const ii_natoms = downstream_pose.residue( ii ).natoms();
		per_res_atom_ind_[ ii ].resize( ii_natoms );
		for ( Size jj = 1; jj <= ii_natoms; ++jj ) {
			per_res_atom_ind_[ ii ][ jj ] = ++count_atoms;
		}
	}
	downstream_atoms_.resize( count_atoms );
	coords_.resize( count_atoms );
	count_atoms = 0;
	for ( Size ii = 1; ii <= downstream_pose.total_residue(); ++ii ) {
		Size const ii_natoms = downstream_pose.residue( ii ).natoms();
		for ( Size jj = 1; jj <= ii_natoms; ++jj ) {
			downstream_atoms_[ ++count_atoms ] = core::id::AtomID( jj, ii );
		}
	}
}

void
UpstreamDownstreamCollisionFilter::set_num_geometric_constraints( Size n_geomcst )
{
	dsbuilders_.resize( n_geomcst );
	us_ds_chemical_bond_.resize( n_geomcst );
	std::fill( us_ds_chemical_bond_.begin(), us_ds_chemical_bond_.end(), false );
}

void
UpstreamDownstreamCollisionFilter::set_chemical_bond_from_upstream_to_downstream( Size geomcst_id )
{
	us_ds_chemical_bond_[ geomcst_id ] = true;
}


void
UpstreamDownstreamCollisionFilter::set_downstream_builder(
	Size geomcst_id,
	downstream::DownstreamBuilderCOP builder
)
{
	runtime_assert( dsbuilders_.size() >= geomcst_id && geomcst_id > 0 );
	dsbuilders_[ geomcst_id ] = builder;
}



bool
UpstreamDownstreamCollisionFilter::passes_filter(
	match const & m
) const
{
	for ( Size ii = 1; ii <= m.size(); ++ii ) {
		if ( ! dsbuilders_[ ii ] ) continue;
		if ( passes_filter( match_dspos1( m, ii ) ) )	return true;
	}
	return false;
}

bool
UpstreamDownstreamCollisionFilter::passes_filter(
	match_dspos1 const & m
) const
{
	if ( filter_by_lj() ) {
		return passes_etable_filter( m );
	} else {
		return passes_hardsphere_filter( m );
	}
}

void UpstreamDownstreamCollisionFilter::set_tolerated_overlap( Real setting )
{
	parent::set_tolerated_overlap( setting );
	max_overlap_dis_ = 0;
	for ( Size ii = 1; ii <= n_probe_radii; ++ii ) {
		for ( Size jj = ii; jj <= n_probe_radii; ++jj ) {
			if ( bump_grid()->required_separation_distance(
					ProbeRadius( ii ), ProbeRadius( jj ) ) > max_overlap_dis_ ) {
				max_overlap_dis_ = bump_grid()->required_separation_distance( ProbeRadius( ii ), ProbeRadius( jj ) );
			}
		}
	}

}

bool UpstreamDownstreamCollisionFilter::passes_etable_filter( match_dspos1 const & m ) const
{
	using namespace core;
	using namespace core::conformation;
	using namespace core::pose;

	runtime_assert( dsbuilders_[ m.originating_geom_cst_for_dspos ] != 0 );
	dsbuilders_[ m.originating_geom_cst_for_dspos ]->coordinates_from_hit(
		full_hit( m ), downstream_atoms_, coords_ );
	for ( Size ii = 1; ii <= downstream_atoms_.size(); ++ii ) {
		downstream_pose_->set_xyz( downstream_atoms_[ ii ], coords_[ ii ] );
	}

	using namespace core::scoring;
	EnergyMap emap;
	for ( Size ii = 1; ii < m.upstream_hits.size(); ++ii ) {
		if ( ii == m.originating_geom_cst_for_dspos ) continue; // don't collision check since we've presumably done so already
		if ( us_ds_chemical_bond_[ ii ] ) continue;
		for ( Size jj = 1; jj <= downstream_pose_->total_residue(); ++jj ) {
			emap[ fa_atr ] = 0; emap[ fa_rep ] = 0; emap[ fa_sol ] = 0;
			etable_energy()->residue_pair_energy(
				*( cacher_->upstream_conformation_for_hit( ii, fake_hit( m.upstream_hits[ ii ] )) ),
				downstream_pose_->residue( jj ),
					*empty_pose(), *empty_sfxn(),
				emap );
			Real energy = wfa_atr() * emap[ fa_atr ] + wfa_rep() * emap[ fa_rep ] + wfa_sol() * emap[ fa_sol ];
			if ( energy > lj_cutoff() ) return false;
		}
	}
	return true;

}

bool UpstreamDownstreamCollisionFilter::passes_hardsphere_filter( match_dspos1 const & m ) const
{
	runtime_assert( dsbuilders_[ m.originating_geom_cst_for_dspos ] != 0 );
	dsbuilders_[ m.originating_geom_cst_for_dspos ]->coordinates_from_hit(
		full_hit( m ), downstream_atoms_, coords_ );

	for ( Size ii = 1; ii <= m.upstream_hits.size(); ++ii ) {
		if ( ii == m.originating_geom_cst_for_dspos )	continue;
		if ( us_ds_chemical_bond_[ ii ] ) continue;

		if( !passes_hardsphere_filter( ii, fake_hit( m.upstream_hits[ ii ] ), coords_ ) ) return false;
	} //loop over all upstream hits
	return true;
}

bool
UpstreamDownstreamCollisionFilter::passes_hardsphere_filter(
	core::Size geomcst_up,
	core::Size geomcst_down,
	Hit const & upstream_hit,
	Hit const & downstream_hit
) const
{
	dsbuilders_[ geomcst_down ]->coordinates_from_hit( downstream_hit, downstream_atoms_, coords_ );
	return passes_hardsphere_filter( geomcst_up, upstream_hit, coords_ );
}

}
}
}
