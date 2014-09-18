// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/downstream/DownstreamAlgorithm.cc
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/downstream/DownstreamAlgorithm.hh>

// Package headers
#include <protocols/match/BumpGrid.hh>
#include <protocols/match/Matcher.hh>
#include <protocols/match/downstream/ActiveSiteGrid.hh>
#include <protocols/match/upstream/UpstreamBuilder.hh>
// AUTO-REMOVED #include <protocols/match/upstream/ScaffoldBuildPoint.hh>
#include <protocols/match/downstream/DownstreamBuilder.hh>
// AUTO-REMOVED #include <protocols/match/downstream/LigandConformer.hh>

//Project header
#include <core/conformation/Residue.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <list>

#include <core/id/AtomID.hh>
#include <protocols/match/Hit.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace downstream {

static thread_local basic::Tracer TR( "protocols.match.downstream.DownstreamAlgorithm" );

DownstreamAlgorithm::DownstreamAlgorithm( Size geom_cst_id ) : geom_cst_id_( geom_cst_id ) {}
DownstreamAlgorithm::DownstreamAlgorithm( DownstreamAlgorithm const & other ) :
	utility::pointer::ReferenceCount(),
	geom_cst_id_( other.geom_cst_id_ ),
	bbgrid_( other.bbgrid_ ),
	active_site_grid_( other.active_site_grid_ )
{}

DownstreamAlgorithm const &
DownstreamAlgorithm::operator = ( DownstreamAlgorithm const & rhs )
{
	if ( this != & rhs ) {
		geom_cst_id_ = rhs.geom_cst_id_;
		bbgrid_ = rhs.bbgrid_;
		active_site_grid_ = rhs.active_site_grid_;
	}
	return *this;
}


DownstreamAlgorithm::~DownstreamAlgorithm() {}

/// @details By initializing local std::list< Hit > variables inside the loop
/// over all of the build points, and then splicing them into a central vector
/// of hit lists, I can avoid expensive list-copy operations while guaranteeing
/// OpenMP thread saftey.
std::list< Hit >
DownstreamAlgorithm::build_hits_at_all_positions(
	Matcher & matcher
)
{
	return default_build_hits_at_all_positions( matcher );
}

/// @details no-op
void
DownstreamAlgorithm::prepare_for_match_enumeration( Matcher const & ) {}



/// @details Noop in base class.
void
DownstreamAlgorithm::respond_to_primary_hitlist_change( Matcher &, Size )
{}

/// @details Noop in base class.
void
DownstreamAlgorithm::respond_to_peripheral_hitlist_change( Matcher & )
{}

DownstreamAlgorithm::Size
DownstreamAlgorithm::geom_cst_id() const {
	return geom_cst_id_;
}

void
DownstreamAlgorithm::set_bb_grid(
	BumpGridCOP bbgrid
)
{
	bbgrid_ = bbgrid;
}

void
DownstreamAlgorithm::set_active_site_grid(
	ActiveSiteGridCOP active_site_grid
)
{
	active_site_grid_ = active_site_grid;
}

void
DownstreamAlgorithm::set_dsbuilder(
	DownstreamBuilderOP dsbuilder )
{
  if (dsbuilder_) {
		dsbuilder_.reset_to_null();
  }
	dsbuilder_ = dsbuilder;
}

DownstreamBuilderOP
DownstreamAlgorithm::get_dsbuilder() const
{
	return dsbuilder_;
}

bool
DownstreamAlgorithm::are_colliding(
	core::conformation::Residue const & us_res /*upstream*/,
	core::conformation::Residue const & ds_res /*downstream*/,
	utility::vector1< core::id::AtomID > const & ds_atoms,
	utility::vector1< core::Size > const & catalytic_atoms
) const
{
	//TR << "RES name" << us_res.name() << "SIZE " << ds_atoms.size() << std::endl;
	//TR << "CAT 1: " << us_res.atom_name( catalytic_atoms[1] ) <<
	//      "CAT 2: " << us_res.atom_name( catalytic_atoms[2] ) <<
	//      "CAT 3: " << ds_res.atom_name( catalytic_atoms[3] ) <<
	//      "CAT 4: " << ds_res.atom_name( catalytic_atoms[4] ) << std::endl;

  for ( Size atomid_ds = 1; atomid_ds <= ds_atoms.size(); ++atomid_ds ) {
     ProbeRadius ds_rad = probe_radius_for_atom_type( ds_res.atom_type_index( ds_atoms[ atomid_ds ].atomno() ) );
		 if ( ! ( catalytic_atoms[3] == ds_atoms[ atomid_ds ].atomno() ||
							catalytic_atoms[4] ==	ds_atoms[ atomid_ds ].atomno() ) ){
  		 //TR << " x "<< ds_res.xyz( ds_atoms[ atomid_ds ].atomno() )[0] << " y " << ds_res.xyz( ds_atoms[ atomid_ds ].atomno()  )[1] << " z " << ds_res.xyz( ds_atoms[ atomid_ds ].atomno()  )[2] << std::endl;

    		for (Size atomid_us = us_res.first_sidechain_atom(); atomid_us <= us_res.nheavyatoms(); ++atomid_us ) {
		    	if ( ! ( catalytic_atoms[1] == atomid_us ||
              		 catalytic_atoms[2] == atomid_us ) ){

						//TR << "RES ds: " << ds_res.atom_name( ds_atoms[ atomid_ds ].atomno() ) << " RES us: " << us_res.atom_name(atomid_us) << std::endl;
      			core::Real dist2( ds_res.xyz( ds_atoms[ atomid_ds ].atomno() ).distance_squared( us_res.xyz( atomid_us ) ) );
      			ProbeRadius us_rad = probe_radius_for_atom_type( us_res.atom_type_index( atomid_us ) );
						core::Real dis = bbgrid().required_separation_distance( us_rad, ds_rad );
						dis = dis*dis;
						//TR << "griddis "  << dis << " dist2 " << dist2 << std::endl;
      		if ( dist2 < dis ) {
						//TR << "Colliding" << std::endl;
  					return true;
      		}
				}
			}
		}
  }
	return false;
}

std::list< Hit >
DownstreamAlgorithm::default_build_hits_at_all_positions(
	Matcher const & matcher
) const
{
	utility::vector1< upstream::ScaffoldBuildPointCOP > const & launch_points
		( matcher.per_constraint_build_points( geom_cst_id_ ) );
	Size n_build_points = launch_points.size();

	std::list< Hit > all_hits;
	utility::vector1< std::list< Hit > > hits( n_build_points );

	/// Generate conformations for the upstream and downstream partners for each of the
	/// possible scaffold build points for this geometric constraint.
	/// This loop will be parallelized.  Everything down stream of this call is const,
	/// in spite of the fact that the matcher is handed as a non-const reference.
#ifdef USE_OPENMP
	#pragma omp parallel for schedule(dynamic,1)
#endif
	for ( Size ii = 1; ii <= n_build_points; ++ii ) {
		std::list< Hit > iihits = matcher.upstream_builder( geom_cst_id_ )->build( * launch_points[ ii ] );
		hits[ ii ].splice( hits[ ii ].end(), iihits );
	}

	for ( Size ii = 1; ii <= n_build_points; ++ii ) {
		all_hits.splice( all_hits.end(), hits[ ii ] );
	}

	return all_hits;
}


}
}
}
