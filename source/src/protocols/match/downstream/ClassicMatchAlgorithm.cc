// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/match/downstream/ClassicMatchAlgorithm.cc
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/downstream/ClassicMatchAlgorithm.hh>

// Package headers
#include <protocols/match/BumpGrid.hh>
#include <protocols/match/Matcher.hh>
#include <protocols/match/OccupiedSpaceHash.hh>
#include <protocols/match/downstream/ActiveSiteGrid.hh>
#include <protocols/match/upstream/UpstreamBuilder.hh>
#include <protocols/match/downstream/DownstreamBuilder.hh>

// Project headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <basic/Tracer.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <list>

#include <protocols/match/Hit.hh>
#include <utility/vector1.hh>
#include <utility/options/StringVectorOption.hh>


namespace protocols {
namespace match {
namespace downstream {

static basic::Tracer TR( "protocols.match.downstream.ClassicMatchAlgorithm" );

ClassicMatchAlgorithm::ClassicMatchAlgorithm( Size geom_cst_id ) :
	parent( geom_cst_id ),
	occspace_rev_id_at_last_update_( 0 ),
	build_round1_hits_twice_( false ),
	completed_first_round1_hit_building_( false )
{}

ClassicMatchAlgorithm::~ClassicMatchAlgorithm() {}

DownstreamAlgorithmOP
ClassicMatchAlgorithm::clone() const {
	return DownstreamAlgorithmOP( new ClassicMatchAlgorithm( *this ) );
}

void
ClassicMatchAlgorithm::set_build_round1_hits_twice()
{
	debug_assert( geom_cst_id() == 1 );
	build_round1_hits_twice_ = true;
}

std::list< Hit >
ClassicMatchAlgorithm::build_hits_at_all_positions(
	Matcher & matcher
)
{
	/// First, make sure that all downstream builders for this geometric constraint
	/// are pointing at the matcher's occupied space grid.

	if ( geom_cst_id() == 1 && build_round1_hits_twice_ && ! completed_first_round1_hit_building_ ) {
		/// set the occupied space grid pointer to NULL: do not discard hits because they do not fall in
		/// the same volume of space as some other hit.  The OccSpaceHash itself contains logic to avoid
		/// discarding such hits only so long as it is completely empty.  However, that logic is not good
		/// enough if the strategy is to build round1 hits twice.
		std::list< DownstreamBuilderOP > const & dsbuilders(
			matcher.downstream_builders( geom_cst_id() ));
		for ( std::list< DownstreamBuilderOP >::const_iterator
				iter = dsbuilders.begin(),
				iter_end = dsbuilders.end();
				iter != iter_end; ++iter ) {
			(*iter)->set_occupied_space_hash( 0 );
		}
	} else {
		std::list< DownstreamBuilderOP > const & dsbuilders(
			matcher.downstream_builders( geom_cst_id() ));
		for ( std::list< DownstreamBuilderOP >::const_iterator
				iter = dsbuilders.begin(),
				iter_end = dsbuilders.end();
				iter != iter_end; ++iter ) {
			(*iter)->set_occupied_space_hash( matcher.occ_space_hash() );
		}
	}

	/// Later, we'll have code right here that looks at the occspace hash and figures out which
	/// regions of 3D each sidechain atom of the upstream residues would have to be in to produce
	/// a viable hit, and then we'll adjust the usptream builder for this residue to filter out
	/// atoms early on in hit generation... so exciting.  For now, we do not.

	if ( geom_cst_id() == 1 && build_round1_hits_twice_ && ! completed_first_round1_hit_building_ ) {
		std::list< Hit > hits = build_and_discard_first_round_hits_at_all_positions( matcher );
		completed_first_round1_hit_building_ = true;
		return hits;
	} else {
		return parent::default_build_hits_at_all_positions( matcher );
	}
}

std::list< Hit >
ClassicMatchAlgorithm::build_and_discard_first_round_hits_at_all_positions(
	Matcher & matcher
)
{
	debug_assert( geom_cst_id() == 1 );

	utility::vector1< upstream::ScaffoldBuildPointCOP > const & launch_points
		( matcher.per_constraint_build_points( geom_cst_id() ) );
	Size n_build_points = launch_points.size();

	std::list< Hit > return_hits; // Only return a single hit from this function
	OccupiedSpaceHashOP occspace = matcher.occ_space_hash();
	for ( Size ii = 1; ii <= n_build_points; ++ii ) {
		// generate hits for build point ii, and insert them into the occspace hash,
		// but throw them out at the end of this iteration since usually there are too
		// many that get generated
		std::list< Hit > iihits = matcher.upstream_builder( geom_cst_id() )->build( * launch_points[ ii ] );
		for ( std::list< Hit >::const_iterator iter = iihits.begin(), iter_end = iihits.end();
				iter != iter_end; ++iter ) {
			occspace->insert_hit_geometry( iter->second() );
		}
		// save one hit so that the Matcher doesn't exit early (it will if we return 0 hits)
		if ( return_hits.empty() && ! iihits.empty() ) { return_hits.push_back( *iihits.begin() ); }
	}
	return return_hits;
}


/// @details Reset the occupied space hash that the matcher uses so that
/// it reflects the hits generated this round; this will cause the invalidation
/// of hits generated in previous rounds.  These invalidated hits will be deleted
/// in calls to respond_to_peripheral_hitlist_change.
void
ClassicMatchAlgorithm::respond_to_primary_hitlist_change(
	Matcher & matcher,
	Size round_just_completed
)
{
	OccupiedSpaceHashOP occspace = matcher.occ_space_hash();

	if ( geom_cst_id() == round_just_completed ) {
		TR << "Finished round " << geom_cst_id() << " with " << matcher.hits( geom_cst_id() ).size();
		TR << " hits." << std::endl;
	}

	if ( geom_cst_id() == 1 && round_just_completed == 1 ) {

		/// the occspace hash has already been updated if we're using the build_round1_hits_twice logic
		if ( build_round1_hits_twice_ ) return;

		/// The first geometric constraint inserts hits into the occupied space grid;
		/// the later geometric constraints merely mark voxels already present in the
		/// occupied space grid with 1.
		for ( Matcher::HitListConstIterator
				iter = matcher.hits( geom_cst_id() ).begin(),
				iter_end = matcher.hits( geom_cst_id() ).end();
				iter != iter_end; ++iter ) {
			occspace->insert_hit_geometry( iter->second() );
		}

	} else {

		if ( geom_cst_id() != round_just_completed ) {
			TR << "Updating the occupied space grid with " << matcher.hits( geom_cst_id() ).size();
			TR << " hits from geometric constraint # " << geom_cst_id() << std::endl;
		}

		occspace->prepare_to_note_hits_for_completed_round();

		/// Prepare to clear the cobwebs.  Note which voxels in the occ_space_hash_
		/// could lead to matches, and which voxels could not possibly lead to matches.
		for ( Matcher::HitListConstIterator
				iter = matcher.hits( geom_cst_id() ).begin(),
				iter_end = matcher.hits( geom_cst_id() ).end();
				iter != iter_end; ++iter ) {
			occspace->note_hit_geometry( iter->second() );
		}

		occspace->drop_unsatisfied_voxels();
	}

	occspace_rev_id_at_last_update_ = occspace->revision_id();
}

/// @details Drop hits that had previously seemed viable after another round completed;
/// during that round, certain previously occupied regions of 6D were not filled
/// with new hits.  Any previously-generated hit that falls into a region of 6D which is
/// no longer occupied should be elminated since it could not ever result in a match;
/// it is inviable.
void
ClassicMatchAlgorithm::respond_to_peripheral_hitlist_change( Matcher & matcher )
{
	OccupiedSpaceHashOP occspace = matcher.occ_space_hash();
	if ( occspace_rev_id_at_last_update_ == occspace->revision_id() ) {
		/// Nothing about the occupied space hash has changed since I last pruned
		/// my non-viable hits.  There are no more non-viable hits that need pruning.
		return;
	}

	Matcher::HitListIterator iter = matcher.hit_list_begin(   geom_cst_id() );
	Matcher::HitListIterator iter_end = matcher.hit_list_end( geom_cst_id() );

	Size drop_count( 0 );

	while ( iter != iter_end  ) {
		Matcher::HitListIterator iter_next = iter;
		++iter_next;
		if ( ! occspace->previous_round_geometry_still_matchable( iter->second() ) ) {
			matcher.erase_hit( *this, geom_cst_id(), iter );
			++drop_count;
		}
		iter = iter_next;
	}

	TR << "Erased " << drop_count << " round " << geom_cst_id();
	TR << " hits with ";
	TR << matcher.hits( geom_cst_id() ).size() << " hits remaining." << std::endl;

	occspace_rev_id_at_last_update_ = occspace->revision_id();
}


std::list< Hit >
ClassicMatchAlgorithm::build(
	Size const scaffold_build_point_id,
	Size const upstream_conf_id,
	core::conformation::Residue const & upstream_residue
) const
{
	std::list< Hit > local_hit_list;

	/// Sanity check: the input residue must match the residue type this MatchAlgorithm is expecting.
	if ( & upstream_residue.type() != restype_.get() ) {
		utility_exit_with_message( "ERROR: Classic Match Algorithm was expecting restype: " + restype_->name() + " but was given " + upstream_residue.name() + ".  Cannot continue" );
	}

	/// Sanity check: all transforms must have been precomputed before reaching this stage.
	for ( auto const & exsampler_vec : external_samplers_ ) {
		for ( auto const & exsampler : exsampler_vec ) {
			if ( ! exsampler.transforms_uptodate() ) {
				utility_exit_with_message( "CRITICAL ERROR in ClassicMatchAlgorithm::build.  ExternalGeomSampler transforms are not up-to-date" );
			}
		}
	}

	for ( Size ii = 1; ii <= n_external_samplers(); ++ii ) {
		std::list< Hit > hits = build_from_three_coords(
			ii,
			scaffold_build_point_id,
			upstream_conf_id,
			upstream_residue );
		local_hit_list.splice( local_hit_list.end(), hits );
	}

	return local_hit_list;
}

bool
ClassicMatchAlgorithm::upstream_only() const
{
	return false;
}

bool
ClassicMatchAlgorithm::generates_primary_hits() const
{
	return true;
}


/// @details If this function is causing an exit, then there is a bug within the Matcher's
/// match-enumeration logic.  There is no meaningful way forward after this function is invoked.
/// It should not be invoked.
HitPtrListCOP
ClassicMatchAlgorithm::hits_to_include_with_partial_match( match_dspos1 const & ) const
{
	HitPtrListCOP empty;
	utility_exit_with_message( "Cannot invoke ClassicMatchAlgorithm::hits_to_include_with_partial_match()" );
	return empty;
}


ClassicMatchAlgorithm::Size
ClassicMatchAlgorithm::n_possible_hits_per_upstream_conformation() const
{
	Size total = 0;
	for ( Size ii = 1; ii <= n_external_samplers(); ++ii ) {
		Size iitotal = 1;
		utility::vector1< toolbox::match_enzdes_util::ExternalGeomSampler > const & exsampler_list( external_samplers_[ ii ] );

		// for ( Size jj = 1; jj <= exsampler.size(); ++jj ) {
		for ( auto const & exsampler : exsampler_list ) {
			iitotal *= exsampler.n_tor_U3D1_samples();
			iitotal *= exsampler.n_ang_U2D1_samples();
			iitotal *= exsampler.n_dis_U1D1_samples();
			iitotal *= exsampler.n_tor_U2D2_samples();
			iitotal *= exsampler.n_ang_U1D2_samples();
			iitotal *= exsampler.n_tor_U1D3_samples();

			iitotal *= dsbuilders_[ ii ]->n_possible_hits_per_at3frame();

			total += iitotal;
		}
	}
	return total;
}


std::list< Hit >
ClassicMatchAlgorithm::build_from_three_coords(
	Size const which_external_sampler,
	Size const scaffold_build_point_id,
	Size const upstream_conf_id,
	core::conformation::Residue const & upstream_residue
) const
{
	using namespace toolbox::match_enzdes_util;
	Vector coord1( upstream_residue.xyz( launch_atom( which_external_sampler, 1 )));
	Vector coord2( upstream_residue.xyz( launch_atom( which_external_sampler, 2 )));
	Vector coord3( upstream_residue.xyz( launch_atom( which_external_sampler, 3 )));


	// this returns a vector of ExternalGeomSamplers
	utility::vector1< toolbox::match_enzdes_util::ExternalGeomSampler > const & exsampler_list( external_samplers_[ which_external_sampler ] );

	DownstreamBuilderCOP dsbuilder( dsbuilders_[ which_external_sampler ] );

	std::list< Hit > hitlist;

	ProbeRadius const
		radD1( dsbuilder->atom1_radius() ),
		radD2( dsbuilder->atom2_radius()),
		radD3( dsbuilder->atom3_radius() );

	bool active_site_check_D1( active_site_grid_set() && dsbuilder->atom1_belongs_in_active_site() );
	bool active_site_check_D2( active_site_grid_set() && dsbuilder->atom2_belongs_in_active_site() );
	bool active_site_check_D3( active_site_grid_set() && dsbuilder->atom3_belongs_in_active_site() );


	/// 1. Define a coordinate frame from the three coordinates of atoms 1 2 and 3 in the
	/// canonical form such that the point is at coord1, z lies along the coord1-coord2 vector,
	/// y lies in the plane of coord1, 2, & 3, and x is the cross product of y and z.
	HTReal ht_start( coord3, coord2, coord1 );


	// This 7th loop (the outer-most loop) iterates over a list of ExternalGeomSamplers associated with the same constraints
	for ( auto const & exsampler : exsampler_list ) {
		/// Six nested for loops to iterate across all the specified geometries that describe how to
		/// build the downstream target from the upstream conformation.  No transcendental function evaulations
		/// necessary here, as the coordinate transformations have already been computed!  Just good old
		/// addition and multiplication.
		for ( Size ii = 1; ii <= exsampler.n_tor_U3D1_samples(); ++ii ) {
			HTReal ht_ii = ht_start * exsampler.transform( HT_tor_U3D1, ii );

			for ( Size jj = 1; jj <= exsampler.n_ang_U2D1_samples(); ++jj ) {
				HTReal ht_jj = ht_ii * exsampler.transform( HT_ang_U2D1, jj );

				for ( Size kk = 1; kk <= exsampler.n_dis_U1D1_samples(); ++kk ) {
					HTReal ht_kk = ht_jj;
					ht_kk.walk_along_z( exsampler.dis_U1D1_samples()[ kk ] );
					Vector pD1 = ht_kk.point();
					if ( radD1 > ZERO && bbgrid().occupied( radD1, pD1 ) ) continue; /// Collision check atom D1
					if ( active_site_check_D1 && ! active_site_grid().occupied( pD1 ) ) continue;

					for ( Size ll = 1; ll <= exsampler.n_tor_U2D2_samples(); ++ll ) {
						HTReal ht_ll = ht_kk * exsampler.transform( HT_tor_U2D2, ll );

						for ( Size mm = 1; mm <= exsampler.n_ang_U1D2_samples(); ++mm ) {
							HTReal ht_mm = ht_ll * exsampler.transform( HT_ang_U1D2, mm );
							Vector pD2 = ht_mm.point();
							if ( radD2 > ZERO && bbgrid().occupied( radD2, pD2 ) ) continue; /// Collision check atom D2
							if ( active_site_check_D2 && ! active_site_grid().occupied( pD2 ) ) continue;

							for ( Size nn = 1; nn <= exsampler.n_tor_U1D3_samples(); ++nn ) {
								HTReal ht_nn = ht_mm * exsampler.transform( HT_tor_U1D3, nn );
								Vector pD3 = ht_nn.point();
								if ( radD3 > ZERO && bbgrid().occupied( radD3, pD3 ) ) continue; /// Collision check atom D3
								if ( active_site_check_D3 && ! active_site_grid().occupied( pD3 ) ) continue;

								std::list< Hit > nn_hits = dsbuilder->build(
									ht_nn,
									scaffold_build_point_id,
									upstream_conf_id,
									exgeom_ids_[ which_external_sampler ],
									upstream_residue );
								hitlist.splice( hitlist.end(), nn_hits );
							}
						}
					}
				}
			}
		}
	}
	return hitlist;

}

////// Debugging code below -- insert this after atom 6's collision check.
/*
std::cout << "Collision free placement of atoms 4 5 and 6: ";
std::cout << "p4: ";
for ( Size oo = 1; oo <= 3; ++oo ) std::cout << p4( oo ) << " ";
std::cout << "p5: ";
for ( Size oo = 1; oo <= 3; ++oo ) std::cout << p5( oo ) << " ";
std::cout << "p6: ";
for ( Size oo = 1; oo <= 3; ++oo ) std::cout << p6( oo ) << " ";
std::cout << std::endl;

std::cout << "tor_U3D1: expected: " << exsampler.tor_U3D1_samples()[ ii ] << " real: "
<< numeric::constants::d::radians_to_degrees * numeric::dihedral_radians(
coord1, coord2, coord3, p4 ) << std::endl;

std::cout << "ang_U2D1: expected: " << exsampler.ang_U2D1_samples()[ jj ] << " real: "
<< numeric::constants::d::radians_to_degrees * numeric::angle_radians(
coord2, coord3, p4 ) << std::endl;

std::cout << "dis_U1D1: expected: " << exsampler.dis_U1D1_samples()[ kk ] << " real: "
<< p4.distance( coord3 ) << std::endl;

std::cout << "tor_U2D2: expected: " << exsampler.tor_U2D2_samples()[ ll ] << " real: "
<< numeric::constants::d::radians_to_degrees * numeric::dihedral_radians(
coord2, coord3, p4, p5 ) << std::endl;

std::cout << "ang_U1D2: expected: " << exsampler.ang_U1D2_samples()[ mm ] << " real: "
<< numeric::constants::d::radians_to_degrees * numeric::angle_radians(
coord3, p4, p5 ) << std::endl;

std::cout << "tor_U1D3: expected: " << exsampler.tor_U1D3_samples()[ nn ] << " real: "
<< numeric::constants::d::radians_to_degrees * numeric::dihedral_radians(
coord3, p4, p5, p6 ) << std::endl;

std::cout << "EXGEOM SAMPLE tor_U3D1: " << exsampler.tor_U3D1_samples()[ ii ];
std::cout << " ang_U2D1: " << exsampler.ang_U2D1_samples()[ jj ];
std::cout << " dis_U1D1: " << exsampler.dis_U1D1_samples()[ kk ];
std::cout << " tor_U2D2: " << exsampler.tor_U2D2_samples()[ ll ];
std::cout << " ang_U1D2: " << exsampler.ang_U1D2_samples()[ mm ];
std::cout << " tor_U1D3: " << exsampler.tor_U1D3_samples()[ nn ] << std::endl;


*/


void ClassicMatchAlgorithm::set_residue_type( core::chemical::ResidueTypeCOP restype )
{
	restype_ = restype;
	external_samplers_.clear();
	launch_points_.clear();
}


/// @details Precompute transforms for the external geom sampler as it is added
/// so that the transforms are ready when build() is called.
void ClassicMatchAlgorithm::add_external_geom_sampler(
	utility::vector1< toolbox::match_enzdes_util::ExternalGeomSampler > const & sampler,
	Size const exgeom_id,
	std::string const & atom1,
	std::string const & atom2,
	std::string const & atom3,
	DownstreamBuilderCOP dsbuilder
)
{

	Size id1( restype().has( atom1 ) ? restype().atom_index( atom1 ) : 0 );
	Size id2( restype().has( atom2 ) ? restype().atom_index( atom2 ) : 0 );
	Size id3( restype().has( atom3 ) ? restype().atom_index( atom3 ) : 0 );

	if ( id1 == 0 ) {
		utility_exit_with_message( "ERROR in adding external geom sampler to ClassicMatchAlgorithm: " + restype().name() + " does not contain requested atom " + atom1 );
	}
	if ( id2 == 0 ) {
		utility_exit_with_message( "ERROR in adding external geom sampler to ClassicMatchAlgorithm: " + restype().name() + " does not contain requested atom " + atom2 );
	}
	if ( id3 == 0 ) {
		utility_exit_with_message( "ERROR in adding external geom sampler to ClassicMatchAlgorithm: " + restype().name() + " does not contain requested atom " + atom3 );
	}

	utility::fixedsizearray1< Size, 3 > atids;
	atids[ 1 ] = id1;
	atids[ 2 ] = id2;
	atids[ 3 ] = id3;

	external_samplers_.push_back( sampler );
	launch_points_.push_back( atids );
	dsbuilders_.push_back( dsbuilder );
	exgeom_ids_.push_back( exgeom_id );

	/// initialize the external sampler transforms
	utility::vector1< toolbox::match_enzdes_util::ExternalGeomSampler > & samp_list( external_samplers_[ external_samplers_.size() ] );


	Real atom1_atom2_distance = dsbuilder->atom1_atom2_distance();
	Real atom2_atom3_distance = dsbuilder->atom2_atom3_distance();
	Real atom1_atom2_atom3_angle = dsbuilder->atom1_atom2_atom3_angle();

	for ( auto & samp : samp_list ) {
		samp.set_dis_D1D2( atom1_atom2_distance );
		samp.set_dis_D2D3( atom2_atom3_distance );
		samp.set_ang_D1D2D3( atom1_atom2_atom3_angle );
		samp.precompute_transforms();
	}
}

void
ClassicMatchAlgorithm::clear_external_geom_samplers()
{
	external_samplers_.clear();
	launch_points_.clear();
	dsbuilders_.clear();
}

}
}
}
