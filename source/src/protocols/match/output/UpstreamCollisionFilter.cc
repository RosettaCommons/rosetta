// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/output/UpstreamCollisionFilter.hh
/// @brief  Implementation for class to filter matches where the upstream residues collide.
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/output/UpstreamCollisionFilter.hh>

// Package headers
#include <protocols/match/BumpGrid.hh>
#include <protocols/match/output/UpstreamHitCacher.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include	<core/pose/Pose.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/scoring/etable/EtableEnergy.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Utility headers
#include <utility/pointer/ReferenceCount.hh>

// C++ headers
#include <iostream>

#include <protocols/match/Hit.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace output {


MatchCollisionFilter::MatchCollisionFilter(
  std::string filter_name,
	UpstreamHitCacherOP coordinate_cacher
) :
	MatchFilter( filter_name ),
	cacher_( coordinate_cacher ),
	filter_by_lj_( false ),
	wfa_atr_( 0.8 ),
	wfa_rep_( 0.44 ),
	wfa_sol_( 0.6 ),
	lj_cutoff_( 10 ),
	tolerated_overlap_( 0.0 ),
	empty_pose_( core::pose::PoseOP( new core::pose::Pose ) ),
	empty_sfxn_( core::scoring::ScoreFunctionOP( new core::scoring::ScoreFunction ) ),
	etable_energy_( /* 0 */ ),
	bump_grid_( BumpGridOP( new BumpGrid ) )
{}

MatchCollisionFilter::~MatchCollisionFilter(){}

void MatchCollisionFilter::set_filter_by_lj( bool setting )
{
	filter_by_lj_ = setting;
	if ( filter_by_lj_ ) {
		using namespace core::scoring;
		using namespace core::scoring::etable;
		using namespace core::scoring::methods;
		EnergyMethodOptions eopts;
		etable_energy_ = core::scoring::methods::ShortRangeTwoBodyEnergyOP( new TableLookupEtableEnergy(
			*ScoringManager::get_instance()->etable( eopts ).lock(), eopts, false /*do_classic_intrares*/ ) ); // FIXME: passing reference of a temporairly locked AP
	}
}

void MatchCollisionFilter::set_lj_cutoff( Real setting )
{
	lj_cutoff_ = setting;
}

void MatchCollisionFilter::set_lj_atr_weight( Real setting )
{
	wfa_atr_ = setting;
}

void MatchCollisionFilter::set_lj_rep_weight( Real setting )
{
	wfa_rep_ = setting;
}

void MatchCollisionFilter::set_lj_sol_weight( Real setting )
{
	wfa_sol_ = setting;
}

void MatchCollisionFilter::set_tolerated_overlap( Real setting )
{
	tolerated_overlap_ = setting;
	bump_grid_->set_general_overlap_tolerance( tolerated_overlap_ );
}

bool
MatchCollisionFilter::filter_by_lj() const {
	return filter_by_lj_;}

MatchCollisionFilter::Real
MatchCollisionFilter::wfa_atr() const {
	return wfa_atr_;}

MatchCollisionFilter::Real
MatchCollisionFilter::wfa_rep() const {
	return wfa_rep_;}

MatchCollisionFilter::Real
MatchCollisionFilter::wfa_sol() const {
	return wfa_sol_;}

MatchCollisionFilter::Real
MatchCollisionFilter::lj_cutoff() const {
	return lj_cutoff_;}

MatchCollisionFilter::Real
MatchCollisionFilter::tolerated_overlap() const {
	return tolerated_overlap_;}

core::pose::PoseCOP
MatchCollisionFilter::empty_pose() const {
	return empty_pose_;}

core::scoring::ScoreFunctionOP
MatchCollisionFilter::empty_sfxn() const{
	return empty_sfxn_;}

core::scoring::methods::ShortRangeTwoBodyEnergyCOP
MatchCollisionFilter::etable_energy() const{
	return etable_energy_; }

BumpGridCOP
MatchCollisionFilter::bump_grid() const{
	return bump_grid_; }


UpstreamCollisionFilter::UpstreamCollisionFilter(
  std::string filter_name,
	UpstreamHitCacherOP coordinate_cacher
) :
	parent( filter_name, coordinate_cacher )
{
	std::cout << "Created UpstreamCollisionFilter" << std::endl;
}

UpstreamCollisionFilter::~UpstreamCollisionFilter()
{}

/// @brief Either sphere-overlap checks the atom pairs in the match residues, or
/// evaluates the Etable energies.  Returns false if any atom pair collides more than
/// the cutoff threshold or if the residue pair energy exceeds the energy cutoff.
/// Returns true otherwise.
bool
UpstreamCollisionFilter::passes_filter(
	match const & m
) const
{
	return passes_filter( match_dspos1( m, 1 ) );
}

bool
UpstreamCollisionFilter::passes_filter(
	match_dspos1 const & m
) const
{
	if ( filter_by_lj() ) {
		using namespace core::scoring;
		EnergyMap emap;
		for ( Size ii = 1; ii < m.upstream_hits.size(); ++ii ) {
			for ( Size jj = ii + 1; jj <= m.upstream_hits.size(); ++jj ) {
				emap[ fa_atr ] = 0; emap[ fa_rep ] = 0; emap[ fa_sol ] = 0;
				etable_energy()->residue_pair_energy(
					*( cacher_->upstream_conformation_for_hit( ii, fake_hit( m.upstream_hits[ ii ] )) ),
					*( cacher_->upstream_conformation_for_hit( jj, fake_hit( m.upstream_hits[ jj ] )) ),
						*empty_pose(), *empty_sfxn(),
					emap );
				Real energy = wfa_atr() * emap[ fa_atr ] + wfa_rep() * emap[ fa_rep ] + wfa_sol() * emap[ fa_sol ];
				if ( energy > lj_cutoff() ) return false;
			}
		}
		return true;
	} else {
		for ( Size ii = 1; ii < m.upstream_hits.size(); ++ii ) {
			core::conformation::ResidueCOP iires = cacher_->upstream_conformation_for_hit( ii, fake_hit( m.upstream_hits[ ii ] ) );
			for ( Size jj = ii + 1; jj <= m.upstream_hits.size(); ++jj ) {
				core::conformation::ResidueCOP jjres = cacher_->upstream_conformation_for_hit( jj, fake_hit( m.upstream_hits[ jj ] ) );
				if( !passes_hardsphere_filter( *iires, *jjres ) ) return false;
			} // jj loop over upstream_hits.size()
		} //ii loop over upstream_hits.size()
		return true;
	}
}

bool
UpstreamCollisionFilter::passes_hardsphere_filter(
	core::Size geomcst_a,
  core::Size geomcst_b,
	Hit const & upstream_hit_a,
	Hit const & upstream_hit_b
) const
{
	core::conformation::ResidueCOP resA = cacher_->upstream_conformation_for_hit( geomcst_a, upstream_hit_a );
	core::conformation::ResidueCOP resB = cacher_->upstream_conformation_for_hit( geomcst_b, upstream_hit_b );
	return passes_hardsphere_filter( *resA, *resB );
}


}
}
}
