// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/match/output/DownstreamRMSEvaluator.cc
/// @brief
/// @author Alex Zanghellini (zanghell@u.washington.edu)
/// @author Andrew Leaver-Fay (aleaverfay@gmail.com), porting to mini

// Unit headers
#include <protocols/match/output/DownstreamRMSEvaluator.hh>

// Package headers
#include <protocols/match/output/MatchEvaluator.hh>
#include <protocols/match/downstream/DownstreamBuilder.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <core/id/AtomID.hh>
#include <protocols/match/Hit.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace match {
namespace output {


DownstreamRMSEvaluator::DownstreamRMSEvaluator() {}

DownstreamRMSEvaluator::~DownstreamRMSEvaluator() {}

void
DownstreamRMSEvaluator::set_downstream_pose( core::pose::PoseCOP dspose )
{
	dspose_ = dspose;

	/// TEMP!
	/// Compare all heavy atoms on all residues
	Size count_atoms = 0;
	for ( Size ii = 1; ii <= dspose_->total_residue(); ++ii ) {
		count_atoms += dspose_->residue(ii).nheavyatoms();
	}

	atoms_to_compare_.resize( count_atoms );
	count_atoms = 1;
	for ( Size ii = 1; ii <= dspose_->total_residue(); ++ii ) {
		for ( Size jj = 1; jj <= dspose_->residue(ii).nheavyatoms(); ++jj ) {
			atoms_to_compare_[ count_atoms ] = core::id::AtomID( jj, ii );
			++count_atoms;
		}
	}
}

void
DownstreamRMSEvaluator::set_n_geometric_constraints( Size setting )
{
	n_geometric_constraints_ = setting;
	ds_builders_.resize( n_geometric_constraints_ );
}

/// @details upstream-only downstream algorithms are handled by having null-pointing
/// downstream builder pointers.  The score() method checks that the ds_builder is not
/// null;
void
DownstreamRMSEvaluator::set_downstream_builder(
	Size which_geom_cst,
	downstream::DownstreamBuilderCOP ds_builder
)
{
	ds_builders_[ which_geom_cst ] = ds_builder;

}

DownstreamRMSEvaluator::Real
DownstreamRMSEvaluator::score( match const & m ) const
{
	utility::vector1< utility::vector1< Vector > > ds_coords( n_geometric_constraints_ );
	for ( Size ii = 1; ii <= n_geometric_constraints_; ++ii ) {
		ds_coords[ ii ].resize( atoms_to_compare_.size() );
	}

	for ( Size ii = 1; ii <= n_geometric_constraints_; ++ii ) {
		if ( ds_builders_[ ii ] ) {
			ds_builders_[ ii ]->coordinates_from_hit( m[ ii ], atoms_to_compare_, ds_coords[ ii ] );
		}
	}

	Real rms_sum = 0.0;
	for ( Size ii = 1; ii <= n_geometric_constraints_; ++ii ) {
		if ( ! ds_builders_[ ii ] ) continue;
		for ( Size jj = ii+1; jj <= n_geometric_constraints_; ++jj ) {
			if ( ! ds_builders_[ jj ] ) continue;

			Real iijj_rms = 0.0;
			Size num_atoms_to_compare=atoms_to_compare_.size();
			for ( Size kk = 1; kk <= num_atoms_to_compare; ++kk ) {
				iijj_rms += ds_coords[ ii ][ kk ].distance_squared( ds_coords[ jj ][ kk ] );
			}
			rms_sum += std::sqrt( iijj_rms/num_atoms_to_compare );
		}
	}
	return rms_sum;
}

DownstreamRMSEvaluator::Real
DownstreamRMSEvaluator::score( match_dspos1 const & ) const
{
	utility_exit_with_message( "DownstreamRMSEvaluator is not able to evaluate a single-downstream-position match" );
	return 0.0;
}

}
}
}

