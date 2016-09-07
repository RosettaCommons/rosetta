// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file BuildPosition.cc
/// @brief Position where motifs are built in MotifSearch
/// @author sthyme (sthyme@gmail.com)

// Unit Headers
#include <protocols/motifs/BuildPosition.hh>

// Package Headers
#include <protocols/motifs/Motif.hh>
#include <protocols/motifs/MotifHit.hh>

// Project Headers
#include <core/conformation/Residue.hh>

#include <utility>
#include <utility/vector1.hh>


// Utility Headers

// C++ Headers

namespace protocols {
namespace motifs {

BuildPosition::BuildPosition(
	Size const seqpos,
	utility::vector1< Size > const & target_positions
) : seqpos_( seqpos ),
	target_positions_( target_positions ),
	allowed_types_(),
	best_rotamers_(/* 0 */),
	best_motifs_(/* 0 */),
	best_motifhits_(/* 0 */)
{}


BuildPosition::BuildPosition(
	Size const seqpos,
	utility::vector1< Size > const & target_positions,
	std::set< std::string >  allowed_types
) : seqpos_( seqpos ),
	target_positions_( target_positions ),
	allowed_types_(std::move( allowed_types )),
	best_rotamers_(/* 0 */),
	best_motifs_(/* 0 */),
	best_motifhits_(/* 0 */)
{}

BuildPosition::~BuildPosition() = default;

BuildPosition::BuildPosition( BuildPosition const & src ) :
	utility::pointer::ReferenceCount( src ),
	seqpos_( src.seqpos() ),
	target_positions_( src.target_positions() ),
	allowed_types_( src.allowed_types() ),
	best_rotamers_( src.best_rotamers() ),
	best_motifs_( src.best_motifs() ),
	best_motifhits_( src.best_motifhits() )
{}

BuildPositionOP
BuildPosition::clone() const
{
	return BuildPositionOP( new BuildPosition(*this) );
}

void
BuildPosition::keep_rotamer(
	core::conformation::Residue const & res
)
{
	best_rotamers_.push_back( res.clone() );
}

void
BuildPosition::keep_motif(
	Motif const & motif
)
{
	best_motifs_.push_back( motif.clone() );
}

void
BuildPosition::keep_motifhit(
	MotifHit const & motifhit
)
{
	best_motifhits_.push_back( motifhit.clone() );
}

void
BuildPosition::clear_data()
{
	best_rotamers_.clear();
	best_motifs_.clear();
}

void
BuildPosition::clear_rots()
{
	best_rotamers_.clear();
}

} // namespace motifs
} // namespace protocols
