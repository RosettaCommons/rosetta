// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file MotifHit.cc
/// @brief Class that holds information about a motif in the context of the search
/// @author sthyme (sthyme@gmail.com)

// Unit Headers
#include <protocols/motifs/MotifHit.hh>

// Package Headers
#include <protocols/motifs/Motif.hh>

// Project Headers
#include <core/conformation/Residue.hh>

#include <utility/vector1.hh>


// Utility Headers

// C++ Headers

namespace protocols {
namespace motifs{

MotifHit::MotifHit(
	Motif const & motif,
	Size const & vbpos,
	bool const passed_automorphism
) : motifcop_( motif.clone() ),
		vbpos_( vbpos ),
		passed_automorphism_( passed_automorphism ),
		build_rotamer_(0),
		target_conformer_(0)
{}

MotifHit::~MotifHit()
{}

MotifHit::MotifHit( MotifHit const & src ) :
	utility::pointer::ReferenceCount( src ),
	motifcop_( src.motifcop()->clone() ),
	vbpos_( src.vbpos_ ),
	passed_automorphism_( src.passed_automorphism_ ),
	build_rotamer_( src.build_rotamer()->clone() ),
	target_conformer_( src.target_conformer()->clone() )
{}

MotifHitOP
MotifHit::clone() const
{
	return new MotifHit(*this);
}

void
MotifHit::final_test(
	core::Real const & final_test
)
{
	final_test_ = final_test;
}

void
MotifHit::build_rotamer(
	core::conformation::Residue const & build_rotamer
)
{
	build_rotamer_ = build_rotamer.clone();
}

void
MotifHit::target_conformer(
	core::conformation::Residue const & target_conformer
)
{
	target_conformer_ = target_conformer.clone();
}

} // namespace motifs
} // namespace protocols
