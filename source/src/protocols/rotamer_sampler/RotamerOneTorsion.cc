// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/rotamer_sampler/RotamerOneTorsion.hh
/// @brief Generate glycosidic chi rotamers for RNA.
/// @author Fang-Chieh Chou

// Unit headers
#include <protocols/rotamer_sampler/RotamerOneTorsion.hh>

// Project headers
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/random/random.hh>

using namespace core;
static basic::Tracer TR( "protocols.rotamer_sampler.RotamerOneTorsion" );
static numeric::random::RandomGenerator RG( 2560199 );  // Magic number

namespace protocols {
namespace rotamer_sampler {
///////////////////////////////////////////////////////////////////////////
RotamerOneTorsion::RotamerOneTorsion():
	RotamerSized(),
	id_( 0 ),
	torsion_id_()
{}

RotamerOneTorsion::RotamerOneTorsion(
		core::id::TorsionID const tor_id,
		TorsionList const allowed_torsions
):
	RotamerSized(),
	id_( 0 ),
	torsions_( allowed_torsions ),
	torsion_id_( tor_id )
{}

RotamerOneTorsion::~RotamerOneTorsion(){}
///////////////////////////////////////////////////////////////////////////
void RotamerOneTorsion::reset() {
	runtime_assert( is_init() );
	if ( random() ) {
		++( *this );
	} else {
		id_ = 1;
	}
}
///////////////////////////////////////////////////////////////////////////
void RotamerOneTorsion::operator++() {
	runtime_assert( not_end() );
	if ( random() ) {
		id_ = RG.random_range( 1, size() );
	} else {
		++id_;
	}
}
///////////////////////////////////////////////////////////////////////////
bool RotamerOneTorsion::not_end() const {
	runtime_assert( is_init() );
	if ( random() ) return true;
	return ( id_ <= size() );
}
///////////////////////////////////////////////////////////////////////////
void RotamerOneTorsion::apply( core::pose::Pose & pose, Size const i ) {
	pose.set_torsion( torsion_id_, value( i ) );
}
///////////////////////////////////////////////////////////////////////////

}
}
