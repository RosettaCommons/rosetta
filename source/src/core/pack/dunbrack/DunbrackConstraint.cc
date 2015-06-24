// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/dunbrack/DunbrackConstraint.cc
/// @brief
/// @author James Thompson

// Unit headers
#include <core/pack/dunbrack/DunbrackConstraint.hh>
#include <core/pack/dunbrack/DunbrackConstraintCreator.hh>


#include <core/id/AtomID.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/pack/dunbrack/RotamerLibrary.hh>
#include <core/pack/rotamers/SingleResidueRotamerLibrary.hh>
#include <core/scoring/func/XYZ_Func.hh>

#include <basic/Tracer.hh>

#include <utility/vector1.hh>


namespace core {
namespace pack {
namespace dunbrack {

static thread_local basic::Tracer TR( "core.pack.dunbrack.DunbrackConstraint" );

DunbrackConstraintCreator::DunbrackConstraintCreator() {}
DunbrackConstraintCreator::~DunbrackConstraintCreator() {}

scoring::constraints::ConstraintOP
DunbrackConstraintCreator::create_constraint() const {
	return scoring::constraints::ConstraintOP( new DunbrackConstraint );
}

std::string DunbrackConstraintCreator::keyname() const
{
	return "Dunbrack";
}

DunbrackConstraint::DunbrackConstraint() :
	Constraint( core::scoring::dunbrack_constraint ),
	bonus_( 0 ),
	seqpos_( 0 ),
	rot_vec_pos_( 0 ),
	rot_bin_( 0 )
{}

DunbrackConstraint::~DunbrackConstraint() {}

std::string
DunbrackConstraint::type() const {
	return "Dunbrack";
}

Size
DunbrackConstraint::natoms() const {
	return atom_ids_.size();
}

id::AtomID const &
DunbrackConstraint::atom( Size const index ) const {
	return atom_ids_[index];
}

scoring::constraints::ConstraintOP
DunbrackConstraint::clone() const {
	return scoring::constraints::ConstraintOP( new DunbrackConstraint( *this ) );
}

// Calculates a score for this constraint using XYZ_Func, and puts the
// UNWEIGHTED score into emap. Although the current set of weights currently is
// provided, Constraint objects should put unweighted scores into emap.
void
DunbrackConstraint::score(
	scoring::func::XYZ_Func const & xyz_func,
	scoring::EnergyMap const & weights,
	scoring::EnergyMap & emap
) const {
	if ( weights[ this->score_type() ] == 0 ) return; // what's the point?

	conformation::Residue const & rsd( xyz_func.residue(seqpos_) );

	pack::dunbrack::RotVector rot;
	pack::dunbrack::rotamer_from_chi( rsd, rot );
	// don't try to restrain angles that don't exist
	if ( rot.size() > rot_vec_pos_ ) return;
	if ( rot[ rot_vec_pos_ ] == rot_bin_ ) {
		emap[ this->score_type() ] += bonus_;
	}
}

void
DunbrackConstraint::fill_f1_f2(
	AtomID const & ,//atom,
	scoring::func::XYZ_Func const &,
	Vector & ,//F1,
	Vector & ,//F2,
	scoring::EnergyMap const & //weights
) const {
	// Do nothing.
	// Derivative of this restraint is effectively zero (because it's constant
	// within each rotamer well), so we just "add zero" to F1 and F2.
	// APL Note: these derivatives will be wrong for phi and psi, since a change in phi or psi
	// will cause a change in the energy of the lowest rotamer in a well. This ought to be fixed.
}

void DunbrackConstraint::show( std::ostream & out ) const {
	out << type() << " " << seqpos_
		<< " " << rot_vec_pos_ << " " << rot_bin_
		<< " " << bonus_;
	out << std::endl;
}

/// @brief Format should look like:
/// Dunbrack seqpos_ rot_vec_pos_ rot_bin_ bonus_
void DunbrackConstraint::read_def(
	std::istream & in,
	pose::Pose const & /* pose */,
	scoring::func::FuncFactory const & /* func_factory */
) {
	in >> seqpos_ >> rot_vec_pos_ >> rot_bin_ >> bonus_;
	if ( seqpos_ == 0 || rot_vec_pos_ == 0 || rot_bin_ == 0 || bonus_ == 0 ) {
		in.setstate( std::ios_base::failbit );
	}
}

} // namespace constraints
} // namespace scoring
} // namespace core
