// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constraints/ResidueTypeConstraint.cc
///
/// @brief
/// @author Sarel Fleishman


#include <core/scoring/constraints/ResidueTypeConstraint.hh>

#include <core/conformation/Residue.hh>
//#include <core/io/pdb/pose_io.hh> -- REALLY?
// AUTO-REMOVED #include <basic/options/option.hh>
// AUTO-REMOVED #include <basic/options/keys/OptionKeys.hh>
#include <core/scoring/ScoreType.hh>
// AUTO-REMOVED #include <core/scoring/constraints/ConstraintSet.hh>
#include <basic/Tracer.hh>

#include <core/id/SequenceMapping.hh>
#include <core/pose/util.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/pose/Pose.hh>

namespace core {
namespace scoring {
namespace constraints {


static basic::Tracer TR("core.scoring.constraints.ResidueTypeConstraint");

ResidueTypeConstraint::ResidueTypeConstraint():
			Constraint( core::scoring::res_type_constraint )
{}

ResidueTypeConstraint::ResidueTypeConstraint(
	core::pose::Pose const & pose,
	Size seqpos,
	core::Real favor_native_bonus
):
	Constraint( core::scoring::res_type_constraint ),
	seqpos_( seqpos ),
	rsd_type_name3_( pose.residue_type( seqpos ).name3() ),
	favor_native_bonus_( favor_native_bonus )
{}


ResidueTypeConstraint::ResidueTypeConstraint(
	core::pose::Pose const &, //pose,
	Size seqpos,
	std::string AAname,
	core::Real favor_native_bonus
):
	Constraint( core::scoring::res_type_constraint ),
	seqpos_( seqpos ),
	rsd_type_name3_( AAname ),
	favor_native_bonus_( favor_native_bonus )
{}

ResidueTypeConstraint::ResidueTypeConstraint(
		Size seqpos,
		std::string aa_in,
		std::string name3_in,
		core::Real bonus_in
):
	Constraint( core::scoring::res_type_constraint ),
	seqpos_( seqpos ),
	AAname( aa_in ),
	rsd_type_name3_( name3_in ),
	favor_native_bonus_( bonus_in )
{}


ResidueTypeConstraint::~ResidueTypeConstraint() {}

ConstraintOP
ResidueTypeConstraint::clone() const
{
	return ConstraintOP( new ResidueTypeConstraint( *this ) );
}

utility::vector1< core::Size >
ResidueTypeConstraint::residues() const {
	utility::vector1< core::Size > pos_list(1, seqpos_); // length 1 containing "all" seqpos_ values
	return pos_list;
}

void
ResidueTypeConstraint::show( std::ostream & out ) const {
	out << "ResidueTypeConstraint; ";
	out << "seqpos: " << seqpos_;
	out << "; AAname: "<< AAname;
	out << "; rsd_type_name3: "<< rsd_type_name3_;
	out << "; favor_native_bonus: "<< favor_native_bonus_;
}

ConstraintOP
ResidueTypeConstraint::remap_resid( core::id::SequenceMapping const &seqmap ) const
{
	core::Size newseqpos = seqmap[ seqpos_ ];
  if ( newseqpos != 0 ) {
		return new ResidueTypeConstraint(	newseqpos, AAname, rsd_type_name3_, favor_native_bonus_ );
  } else {
    return NULL;
  }
}

bool
ResidueTypeConstraint::operator == ( Constraint const & other_cst ) const
{
	if( !dynamic_cast< ResidueTypeConstraint const * > ( &other_cst ) ) return false;

	ResidueTypeConstraint const & other( static_cast< ResidueTypeConstraint const & > (other_cst) );

	if( seqpos_ != other.seqpos_ ) return false;
	if( AAname != other.AAname ) return false;
	if( rsd_type_name3_ != other.rsd_type_name3_ ) return false;
	if( favor_native_bonus_ != other.favor_native_bonus_ ) return false;
	if( this->score_type() != other.score_type() ) return false;

	return true;
}

ConstraintOP
ResidueTypeConstraint::remapped_clone( pose::Pose const&, pose::Pose const&, id::SequenceMappingCOP smap ) const {

	core::Size newseqpos = seqpos_;
	if ( smap ) {
		newseqpos = (*smap)[ seqpos_ ];
		if( newseqpos == 0 ) return NULL;
	}

	return new ResidueTypeConstraint(newseqpos, AAname, rsd_type_name3_, favor_native_bonus_);
}


// Calculates a score for this constraint using XYZ_Func, and puts the UNWEIGHTED score into
// emap. Although the current set of weights currently is provided, Constraint objects
// should put unweighted scores into emap.
void
ResidueTypeConstraint::score( func::XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const
{
	Real const weight(weights[ this->score_type() ] );
	if( weight == 0 ) return; // what's the point?

	conformation::Residue const & rsd( xyz_func.residue(seqpos_) );
	if( rsd.type().name3() == rsd_type_name3_ )
		emap[ this->score_type() ] -= favor_native_bonus_;
	// no match, don't adjust score
}


void
ResidueTypeConstraint::fill_f1_f2(
	AtomID const & ,//atom,
	func::XYZ_Func const &,
	Vector & ,//F1,
	Vector & ,//F2,
	EnergyMap const & //weights
) const
{
	// Do nothing.
	// Derivative of this restraint is effectively zero
	// so we just "add zero" to F1 and F2.
}

} // namespace constraints
} // namespace scoring
} // namespace core
