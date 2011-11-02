// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constraints/ResidueTypeLinkingConstraint.cc
///
/// @brief
/// @author Sarel Fleishman


#include <core/scoring/constraints/ResidueTypeLinkingConstraint.hh>

#include <core/conformation/Residue.hh>
#include <core/scoring/ScoreType.hh>
#include <basic/Tracer.hh>

#include <core/id/AtomID.hh>
#include <core/kinematics/Jump.hh>
#include <core/pose/Pose.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/constraints/XYZ_Func.hh>
#include <utility/vector1.hh>


namespace core {
namespace scoring {
namespace constraints {


static basic::Tracer TR("core.scoring.constraints.ResidueTypeLinkingConstraint");

ResidueTypeLinkingConstraint::ResidueTypeLinkingConstraint():
			Constraint( core::scoring::res_type_linking_constraint )
{}

ResidueTypeLinkingConstraint::ResidueTypeLinkingConstraint(
	core::pose::Pose const & pose,
	Size seqpos1,
	Size seqpos2,
	core::Real bonus
):
	Constraint( core::scoring::res_type_linking_constraint ),
	seqpos1_( seqpos1 ),
	seqpos2_( seqpos2 ),
	rsd1_type_name3_( pose.residue_type( seqpos1 ).name3() ),
	rsd2_type_name3_( pose.residue_type( seqpos2 ).name3() ),
	bonus_( bonus )
{
		//conformation::Residue const & rsd( pose.residue(seqpos_) );
		//flo feb 11 by default this constraint now only has
		// one atom, good enough since we're not doing minimization
		//anyway. having all heavyatoms fucks up remapping this constraint
		//from one residue type to another
		atom_ids_.push_back(AtomID(1, seqpos1_));
		atom_ids_.push_back(AtomID(1, seqpos2_));
		//for(Size i = 1, i_end = rsd.nheavyatoms(); i <= i_end; ++i) {
		//	atom_ids_.push_back(AtomID( i, seqpos_ )); // atom, rsd
		//}
}


ResidueTypeLinkingConstraint::ResidueTypeLinkingConstraint(
	core::pose::Pose const &, //pose,
	Size seqpos1,
	Size seqpos2,
	std::string AA1name,
    std::string AA2name,
	core::Real bonus
):
	Constraint( core::scoring::res_type_linking_constraint ),
	seqpos1_( seqpos1 ),
	seqpos2_( seqpos2 ),
	rsd1_type_name3_( AA1name ),
	rsd2_type_name3_( AA2name ),
	bonus_( bonus )
{
		atom_ids_.push_back(AtomID(1, seqpos1_));
		atom_ids_.push_back(AtomID(1, seqpos2_));
		//conformation::Residue const & rsd( pose.residue(seqpos_) );
		//for(Size i = 1, i_end = rsd.nheavyatoms(); i <= i_end; ++i) {
		//	atom_ids_.push_back(AtomID( i, seqpos_ )); // atom, rsd
		//}
}

ResidueTypeLinkingConstraint::~ResidueTypeLinkingConstraint() {}

ConstraintOP
ResidueTypeLinkingConstraint::clone() const
{
	return ConstraintOP( new ResidueTypeLinkingConstraint( *this ) );
}

Size
ResidueTypeLinkingConstraint::natoms() const
{
	return atom_ids_.size();
}

void
ResidueTypeLinkingConstraint::show( std::ostream & out ) const {
	out << "ResidueTypeLinkingConstraint; ";
	out << "seqpos1: " << seqpos1_;
	out << "seqpos2: " << seqpos2_;
	out << "; AA1name: "<< AA1name;
	out << "; AA2name: "<< AA2name;
	out << "; rsd1_type_name3: "<< rsd1_type_name3_;
	out << "; rsd2_type_name3: "<< rsd2_type_name3_;
	out << "; favor_native_bonus: "<< bonus_;
}

id::AtomID const &
ResidueTypeLinkingConstraint::atom( Size const index ) const
{
	return atom_ids_[index];
}

/*
ConstraintOP
ResidueTypeLinkingConstraint::remap_resid( core::id::SequenceMapping const &seqmap ) const
{
	core::Size newseqpos = seqmap[ seqpos_ ];
  if ( newseqpos != 0 ) {

		utility::vector1< id::AtomID > new_atomids;
		for( utility::vector1< id::AtomID >::const_iterator at_it = atom_ids_.begin(); at_it != atom_ids_.end(); ++at_it ){
			if( seqmap[ at_it->rsd() ] != 0 ){
				new_atomids.push_back( id::AtomID( at_it->atomno(), seqmap[ at_it->rsd() ] ) );
			}
		}
		return ConstraintOP( new ResidueTypeLinkingConstraint(	newseqpos, AAname, rsd_type_name3_, favor_native_bonus_, new_atomids ) );
  } else {
    return NULL;
  }
}
*/
bool
ResidueTypeLinkingConstraint::operator == ( Constraint const & other_cst ) const
{
	if( !dynamic_cast< ResidueTypeLinkingConstraint const * > ( &other_cst ) ) return false;

	ResidueTypeLinkingConstraint const & other( static_cast< ResidueTypeLinkingConstraint const & > (other_cst) );

	if( seqpos1_ != other.seqpos1_ ) return false;
	if( seqpos2_ != other.seqpos2_ ) return false;
	if( AA1name != other.AA1name ) return false;
	if( AA2name != other.AA2name ) return false;
	if( rsd1_type_name3_ != other.rsd1_type_name3_ ) return false;
	if( rsd2_type_name3_ != other.rsd2_type_name3_ ) return false;
	if( bonus_ != other.bonus_ ) return false;
	core::Size natoms = atom_ids_.size();
	if( natoms != other.atom_ids_.size() ) return false;
	for(core::Size i =1; i <= natoms; ++i){
		if( atom_ids_[i] != other.atom_ids_[i] ) return false;
	}
	if( this->score_type() != other.score_type() ) return false;

	return true;
}
/*
ConstraintOP
ResidueTypeLinkingConstraint::remapped_clone( pose::Pose const& src, pose::Pose const& dest, id::SequenceMappingCOP smap ) const {

	core::Size newseqpos = seqpos_;
	if ( smap ) {
		newseqpos = (*smap)[ seqpos_ ];
		if( newseqpos == 0 ) return NULL;
	}

	utility::vector1< AtomID > new_atomids;
	for( core::Size  i =1; i <= atom_ids_.size(); ++i ){
		id::NamedAtomID atom( atom_id_to_named_atom_id( atom_ids_[i], src ) );
		if( smap ){
			atom.rsd() = (*smap)[atom_ids_[i].rsd() ];
		}
		new_atomids.push_back( named_atom_id_to_atom_id( atom, dest ) );
	}
	return new ResidueTypeLinkingConstraint(newseqpos, AAname, rsd_type_name3_, favor_native_bonus_, new_atomids);
}
*/

// Calculates a score for this constraint using XYZ_Func, and puts the UNWEIGHTED score into
// emap. Although the current set of weights currently is provided, Constraint objects
// should put unweighted scores into emap.
void
ResidueTypeLinkingConstraint::score( XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const
{
	Real const weight(weights[ this->score_type() ] );
	if( weight == 0 ) return; // what's the point?

	conformation::Residue const & rsd1( xyz_func.residue(seqpos1_) );
	conformation::Residue const & rsd2( xyz_func.residue(seqpos2_) );
	if( rsd1.type().name3() == rsd2.type().name3() )
		emap[ this->score_type() ] -= bonus_;
	// no match, don't adjust score
}


void
ResidueTypeLinkingConstraint::fill_f1_f2(
	AtomID const & ,//atom,
	XYZ_Func const &,
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
