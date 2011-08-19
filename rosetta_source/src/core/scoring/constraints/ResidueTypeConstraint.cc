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

//Auto Headers
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/constraints/XYZ_Func.hh>
#include <core/id/SequenceMapping.hh>
#include <utility/options/keys/BooleanOptionKey.hh>


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
{
		//conformation::Residue const & rsd( pose.residue(seqpos_) );
		//flo feb 11 by default this constraint now only has
		// one atom, good enough since we're not doing minimization
		//anyway. having all heavyatoms fucks up remapping this constraint
		//from one residue type to another
		atom_ids_.push_back(AtomID(1, seqpos_));
		//for(Size i = 1, i_end = rsd.nheavyatoms(); i <= i_end; ++i) {
		//	atom_ids_.push_back(AtomID( i, seqpos_ )); // atom, rsd
		//}
}


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
{
		atom_ids_.push_back(AtomID(1, seqpos_));
		//conformation::Residue const & rsd( pose.residue(seqpos_) );
		//for(Size i = 1, i_end = rsd.nheavyatoms(); i <= i_end; ++i) {
		//	atom_ids_.push_back(AtomID( i, seqpos_ )); // atom, rsd
		//}
}

ResidueTypeConstraint::ResidueTypeConstraint(
		Size seqpos,
		std::string aa_in,
		std::string name3_in,
		core::Real bonus_in,
		utility::vector1< id::AtomID > const & atoms_in
):
	Constraint( core::scoring::res_type_constraint ),
	seqpos_( seqpos ),
	AAname( aa_in ),
	rsd_type_name3_( name3_in ),
	favor_native_bonus_( bonus_in )
{
	for( utility::vector1< id::AtomID >::const_iterator at_it = atoms_in.begin(); at_it != atoms_in.end(); ++at_it ){
		atom_ids_.push_back( *at_it );
	}
}


ResidueTypeConstraint::~ResidueTypeConstraint() {}

ConstraintOP
ResidueTypeConstraint::clone() const
{
	return ConstraintOP( new ResidueTypeConstraint( *this ) );
}

Size
ResidueTypeConstraint::natoms() const
{
	return atom_ids_.size();
}

void
ResidueTypeConstraint::show( std::ostream & out ) const {
	out << "ResidueTypeConstraint; ";
	out << "seqpos: " << seqpos_;
	out << "; AAname: "<< AAname;
	out << "; rsd_type_name3: "<< rsd_type_name3_;
	out << "; favor_native_bonus: "<< favor_native_bonus_;
}

id::AtomID const &
ResidueTypeConstraint::atom( Size const index ) const
{
	return atom_ids_[index];
}


ConstraintOP
ResidueTypeConstraint::remap_resid( core::id::SequenceMapping const &seqmap ) const
{
	core::Size newseqpos = seqmap[ seqpos_ ];
  if ( newseqpos != 0 ) {

		utility::vector1< id::AtomID > new_atomids;
		for( utility::vector1< id::AtomID >::const_iterator at_it = atom_ids_.begin(); at_it != atom_ids_.end(); ++at_it ){
			if( seqmap[ at_it->rsd() ] != 0 ){
				new_atomids.push_back( id::AtomID( at_it->atomno(), seqmap[ at_it->rsd() ] ) );
			}
		}
		return ConstraintOP( new ResidueTypeConstraint(	newseqpos, AAname, rsd_type_name3_, favor_native_bonus_, new_atomids ) );
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
	core::Size natoms = atom_ids_.size();
	if( natoms != other.atom_ids_.size() ) return false;
	for(core::Size i =1; i <= natoms; ++i){
		if( atom_ids_[i] != other.atom_ids_[i] ) return false;
	}
	if( this->score_type() != other.score_type() ) return false;

	return true;
}

ConstraintOP
ResidueTypeConstraint::remapped_clone( pose::Pose const& src, pose::Pose const& dest, id::SequenceMappingCOP smap ) const {

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
	return new ResidueTypeConstraint(newseqpos, AAname, rsd_type_name3_, favor_native_bonus_, new_atomids);
}


// Calculates a score for this constraint using XYZ_Func, and puts the UNWEIGHTED score into
// emap. Although the current set of weights currently is provided, Constraint objects
// should put unweighted scores into emap.
void
ResidueTypeConstraint::score( XYZ_Func const & xyz_func, EnergyMap const & weights, EnergyMap & emap ) const
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
