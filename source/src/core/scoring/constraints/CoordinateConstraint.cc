// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief

// Unit headers
#include <core/scoring/constraints/CoordinateConstraint.hh>

// Package Headers
#include <core/scoring/func/FuncFactory.hh>
#include <core/scoring/constraints/Constraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <core/scoring/ScoreType.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/id/AtomID.hh>

#include <core/conformation/Conformation.hh>

// Utility Headers
#include <basic/Tracer.hh>


// C++ Headers
#include <cstdlib>
#include <iostream>

#include <core/id/NamedAtomID.hh>
#include <core/id/SequenceMapping.hh>
#include <utility/vector1.hh>


static thread_local basic::Tracer tr( "core.io.constraints" );

namespace core {
namespace scoring {
namespace constraints {

using pose::named_atom_id_to_atom_id;
using pose::atom_id_to_named_atom_id;

CoordinateConstraint::CoordinateConstraint() :
	Constraint( coordinate_constraint ),
	atom_( id::BOGUS_ATOM_ID ),
	fixed_atom_( id::BOGUS_ATOM_ID ),
	func_( /* NULL */ ) {}

///c-tor
CoordinateConstraint::CoordinateConstraint(
	AtomID const & a1,
	AtomID const & fixed_atom_in,
	Vector const & xyz_target_in,
	func::FuncOP func,
	ScoreType scotype
):
	Constraint( scotype ),
	atom_(a1),
	fixed_atom_( fixed_atom_in ),
	xyz_target_( xyz_target_in ),
	func_( func )
{}

CoordinateConstraint::~CoordinateConstraint() {}

std::string
CoordinateConstraint::type() const {
	return "CoordinateConstraint";
}

ConstraintOP
CoordinateConstraint::clone() const {
	return ConstraintOP( new CoordinateConstraint( *this ) );
}

ConstraintOP
CoordinateConstraint::clone( core::scoring::func::FuncOP newfunc ) const
{
	return ConstraintOP( new CoordinateConstraint( atom_, fixed_atom_, xyz_target_, newfunc, score_type() ));
}



void CoordinateConstraint::show( std::ostream& out ) const
{
	out << "CoordinateConstraint ("
		<< atom_.atomno() << "," << atom_.rsd() << "-"
		<< fixed_atom_.atomno() << "," << fixed_atom_.rsd() << ")" << std::endl;
	func_->show( out );
}

void CoordinateConstraint::show_def( std::ostream& out, pose::Pose const& pose ) const
{
	out << type() << " " << pose::atom_id_to_named_atom_id( atom_, pose ) << " " << pose::atom_id_to_named_atom_id(fixed_atom_, pose ) << " ";
	out << xyz_target_.x() << " " << xyz_target_.y() << " " << xyz_target_.z() << " " ;
	func_->show_definition( out );
}


Size
CoordinateConstraint::show_violations(
	std::ostream & out,
	pose::Pose const & pose,
	Size verbose_level,
	Real threshold
) const{
	if ( verbose_level > 80 ) {
		out << "CoordConstraint ("
			<< pose.residue_type(atom_.rsd() ).atom_name( atom_.atomno() )
			<< ":" << atom_.atomno() << "," << atom_.rsd() << " ) ";
	}
	return func_->show_violations( out, dist( pose ), verbose_level, threshold );
}


/// @brief Copies the data from this Constraint into a new object and returns an OP
/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
/// to the new object. Intended to be implemented by derived classes.
ConstraintOP CoordinateConstraint::remapped_clone( pose::Pose const& src, pose::Pose const& dest, id::SequenceMappingCOP smap ) const {
	id::NamedAtomID atom1( atom_id_to_named_atom_id( atom(1), src ) );
	id::NamedAtomID atom2( atom_id_to_named_atom_id( atom(2), src ) );

	if ( smap ) {
		atom1.rsd() = (*smap)[ atom(1).rsd() ];
		atom2.rsd() = (*smap)[ atom(2).rsd() ];
	}

	//get AtomIDs for target pose
	id::AtomID id1( named_atom_id_to_atom_id( atom1, dest, false /*raise exception*/ ) );
	id::AtomID id2( named_atom_id_to_atom_id( atom2, dest, false /*raise exception*/ ) );

	// if ( atom(1) != id1 ) tr << "REMAPPING: " << atom(1) << " to " << id1 << std::endl;

	if ( id1.valid() && id2.valid() ) {
		return ConstraintOP( new CoordinateConstraint( id1, id2, xyz_target_, func_, score_type() ) );
	} else {
		return NULL;
	}
}


Real
CoordinateConstraint::non_virtual_score(
	Vector const & xyz
) const
{
	return func( xyz.distance( xyz_target_ ) );
}


void
CoordinateConstraint::score( func::XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const
{
	emap[ this->score_type() ] += non_virtual_score( xyz( atom_ ) );
}

// atom deriv
void
CoordinateConstraint::fill_f1_f2(
	AtomID const & atom,
	func::XYZ_Func const & xyz,
	Vector & F1,
	Vector & F2,
	EnergyMap const & weights
) const
{
	if ( atom != atom_ ) return;

	Vector const & xyz1( xyz( atom_ ) ), xyz2( xyz_target_ );

	Vector const f2( xyz1 - xyz2 );
	Real const dist( f2.length() ), deriv( dfunc( dist ) );
	if ( deriv != 0.0 && dist != 0.0 ) {
		Vector const f1( xyz1.cross( xyz2 ) );
		// jk: double F1 and F2 because the target is fixed
		// (matches deriv_check, and minimizes faster)
		// rhiju: No, JK, this isn't working...
		F1 += ( ( deriv / dist ) * f1 ) * weights[ this->score_type() ];
		F2 += ( ( deriv / dist ) * f2 ) * weights[ this->score_type() ];
	}
}


Size
CoordinateConstraint::natoms() const
{
	return 2;
}


core::id::AtomID const &
CoordinateConstraint::atom( Size const n ) const
{
	switch( n ) {
	case 1 :
		return atom_;
	case 2 :
		return fixed_atom_;
	default :
		utility_exit_with_message( "CoordinateConstraint::atom() bad argument" );
	}
	return atom_;
}


Real
CoordinateConstraint::dist( pose::Pose const & pose ) const {
	conformation::Conformation const & conformation( pose.conformation() );
	Vector const & xyz( conformation.xyz( atom_ ) );
	return xyz.distance( xyz_target_ );
}

ConstraintOP
CoordinateConstraint::remap_resid( core::id::SequenceMapping const &seqmap ) const
{

	if ( seqmap[atom_.rsd()] != 0 && seqmap[fixed_atom_.rsd()] != 0 ) {
		AtomID remap_a( atom_.atomno(), seqmap[atom_.rsd()] ),
			remap_fa( fixed_atom_.atomno(), seqmap[fixed_atom_.rsd()] );
		return ConstraintOP( new CoordinateConstraint( remap_a, remap_fa, xyz_target_, this->func_, score_type() ) );
	} else {
		return NULL;
	}

}

void
CoordinateConstraint::steal_def( pose::Pose const& pose ) {
	conformation::Conformation const & conformation( pose.conformation() );
	xyz_target_ = conformation.xyz( atom_ );
}

/// @details one line definition "CoordinateConstraint Atom1_Name Atom1_ResNum Atom2_Name Atom2_ResNum Atom1_target_X_coordinate Atom1_target_Y_coordinate Atom1_target_Z_coordinate func::Func_Type func::Func_Def"
void
CoordinateConstraint::read_def(
	std::istream& data,
	core::pose::Pose const& pose,
	func::FuncFactory const& func_factory
) {
	Size res1, res2;
	std::string tempres1, tempres2;
	std::string name1, name2;
	std::string func_type;

	data
		>> name1 >> tempres1
		>> name2 >> tempres2
		>> xyz_target_.x()
		>> xyz_target_.y()
		>> xyz_target_.z()
		>> func_type;

	ConstraintIO::parse_residue( pose, tempres1, res1 );
	ConstraintIO::parse_residue( pose, tempres2, res2 );

	tr.Debug << "read: " << name1 << " " << name2 << " " << res1 << " " << res2 << " func: " << func_type << std::endl;
	if ( res1 > pose.total_residue() || res2 > pose.total_residue() ) {
		tr.Warning  << "ignored constraint (no such atom in pose!)"
			<< name1 << " " << name2 << " " << res1 << " " << res2 << std::endl;
		data.setstate( std::ios_base::failbit );
		return;
	}

	atom_ = named_atom_id_to_atom_id( id::NamedAtomID( name1, res1 ), pose );
	fixed_atom_ = named_atom_id_to_atom_id( id::NamedAtomID( name2, res2 ), pose );
	if ( atom_.atomno() == 0 || fixed_atom_.atomno() == 0 ) {
		tr.Warning << "Error reading atoms: read in atom names("
			<< name1 << "," << name2 << "), "
			<< "and found AtomIDs (" << atom_ << "," << fixed_atom_ << ")" << std::endl;
		data.setstate( std::ios_base::failbit );
		return;
	}

	func_ = func_factory.new_func( func_type );
	func_->read_data( data );

	if ( data.good() ) {
		//chu skip the rest of line since this is a single line defintion.
		while ( data.good() && (data.get() != '\n') ) {}
		if ( !data.good() ) data.setstate( std::ios_base::eofbit );
	}

	if ( tr.Debug.visible() ) {
		func_->show_definition( std::cout );
		std::cout << std::endl;
	}
} // read_def

bool
CoordinateConstraint::operator == ( Constraint const & other_cst ) const
{
	if ( !dynamic_cast< CoordinateConstraint const * > ( &other_cst ) ) return false;

	CoordinateConstraint const & other( static_cast< CoordinateConstraint const & > (other_cst) );

	if ( atom_ != other.atom_ ) return false;
	if ( fixed_atom_ != other.fixed_atom_ ) return false;
	if ( func_ != other.func_ ) return false;
	if ( xyz_target_ != other.xyz_target_ ) return false;
	if ( this->score_type() != other.score_type() ) return false;

	return true;
}


// functions
Real
CoordinateConstraint::func( Real const theta ) const
{
	return func_->func( theta );
}

// deriv
Real
CoordinateConstraint::dfunc( Real const theta ) const
{
	return func_->dfunc( theta );
}


} // constraints
} // scoring
} // core
