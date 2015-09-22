// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/constraints/NamedAngleConstraint.cc
///
/// @brief Provides an angle constraint based on atom names, rather than numbers. Useful for when atom numbers change
/// @author Tom Linsky (tlinsky at uw dot edu)

// Unit Headers
#include <core/scoring/constraints/NamedAngleConstraint.hh>

// Package Headers
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/FuncFactory.hh>

// Project Headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Residue.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/id/Exceptions.hh>

// Utility Headers
#include <basic/Tracer.hh>

#include <core/id/SequenceMapping.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <utility/vector1.hh>

//Auto Headers
#include <core/id/types.hh>

static THREAD_LOCAL basic::Tracer TR( "core.scoring.constraints.NamedAngleConstraint" );

namespace core {
namespace scoring {
namespace constraints {

NamedAngleConstraint::NamedAngleConstraint(
	id::NamedAtomID const & a1,
	id::NamedAtomID const & a2,
	id::NamedAtomID const & a3,
	func::FuncOP func,
	ScoreType scoretype ):
	AngleConstraint( id::AtomID( 0, a1.rsd() ), id::AtomID( 0, a2.rsd() ), id::AtomID( 0, a3.rsd() ), func, scoretype ),
	named_atom1_( a1 ),
	named_atom2_( a2 ),
	named_atom3_( a3 )
{
}

std::string
NamedAngleConstraint::type() const
{
	return "NamedAngleConstraint";
}

ConstraintOP
NamedAngleConstraint::clone() const
{
	return ConstraintOP( new NamedAngleConstraint( *this ) );
}

/// @brief Copies the data from this Constraint into a new object and returns an OP
/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
/// to the new object. Intended to be implemented by derived classes.
ConstraintOP
NamedAngleConstraint::remapped_clone(
	pose::Pose const &,
	pose::Pose const& dest,
	id::SequenceMappingCOP smap ) const
{
	id::NamedAtomID atom1( named_atom1_ );
	id::NamedAtomID atom2( named_atom2_ );
	id::NamedAtomID atom3( named_atom3_ );

	if ( smap ) {
		atom1.rsd() = (*smap)[ atom1.rsd() ];
		atom2.rsd() = (*smap)[ atom2.rsd() ];
		atom3.rsd() = (*smap)[ atom3.rsd() ];
	}

	//get AtomIDs for target pose
	id::AtomID id1( core::pose::named_atom_id_to_atom_id( atom1, dest ) );
	id::AtomID id2( core::pose::named_atom_id_to_atom_id( atom2, dest ) );
	id::AtomID id3( core::pose::named_atom_id_to_atom_id( atom3, dest ) );
	if ( id1.valid() && id2.valid() && id3.valid() ) {
		return ConstraintOP( new NamedAngleConstraint( atom1, atom2, atom3, func()->clone(), score_type() ) );
	} else {
		return NULL;
	}
}

void
NamedAngleConstraint::score( func::XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const
{
	emap[ this->score_type() ] += AngleConstraint::score(
		xyz.residue( named_atom1_.rsd() ).xyz( named_atom1_.atom() ),
		xyz.residue( named_atom2_.rsd() ).xyz( named_atom2_.atom() ),
		xyz.residue( named_atom3_.rsd() ).xyz( named_atom3_.atom() ) );
}

void
NamedAngleConstraint::show_def( std::ostream & out, pose::Pose const & ) const
{
	show_def_nopose( out );
}

void
NamedAngleConstraint::show_def_nopose( std::ostream & out ) const
{
	out << type() << " " << named_atom1_ << " " << named_atom2_ << " " << named_atom3_ << " ";
	func()->show_definition( out );
}

/// @details one line definition "AngleConstraint atom1 res1 atom2 res2 atom3 res3 function_type function_definition"
void
NamedAngleConstraint::read_def(
	std::istream & data,
	core::pose::Pose const & pose,
	func::FuncFactory const & func_factory )
{
	Size res1, res2, res3;
	std::string tempres1, tempres2, tempres3;
	std::string name1, name2, name3;
	std::string func_type;

	data >> name1 >> tempres1
		>> name2 >> tempres2
		>> name3 >> tempres3
		>> func_type;

	ConstraintIO::parse_residue( pose, tempres1, res1 );
	ConstraintIO::parse_residue( pose, tempres2, res2 );
	ConstraintIO::parse_residue( pose, tempres3, res3 );

	TR.Debug << "read: " << name1 << " " << name2 << " " << name3 << " "
		<< res1 << " " << res2 << " " << res3 << " func: " << func_type << std::endl;
	if ( (res1 > pose.total_residue()) || (res2 > pose.total_residue()) || (res3 > pose.total_residue()) ) {
		TR.Warning  << "ignored constraint (no such atom in pose!)"
			<< name1 << " " << name2 << " " << name3 << " " << res1 << " " << res2 << " " << res3 << std::endl;
		data.setstate( std::ios_base::failbit );
		return;
	}

	named_atom1_ = id::NamedAtomID( name1, res1 );
	named_atom2_ = id::NamedAtomID( name2, res2 );
	named_atom3_ = id::NamedAtomID( name3, res3 );

	// check atom ids
	id::AtomID a1 = id::AtomID( core::pose::named_atom_id_to_atom_id(named_atom1_, pose ) );
	id::AtomID a2 = id::AtomID( core::pose::named_atom_id_to_atom_id(named_atom2_, pose ) );
	id::AtomID a3 = id::AtomID( core::pose::named_atom_id_to_atom_id(named_atom3_, pose ) );
	if ( (a1.atomno() == 0) || (a2.atomno() == 0) || (a3.atomno() == 0) ) {
		TR.Warning << "Error reading atoms: read in atom names("
			<< name1 << "," << name2 << "," << name3 << "), "
			<< "and found AtomIDs (" << a1 << "," << a2 << "," << a3 << ")" << std::endl;
		data.setstate( std::ios_base::failbit );
		return;
	}

	func::FuncOP f = func_factory.new_func( func_type );
	f->read_data( data );
	set_func( f );

	while ( data.good() && (data.get() != '\n') ) {}

	if ( TR.Debug.visible() ) {
		func()->show_definition( TR.Debug );
		TR.Debug << std::endl;
	}
} // read_def

}
}
}
