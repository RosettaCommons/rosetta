// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/NamedDihedralConstraint.cc
///
/// @brief Provides a dihedral constraint based on atom names, rather than numbers. Useful for when atom numbers change
/// @author Tom Linsky (tlinsky at uw dot edu)
/// @author Andy Wtakins (amw579@nyu.edu)

// Unit Headers
#include <core/scoring/constraints/NamedDihedralConstraint.hh>

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

static THREAD_LOCAL basic::Tracer TR( "core.scoring.constraints.NamedDihedralConstraint" );

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/access.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

NamedDihedralConstraint::NamedDihedralConstraint(
	id::NamedAtomID const & a1,
	id::NamedAtomID const & a2,
	id::NamedAtomID const & a3,
	id::NamedAtomID const & a4,
	func::FuncOP func,
	ScoreType scoretype
) :
	DihedralConstraint( id::AtomID( 0, a1.rsd() ), id::AtomID( 0, a2.rsd() ), id::AtomID( 0, a3.rsd() ), id::AtomID( 0, a4.rsd() ), func, scoretype ),
	named_atom1_( a1 ),
	named_atom2_( a2 ),
	named_atom3_( a3 )
{
}

std::string
NamedDihedralConstraint::type() const
{
	return "NamedDihedralConstraint";
}

ConstraintOP
NamedDihedralConstraint::clone() const
{
	return ConstraintOP( new NamedDihedralConstraint( *this ) );
}

/// @brief Copies the data from this Constraint into a new object and returns an OP
/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
/// to the new object. Intended to be implemented by derived classes.
ConstraintOP
NamedDihedralConstraint::remapped_clone(
	pose::Pose const &,
	pose::Pose const& dest,
	id::SequenceMappingCOP smap ) const
{
	id::NamedAtomID atom1( named_atom1_ );
	id::NamedAtomID atom2( named_atom2_ );
	id::NamedAtomID atom3( named_atom3_ );
	id::NamedAtomID atom4( named_atom4_ );


	if ( smap ) {
		atom1.rsd() = (*smap)[ atom1.rsd() ];
		atom2.rsd() = (*smap)[ atom2.rsd() ];
		atom3.rsd() = (*smap)[ atom3.rsd() ];
		atom4.rsd() = (*smap)[ atom4.rsd() ];
	}

	//get AtomIDs for target pose
	id::AtomID id1( core::pose::named_atom_id_to_atom_id( atom1, dest ) );
	id::AtomID id2( core::pose::named_atom_id_to_atom_id( atom2, dest ) );
	id::AtomID id3( core::pose::named_atom_id_to_atom_id( atom3, dest ) );
	id::AtomID id4( core::pose::named_atom_id_to_atom_id( atom4, dest ) );
	if ( id1.valid() && id2.valid() && id3.valid() && id4.valid() ) {
		return ConstraintOP( new NamedDihedralConstraint( atom1, atom2, atom3, atom4, func()->clone(), score_type() ) );
	} else {
		return NULL;
	}
}

bool NamedDihedralConstraint::operator == ( Constraint const & rhs ) const {
	// base class operator== ensures that both classes are of type NamedDihedralConstraint
	// through the mutual invocation of same_type_as_me
	if ( ! DihedralConstraint::operator == ( rhs ) ) return false;

	NamedDihedralConstraint const & rhs_napc( static_cast< NamedDihedralConstraint const & > ( rhs ) );
	if ( named_atom1_ != rhs_napc.named_atom1_ ) return false;
	if ( named_atom2_ != rhs_napc.named_atom2_ ) return false;
	if ( named_atom3_ != rhs_napc.named_atom3_ ) return false;
	if ( named_atom4_ != rhs_napc.named_atom4_ ) return false;

	return true;
}

bool NamedDihedralConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< NamedDihedralConstraint const * > ( &other );
}


void
NamedDihedralConstraint::score( func::XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const
{
	emap[ this->score_type() ] += DihedralConstraint::score(
		xyz.residue( named_atom1_.rsd() ).xyz( named_atom1_.atom() ),
		xyz.residue( named_atom2_.rsd() ).xyz( named_atom2_.atom() ),
		xyz.residue( named_atom4_.rsd() ).xyz( named_atom4_.atom() ),
		xyz.residue( named_atom4_.rsd() ).xyz( named_atom4_.atom() ) );
}

void
NamedDihedralConstraint::show_def( std::ostream & out, pose::Pose const & ) const
{
	show_def_nopose( out );
}

void
NamedDihedralConstraint::show_def_nopose( std::ostream & out ) const
{
	out << type() << " " << named_atom1_ << " " << named_atom2_ << " " << named_atom3_ << " ";
	func()->show_definition( out );
}

/// @details one line definition "DihedralConstraint atom1 res1 atom2 res2 atom3 res3 function_type function_definition"
void
NamedDihedralConstraint::read_def(
	std::istream & data,
	core::pose::Pose const & pose,
	func::FuncFactory const & func_factory )
{
	Size res1, res2, res3, res4;
	std::string tempres1, tempres2, tempres3, tempres4;
	std::string name1, name2, name3, name4;
	std::string func_type;

	data >> name1 >> tempres1
		>> name2 >> tempres2
		>> name3 >> tempres3
		>> name4 >> tempres4
		>> func_type;

	ConstraintIO::parse_residue( pose, tempres1, res1 );
	ConstraintIO::parse_residue( pose, tempres2, res2 );
	ConstraintIO::parse_residue( pose, tempres3, res3 );
	ConstraintIO::parse_residue( pose, tempres4, res4 );

	TR.Debug << "read: " << name1 << " " << name2 << " " << name3 << " " << name4 << " "
		<< res1 << " " << res2 << " " << res3 << " func: " << res4 << " func: " << func_type << std::endl;
	if ( (res1 > pose.size()) || (res2 > pose.size()) || (res3 > pose.size()) || (res4 > pose.size()) ) {
		TR.Warning  << "ignored constraint (no such atom in pose!)"
			<< name1 << " " << name2 << " " << name3 << " " << name4 << " " << res1 << " " << res2 << " " << res3 << " " << res4 << std::endl;
		data.setstate( std::ios_base::failbit );
		return;
	}

	named_atom1_ = id::NamedAtomID( name1, res1 );
	named_atom2_ = id::NamedAtomID( name2, res2 );
	named_atom3_ = id::NamedAtomID( name3, res3 );
	named_atom4_ = id::NamedAtomID( name4, res4 );

	// check atom ids
	id::AtomID a1 = id::AtomID( core::pose::named_atom_id_to_atom_id(named_atom1_, pose ) );
	id::AtomID a2 = id::AtomID( core::pose::named_atom_id_to_atom_id(named_atom2_, pose ) );
	id::AtomID a3 = id::AtomID( core::pose::named_atom_id_to_atom_id(named_atom3_, pose ) );
	id::AtomID a4 = id::AtomID( core::pose::named_atom_id_to_atom_id(named_atom4_, pose ) );
	if ( (a1.atomno() == 0) || (a2.atomno() == 0) || (a3.atomno() == 0) || (a4.atomno() == 0) ) {
		TR.Warning << "Error reading atoms: read in atom names("
			<< name1 << "," << name2 << "," << name3 << "," << name4 << "), "
			<< "and found AtomIDs (" << a1 << "," << a2 << "," << a3 << "," << a4 << ")" << std::endl;
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

NamedDihedralConstraint::NamedDihedralConstraint( NamedDihedralConstraint const & src ) :
	DihedralConstraint( src ),
	named_atom1_( src.named_atom1_ ),
	named_atom2_( src.named_atom2_ ),
	named_atom3_( src.named_atom3_ ),
	named_atom4_( src.named_atom4_ )
{}


}
}
}

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::constraints::NamedDihedralConstraint::NamedDihedralConstraint() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::NamedDihedralConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< DihedralConstraint >( this ) );
	arc( CEREAL_NVP( named_atom1_ ) ); // id::NamedAtomID
	arc( CEREAL_NVP( named_atom2_ ) ); // id::NamedAtomID
	arc( CEREAL_NVP( named_atom3_ ) ); // id::NamedAtomID
	arc( CEREAL_NVP( named_atom4_ ) ); // id::NamedAtomID
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::NamedDihedralConstraint::load( Archive & arc ) {
	arc( cereal::base_class< DihedralConstraint >( this ) );
	arc( named_atom1_ ); // id::NamedAtomID
	arc( named_atom2_ ); // id::NamedAtomID
	arc( named_atom3_ ); // id::NamedAtomID
	arc( named_atom4_ ); // id::NamedAtomID
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::NamedDihedralConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::NamedDihedralConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_NamedDihedralConstraint )
#endif // SERIALIZATION
