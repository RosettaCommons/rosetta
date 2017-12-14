// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

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

static basic::Tracer TR( "core.scoring.constraints.NamedAngleConstraint" );

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

NamedAngleConstraint::NamedAngleConstraint(
	id::NamedAtomID const & a1,
	id::NamedAtomID const & a2,
	id::NamedAtomID const & a3,
	func::FuncOP func,
	ScoreType scoretype
) :
	AngleConstraint( id::AtomID( 0, a1.rsd() ), id::AtomID( 0, a2.rsd() ), id::AtomID( 0, a3.rsd() ), func, scoretype ),
	named_atom1_( a1 ),
	named_atom2_( a2 ),
	named_atom3_( a3 ),
	type1_id_( 0 ),
	type2_id_( 0 ),
	type3_id_( 0 )
{
}

NamedAngleConstraint::NamedAngleConstraint(NamedAngleConstraint const & /*other*/) = default;



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
		return nullptr;
	}
}

/// @brief This overrride updates the sequence numbering but not the atom names.
/// @author Vikram K. Mulligan (vmullig@uw.edu).
ConstraintOP
NamedAngleConstraint::remap_resid(
	core::id::SequenceMapping const &seqmap
) const {
	if ( seqmap[named_atom1_.rsd()] !=0 && seqmap[named_atom2_.rsd()] !=0 && seqmap[named_atom3_.rsd()] !=0 ) {
		core::id::NamedAtomID newat1( named_atom1_.atom(), seqmap[named_atom1_.rsd()] );
		core::id::NamedAtomID newat2( named_atom2_.atom(), seqmap[named_atom2_.rsd()] );
		core::id::NamedAtomID newat3( named_atom3_.atom(), seqmap[named_atom3_.rsd()] );
		NamedAngleConstraintOP new_cst( new NamedAngleConstraint( newat1, newat2, newat3, func()->clone(), score_type() ) );
		debug_assert( new_cst );
		return ConstraintOP( new_cst );
	}

	return ConstraintOP( /*NULL*/ );
}

bool NamedAngleConstraint::operator == ( Constraint const & rhs ) const {
	// base class operator== ensures that both classes are of type NamedAngleConstraint
	// through the mutual invocation of same_type_as_me
	if ( ! AngleConstraint::operator == ( rhs ) ) return false;

	auto const & rhs_napc( static_cast< NamedAngleConstraint const & > ( rhs ) );
	if ( named_atom1_ != rhs_napc.named_atom1_ ) return false;
	if ( named_atom2_ != rhs_napc.named_atom2_ ) return false;
	if ( named_atom3_ != rhs_napc.named_atom3_ ) return false;

	if ( type1_id_ != rhs_napc.type1_id_ ) return false;
	if ( type2_id_ != rhs_napc.type2_id_ ) return false;
	if ( type3_id_ != rhs_napc.type3_id_ ) return false;

	return true;
}

bool NamedAngleConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< NamedAngleConstraint const * > ( &other );
}

void
NamedAngleConstraint::setup_for_scoring(  func::XYZ_Func const & xyz, ScoreFunction const& ) const {
	// if ( pose_chemical_checksum_ == pose.get_current_chemical_checksum() )
	auto type1_id_now = (core::Size)&( xyz.residue( named_atom1_.rsd() ).type() );
	auto type2_id_now = (core::Size)&( xyz.residue( named_atom2_.rsd() ).type() );
	auto type3_id_now = (core::Size)&( xyz.residue( named_atom3_.rsd() ).type() );
	if ( type1_id_ != type1_id_now || type2_id_ != type2_id_now  || type3_id_ != type3_id_now ) {
		atom1( id::AtomID( xyz.residue( named_atom1_.rsd() ).atom_index( named_atom1_.atom() ), named_atom1_.rsd() ));
		atom2( id::AtomID( xyz.residue( named_atom2_.rsd() ).atom_index( named_atom2_.atom() ), named_atom2_.rsd() ));
		atom3( id::AtomID( xyz.residue( named_atom3_.rsd() ).atom_index( named_atom3_.atom() ), named_atom3_.rsd() ));
		if ( !atom1().valid() || !atom2().valid() || !atom3().valid() ) {
			TR.Warning << "can't find atom for constraint"; show_def_nopose( TR.Warning );
			TR.Warning << std::endl;
		}
		if ( !atom1().valid() ) {
			throw CREATE_EXCEPTION(core::id::EXCN_AtomNotFound, named_atom1_ );
		}
		if ( !atom2().valid() ) {
			throw CREATE_EXCEPTION(core::id::EXCN_AtomNotFound, named_atom2_ );
		}
		if ( !atom3().valid() ) {
			throw CREATE_EXCEPTION(core::id::EXCN_AtomNotFound, named_atom3_ );
		}
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

core::Real
NamedAngleConstraint::dist( core::scoring::func::XYZ_Func const & xyz ) const {
	return angle(
		xyz.residue( named_atom1_.rsd() ).xyz( named_atom1_.atom() ),
		xyz.residue( named_atom2_.rsd() ).xyz( named_atom2_.atom() ),
		xyz.residue( named_atom3_.rsd() ).xyz( named_atom3_.atom() )
	);
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
	if ( (res1 > pose.size()) || (res2 > pose.size()) || (res3 > pose.size()) ) {
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

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::constraints::NamedAngleConstraint::NamedAngleConstraint() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::NamedAngleConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< AngleConstraint >( this ) );
	arc( CEREAL_NVP( named_atom1_ ) ); // id::NamedAtomID
	arc( CEREAL_NVP( named_atom2_ ) ); // id::NamedAtomID
	arc( CEREAL_NVP( named_atom3_ ) ); // id::NamedAtomID
	arc( CEREAL_NVP( type1_id_ ) ); // core::Size
	arc( CEREAL_NVP( type2_id_ ) ); // core::Size
	arc( CEREAL_NVP( type3_id_ ) ); // core::Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::NamedAngleConstraint::load( Archive & arc ) {
	arc( cereal::base_class< AngleConstraint >( this ) );
	arc( named_atom1_ ); // id::NamedAtomID
	arc( named_atom2_ ); // id::NamedAtomID
	arc( named_atom3_ ); // id::NamedAtomID
	arc( type1_id_ ); // core::Size
	arc( type2_id_ ); // core::Size
	arc( type3_id_ ); // core::Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::NamedAngleConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::NamedAngleConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_NamedAngleConstraint )
#endif // SERIALIZATION
