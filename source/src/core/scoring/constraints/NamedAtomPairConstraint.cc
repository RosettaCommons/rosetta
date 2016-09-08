// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/NamedAtomPairConstraint.cc
///
/// @brief
/// @author Oliver Lange

// Unit Headers
#include <core/scoring/constraints/NamedAtomPairConstraint.hh>

// Package Headers
#include <core/scoring/constraints/ConstraintIO.hh>
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

static THREAD_LOCAL basic::Tracer tr( "core.io.constraints" );


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

/// @details Auto-generated virtual destructor
Obsolet_NamedAtomPairConstraint::~Obsolet_NamedAtomPairConstraint() {}

NamedAtomPairConstraint::NamedAtomPairConstraint(
	id::NamedAtomID const& a1,
	id::NamedAtomID const& a2,
	func::FuncOP func,
	ScoreType scoretype /* = atom_pair_constraint */
) :
	AtomPairConstraint( id::AtomID( 0, a1.rsd() ), id::AtomID( 0, a2.rsd() ), func, scoretype ),
	named_atom1_( a1 ),
	named_atom2_( a2 ),
	type1_id_( 0 ),
	type2_id_( 0 )
{}

ConstraintOP
NamedAtomPairConstraint::clone() const {
	/// invokes the parent copy-ctor, which performs a deep copy of the Func object
	return ConstraintOP( new NamedAtomPairConstraint( *this ));
}


/// @brief Copies the data from this Constraint into a new object and returns an OP
/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
/// to the new object. Intended to be implemented by derived classes.
ConstraintOP NamedAtomPairConstraint::remapped_clone( pose::Pose const&, pose::Pose const& dest, id::SequenceMappingCOP smap ) const {
	id::NamedAtomID named_atom1( named_atom1_ );
	id::NamedAtomID named_atom2( named_atom2_ );

	if ( smap ) {
		named_atom1.rsd() = (*smap)[ atom1().rsd() ];
		named_atom2.rsd() = (*smap)[ atom2().rsd() ];
	}

	//get AtomIDs for target pose
	id::AtomID id1( core::pose::named_atom_id_to_atom_id( named_atom1, dest ) );
	id::AtomID id2( core::pose::named_atom_id_to_atom_id( named_atom2, dest ) );
	if ( id1.valid() && id2.valid() ) {
		return ConstraintOP( new NamedAtomPairConstraint( named_atom1, named_atom2, get_func().clone(), score_type() ) );
	} else {
		return NULL;
	}
}

bool NamedAtomPairConstraint::operator == ( Constraint const & rhs ) const {
	// base class operator== ensures that both classes are of type NamedAtomPairConstraint
	// through the mutual invocation of same_type_as_me
	if ( ! AtomPairConstraint::operator == ( rhs ) ) return false;

	NamedAtomPairConstraint const & rhs_napc( static_cast< NamedAtomPairConstraint const & > ( rhs ) );
	if ( named_atom1_ != rhs_napc.named_atom1_ ) return false;
	if ( named_atom2_ != rhs_napc.named_atom2_ ) return false;
	if ( type1_id_ != rhs_napc.type1_id_ ) return false;
	if ( type2_id_ != rhs_napc.type2_id_ ) return false;

	return true;
}

bool NamedAtomPairConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< NamedAtomPairConstraint const * > ( & other );
}


void
NamedAtomPairConstraint::setup_for_scoring(  func::XYZ_Func const & xyz, ScoreFunction const& ) const {
	// if ( pose_chemical_checksum_ == pose.get_current_chemical_checksum() )
	core::Size type1_id_now = (core::Size)&( xyz.residue( named_atom1_.rsd() ).type() );
	core::Size type2_id_now = (core::Size)&( xyz.residue( named_atom2_.rsd() ).type() );
	if ( type1_id_ != type1_id_now || type2_id_ != type2_id_now ) {
		atom1( id::AtomID( xyz.residue( named_atom1_.rsd() ).atom_index( named_atom1_.atom() ), named_atom1_.rsd() ));
		atom2( id::AtomID( xyz.residue( named_atom2_.rsd() ).atom_index( named_atom2_.atom() ), named_atom2_.rsd() ));
		if ( !atom1().valid() || !atom2().valid() ) {
			tr.Warning << "[WARNING] can't find atom for constraint"; show_def_nopose( tr.Warning );
			tr.Warning << std::endl;
		}
		if ( !atom1().valid() ) {
			throw core::id::EXCN_AtomNotFound( named_atom1_ );
		}
		if ( !atom2().valid() ) {
			throw core::id::EXCN_AtomNotFound( named_atom2_ );
		}
	}
}

void NamedAtomPairConstraint::show_def( std::ostream& out, pose::Pose const& ) const {
	show_def_nopose( out );
}

void NamedAtomPairConstraint::show_def_nopose( std::ostream& out ) const {
	out << type() << " " << named_atom1_ << " " << named_atom2_ << " ";
	get_func().show_definition( out );
}

/// @details one line definition "AtomPairs atom1 res1 atom2 res2 function_type function_definition"
void
NamedAtomPairConstraint::read_def(
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
		>> func_type;

	ConstraintIO::parse_residue( pose, tempres1, res1 );
	ConstraintIO::parse_residue( pose, tempres2, res2 );

	tr.Debug << "read: " << name1 << " " << name2 << " " << res1 << " " << res2 << " func: " << func_type << std::endl;
	if ( res1 > pose.size() || res2 > pose.size() ) {
		tr.Warning  << "ignored constraint (no such atom in pose!)"
			<< name1 << " " << name2 << " " << res1 << " " << res2 << std::endl;
		data.setstate( std::ios_base::failbit );
		return;
	}

	named_atom1_ =  id::NamedAtomID( name1, res1 );
	named_atom2_ =  id::NamedAtomID( name2, res2 );
	atom1( id::AtomID( core::pose::named_atom_id_to_atom_id(named_atom1_, pose )));
	atom2( id::AtomID( core::pose::named_atom_id_to_atom_id(named_atom2_, pose )));
	if ( atom1().atomno() == 0 || atom2().atomno() == 0 ) {
		tr.Warning << "Error reading atoms: read in atom names("
			<< name1 << "," << name2 << "), "
			<< "and found AtomIDs (" << atom1() << "," << atom2() << ")" << std::endl;
		data.setstate( std::ios_base::failbit );
		return;
	}

	func::FuncOP new_func = func_factory.new_func( func_type );
	if ( ! new_func ) {
		throw utility::excn::EXCN_Msg_Exception( "The function type '" + func_type + "' is not valid" );
	}
	new_func->read_data( data );
	set_func( new_func );

	//chu skip the rest of line since this is a single line defintion.
	while ( data.good() && (data.get() != '\n') ) {}

	if ( tr.Debug.visible() ) {
		get_func().show_definition( tr.Debug );
		tr.Debug << std::endl;
	}
} // read_def


Obsolet_NamedAtomPairConstraint::Obsolet_NamedAtomPairConstraint(
	AtomPairConstraintOP cst_in,
	pose::Pose const& ref_pose ) :
	atom1_( core::pose::atom_id_to_named_atom_id(cst_in->atom( 1 ), ref_pose )),
	atom2_( core::pose::atom_id_to_named_atom_id(cst_in->atom( 2 ), ref_pose )),
	cst_( cst_in )
{}

Obsolet_NamedAtomPairConstraint::Obsolet_NamedAtomPairConstraint(
	NamedAtomID const& atom1,
	NamedAtomID const& atom2,
	AtomPairConstraintOP cst ) :
	atom1_( atom1 ),
	atom2_( atom2 ),
	cst_( cst )
{}

AtomPairConstraintOP
Obsolet_NamedAtomPairConstraint::mapto(
	id::SequenceMapping const& map,
	pose::Pose const& pose ) const
{
	//get AtomIDs for target sequence and target pose
	id::AtomID id1( core::pose::named_atom_id_to_atom_id( NamedAtomID( atom1_.atom(), map[ atom1_.rsd() ] ), pose ));
	id::AtomID id2( core::pose::named_atom_id_to_atom_id( NamedAtomID( atom2_.atom(), map[ atom2_.rsd() ] ), pose ));
	if ( id1.valid() && id2.valid() ) {
		return AtomPairConstraintOP( new AtomPairConstraint( id1, id2, cst_->get_func().clone(), cst_->score_type() ) );
	}
	return NULL; // if translation not possible
}

Obsolet_NamedAtomPairConstraintOP
Obsolet_NamedAtomPairConstraint::mapto( id::SequenceMapping const& map ) const
{
	NamedAtomID id1 ( atom1_ );
	NamedAtomID id2 ( atom2_ );
	id1.rsd() = map[ id1.rsd() ];
	id2.rsd() = map[ id2.rsd() ];
	if ( id1.valid() && id2.valid() ) {
		return Obsolet_NamedAtomPairConstraintOP( new Obsolet_NamedAtomPairConstraint( id1, id2, cst_ ) );
	}
	return NULL;
}

AtomPairConstraintOP
Obsolet_NamedAtomPairConstraint::mapto( core::pose::Pose const& pose ) const {
	id::AtomID id1( core::pose::named_atom_id_to_atom_id(atom1_, pose ));
	id::AtomID id2( core::pose::named_atom_id_to_atom_id(atom2_, pose ));
	if ( id1.valid() && id2.valid() ) {
		return AtomPairConstraintOP( new AtomPairConstraint( id1, id2, cst_->get_func().clone(), cst_->score_type() ) );
	}
	return NULL;
}

std::ostream& operator<< ( std::ostream& out, Obsolet_NamedAtomPairConstraint const& cst ) {
	return out << cst.atom1_ << " " << cst.atom2_ << " " << cst.cst_;
}


}
}
}

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::constraints::NamedAtomPairConstraint::NamedAtomPairConstraint() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::NamedAtomPairConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< AtomPairConstraint >( this ) );
	arc( CEREAL_NVP( named_atom1_ ) ); // id::NamedAtomID
	arc( CEREAL_NVP( named_atom2_ ) ); // id::NamedAtomID
	arc( CEREAL_NVP( type1_id_ ) ); // core::Size
	arc( CEREAL_NVP( type2_id_ ) ); // core::Size
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::NamedAtomPairConstraint::load( Archive & arc ) {
	arc( cereal::base_class< AtomPairConstraint >( this ) );
	arc( named_atom1_ ); // id::NamedAtomID
	arc( named_atom2_ ); // id::NamedAtomID
	arc( type1_id_ ); // core::Size
	arc( type2_id_ ); // core::Size
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::NamedAtomPairConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::NamedAtomPairConstraint )


/// @brief Default constructor required by cereal to deserialize this class
core::scoring::constraints::Obsolet_NamedAtomPairConstraint::Obsolet_NamedAtomPairConstraint() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::Obsolet_NamedAtomPairConstraint::save( Archive & arc ) const {
	arc( CEREAL_NVP( atom1_ ) ); // id::NamedAtomID
	arc( CEREAL_NVP( atom2_ ) ); // id::NamedAtomID
	arc( CEREAL_NVP( cst_ ) ); // AtomPairConstraintOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::Obsolet_NamedAtomPairConstraint::load( Archive & arc ) {
	arc( atom1_ ); // id::NamedAtomID
	arc( atom2_ ); // id::NamedAtomID
	arc( cst_ ); // AtomPairConstraintOP
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::Obsolet_NamedAtomPairConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::Obsolet_NamedAtomPairConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_NamedAtomPairConstraint )
#endif // SERIALIZATION
