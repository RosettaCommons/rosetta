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

// Package Headers
#include <core/scoring/func/FuncFactory.hh>
#include <core/scoring/constraints/LocalCoordinateConstraint.hh>
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

#include <core/id/NamedStubID.hh>
#include <core/id/SequenceMapping.hh>
#include <utility/vector1.hh>
#include <numeric/xyz.functions.hh>


static THREAD_LOCAL basic::Tracer tr( "core.LocalCoordinateConstraint" );

#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Numeric serialization headers
#include <numeric/xyz.serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

LocalCoordinateConstraint::LocalCoordinateConstraint() :
	Constraint( coordinate_constraint ),
	atom_( id::BOGUS_ATOM_ID ),
	fixed_stub_( id::BOGUS_STUB_ID ),
	func_( /* NULL */ ) {}

///c-tor
LocalCoordinateConstraint::LocalCoordinateConstraint(
	id::AtomID const & a1,
	id::StubID const & fixed_stub_in,
	Vector const & xyz_target_in,
	func::FuncOP func,
	ScoreType scotype /* = coordinate_constraint */
):
	Constraint( scotype ),
	atom_(a1),
	fixed_stub_( fixed_stub_in ),
	xyz_target_( xyz_target_in ),
	func_( func )
{
	runtime_assert( fixed_stub_.atom( 1 ) == fixed_stub_.center() || !fixed_stub_.center().valid() );
	//don't allow 4-atom stubs, because that changes other functions
}

LocalCoordinateConstraint::~LocalCoordinateConstraint() {}

std::string
LocalCoordinateConstraint::type() const {
	return "LocalCoordinateConstraint";
}

ConstraintOP
LocalCoordinateConstraint::clone() const {
	return ConstraintOP( new LocalCoordinateConstraint( *this ) );
}

/// @brief Copies the data from this Constraint into a new object and returns an OP
/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
/// to the new object. Intended to be implemented by derived classes.
ConstraintOP
LocalCoordinateConstraint::remapped_clone( pose::Pose const& src, pose::Pose const& dest, id::SequenceMappingCOP smap ) const {
	id::NamedAtomID atom1( core::pose::atom_id_to_named_atom_id(atom(1), src ));
	id::NamedAtomID atom2( core::pose::atom_id_to_named_atom_id(atom(2), src ));
	id::NamedAtomID atom3( core::pose::atom_id_to_named_atom_id(atom(3), src ));
	id::NamedAtomID atom4( core::pose::atom_id_to_named_atom_id(atom(4), src ));

	if ( smap ) {
		atom1.rsd() = (*smap)[ atom(1).rsd() ];
		atom2.rsd() = (*smap)[ atom(2).rsd() ];
		atom3.rsd() = (*smap)[ atom(3).rsd() ];
		atom4.rsd() = (*smap)[ atom(4).rsd() ];
	}

	//get AtomIDs for target pose
	id::AtomID id1( core::pose::named_atom_id_to_atom_id(atom1, dest ));
	id::AtomID id2( core::pose::named_atom_id_to_atom_id(atom2, dest ));
	id::AtomID id3( core::pose::named_atom_id_to_atom_id(atom3, dest ));
	id::AtomID id4( core::pose::named_atom_id_to_atom_id(atom4, dest ));
	if ( id1.valid() && id2.valid() && id3.valid() && id4.valid() ) {
		return ConstraintOP( new LocalCoordinateConstraint( id1, id::StubID( id2, id3, id4) , xyz_target_, func_ ? func_->clone() : func_, score_type() ) );
	} else {
		return NULL;
	}
}

ConstraintOP
LocalCoordinateConstraint::remap_resid( core::id::SequenceMapping const &seqmap ) const
{
	if ( seqmap[atom_.rsd()] !=0
			&& seqmap[ fixed_stub_.atom( 1 ).rsd() ] != 0
			&& seqmap[ fixed_stub_.atom( 2 ).rsd() ] != 0
			&& seqmap[ fixed_stub_.atom( 3 ).rsd() ] != 0 ) {
		id::AtomID remap_a( atom_.atomno(), seqmap[ atom_.rsd() ] );
		id::AtomID remap_s1( fixed_stub_.atom( 1 ).atomno(), seqmap[ fixed_stub_.atom( 1 ).rsd() ] );
		id::AtomID remap_s2( fixed_stub_.atom( 2 ).atomno(), seqmap[ fixed_stub_.atom( 2 ).rsd() ] );
		id::AtomID remap_s3( fixed_stub_.atom( 3 ).atomno(), seqmap[ fixed_stub_.atom( 3 ).rsd() ] );
		id::StubID remap_stub( remap_s1, remap_s2, remap_s3 );
		return ConstraintOP( new LocalCoordinateConstraint( remap_a, remap_stub, xyz_target_, this->func_, score_type() ) );
	} else {
		return NULL;
	}
}

bool LocalCoordinateConstraint::operator == ( Constraint const & rhs ) const
{
	if ( !     same_type_as_me(   rhs ) ) return false;
	if ( ! rhs.same_type_as_me( *this ) ) return false;

	LocalCoordinateConstraint const & rhs_lcc( static_cast< LocalCoordinateConstraint const & > ( rhs ));
	if ( atom_ != rhs_lcc.atom_ ) return false;
	if ( fixed_stub_ != rhs_lcc.fixed_stub_ ) return false;
	if ( xyz_target_ != rhs_lcc.xyz_target_ ) return false;
	if ( score_type() != rhs.score_type() ) return false;

	return func_ == rhs_lcc.func_ || ( func_ && rhs_lcc.func_ && *func_ == *rhs_lcc.func_ );
}

bool LocalCoordinateConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< LocalCoordinateConstraint const * > ( &other );
}



/// @details one line definition "LocalCoordinateConstraint Atom1_Name Atom1_ResNum Atom2_Name Atom3_Name Atom4_Name Atom234_ResNum Atom1_target_X_coordinate Atom1_target_Y_coordinate Atom1_target_Z_coordinate func::Func_Type func::Func_Def"
void
LocalCoordinateConstraint::read_def(
	std::istream& data,
	core::pose::Pose const& pose,
	func::FuncFactory const& func_factory
) {
	Size res1, res2;
	std::string tempres1, tempres2;
	std::string name1, name2, name3, name4;
	std::string func_type;

	data
		>> name1 >> tempres1
		>> name2 >> name3 >> name4 >> tempres2
		>> xyz_target_.x()
		>> xyz_target_.y()
		>> xyz_target_.z()
		>> func_type;

	ConstraintIO::parse_residue( pose, tempres1, res1 );
	ConstraintIO::parse_residue( pose, tempres2, res2 );

	tr.Debug << "read: " << name1 << " " << res1 << " " << name2 << " " << name3
		<< " " << name4 << " " << res2 << " func: " << func_type << std::endl;
	if ( res1 > pose.total_residue() || res2 > pose.total_residue() ) {
		tr.Warning  << "ignored constraint (no such atom in pose!)"
			<< name1 << " " << name2 << " " << res1 << " " << res2 << std::endl;
		data.setstate( std::ios_base::failbit );
		return;
	}

	atom_ = id::AtomID( core::pose::named_atom_id_to_atom_id( id::NamedAtomID( name1, res1 ), pose ));
	fixed_stub_ = id::StubID( core::pose::named_stub_id_to_stub_id( id::NamedStubID( name2, name3, name4, res2 ), pose ));
	if ( !atom_.valid() || !fixed_stub_.valid() ) {
		tr.Warning << "Error reading atoms: read in atom names("
			<< name1 << "," << name2 << "), "
			<< "and found AtomIDs (" << atom_ << "," << fixed_stub_ << ")" << std::endl;
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

void
LocalCoordinateConstraint::show_def(
	std::ostream& out,
	core::pose::Pose const& pose ) const
{
	out << type() << " " << core::pose::atom_id_to_named_atom_id( atom_, pose ) << " ";
	id::NamedStubID stubid( core::pose::stub_id_to_named_stub_id(fixed_stub_, pose ));
	out << stubid.atom( 1 ).atom() << " " <<
		stubid.atom( 2 ).atom() << " " <<
		stubid.atom( 3 ).atom() << " " << stubid.atom(1).rsd() << " ";

	out << xyz_target_.x() << " " << xyz_target_.y() << " " << xyz_target_.z() << " " ;
	func_->show_definition( out );
}

void
LocalCoordinateConstraint::steal_def( pose::Pose const& pose ) {
	conformation::Conformation const & conformation( pose.conformation() );
	kinematics::Stub my_stub( conformation.stub_from_id( fixed_stub_ ) );
	xyz_target_ = my_stub.global2local( conformation.xyz( atom_ ) );
	//tr.Trace << "get local xyz for " << atom_ << " " << xyz_target_.x() << " " << xyz_target_.y() << " " << xyz_target_.z() << " vs global: " // Vector bla( conformation.xyz( atom_ ) );
	//tr.Trace << bla.x() << " " << bla.y() << " " << bla.z() << " " << std::endl;
}


///
void LocalCoordinateConstraint::show( std::ostream& out ) const
{
	out << "LocalCoordinateConstraint ("
		<< atom_.atomno() << "," << atom_.rsd() << "-"
		<< fixed_stub_ << ")" << std::endl;
	func_->show( out );
}



///
Real
LocalCoordinateConstraint::score(
	Vector const & xyz, //target
	Vector const & s1, //fixed_stub.a
	Vector const & s2, //fixed_stub.b
	Vector const & s3 //fixed_stub.c
) const
{

	kinematics::Stub my_stub( s1, s2, s3 );
	//tr.Trace << "score: global "<< xyz.x() << " " << xyz.y() << " " << xyz.z() << std::endl;
	Vector xyz_local( my_stub.global2local( xyz ) );
	//tr.Trace << xyz_local.x() << " " << xyz_local.y() << " " << xyz_local.z() << " " << std::endl;
	return func( xyz_target_.distance( my_stub.global2local( xyz ) ) );
}

///
void
LocalCoordinateConstraint::score( func::XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const
{
	emap[ this->score_type() ] += score( xyz( atom_ ),
		xyz( fixed_stub_.atom( 1 ) ),
		xyz( fixed_stub_.atom( 2 ) ),
		xyz( fixed_stub_.atom( 3 ) )
	);
}

// atom deriv
void
LocalCoordinateConstraint::fill_f1_f2(
	AtomID const & atom,
	func::XYZ_Func const & xyz,
	Vector & F1,
	Vector & F2,
	EnergyMap const & weights
) const
{
	utility_exit_with_message( " derivative of LocalCoordinateConstraint not supported yet " );
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



///
Size
LocalCoordinateConstraint::natoms() const
{
	return 4;
}

///
id::AtomID const &
LocalCoordinateConstraint::atom( Size const n ) const
{
	if ( n == 1 ) {
		return atom_;
	} else if ( n <= 4 ) {
		return fixed_stub_.atom( n - 1 );
	} else {
		utility_exit_with_message( "LocalCoordinateConstraint::atom() bad argument" );
	}
	return atom_;
}


Real
LocalCoordinateConstraint::dist( pose::Pose const & pose ) const {
	conformation::Conformation const & conformation( pose.conformation() );
	Vector const & xyz( conformation.xyz( atom_ ) );
	kinematics::Stub my_stub( conformation.stub_from_id( fixed_stub_ ) );
	return xyz_target_.distance( my_stub.global2local( xyz ) );
}

Size
LocalCoordinateConstraint::show_violations(
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

void LocalCoordinateConstraint::set_fixed_stub( id::StubID new_stub ) {
	fixed_stub_ = new_stub;
}


Vector LocalCoordinateConstraint::xyz_target( core::pose::Pose const& local_frame_pose ) const {
	tr.Trace << "local target: " <<  xyz_target_.x() << " " << xyz_target_.y() << " " << xyz_target_.z() << std::endl;
	conformation::Conformation const & conformation( local_frame_pose.conformation() );
	kinematics::Stub my_stub( conformation.stub_from_id( fixed_stub_ ) );
	tr.Trace << "stub: " << my_stub << std::endl;
	return my_stub.local2global( xyz_target_ );
}

void LocalCoordinateConstraint::set_xyz_target( Vector const& xyz_in, core::pose::Pose const& local_frame_pose ) {
	conformation::Conformation const & conformation( local_frame_pose.conformation() );
	kinematics::Stub my_stub( conformation.stub_from_id( fixed_stub_ ) );
	xyz_target_=my_stub.global2local( xyz_in );
}

Real
LocalCoordinateConstraint::func( Real const theta ) const
{
	return func_->func( theta );
}

Real
LocalCoordinateConstraint::dfunc( Real const theta ) const
{
	return func_->dfunc( theta );
}

LocalCoordinateConstraint::LocalCoordinateConstraint( LocalCoordinateConstraint const & src ) :
	Constraint( src ),
	atom_( src.atom_ ),
	fixed_stub_( src.fixed_stub_ ),
	xyz_target_( src.xyz_target_ ),
	func_( src.func_ ? src.func_->clone() : src.func_ )
{}

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::LocalCoordinateConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< Constraint >( this ) );
	arc( CEREAL_NVP( atom_ ) ); // id::AtomID
	arc( CEREAL_NVP( fixed_stub_ ) ); // id::StubID
	arc( CEREAL_NVP( xyz_target_ ) ); // Vector
	arc( CEREAL_NVP( func_ ) ); // func::FuncOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::LocalCoordinateConstraint::load( Archive & arc ) {
	arc( cereal::base_class< Constraint >( this ) );
	arc( atom_ ); // id::AtomID
	arc( fixed_stub_ ); // id::StubID
	arc( xyz_target_ ); // Vector
	arc( func_ ); // func::FuncOP
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::LocalCoordinateConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::LocalCoordinateConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_LocalCoordinateConstraint )
#endif // SERIALIZATION
