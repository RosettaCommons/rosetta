// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/AtomPairConstraint.cc
///
/// @brief
/// @author Oliver Lange

// Unit Headers
#include <core/scoring/constraints/AtomPairConstraint.hh>

// Package Headers
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/Func.hh>
#include <core/scoring/func/FuncFactory.hh>

// Project Headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/NamedAtomID.hh>
#include <core/id/AtomID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/kinematics/ShortestPathInFoldTree.hh>
// Utility Headers
#include <basic/Tracer.hh>

// Numeric Headers
#include <numeric/deriv/distance_deriv.hh>

//Auto Headers
#include <core/id/SequenceMapping.hh>
#include <core/kinematics/AtomTree.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <utility/vector1.hh>


static THREAD_LOCAL basic::Tracer tr( "core.io.constraints" );


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

using core::pose::atom_id_to_named_atom_id;
using core::pose::named_atom_id_to_atom_id;

// default c-tor
AtomPairConstraint::AtomPairConstraint() : Constraint( atom_pair_constraint ) {}

///c-tor
AtomPairConstraint::AtomPairConstraint(
	AtomID const & a1,
	AtomID const & a2,
	core::scoring::func::FuncOP func,
	ScoreType scoretype /* = atom_pair_constraint*/
):
	Constraint( scoretype ),
	atom1_(a1),
	atom2_(a2),
	func_( func )
{}


ConstraintOP
AtomPairConstraint::clone() const {
	return ConstraintOP( new AtomPairConstraint( *this ));
}

ConstraintOP
AtomPairConstraint::clone( func::FuncOP func ) const {
	return ConstraintOP( new AtomPairConstraint( atom1_, atom2_, func, score_type() ) );
}

///
void
AtomPairConstraint::score( func::XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const
{
	Real const score_val =  score( xyz( atom1_ ), xyz( atom2_ ) );
	emap[ this->score_type() ] += score_val;
}

Real AtomPairConstraint::score( pose::Pose const& pose ) const {
	return func_->func( dist( pose ) );
}


Real
AtomPairConstraint::score(
	Vector const & xyz1,
	Vector const & xyz2
) const
{
	// std::cout << "score " <<  atom1_ << ' ' << atom2_ << ' ' << xyz1.distance( xyz2 ) << " --> " << func( xyz1.distance( xyz2 ) ) << std::endl;;
	return func( xyz1.distance( xyz2 ) );
}

/// @brief Copies the data from this Constraint into a new object and returns an OP
/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
/// to the new object. Intended to be implemented by derived classes.
ConstraintOP AtomPairConstraint::remapped_clone( pose::Pose const& src, pose::Pose const& dest, id::SequenceMappingCOP smap ) const {
	id::NamedAtomID atom1( atom_id_to_named_atom_id( atom(1), src ) );
	id::NamedAtomID atom2( atom_id_to_named_atom_id( atom(2), src ) );

	if ( smap ) {
		atom1.rsd() = (*smap)[ atom1_.rsd() ];
		atom2.rsd() = (*smap)[ atom2_.rsd() ];
	}

	if ( atom1.rsd() == 0 || atom2.rsd() == 0 ) return NULL;

	//get AtomIDs for target pose
	id::AtomID id1( named_atom_id_to_atom_id( atom1, dest ) );
	id::AtomID id2( named_atom_id_to_atom_id( atom2, dest ) );
	if ( id1.valid() && id2.valid() ) {
		return ConstraintOP( new AtomPairConstraint( id1, id2, func_, score_type() ) );
	} else {
		return NULL;
	}
}

bool
AtomPairConstraint::operator == ( Constraint const & other_cst ) const
{
	if ( !           same_type_as_me( other_cst ) ) return false;
	if ( ! other_cst.same_type_as_me(     *this ) ) return false;

	AtomPairConstraint const & other( static_cast< AtomPairConstraint const & > (other_cst) );
	if ( atom1_ != other.atom1_ ) return false;
	if ( atom2_ != other.atom2_ ) return false;
	if ( this->score_type() != other.score_type() ) return false;

	return func_ == other.func_ || ( func_ && other.func_ && *func_ == *other.func_ );
}

bool
AtomPairConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< AtomPairConstraint const * > (&other);
}

void AtomPairConstraint::show( std::ostream& out ) const {
	out << "AtomPairConstraint ("
		<< atom1_.atomno() << "," << atom1_.rsd() << "-"
		<< atom2_.atomno() << "," << atom2_.rsd() << ")" << std::endl;
	func_->show( out );
}

void AtomPairConstraint::show_def( std::ostream& out, pose::Pose const& pose ) const {
	out << type() << " " << atom_id_to_named_atom_id( atom1_, pose ) << " " << atom_id_to_named_atom_id( atom2_, pose ) << " ";
	func_->show_definition( out );
}

Real
AtomPairConstraint::dist( pose::Pose const & pose ) const {
	return dist( pose.conformation() );
}

Real AtomPairConstraint::dist( conformation::Conformation  const& conformation ) const {
	conformation::Residue const& res1( conformation.residue( atom1_.rsd() ) );
	conformation::Residue const& res2( conformation.residue( atom2_.rsd() ) );
#ifndef NDEBUG
	bool fail( false );
#endif
	if ( atom1_.atomno() == 0 || atom1_.atomno() > res1.natoms() ) {
		std::cerr << "AtomPairConstraint: atom1 out of bounds" << atom1_ << std::endl;
#ifndef NDEBUG
		fail = true;
#endif
	}
	if ( atom2_.atomno() == 0 || atom2_.atomno() > res2.natoms() ) {
		std::cerr << "AtomPairConstraint: atom2 out of bounds" << atom2_ << std::endl;
#ifndef NDEBUG
		fail = true;
#endif
	}
	debug_assert( conformation.atom_tree().has( atom1_ ) );
	debug_assert( conformation.atom_tree().has( atom2_ ) );
	if ( !conformation.atom_tree().has( atom1_ ) ) {
		std::cerr << "AtomPairConstraint: cannot find atom " << atom1_ << std::endl;
	}
	if ( !conformation.atom_tree().has( atom2_ ) ) {
		std::cerr << "AtomPairConstraint: cannot find atom " << atom2_ << std::endl;
	}
	debug_assert( !fail );
	Vector const & xyz1( conformation.xyz( atom1_ ) ), xyz2( conformation.xyz( atom2_ ) );
	Vector const f2( xyz1 - xyz2 );
	Real const dist( f2.length() );
	return dist;
}

Real AtomPairConstraint::dist( func::XYZ_Func const & xyz ) const {
	debug_assert( atom1_.atomno() );
	debug_assert( atom2_.atomno() );
	return xyz( atom1_ ).distance( xyz( atom2_ ) );
}

func::Func const & AtomPairConstraint::get_func() const {
	return *func_;
}

Size AtomPairConstraint::show_violations(
	std::ostream& out,
	pose::Pose const& pose,
	Size verbose_level,
	Real threshold
) const {

	if ( verbose_level > 80 ) {
		out << "AtomPairConstraint ("
			<< pose.residue_type(atom1_.rsd() ).atom_name( atom1_.atomno() ) << ":"
			<< atom1_.atomno() << "," << atom1_.rsd() << "-"
			<< pose.residue_type(atom2_.rsd() ).atom_name( atom2_.atomno() ) << ":"
			<< atom2_.atomno() << "," << atom2_.rsd() << ") ";
	}
	if ( verbose_level > 120 ) { //don't ask but I had a really weird bug to track down!
		conformation::Conformation const & conformation( pose.conformation() );
		Vector const & xyz1( conformation.xyz( atom1_ ) ), xyz2( conformation.xyz( atom2_ ) );
		out << "\ncoords1: " << xyz1[ 1 ] << " " << xyz1[ 2 ] << " " << xyz1[ 3 ] << " --- ";
		out << "coords1: " << xyz2[ 1 ] << " " << xyz2[ 2 ] << " " << xyz2[ 3 ] << "\n";
	}

	return func_->show_violations( out, dist( pose ), verbose_level, threshold );
}

// atom deriv
void
AtomPairConstraint::fill_f1_f2(
	AtomID const & atom,
	func::XYZ_Func const & xyz,
	Vector & F1,
	Vector & F2,
	EnergyMap const & weights
) const
{
	AtomID other_atom;
	if ( atom == atom1_ ) other_atom = atom2_;
	else if ( atom == atom2_ ) other_atom = atom1_;
	else {
		// std::cout << "Error in AtomPairConstraint::fill_f1_f2()" << std::endl;
		// std::cout << "Bad AtomID: (" << atom.rsd() << ", " << atom.atomno() << ") -- options: (";
		// std::cout << atom1_.rsd() << ", " << atom1_.atomno() << ") and (";
		// std::cout << atom2_.rsd() << ", " << atom2_.atomno() << ");" << std::endl;
		return;
	}

	//Vector const & xyz1( conformation.xyz( atom ) ), xyz2( conformation.xyz( other_atom ) );

	/*
	Vector const f2( xyz1 - xyz2 );
	Real const dist( f2.length() ), deriv( dfunc( dist ) );
	if ( deriv != 0.0 && dist != 0.0 ) {
	Vector const f1( xyz1.cross( xyz2 ) );
	F1 += ( ( deriv / dist ) * f1 ) * weights[ this->score_type() ];
	F2 += ( ( deriv / dist ) * f2 ) * weights[ this->score_type() ];
	}*/

	Real dist(0.0);
	Vector f1(0.0), f2(0.0);
	numeric::deriv::distance_f1_f2_deriv( xyz( atom ), xyz( other_atom ), dist, f1, f2 );
	Real wderiv( weights[ this->score_type() ] * dfunc( dist ));
	F1 += wderiv * f1;
	F2 += wderiv * f2;

}

ConstraintOP
AtomPairConstraint::remap_resid( core::id::SequenceMapping const &seqmap ) const
{
	if ( seqmap[atom1_.rsd()] != 0 && seqmap[atom2_.rsd()] != 0 ) {
		AtomID remap_a1( atom1_.atomno(), seqmap[atom1_.rsd()] ),
			remap_a2( atom2_.atomno(), seqmap[atom2_.rsd()] );
		return ConstraintOP( new AtomPairConstraint( remap_a1, remap_a2, this->func_ ) );
	} else {
		return NULL;
	}
}

/// @details one line definition "AtomPairs atom1 res1 atom2 res2 function_type function_definition"
void
AtomPairConstraint::read_def(
	std::istream & data,
	core::pose::Pose const & pose,
	func::FuncFactory const & func_factory
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
	if ( res1 > pose.total_residue() || res2 > pose.total_residue() ) {
		tr.Warning  << "ignored constraint (requested residue numbers exceed numbers of residues in pose): " << "Total in Pose: " << pose.total_residue() << " "
			<< name1 << " " << name2 << " " << res1 << " " << res2 << std::endl;
		data.setstate( std::ios_base::failbit );
		return;
	}

	if ( name1 == "H" && res1 == 1 && pose.is_fullatom() ) name1 = "1H";
	if ( name2 == "H" && res2 == 1 && pose.is_fullatom() ) name2 = "1H";

	atom1_ = named_atom_id_to_atom_id( id::NamedAtomID( name1, res1 ), pose );
	atom2_ = named_atom_id_to_atom_id( id::NamedAtomID( name2, res2 ), pose );

	if ( atom1_.atomno() == 0 || atom2_.atomno() == 0 ) {
		tr.Warning << "Error reading atoms: read in atom names("
			<< name1 << "," << name2 << "), "
			<< "and found AtomIDs (" << atom1_ << "," << atom2_ << ")" << std::endl;
		data.setstate( std::ios_base::failbit );
		runtime_assert( false );
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
		func_->show_definition( tr.Debug );
		tr.Debug << std::endl;
	}
	runtime_assert( atom1_.valid() && atom2_.valid() );
} // read_def

core::Size
AtomPairConstraint::effective_sequence_separation( core::kinematics::ShortestPathInFoldTree const& sp ) const {
	return sp.dist( atom1_.rsd(), atom2_.rsd() );
}

// functions
Real
AtomPairConstraint::func( Real const theta ) const
{
	return func_->func( theta );
}

// deriv
Real
AtomPairConstraint::dfunc( Real const theta ) const
{
	return func_->dfunc( theta );
}

void
AtomPairConstraint::set_func( func::FuncOP setting ) {
	func_ = setting;
}

void
AtomPairConstraint::atom1( AtomID setting ) const {
	atom1_ = setting;
}

void
AtomPairConstraint::atom2( AtomID setting ) const {
	atom2_ = setting;
}

AtomPairConstraint::AtomPairConstraint( AtomPairConstraint const & src ) :
	Constraint( src ),
	atom1_( src.atom1_ ),
	atom2_( src.atom2_ ),
	func_( src.func_ ? src.func_->clone() : src.func_ )
{}

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::AtomPairConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< Constraint >( this ) );
	arc( CEREAL_NVP( atom1_ ) ); // AtomID
	arc( CEREAL_NVP( atom2_ ) ); // AtomID
	arc( CEREAL_NVP( func_ ) ); // core::scoring::func::FuncOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::AtomPairConstraint::load( Archive & arc ) {
	arc( cereal::base_class< Constraint >( this ) );
	arc( atom1_ ); // AtomID
	arc( atom2_ ); // AtomID
	arc( func_ ); // core::scoring::func::FuncOP
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::AtomPairConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::AtomPairConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_AtomPairConstraint )
#endif // SERIALIZATION
