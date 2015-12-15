// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file src/core/scoring/DistancePairConstraint.cc
/// @brief Restrain a pair of residues to take the same AtomPair constraint of another pair
/// @author Frank DiMaio, Fabio Parmeggiani

#include <core/scoring/constraints/DistancePairConstraint.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/FuncFactory.hh>
#include <core/scoring/func/Func.hh>

#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/conformation/Conformation.hh>
#include <basic/Tracer.hh>

#include <numeric/xyz.functions.hh>
#include <numeric/trig.functions.hh>
#include <numeric/deriv/distance_deriv.hh>

#include <utility/exit.hh>

#include <core/id/NamedAtomID.hh>
#include <core/id/SequenceMapping.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <utility/vector1.hh>

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

static THREAD_LOCAL basic::Tracer TR( "core.io.constraints" );


/////////////////////////////////////////////////////////////////////////////

ConstraintOP
DistancePairConstraint::remap_resid(
	core::id::SequenceMapping const & seqmap
) const {
	if (   seqmap[atomA1_.rsd()] != 0 && seqmap[atomA2_.rsd()] != 0
			&& seqmap[atomB1_.rsd()] != 0 && seqmap[atomB2_.rsd()] != 0 ) {
		AtomID remap_a1( atomA1_.atomno(), seqmap[atomA1_.rsd()] ),
			remap_a2( atomA2_.atomno(), seqmap[atomA2_.rsd()] );
		AtomID remap_b1( atomB1_.atomno(), seqmap[atomB1_.rsd()] ),
			remap_b2( atomB2_.atomno(), seqmap[atomB2_.rsd()] );
		return ConstraintOP( new DistancePairConstraint(
			remap_a1, remap_a2,
			remap_b1, remap_b2,
			this->func_ ) );
	} else {
		return NULL;
	}
}


/// @brief Copies the data from this Constraint into a new object and returns an OP
/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
/// to the new object. Intended to be implemented by derived classes.
ConstraintOP DistancePairConstraint::remapped_clone( pose::Pose const& src, pose::Pose const& dest, id::SequenceMappingCOP smap ) const {
	id::NamedAtomID atomA1( core::pose::atom_id_to_named_atom_id(atom(1), src ) );
	id::NamedAtomID atomA2( core::pose::atom_id_to_named_atom_id(atom(2), src ) );
	id::NamedAtomID atomB1( core::pose::atom_id_to_named_atom_id(atom(3), src ) );
	id::NamedAtomID atomB2( core::pose::atom_id_to_named_atom_id(atom(4), src ) );
	if ( smap ) {
		atomA1.rsd() = (*smap)[ atomA1_.rsd() ];
		atomA2.rsd() = (*smap)[ atomA2_.rsd() ];
		atomB1.rsd() = (*smap)[ atomB1_.rsd() ];
		atomB2.rsd() = (*smap)[ atomB2_.rsd() ];
	}

	//get AtomIDs for target pose
	id::AtomID id1( core::pose::named_atom_id_to_atom_id(atomA1, dest ));
	id::AtomID id2( core::pose::named_atom_id_to_atom_id(atomA2, dest ));
	id::AtomID id3( core::pose::named_atom_id_to_atom_id(atomB1, dest ));
	id::AtomID id4( core::pose::named_atom_id_to_atom_id(atomB2, dest ));
	if (    id1.valid() && id2.valid() &&  id3.valid() && id4.valid() ) {
		return ConstraintOP( new DistancePairConstraint( id1, id2, id3, id4, func_, score_type() ) );
	} else {
		return NULL;
	}
}


id::AtomID const &
DistancePairConstraint::atom( Size const n ) const {
	switch( n ) {
	case 1 : return atomA1_;
	case 2 : return atomA2_;
	case 3 : return atomB1_;
	case 4 : return atomB2_;
	default :
		utility_exit_with_message( "DistancePairConstraint::atom() bad argument" );
	}
	return atomA1_;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/// @details one line definition "DistancePair atom1 res1 atom2 res2 atom3 res3 atom4 res4 function_type function_definition"
void
DistancePairConstraint::read_def(
	std::istream & in,
	pose::Pose const & pose,
	func::FuncFactory const & func_factory
) {
	Size res1, res2, res3, res4;
	std::string tempres1, tempres2, tempres3, tempres4;
	std::string name1, name2, name3, name4;
	std::string func_type;

	in
		>> name1 >> tempres1
		>> name2 >> tempres2
		>> name3 >> tempres3
		>> name4 >> tempres4
		>> func_type;

	ConstraintIO::parse_residue( pose, tempres1, res1 );
	ConstraintIO::parse_residue( pose, tempres2, res2 );
	ConstraintIO::parse_residue( pose, tempres3, res3 );
	ConstraintIO::parse_residue( pose, tempres4, res4 );

	TR.Debug  << "read: " << name1 << " " << name2 << " "
		<< res1 << " " << res2 << " func: " << func_type
		<< std::endl;
	if (    res1 > pose.total_residue() || res2 > pose.total_residue()
			|| res3 > pose.total_residue() || res4 > pose.total_residue()  ) {
		TR.Warning  << "ignored constraint (no such atom in pose!)"
			<< name1 << " " << name2 << " "
			<< res1 << " " << res2 << " func: " << func_type
			<< std::endl;
		in.setstate( std::ios_base::failbit );
		return;
	}

	atomA1_ = id::AtomID( core::pose::named_atom_id_to_atom_id( id::NamedAtomID( name1, res1 ), pose ) );
	atomA2_ = id::AtomID( core::pose::named_atom_id_to_atom_id( id::NamedAtomID( name2, res2 ), pose ) );
	atomB1_ = id::AtomID( core::pose::named_atom_id_to_atom_id( id::NamedAtomID( name3, res3 ), pose ) );
	atomB2_ = id::AtomID( core::pose::named_atom_id_to_atom_id( id::NamedAtomID( name4, res4 ), pose ) );


	if (    atomA1_.atomno() == 0 || atomA2_.atomno() == 0
			|| atomB1_.atomno() == 0 || atomB2_.atomno() == 0 ) {
		TR.Warning << "Error reading atoms: read in atom names("
			<< name1 << "," << name2 << "," << name3 << "," << name4 << "), "
			<< "and found AtomIDs (" << atomA1_ << "," << atomA2_ << " -- "
			<< atomB1_ << "," << atomB2_ << ")" << std::endl;
		in.setstate( std::ios_base::failbit );
		return;
	}

	func_ = func_factory.new_func( func_type );
	func_->read_data( in );

	while ( in.good() && (in.get() != '\n') ) {}

	if ( TR.Debug.visible() ) {
		func_->show_definition( std::cout );
		std::cout << std::endl;
	}
}

/////////////////////////////////////////////////////////////////////////////

Real
DistancePairConstraint::score(
	conformation::Conformation const & conformation
) const {
	return score(
		conformation.xyz( atomA1_ ), conformation.xyz( atomA2_ ),
		conformation.xyz( atomB1_ ), conformation.xyz( atomB2_ ));
}

void
DistancePairConstraint::score(
	func::XYZ_Func const & xyz,
	EnergyMap const &,
	EnergyMap & emap
) const {
	emap[ this->score_type() ] += score(
		xyz( atomA1_ ), xyz( atomA2_ ),
		xyz( atomB1_ ), xyz( atomB2_ )
	);
}


Real
DistancePairConstraint::score(
	Vector const & p1, Vector const & p2, Vector const & p3, Vector const & p4
) const {
	core::Real dis1 = p1.distance(p2) ;
	core::Real dis2 = p3.distance(p4) ;
	core::Real difference = dis2 - dis1 ;

	TR.Debug << "difference of " << difference << " Angstroms" << std::endl;
	return func(difference);
}


void
DistancePairConstraint::fill_f1_f2(
	AtomID const & /* atom */,
	func::XYZ_Func const & xyz,
	Vector & F1,
	Vector & F2,
	EnergyMap const & weights
) const {
	using namespace numeric::deriv;

	Vector f1(0.0) ,f2(0.0);

	Real dist(0.0), dist0(0.0);

	numeric::deriv::distance_f1_f2_deriv( xyz( atomA1_ ), xyz( atomA2_ ), dist0, f1, f2 );
	numeric::deriv::distance_f1_f2_deriv( xyz( atomB1_ ), xyz( atomB2_ ), dist, f1, f2 );

	core::Real difference = dist - dist0;
	Real const dE_dist( dfunc( difference ) );

	F1 += dE_dist * weights[ this->score_type() ] * f1;
	F2 += dE_dist * weights[ this->score_type() ] * f2;
}

void DistancePairConstraint::show( std::ostream & out ) const {
	out << "DistancePairConstraint";
	for ( Size i = 1; i <= natoms(); ++i ) {
		AtomID const & id = atom(i);
		out << ' ' << id.rsd() << ' ' << id.atomno();
	}
	out << ' ';
	func_->show_definition(out);
}

std::string DistancePairConstraint::type() const {
	return "DistancePair";
}

ConstraintOP DistancePairConstraint::clone() const {
	return ConstraintOP( new DistancePairConstraint( *this ) );
}

Size DistancePairConstraint::show_violations(
	std::ostream& out,
	pose::Pose const& pose,
	Size verbose_level,
	Real threshold
) const {
	if ( verbose_level > 80 ) {
		out << "DistancePair ("
			<< pose.residue_type(atomA1_.rsd() ).atom_name( atomA1_.atomno() ) << ":"
			<< atomA1_.atomno() << "," << atomA1_.rsd() << "-"
			<< pose.residue_type(atomA2_.rsd() ).atom_name( atomA2_.atomno() ) << ":"
			<< atomA2_.atomno() << "," << atomA2_.rsd() << " === "
			<< pose.residue_type(atomB1_.rsd() ).atom_name( atomB1_.atomno() ) << ":"
			<< atomB1_.atomno() << "," << atomB1_.rsd() << "-"
			<< pose.residue_type(atomB2_.rsd() ).atom_name( atomB2_.atomno() ) << ":"
			<< atomB2_.atomno() << "," << atomB2_.rsd() << "-";
	}
	return func_->show_violations( out, 0.0, verbose_level, threshold );
}

bool DistancePairConstraint::operator == ( Constraint const & rhs ) const
{
	if ( !     same_type_as_me(   rhs ) ) return false;
	if ( ! rhs.same_type_as_me( *this ) ) return false;

	DistancePairConstraint const & rhs_dpc( static_cast< DistancePairConstraint const & > ( rhs ));
	if ( atomA1_      != rhs_dpc.atomA1_ ) return false;
	if ( atomA2_      != rhs_dpc.atomA2_ ) return false;
	if ( atomB1_      != rhs_dpc.atomB1_ ) return false;
	if ( atomB2_      != rhs_dpc.atomB2_ ) return false;
	if ( score_type() != rhs_dpc.score_type() ) return false;

	return func_ == rhs_dpc.func_ || ( func_ && rhs_dpc.func_ && *func_ == *rhs_dpc.func_ );
}

bool DistancePairConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< DistancePairConstraint const * > ( &other );
}

Real
DistancePairConstraint::func( Real const theta ) const {
	return func_->func( theta );
}

Real
DistancePairConstraint::dfunc( Real const theta ) const {
	return func_->dfunc( theta );
}

DistancePairConstraint::DistancePairConstraint( DistancePairConstraint const & src ) :
	Constraint( src ),
	atomA1_( src.atomA1_ ),
	atomA2_( src.atomA2_ ),
	atomB1_( src.atomB1_ ),
	atomB2_( src.atomB2_ ),
	func_( src.func_ ? src.func_->clone() : src.func_ )
{}

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::constraints::DistancePairConstraint::DistancePairConstraint() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::DistancePairConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< Constraint >( this ) );
	arc( CEREAL_NVP( atomA1_ ) ); // AtomID
	arc( CEREAL_NVP( atomA2_ ) ); // AtomID
	arc( CEREAL_NVP( atomB1_ ) ); // AtomID
	arc( CEREAL_NVP( atomB2_ ) ); // AtomID
	arc( CEREAL_NVP( func_ ) ); // func::FuncOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::DistancePairConstraint::load( Archive & arc ) {
	arc( cereal::base_class< Constraint >( this ) );
	arc( atomA1_ ); // AtomID
	arc( atomA2_ ); // AtomID
	arc( atomB1_ ); // AtomID
	arc( atomB2_ ); // AtomID
	arc( func_ ); // func::FuncOP
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::DistancePairConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::DistancePairConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_DistancePairConstraint )
#endif // SERIALIZATION
