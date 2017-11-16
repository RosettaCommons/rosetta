// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief Definition of AngularConstraints

// Unit headers
#include <core/scoring/constraints/AngleConstraint.hh>

// Package headers
#include <core/scoring/func/Func.hh>

// Project headers
#include <core/chemical/ResidueType.hh>
#include <core/conformation/Conformation.hh>
#include <core/id/SequenceMapping.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/func/FuncFactory.hh>
#include <core/scoring/func/XYZ_Func.hh>
#include <core/scoring/EnergyMap.hh>
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <basic/Tracer.hh>
#include <core/id/AtomID.hh>
#include <core/id/NamedAtomID.hh>

// Numeric headers
#include <numeric/trig.functions.hh>
#include <numeric/deriv/angle_deriv.hh>

#include <utility/vector1.hh>
#include <utility/assert.hh>

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

static basic::Tracer TR( "core.scoring.constraints.AngleConstraint" );

AngleConstraint::AngleConstraint(
	AtomID const & a1,
	AtomID const & a2,
	AtomID const & a3,
	core::scoring::func::FuncOP func_in, // we take ownership of this guy
	ScoreType scotype /* = angle_constraint */
):
	Constraint( scotype ),
	atom1_(a1),
	atom2_(a2),
	atom3_(a3),
	func_( func_in )
{}

AngleConstraint::AngleConstraint(
	core::scoring::func::FuncOP func_in,
	ScoreType scoretype /* = angle_constraint */
):
	Constraint( scoretype ),
	func_( func_in )
{}

std::string AngleConstraint::type() const {
	return "Angle";
}

ConstraintOP
AngleConstraint::clone() const {
	return ConstraintOP( new AngleConstraint( *this ));
}


bool AngleConstraint::operator == ( Constraint const & other ) const
{
	if ( ! same_type_as_me( other ) ) return false;
	if ( ! other.same_type_as_me( *this ) ) return false;

	AngleConstraint const & other_ang( static_cast< AngleConstraint const & > (other));
	if ( atom1_ != other_ang.atom1_ ) return false;
	if ( atom2_ != other_ang.atom2_ ) return false;
	if ( atom3_ != other_ang.atom3_ ) return false;
	if ( score_type() != other_ang.score_type() ) return false;

	return func_ == other_ang.func_ || ( func_ && other_ang.func_ && *func_ == *other_ang.func_ );
}

bool AngleConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< AngleConstraint const * > (&other);
}

void AngleConstraint::show( std::ostream & out ) const {
	out << type();
	for ( Size i = 1; i <= natoms(); ++i ) {
		AtomID const & id = atom(i);
		out << ' ' << id.rsd() << ' ' << id.atomno();
	}
	out << ' ';
	func_->show_definition(out);
}

void AngleConstraint::show_def( std::ostream& out, pose::Pose const& pose ) const {
	out << type();
	for ( Size i = 1; i <= natoms(); ++i ) {
		AtomID const & id = atom(i);
		out << ' ' <<  atom_id_to_named_atom_id( id, pose );
	}
	out << ' ';
	func_->show_definition( out );
}

/////////////////////////////////////////////////////////////////////////////
/// @details one line definition "Angle atom1 res1 atom2 res2 atom3 res3 function_type function_definition"
///SML: It appears to be reading the angle in radians, because score ultimately uses std:acos, which returns radians.
void
AngleConstraint::read_def(
	std::istream & data,
	pose::Pose const & pose,
	func::FuncFactory const & func_factory
)
{
	Size res1, res2, res3;
	std::string tempres1, tempres2, tempres3;
	std::string name1, name2, name3;
	std::string func_type;

	data
		>> name1 >> tempres1
		>> name2 >> tempres2
		>> name3 >> tempres3
		>> func_type;

	ConstraintIO::parse_residue( pose, tempres1, res1 );
	ConstraintIO::parse_residue( pose, tempres2, res2 );
	ConstraintIO::parse_residue( pose, tempres3, res3 );

	TR.Debug  << "read: " << name1 << " " << name2 << " " << name3 << " "
		<< res1 << " " << res2 << " " << res3 << " func: " << func_type
		<< std::endl;
	if ( res1 > pose.size() || res2 > pose.size() || res3 > pose.size() ) {
		TR.Warning  << "ignored constraint (no such atom in pose!)"
			<< name1 << " " << name2 << " " << res1 << " " << res2 << std::endl;
		data.setstate( std::ios_base::failbit );
		return;
	}

	atom1_ = id::AtomID( pose.residue(res1).atom_index(name1), res1 );
	atom2_ = id::AtomID( pose.residue(res2).atom_index(name2), res2 );
	atom3_ = id::AtomID( pose.residue(res3).atom_index(name3), res3 );
	if ( atom1_.atomno() == 0 || atom2_.atomno() == 0 || atom3_.atomno() == 0 ) {
		TR.Warning << "Error reading atoms: read in atom names("
			<< name1 << "," << name2 << "," << name3 << "), "
			<< "and found AtomIDs (" << atom1_ << "," << atom2_ << "," << atom3_ << ")"
			<< std::endl;
		data.setstate( std::ios_base::failbit );
		return;
	}

	func_ = func_factory.new_func( func_type );
	func_->read_data( data );

	//chu skip the rest of line since this is a single line defintion.
	while ( data.good() && (data.get() != '\n') ) {}

	if ( TR.Debug.visible() ) {
		func_->show_definition( std::cout );
		std::cout << std::endl;
	}

} // read_def

/// @brief Copies the data from this Constraint into a new object and returns an OP
/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
/// to the new object. Intended to be implemented by derived classes.
ConstraintOP
AngleConstraint::remapped_clone( pose::Pose const& src, pose::Pose const& dest, id::SequenceMappingCOP smap ) const {
	id::NamedAtomID atom1( atom_id_to_named_atom_id( atom(1), src ) );
	id::NamedAtomID atom2( atom_id_to_named_atom_id( atom(2), src ) );
	id::NamedAtomID atom3( atom_id_to_named_atom_id( atom(3), src ) );

	if ( smap ) {
		atom1.rsd() = (*smap)[ atom1_.rsd() ];
		atom2.rsd() = (*smap)[ atom2_.rsd() ];
		atom3.rsd() = (*smap)[ atom3_.rsd() ];
	}

	//get AtomIDs for target pose
	id::AtomID id1( named_atom_id_to_atom_id( atom1, dest ) );
	id::AtomID id2( named_atom_id_to_atom_id( atom2, dest ) );
	id::AtomID id3( named_atom_id_to_atom_id( atom3, dest ) );
	if ( id1.valid() && id2.valid() && id3.valid() ) {
		return ConstraintOP( new AngleConstraint( id1, id2, id3, func_, score_type() ) );
	} else {
		return NULL;
	}
}

/////////////////////////////////////////////////////////////////////////////
Real
AngleConstraint::score(
	Vector const & p1,
	Vector const & p2,
	Vector const & p3
) const
{
	return func( angle( p1, p2, p3 ) );
}

void
AngleConstraint::score( func::XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const
{
	emap[ this->score_type() ] += score( xyz( atom1_ ), xyz( atom2_ ), xyz( atom3_ ) );
}

core::Real
AngleConstraint::score( conformation::Conformation const & conf ) const {
	return score( conf.xyz( atom1_ ), conf.xyz( atom2_ ), conf.xyz( atom3_ ) );
}

/// @brief compute
Real
AngleConstraint::angle(
	Vector const & p1,
	Vector const & p2,
	Vector const & p3
) const {
	Vector u1( p1 - p2 );
	Vector u2( p3 - p2 );
	Real const n1( u1.length() );
	Real const n2( u2.length() );
	if ( n1 > 1e-12 && n2 > 1e-12 ) {
		return numeric::arccos( dot( u1,u2 ) / ( n1 * n2 ) );
	}
	TR.Error << "AngleConstraint::score: warning: 0-length bonds!" << std::endl;
	utility_exit_with_message( "0-length bonds in AngleConstraint" );
	return 0.0;
}

core::Real
AngleConstraint::dist( core::scoring::func::XYZ_Func const & xyz ) const {
	return angle( xyz( atom1_ ), xyz( atom2_ ), xyz( atom3_ ) );
}

void
AngleConstraint::setup_for_scoring( func::XYZ_Func const &, ScoreFunction const & ) const {} //Do nothing.

/////////////////////////////////////////////////////////////////////////////
// accumulate F1, F2 contributions from terms that look like
//
//       d(M-F)
// dot( -------- , w )
//       d phi
//
// where phi is moving M and F is fixed.
//
// F1 collects terms f1 that look like: (-u_phi) dot f1
// F2 collects terms f2 that look like: (-u_phi x R_phi) dot f2
//
// where u_phi is the unit vector axis of phi and R_phi is a point
// on the axis
//
// basic id: d(M-F)/dphi = u_phi x ( M - R_phi )
//
// and dot( a, b x c ) = dot( b, c x a )
//
//
void
AngleConstraint::helper(
	Vector const & M,
	Vector const & w,
	Vector & F1,
	Vector & F2
)
{
	F2 += w;
	F1 += cross( w, M );
}

/// @brief evaluate func at theta
Real
AngleConstraint::func( Real const theta ) const
{
	return func_->func( theta );
}

/// @brief evaluate dfunc at theta
Real
AngleConstraint::dfunc( Real const theta ) const
{
	return func_->dfunc( theta );
}

/////////////////////////////////////////////////////////////////////////////
// calculates f1,f2 contributions for dtheta_dphi
//
// where phi is a torsion angle moving p1 while p2 and p3 are fixed
void
AngleConstraint::p1_theta_deriv(
	Vector const & p1,
	Vector const & p2,
	Vector const & p3,
	Vector & F1,
	Vector & F2
)
{
	using numeric::conversions::radians;

	// to avoid problems with dtheta/dx around 0 and 180 degrees
	// truncate x a bit in the calculation of the derivative
	static Real const small_angle( radians( Real(0.1) ) );
	static Real const big_angle( radians( Real(179.9) ) );
	static Real const max_x( std::cos( small_angle ));
	static Real const min_x( std::cos( big_angle ));
	// dtheta_dx has a value of ~ 572.96 for min_x and max_x
	// this goes to infinity as x goes to -1 or 1

	Vector v1( p1 - p2 );
	Vector v2( p3 - p2 );
	Real const n1( v1.length() );
	Real const n2( v2.length() );
	if ( n1 < 1e-9 || n2 < 1e-9 ) {
		return;
	}

	// calculate dx/dphi where x = cos theta = dot( v1,v2) / (n1*n2)

	Real x = dot(v1,v2)/(n1*n2);

	// only v1 depends on phi, not v2 (we are assuming only p1 is moving)
	// so we get two terms, one from v1 on top and one from n1 = |v1| on
	// the bottom

	Vector f1(0.0),f2(0.0);

	{ // first term
		Real const f = 1.0f / ( n1 * n2 );
		helper( p1, f * v2, f1, f2 );
	}

	{ // second term
		Real const f = -1.0f * x / ( n1 * n1 );
		helper( p1, f * v1, f1, f2 );
	}

	x = std::min( std::max( min_x, x ), max_x );
	Real const dtheta_dx = -1.0f / sqrt( 1.0f - x*x );
	f1 *= dtheta_dx;
	f2 *= dtheta_dx;

	// translation of p1 along v1 or perpendicular to v1 and v2 ==> deriv=0
	debug_assert( f1.distance( cross(f2,p1) ) < 1e-3 && // see helper fcn
		std::abs( dot( f2, v1 ) ) < 1e-3 &&
		std::abs( dot( f2, cross( v1, v2 ) ) ) < 1e-3 );


	{ // more debugging
		// pretend axis = u2, R_phi = p2
		ASSERT_ONLY(Vector const u_phi( v2.normalized() );)
		ASSERT_ONLY(Vector const R_phi( p2 );)
		ASSERT_ONLY(Real const deriv = - dot( u_phi, f1 ) - dot( cross( u_phi, R_phi ), f2);)
		debug_assert( std::abs( deriv ) < 1e-3 );
		//std::cout << "deriv: " << deriv<< ' ' <<
		// F(9,3,u_phi(1)) << F(9,3,u_phi(2)) << F(9,3,u_phi(3)) << ' ' <<
		// F(9,3,R_phi(1)) << F(9,3,R_phi(2)) << F(9,3,R_phi(3)) << "\nF1,F2: " <<
		// F(9,3,f1(1)) << F(9,3,f1(2)) << F(9,3,f1(3)) << ' ' <<
		// F(9,3,f2(1)) << F(9,3,f2(2)) << F(9,3,f2(3)) << std::endl;
	}


	F1 += f1;
	F2 += f2;
}

/////////////////////////////////////////////////////////////////////////////
void
AngleConstraint::p1_deriv(
	Vector const & p1,
	Vector const & p2,
	Vector const & p3,
	Vector & F1,
	Vector & F2
) const
{
	Vector u1( p1 - p2 );
	Vector u2( p3 - p2 );
	Real const n1_n2( u1.length() * u2.length() );
	if ( n1_n2 < 1e-12 ) {
		std::cout << "AngleConstraint::p1_deriv: short bonds: " << n1_n2 <<
			std::endl;
		return;
	}

	Vector f1(0.0),f2(0.0);
	p1_theta_deriv( p1, p2, p3, f1, f2 );

	Real d( dot(u1,u2) / n1_n2 );
	Real const tol(0.001);
	if ( d <= -1.0 + tol ) {
		//std::cout << "out-of-tol: " << d << ' ' << std::endl;
		d = -1.0+tol;
	} else if ( d >= 1.0 - tol ) {
		//std::cout << "out-of-tol: " << d << std::endl;
		d = 1.0-tol;
	}
	Real const theta = numeric::arccos( d );


	Real const dE_dtheta = dfunc( theta );

	//std::cout << "dE_dtheta_p1_deriv: " << dE_dtheta << ' ' <<
	// f1(1) << ' ' << f1(2) << ' ' << f1(3) << std::endl;

	F1 += dE_dtheta * f1;
	F2 += dE_dtheta * f2;

}

/////////////////////////////////////////////////////////////////////////////
void
AngleConstraint::p2_deriv(
	Vector const & p1,
	Vector const & p2,
	Vector const & p3,
	Vector & F1,
	Vector & F2
) const {
	Vector v1( p1 - p2 );
	Vector v2( p3 - p2 );
	Real const v12( v1.length() * v2.length() );
	if ( v12 < 1e-12 ) return;

	// here we use the trick that theta = pi - alpha - beta
	//
	// where alpha and beta are the other angles in the p1,p2,p3 triangle
	// so dtheta_dphi  = -dalpha_dphi - dbeta_dphi
	//
	// for these we can use the p1_theta_deriv formula
	//
	Vector f1(0.0),f2(0.0);
	p1_theta_deriv( p2, p1, p3, f1, f2 ); // alpha deriv
	p1_theta_deriv( p2, p3, p1, f1, f2 ); // beta deriv

	// translation of p2 atom perpendicular to plane ==> deriv = 0
	//std::cout << "p2 deriv check: " << std::endl;
	debug_assert( std::abs( dot( f2, cross( v1,v2) ) ) < 1e-3 );

	Real d( dot(v1,v2) / v12 );
	Real const tol(0.001);
	if ( d <= -1.0 + tol ) {
		//std::cout << "out-of-tol: " << d << ' ' << std::endl;
		d = -1.0+tol;
	} else if ( d >= 1.0 - tol ) {
		//std::cout << "out-of-tol: " << d << std::endl;
		d = 1.0-tol;
	}
	Real const theta = numeric::arccos( d );
	//Real const theta = arccos( dot( v1,v2 ) / v12 );
	Real const dE_dtheta = dfunc( theta );


	//std::cout << "dE_dtheta_p2_deriv: " << dE_dtheta << ' ' <<
	// f1(1) << ' ' << f1(2) << ' ' << f1(3) << std::endl;

	F1 += -1.0f * dE_dtheta * f1;
	F2 += -1.0f * dE_dtheta * f2;
}


/////////////////////////////////////////////////////////////////////////////
void
AngleConstraint::fill_f1_f2(
	AtomID const & atom,
	func::XYZ_Func const & xyz,
	Vector & F1,
	Vector & F2,
	EnergyMap const & weights
) const {

	Vector unweighted_F1(0.);
	Vector unweighted_F2(0.);
	using numeric::deriv::angle_p1_deriv;
	using numeric::deriv::angle_p2_deriv;
	Real theta( 0.0 );
	if ( atom == atom1_ ) {
		angle_p1_deriv( xyz( atom1_ ), xyz( atom2_ ), xyz( atom3_ ), theta, unweighted_F1, unweighted_F2 );
	} else if ( atom == atom2_ ) {
		angle_p2_deriv( xyz( atom1_ ), xyz( atom2_ ), xyz( atom3_ ), theta, unweighted_F1, unweighted_F2 );
	} else if ( atom == atom3_ ) {
		angle_p1_deriv( xyz( atom3_ ), xyz( atom2_ ), xyz( atom1_ ), theta, unweighted_F1, unweighted_F2 );
	} else {
		return;
	}

	Real dfunc_dtheta( dfunc( theta ));

	F1 += weights[ this->score_type() ] * dfunc_dtheta * unweighted_F1;
	F2 += weights[ this->score_type() ] * dfunc_dtheta * unweighted_F2;
}

ConstraintOP
AngleConstraint::remap_resid(
	core::id::SequenceMapping const & seqmap
) const {
	if ( seqmap[atom1_.rsd()] != 0 && seqmap[atom2_.rsd()] != 0 && seqmap[atom3_.rsd()] != 0 ) {
		AtomID remap_a1( atom1_.atomno(), seqmap[atom1_.rsd()] ),
			remap_a2( atom2_.atomno(), seqmap[atom2_.rsd()] ),
			remap_a3( atom3_.atomno(), seqmap[atom3_.rsd()] );
		return ConstraintOP( new AngleConstraint( remap_a1, remap_a2, remap_a3, this->func_ ) );
	} else {
		return NULL;
	}
}

id::AtomID const &
AngleConstraint::atom( Size const n ) const
{
	switch( n ) {
	case 1 :
		return atom1_;
	case 2 :
		return atom2_;
	case 3 :
		return atom3_;
	default :
		utility_exit_with_message( "AngleConstraint::atom() bad argument" );
	}
	return atom1_;
}

/// @details show violations of angular constraint. Control verbosity level between 1 .. 100
/// >80 : write MET 4 CA
/// otherwise: compute angle and call func_->show_violations
Size AngleConstraint::show_violations(
	std::ostream & out,
	pose::Pose const & pose,
	Size verbose_level,
	Real threshold
) const {
	conformation::Conformation const & conformation( pose.conformation() );
	if ( verbose_level > 80 ) {
		out << "AngleConstraint ("
			<< pose.residue_type(atom1_.rsd() ).atom_name( atom1_.atomno() ) << ":" << atom1_.atomno() << "," << atom1_.rsd() << "-"
			<< pose.residue_type(atom2_.rsd() ).atom_name( atom2_.atomno() ) << ":" << atom2_.atomno() << "," << atom2_.rsd() << "-"
			<< pose.residue_type(atom3_.rsd() ).atom_name( atom3_.atomno() ) << ":" << atom3_.atomno() << "," << atom3_.rsd() << ") ";
	};

	// compute angle
	Vector const & p1( conformation.xyz( atom1_ ) ), p2( conformation.xyz( atom2_ ) ), p3( conformation.xyz( atom3_ ) );
	Vector u1( p1 - p2 );
	Vector u2( p3 - p2 );
	Real const n1( u1.length() );
	Real const n2( u2.length() );
	if ( n1 > 1e-12 && n2 > 1e-12 ) {
		Real const theta = numeric::arccos( dot( u1,u2 ) / ( n1 * n2 ) );

		// angle is meaningful, call show_violation of func_
		return func_->show_violations( out, theta, verbose_level, threshold );
	};
	std::cout << "AngleConstraint::show_violations: error: 0-length bonds!"
		<< std::endl;
	return 0;
}

AngleConstraint::AngleConstraint( AngleConstraint const & src ) :
	Constraint( src ),
	atom1_( src.atom1_ ),
	atom2_( src.atom2_ ),
	atom3_( src.atom3_ ),
	func_( src.func_->clone() )
{}

/// @brief const access to func
func::FuncCOP AngleConstraint::func() const { return func_; }

/// @brief set func
void AngleConstraint::set_func( func::FuncOP f ) { func_ = f; }

void
AngleConstraint::atom1( AtomID setting ) const {
	atom1_ = setting;
}

void
AngleConstraint::atom2( AtomID setting ) const {
	atom2_ = setting;
}

void
AngleConstraint::atom3( AtomID setting ) const {
	atom3_ = setting;
}

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION

/// @brief Default constructor required by cereal to deserialize this class
core::scoring::constraints::AngleConstraint::AngleConstraint() {}

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::AngleConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< Constraint >( this ) );
	arc( CEREAL_NVP( atom1_ ) ); // AtomID
	arc( CEREAL_NVP( atom2_ ) ); // AtomID
	arc( CEREAL_NVP( atom3_ ) ); // AtomID
	arc( CEREAL_NVP( func_ ) ); // core::scoring::func::FuncOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::AngleConstraint::load( Archive & arc ) {
	arc( cereal::base_class< Constraint >( this ) );
	arc( atom1_ ); // AtomID
	arc( atom2_ ); // AtomID
	arc( atom3_ ); // AtomID
	arc( func_ ); // core::scoring::func::FuncOP
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::AngleConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::AngleConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_AngleConstraint )
#endif // SERIALIZATION
