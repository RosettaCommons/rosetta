// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/scoring/constraints/AtomToAxisConstraint.cc
/// @author Jack Maguire, jackmaguire1444@gmail.com

// Unit Headers
#include <core/scoring/constraints/AtomToAxisConstraint.hh>

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
#include <utility>
#include <utility/vector1.hh>


static basic::Tracer TR( "core.scoring.constraints.AtomToAxisConstraint" );


#ifdef SERIALIZATION
// Utility serialization headers
#include <utility/serialization/serialization.hh>

// Cereal headers
#include <cereal/types/polymorphic.hpp>
#endif // SERIALIZATION


namespace core {
namespace scoring {
namespace constraints {

/// @details This can be sped up by using distances-squared
core::Real
distance_from_line(
	XYZ const & p,
	XYZ const & x1,
	XYZ const & x2
) {
	//https://mathworld.wolfram.com/Point-LineDistance3-Dimensional.html
	// d = |(x2-x1)*(x1-p)| / |x2-x1|
	// where
	//   x1 and x2 are on the line
	//   p is the point
	//   * is a cross product

	XYZ const x2mx1 = x2 - x1;
	XYZ const x1mp = x1 - p;
	core::Real const numerator = x2mx1.cross( x1mp ).length();
	core::Real const denominator = x2mx1.length();
	return numerator / denominator;
}

XYZ
closest_point_on_line(
	XYZ const & p,
	XYZ const & x1,
	XYZ const & x2
) {
	//https://stackoverflow.com/questions/47481774/getting-point-on-line-segment-that-is-closest-to-another-point
	XYZ const x2mx1 = x2 - x1; //TODO cache this
	XYZ const x1mp = p - x1;
	core::Real const x2mx1_l2 = x2mx1.length_squared();
	core::Real const interpolation = x2mx1.dot( x1mp ) / x2mx1_l2;
	return x1 + interpolation * x2mx1;
}

XYZ
closest_point_on_line(
	XYZs const & xyzs
) {
	return closest_point_on_line( xyzs.atom, xyzs.axis1, xyzs.axis2 );
}

using core::pose::atom_id_to_named_atom_id;
using core::pose::named_atom_id_to_atom_id;

// default c-tor
AtomToAxisConstraint::AtomToAxisConstraint() : Constraint( atom_pair_constraint ) {}

///c-tor
AtomToAxisConstraint::AtomToAxisConstraint(
	AtomID const & atom,
	utility::vector1< AtomID > const & axis1,
	utility::vector1< AtomID > const & axis2,
	core::scoring::func::FuncOP const & func,
	ScoreType scoretype /* = atom_pair_constraint*/
):
	Constraint( scoretype ),
	atom_(atom),
	axis1_(axis1),
	axis2_(axis2),
	func_( func )
{}


ConstraintOP
AtomToAxisConstraint::clone() const {
	return utility::pointer::make_shared< AtomToAxisConstraint >( *this );
}

ConstraintOP
AtomToAxisConstraint::clone( func::FuncOP func ) const {
	return utility::pointer::make_shared< AtomToAxisConstraint >( atom_, axis1_, axis2_, func, score_type() );
}

///
void
AtomToAxisConstraint::score( func::XYZ_Func const & xyz, EnergyMap const &, EnergyMap & emap ) const
{
	XYZs const xyzs = calc_xyzs( xyz );
	XYZ const closest = closest_point_on_line( xyzs );
	Real const score_val = score( xyz( atom_ ), closest );
	emap[ this->score_type() ] += score_val;
}

Real
AtomToAxisConstraint::score( pose::Pose const& pose ) const {
	return func_->func( dist( pose ) );
}

Real
AtomToAxisConstraint::score(
	Vector const & xyz1,
	Vector const & xyz2
) const
{
	return func( xyz1.distance( xyz2 ) );
}

XYZs
AtomToAxisConstraint::calc_xyzs(
	core::pose::Pose const & pose
) const {
	XYZs result;
	result.atom = pose.xyz( atom_ );

	debug_assert( ! axis1_.empty() );
	debug_assert( ! axis2_.empty() );

	result.axis1 = XYZ( 0, 0, 0 );
	for ( AtomID const & atomid : axis1_ ) {
		result.axis1 += pose.xyz( atomid );
	}
	result.axis1 /= axis1_.size();

	result.axis2 = XYZ( 0, 0, 0 );
	for ( AtomID const & atomid : axis2_ ) {
		result.axis2 += pose.xyz( atomid );
	}
	result.axis2 /= axis2_.size();

	return result;
}

Real
AtomToAxisConstraint::dist( pose::Pose const & pose ) const {
	auto && assert_atomid = [&]( AtomID const & atom ){
		if ( ! pose.atom_tree().has( atom ) ) {
			TR << "AtomToAxisConstraint: cannot find atom " << atom << std::endl;
			utility_exit_with_message( "AtomToAxisConstraint: cannot find atom" );
		}
		conformation::Residue const& res = pose.residue( atom.rsd() );
		//Is this redundant? I'm just copying from AtomPairConstraint
		if ( atom.atomno() == 0 || atom.atomno() > res.natoms() ) {
			TR << "AtomToAxisConstraint: atom1 out of bounds " << atom << std::endl;
			utility_exit_with_message( "AtomToAxisConstraint: atom out of bounds" );
		}
	};

	assert_atomid( atom_ );
	for ( AtomID const & atom : axis1_ ) assert_atomid( atom );
	for ( AtomID const & atom : axis2_ ) assert_atomid( atom );

	XYZs const xyzs = calc_xyzs( pose );
	return pose.xyz( atom_ ).distance( closest_point_on_line( xyzs ) );
}

XYZs
AtomToAxisConstraint::calc_xyzs(
	func::XYZ_Func const & xyz
) const {
	XYZs result;
	result.atom = xyz( atom_ );

	debug_assert( ! axis1_.empty() );
	debug_assert( ! axis2_.empty() );

	result.axis1 = XYZ( 0, 0, 0 );
	for ( AtomID const & atomid : axis1_ ) {
		result.axis1 += xyz( atomid );
	}
	result.axis1 /= axis1_.size();

	result.axis2 = XYZ( 0, 0, 0 );
	for ( AtomID const & atomid : axis2_ ) {
		result.axis2 += xyz( atomid );
	}
	result.axis2 /= axis2_.size();

	return result;
}

Real
AtomToAxisConstraint::dist( func::XYZ_Func const & xyz ) const {
	XYZs const xyzs = calc_xyzs( xyz );
	return xyz( atom_ ).distance( closest_point_on_line( xyzs ) );
}

func::Func const & AtomToAxisConstraint::get_func() const {
	return *func_;
}

Size AtomToAxisConstraint::show_violations(
	std::ostream& out,
	pose::Pose const& pose,
	Size verbose_level,
	Real threshold
) const {
	return func_->show_violations( out, dist( pose ), verbose_level, threshold );
}

/// @brief Copies the data from this Constraint into a new object and returns an OP
/// atoms are mapped to atoms with the same name in dest pose ( e.g. for switch from centroid to fullatom )
/// if a sequence_mapping is present it is used to map residue numbers .. NULL = identity mapping
/// to the new object. Intended to be implemented by derived classes.
ConstraintOP
AtomToAxisConstraint::remapped_clone(
	pose::Pose const & src,
	pose::Pose const& dest,
	id::SequenceMappingCOP smap
) const {
	auto && update = [&](AtomID const & atom){
		id::NamedAtomID named_atom( atom_id_to_named_atom_id( atom, src ) );
		if ( smap != nullptr ) {
			named_atom.rsd() = (*smap)[ atom.rsd() ];
		}
		if ( named_atom.rsd() == 0 ) {
			TR << "Can't find new residue for AtomID: " << atom << std::endl;
			runtime_assert( false );
		}
		id::AtomID new_id( named_atom_id_to_atom_id( named_atom, dest ) );
		runtime_assert( new_id.valid() );
		return new_id;
	};

	id::AtomID const new_atom = update( atom_ );

	utility::vector1< AtomID > new_axis1;
	for ( auto const & atom : axis1_ ) new_axis1.push_back( update( atom ) );

	utility::vector1< AtomID > new_axis2;
	for ( auto const & atom : axis2_ ) new_axis2.push_back( update( atom ) );

	return utility::pointer::make_shared< AtomToAxisConstraint >( new_atom, new_axis1, new_axis2, func_, score_type() );
}

bool
AtomToAxisConstraint::operator == ( Constraint const & other_cst ) const
{
	if ( !           same_type_as_me( other_cst ) ) return false;
	if ( ! other_cst.same_type_as_me(     *this ) ) return false;

	auto const & other( static_cast< AtomToAxisConstraint const & > (other_cst) );
	if ( atom_ != other.atom_ ) return false;
	if ( axis1_.size() != other.axis1_.size() ) return false;
	if ( axis2_.size() != other.axis2_.size() ) return false;

	for ( core::Size ii = 1; ii <= axis1_.size(); ++ii ) {
		if ( axis1_[ii] != other.axis1_[ii] ) return false;
	}
	for ( core::Size ii = 1; ii <= axis2_.size(); ++ii ) {
		if ( axis2_[ii] != other.axis2_[ii] ) return false;
	}

	if ( this->score_type() != other.score_type() ) return false;

	return func_ == other.func_ || ( func_ && other.func_ && *func_ == *other.func_ );
}

bool
AtomToAxisConstraint::same_type_as_me( Constraint const & other ) const
{
	return dynamic_cast< AtomToAxisConstraint const * > (&other);
}

void AtomToAxisConstraint::show( std::ostream& out ) const {
	out << "AtomToAxisConstraint ("
		<< atom_.atomno() << "," << atom_.rsd() << ")\n";
	func_->show( out );
}

void AtomToAxisConstraint::show_def( std::ostream& out, pose::Pose const& pose ) const {
	out << type() << " " << atom_id_to_named_atom_id( atom_, pose ) << " ";
	func_->show_definition( out );
}

// atom deriv
void
AtomToAxisConstraint::fill_f1_f2(
	AtomID const & atom,
	func::XYZ_Func const & xyz,
	Vector & F1,
	Vector & F2,
	EnergyMap const & weights
) const {
	if ( atom != atom_ ) return;

	XYZs const xyzs = calc_xyzs( xyz );
	XYZ const closest = closest_point_on_line( xyzs );

	Real dist(0.0);
	Vector f1(0.0), f2(0.0);
	numeric::deriv::distance_f1_f2_deriv( xyz( atom ), closest, dist, f1, f2 );
	Real wderiv( weights[ this->score_type() ] * dfunc( dist ));
	F1 += wderiv * f1;
	F2 += wderiv * f2;
}

ConstraintOP
AtomToAxisConstraint::remap_resid( core::id::SequenceMapping const & ) const
{
	utility_exit_with_message( "AtomToAxisConstraint does not support remapping yet" );
	return nullptr;
}

/// @details one line definition "AtomPairs atom1 res1 atom2 res2 function_type function_definition"
void
AtomToAxisConstraint::read_def(
	std::istream &,
	core::pose::Pose const &,
	func::FuncFactory const &
) {
	utility_exit_with_message( "AtomToAxisConstraint does not support read_def yet" );
} // read_def

void
AtomToAxisConstraint::setup_for_scoring(
	func::XYZ_Func const &,
	ScoreFunction const &
) const {} //Do nothing.

// functions
Real
AtomToAxisConstraint::func( Real const theta ) const
{
	return func_->func( theta );
}

// deriv
Real
AtomToAxisConstraint::dfunc( Real const theta ) const
{
	return func_->dfunc( theta );
}

void
AtomToAxisConstraint::set_func( func::FuncOP setting ) {
	func_ = setting;
}

} // constraints
} // scoring
} // core

#ifdef    SERIALIZATION

/// @brief Automatically generated serialization method
template< class Archive >
void
core::scoring::constraints::AtomToAxisConstraint::save( Archive & arc ) const {
	arc( cereal::base_class< Constraint >( this ) );
	arc( CEREAL_NVP( atom_ ) ); // AtomID
	arc( CEREAL_NVP( axis1_ ) ); // AtomIDs
	arc( CEREAL_NVP( axis2_ ) ); // AtomIDs
	arc( CEREAL_NVP( func_ ) ); // core::scoring::func::FuncOP
}

/// @brief Automatically generated deserialization method
template< class Archive >
void
core::scoring::constraints::AtomToAxisConstraint::load( Archive & arc ) {
	arc( cereal::base_class< Constraint >( this ) );
	arc( atom_ ); // AtomID
	arc( axis1_ ); // AtomIDs
	arc( axis2_ ); // AtomIDs
	arc( func_ ); // core::scoring::func::FuncOP
}

SAVE_AND_LOAD_SERIALIZABLE( core::scoring::constraints::AtomToAxisConstraint );
CEREAL_REGISTER_TYPE( core::scoring::constraints::AtomToAxisConstraint )

CEREAL_REGISTER_DYNAMIC_INIT( core_scoring_constraints_AtomToAxisConstraint )
#endif // SERIALIZATION
