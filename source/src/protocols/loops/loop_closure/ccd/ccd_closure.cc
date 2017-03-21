// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/loops/loop_closure/ccd.ccd_closure.cc
/// @brief  Method definitions for cyclic coordinate descent loop closure.
/// @author Phil Bradley
/// @author Brian Weitzner
/// @author Labonte <JWLabonte@jhu.edu>


// Unit Header
#include <protocols/loops/loop_closure/ccd/ccd_closure.hh>

// Project Headers
#include <core/types.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/kinematics/Stub.hh>
#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

//Utility Headers
#include <utility/assert.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/vector1.hh>

// Numeric Headers
#include <numeric/xyz.functions.hh>
#include <numeric/numeric.functions.hh>

// Basic Headers
#include <basic/prof.hh>
#include <basic/basic.hh>
#include <basic/Tracer.hh>

// External Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>


using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

using namespace core;
using utility::vector1;
typedef numeric::xyzMatrix< Real > Matrix;

static THREAD_LOCAL basic::Tracer TR( "protocols.loops.loop_closure.ccd.ccd_closure" );

// FIXME: This is duplicated code currently with that in the Mover.
/// @param <coords>: an array of coordinates of the main-chain atoms for every residue in the Pose
/// @param <torsions>: an array of values for the main-chain torsion angles for every residue in the Pose
void
load_coords_and_torsions(
	pose::Pose const & pose,
	utility::vector1< utility::vector1< core::Vector > > & coords,
	utility::vector1< utility::vector1< core::Angle > > & torsions )
{
	// TODO: This is a useful utility function that should probably be moved to core
	//
	// I would also toss in a vector1< bool > (we would test its size()) or an atom ID map that can be used as a mask...
	// Of course, it would then be the protocol writer's responsibility to keep track of the vector index to residue
	// number mapping.
	//
	// In the case of loop modeling, this is really easy -- index 1 is the first residue of the loop, we can
	// test that the size() of the vector is the same as the loop length, etc.
	Size const nres( pose.size() );
	coords.resize( nres );
	torsions.resize( nres );

	for ( core::uint i = 1; i <= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue( i ) );

		// This needs to be private data of the new class. ~Labonte
		// Also, here, we define this for each residue, where higher in the code, we do not. ~Labonte
		Size const n_mainchain_atoms( rsd.mainchain_atoms().size() );
		if ( n_mainchain_atoms ) {
			coords[i].resize( n_mainchain_atoms );
			torsions[i].resize( n_mainchain_atoms );
			for ( core::uint j = 1; j <= n_mainchain_atoms; ++j ) {
				coords[i][j] = rsd.atom( rsd.mainchain_atoms()[ j ] ).xyz();
				torsions[i][j] = rsd.mainchain_torsion( j );
			}
		}
	}
}

// FIXME: This is duplicated code currently with that in the Mover.
void
get_overlap_pos(
	utility::vector1< utility::vector1< core::Vector > > const & coords,
	utility::vector1< utility::vector1< core::Angle > > const & torsions,
	core::uint const cutpoint,
	core::uint const direction,
	core::Angle const bond_angle1,
	core::Length const bond_length,
	core::Angle const bond_angle2,
	Matrix & M )
{
	using numeric::conversions::radians;
	using namespace utility;
	using kinematics::Stub;

	bool const local_debug( false/*true*/ );

	Size const n_mainchain_atoms( coords[cutpoint].size() );
	debug_assert( n_mainchain_atoms >= 3 ); // only true for proteins. what if you want to close a sugar loop?

	// These  magic numbers are correct for proteins -- what should they be for sugars?
	// pack the coords, torsions, and angles
	vector1< Vector > atoms( 6 );    // should this be vector1< xyzVecotr >?
	vector1< Angle > dihedrals( 3 ); // psi[i], omega[i], phi[i + 1]
	vector1< Angle > angles( 2 );    // theta_lower, theta_upper

	if ( direction == 1 ) { // eg N->CA->C->N'->CA'->C'
		for ( int i=1; i<= 3; ++i ) atoms[i  ] = coords[ cutpoint   ][ n_mainchain_atoms-3+i ]; // 1,2,3: N, CA, C
		for ( int i=1; i<= 3; ++i ) atoms[i+3] = coords[ cutpoint+1 ][ i       ]; // 4,5,6: N',CA',C'
		dihedrals[1] = radians( torsions[cutpoint  ][ n_mainchain_atoms-1 ] ); // N-CA-C-N'
		dihedrals[2] = radians( torsions[cutpoint  ][ n_mainchain_atoms   ] ); // CA-C-N'-CA'
		dihedrals[3] = radians( torsions[cutpoint+1][ 1     ] ); // C-N'-CA'-C'
		angles[1] = bond_angle1; // CA-C-N' already in radians
		angles[2] = bond_angle2; // C-N'-CA' already in radian
	} else { // eg N<-CA<-C<-N'<-CA'<-C'
		for ( int i=1; i<= 3; ++i ) atoms[7-i  ] = coords[ cutpoint   ][ n_mainchain_atoms-3+i ]; // 6,5,4: N, CA, C
		for ( int i=1; i<= 3; ++i ) atoms[7-i-3] = coords[ cutpoint+1 ][ i       ]; // 3,2,1: N',CA',C'
		dihedrals[3] = radians( torsions[cutpoint  ][ n_mainchain_atoms-1 ] ); // N-CA-C-N'
		dihedrals[2] = radians( torsions[cutpoint  ][ n_mainchain_atoms   ] ); // CA-C-N'-CA'
		dihedrals[1] = radians( torsions[cutpoint+1][ 1     ] ); // C-N'-CA'-C'
		angles[2] = bond_angle1; // CA-C-N' already in radians
		angles[1] = bond_angle2; // C-N'-CA' already in radians
	}

	//build the two overlap atoms using bond_angle...
	Vector atom3p, atom4p;
	{
		// atom4p is build forward, eg given N-CA-C, how to build N'
		Stub const stub321( atoms[3], atoms[2], atoms[1] );
		atom4p = stub321.spherical( dihedrals[1], angles[1], bond_length );

		// atom3p is built backward, eg given N'-CA'-C, how to build C // how to build C', right? ~BDW
		Stub const stub456( atoms[4], atoms[5], atoms[6] );
		atom3p = stub456.spherical( dihedrals[3], angles[2], bond_length );
	}

	Stub const stub1( atom4p, atoms[3], atoms[2] );
	Stub const stub2( atoms[4], atom3p, atoms[5] );

	for ( int i=1; i<= 3; ++i ) {
		M.col( i, stub1.local2global( stub2.global2local( atoms[3+i] ) ) );
	}

	// M now corresponds to a torsion at the inter-residue bond of 0.0, need to adjust to proper torsion
	Matrix const R( numeric::rotation_matrix_radians( ( atom4p - atoms[3] ).normalized(), dihedrals[2] ) );
	M.col( 2, R * ( M.col(2) - atom4p ) + atom4p );
	M.col( 3, R * ( M.col(3) - atom4p ) + atom4p );

	// confirm that the bond length, angle, and torsions are satisfied
	if ( local_debug ) {
		using basic::subtract_radian_angles;
		ASSERT_ONLY( Real const dihedral1( dihedral_radians( atoms[1], atoms[2], atoms[3], M.col(1) ) ););
		ASSERT_ONLY(Real const dihedral2( dihedral_radians( atoms[2], atoms[3], M.col(1), M.col(2) ) ););
		ASSERT_ONLY(Real const dihedral3( dihedral_radians( atoms[3], M.col(1), M.col(2), M.col(3) ) ););

		ASSERT_ONLY(Real const angle1( std::acos( dot( ( atoms[3] - atoms[2] ).normalized(), ( M.col(1) - atoms[3] ).normalized() ) ) ););
		ASSERT_ONLY(Real const angle2( std::acos( dot( ( M.col(2) - M.col(1) ).normalized(), ( M.col(1) - atoms[3] ).normalized() ) ) ););
		ASSERT_ONLY(Real const length( atoms[3].distance( M.col(1) ) ););

		debug_assert( std::abs( subtract_radian_angles( dihedral1, dihedrals[1] ) ) < 1e-3 );
		debug_assert( std::abs( subtract_radian_angles( dihedral2, dihedrals[2] ) ) < 1e-3 );
		debug_assert( std::abs( subtract_radian_angles( dihedral3, dihedrals[3] ) ) < 1e-3 );
		debug_assert( std::abs( subtract_radian_angles( angle1, angles[1] ) ) < 1e-3 );
		debug_assert( std::abs( subtract_radian_angles( angle2, angles[2] ) ) < 1e-3 );
		debug_assert( std::abs( length - bond_length ) < 1e-3 );
	}

	// restore proper chain order for M if folding backward
	if ( direction == 2 ) {
		Matrix const tmp( M );
		for ( int i=1; i<= 3; ++i ) M.col( i, tmp.col(4-i) );
	}
}

Distance
compute_single_direction_deviation(
	vector1< vector1< Vector > > const & coords,
	vector1< vector1< Angle > > const & torsions,
	int const cutpoint,
	core::uint const direction,
	Angle const ideal_theta_lower,
	Length const bond_length,
	Angle const bond_angle2,
	Matrix const & F
) {
	// Matrix F(0);
	Matrix M(0);

	get_overlap_pos( coords, torsions, cutpoint, direction, ideal_theta_lower, bond_length, bond_angle2, M );

	DistanceSquared deviation = 0.0;
	for ( int i = 1; i <= 3; ++i ) {
		for ( int j = 1; j <= 3; ++j ) {
			deviation += numeric::square( M(j,i) - F(j,i) );
		}
	}
	return sqrt( deviation / 3 );
}

// TODO: Move to a util file.
//Pulling out the deviation calculations from fast_ccd_closure, returns a pair of Reals containing forward and backward deviations
std::pair<core::Real, core::Real>
get_deviation( core::pose::Pose const & pose, core::uint const cutpoint )
{
	// Store the ideal geometry of the residues across the cutpoint
	Angle const ideal_theta_lower( pose.residue( cutpoint ).upper_connect().icoor().theta() );// CA-C=N bond angle
	Angle const bond_angle2( pose.residue( cutpoint+1 ).lower_connect().icoor().theta() ); // C=N-CA bond angle
	Length const bond_length( pose.residue( cutpoint+1 ).lower_connect().icoor().d() ); // C=N distance

	Size const n_mainchain_atoms( pose.residue( cutpoint ).mainchain_atoms().size() );

	vector1< vector1< Vector > > coords;
	vector1< vector1< Real > > torsions;
	load_coords_and_torsions( pose, coords, torsions );

	// Store relevant coordinates for deviation calculation:
	Matrix F(0);

	// forward_deviation:
	for ( int i = 1; i <= 3; ++i ) {
		F.col( i, coords[ cutpoint + 1 ][ i ] );
	}
	Distance forward_deviation = compute_single_direction_deviation( coords, torsions, cutpoint, 1,
		ideal_theta_lower, bond_length, bond_angle2, F);

	// backward_deviation:
	// Overwrite exisiting values in F -- we don't need them anymore!
	for ( int i = 1; i <= 3; ++i ) {
		F.col( i, coords[ cutpoint ][ n_mainchain_atoms - 3 + i ] );
	}
	Distance backward_deviation = compute_single_direction_deviation( coords, torsions, cutpoint, 2,
		ideal_theta_lower, bond_length, bond_angle2, F);
	return std::make_pair(forward_deviation, backward_deviation);
}

} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols
