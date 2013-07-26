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
/// @author

#include <core/scoring/dna/base_geometry.hh>
#include <core/scoring/dna/BasePartner.hh>
#include <core/scoring/dna/DNA_BasePotential.hh>

#include <core/scoring/ScoringManager.hh>

#include <core/types.hh>
#include <basic/Tracer.hh>
#include <basic/basic.hh>
#include <core/conformation/Residue.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/kinematics/Stub.hh>
#include <core/pose/Pose.hh>

// AUTO-REMOVED #include <utility/vector1.hh>
#include <numeric/conversions.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyz.functions.hh>

#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

#include <utility/vector1.hh>
#include <cmath>


namespace core {
namespace scoring {
namespace dna {

using namespace ObjexxFCL;
using namespace ObjexxFCL::fmt;

typedef numeric::xyzMatrix< Real > Matrix;
	//using kinematics::Stub::Matrix;
typedef utility::vector1< Real > Params;

////////////////////////////////////////////////////////////////////////////
/// FWD

Vector
lsf_normal(
	utility::vector1< Vector > const & atoms_in
);

////////////////////////////////////////////////////////////////////////////
void
get_base_pucker(
	conformation::Residue const & rsd,
	std::pair< std::string, int > & pucker
)
{

	utility::vector1< std::string > names;
	names.push_back( "C1'" );
	names.push_back( "C2'" );
	names.push_back( "C3'" );
	names.push_back( "C4'" );
	names.push_back( "O4'" );

	utility::vector1< Vector > atoms;
	for ( int i=1; i<= 5; ++i ) {
		atoms.push_back( rsd.xyz( names[i] ) );
	}

	Real mindot = 1000.0;
	bool exxo( false );
	for ( int ii=1; ii<= 5; ++ii ) {

		Vector n12 = (( atoms[2]-atoms[1] ).cross( atoms[3]-atoms[2] ) ).normalized();
		Real dot = std::fabs( n12.dot( ( atoms[4]-atoms[3] ).normalized() ) );
		if ( dot < mindot ) {
			// get pucker
			//Real pucker_dot = n12.dot( ( atoms[5] - Real(0.5) * ( atoms[4] + atoms[1] ) ).normalized() );

			mindot = dot;
			pucker.first = names[5];
			exxo = ( n12.dot( ( atoms[5] - Real(0.5) * ( atoms[4] + atoms[1] ) ).normalized() ) > 0.0 );
		}

		atoms.push_back( atoms[1] );
		atoms.erase( atoms.begin() );

		names.push_back( names[1] );
		names.erase( names.begin() );

	}


	// additional integer for scannability
	{
		int const atom_index( std::find( names.begin(), names.end(), pucker.first ) - names.begin() );
		int const sign_index( exxo ? 0 : 1 );
		if ( atom_index%2 == sign_index ) pucker.second = atom_index+1;
		else                              pucker.second = atom_index-4;
	}

	if ( exxo ) pucker.first += " exxo";
	else pucker.first += " endo";
}


///////////////////////////////////////////////////////////////////////////////
std::string
dihedral_bin( Real phi )
{
	phi = basic::periodic_range( phi, 360.0 );
	if      ( -120 <= phi && phi <   0 ) return "g-";
	else if (    0 <= phi && phi < 120 ) return "g+";
	else return "t ";
}

///////////////////////////////////////////////////////////////////////////////
std::string
get_DNA_backbone_bin( conformation::Residue const & rsd )
{
	std::string ag_bin, b12_bin;
	Real const alpha  ( rsd.mainchain_torsion( 1 ) );
	Real const gamma  ( rsd.mainchain_torsion( 3 ) );
	Real const epsilon( rsd.mainchain_torsion( 5 ) );
	Real const zeta   ( rsd.mainchain_torsion( 6 ) );

	if ( rsd.is_lower_terminus() ) ag_bin = "-- --";
	else ag_bin = dihedral_bin( alpha ) + " " + dihedral_bin( gamma );

	if ( rsd.is_upper_terminus() ) b12_bin = "--   0";
	else {
		Real const dev( basic::subtract_degree_angles( epsilon, zeta ) );
		if ( dev < 0 ) {
			b12_bin = "B1";
		} else {
			b12_bin = "B2";
		}
		b12_bin += ObjexxFCL::right_string_of( static_cast<int>( dev/10.0 ), 4 );
	}
	return ag_bin + " " + b12_bin;
}


///////////////////////////////////////////////////////////////////////////////
// y-axis goes from a1 to a2
//
void
get_y_axis_atoms(
	chemical::ResidueType const & rsd_type,
	int const strand, // 1 or 2
	std::string & a1,
	std::string & a2
)
{
	using namespace chemical;
	if ( rsd_type.aa() == na_ade || rsd_type.aa() == na_gua || rsd_type.aa() == na_rgu || rsd_type.aa() == na_rad ) {
		if ( strand == 1 ) {
			a1 = "N1";
			a2 = "C4";
		} else {
			a1 = "C4";
			a2 = "N1";
		}
	} else if ( rsd_type.aa() == na_cyt || rsd_type.aa() == na_thy || rsd_type.aa() == na_rcy || rsd_type.aa() == na_ura ) {
		if ( strand == 1 ) {
			a1 = "N3";
			a2 = "C6";
		} else {
			a1 = "C6";
			a2 = "N3";
		}
	} else {
		std::cout << rsd_type.aa() << " is unknown to me, tovarisch!\n"; // AM
		utility_exit();
	}
}
///////////////////////////////////////////////////////////////////////////////
Vector
get_y_axis(
	conformation::Residue const & rsd,
	int const strand
)
{
	std::string a1,a2;
	get_y_axis_atoms( rsd.type(), strand, a1, a2 );
	return ( rsd.xyz( a2 ) - rsd.xyz( a1 ) ).normalized();
}



///////////////////////////////////////////////////////////////////////////////
Vector
get_base_pair_y_axis_atom_xyz(
	conformation::Residue const & rsd
)
{
	using namespace chemical;
	if ( rsd.aa() == na_ade || rsd.aa() == na_gua || rsd.aa() == na_rgu || rsd.aa() == na_rad ) {
		return rsd.xyz("C8");

	} else if ( rsd.aa() == na_thy || rsd.aa() == na_cyt || rsd.aa() == na_rcy || rsd.aa() == na_ura ) {
		return rsd.xyz("C6");

	}
	utility_exit_with_message("get_base_pair_y_axis_atom_xyz: bad aa");
	return Vector(0.0);
}

///////////////////////////////////////////////////////////////////////////////
// only uses 6-membered ring
//
Vector
get_z_axis(
	conformation::Residue const & rsd,
	Vector const & y_axis
)
{
	using namespace chemical;
	assert( rsd.is_DNA() || rsd.is_RNA() );
	Vector xx(0); // approximate x-axis direction
	if ( rsd.aa() == na_ade || rsd.aa() == na_gua || rsd.aa() == na_rgu || rsd.aa() == na_rad ) {
		xx = rsd.xyz("C5") + rsd.xyz("C6") - rsd.xyz("N3") - rsd.xyz("C2");
	} else if ( rsd.aa() == na_thy || rsd.aa() == na_cyt || rsd.aa() == na_rcy || rsd.aa() == na_ura ) {
		xx = rsd.xyz("C5") + rsd.xyz("C4") - rsd.xyz("N1") - rsd.xyz("C2");
	} else {
		utility_exit_with_message("get_z_axis: bad aa");
	}

	return xx.cross( y_axis ).normalized();
}


///////////////////////////////////////////////////////////////////////////////
Vector
strand_orientation_vector( conformation::Residue const & rsd, int const strand )
{
	Vector orient( rsd.xyz("C3'") - rsd.xyz("C4'" ) );
	if ( strand == 2 ) orient *= -1.0f;
	orient.normalize();
	return orient;
}

///////////////////////////////////////////////////////////////////////////////
// only uses 6-membered ring
//
Vector
get_z_axis(
	conformation::Residue const & rsd,
	Vector const & y_axis,
	int const strand,
	bool & flipped
)
{
	Vector z_axis( get_z_axis( rsd, y_axis ) );

	// flip z-axis if necessary
	Vector const orient( strand_orientation_vector( rsd, strand ) );
	if ( dot( z_axis, orient ) < 0.0 ) {
		flipped = true;
		z_axis *= -1.0f;
	} else {
		flipped = false;
	}
	assert( std::fabs( z_axis.dot( y_axis ) ) < 1e-3 );

	return z_axis;
}

///////////////////////////////////////////////////////////////////////////////
// only uses 6-membered ring
//
// dont care about flips??
Vector
get_z_axis(
	conformation::Residue const & rsd,
	Vector const & y_axis,
	int const strand
)
{
	bool flipped( false );
	return get_z_axis( rsd, y_axis, strand, flipped );
}



	/// helper fxn
bool
is_orthonormal(
	numeric::xyzMatrix< Real > const & M,
	//Matrix const & M,
	Real const tol
)
{
	Matrix X( M*M.transposed() );
	float dev( 0.0 );
	for( int i=1; i<=3; ++i ) {
		for( int j=1; j<=3; ++j ) {
			dev += std::fabs( X(i,j) - ( i == j ? 1.0 : 0.0 ) );
		}
	}
	return dev < tol;
}


///////////////////////////////////////////////////////////////////////////////
// stolen/borrowed from Alex Morozov!
//

kinematics::Stub
get_base_stub(
	conformation::Residue const & rsd,
	int const strand
)
{
	using numeric::conversions::degrees;

	std::string a1,a2;
	get_y_axis_atoms( rsd.type(), strand, a1, a2 );

	Vector const & v1( rsd.xyz( a1 ) );
	Vector const & v2( rsd.xyz( a2 ) );

	Vector x_axis, y_axis, z_axis, origin;

	y_axis = v2 - v1;
	y_axis.normalize();

	origin = 0.5f * (v1 + v2 );

	bool flipped( false );
	z_axis = get_z_axis( rsd, y_axis, strand, flipped );
	if ( flipped ) {
		basic::T( "core.scoring.dna.base_geometry", basic::t_warning ) << "base flip in get_base_stub!!!" << '\n';
	}
	assert( std::fabs( dot(y_axis, z_axis) ) < 1e-3 );

	x_axis = cross( y_axis, z_axis );
	x_axis.normalize();

	return kinematics::Stub( kinematics::Stub::Matrix::cols( x_axis, y_axis, z_axis ), origin );
}

///////////////////////////////////////////////////////////////////////////////

kinematics::Stub
get_base_pair_stub(
	conformation::Residue const & rsd1, // on strand I
	conformation::Residue const & rsd2  // on strand II
)
{
	using numeric::conversions::degrees;

	Vector const y1( get_base_pair_y_axis_atom_xyz( rsd1 ) );
	Vector const y2( get_base_pair_y_axis_atom_xyz( rsd2 ) );
	Vector const origin( Real( 0.5 )* ( y1 + y2 ) );
	Vector const y_axis( ( y1 - y2 ).normalized() );

	Vector const z1_axis( get_z_axis( rsd1, y_axis, 1 ) );
	Vector const z2_axis( get_z_axis( rsd2, y_axis, 2 ) );
	Vector z_axis;
	if ( z1_axis.dot( z2_axis ) < 0.0 ) {
		basic::T( "core.scoring.dna.base_geometry", basic::t_warning ) << "wacky base flip in get_base_pair_stub!!!" << '\n';
		z_axis = z1_axis;
	} else {
		z_axis = ( z1_axis + z2_axis ).normalized();
	}
	assert( std::fabs( y_axis.dot( z_axis ) ) <1e-3 );
	Vector x_axis( cross( y_axis, z_axis ) );
	x_axis.normalize(); // prob unnecessary

	return kinematics::Stub( kinematics::Stub::Matrix::cols( x_axis, y_axis, z_axis ), origin );
}

///////////////////////////////////////////////////////////////////////////////

kinematics::Stub
get_base_pair_stub_slow(
	conformation::Residue const & rsd1, // on strand I
	conformation::Residue const & rsd2  // on strand II
)
{
	using numeric::conversions::degrees;

	Vector const y1( get_base_pair_y_axis_atom_xyz( rsd1 ) );
	Vector const y2( get_base_pair_y_axis_atom_xyz( rsd2 ) );
	Vector const origin( Real( 0.5 )* ( y1 + y2 ) );
	Vector const y_axis( ( y1 - y2 ).normalized() );

	assert( rsd1.atom_is_backbone( rsd1.chi_atoms(1)[2] ) && !rsd1.atom_is_backbone( rsd1.chi_atoms(1)[3] ) &&
					rsd2.atom_is_backbone( rsd2.chi_atoms(1)[2] ) && !rsd2.atom_is_backbone( rsd2.chi_atoms(1)[3] ) );

	utility::vector1< Vector > basepair_atoms;

	Size first_base_sidechain_atom = rsd1.first_sidechain_atom();
	if ( rsd1.is_RNA() ) {
		++first_base_sidechain_atom;
	}
	for ( Size i = first_base_sidechain_atom; i<= rsd1.nheavyatoms(); ++i ) {
		basepair_atoms.push_back( rsd1.xyz(i) );
	}

	first_base_sidechain_atom = rsd2.first_sidechain_atom();
	if ( rsd1.is_RNA() ) {
		++first_base_sidechain_atom;
	}
	for ( Size i = first_base_sidechain_atom; i<= rsd2.nheavyatoms(); ++i ) {
		basepair_atoms.push_back( rsd2.xyz(i) );
	}

	Vector z_axis( lsf_normal( basepair_atoms ) );
	z_axis = ( z_axis - y_axis.dot( z_axis ) * y_axis ).normalized();
	assert( z_axis.is_normalized( 1e-3 ) && z_axis.dot( y_axis ) < 1e-3 );
	if ( z_axis.dot( strand_orientation_vector( rsd1, 1 ) ) < 0.0 ) z_axis *= Real(-1.0);

	Vector x_axis( cross( y_axis, z_axis ) );
	x_axis.normalize(); // prob unnecessary

	return kinematics::Stub( kinematics::Stub::Matrix::cols( x_axis, y_axis, z_axis ), origin );
}

///////////////////////////////////////////////////////////////////////////////
void
show_base_pair_params_with_z_scores(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	std::ostream & out
)
{
	DNA_BasePotential const & potential( ScoringManager::get_instance()->get_DNA_BasePotential() );

	using namespace fmt;
	Params params(6);
	get_base_pair_params(rsd1,rsd2,params);
	utility::vector1< Real > z_scores(6,0.0);

	potential.eval_base_pair_Z_scores( rsd1, rsd2, z_scores );
	Real dev(0.0);
	for ( Size i=1; i<= 6; ++i ) dev += z_scores[i]*z_scores[i];
	out << "BP_PARAMS " <<
		F(7,1,potential.base_pair_score( rsd1, rsd2 ) ) << F(7,1,dev) <<
		I(4,rsd1.seqpos()) << I(4,rsd2.seqpos()) << ' ' << rsd1.name1() << rsd2.name1() <<
		" Prop: " << F(7,1,params[1]) << F(6,1,z_scores[1]) <<
		" Buck: " << F(7,1,params[2]) << F(6,1,z_scores[2]) <<
		" Open: " << F(7,1,params[3]) << F(6,1,z_scores[3]) <<
		" Sher: " << F(7,2,params[4]) << F(6,1,z_scores[4]) <<
		" Strc: " << F(7,2,params[5]) << F(6,1,z_scores[5]) <<
		" Stag: " << F(7,2,params[6]) << F(6,1,z_scores[6]) << '\n';
}

///////////////////////////////////////////////////////////////////////////////
void
show_base_pair_params(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	std::ostream & out
)
{
	//DNA_BasePotential const & potential( ScoringManager::get_instance()->get_DNA_BasePotential() );

	using namespace fmt;
	Params params(6);
	get_base_pair_params(rsd1,rsd2,params);
	out << "BP_PARAMS " << I(4,rsd1.seqpos()) << I(4,rsd2.seqpos()) << ' ' << rsd1.name1() << rsd2.name1() <<
		" Prop:" << F(6,1,params[1]) <<
		" Buck:" << F(6,1,params[2]) <<
		" Open:" << F(6,1,params[3]) <<
		" Sher:" << F(6,2,params[4]) <<
		" Strc:" << F(6,2,params[5]) <<
		" Stag:" << F(6,2,params[6]) << '\n';
}

///////////////////////////////////////////////////////////////////////////////
void
show_base_pair_params(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
)
{
	show_base_pair_params( rsd1, rsd2, std::cout );
}
///////////////////////////////////////////////////////////////////////////////
/**
	 works as-is for base pair with stub1 ~ strand I, stub2 ~ strand II
	 for base-step we use the mapping:

	 step    pair
	 ------------
	 z    ->    y
	 y    ->    x
	 x    ->    z
	 ------------
	 i+1  ->    I  (stub1)
	 i    ->   II  (stub2)
	 ------------
	 ------------


	 mid-stub coordsys generated by aligning the y-axes.

	 params[1] = oriented angle from the stub2 z-axis to the stub1 z-axis after aligning y-axes
	           = ( same as oriented angle from the stub2 x-axis to the stub1 x-axis after aligning y-axes )
						 = atan2( dot( stub1-z, stub2-x ), dot( stub1-z, stub2-z ) )
						 = atan2( y-coord( stub1-z ), x-coord( stub1-z ) ) in z-x coordinate system.
						 = propeller (twist)

	 params[2] = x-component of the rotation axis from stub2-y to stub1-y, weighted by rotation angle
	           = buckle (roll)

	 params[3] = z-component of the rotation axis from stub2-y to stub1-y, weighted by rotation angle
	           = opening (tilt)

	 params[4] = x-coordinate of the vector from stub2-origin to stub1-origin
	           = shear (slide)

	 params[5] = y-coordinate of the vector from stub2-origin to stub1-origin
	           = stretch (rise)

	 params[6] = z-coordinate of the vector from stub2-origin to stub1-origin
	           = stagger (shift)


	 NAMES       pair           step
	 ------------------------------
	 params[1] = propeller      twist
	 params[2] = buckle         roll
	 params[3] = opening        tilt
	 params[4] = shear          slide
	 params[5] = stretch        rise
	 params[6] = stagger        shift

**/

void
get_stub_stub_params(
	kinematics::Stub const & stub1,
	kinematics::Stub const & stub2,
	Params & params
)
{
	using namespace std;
	using numeric::conversions::degrees;
	using numeric::arccos;

	bool const local_debug( true ); // PBHACK!!!!!!!!!!!!!!!!!!!!!!
	params.resize(6);

	// copy matrices
	Matrix M1( stub1.M ), M2( stub2.M );

	assert( is_orthonormal( M1, 1e-3 ) );
	assert( is_orthonormal( M2, 1e-3 ) );

#ifndef NDEBUG
	bool base_flipped = false;
#endif
	if ( dot( M1.col_z(), M2.col_z() ) < 0.0 ) {
#ifndef NDEBUG
		base_flipped = true;
#endif
		basic::T("core.scoring.base_geometry") << "get_stub_stub_params: base flip!!!\n";
		for (Size i = 1; i <= 6; ++i) {
			params[i] = -9999;
		}
		return;
		//utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	}

	// get angle between the y-axes
	Real const gamma( arccos( dot( M1.col_y(), M2.col_y() ) ) );

	Vector const bo( ( cross( M2.col_y(), M1.col_y() ) ).normalized() );

	Matrix R_gamma_2( rotation_matrix( bo, gamma/2.0f ) );

	M2 = R_gamma_2 * M2;
	M1 = R_gamma_2.transposed() * M1;

	assert( is_orthonormal( M1, 1e-3 ) );
	assert( is_orthonormal( M2, 1e-3 ) );

	// build mid-stub triad
	assert( M1.col_y().distance( M2.col_y() ) < 1e-3 );
	assert( std::fabs( dot( bo, M1.col_y() ) ) < 1e-3 );

	Matrix MBT;
	MBT.col_y( M1.col_y() );

	assert( std::fabs( dot( M1.col_z(), MBT.col_y() ) ) < 1e-3 );
	assert( std::fabs( dot( M2.col_z(), MBT.col_y() ) ) < 1e-3 );
	assert( std::fabs( dot( M1.col_x(), MBT.col_y() ) ) < 1e-3 );
	assert( std::fabs( dot( M2.col_x(), MBT.col_y() ) ) < 1e-3 );

	// get
	MBT.col_x( ( 0.5f * ( M1.col_x() + M2.col_x() ) ).normalized() );
	MBT.col_z( ( 0.5f * ( M1.col_z() + M2.col_z() ) ).normalized() );

	assert( is_orthonormal( MBT, 1e-3 ) );

	// angular params

	// propellor
	// z,x,y make rh coord system
	params[1] = std::atan2( dot( M1.col_z(), M2.col_x() ),
													dot( M1.col_z(), M2.col_z() ) );

	assert( !local_debug || ( std::fabs( std::fabs( params[1] ) - arccos( dot( M1.col_z(), M2.col_z() ) ) ) < 1e-2 ) );

	// buckle:
	params[2] = gamma * dot( bo, MBT.col_x() );

	// opening
	params[3] = gamma * dot( bo, MBT.col_z() );

	// translational params
	Vector const displacement( stub1.v - stub2.v );

	params[4] = dot( displacement, MBT.col_x() );
	params[5] = dot( displacement, MBT.col_y() );
	params[6] = dot( displacement, MBT.col_z() );

	/////////////
	// remove this debugging stuff, preserved in old code at the end of the file
	if ( local_debug ) {
		{ // sin gamma version of params[2] is a simple dot product:
			//Real const tmp1 = (sin( gamma ) * params[2] / gamma);
			//Real const tmp2 = dot( Vector( cross( stub2.M.col_y(), stub1.M.col_y() ) ), MBT.col_x() );
			assert( std::fabs( (sin( gamma ) * params[2] / gamma) -
						  dot( Vector( cross( stub2.M.col_y(), stub1.M.col_y() ) ), MBT.col_x() )) < 1e-2 );
		}

		{ // sin gamma version of params[3] is a simple dot product:
			//Real const tmp1( sin( gamma ) * params[3] / gamma );
			//Real const tmp2( dot( Vector( cross( stub2.M.col_y(), stub1.M.col_y() )), MBT.col_z() ) );
			assert( std::fabs( ( sin( gamma ) * params[3] / gamma ) -
						   dot( Vector( cross( stub2.M.col_y(), stub1.M.col_y() )), MBT.col_z() ) ) < 1e-2 );
		}

		// check sign conventions
		Real phi_prime( arccos( dot( bo, MBT.col_x() ) ) );
		if ( dot( cross( bo, MBT.col_x() ), MBT.col_y() ) < 0.0f ) {
			phi_prime *= -1.0f;
		}

		Vector tmp( cross( M2.col_z(), M1.col_z() ) );
		assert( cross( tmp, MBT.col_y() ).length() <1e-2 );

		//Real const p1x = std::asin( dot( MBT.col_y(), cross( M2.col_x(), M1.col_x() ) ) );
		//Real const p1z = std::asin( dot( MBT.col_y(), cross( M2.col_z(), M1.col_z() ) ) );
		assert( ( base_flipped ) ||
				( std::fabs( params[1] - asin( dot( MBT.col_y(), cross( M2.col_x(), M1.col_x() ) ) ) ) +
				  std::fabs( params[1] - asin( dot( MBT.col_y(), cross( M2.col_z(), M1.col_z() ) ) ) ) < 1e-2 ) );
		//std::cout << "equal? p1: " << params[1] << ' ' << p1x << ' ' << p1z <<
 		//	std::endl;

		//Real const p2 = gamma * cos( phi_prime );
		//Real const p3 = gamma * sin( phi_prime );
		//Real const dev( std::fabs( p2 - params[2] ) + std::fabs( p3 - params[3] ) );
		//std::cout << "dev: " << dev << std::endl;
		assert( std::fabs( gamma * cos( phi_prime ) - params[2] ) + std::fabs( gamma * sin( phi_prime ) - params[3] ) < 1e-2 );

		// check sign conventions
		assert( params[1] * dot( MBT.col_y(), cross( M2.col_x(), M1.col_x() ) ) > 0);
	}

	// convert to degrees
	params[1] = degrees( params[1] );
	params[2] = degrees( params[2] );
	params[3] = degrees( params[3] );

}


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
// checking these with X3DNA
//
// use find_pair followed by cehs:
//
// 1111  ~/X3DNA/bin/find_pair -t pap1_subset.pdb pap1_subset.pdb.inp
// 1112  ~/X3DNA/bin/cehs pap1_subset.pdb.inp
// 1113  new
// 1114  more pap1_subset.outc
//
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//
// aa1 is in strand1, aa2 is in strand2
//
void
get_base_pair_params(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	Params & params // output
)
{
	get_stub_stub_params( get_base_stub( rsd1, 1 /*strand*/ ),
												get_base_stub( rsd2, 2 /*strand*/ ), params );
}



///////////////////////////////////////////////////////////////////////////////
//
// aa1 is in strand1, aa2 is in strand2
//
void
get_base_step_params(
	conformation::Residue const & rsd11, // pair1 strand I
	conformation::Residue const & rsd12, // pair1 strand II
	conformation::Residue const & rsd21, // pair2 strand I
	conformation::Residue const & rsd22, // pair2 strand II
	Params & params // output
)
{
	using kinematics::Stub;

	assert( rsd21.seqpos() == rsd11.seqpos() + 1 && rsd12.seqpos() == rsd22.seqpos() + 1 );

	//Stub const stub1( get_base_pair_stub( rsd11, rsd12)), stub2( get_base_pair_stub( rsd21, rsd22 ) );
	Stub const stub1( get_base_pair_stub_slow( rsd11, rsd12)), stub2( get_base_pair_stub_slow( rsd21, rsd22 ) );

	get_stub_stub_params( Stub( Stub::Matrix::cols( stub2.M.col_y(), stub2.M.col_z(), stub2.M.col_x() ), stub2.v ),
												Stub( Stub::Matrix::cols( stub1.M.col_y(), stub1.M.col_z(), stub1.M.col_x() ), stub1.v ),
												params );

}



///////////////////////////////////////////////////////////////////////////////
bool
seqpos_is_base_step_anchor(
	Size const seqpos,
	pose::Pose const & pose
)
{
	BasePartner const & partner( retrieve_base_partner_from_pose( pose ) );
	conformation::Residue const & rsd( pose.residue( seqpos ) );

	return ( seqpos < pose.total_residue() && ( rsd.is_DNA()  || rsd.is_RNA() )&& !rsd.is_lower_terminus() && partner[ seqpos ] &&
	partner[ seqpos+1 ] && partner[seqpos] == partner[seqpos+1]+1 && partner[seqpos] != seqpos+1 );
}


///////////////////////////////////////////////////////////////////////////////

void
show_base_step_params(
	Size const seqpos,
	pose::Pose const & pose,
	std::ostream & out
)
{

	using namespace fmt;
	Params params(6);

	BasePartner const & partner( retrieve_base_partner_from_pose( pose ) );

	if ( !seqpos_is_base_step_anchor( seqpos, pose ) ) {
		out << "BS_PARAMS " << seqpos << " N/A\n";
		return;
	}

	conformation::Residue const & rsd11( pose.residue( seqpos   ) );
	conformation::Residue const & rsd12( pose.residue( partner[ seqpos ] ) );
	conformation::Residue const & rsd21( pose.residue( seqpos+1   ) );
	conformation::Residue const & rsd22( pose.residue( partner[ seqpos+1 ] ) );
	get_base_step_params( rsd11, rsd12, rsd21, rsd22, params );

	out << "BS_PARAMS " <<
		I(4,seqpos  ) << I(4,partner[seqpos  ]) << ' ' << rsd11.name1() << rsd12.name1() << " to " <<
		I(4,seqpos+1) << I(4,partner[seqpos+1]) << ' ' << rsd21.name1() << rsd22.name1() <<
		" Twst:" << F(6,1,params[1]) <<
		" Roll:" << F(6,1,params[2]) <<
		" Tilt:" << F(6,1,params[3]) <<
		" Slid:" << F(6,2,params[4]) <<
		" Rise:" << F(6,2,params[5]) <<
		" Shft:" << F(6,2,params[6]) << '\n';
}

void
show_new_base_step_params(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2
)
{
	DNA_BasePotential const & potential( ScoringManager::get_instance()->get_DNA_BasePotential() );

	using namespace fmt;
	Params params(6);
	get_base_step_params(rsd1,rsd2,params);
	utility::vector1< Real > z_scores(6,0.0);
	potential.eval_base_step_Z_scores( rsd1, rsd2, z_scores );
	Real dev(0.0);
	for ( Size i=1; i<= 6; ++i ) dev += z_scores[i]*z_scores[i];
	std::cout << "BS-params: " <<
		F(7,1,potential.base_step_score( rsd1, rsd2 ) ) << F(7,1,dev) <<
		I(4,rsd1.seqpos()) << I(4,rsd2.seqpos()) << ' ' << rsd1.name1() << rsd2.name1() <<
		" Twst: " << F(7,1,params[1]) << F(6,1,z_scores[1]) <<
		" Roll: " << F(7,1,params[2]) << F(6,1,z_scores[2]) <<
		" Tilt: " << F(7,1,params[3]) << F(6,1,z_scores[3]) <<
		" Shft: " << F(7,2,params[4]) << F(6,1,z_scores[4]) <<
		" Slid: " << F(7,2,params[5]) << F(6,1,z_scores[5]) <<
		" Rise: " << F(7,2,params[6]) << F(6,1,z_scores[6]) << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
void
show_base_step_params(
											pose::Pose const & pose,
											std::ostream & out
											)
{

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( seqpos_is_base_step_anchor( i, pose ) ) {
			show_base_step_params( i, pose, out );
		}
	}

}

///////////////////////////////////////////////////////////////////////////////
void
show_base_pair_params(
											pose::Pose const & pose,
											std::ostream & out
											)
{
	BasePartner const & partner( retrieve_base_partner_from_pose( pose ) );

	for ( Size i=1; i<= pose.total_residue(); ++i ) {
		if ( partner[i] > i ) {
			show_base_pair_params( pose.residue(i), pose.residue(partner[i]), out );
		}
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////
// might be useful at some point to get them from coordinates. in case atomtree has breaks...
//
// void
// get_dna_dihedrals(
// 									Size const seqpos,
// 									pose::Pose const & pose
// 									vector1< Real > & dihedrals
// 									)
// {
// 	dihedrals.resize( 7 );
// 	Residue const & rsd( pose.residue(seqpos) );


// }


/////////////////////////////////////////////////////////////////////////////////////////////
void
show_dna_geometry(
									pose::Pose const & pose,
									std::ostream & out
									)
{

	// dihedrals + a/g bin + typeI,II + pucker
	Size const nres( pose.total_residue() );

	for ( Size i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		if ( !rsd.is_DNA() || rsd.is_RNA() ) continue;
		std::pair< std::string, int > pucker;
		get_base_pucker( rsd, pucker );

		out << "DNA_DIHEDRALS " << I(4,i) << ' ' << rsd.name1() << ' ' <<
			pucker.first << right_string_of( pucker.second, 3 ) << ' ' << get_DNA_backbone_bin( rsd ) <<
			F(7,1,rsd.mainchain_torsion(1)) <<
			F(7,1,rsd.mainchain_torsion(2)) <<
			F(7,1,rsd.mainchain_torsion(3)) <<
			F(7,1,rsd.mainchain_torsion(4)) <<
			F(7,1,rsd.mainchain_torsion(5)) <<
			F(7,1,rsd.mainchain_torsion(6)) <<
			F(7,1,rsd.chi(1)) << '\n';
	}

	// base-pair params
	show_base_pair_params( pose, out );

	// base-step params
	show_base_step_params( pose, out );

}

///////////////////////////////////////////////////////////////////
/// if you really want the least-squares plane:
/// Thanks to Alex Morozov for implementing this:
///
/// V.Schomaker et al. Acta Cryst. (1959) 12, 600-604
/// D.M.Blow Acta Cryst. (1960) 13, 168
///
/// Note: no guarantee about the orientation of the returned vector...
///

Vector
lsf_normal(
	utility::vector1< Vector > const & atoms_in
)
{
	// translate the atoms so center of mass is at the origin
	Vector cm(0.0);
	int const natoms( atoms_in.size() );
	for ( int i=1; i<= natoms; ++i ) {
		cm += atoms_in[i];
	}

	cm /= natoms;

	utility::vector1< Vector > atoms; atoms.reserve( natoms );
	for ( int i=1; i<= natoms; ++i ) {
		atoms.push_back( atoms_in[i] - cm );
	}

	// Create a matrix which provides coeffs for the cubic equation.
	// Note that A(i,j) = A(j,i).
	//FArray2D_Real A( 3, 3, 0.0f );
	Matrix A( 0.0 );
	for ( int i = 1; i <= 3; ++i ) {
		for ( int j = i; j <= 3; ++j ) {
			for ( int a = 1; a <= natoms; ++a ) {
				A(i,j) += atoms[a](i) * atoms[a](j);//coords2fit(i,a)*coords2fit(j,a);
			}
		}
	}

	// Compute cubic eqn. coeffs:
	Real a = A(1,1) + A(2,2) + A(3,3); // alpha
	Real b = A(1,3)*A(1,3) + A(1,2)*A(1,2) + A(2,3)*A(2,3)
		- A(2,2)*A(3,3) - A(1,1)*A(3,3) - A(1,1)*A(2,2); // beta
	Real g = A(1,1)*A(2,2)*A(3,3) + 2*A(1,2)*A(2,3)*A(1,3)
		- A(1,1)*A(2,3)*A(2,3) - A(2,2)*A(1,3)*A(1,3)
		- A(3,3)*A(1,2)*A(1,2); // gamma

	// The solution to the cubic eqn. is obtained through neglecting the
	// lambda**3 term (lambda is small for nearly planar atom sets):
	// WARNING: bound to fail for REALLY nonplanar atom sets!
	Real lambda;
	if (b*b - 4*a*g > 0) {
		lambda = (-b - std::sqrt(b*b - 4*a*g) )/(2*a);
	} else {
		lambda = 0;
	}

	// debug
	//std::cout << "lambda = " << lambda << std::endl;
	//std::cout << "D = " << (b*b - 4*a*g) << std::endl;
	// Finally, compute the eigenvector corresponding to the least squares plane:

	Vector normal_f;
	normal_f(1) = (A(2,2)-lambda)*A(1,3) - A(1,2)*A(2,3);
	normal_f(2) = (A(1,1)-lambda)*A(2,3) - A(1,2)*A(1,3);
	normal_f(3) = A(1,2)*A(1,2) - (A(2,2)-lambda)*(A(1,1)-lambda);

	normal_f.normalize();
	return normal_f;
}
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////
/// saved for comparison and debugging purposes
/////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////

void
get_base_pair_params_old(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	Params & params // output
)
{
	using numeric::conversions::degrees;
	using numeric::arccos;

	bool const local_debug( false );

	params.resize(6);

	kinematics::Stub const stub1( get_base_stub( rsd1, 1 /*strand*/ )), stub2( get_base_stub( rsd2, 2 /*strand*/ ) );

	// copy matrices
	Matrix M1( stub1.M ), M2( stub2.M );

	assert( is_orthonormal( M1, 1e-3 ) );
	assert( is_orthonormal( M2, 1e-3 ) );

#ifndef NDEBUG
	bool base_flipped = false;
#endif
	if ( dot( M1.col_z(), M2.col_z() ) < 0.0 ) {
#ifndef NDEBUG
		base_flipped = true;
#endif
		basic::T("core.scoring.base_geometry") << "base_pair_params: base flip!!!\n";
		//utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	}

	// get angle between the y-axes
	Real const gamma( arccos( dot( M1.col_y(), M2.col_y() ) ) );

	Vector const bo( ( cross( M2.col_y(), M1.col_y() ) ).normalized() );

	Matrix R_gamma_2( rotation_matrix( bo, gamma/2.0f ) );

	M2 = R_gamma_2 * M2;
	M1 = R_gamma_2.transposed() * M1;

	assert( is_orthonormal( M1, 1e-3 ) );
	assert( is_orthonormal( M2, 1e-3 ) );

	// build mid-base-pair triad
	assert( M1.col_y().distance( M2.col_y() ) < 1e-3 );
	assert( std::fabs( dot( bo, M1.col_y() ) ) < 1e-3 );

	Matrix MBT;
	MBT.col_y( M1.col_y() );

	assert( std::fabs( dot( M1.col_z(), MBT.col_y() ) ) < 1e-3 );
	assert( std::fabs( dot( M2.col_z(), MBT.col_y() ) ) < 1e-3 );
	assert( std::fabs( dot( M1.col_x(), MBT.col_y() ) ) < 1e-3 );
	assert( std::fabs( dot( M2.col_x(), MBT.col_y() ) ) < 1e-3 );

	// get
	MBT.col_x( ( 0.5f * ( M1.col_x() + M2.col_x() ) ).normalized() );
	MBT.col_z( ( 0.5f * ( M1.col_z() + M2.col_z() ) ).normalized() );

	assert( is_orthonormal( MBT, 1e-3 ) );

	// angular params

	// propellor
	// z,x,y make rh coord system
	params[1] = std::atan2( dot( M1.col_z(), M2.col_x() ),
													dot( M1.col_z(), M2.col_z() ) );

	if ( local_debug ) {
		assert( std::fabs( std::fabs( params[1] ) - arccos( dot( M1.col_z(), M2.col_z() ) ) ) < 1e-2 );
	}

	// buckle:
	params[2] = gamma * dot( bo, MBT.col_x() );

	// opening
	params[3] = gamma * dot( bo, MBT.col_z() );

	// translational params
	Vector const displacement( stub1.v - stub2.v );

	params[4] = dot( displacement, MBT.col_x() );
	params[5] = dot( displacement, MBT.col_y() );
	params[6] = dot( displacement, MBT.col_z() );

	/////////////
	// debugging:
	if ( local_debug ) {
		{ // sin gamma version of params[2] is a simple dot product:
			//Real const tmp1 = sin( gamma ) * params[2] / gamma;
			//Real const tmp2 = dot( Vector( cross( stub2.M.col_y(), stub1.M.col_y() ) ), MBT.col_x() );
			assert( std::fabs(sin( gamma ) * params[2] / gamma -
						dot( Vector( cross( stub2.M.col_y(), stub1.M.col_y() ) ), MBT.col_x() )) < 1e-2 );
		}

		{ // sin gamma version of params[3] is a simple dot product:
			//Real const tmp1( sin( gamma ) * params[3] / gamma );
			//Real const tmp2( dot( Vector( cross( stub2.M.col_y(), stub1.M.col_y() )), MBT.col_z() ) );
			assert( std::fabs(sin( gamma ) * params[3] / gamma -
					dot( Vector( cross( stub2.M.col_y(), stub1.M.col_y() )), MBT.col_z() )) < 1e-2 );
		}

		// check sign conventions
		Real phi_prime( arccos( dot( bo, MBT.col_x() ) ) );
		if ( dot( cross( bo, MBT.col_x() ), MBT.col_y() ) < 0.0f ) {
			phi_prime *= -1.0f;
		}

		Vector tmp( cross( M2.col_z(), M1.col_z() ) );
		assert( cross( tmp, MBT.col_y() ).length() <1e-2 );

		//Real const p1x = std::asin( dot( MBT.col_y(), cross( M2.col_x(), M1.col_x() ) ) );
		//Real const p1z = std::asin( dot( MBT.col_y(), cross( M2.col_z(), M1.col_z() ) ) );
		assert( ( base_flipped ) ||
				( std::fabs( params[1] - asin( dot( MBT.col_y(), cross( M2.col_x(), M1.col_x() ) ) ) ) +
				  std::fabs( params[1] - asin( dot( MBT.col_y(), cross( M2.col_z(), M1.col_z() ) ) ) ) < 1e-2 ) );
		//std::cout << "equal? p1: " << params[1] << ' ' << p1x << ' ' << p1z <<
 		//	std::endl;

		//Real const p2 = gamma * cos( phi_prime );
		//Real const p3 = gamma * sin( phi_prime );
		//Real const dev( std::fabs( p2 - params[2] ) + std::fabs( p3 - params[3] ) );
		//std::cout << "dev: " << dev << std::endl;
		assert( std::fabs( gamma * cos( phi_prime ) - params[2] ) + std::fabs( gamma * sin( phi_prime ) - params[3] ) < 1e-2 );

		// check sign conventions
		assert( params[1] * dot( MBT.col_y(), cross( M2.col_x(), M1.col_x() ) ) > 0);
	}

	// convert to degrees
	params[1] = degrees( params[1] );
	params[2] = degrees( params[2] );
	params[3] = degrees( params[3] );

}

///////////////////////////////////////////////////////////////////////////////
//
//// saved for comparison and debugging purposes
void
get_base_step_params(
	conformation::Residue const & rsd1,
	conformation::Residue const & rsd2,
	Params & params // output
)
{
	using numeric::conversions::degrees;
	using numeric::arccos;
	using numeric::cross;

	bool const local_debug( false );

	params.resize(6);

	kinematics::Stub const stub1( get_base_stub( rsd1, 1 ) ), stub2( get_base_stub( rsd2, 1 ) ); // strand = 2

	// copy matrices
	Matrix M1( stub1.M ), M2( stub2.M );

	assert( is_orthonormal( M1, 1e-3 ) );
	assert( is_orthonormal( M2, 1e-3 ) );

	if ( dot( M1.col_z(), M2.col_z() ) < 0.0 ) {
		// BASE FLIP !!!!!!!!!!!!!!!!!!!!!!!!!!!!
		basic::T("core.scoring.base_geometry") << "base_pair_params: base flip!!!\n";
		//std::cout << "new_base_step_params: base flip!" << std::endl;
		//utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	}

	// get angle between the z-axes
	Real const gamma( arccos( dot( M1.col_z(), M2.col_z() ) ) );

	Vector const rt( ( cross( M2.col_z(), M1.col_z() ) ).normalized() );

	Matrix R_gamma_2( rotation_matrix( rt, gamma/2.0f ) );

	M2 = R_gamma_2 * M2;
	M1 = R_gamma_2.transposed() * M1;

	assert( is_orthonormal( M1, 1e-3 ) );
	assert( is_orthonormal( M2, 1e-3 ) );

	// build mid-base-pair triad
	assert( M1.col_z().distance( M2.col_z() ) < 1e-3 );
	assert( std::fabs( dot( rt, M1.col_z() ) ) < 1e-3 );

	Matrix MBT;
	MBT.col_z( M1.col_z() );

	assert( std::fabs( dot( M1.col_x(), MBT.col_z() ) ) < 1e-3 );
	assert( std::fabs( dot( M2.col_x(), MBT.col_z() ) ) < 1e-3 );
	assert( std::fabs( dot( M1.col_y(), MBT.col_z() ) ) < 1e-3 );
	assert( std::fabs( dot( M2.col_y(), MBT.col_z() ) ) < 1e-3 );

	// get
	MBT.col_y( ( 0.5f * ( M1.col_y() + M2.col_y() ) ).normalized() );
	MBT.col_x( ( 0.5f * ( M1.col_x() + M2.col_x() ) ).normalized() );

	assert( is_orthonormal( MBT, 1e-3 ) );

	// angular params

	// TWIST
	// x,y,z make rh coord system
	params[1] = atan2( dot( M1.col_x(), M2.col_y() ), dot( M1.col_x(), M2.col_x() ) );

	if ( local_debug ) {
		assert( std::fabs( std::fabs( params[1] ) - arccos( dot( M1.col_x(), M2.col_x() ) ) ) < 1e-2 );
	}


	// ROLL:
	params[2] = gamma * dot( rt, MBT.col_y() );

	// TILT
	params[3] = gamma * dot( rt, MBT.col_x() );

	// translational params
	Vector const displacement( stub1.v - stub2.v );

	params[4] = dot( displacement, MBT.col_x() ); // SHIFT
	params[5] = dot( displacement, MBT.col_y() ); // SLIDE
	params[6] = dot( displacement, MBT.col_z() ); // RISE

	/////////////
	// debugging:
	if ( local_debug ) {
		{ // sin gamma version of params[2] (roll) is a simple dot product:
			//Real const tmp1 = sin( gamma ) * params[2] / gamma;
			//Real const tmp2 = dot( Vector( cross( stub2.M.col_z(), stub1.M.col_z() ) ), MBT.col_y() );
			assert( std::fabs( sin( gamma ) * params[2] / gamma -
						 dot( Vector( cross( stub2.M.col_z(), stub1.M.col_z() ) ), MBT.col_y() )) < 1e-2 );
		}

		{ // sin gamma version of params[3] (tilt) is a simple dot product:
			//Real const tmp1( sin( gamma ) * params[3] / gamma );
			//Real const tmp2( dot( Vector( cross( stub2.M.col_z(), stub1.M.col_z() )), MBT.col_x() ) );
			assert( std::fabs( sin( gamma ) * params[3] / gamma -
						 dot( Vector( cross( stub2.M.col_z(), stub1.M.col_z() )), MBT.col_x() )) < 1e-2 );
		}

		// check sign conventions
		Real phi_prime( arccos( dot( rt, MBT.col_y() ) ) );
		if ( dot( cross( rt, MBT.col_y() ), MBT.col_z() ) < 0.0f ) {
			phi_prime *= -1.0f;
		}

		Vector tmp( cross( M2.col_x(), M1.col_x() ) );
		assert( cross( tmp, MBT.col_z() ).length() <1e-2 );

		//Real const p1x = std::asin( dot( MBT.col_z(), cross( M2.col_y(), M1.col_y() ) ) );
		//Real const p1z = std::asin( dot( MBT.col_z(), cross( M2.col_x(), M1.col_x() ) ) );
		//std::cout << "equal? p1: " << params[1] << ' ' << p1x << ' ' << p1z <<
 		//	std::endl;
		assert( std::fabs( params[1] - asin( dot( MBT.col_z(), cross( M2.col_y(), M1.col_y() ) ) ) ) +
				std::fabs( params[1] - asin( dot( MBT.col_z(), cross( M2.col_x(), M1.col_x() ) ) ) ) < 1e-2 );

		//Real const p2 = gamma * cos( phi_prime );
		//Real const p3 = gamma * sin( phi_prime );
		//Real const dev( std::fabs( p2 - params[2] ) + std::fabs( p3 - params[3] ) );
		//std::cout << "dev: " << dev << std::endl;
		assert( std::fabs( gamma * cos( phi_prime ) - params[2] ) + std::fabs( gamma * sin( phi_prime ) - params[3] ) < 1e-2 );

		// check sign conventions
		assert( params[1] * dot( MBT.col_z(), cross( M2.col_y(), M1.col_y() ) ) > 0);
	}

	// convert to degrees
	params[1] = degrees( params[1] );
	params[2] = degrees( params[2] );
	params[3] = degrees( params[3] );
}


////////////////////////////////////////////////////////////////////////////
void
get_base_pucker(
	conformation::Residue const & rsd,
	std::pair< std::string, int > & pucker,
	Real & pseudorotation,
	Real & amplitude
)
{
	using numeric::conversions::radians;
	using numeric::conversions::degrees;

	utility::vector1< std::string > names;
	names.push_back( "C1'" );
	names.push_back( "C2'" );
	names.push_back( "C3'" );
	names.push_back( "C4'" );
	names.push_back( "O4'" );

	utility::vector1< Vector > atoms;
	for ( int i=1; i<= 5; ++i ) {
		atoms.push_back( rsd.xyz( names[i] ) );
	}

	Real mindot = 1000.0;
	bool exxo( false );
	utility::vector1< Real > torsions(5, 0.0);
	for ( int ii=1; ii<= 5; ++ii ) {

		torsions[ ii ] = dihedral_radians(
																			atoms[ 4 ],
																			atoms[ 5 ],
																			atoms[ 1 ],
																			atoms[ 2 ]
																			);

		Vector n12 = (( atoms[2]-atoms[1] ).cross( atoms[3]-atoms[2] ) ).normalized();
		Real dot = std::fabs( n12.dot( ( atoms[4]-atoms[3] ).normalized() ) );
		if ( dot < mindot ) {
			// get pucker
			//Real pucker_dot = n12.dot( ( atoms[5] - Real(0.5) * ( atoms[4] + atoms[1] ) ).normalized() );

			mindot = dot;
			pucker.first = names[5];
			exxo = ( n12.dot( ( atoms[5] - Real(0.5) * ( atoms[4] + atoms[1] ) ).normalized() ) > 0.0 );
		}

		atoms.push_back( atoms[1] );
		atoms.erase( atoms.begin() );

		names.push_back( names[1] );
		names.erase( names.begin() );

	}

	pseudorotation = atan(
													( ( torsions[2] + torsions[5] ) - ( torsions[1] + torsions[4] ) ) /
													( 2.0 * torsions[3] * ( sin( radians(36.0) ) + sin( radians(72.0)) ) )
												);

	pseudorotation = degrees( pseudorotation );
	if( torsions[3] < 0.0 ) pseudorotation += 180.0;

	pseudorotation = basic::periodic_range( pseudorotation, 360.0 );

	amplitude = degrees( torsions[3] / ( cos( radians( pseudorotation ) ) + 1.0e-20 ) );

	// additional integer for scannability
	{
		int const atom_index( std::find( names.begin(), names.end(), pucker.first ) - names.begin() );
		int const sign_index( exxo ? 0 : 1 );
		if ( atom_index%2 == sign_index ) pucker.second = atom_index+1;
		else                              pucker.second = atom_index-4;
	}

	if ( exxo ) pucker.first += " exxo";
	else pucker.first += " endo";
}

///////////////////////////////////////////////////////////////////////////////
kinematics::Stub
get_midstep_stub(
	kinematics::Stub const & in_stub1,
	kinematics::Stub const & in_stub2
)
{
	using numeric::conversions::degrees;
	using numeric::arccos;

	// Reordering as Phil did
	kinematics::Stub stub1( kinematics::Stub::Matrix::cols( in_stub2.M.col_y(), in_stub2.M.col_z(), in_stub2.M.col_x() ), in_stub2.v );
	kinematics::Stub stub2( kinematics::Stub::Matrix::cols( in_stub1.M.col_y(), in_stub1.M.col_z(), in_stub1.M.col_x() ), in_stub1.v );

	// copy matrices
	Matrix M1( stub1.M ), M2( stub2.M );

	assert( is_orthonormal( M1, 1e-3 ) );
	assert( is_orthonormal( M2, 1e-3 ) );

	//bool base_flipped = false;  // unused ~Labonte
	if ( dot( M1.col_z(), M2.col_z() ) < 0.0 ) {
		//base_flipped = true;  // unused ~Labonte
		basic::T("core.scoring.base_geometry") << "get_midstep_stub: base flip!!!\n";
		//utility::exit( EXIT_FAILURE, __FILE__, __LINE__);
	}

	// The nomenclature for the angles is the base pair convention (buckle, opening)
	// but the math is the same as for base steps, and the code in get_stub_stub_params works
	// for both.
	//
	// get angle between the y-axes
	Real const gamma( arccos( dot( M1.col_y(), M2.col_y() ) ) );

	Vector const bo( ( cross( M2.col_y(), M1.col_y() ) ).normalized() );

	Matrix R_gamma_2( rotation_matrix( bo, gamma/2.0f ) );

	M2 = R_gamma_2 * M2;
	M1 = R_gamma_2.transposed() * M1;

	assert( is_orthonormal( M1, 1e-3 ) );
	assert( is_orthonormal( M2, 1e-3 ) );

	// build mid-stub triad
	assert( M1.col_y().distance( M2.col_y() ) < 1e-3 );
	assert( std::fabs( dot( bo, M1.col_y() ) ) < 1e-3 );

	Matrix MBT;
	MBT.col_y( M1.col_y() );

	assert( std::fabs( dot( M1.col_z(), MBT.col_y() ) ) < 1e-3 );
	assert( std::fabs( dot( M2.col_z(), MBT.col_y() ) ) < 1e-3 );
	assert( std::fabs( dot( M1.col_x(), MBT.col_y() ) ) < 1e-3 );
	assert( std::fabs( dot( M2.col_x(), MBT.col_y() ) ) < 1e-3 );

	// get
	MBT.col_x( ( 0.5f * ( M1.col_x() + M2.col_x() ) ).normalized() );
	MBT.col_z( ( 0.5f * ( M1.col_z() + M2.col_z() ) ).normalized() );

	assert( is_orthonormal( MBT, 1e-3 ) );

	return kinematics::Stub( MBT, 0.5f*( stub1.v + stub2.v ) );

}



} // namespace dna
}} // scoring core
