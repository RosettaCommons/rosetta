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
/// @author Phil Bradley


#include <protocols/loops/loop_closure/ccd/ccd_closure.hh>

// Rosetta Headers
#include <core/types.hh>
#include <core/pose/Pose.hh>
#include <core/chemical/ResidueConnection.hh>
#include <core/conformation/Residue.hh>
#include <core/scoring/Ramachandran.hh>
#include <core/scoring/ScoringManager.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/kinematics/Stub.hh>
#include <basic/prof.hh>
#include <basic/basic.hh>
#include <basic/Tracer.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/format.hh>
#include <ObjexxFCL/string.functions.hh>

// C++ Headers

//Utility Headers
#include <numeric/xyz.functions.hh>
#include <numeric/numeric.functions.hh>
#include <numeric/random/random.fwd.hh>
#include <utility/assert.hh>

#include <utility/vector1.hh>

//Auto Headers



using namespace ObjexxFCL;
using namespace ObjexxFCL::format;

namespace protocols {
namespace loops {
namespace loop_closure {
namespace ccd {

using namespace core;
using utility::vector1;
typedef numeric::xyzMatrix< Real > Matrix;

basic::Tracer tccd( "protocols.loops.loop_closure.ccd.ccd_closure" );

/////////////////////////////////////////////////////////////////////////////////
/// @brief copy mainchain atoms xyz and torsions from pose to coords and torsions
void
load_coords_and_torsions(
	pose::Pose const & pose,
	vector1< vector1< Vector > > & coords,
	vector1< vector1< Real > > & torsions
)
{
	Size const nres( pose.total_residue() );
	coords.resize( nres );
	torsions.resize( nres );

	for ( Size i=1; i<= nres; ++i ) {
		conformation::Residue const & rsd( pose.residue(i) );
		Size const nbb( rsd.mainchain_atoms().size() );
		if ( nbb ) {
			coords[i].resize( nbb );
			torsions[i].resize( nbb );
			for ( Size j=1; j <= nbb; ++j ) {
				coords[i][j] = rsd.atom( rsd.mainchain_atoms()[j] ).xyz();
				torsions[i][j] = rsd.mainchain_torsion(j);
			}
		}
	}

}
/////////////////////////////////////////////////////////////////////////////
///@brief copy torsions into pose between residue loop_begin and loop_end
void
copy_torsions_to_pose(
	pose::Pose & pose,
	int const loop_begin,
	int const loop_end,
	vector1< vector1< Real > > const & torsions
)
{
	for ( int i=loop_begin; i<= loop_end; ++i ) {
		for ( Size j=1; j<= torsions[i].size(); ++j ) {
			//TorsionID id( i, id::BB, int(j) ), torsions[i][j] );
			pose.set_torsion( id::TorsionID( i, id::BB, int(j) ), torsions[i][j] );
		}
	}
}

///////////////////////////////////////////////////////////////////////////

void
get_overlap_pos(
	vector1< vector1< Vector > > const & coords,
	vector1< vector1< Real > > const & torsions,
	int const cutpoint,
	int const direction,
	Real const bond_angle1,
	Real const bond_length,
	Real const bond_angle2,
	Matrix & M
)
{
	using numeric::conversions::radians;
	using kinematics::Stub;

	bool const local_debug( false/*true*/ );

	Size const nbb( coords[cutpoint].size() );
	assert( nbb >= 3 );

	// pack the coords, torsions, and angles
	vector1< Vector > atoms( 6 );
	vector1< Real > dihedrals( 3 );
	vector1< Real > angles( 2 );

	if ( direction == 1 ) { // eg N->CA->C->N'->CA'->C'
		for ( int i=1; i<= 3; ++i ) atoms[i  ] = coords[ cutpoint   ][ nbb-3+i ]; // 1,2,3: N, CA, C
		for ( int i=1; i<= 3; ++i ) atoms[i+3] = coords[ cutpoint+1 ][ i       ]; // 4,5,6: N',CA',C'
		dihedrals[1] = radians( torsions[cutpoint  ][ nbb-1 ] ); // N-CA-C-N'
		dihedrals[2] = radians( torsions[cutpoint  ][ nbb   ] ); // CA-C-N'-CA'
		dihedrals[3] = radians( torsions[cutpoint+1][ 1     ] ); // C-N'-CA'-C'
		angles[1] = bond_angle1; // CA-C-N' already in radians
		angles[2] = bond_angle2; // C-N'-CA' already in radian
	} else { // eg N<-CA<-C<-N'<-CA'<-C'
		for ( int i=1; i<= 3; ++i ) atoms[7-i  ] = coords[ cutpoint   ][ nbb-3+i ]; // 6,5,4: N, CA, C
		for ( int i=1; i<= 3; ++i ) atoms[7-i-3] = coords[ cutpoint+1 ][ i       ]; // 3,2,1: N',CA',C'
		dihedrals[3] = radians( torsions[cutpoint  ][ nbb-1 ] ); // N-CA-C-N'
		dihedrals[2] = radians( torsions[cutpoint  ][ nbb   ] ); // CA-C-N'-CA'
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

		// atom3p is built backward, eg given N'-CA'-C, how to build C
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
		ASSERT_ONLY( Real const dihedral1
			( dihedral_radians( atoms[1], atoms[2], atoms[3], M.col(1) ) );)
		ASSERT_ONLY(Real const dihedral2
			( dihedral_radians( atoms[2], atoms[3], M.col(1), M.col(2) ) );)
		ASSERT_ONLY(Real const dihedral3
			( dihedral_radians( atoms[3], M.col(1), M.col(2), M.col(3) ) );)

		ASSERT_ONLY(Real const angle1( std::acos( dot( ( atoms[3] - atoms[2] ).normalized(),
																									 ( M.col(1) - atoms[3] ).normalized() ) ) );)
		ASSERT_ONLY(Real const angle2( std::acos( dot( ( M.col(2) - M.col(1) ).normalized(),
																									 ( M.col(1) - atoms[3] ).normalized() ) ) );)
		ASSERT_ONLY(Real const length( atoms[3].distance( M.col(1) ) );)

		assert( std::abs( subtract_radian_angles( dihedral1, dihedrals[1] ) ) < 1e-3 );
		assert( std::abs( subtract_radian_angles( dihedral2, dihedrals[2] ) ) < 1e-3 );
		assert( std::abs( subtract_radian_angles( dihedral3, dihedrals[3] ) ) < 1e-3 );
		assert( std::abs( subtract_radian_angles( angle1, angles[1] ) ) < 1e-3 );
		assert( std::abs( subtract_radian_angles( angle2, angles[2] ) ) < 1e-3 );
		assert( std::abs( length - bond_length ) < 1e-3 );
	}

	// restore proper chain order for M if folding backward
	if ( direction == -1 ) {
		Matrix const tmp( M );
		for ( int i=1; i<= 3; ++i ) M.col( i, tmp.col(4-i) );
	}
}


///////////////////////////////////////////////////////////////////////////
/// @brief check whether xyz coords of overlap position match internal coords
void
check_overlap_pos(
	vector1< vector1< Vector > > const & coords,
	vector1< vector1< Real > > const & torsions,
	int const cutpoint,
	int const direction,
	Real const bond_angle1,
	Real const ASSERT_ONLY(bond_length),
	Real const bond_angle2,
	Matrix const & M_in
)
{
	using numeric::conversions::radians;
	using kinematics::Stub;

	Size const nbb( coords[cutpoint].size() );
	assert( nbb >= 3 );

	// pack the coords, torsions, and angles
	vector1< Vector > atoms( 6 );
	vector1< Real > dihedrals( 3 );
	vector1< Real > angles( 3 );

	if ( direction == 1 ) {
		for ( int i=1; i<= 3; ++i ) atoms[i  ] = coords[ cutpoint   ][ nbb-3+i ];
		for ( int i=1; i<= 3; ++i ) atoms[i+3] = coords[ cutpoint+1 ][ i       ];
		dihedrals[1] = radians( torsions[cutpoint  ][ nbb-1 ] );
		dihedrals[2] = radians( torsions[cutpoint  ][ nbb   ] );
		dihedrals[3] = radians( torsions[cutpoint+1][ 1     ] );
		angles[1] = bond_angle1; // already in radians
		angles[2] = bond_angle2;
	} else {
		for ( int i=1; i<= 3; ++i ) atoms[7-i  ] = coords[ cutpoint   ][ nbb-3+i ];
		for ( int i=1; i<= 3; ++i ) atoms[7-i-3] = coords[ cutpoint+1 ][ i       ];
		dihedrals[3] = radians( torsions[cutpoint  ][ nbb-1 ] );
		dihedrals[2] = radians( torsions[cutpoint  ][ nbb   ] );
		dihedrals[1] = radians( torsions[cutpoint+1][ 1     ] );
		angles[2] = bond_angle1; // already in radians
		angles[1] = bond_angle2;
	}

	Matrix M( M_in );

	if ( direction == -1 ) {
		for ( int i=1; i<= 3; ++i ) M.col( i, M_in.col(4-i) );
	}

	{
		using basic::subtract_radian_angles;
		ASSERT_ONLY(Real const dihedral1
			( dihedral_radians( atoms[1], atoms[2], atoms[3], M.col(1) ) );)
		ASSERT_ONLY(Real const dihedral2
			( dihedral_radians( atoms[2], atoms[3], M.col(1), M.col(2) ) );)
		ASSERT_ONLY(Real const dihedral3
			( dihedral_radians( atoms[3], M.col(1), M.col(2), M.col(3) ) );)

		ASSERT_ONLY(Real const angle1( std::acos( dot( ( atoms[3] - atoms[2] ).normalized(),
																									 ( M.col(1) - atoms[3] ).normalized() ) ) );)
		ASSERT_ONLY(Real const angle2( std::acos( dot( ( M.col(2) - M.col(1) ).normalized(),
																									 ( M.col(1) - atoms[3] ).normalized() ) ) );)
		ASSERT_ONLY(Real const length( atoms[3].distance( M.col(1) ) );)

		assert( std::abs( subtract_radian_angles( dihedral1, dihedrals[1] ) ) < 1e-3 );
		assert( std::abs( subtract_radian_angles( dihedral2, dihedrals[2] ) ) < 1e-3 );
		assert( std::abs( subtract_radian_angles( dihedral3, dihedrals[3] ) ) < 1e-3 );
		assert( std::abs( subtract_radian_angles( angle1, angles[1] ) ) < 1e-3 );
		assert( std::abs( subtract_radian_angles( angle2, angles[2] ) ) < 1e-3 );
		assert( std::abs( length - bond_length ) < 1e-3 );
	}

}


///////////////////////////////////////////////////////////////////////////
void
index_pair_in_range( int & pos, int & atom, int const nbb ) {
	while ( atom > nbb ) {
		atom -= nbb;
		pos += 1;
	}
	while ( atom < 1 ) {
		atom += nbb;
		pos -= 1;
	}
}

//////////////////////////////////////////////////////////////////////////
void
get_torsion_axis(
	vector1< vector1< Vector > > const & coords,
	int const seqpos,
	int const torsion,
	Vector & axis_atom, // output
	Vector & axis_vector // output
)
{
	int const nbb( coords[seqpos].size() );
	int pos1( seqpos ), pos2( seqpos ),
		atom1( torsion ), atom2( torsion+1 );
	index_pair_in_range( pos2, atom2, nbb ); // in case torsion == nbb

	axis_atom = coords[pos2][atom2];
	axis_vector = ( axis_atom - coords[pos1][atom1] ).normalized();
}


//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
// this routine is used by ccd loop closure

void
refold_loop_torsion(
	Real const alpha, // degrees
	int const pos,
	int const torsion,
	int const dir,
	int const cutpoint,
	vector1< vector1< Vector > > & coords,
	Matrix & M
)
{
	using numeric::conversions::radians;
	Size const nbb( coords[pos].size() );
	int const fold_end( ( dir == 1 ) ? cutpoint : cutpoint + 1 );

	assert( ( fold_end - pos ) * dir >= 0 );

	// get the rotation matrix, offset
	Vector axis_atom, axis_vector;
	get_torsion_axis( coords, pos, torsion, axis_atom, axis_vector );

	Matrix const R( numeric::rotation_matrix_degrees( axis_vector, Real( dir * alpha) ) );

	Vector const v( axis_atom - R * axis_atom );
	// get the rotation matrix about this axis_vector
	// A is an atom on the torsion axis (axis_atom)
	//
	// transformation is x --> R * ( x - A ) + A = R * x - R*A + A = R * x + v where v = ( A - R * A )
	//

	int fold_begin( pos );
	int fold_begin_atom( ( dir == 1 ) ? torsion + 2 : torsion - 1 );
	index_pair_in_range( fold_begin, fold_begin_atom, nbb );

	//
	for ( int i = fold_begin; ( fold_end - i ) * dir >= 0; i += dir ) {
		int j_begin( ( i == fold_begin ) ? fold_begin_atom : ( dir == 1 ? 1 : nbb ) );

		for ( Size j = j_begin; j >= 1 && j <= nbb; j += dir ) {
			coords[i][j] = R * coords[i][j] + v;
		}
	}

	// move the overlap pos
	for ( int i=1; i<= 3; ++i ) {
		M.col(i, R*M.col(i) + v);
	}
}

/////////////////////////////////////////////////////////////////////////////
void
calculate_ccd_angle(
	Matrix const & F,
	Matrix const & M,
	vector1< vector1< Vector > > const & coords,
	int const pos,
	int const torsion,
	int const direction,
	Real & angle,
	Real & dev
)
{
	using numeric::conversions::radians;
	using numeric::conversions::degrees;


// !!!!!!!!!! REAL PRECISION !!!!!!!!!!!!!!!!!!!
// Reald letters are scalars: aa,bb,cc,dd
	//Real r_length, f_length, alpha;

	int const n2c = { 1 };
	int const c2n = { -1 };

	bool const local_debug = { false };//{ true };

	Vector A, axis_vector;
	//Vector A, axis_vector, O, r, r_hat, s_hat, f_vector;
	get_torsion_axis( coords, pos, torsion, A, axis_vector );

// accumulate the components of the ccd equation
	Real aa = 0.0;
	Real bb = 0.0;
	Real cc = 0.0;

	for ( int i = 1; i <= 3; ++i ) {
		// project M onto the axis ==> O

		//AM = M.col(i) - A;
		//dd = dot( AM, axis_vector );
		//v = dd * axis_vector;
		Vector const O( A + dot( M.col(i) - A, axis_vector ) * axis_vector );

		// calculate r
		Vector const r  = M.col(i) - O;
		Real const r_length = r.length();
		// calculate r_hat
		Vector const r_hat = r / r_length;

		// calculate s_hat; orientation choice comes here!!
		Vector const s_hat = cross( axis_vector, r_hat );

		// calculate f_vector
		Vector const f_vector = F.col(i) - O;
		Real const f_length = f_vector.length();

		aa += r_length * r_length + f_length * f_length;

		bb += 2 * r_length * dot( f_vector, r_hat );

		cc += 2 * r_length * dot( f_vector, s_hat );

		if ( local_debug ) {
			// some checks on orthonormality
			assert( std::abs( dot( axis_vector, r_hat ) ) < 1e-3 );
			assert( std::abs( dot( axis_vector, s_hat ) ) < 1e-3 );
			assert( std::abs( dot( s_hat, r_hat ) ) < 1e-3 );
			assert( axis_vector.is_unit( 1e-3 ) );
			assert( s_hat.is_unit( 1e-3 ) );
			assert( r_hat.is_unit( 1e-3 ) );

			assert( std::abs( dot( M.col(i), axis_vector ) - dot( O, axis_vector ) ) < 1e-3 );
			assert( std::abs( dot( r, axis_vector ) ) < 1e-3 );

		}
	}

	Real const alpha = std::atan2( cc, bb );

	dev = aa - bb * std::cos( alpha ) - cc * std::sin( alpha );

	if ( local_debug ) {
		Real const one_degree = radians( 1.0 ); // one degree expressed in radians
		if ( ( dev > aa - bb * std::cos( alpha - one_degree ) - cc * std::sin( alpha - one_degree ) ) ||
				 ( dev > aa - bb * std::cos( alpha + one_degree ) - cc * std::sin( alpha + one_degree ) ) ) {
			utility_exit_with_message( "ccd problemo!" );
		}
	}

	if ( direction == n2c ) {
		angle = degrees(alpha);
	} else if ( direction == c2n ) {
		angle = -degrees(alpha);
	} else {
		utility_exit_with_message( "ccd_angle: unrecognized direction: "+string_of(direction) );
	}

}

///////////////////////////////////////////////////////////////////////////////////
void
check_torsions(
							 int const loop_begin,
							 int const loop_end,
							 int const cutpoint,
							 vector1< vector1< Real > > const & torsions,
							 vector1< vector1< Vector > > const & coords
							 )
{
	int const nbb( torsions[ loop_begin ].size() );
	for ( int i=loop_begin; i<= loop_end; ++i ) {
		for ( int torsion=1; torsion<= nbb; ++torsion ) {
			vector1< Vector > xyz;
			for ( int k=1; k<= 4; ++k ) {
				int seqpos(i), atomno(torsion+k-2);
				index_pair_in_range( seqpos, atomno, nbb );
				xyz.push_back( coords[seqpos][atomno] );
			}
			Real const xyz_dihedral( numeric::dihedral_degrees( xyz[1], xyz[2], xyz[3], xyz[4] ) );
			if ( ( i == cutpoint && ( torsion == nbb -1 || torsion == nbb ) ) ||
					 ( i == cutpoint+1 && ( torsion == 1 ) ) ) continue;
			assert( std::abs( basic::subtract_degree_angles( xyz_dihedral, torsions[i][torsion] ) ) < 1e-3 );
			tccd.Debug << "torsion-check: " << i << ' ' << torsion << ' ' << torsions[i][torsion] << ' ' << xyz_dihedral << std::endl;
		}

	}
}


//------------------------------------------------------------------------------
//
/// @details tolerance in Angstroms, forward and backward splice RMS over N,CA,C must
/// be less than tolerance for an early return. otherwise goes through the
/// loop ii_cycles, each time both forward and backward
/// returns the number of cycles actually taken

int
fast_ccd_loop_closure(
	pose::Pose & pose,
	kinematics::MoveMap const & mm,
	int const loop_begin,
	int const loop_end,
	int const cutpoint,
	int const ii_cycles,
	Real const tolerance,
	bool const rama_check,
	Real const max_rama_score_increase,
	Real const max_total_delta_helix,
	Real const max_total_delta_strand,
	Real const max_total_delta_loop,
	Real & forward_deviation, // output
	Real & backward_deviation, // output
	Real & torsion_delta,
	Real & rama_delta
)
{
	PROF_START( basic::CCD_CLOSE );
	using namespace id;

	/////////
	// params
	/////////
	bool const local_debug = { false }; //true };
	bool const local_verbose = { false }; //true };
	int const n2c = { 1 }; // must be 1 and -1 (below) for proper incrementing in loops
	int const c2n = { -1 };
	// per-cycle max move params
	//Real const max_angle_delta_helix = { 100.0 };
	//Real const max_angle_delta_strand = { 100.0 };
	//Real const max_angle_delta_loop = { 100.0 };
 	Real const max_angle_delta_helix = { 1.0 };
 	Real const max_angle_delta_strand = { 5.0 };
 	Real const max_angle_delta_loop = { 10.0 };
	// for rama-checking with Boltzman criterion (should be higher?)
	Real const ccd_temperature = { 0.25 };
	Size const nres( pose.total_residue() );

	Real const bond_angle1( pose.residue( cutpoint ).upper_connect().icoor().theta() );// CA-C=N bond angle
	Real const bond_angle2( pose.residue( cutpoint+1 ).lower_connect().icoor().theta() ); // C=N-CA bond angle
	Real const bond_length( pose.residue( cutpoint+1 ).lower_connect().icoor().d() ); // C=N distance

	// pack mainchain coordinates, torsions into local arrays
	Size const nbb( pose.residue( loop_begin ).mainchain_atoms().size() );
	vector1< vector1< Vector > > coords;
	vector1< vector1< Real > > torsions;

	load_coords_and_torsions( pose, coords, torsions );
	assert( coords[loop_begin].size() == nbb && torsions[loop_begin].size() == nbb );

	// save the starting torsions for comparison, starting rama score if rama_check == true
	vector1< vector1< Real > > start_torsions( torsions );
	vector1< Real > start_rama_score( nres, 0.0 );
	scoring::RamachandranCOP rama(0);
	if ( rama_check ) {
		rama = &(scoring::ScoringManager::get_instance()->get_Ramachandran());
		for ( int i=loop_begin; i<= loop_end; ++i ) {
			start_rama_score[i] = rama->eval_rama_score_residue( pose.aa(i), torsions[i][1], torsions[i][2] );
		}
	}

	// reuse this for saving coords in case of move-rejection to prevent unnecessary refolds
	vector1< vector1< Vector > > save_coords( coords );

	// for reporting how many cycles we actually took
	int actual_cycles = ii_cycles;


	////////////////////////////
	// now start cycling through
	////////////////////////////
	for ( int ii = 1; ii <= ii_cycles; ++ii ) {

		if ( ii%10 == 0 ) {
			// avoid coordinate errors due to roundoff
			copy_torsions_to_pose( pose, loop_begin, loop_end, torsions );
			load_coords_and_torsions( pose, coords, torsions );
		}

		// first forward, then backward
		for ( int repeat = 1; repeat <= 2; ++repeat ) {
			int direction,start_pos,stop_pos,increment;
			if ( repeat == 1 ) {
				direction = n2c;
				start_pos = loop_begin;
				stop_pos = std::min(loop_end,cutpoint);
				increment = 1;
			} else {
				direction = c2n;
				start_pos = loop_end;
				stop_pos = std::max(loop_begin,cutpoint+1);
				increment = -1;
			}

			// get fixed and moving positions for this direction:
			Matrix F,M;
			// F is the fixed landing position,
			// for n->c direction, that would be the first three backbone atom of the residue after cutpoint
			// for c->n direction, that would be the last three backbone atom of the residue of cutpoint
			for ( int i = 1; i <= 3; ++i ) {
				if ( direction == n2c ) {
					F.col( i, coords[ cutpoint + 1 ][ i           ] );
				} else {
					F.col( i, coords[ cutpoint     ][ nbb - 3 + i ] );
				}
			}
			// M is the actual landing postion
			get_overlap_pos( coords, torsions, cutpoint, direction, bond_angle1, bond_length, bond_angle2, M );

			if ( local_debug )
				check_overlap_pos( coords, torsions, cutpoint, direction, bond_angle1, bond_length, bond_angle2, M );



			// as we change torsion angles through the loop the moving
			// atoms in M will be transformed along with the downstream segment
			// of the loop, by the routine refold_loop_torsions(...)

			for ( int pos = start_pos; pos*increment <= stop_pos*increment;
						pos += increment ) {

				if ( !mm.get( TorsionID( pos, BB, phi_torsion ) ) && !mm.get( TorsionID( pos, BB, psi_torsion ) ) ) continue;


				Real max_angle_delta, max_total_delta;
				if ( pose.secstruct(pos) == 'H' ) {
					max_angle_delta = max_angle_delta_helix;
					max_total_delta = max_total_delta_helix;
				} else if ( pose.secstruct(pos) == 'E' ) {
					max_angle_delta = max_angle_delta_strand;
					max_total_delta = max_total_delta_strand;
				} else {
					max_angle_delta = max_angle_delta_loop;
					max_total_delta = max_total_delta_loop;
				}

				if ( max_total_delta <= 0.01 ) {
					tccd.Debug << "cant move this residue " << pos << std::endl;
					continue;
				}

				// save values for rama score eval, and in case we rama-reject
				vector1< Real > const save_torsions( torsions[ pos ] );

				// save in case we rama-reject!
				save_coords.resize( nres );
				for ( int i=loop_begin; i<= loop_end; ++i ) {
					save_coords[i].resize( nbb );
					for ( Size j=1; j<= nbb; ++j ) {
						save_coords[i][j] = coords[i][j];
					}
				}

				for ( Size torsion = 1; torsion <= nbb; ++torsion ) {

					if ( !mm.get( TorsionID( pos, BB, torsion ) ) ) continue;

					if ( pos == cutpoint && torsion == nbb ) continue;

					if ( local_debug )
						check_overlap_pos( coords, torsions, cutpoint, direction, bond_angle1, bond_length, bond_angle2, M );

					Real alpha,dev;
					calculate_ccd_angle( F, M, coords, pos, torsion, direction, alpha, dev );
					Real const alpha_orig( alpha );

					// impose per-move deltas
					if ( alpha >  max_angle_delta ) alpha =  max_angle_delta;
					if ( alpha < -max_angle_delta ) alpha = -max_angle_delta;

					// check for total movement during closure run:
					Real const total_delta
						( basic::subtract_degree_angles( start_torsions[pos][torsion], torsions[pos][torsion] + alpha ) );

					if ( total_delta > max_total_delta ) {
						// this logic is a little tricky: if adding alpha to the previous
						// delta pushes us past 180 from start, then it wont work, so check for
						// that (note that if max_total_delta > 180 we wont even get here ):
						assert( alpha + max_total_delta < 180 );
						if ( alpha > 0 ) {
							alpha -= ( total_delta - max_total_delta + 0.01 );
						} else {
							alpha += ( total_delta - max_total_delta + 0.01 );
						}
					}

					// update the coordinates to reflect the torsion angle change
					refold_loop_torsion( alpha, pos, torsion, direction, cutpoint,  coords, M );

					// update torsions
					torsions[ pos ][ torsion ] = basic::periodic_range( torsions[pos][torsion] + alpha, 360.0 );

					if ( local_debug )
						check_overlap_pos( coords, torsions, cutpoint, direction, bond_angle1, bond_length, bond_angle2, M );

					// debugging: ///////////////////////
					if ( local_debug ) {
						// this assumes that the tree is not rooted at cutpoint or
						// cutpoint+1
						check_torsions( loop_begin, loop_end, cutpoint, torsions, coords );

						copy_torsions_to_pose( pose, loop_begin, loop_end, torsions );

						vector1< vector1< Vector > > tmp_coords;
						vector1< vector1< Real > > tmp_torsions;

						load_coords_and_torsions( pose, tmp_coords, tmp_torsions );

						// calculate coordinate dev
						Real coord_dev = 0.0;
						for ( int i = loop_begin; i <= loop_end; ++i ) {
							for ( Size j = 1; j <= nbb; ++j ) {
								coord_dev += coords[i][j].distance( tmp_coords[i][j] );
							}
						}

						if ( coord_dev > 0.1 ) {
							utility_exit_with_message( "fold_loop_torsion dev:: " + string_of( coord_dev ) );
						}

						Real tmp_dev = 0.0;
						for ( int i = 1; i <= 3; ++i ) {
							tmp_dev += M.col(i).distance_squared( F.col(i) );
						}

						if ( alpha == alpha_orig && std::abs( dev - tmp_dev) > 0.1 ) {
							tccd.Warning << "WARNING:: ccd-dev error: " << dev << ' ' <<
								tmp_dev << ' ' << std::abs( dev-tmp_dev) << std::endl;
						}

						if ( local_verbose )
							tccd.Warning << "pos: " << I( 3, pos ) << I( 3, torsion) <<
								" alpha-orig: " << format::F( 7, 3, alpha_orig ) <<
								" alpha: " << format::F( 7, 3, alpha ) <<
								" dev1: " << format::F( 13, 9, dev ) <<
								" dev2: " << format::F( 13, 9, tmp_dev ) << std::endl;

					} // if ( local_debug ) ///////////////////////////////////
				} // torsion = 1,nbb


				if ( rama_check ) {
					assert( nbb == 3 );
					//////////////////////////////////////
					// evaluate the rama score difference:
					Real const old_rama_score( rama->eval_rama_score_residue( pose.aa(pos), torsions[pos][1], torsions[pos][2]));
					Real const new_rama_score( rama->eval_rama_score_residue( pose.aa(pos), save_torsions[1], save_torsions[2]));

					if (local_verbose)
						tccd.Warning  << "rama_check: " << pos << ' ' << old_rama_score << ' ' << new_rama_score << std::endl;

					if ( new_rama_score > old_rama_score ) {
						Real const boltz_factor ( (old_rama_score-new_rama_score)/ccd_temperature );
						Real const negative_forty( -40.0f );
						Real const probability ( std::exp(std::max(negative_forty,boltz_factor) ) );
						if ( new_rama_score - start_rama_score[ pos ] > max_rama_score_increase ||
							numeric::random::uniform() >= probability ) {
							// undo the change:
							torsions[pos] = save_torsions;

							// recover saved coords
							for ( int i=loop_begin; i<= loop_end; ++i ) {
								for ( Size j=1; j<= nbb; ++j ) {
									coords[i][j] = save_coords[i][j];
								}
							}

							get_overlap_pos( coords, torsions, cutpoint, direction, bond_angle1, bond_length, bond_angle2, M );

						} // rama reject
					} // rama-score got worse
				} // if ( rama_check )
			} // pos = start_pos,stop_pos
		} // repeat = 1,2   1=n2c; 2=c2n

		if ( ii%5 == 0 || ii == ii_cycles ) {
			// check overlap deviations to see if loop is closed
			// every 5 cycles or on the last one.
			//
			// forward_deviation:
			Matrix F,M;
			int direction = n2c;
			for ( int i = 1; i <= 3; ++i ) {
				F.col( i, coords[ cutpoint + 1 ][ i ] );
			}

			get_overlap_pos( coords, torsions, cutpoint, direction, bond_angle1, bond_length, bond_angle2, M );

			forward_deviation = 0.0;
			for ( int i = 1; i <= 3; ++i ) {
				for ( int j = 1; j <= 3; ++j ) {
					forward_deviation += numeric::square( M(j,i) - F(j,i) );
				}
			}
			forward_deviation = sqrt( forward_deviation / 3 );

			// backward_deviation:
			direction = c2n;
			for ( int i = 1; i <= 3; ++i ) {
				F.col( i, coords[ cutpoint     ][ nbb - 3 + i ] );
			}
			get_overlap_pos( coords, torsions, cutpoint, direction, bond_angle1, bond_length, bond_angle2, M );

			backward_deviation = 0.0;
			for ( int i = 1; i <= 3; ++i ) {
				for ( int j = 1; j <= 3; ++j ) {
					backward_deviation += numeric::square( M(j,i) - F(j,i) );
				}
			}
			backward_deviation = sqrt( backward_deviation / 3 );

			if (local_verbose) tccd.Debug << "cycle/forward_dev/backward_dev: " << ii << " " << forward_deviation << " "
																	 << backward_deviation << std::endl;

			if ( forward_deviation < tolerance && backward_deviation < tolerance ) {
				if ( local_verbose )
					tccd.Debug << "closed early: ii= " << ii << ' ' << tolerance	<< std::endl;
				actual_cycles = ii;
				break;
			}
		}                  // check deviations
	}                     // ii=1,ii_cycles

	// DONE!
	copy_torsions_to_pose( pose, loop_begin, loop_end, torsions );

	// calculate torsion dev, rama delta
	torsion_delta = 0.0;
	rama_delta = 0.0;

	for ( int i=loop_begin; i<= loop_end; ++i ) {
		for ( Size j=1; j<= nbb; ++j ) {
			torsion_delta += std::abs( basic::subtract_degree_angles( start_torsions[i][j], torsions[i][j] ) );
		}
		if ( rama_check ) {
			Real const final_rama_score( rama->eval_rama_score_residue( pose.aa(i), torsions[i][1], torsions[i][2] ) );
			rama_delta += ( final_rama_score - start_rama_score[i] );
		}
	}
	torsion_delta /= ( loop_end - loop_begin + 1);
	rama_delta /= ( loop_end - loop_begin + 1);

	PROF_STOP( basic::CCD_CLOSE );
	return actual_cycles;
}

/////////////////////////////////////////////////////////
void
ccd_moves(
	int const total_moves,
	core::pose::Pose & pose,
	core::kinematics::MoveMap const & mm,
	int const loop_begin,
	int const loop_end,
	int const cutpoint
)
{


	using namespace id;

	/////////
	// params
	/////////
	bool const local_debug = { false }; //true };
	bool const local_verbose = { false }; //true };
	bool const rama_check = { true };

	int const n2c = { 1 }; // must be 1 and -1 (below) for proper incrementing in loops
	int const c2n = { -1 };

	// per-cycle max move params
	//Real const max_angle_delta_helix = { 100.0 };
	//Real const max_angle_delta_strand = { 100.0 };
	//Real const max_angle_delta_loop = { 100.0 };
 	Real const max_angle_delta_helix = { 1.0 };
 	Real const max_angle_delta_strand = { 5.0 };
 	Real const max_angle_delta_loop = { 10.0 };
	// for rama-checking with Boltzman criterion (should be higher?)
	Real const ccd_temperature = { 0.25 };
	Size const nres( pose.total_residue() );

	Real const bond_angle1( pose.residue( cutpoint ).upper_connect().icoor().theta() );// CA-C=N bond angle
	Real const bond_angle2( pose.residue( cutpoint+1 ).lower_connect().icoor().theta() ); // C=N-CA bond angle
	Real const bond_length( pose.residue( cutpoint+1 ).lower_connect().icoor().d() ); // C=N distance

	// pack mainchain coordinates, torsions into local arrays
	Size const nbb( pose.residue( loop_begin ).mainchain_atoms().size() );
	vector1< vector1< Vector > > coords;
	vector1< vector1< Real > > torsions;

	load_coords_and_torsions( pose, coords, torsions );
	assert( coords[loop_begin].size() == nbb && torsions[loop_begin].size() == nbb );

	// save the starting torsions for comparison,
	vector1< vector1< Real > > start_torsions( torsions );

	// reuse this for saving coords in case of move-rejection to prevent unnecessary refolds
	vector1< vector1< Vector > > save_coords( coords );

	Size total_insert( loop_end - loop_begin + 1 );
	Real const H_weight ( 0.5 );
	Real const E_weight ( 1.0 );
	Real const L_weight ( 8.5 );
	Real total_weight(0.0);
	vector1< Real > weight_map( total_insert, 0.0 );

	for ( Size i = 1; i <= total_insert; ++i ) {

			char ss( pose.secstruct( loop_begin - 1 + i ) );

			if ( ss == 'H' ) {
				total_weight += H_weight;
			} else if (  ss == 'E' ) {
				total_weight += E_weight;
			} else {
				assert( ss == 'L' );
				total_weight += L_weight;
			}
			weight_map[ i ] = total_weight;

	}


	Size nmoves( 0 );
	Size ntries( 0 );

	while ( nmoves < Size ( total_moves ) ) {
		++ntries;
		if ( ntries > Size ( 5*total_moves ) ) {
			break;
		}

		Real const weight( numeric::random::uniform()*total_weight );
		Size ipos;
		for( ipos = 1; ipos <= total_insert; ++ipos ) {
			if( weight < weight_map[ ipos ] ) break;
		}
		if( ipos > total_insert ) {
			tccd.Debug << "ccd_moves: bad_weight " << weight << ' ' <<
				total_weight << ' ' << weight_map[ total_insert ] << std::endl;
			ipos = total_insert;
		}

		Size const pos = ipos + loop_begin - 1;
		int const direction( pos <= Size( cutpoint ) ? n2c : c2n );


		// get fixed and moving positions for this direction:
		Matrix F,M;
		// F is the fixed landing position,
		// for n->c direction, that would be the first three backbone atom of the residue after cutpoint
		// for c->n direction, that would be the last three backbone atom of the residue of cutpoint
		for ( int i = 1; i <= 3; ++i ) {
			if ( direction == n2c ) {
				F.col( i, coords[ cutpoint + 1 ][ i           ] );
			} else {
				F.col( i, coords[ cutpoint     ][ nbb - 3 + i ] );
			}
		}
		// M is the actual landing postion
			get_overlap_pos( coords, torsions, cutpoint, direction, bond_angle1, bond_length, bond_angle2, M );

			if ( local_debug )
				check_overlap_pos( coords, torsions, cutpoint, direction, bond_angle1, bond_length, bond_angle2, M );


			Real max_angle_delta;
			if ( pose.secstruct(pos) == 'H' ) {
				max_angle_delta = max_angle_delta_helix;
			} else if ( pose.secstruct(pos) == 'E' ) {
				max_angle_delta = max_angle_delta_strand;
			} else {
				max_angle_delta = max_angle_delta_loop;
			}

			// save values for rama score eval, and in case we rama-reject
			vector1< Real > const save_torsions( torsions[ pos ] );

			// save in case we rama-reject!
			save_coords.resize( nres );
			for ( int i=loop_begin; i<= loop_end; ++i ) {
				save_coords[i].resize( nbb );
				for ( Size j=1; j<= nbb; ++j ) {
					save_coords[i][j] = coords[i][j];
				}
			}

			for ( Size torsion = 1; torsion <= nbb; ++torsion ) {

				if ( !mm.get( TorsionID( pos, BB, torsion ) ) ) continue;

				if ( pos == Size( cutpoint ) && torsion == nbb ) continue;

				if ( local_debug )
					check_overlap_pos( coords, torsions, cutpoint, direction, bond_angle1, bond_length, bond_angle2, M );

				Real alpha, dev;
				calculate_ccd_angle( F, M, coords, pos, torsion, direction, alpha, dev );
				//v				Real const alpha_orig( alpha );

				// impose per-move deltas
				if ( alpha >  max_angle_delta ) alpha =  max_angle_delta;
				if ( alpha < -max_angle_delta ) alpha = -max_angle_delta;

				// update the coordinates to reflect the torsion angle change
				refold_loop_torsion( alpha, pos, torsion, direction, cutpoint,  coords, M );

				// update torsions
				torsions[ pos ][ torsion ] = basic::periodic_range( torsions[pos][torsion] + alpha, 360.0 );

				if ( local_debug )
					check_overlap_pos( coords, torsions, cutpoint, direction, bond_angle1, bond_length, bond_angle2, M );


			} // loop over nbb


			scoring::RamachandranCOP rama(0);
			if ( rama_check ) {
				rama = &(scoring::ScoringManager::get_instance()->get_Ramachandran());
				assert( nbb == 3 );
				//////////////////////////////////////
				// evaluate the rama score difference:
				Real const old_rama_score( rama->eval_rama_score_residue( pose.aa(pos), torsions[pos][1], torsions[pos][2]));
				Real const new_rama_score( rama->eval_rama_score_residue( pose.aa(pos), save_torsions[1], save_torsions[2]));

				if (local_verbose)
					tccd.Debug << "rama_check: " << pos << ' ' << old_rama_score << ' ' << new_rama_score << std::endl;

				if ( new_rama_score > old_rama_score ) {
					Real const boltz_factor ( (old_rama_score-new_rama_score)/ccd_temperature );
					Real const negative_forty( -40.0f );
					Real const probability ( std::exp(std::max(negative_forty,boltz_factor) ) );
					if ( numeric::random::uniform() >= probability ) {
						// undo the change:
						torsions[pos] = save_torsions;

						// recover saved coords
						for ( int i=loop_begin; i<= loop_end; ++i ) {
							for ( Size j=1; j<= nbb; ++j ) {
								coords[i][j] = save_coords[i][j];
							}
						}

						get_overlap_pos( coords, torsions, cutpoint, direction, bond_angle1, bond_length, bond_angle2, M );

					} // rama reject
				} // rama-score got worse
			} // if ( rama_check )


			++nmoves;
	} //end while( nmoves < total_moves)
	// DONE!
	copy_torsions_to_pose( pose, loop_begin, loop_end, torsions );

}

//Pulling out the deviation calculations from fast_ccd_closure, returns a pair of Reals containing forward and backward deviations
std::pair<core::Real, core::Real>
get_deviation(
	core::pose::Pose const & pose,
	int const cutpoint
){
	Real const bond_angle1( pose.residue( cutpoint ).upper_connect().icoor().theta() );// CA-C=N bond angle
	Real const bond_angle2( pose.residue( cutpoint+1 ).lower_connect().icoor().theta() ); // C=N-CA bond angle
	Real const bond_length( pose.residue( cutpoint+1 ).lower_connect().icoor().d() ); // C=N distance
	
	int const n2c = { 1 }; // must be 1 and -1 (below) for proper incrementing in loops
	int const c2n = { -1 };
	
	Size const nbb( pose.residue( cutpoint ).mainchain_atoms().size() );
	
	vector1< vector1< Vector > > coords;
	vector1< vector1< Real > > torsions;
	load_coords_and_torsions( pose, coords, torsions );
	
	Matrix F,M;
	
	// forward_deviation:
	int direction = n2c;
	for ( int i = 1; i <= 3; ++i ) {
		F.col( i, coords[ cutpoint + 1 ][ i ] );
	}
	get_overlap_pos( coords, torsions, cutpoint, direction, bond_angle1, bond_length, bond_angle2, M );

	core::Real forward_deviation = 0.0;
	for ( int i = 1; i <= 3; ++i ) {
		for ( int j = 1; j <= 3; ++j ) {
			forward_deviation += numeric::square( M(j,i) - F(j,i) );
		}
	}
	forward_deviation = sqrt( forward_deviation / 3 );

	// backward_deviation:
	direction = c2n;
	for ( int i = 1; i <= 3; ++i ) {
		F.col( i, coords[ cutpoint     ][ nbb - 3 + i ] );
	}
	get_overlap_pos( coords, torsions, cutpoint, direction, bond_angle1, bond_length, bond_angle2, M );

	core::Real backward_deviation = 0.0;
	for ( int i = 1; i <= 3; ++i ) {
		for ( int j = 1; j <= 3; ++j ) {
			backward_deviation += numeric::square( M(j,i) - F(j,i) );
		}
	}
	backward_deviation = sqrt( backward_deviation / 3 );
	
	return std::make_pair(forward_deviation, backward_deviation);
}

} // namespace ccd
} // namespace loop_closure
} // namespace loops
} // namespace protocols
