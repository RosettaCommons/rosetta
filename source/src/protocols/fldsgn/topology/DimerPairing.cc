// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )


#include <protocols/fldsgn/topology/DimerPairing.hh>
#include <protocols/fldsgn/topology/SS_Info2.hh>
#include <ObjexxFCL/format.hh>

#include <utility/vector1.hh>


namespace protocols {
namespace fldsgn {
namespace topology {

/// @details Auto-generated virtual destructor
DimerPairing::~DimerPairing() = default;


/// @brief
DimerPairing::DimerPairing(
	Size const res1,
	Size const res2,
	Real const dist,
	Real const phi,
	Real const theta,
	Real const sigma,
	Real const dp,
	Size const sign1,
	Size const sign2,
	Real const score
): res1_( res1 ),
	res2_( res2 ),
	dist_( dist ),
	phi_( phi ),
	theta_( theta ),
	sigma_( sigma ),
	dp_( dp ),
	sign1_( sign1 ),
	sign2_( sign2 ),
	score_( score ),
	orient_( 'N' ),
	valid_( true )
{
	if ( is_parallel( phi, theta ) ) {
		orient_ = 'P';
	} else {
		orient_ = 'A';
	}
}


/// @brief the pairing is parallel ?
bool
DimerPairing::is_parallel( Real const phi, Real const theta )
{
	if ( phi <= 0.0 ) {
		if ( theta <= 70.0 ) {
			return true;
		} else {
			return false;
		}
	} else {
		if ( theta <= 110 ) {
			return true;
		} else {
			return false;
		}
	}
	// Should never get here
	runtime_assert( false );
} // is_parallel

std::ostream& operator<<( std::ostream & out, const DimerPairing & dp )
{
	using ObjexxFCL::format::I;
	using ObjexxFCL::format::F;
	using ObjexxFCL::format::A;
	out << I( 4, dp.res1() )     << ' ' << I( 4, dp.res2() )     << ' '
		<< F( 9, 3, dp.dist() )  << ' ' << F( 9, 3, dp.sigma() ) << ' '
		<< F( 9, 3, dp.phi() )   << ' ' << F( 9, 3, dp.theta() ) << ' ' << F( 9, 3, dp.dp() )       << ' '
		<< I( 2,    dp.sign1() ) << ' ' << I( 2, dp.sign2() )    << ' ' << A( 2,    dp.orient() ) << ' '
		<< F( 9, 3, dp.score() ) << std::endl;
	return out;
}

/// @brief
bool pointer_sorter( DimerPairingCOP const a, DimerPairingCOP const b )
{
	return ( a->score() < b->score() );
}

/// @brief
void
DimerPairings::finalize( SS_Info2 const & ss_info )
{
	std::sort( begin(), end(), pointer_sorter );

	for ( auto it=begin(), ite=end(); it!=ite ; ++it ) {

		DimerPairing const & pairing( **it );
		if ( !pairing.valid() ) continue;

		Size const res1( pairing.res1() );
		Size const res2( pairing.res2() );
		runtime_assert( res2 > res1 );

		Size const sign1( pairing.sign1() );
		Size const sign2( pairing.sign2() );

		Size const strand1( ss_info.strand_id( pairing.res1() ) );
		Size const strand2( ss_info.strand_id( pairing.res2() ) );
		runtime_assert( strand2 > strand1 );

		auto it2( it );
		++it2;
		for ( ; it2 != ite; ++it2 ) {

			DimerPairing & other( **it2 );
			if ( !other.valid() ) continue;

			Size const other_strand1( ss_info.strand_id( other.res1() ) );
			Size const other_strand2( ss_info.strand_id( other.res2() ) );

			//car if dimer1 and ss2 interact favorably and ss2 and dimer2 are in different strands
			//car and ss2 is on the same side of dimer1 as dimer2 is... then mark this pair as dissallowed
			if ( ( other.res1() == res1 && other_strand2 != strand2 && other.sign1() == sign1 ) ||
					( other.res2() == res1 && other_strand1 != strand2 && other.sign2() == sign1 ) ||
					( other.res1() == res2 && other_strand2 != strand1 && other.sign1() == sign2 ) ||
					( other.res2() == res2 && other_strand1 != strand1 && other.sign2() == sign2 ) ) {
				other.valid( false );
			}

		} // it2
	} // it

} // finalize


std::ostream& operator<<( std::ostream & out, const DimerPairings &dps )
{
	using ObjexxFCL::format::I;
	using ObjexxFCL::format::LJ;
	using ObjexxFCL::format::RJ;
	out  << LJ(4, "#" ) << ' '
		<< LJ(4, "res1" ) << ' ' << LJ(4, "res2"  ) << ' '
		<< RJ(9, "dist" ) << ' ' << RJ(9, "sigma" ) << ' '
		<< RJ(9,  "phi" ) << ' ' << RJ(9, "theta" ) << ' ' << RJ(9, "dp" ) << ' '
		<< LJ(2,   "s1" ) << ' ' << LJ(2,    "s2" ) << ' ' << LJ(2, "pr" ) << ' '
		<< RJ(9, "score") << std::endl;

	core::Size count( 0 );
	for ( auto const & dp : dps ) {
		if ( !( *dp).valid() ) continue;
		out << I( 4, count++ ) << ' ' ;
		out << (*dp);
	}
	return out;
}


} // ns topology
} // ns fldsgn
} // ns protocols
