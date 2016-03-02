// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/sequence/ABEGOManager.cc
/// @brief class for ABEGO
/// @author Nobuyasu Koga ( nobuyasu@uw.edu )

#include <core/chemical/ResidueType.hh>
#include <core/sequence/ABEGOManager.hh>
#include <core/pose/Pose.hh>
#include <utility/exit.hh>
#include <basic/Tracer.hh>
#include <cmath>
#include <utility/vector1.hh>

static THREAD_LOCAL basic::Tracer TR( "core.sequence.ABEGOManager" );

namespace core {
namespace sequence {

// @brief Auto-generated virtual destructor
ABEGOManager::~ABEGOManager() {}

/// @brief add line to specify abego region
void
ABEGO::add_line( Real const slope, Real const intercept, bool const region )
{
	lines_.push_back( Line( slope, intercept, region ) );
}

/// @brief  check input torsion angles is compatible with defined abego
bool
ABEGO::check_rama2( Real const & phi, Real const & psi )
{
	for ( Size ii=1; ii<=lines_.size(); ii++ ) {
		Real sign( psi - lines_[ ii ].slope * phi - lines_[ ii ].intercept );
		if ( (sign >= 0 && lines_[ ii ].region) || (sign < 0 && ! lines_[ ii ].region) ) {}
		else {
			return false;
		}
	}
	return true;
}

/// @brief check input torsion angles is compatible with defined abego
bool
ABEGO::check_rama( Real const & phi, Real const & psi, Real const & omega )
{
	if ( name() == 'X' ) {
		return true;
	}

	if ( cis_omega_ ) {

		return ( fabs( omega ) < 90.0 );

	} else {

		if ( fabs( omega ) < 90.0 ) return false;

		Real psii = psi;
		if ( (psi < -75.0 && phi < 0) || (psi < -100.0 && phi >= 0) ) {
			psii = psi + 360.0;
		}

		if ( phi >= phi_min_ && phi < phi_max_  &&
				psii >= psi_min_ && psii < psi_max_ ) {

			if ( lines_.size() < 1 ) {
				return true;
			} else {
				return check_rama2( phi, psii );
			}

		} else {
			return false;
		}
	}
}

/// @brief default constructor
ABEGOManager::ABEGOManager() : utility::pointer::ReferenceCount()
{
	initialize();
}

/// @brief copy constructor
ABEGOManager::ABEGOManager( ABEGOManager const & rval ) : utility::pointer::ReferenceCount(),
	totnum_abego_( rval.totnum_abego_ ),
	name2abego_( rval.name2abego_ )
{}

/// @brief intialize ABEGOManager
void
ABEGOManager::initialize()
{
	ABEGO A( 'A', -180.0,   0.0,  -75.0,  50.0, false );
	ABEGO B( 'B', -180.0,   0.0,   50.0, 285.5, false );
	ABEGO E( 'E',    0.0, 180.5,  100.0, 260.5, false );
	ABEGO G( 'G',    0.0, 180.5, -100.0, 100.0, false );
	ABEGO O( 'O',    0.0,   0.0,   0.0,    0.0,  true );

	ABEGO S( 'S', -180.0,   0.0,  100.0, 195.0, false );
	S.add_line( -1.6, 4, false );

	ABEGO P( 'P', -180.0,   0.0,  100.0, 195.0, false );
	P.add_line( -1.6, 4, true  );

	ABEGO D( 'D', -180.0,    0.0, 195.0, 285.5, false );
	ABEGO Z( 'Z', -180.0, -100.0,  50.0, 100.0, false );
	ABEGO Y( 'Y', -100.0,    0.0,  50.0, 100.0, false );

	ABEGO M( 'M', -180.0,  -90.0,  -75.0,  50.0, false );
	ABEGO N( 'N',  -90.0,    0.0,  -75.0,  50.0, false );

	// X represents all torsion space ( the following parameters are no meaning. )
	ABEGO X( 'X', 0.0, 0.0, 0.0, 0.0, false );

	name2abego_[ 1 ] = A;
	name2abego_[ 2 ] = B;
	name2abego_[ 3 ] = E;
	name2abego_[ 4 ] = G;
	name2abego_[ 5 ] = O;
	name2abego_[ 6 ] = S;
	name2abego_[ 7 ] = P;
	name2abego_[ 8 ] = D;
	name2abego_[ 9 ] = Z;
	name2abego_[ 10 ] = Y;
	name2abego_[ 11 ] = M;
	name2abego_[ 12 ] = N;
	name2abego_[ 13 ] = X;

	totnum_abego_ = 13;

}

/// @brief check input torsion angle are in a given abego region
bool
ABEGOManager::check_rama( char const & symbol, Real const & phi, Real const & psi, Real const & omega )
{
	Size idx = symbol2index( symbol );
	return name2abego_[ idx ].check_rama( phi, psi, omega );
}


/// @brief get abego
Size
ABEGOManager::torsion2index( Real const phi, Real const psi, Real const omega, Size const level )
{
	switch( level ) {
	case 1 :
		return torsion2index_level1( phi, psi, omega );
	case 2 :
		return torsion2index_level2( phi, psi, omega );
	case 3 :
		return torsion2index_level3( phi, psi, omega );
	case 4 :
		return torsion2index_level4( phi, psi, omega );
	default :
		TR << " [ERROR] Unrecognized level  " << level << std::endl;
		runtime_assert( false );
		return 0;
	}
}


/// @brief get abegeo index from torsion angles: ABEGO
Size
ABEGOManager::torsion2index_level1( Real const phi, Real const psi, Real const omega )
{
	if ( fabs( omega ) < 90.0 ) {
		return 5;   // cis-omega
	} else if ( phi >= 0.0 ) {
		if ( -100.0 <= psi && psi < 100.0 ) {
			return 4; // alpha-L
		} else {
			return 3; // E
		}
	} else {
		if ( -75.0 <= psi && psi < 50.0 ) {
			return 1; // helical
		} else {
			return 2; // beta
		}
	}
	return 0;
}


/// @brief get abego index from torsion angles: ABEGOD
Size
ABEGOManager::torsion2index_level2( Real const phi, Real const psi, Real const omega )
{
	if ( fabs( omega ) < 90.0 ) {
		return 5; // cis-omega
	} else if ( phi >= 0.0 ) {
		if ( -100.0 <= psi && psi < 100.0 ) {
			return 4; // alpha-L
		} else {
			return 3; // E
		}
	} else {

		if ( -75.0 <= psi && psi < 50.0 ) {
			return 1; //helical
		} else {

			Real ppsi;
			if ( psi < -75.0 ) {
				ppsi = psi + 360.0;
			} else {
				ppsi = psi;
			}
			if ( ppsi >= 195.0 ) {
				return 8; // D
			} else {
				return 2; // B
			}

		}

	}
	return 0;
}

/// @brief get abego index from torsion angles: ASPZYD
Size
ABEGOManager::torsion2index_level3( Real const phi, Real const psi, Real const omega )
{
	if ( fabs( omega ) < 90.0 ) {
		return 5; // cis-omega
	} else if ( phi >= 0.0 ) {
		if ( -100.0 <= psi && psi < 100.0 ) {
			return 4; // alpha-L
		} else {
			return 3; // E
		}
	} else {
		if ( -75.0 <= psi && psi < 50.0 ) {
			return 1; //helical
		} else if ( 50.0 <= psi && psi < 100.0 ) {
			if ( phi >= -100 ) {
				return 10; // Y
			} else {
				return 9; // Z
			}
		} else {
			Real ppsi;
			if ( psi < -75.0 ) {
				ppsi = psi + 360.0;
			} else {
				ppsi = psi;
			}
			if ( ppsi >= 195.0 ) {
				return 8; // D
			} else {
				Real sign( ppsi - ( -1.6*phi + 4.0 ) );
				if ( sign >= 0 ) {
					return 7; // P
				} else {
					return 6; // S
				}
			}
		}
	}
	return 0;
}

/// @brief get abego index from torsion angles: MNSPZYD
Size
ABEGOManager::torsion2index_level4( Real const phi, Real const psi, Real const omega )
{
	if ( fabs( omega ) < 90.0 ) {
		return 5; // cis-omega
	} else if ( phi >= 0.0 ) {
		if ( -100.0 <= psi && psi < 100.0 ) {
			return 4; // alpha-L
		} else {
			return 3; // E
		}
	} else {
		if ( -75.0 <= psi && psi < 50.0 ) {
			if ( phi >= -90.0 ) {
				return 12; // N
			} else {
				return 11; // M
			}
		} else if ( 50.0 <= psi && psi < 100.0 ) {
			if ( phi >= -100 ) {
				return 10; // Y
			} else {
				return 9; // Z
			}
		} else {
			Real ppsi;
			if ( psi < -75.0 ) {
				ppsi = psi + 360.0;
			} else {
				ppsi = psi;
			}

			if ( ppsi >= 195.0 ) {
				return 8; // D
			} else {
				Real sign( ppsi - ( -1.6*phi + 4.0 ) );
				if ( sign >= 0 ) {
					return 7; // P
				} else {
					return 6; // S
				}
			}
		}
	}
	return 0;
}

/// @brief transform abego index to symbol
char
ABEGOManager::index2symbol( Size const & idx )
{
	switch( idx ) {
	case 1 :
		return 'A';
	case 2 :
		return 'B';
	case 3 :
		return 'E';
	case 4 :
		return 'G';
	case 5 :
		return 'O';
	case 6 :
		return 'S';
	case 7 :
		return 'P';
	case 8 :
		return 'D';
	case 9 :
		return 'Z';
	case 10 :
		return 'Y';
	case 11 :
		return 'M';
	case 12 :
		return 'N';
	case 13 :
		return 'X';
	default :
		TR << " [ERROR] Unrecognized abego index: " << idx << std::endl;
		runtime_assert( false );
		return 0;
	}
}


/// @brief transform abego symbol to index
Size
ABEGOManager::symbol2index( char const & symbol )
{
	switch( symbol ) {
	case 'A' :
	case 'a' :
		return 1;
	case 'B' :
	case 'b' :
		return 2;
	case 'E' :
	case 'e' :
		return 3;
	case 'G' :
	case 'g' :
		return 4;
	case 'O' :
	case 'o' :
		return 5;
	case 'S' :
	case 's' :
		return 6;
	case 'P' :
	case 'p' :
		return 7;
	case 'D' :
	case 'd' :
		return 8;
	case 'Z' :
	case 'z' :
		return 9;
	case 'Y' :
	case 'y' :
		return 10;
	case 'M' :
	case 'm' :
		return 11;
	case 'N' :
	case 'n' :
		return 12;
	case 'X' :
	case 'x' :
		return 13;
	default :
		TR << " [ERROR] Unrecognized abego index: " << symbol << std::endl;
		runtime_assert( false );
		return 0;
	}
}

/// @brief transform abego symbol string to base5 index. This is used to quickly pool the abego from Alex's hd5 database
Size ABEGOManager::symbolString2base5index( std::string symbolString){
	std::string allowedChar = "ABEGO";
	for ( Size ii=0; ii<symbolString.size(); ++ii ) { //check for only allowed characters
		if ( allowedChar.find(symbolString.substr(ii,1)) == std::string::npos ) {
			utility_exit_with_message("Looking for " + symbolString.substr(ii,1) + " which doesn't exist in the hdf5 database");
		}
	}
	Size base5index = 0;
	for ( Size ii=0; ii<symbolString.size(); ++ii ) {
		Size symbolValue = symbol2index( symbolString[ii]);
		Size abego_index_radix = pow(5,ii);
		base5index += abego_index_radix*(symbolValue-1);
	}
	return(base5index);
}

/// @brief transform abego symbol string to base5 index. This is used to quickly pool the abego from Alex's hd5 database
std::string ABEGOManager::base5index2symbolString( Size base5index,Size length){
	Size tmp_base5index=base5index;
	std::string symbolString ="";
	for ( Size ii=0; ii<length; ++ii ) {
		Size index = tmp_base5index % 5;
		tmp_base5index = tmp_base5index/5;
		symbolString+=index2symbol(index+1);
	}
	return(symbolString);
}

/// @brief get abego sequence from pose
utility::vector1< std::string >
ABEGOManager::get_symbols( Pose const & pose, Size const begin, Size const end, Size const level )
{
	runtime_assert( begin >= 1 && end <= pose.total_residue() );

	utility::vector1< String > symbols;
	for ( Size i=begin; i<=end; i++ ) {
		if ( pose.residue_type( i ).is_protein() ) {
			Size idx = torsion2index( pose.phi( i ), pose.psi( i ), pose.omega( i ), level );
			std::ostringstream symbol;
			symbol << index2symbol( idx );
			symbols.push_back( symbol.str() );
		} else symbols.push_back("-");
	}
	return symbols;
}

/// @brief get abego sequence from pose
utility::vector1< std::string >
ABEGOManager::get_symbols( Pose const & pose, Size const level )
{
	return get_symbols( pose, 1, pose.total_residue(), level );
}


/// @brief get abego string
std::string
ABEGOManager::get_abego_string( utility::vector1< std::string > abego )
{
	std::ostringstream output;
	for ( Size ii=1; ii<=abego.size(); ++ii ) {
		Size length = abego[ ii ].length();
		if ( length > 1 ) {
			std::ostringstream multi;
			multi << "[";
			for ( Size jj=0; jj<abego[ ii ].length(); ++jj ) {
				multi << abego[ ii ].at( jj );
			}
			multi << "]";
			output << multi.str();
		} else {
			output << abego[ ii ].at( 0 );
		}
	}
	return output.str();
}


/// @brief utility for getting abego
utility::vector1< std::string >
get_abego( core::pose::Pose const & pose, core::Size const level )
{
	ABEGOManager am;
	return am.get_symbols( pose, level );
}

} // namespace sequence
} // namespace core

