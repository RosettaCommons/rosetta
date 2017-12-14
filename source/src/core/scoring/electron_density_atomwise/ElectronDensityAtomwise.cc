
// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/scoring/electron_density_atomwise/ElectronDensityAtomwise.cc
/// @brief  elec_dens_atomwise scoring method implementation
/// @author Fang-Chieh Chou

// Unit headers
#include <core/scoring/electron_density_atomwise/ElectronDensityAtomwise.hh>

// Package headers
#include <core/scoring/EnergyMap.hh>
#include <core/scoring/methods/EnergyMethodOptions.hh>

// Project headers
#include <core/pose/Pose.hh>
#include <core/scoring/Energies.hh>
#include <core/conformation/Residue.hh>
#include <core/kinematics/Edge.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <basic/options/option.hh>
#include <basic/Tracer.hh>

#include <numeric/constants.hh>
#include <numeric/fourier/FFT.hh>
#include <numeric/xyzMatrix.hh>
#include <numeric/xyzVector.hh>
#include <numeric/xyz.functions.hh>
#include <numeric/xyzVector.io.hh>
#include <numeric/statistics/functions.hh>

#include <core/scoring/electron_density/SplineInterp.hh>


#include <basic/options/keys/edensity.OptionKeys.gen.hh>

// Utility headers
#include <utility/string_util.hh>

// C++ headers
#include <fstream>
#include <limits>


//Auto Headers
#include <core/chemical/AtomType.hh>
#include <core/chemical/AtomTypeSet.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/RT.hh>


#ifdef WIN32
#define _USE_MATH_DEFINES
#include <math.h>
#endif


namespace core {
namespace scoring {
namespace electron_density_atomwise {

using namespace core;
using namespace basic::options;

static basic::Tracer TR( "core.scoring.electron_density_atomwise.ElectronDensityAtomwise" );


/// null constructor
ElectronDensityAtomwise::ElectronDensityAtomwise() {
	is_map_loaded_ = false;
	is_score_precomputed_ = false;
	grid = numeric::xyzVector< int >( 0, 0, 0 );
	orig = numeric::xyzVector< int >( 0, 0, 0 );
	cell_dimensions = numeric::xyzVector< float >( 1, 1, 1 );
	cell_angles = numeric::xyzVector< float >( 90, 90, 90 );
	// resolution ... set in readMRCandResize
	map_reso = -1.0;
}

////////
//Functions used for map loading

const int CCP4HDSIZE = 1024;  // size of CCP4/MRC header

inline float d2r( float d ) {
	return ( d * M_PI / 180.0 );
}
inline double d2r( double d ) {
	return ( d * M_PI / 180.0 );
}
inline float square( float x ) {
	return ( x * x );
}
inline double square( double x ) {
	return ( x * x );
}

// x mod y, returns z in [0,y-1]
inline int pos_mod( int x, int y ) {
	int r = x % y;

	if ( r < 0 ) r += y;

	return r;
}
inline float pos_mod( float x, float y ) {
	float r = std::fmod( x, y );

	if ( r < 0 ) r += y;

	return r;
}
inline double pos_mod( double x, double y ) {
	double r = std::fmod( x, y );

	if ( r < 0 ) r += y;

	return r;
}

// Endianness swap
// Only works with aligned 4-byte quantities
static void swap4_aligned( void *v, long ndata ) {
	auto *data = ( int * ) v;
	long i;
	int *N;
	// AMW: again cppcheck says scope of N can be reduced,
	// but I want other eyes on this.
	for ( i = 0; i < ndata; i++ ) {
		N = data + i;
		*N = ( ( ( *N >> 24 ) & 0xff ) | ( ( *N & 0xff ) << 24 ) | ( ( *N >> 8 ) & 0xff00 ) | ( ( *N & 0xff00 ) << 8 ) );
	}
}

//////
void ElectronDensityAtomwise::computeCrystParams() {
	// recompute reciprocal cell
	// f2c, c2f
	core::Real ca = cos( d2r( cell_angles[0] ) ), cb = cos( d2r( cell_angles[1] ) ), cg = cos( d2r( cell_angles[2] ) );
	core::Real sa = sin( d2r( cell_angles[0] ) ), sb = sin( d2r( cell_angles[1] ) ), sg = sin( d2r( cell_angles[2] ) );
	// conversion from fractional cell coords to cartesian coords
	f2c = numeric::xyzMatrix<core::Real>::rows(
		cell_dimensions[0], cell_dimensions[1] * cg, cell_dimensions[2] * cb,
		0.0, cell_dimensions[1] * sg, cell_dimensions[2] * ( ca - cb * cg ) / sg,
		0.0, 0.0, cell_dimensions[2] * sb * sqrt( 1.0 - square( ( cb * cg - ca ) / ( sb * sg ) ) ) );
	core::Real D = f2c.det();

	if ( D == 0 ) {
		TR.Warning << "Invalid crystal cell dimensions." << std::endl;
		return;
	}

	// c2f is inverse of f2c
	c2f = numeric::xyzMatrix<core::Real>::rows(
		( f2c( 2, 2 ) * f2c( 3, 3 ) - f2c( 2, 3 ) * f2c( 3, 2 ) ) / D,
		- ( f2c( 1, 2 ) * f2c( 3, 3 ) - f2c( 1, 3 ) * f2c( 3, 2 ) ) / D,
		( f2c( 1, 2 ) * f2c( 2, 3 ) - f2c( 1, 3 ) * f2c( 2, 2 ) ) / D,
		- ( f2c( 2, 1 ) * f2c( 3, 3 ) - f2c( 3, 1 ) * f2c( 2, 3 ) ) / D,
		( f2c( 1, 1 ) * f2c( 3, 3 ) - f2c( 1, 3 ) * f2c( 3, 1 ) ) / D,
		- ( f2c( 1, 1 ) * f2c( 2, 3 ) - f2c( 1, 3 ) * f2c( 2, 1 ) ) / D,
		( f2c( 2, 1 ) * f2c( 3, 2 ) - f2c( 3, 1 ) * f2c( 2, 2 ) ) / D,
		- ( f2c( 1, 1 ) * f2c( 3, 2 ) - f2c( 1, 2 ) * f2c( 3, 1 ) ) / D,
		( f2c( 1, 1 ) * f2c( 2, 2 ) - f2c( 1, 2 ) * f2c( 2, 1 ) ) / D );
	cell_volume = cell_dimensions[0] * cell_dimensions[1] * cell_dimensions[2] * sqrt( 1 - square( ca ) - square( cb ) - square( cg ) + 2 * ca * cb * cg );
	// reciprocal space cell dimensions
	r_cell_dimensions[0] = cell_dimensions[1] * cell_dimensions[2] * sa / cell_volume;
	r_cell_dimensions[1] = cell_dimensions[0] * cell_dimensions[2] * sb / cell_volume;
	r_cell_dimensions[2] = cell_dimensions[0] * cell_dimensions[1] * sg / cell_volume;
	cos_r_cell_angles[0] = cos( asin( std::min( std::max( cell_volume / ( cell_dimensions[0] * cell_dimensions[1] * cell_dimensions[2] * sb * sg ) , -1.0 ) , 1.0 ) ) );
	cos_r_cell_angles[1] = cos( asin( std::min( std::max( cell_volume / ( cell_dimensions[0] * cell_dimensions[1] * cell_dimensions[2] * sa * sg ) , -1.0 ) , 1.0 ) ) );
	cos_r_cell_angles[2] = cos( asin( std::min( std::max( cell_volume / ( cell_dimensions[0] * cell_dimensions[1] * cell_dimensions[2] * sa * sb ) , -1.0 ) , 1.0 ) ) );
	r_cell_volume = 1.0 / cell_volume;
}

/////////////////////////////////////
// parse symmops from ccp4 map header
void ElectronDensityAtomwise::initializeSymmOps ( utility::vector1< std::string > const & symList ) {
	using core::kinematics::RT;
	symmOps.clear();

	if ( symList.empty() ) { // no symminfo in header, assume P 1
		symmOps.push_back( RT( numeric::xyzMatrix< core::Real >::identity(),
			numeric::xyzVector< core::Real >( 0.0, 0.0, 0.0 ) ) );
	}

	for ( std::string const & line : symList ) {
		utility::vector1< std::string > rows = utility::string_split( line, ',' );

		if ( rows.size() != 3 ) {
			TR.Error << "invalid symmop in map file" << std::endl;
			TR.Error << line << std::endl;
			TR.Error << "Setting symmetry to P1 and continuing!" << line << std::endl;
			// should we throw an exception here????  nah, just set symm to P1 and continue
			symmOps.clear();
			symmOps.push_back( RT( numeric::xyzMatrix< core::Real >::identity(),
				numeric::xyzVector< core::Real >( 0.0, 0.0, 0.0 ) ) );
			return;
		}

		// _REALLY_ simple parser
		numeric::xyzMatrix< core::Real > rot( 0 );
		numeric::xyzVector< core::Real > trans( 0, 0, 0 );
		//int k;

		for ( int j = 0; j < 3; ++j ) {
			int k = rows[j+1].find( '/' );

			if ( k != ( int ) std::string::npos ) {
				// numerator
				int startNum = rows[j+1].find_last_not_of( "0123456789", k - 1 ) + 1;
				int startDenom = k + 1;
				float denom = std::atof( &rows[j+1][startDenom] );
				// make sure this shift corresponds to a point in the map
				trans[j] = std::atof( &rows[j+1][startNum] ) / denom;
			} else {
				trans[j] = 0;
			}

			if ( rows[j+1].find( "-X" ) != std::string::npos ) {
				rot( j + 1, 1 ) = -1;
			} else if ( rows[j+1].find( "X" ) != std::string::npos ) {
				rot( j + 1, 1 ) = 1;
			}

			if ( rows[j+1].find( "-Y" ) != std::string::npos ) {
				rot( j + 1, 2 ) = -1;
			} else if ( rows[j+1].find( "Y" ) != std::string::npos ) {
				rot( j + 1, 2 ) = 1;
			}

			if ( rows[j+1].find( "-Z" ) != std::string::npos ) {
				rot( j + 1, 3 ) = -1;
			} else if ( rows[j+1].find( "Z" ) != std::string::npos ) {
				rot( j + 1, 3 ) = 1;
			}
		}

		symmOps.push_back( RT( rot, trans ) );
	}
}


/////////////////////////////////////
// expand density to cover unit cell
// maintain origin
void ElectronDensityAtomwise::expandToUnitCell() {
	numeric::xyzVector< int > extent( density.u1(), density.u2(), density.u3() );

	// if it already covers unit cell do nothing
	if ( grid[0] == extent[0] && grid[1] == extent[1] && grid[2] == extent[2] ) {
		return;
	}

	ObjexxFCL::FArray3D< float > newDensity( grid[0], grid[1], grid[2], 0.0 );
	// copy the block
	int limX = std::min( extent[0], grid[0] ),
		limY = std::min( extent[1], grid[1] ),
		limZ = std::min( extent[2], grid[2] );

	for ( int x = 1; x <= limX; ++x ) {
		for ( int y = 1; y <= limY; ++y ) {
			for ( int z = 1; z <= limZ; ++z ) {
				newDensity( x, y, z ) = density( x, y, z );
			}
		}
	}

	// apply symmetry
	// why backwards? it is a mystery
	for ( int x = grid[0]; x >= 1; --x ) {
		for ( int y = grid[1]; y >= 1; --y ) {
			for ( int z = grid[2]; z >= 1; --z ) {
				if ( x <= limX && y <= limY && z <= limZ )  continue;

				numeric::xyzVector<core::Real> fracX(
					( ( core::Real ) x + orig[0] - 1 ) / grid[0],
					( ( core::Real ) y + orig[1] - 1 ) / grid[1],
					( ( core::Real ) z + orig[2] - 1 ) / grid[2] );

				for ( auto const & symmOp : symmOps ) {
					numeric::xyzVector<core::Real> SfracX =
						symmOp.get_rotation() * fracX + symmOp.get_translation();
					// indices of symm copy
					int Sx = pos_mod( ( int ) floor( SfracX[0] * grid[0] + 0.5 - orig[0] ) , grid[0] ) + 1;
					int Sy = pos_mod( ( int ) floor( SfracX[1] * grid[1] + 0.5 - orig[1] ) , grid[1] ) + 1 ;
					int Sz = pos_mod( ( int ) floor( SfracX[2] * grid[2] + 0.5 - orig[2] ) , grid[2] ) + 1 ;

					if ( Sx <= limX && Sy <= limY && Sz <= limZ ) {
						newDensity( x, y, z ) = density( Sx, Sy, Sz );
					}
				}
			}
		}
	}

	// new map!
	density = newDensity;
}

/////////////////////////////////////
// resize a map (using FFT-interpolation)
void ElectronDensityAtomwise::resize( core::Real approxGridSpacing ) {
	// potentially expand map to cover entire unit cell
	if ( grid[0] != density.u1() || grid[1] != density.u2() || grid[2] != density.u3() ) {
		TR.Error << "resize() not supported for maps not covering the entire unit cell." << std::endl;
		TR.Error << "   " << grid[0] << " != " << density.u1()
			<< " || " << grid[1] << " != " << density.u2()
			<< " || " << grid[2] << " != " << density.u3() << std::endl;
		exit( 1 );
	}

	// compute new dimensions & origin
	numeric::xyzVector<int> newDims,  newGrid;
	numeric::xyzVector<double> newOri;
	newDims[0] = ( int ) floor( cell_dimensions[0] / approxGridSpacing + 0.5 );
	newDims[1] = ( int ) floor( cell_dimensions[1] / approxGridSpacing + 0.5 );
	newDims[2] = ( int ) floor( cell_dimensions[2] / approxGridSpacing + 0.5 );
	newOri[0] = newDims[0] * orig[0] / ( ( core::Real ) grid[0] );
	newOri[1] = newDims[1] * orig[1] / ( ( core::Real ) grid[1] );
	newOri[2] = newDims[2] * orig[2] / ( ( core::Real ) grid[2] );
	newGrid = newDims;
	ObjexxFCL::FArray3D< std::complex<double> > newDensity;
	newDensity.dimension( newDims[0], newDims[1], newDims[2] );
	TR << "Resizing " << density.u1() << "x" << density.u2() << "x" << density.u3() << " to "
		<< newDensity.u1() << "x" << newDensity.u2() << "x" << newDensity.u3() << std::endl;
	// convert map to complex<double>
	ObjexxFCL::FArray3D< std::complex<double> > Foldmap, Fnewmap;
	Fnewmap.dimension( newDims[0], newDims[1], newDims[2] );
	// fft
	ObjexxFCL::FArray3D< std::complex<double> > cplx_density;
	cplx_density.dimension( density.u1() , density.u2() , density.u3() );

	for ( int i = 0; i < density.u1() *density.u2() *density.u3(); ++i ) cplx_density[i] = ( double ) density[i];

	numeric::fourier::fft3( cplx_density, Foldmap );

	// reshape (handles both shrinking and growing in each dimension)
	for ( int i = 0; i < Fnewmap.u1() *Fnewmap.u2() *Fnewmap.u3(); ++i ) Fnewmap[i] = std::complex<double> ( 0, 0 );

	numeric::xyzVector<int> nyq( std::min( Foldmap.u1(), Fnewmap.u1() ) / 2,
		std::min( Foldmap.u2(), Fnewmap.u2() ) / 2,
		std::min( Foldmap.u3(), Fnewmap.u3() ) / 2 );
	numeric::xyzVector<int> nyqplus1_old( std::max( Foldmap.u1() - ( std::min( Foldmap.u1(), Fnewmap.u1() ) - nyq[0] ) + 1 , nyq[0] + 1 ) ,
		std::max( Foldmap.u2() - ( std::min( Foldmap.u2(), Fnewmap.u2() ) - nyq[1] ) + 1 , nyq[1] + 1 ) ,
		std::max( Foldmap.u3() - ( std::min( Foldmap.u3(), Fnewmap.u3() ) - nyq[2] ) + 1 , nyq[2] + 1 ) );
	numeric::xyzVector<int> nyqplus1_new( std::max( Fnewmap.u1() - ( std::min( Foldmap.u1(), Fnewmap.u1() ) - nyq[0] ) + 1 , nyq[0] + 1 ) ,
		std::max( Fnewmap.u2() - ( std::min( Foldmap.u2(), Fnewmap.u2() ) - nyq[1] ) + 1 , nyq[1] + 1 ) ,
		std::max( Fnewmap.u3() - ( std::min( Foldmap.u3(), Fnewmap.u3() ) - nyq[2] ) + 1 , nyq[2] + 1 ) );

	for ( int i = 1; i <= Fnewmap.u1(); i++ ) {
		for ( int j = 1; j <= Fnewmap.u2(); j++ ) {
			for ( int k = 1; k <= Fnewmap.u3(); k++ ) {
				if ( i - 1 <= nyq[0] ) {
					if ( j - 1 <= nyq[1] ) {
						if ( k - 1 <= nyq[2] ) {
							Fnewmap( i, j, k ) = Foldmap( i, j, k );
						} else if ( k - 1 >= nyqplus1_new[2] ) {
							Fnewmap( i, j, k ) = Foldmap( i, j, k - nyqplus1_new[2] + nyqplus1_old[2] );
						}
					} else if ( j - 1 >= nyqplus1_new[1] ) {
						if ( k - 1 <= nyq[2] ) {
							Fnewmap( i, j, k ) = Foldmap( i, j - nyqplus1_new[1] + nyqplus1_old[1], k );
						} else if ( k - 1 >= nyqplus1_new[2] ) {
							Fnewmap( i, j, k ) = Foldmap( i, j - nyqplus1_new[1] + nyqplus1_old[1], k - nyqplus1_new[2] + nyqplus1_old[2] );
						}
					}
				} else if ( i - 1 >= nyqplus1_new[0] ) {
					if ( j - 1 <= nyq[1] ) {
						if ( k - 1 <= nyq[2] ) {
							Fnewmap( i, j, k ) = Foldmap( i - nyqplus1_new[0] + nyqplus1_old[0], j, k );
						} else if ( k - 1 >= nyqplus1_new[2] ) {
							Fnewmap( i, j, k ) = Foldmap( i - nyqplus1_new[0] + nyqplus1_old[0], j, k - nyqplus1_new[2] + nyqplus1_old[2] );
						}
					} else if ( j - 1 >= nyqplus1_new[1] ) {
						if ( k - 1 <= nyq[2] ) {
							Fnewmap( i, j, k ) = Foldmap( i - nyqplus1_new[0] + nyqplus1_old[0], j - nyqplus1_new[1] + nyqplus1_old[1], k );
						} else if ( k - 1 >= nyqplus1_new[2] ) {
							Fnewmap( i, j, k ) = Foldmap( i - nyqplus1_new[0] + nyqplus1_old[0],
								j - nyqplus1_new[1] + nyqplus1_old[1],
								k - nyqplus1_new[2] + nyqplus1_old[2] );
						}
					}
				}
			}
		}
	}

	// ifft
	numeric::fourier::ifft3( Fnewmap, newDensity );
	// update density
	density.dimension( newDims[0], newDims[1], newDims[2] );

	for ( int i = 0; i < newDims[0]*newDims[1]*newDims[2] ; ++i ) {
		density[i] = ( float ) newDensity[i].real();
	}

	grid = newGrid;
	orig = newOri;
	TR << " new extent: " << density.u1() << " x " << density.u2() << " x " << density.u3() << std::endl;
	TR << " new origin: " << orig[0] << " x " << orig[1] << " x " << orig[2] << std::endl;
	TR << "   new grid: " << grid[0] << " x " << grid[1] << " x " << grid[2] << std::endl;
}

///////////////////////////////////////////////////////////////////////////
//Load in the CCP4 map
void
ElectronDensityAtomwise::readMRCandResize() {
	TR << "Loading Density Map" << std::endl;

	if ( !basic::options::option[ basic::options::OptionKeys::edensity::mapfile ].user() ) {
		utility_exit_with_message( "No density map specified for electron density scoring." );
	}

	std::string map_file = basic::options::option[ basic::options::OptionKeys::edensity::mapfile ]();
	map_reso = basic::options::option[ basic::options::OptionKeys::edensity::mapreso ]();
	grid_spacing = basic::options::option[ basic::options::OptionKeys::edensity::grid_spacing ]();
	std::fstream mapin( map_file.c_str() , std::ios::binary | std::ios::in );
	char mapString[4], symData[81];
	int  crs2xyz[3], extent[3], mode, symBytes, origin_xyz[3];
	int  xyz2crs[3], vol_xsize, vol_ysize, vol_zsize;
	int xIndex, yIndex, zIndex, vol_xySize, coord[3];
	long dataOffset, filesize;
	float *rowdata;
	bool swap = false;

	if ( !mapin ) {
		TR.Error << "Error opening MRC map " << map_file << ".  Not loading map." << std::endl;
		utility_exit_with_message( "Fail to load the density map." );
	}

	if ( !mapin.read ( reinterpret_cast <char*> ( extent ), 3 * sizeof( int ) )
			|| !mapin.read ( reinterpret_cast <char*> ( &mode ), 1 * sizeof( int ) )
			|| !mapin.read ( reinterpret_cast <char*> ( &origin_xyz[0] ), 3 * sizeof( int ) )
			|| !mapin.read ( reinterpret_cast <char*> ( &grid[0] ), 3 * sizeof( int ) )
			|| !mapin.read ( reinterpret_cast <char*> ( &cell_dimensions[0] ), 3 * sizeof( float ) )
			|| !mapin.read ( reinterpret_cast <char*> ( &cell_angles[0] ), 3 * sizeof( float ) )
			|| !mapin.read ( reinterpret_cast <char*> ( crs2xyz ), 3 * sizeof( int ) ) )  {
		TR.Error << "Improperly formatted line in MRC map.  Not loading map." << std::endl;
		utility_exit_with_message( "Fail to load the density map." );
	}

	// Check the number of bytes used for storing symmetry operators
	mapin.seekg( 92, std::ios::beg );

	if ( !mapin.read( reinterpret_cast <char*>( &symBytes ), 1 * sizeof( int ) ) ) {
		TR.Error << "Failed reading symmetry bytes record.  Not loading map." << "\n";
		utility_exit_with_message( "Fail to load the density map." );
	}

	// alt: MRC files have floating-point origin info at byte 196
	// read this and try to figure out if it is used
	float altorigin[3];
	mapin.seekg( 196, std::ios::beg );

	if ( !mapin.read( reinterpret_cast <char*>( altorigin ), 3 * sizeof( float ) ) ) {
		TR.Error << "Improperly formatted line in MRC map.  Not loading map." << std::endl;
		utility_exit_with_message( "Fail to load the density map." );
	}

	// Check for the string "MAP" at byte 208, indicating a CCP4 file.
	mapin.seekg( 208, std::ios::beg );
	mapString[3] = '\0';

	if ( !mapin.read( mapString, 3 ) || ( std::string( mapString ) != "MAP" ) ) {
		TR.Error << "'MAP' string missing, not a valid MRC map.  Not loading map." << std::endl;
		utility_exit_with_message( "Fail to load the density map." );
	}

	// Check the file endianness
	if ( mode != 2 ) {
		swap4_aligned( &mode, 1 );

		if ( mode != 2 ) {
			TR.Error << "Non-real (32-bit float) data types are unsupported.  Not loading map." << std::endl;
			utility_exit_with_message( "Fail to load the density map." );
		} else {
			swap = true; // enable byte swapping
		}
	}

	// Swap all the information obtained from the header
	if ( swap ) {
		swap4_aligned( extent, 3 );
		swap4_aligned( &origin_xyz[0], 3 );
		swap4_aligned( &altorigin[0], 3 );
		swap4_aligned( &grid[0], 3 );
		swap4_aligned( &cell_dimensions[0], 3 );
		swap4_aligned( &cell_angles[0], 3 );
		swap4_aligned( crs2xyz, 3 );
		swap4_aligned( &symBytes, 1 );
	}

	TR << " Setting resolution to " << map_reso << "A" << std::endl;
	TR << " Read density map'" << map_file << "'" << std::endl;
	TR << "     extent: " << extent[0] << " x " << extent[1] << " x " << extent[2] << std::endl;
	TR << "     origin: " << origin_xyz[0] << " x " << origin_xyz[1] << " x" << origin_xyz[2] << std::endl;
	TR << "  altorigin: " << altorigin[0] << " x " << altorigin[1] << " x " << altorigin[2] << std::endl;
	TR << "       grid: " << grid[0] << " x " << grid[1] << " x " << grid[2] << std::endl;
	TR << "    celldim: " << cell_dimensions[0] << " x " << cell_dimensions[1] << " x " << cell_dimensions[2] << std::endl;
	TR << " cellangles: " << cell_angles[0] << " x " << cell_angles[1] << " x " << cell_angles[2] << std::endl;
	TR << "    crs2xyz: " << crs2xyz[0] << " x " << crs2xyz[1] << " x " << crs2xyz[2] << std::endl;
	TR << "   symBytes: " << symBytes << std::endl;;
	// Check the dataOffset: this fixes the problem caused by files claiming
	// to have symmetry records when they do not.
	mapin.seekg( 0, std::ios::end );
	filesize = mapin.tellg();
	dataOffset = filesize - 4 * ( extent[0] * extent[1] * extent[2] );

	if ( dataOffset != ( CCP4HDSIZE + symBytes ) ) {
		if ( dataOffset == CCP4HDSIZE ) {
			// Bogus symmetry record information
			TR.Warning << "File contains bogus symmetry record.  Continuing." << std::endl;
			symBytes = 0;
		} else if ( dataOffset < CCP4HDSIZE ) {
			TR.Error << "File appears truncated and doesn't match header.  Not loading map." << std::endl;
			utility_exit_with_message( "Fail to load the density map." );
		} else if ( ( dataOffset > CCP4HDSIZE ) && ( dataOffset < ( 1024 * 1024 ) ) ) {
			// Fix for loading SPIDER files which are larger than usual
			// In this specific case, we must absolutely trust the symBytes record
			dataOffset = CCP4HDSIZE + symBytes;
			TR.Warning << "File is larger than expected and doesn't match header.  Reading anyway." << std::endl;
		} else {
			TR.Error << "File is MUCH larger than expected and doesn't match header.  Not loading map." << std::endl;
			utility_exit_with_message( "Fail to load the density map." );
		}
	}

	// Read symmetry records -- organized as 80-byte lines of text.
	utility::vector1< std::string > symList;
	symData[80] = '\0';

	if ( symBytes != 0 ) {
		TR << "Symmetry records found:" << std::endl;
		mapin.seekg( CCP4HDSIZE, std::ios::beg );

		for ( int i = 0; i < symBytes / 80; i++ ) {
			mapin.read( symData, 80 );
			symList.push_back( symData );
			TR << symData << std::endl;
		}
	} else {
		// no symm info; assume P 1
		symList.push_back( "X,  Y,  Z" );
	}

	initializeSymmOps( symList );

	// check extent and grid interval counts
	if ( grid[0] == 0 && extent[0] > 0 ) {
		grid[0] = extent[0] - 1;
		TR.Warning << "Fixed X interval count.  Continuing." << std::endl;
	}

	if ( grid[1] == 0 && extent[1] > 0 ) {
		grid[1] = extent[1] - 1;
		TR.Warning << "Fixed Y interval count.  Continuing." << std::endl;
	}

	if ( grid[2] == 0 && extent[2] > 0 ) {
		grid[2] = extent[2] - 1;
		TR.Warning << "Fixed Z interval count.  Continuing." << std::endl;
	}

	// Mapping between CCP4 column, row, section and Cartesian x, y, z.
	if ( crs2xyz[0] == 0 && crs2xyz[1] == 0 && crs2xyz[2] == 0 ) {
		TR.Warning << "All crs2xyz records are zero." << std::endl;
		TR.Warning << "Setting crs2xyz to 1, 2, 3 and continuing." << std::endl;
		crs2xyz[0] = 1;
		crs2xyz[1] = 2;
		crs2xyz[2] = 3;
	}

	xyz2crs[crs2xyz[0] - 1] = 0;
	xyz2crs[crs2xyz[1] - 1] = 1;
	xyz2crs[crs2xyz[2] - 1] = 2;
	xIndex = xyz2crs[0];
	yIndex = xyz2crs[1];
	zIndex = xyz2crs[2];
	vol_xsize = extent[xIndex];
	vol_ysize = extent[yIndex];
	vol_zsize = extent[zIndex];
	vol_xySize = vol_xsize * vol_ysize;
	// coord = <col, row, sec>
	// extent = <colSize, rowSize, secSize>
	rowdata = new float[extent[0]];
	mapin.seekg( dataOffset, std::ios::beg );
	// 'alloc' changes ordering of "extent"
	density.dimension( vol_xsize, vol_ysize, vol_zsize );

	for ( coord[2] = 1; coord[2] <= extent[2]; coord[2]++ ) {
		for ( coord[1] = 1; coord[1] <= extent[1]; coord[1]++ ) {
			// Read an entire row of data from the file, then write it into the
			// datablock with the correct slice ordering.
			if ( mapin.eof() ) {
				TR.Error << "Unexpected end-of-file. Not loading map." << std::endl;
				utility_exit_with_message( "Fail to load the density map." );
			}

			if ( mapin.fail() ) {
				TR.Error << "Problem reading the file. Not loading map." << std::endl;
				utility_exit_with_message( "Fail to load the density map." );
			}

			if ( !mapin.read( reinterpret_cast< char* >( rowdata ), sizeof( float ) *extent[0] ) ) {
				TR.Error << "Error reading data row. Not loading map." << std::endl;
				utility_exit_with_message( "Fail to load the density map." );
			}

			for ( coord[0] = 1; coord[0] <= extent[0]; coord[0]++ ) {
				density( coord[xyz2crs[0]], coord[xyz2crs[1]], coord[xyz2crs[2]] ) = rowdata[coord[0] - 1];
			}
		}
	}

	if ( swap == 1 ) {
		swap4_aligned( &density[0], vol_xySize * vol_zsize );
	}

	delete [] rowdata;
	mapin.close();
	orig[0] = origin_xyz[xyz2crs[0]];
	orig[1] = origin_xyz[xyz2crs[1]];
	orig[2] = origin_xyz[xyz2crs[2]];
	// grid doesnt seemed to get remapped in ccp4 maps
	///////////////////////////////////
	/// POST PROCESSING
	// expand to unit cell
	computeCrystParams();

	//fpd  change this so if the alt origin is non-zero use it
	if ( altorigin[0] != 0 &&  altorigin[1] != 0 &&  altorigin[2] != 0 &&
			( altorigin[0] > -10000 && altorigin[0] < 10000 ) &&
			( altorigin[1] > -10000 && altorigin[1] < 10000 ) &&
			( altorigin[2] > -10000 && altorigin[2] < 10000 )
			) {
		orig[0] = altorigin[xyz2crs[0]];
		orig[1] = altorigin[xyz2crs[1]];
		orig[2] = altorigin[xyz2crs[2]];
		numeric::xyzVector<core::Real> fracX = c2f * ( orig );
		orig = numeric::xyzVector<core::Real>( fracX[0] * grid[0] , fracX[1] * grid[1] , fracX[2] * grid[2] );
		TR << "Using ALTERNATE origin\n";
		TR << "     origin =" << orig[0] << " x " << orig[1] << " x " << orig[2] << std::endl;
	}

	expandToUnitCell();

	// resample the map
	if ( grid_spacing > 0 ) resize( grid_spacing );

	// grid spacing in each dim
	max_del_grid = std::max( cell_dimensions[0] / ( ( double ) grid[0] ) , cell_dimensions[1] / ( ( double ) grid[1] ) );
	max_del_grid = std::max( max_del_grid , cell_dimensions[2] / ( ( double ) grid[2] ) );
	calculate_index2cart();
	is_map_loaded_ = true;
}

//////////////////////////////////////////////////////////////////////////
//compute index to cartesian transformation
void
ElectronDensityAtomwise::calculate_index2cart() {
	auto i2f = numeric::xyzMatrix<core::Real>::rows(
		1 / double ( grid[0] ), 0, 0,
		0, 1 / double ( grid[1] ), 0,
		0, 0, 1 / double ( grid[2] ) );
	auto f2i = numeric::xyzMatrix<core::Real>::rows(
		grid[0], 0, 0,
		0, grid[1], 0,
		0, 0, grid[2] );
	i2c = f2c * i2f;
	c2i = f2i * c2f;
}

///////////////////////////////////////////////////////////////////////
//Convert a vector from fractional coordinate to index coordinate
numeric::xyzVector< core::Real >
ElectronDensityAtomwise::frac2index( numeric::xyzVector< core::Real > const & frac_vector ) {
	numeric::xyzVector< core::Real > index_vector;
	index_vector[0] = frac_vector[0] * grid[0];
	index_vector[1] = frac_vector[1] * grid[1];
	index_vector[2] = frac_vector[2] * grid[2];
	return index_vector;
}

//Convert a vector from index coordinate to fractional coordinate
numeric::xyzVector< core::Real >
ElectronDensityAtomwise::index2frac( numeric::xyzVector< core::Real > const &index_vector ) {
	numeric::xyzVector< core::Real > frac_vector;
	frac_vector[0] = index_vector[0] / grid[0];
	frac_vector[1] = index_vector[1] / grid[1];
	frac_vector[2] = index_vector[2] / grid[2];
	return frac_vector;
}
numeric::xyzVector< core::Real >
ElectronDensityAtomwise::index2frac( numeric::xyzVector< int > const & index_vector ) {
	numeric::xyzVector< core::Real > frac_vector;
	frac_vector[0] = ( double )  index_vector[0] / grid[0];
	frac_vector[1] = ( double )  index_vector[1] / grid[1];
	frac_vector[2] = ( double )  index_vector[2] / grid[2];
	return frac_vector;
}

//Convert a vector from xyz coordinate to index coordinate, shift w/ respect to the origin and fold into the unit cell
numeric::xyzVector< core::Real >
ElectronDensityAtomwise::xyz2index_in_cell( numeric::xyzVector< core::Real > const & xyz_vector ) {
	numeric::xyzVector< Real > index_vector;
	index_vector = c2i * xyz_vector;
	index_vector[0] = pos_mod( index_vector[0] - orig[0], ( double ) grid[0] );
	index_vector[1] = pos_mod( index_vector[1] - orig[1], ( double ) grid[1] );
	index_vector[2] = pos_mod( index_vector[2] - orig[2], ( double ) grid[2] );
	return index_vector;
}

//////////////////////////////////////////////////////////////////////////
//generate 1d gaussian function and store it.
void
ElectronDensityAtomwise::generate_gaussian_1d( core::Real const & sigma ) {
	const Real PI = numeric::constants::f::pi;
	//Calculate up to 4*sigma
	gaussian_max_d = 4 * sigma;
	//Using interval of 0.01 Angstong
	auto gaussian_size = Size( std::ceil( gaussian_max_d / 0.01 ) );

	for ( Size i = 0; i <= gaussian_size; ++i ) {
		Real calc_value = exp( -( i * 0.01 ) * ( i * 0.01 ) / ( 2.0 * sigma * sigma ) ) / ( sqrt( 2.0 * PI ) * sigma );
		atom_gaussian_value.push_back ( calc_value );
	}
}

//Return the value of 1D gaussian given the distance using stored values
core::Real
ElectronDensityAtomwise::gaussian_1d( core::Real const & dist ) {
	if ( dist >= gaussian_max_d ) return 0.0;

	return atom_gaussian_value[ ( int ) std::floor( dist*100.0 + 0.5 ) +1];
}
///////////////////////////////////////////////////
//return the weight of atom given its element type
core::Size
ElectronDensityAtomwise::get_atom_weight( std::string const & elt ) {
	if ( elt == "C" ) return 6;
	if ( elt == "N" ) return 7;
	if ( elt == "O" ) return 8;
	if ( elt == "P" ) return 15;
	if ( elt == "S" ) return 16;
	if ( elt == "X" ) return 0; // centroid

	// default to C
	TR.Warning << "Unknown atom " << elt << std::endl;
	return 6;
}
//////////////////////////////////////////////////////////////////
//Trilinear Interpolation
core::Real
ElectronDensityAtomwise::trilinear_interpolation(
	ObjexxFCL::FArray3D< double > & score,
	numeric::xyzVector< core::Real > const & index ) {
	int lower_bound[3], upper_bound[3];
	double residual[3];

	// find bounding grid points
	for ( int i = 0; i < 3; ++i ) {
		lower_bound[i] = ( int )( std::ceil( index[i] ) );
		upper_bound[i] = lower_bound[i] + 1;

		if ( upper_bound[i] > grid[i] ) upper_bound[i] = 1;

		residual[i] = index[i] + 1 - lower_bound[i];
	}

	double &c000 = score( lower_bound[0], lower_bound[1], lower_bound[2] );
	double &c100 = score( upper_bound[0], lower_bound[1], lower_bound[2] );
	double &c010 = score( lower_bound[0], upper_bound[1], lower_bound[2] );
	double &c001 = score( lower_bound[0], lower_bound[1], upper_bound[2] );
	double &c110 = score( upper_bound[0], upper_bound[1], lower_bound[2] );
	double &c101 = score( upper_bound[0], lower_bound[1], upper_bound[2] );
	double &c011 = score( lower_bound[0], upper_bound[1], upper_bound[2] );
	double &c111 = score( upper_bound[0], upper_bound[1], upper_bound[2] );
	//interpolate at x dimension (R is the residual coord)
	Real cR00 = c000 + residual[0] * ( c100 - c000 );
	Real cR10 = c010 + residual[0] * ( c110 - c010 );
	Real cR01 = c001 + residual[0] * ( c101 - c001 );
	Real cR11 = c011 + residual[0] * ( c111 - c011 );
	//interpolate at y dimension
	Real cRR0 = cR00 + residual[1] * ( cR10 - cR00 );
	Real cRR1 = cR01 + residual[1] * ( cR11 - cR01 );
	//interpolate at z dimension
	Real cRRR = cRR0 + residual[2] * ( cRR1 - cRR0 );
	return cRRR;
}

numeric::xyzVector<core::Real>
ElectronDensityAtomwise::trilinear_gradient(
	ObjexxFCL::FArray3D< double > & score,
	numeric::xyzVector< core::Real > const & index ) {
	int lower_bound[3], upper_bound[3];
	double residual[3];
	numeric::xyzVector<core::Real> grad;

	// find bounding grid points
	for ( int i = 0; i < 3; ++i ) {
		lower_bound[i] = ( int )( floor( index[i] + 1 ) );
		upper_bound[i] = lower_bound[i] + 1;

		if ( upper_bound[i] > grid[i] ) upper_bound[i] = 1;

		residual[i] = index[i] + 1 - lower_bound[i];
	}

	double &c000 = score( lower_bound[0], lower_bound[1], lower_bound[2] );
	double &c100 = score( upper_bound[0], lower_bound[1], lower_bound[2] );
	double &c010 = score( lower_bound[0], upper_bound[1], lower_bound[2] );
	double &c001 = score( lower_bound[0], lower_bound[1], upper_bound[2] );
	double &c110 = score( upper_bound[0], upper_bound[1], lower_bound[2] );
	double &c101 = score( upper_bound[0], lower_bound[1], upper_bound[2] );
	double &c011 = score( lower_bound[0], upper_bound[1], upper_bound[2] );
	double &c111 = score( upper_bound[0], upper_bound[1], upper_bound[2] );
	//interpolate at x dimension (x)  (R is the residual coord)
	Real cR00 = c000 + residual[0] * ( c100 - c000 );
	Real cR10 = c010 + residual[0] * ( c110 - c010 );
	Real cR01 = c001 + residual[0] * ( c101 - c001 );
	Real cR11 = c011 + residual[0] * ( c111 - c011 );
	//interpolate at y dimension (x->y)
	Real cRR0 = cR00 + residual[1] * ( cR10 - cR00 );
	Real cRR1 = cR01 + residual[1] * ( cR11 - cR01 );
	//Get the z gradient (x->y->z)
	grad[2] = cRR1 - cRR0;
	//interpolate at z dimension (x->z)
	Real cR0R = cR00 + residual[2] * ( cR01 - cR00 );
	Real cR1R = cR10 + residual[2] * ( cR11 - cR10 );
	//Get the y gradient (x->z->y)
	grad[1] = cR1R - cR0R;
	//interpolate at y dimension (y)
	Real c0R0 = c000 + residual[1] * ( c010 - c000 );
	Real c1R0 = c100 + residual[1] * ( c110 - c100 );
	Real c0R1 = c001 + residual[1] * ( c011 - c001 );
	Real c1R1 = c101 + residual[1] * ( c111 - c101 );
	//interpolate at z dimension (y->z)
	Real c0RR = c0R0 + residual[2] * ( c0R1 - c0R0 );
	Real c1RR = c1R0 + residual[2] * ( c1R1 - c1R0 );
	//get the x gradient (y->z->x)
	grad[0] = c1RR - c0RR;
	return grad;
}

///////////////////////////////////////////////////////////////////
//Spline interpolation
core::Real
ElectronDensityAtomwise::spline_interpolation(
	ObjexxFCL::FArray3D< double > & coeffs ,
	numeric::xyzVector< core::Real > const & idxX
) const {
	int dims[3] = { coeffs.u3(), coeffs.u2(), coeffs.u1() };
	core::Real pt[3] = {idxX[2], idxX[1], idxX[0]};
	core::Real retval = core::scoring::electron_density::SplineInterp::interp3( &coeffs[0], dims, pt );
	return retval;
}

void
ElectronDensityAtomwise::spline_coeffs( ObjexxFCL::FArray3D< double > & data ,
	ObjexxFCL::FArray3D< double > & coeffs ) {
	int dims[3] = { data.u3(), data.u2(), data.u1() };
	coeffs = data;
	core::scoring::electron_density::SplineInterp::compute_coefficients3( &coeffs[0] , dims );
}
//////////////////////////////////////////////////////////////////
//Compute Normalization factor given a pose
void
ElectronDensityAtomwise::compute_normalization ( pose::Pose const & pose ) {
	if ( is_score_precomputed_ ) return;
	atom_weight_stored.clear();

	//rho_calc_array initialization
	ObjexxFCL::FArray3D< double >  rho_calc_array;
	rho_calc_array.dimension( density.u1() , density.u2() , density.u3() );

	for ( int i = 0; i < density.u1() *density.u2() *density.u3(); ++i ) {
		rho_calc_array[i] = 0.0;
	}

	//Generate gaussian function for the electron density of the atoms.
	Real map_reso = basic::options::option[ basic::options::OptionKeys::edensity::mapreso ]();
	Real atom_gaussian_sigma = 0.30 + 0.18 * map_reso;
	generate_gaussian_1d( atom_gaussian_sigma );
	//Calculate the weight and the fractional coordinates of the atoms
	Size n_res = pose.size();
	Size sum_weight = 0;

	numeric::xyzVector< core::Real > grid_half( grid[0] * 0.5, grid[1] * 0.5, grid[2] * 0.5 );

	atom_weight_stored.reserve( n_res );
	for ( Size i = 1; i <= n_res; ++i ) { //rsd
		conformation::Residue const & rsd( pose.residue( i ) );
		Size n_heavyatom_rsd = rsd.nheavyatoms();
		utility::vector1< Size > weight_rsd;

		// skip virtual residues
		if ( rsd.aa() == core::chemical::aa_vrt ) {
			atom_weight_stored.push_back( weight_rsd );
			continue;
		}

		weight_rsd.reserve( n_heavyatom_rsd );
		for ( Size j = 1; j <= n_heavyatom_rsd; ++j ) { //atom
			chemical::AtomTypeSet const & atom_type_set( rsd.atom_type_set() );
			//get the weight of each atom (element based)
			std::string const & element = atom_type_set[ rsd.atom_type_index( j ) ].element();
			Size weight = get_atom_weight( element );
			sum_weight += weight;
			weight_rsd.push_back( weight );
			//Compute the atom position in unit cell (index coord)
			numeric::xyzVector< core::Real > const & coord_index = xyz2index_in_cell( rsd.xyz( j ) );
			numeric::xyzVector< core::Real > dist_index;

			//Compute rho_calc
			for ( int x = 1; x <= density.u1(); ++x ) {
				//distance from the atom to point in the map (x)
				dist_index[0] = coord_index[0] - ( x - 1 );

				//fold to make it in the range [-0.5*grid, 0.5*grid]
				if ( dist_index[0] > grid_half[0] ) dist_index[0] -= grid[0];
				if ( dist_index[0] <= -grid_half[0] ) dist_index[0] += grid[0];

				//skip if the dist is too large
				numeric::xyzVector<Real> dist_x = i2c * numeric::xyzVector<Real>( dist_index[0], 0.0, 0.0 );
				if ( dist_x.length() > gaussian_max_d ) continue;

				for ( int y = 1; y <= density.u2(); ++y ) {
					//distance from the atom to point in the map (y)
					dist_index[1] = coord_index[1] - ( y - 1 );

					//fold to make it in the range [-0.5*grid, 0.5*grid]
					if ( dist_index[1] > grid_half[1] ) dist_index[1] -= grid[1];
					if ( dist_index[1] <= - grid_half[1] ) dist_index[1] += grid[1];

					//skip if the dist is too large
					numeric::xyzVector<Real> dist_xy = i2c * numeric::xyzVector<Real>( dist_index[0], dist_index[1], 0.0 );
					if ( dist_xy.length() > gaussian_max_d ) continue;

					for ( int z = 1; z <= density.u3(); ++z ) {
						//distance from the atom to point in the map (z)
						dist_index[2] = coord_index[2] - ( z - 1 );

						//fold to make it in the range [-0.5*grid, 0.5*grid]
						if ( dist_index[2] > grid_half[2] ) dist_index[2] -= grid[2];
						if ( dist_index[2] <= - grid_half[2] ) dist_index[2] += grid[2];

						//skip if the dist is too large
						numeric::xyzVector<Real> dist_xyz = i2c * dist_index;
						Real dist = dist_xyz.length();

						if ( dist > gaussian_max_d ) continue;

						rho_calc_array( x, y, z ) += gaussian_1d( dist ) * weight;
					}
				}
			}
		}

		atom_weight_stored.push_back( weight_rsd );
	}

	Real sum_rho_calc = 0.0;
	Real sum_rho_calc_square = 0.0;
	Real sum_rho_obs = 0.0;
	Real sum_rho_obs_square = 0.0;
	Real sum_rho_calc_times_rho_obs = 0.0;
	Size n_grid_points = 0;

	for ( int i = 0; i < density.u1() *density.u2() *density.u3(); ++i ) {
		float rho_obs = density[i];
		Real rho_calc = rho_calc_array[i];

		if ( rho_calc < 1.0e-2 ) continue;

		//compute all factors required for normalization
		n_grid_points++;
		sum_rho_calc += rho_calc;
		sum_rho_calc_square += rho_calc * rho_calc;
		sum_rho_obs += rho_obs;
		sum_rho_obs_square += rho_obs * rho_obs;
		sum_rho_calc_times_rho_obs += rho_calc * rho_obs;
	}

	//calculate the real-space correlation coefficient (rscc)
	//rscc = (<xy> - <x><y>) / sqrt( (<x^2> -<x>^2) * (<y^2> -<y>^2) )
	//     = <dx dy> / sqrt( <dx dx> + <dy dy> )
	Real rscc = ( sum_rho_calc_times_rho_obs - sum_rho_calc * sum_rho_obs / n_grid_points ) / ( sqrt( sum_rho_calc_square - sum_rho_calc * sum_rho_calc / n_grid_points ) * sqrt( sum_rho_obs_square - sum_rho_obs * sum_rho_obs / n_grid_points ) );
	//Compute Normalization, Score = -sum_weight * rscc = normalization * <dx dy>
	//Per atom energy = normalization * ( sum(rho_atom*rho_obs) - sum(rho_atom) * avg_rho_obs )
	normalization = - ( ( double ) sum_weight ) / ( sqrt( sum_rho_calc_square - sum_rho_calc * sum_rho_calc / n_grid_points ) * sqrt( sum_rho_obs_square - sum_rho_obs * sum_rho_obs / n_grid_points ) );
	avg_rho_obs = sum_rho_obs / n_grid_points;

	TR.Debug << "RSCC of starting pose = " << rscc << std::endl;
	TR.Debug << "Score of starting pose = " << -( rscc * sum_weight ) << std::endl;
	TR.Debug << "normalization factor = " << normalization << std::endl;
	TR.Debug << "avg_rho_obs = " << avg_rho_obs << std::endl;
}

//Pre-compute the unweighted score of atom in a grid
void
ElectronDensityAtomwise::precompute_unweighted_score() {
	if ( is_score_precomputed_ ) return;

	ObjexxFCL::FArray3D< double > unweighted_score;
	ObjexxFCL::FArray3D< std::complex<double> > density_transformed,
		atom_dens_transformed, atom_dens, cplx_score, cplx_density;
	atom_dens.dimension( density.u1() , density.u2() , density.u3() );

	for ( int i = 0; i < density.u1() *density.u2() *density.u3(); ++i ) {
		atom_dens[i] = 0.0;
	}

	//compute the gaussian density of an atom
	Real sum_rho_calc = 0.0;
	numeric::xyzVector< core::Real > dist_index;
	numeric::xyzVector< core::Real > grid_half( grid[0] * 0.5, grid[1] * 0.5, grid[2] * 0.5 );

	for ( int x = 1; x <= atom_dens.u1(); ++x ) {
		//distance from the atom to point in the map (x)
		dist_index[0] = x - 1;

		//fold to make it in the range [-0.5*grid, 0.5*grid]
		if ( dist_index[0] > grid_half[0] ) dist_index[0] -= grid[0];
		if ( dist_index[0] <= - grid_half[0] ) dist_index[0] += grid[0];

		//skip if the dist is too large
		numeric::xyzVector<Real> dist_x = i2c * numeric::xyzVector<Real>( dist_index[0], 0.0, 0.0 );

		if ( dist_x.length() > gaussian_max_d ) continue;

		for ( int y = 1; y <= atom_dens.u2(); ++y ) {
			//distance from the atom to point in the map (y)
			dist_index[1] = y - 1;

			//fold to make it in the range [-0.5*grid, 0.5*grid]
			if ( dist_index[1] > grid_half[1] ) dist_index[1] -= grid[1];
			if ( dist_index[1] <= - grid_half[1] ) dist_index[1] += grid[1];

			//skip if the dist is too large
			auto dist_xy = i2c * numeric::xyzVector<Real>( dist_index[0], dist_index[1], 0.0 );

			if ( dist_xy.length() > gaussian_max_d ) continue;

			for ( int z = 1; z <= atom_dens.u3(); ++z ) {
				//distance from the atom to point in the map (y)
				dist_index[2] = z - 1;

				//fold to make it in the range [-0.5*grid, 0.5*grid]
				if ( dist_index[2] > grid_half[2] ) dist_index[2] -= grid[2];
				if ( dist_index[2] <= - grid_half[2] ) dist_index[2] += grid[2];

				numeric::xyzVector<Real> dist_xyz = i2c * dist_index;
				Real dist = dist_xyz.length();
				atom_dens( x, y, z ) = gaussian_1d( dist );
				sum_rho_calc += gaussian_1d( dist );
			}
		}
	}

	//Compute convolution using fft
	numeric::fourier::fft3_dynamic( atom_dens, atom_dens_transformed );
	atom_dens.clear();
	cplx_density.dimension( density.u1() , density.u2() , density.u3() );

	for ( int i = 0; i < density.u1() *density.u2() *density.u3(); ++i ) {
		cplx_density[i] = ( double ) density[i];
	}

	// calebgeniesse: cannot recompute normalization if this is cleared,
	//                and no obvious reason that it needs to be cleared.
	//density.clear();

	numeric::fourier::fft3_dynamic( cplx_density, density_transformed );
	cplx_density.clear();

	for ( int i = 0; i < density_transformed.u1() *density_transformed.u2() *density_transformed.u3(); ++i ) {
		density_transformed[i] *= atom_dens_transformed[i];
	}

	atom_dens_transformed.clear();
	numeric::fourier::ifft3_dynamic( density_transformed, cplx_score );
	density_transformed.clear();
	unweighted_score.dimension( cplx_score.u1() , cplx_score.u2() , cplx_score.u3() );

	for ( int i = 0; i < unweighted_score.u1() *unweighted_score.u2() *unweighted_score.u3(); ++i ) {
		unweighted_score[i] = cplx_score[i].real();
		unweighted_score[i] -=  sum_rho_calc * avg_rho_obs;
		unweighted_score[i] *=  normalization;
	}

	cplx_score.clear();
	//compute spline coefficient for interpolation
	spline_coeffs( unweighted_score, unweighted_score_coeff );
	unweighted_score.clear();
	is_score_precomputed_ = true;
}

///////////////////////////////////////////////////////////////////////
//Return the score of a given residue
core::Real
ElectronDensityAtomwise::residue_score( core::conformation::Residue const & rsd ) {
	// skip virtual residues
	if ( rsd.aa() == core::chemical::aa_vrt ) return 0.0;

	Size n_heavyatom_rsd = rsd.nheavyatoms();
	Size rsd_id = rsd.seqpos();
	Real total_score = 0.0;

	for ( Size j = 1 ; j <= n_heavyatom_rsd; ++j ) {
		//Skip virtual atoms
		if ( rsd.is_virtual(j) ) continue;
		//get the weight of each atom (element based)
		Size weight = atom_weight_stored[rsd_id][j];
		//Compute the atom position in unit cell (index coord)
		numeric::xyzVector< Real > coord_index = xyz2index_in_cell( rsd.xyz( j ) );
		Real score = weight * spline_interpolation( unweighted_score_coeff, coord_index );
		total_score += score;
	}

	return total_score;
}

//Return the gradient of an atom
numeric::xyzVector< core::Real >
ElectronDensityAtomwise::atom_gradient( core::pose::Pose const & pose, core::Size
	const & rsd_id, core::Size const & atm_id ) {

	//Skip virtual atoms
	if ( pose.residue( rsd_id ).is_virtual( atm_id ) ) {
		return numeric::xyzVector< core::Real >( 0, 0, 0 );
	}

	//get the weight of the atom (element based)
	Size &weight = atom_weight_stored[rsd_id][atm_id];
	Real incre = 0.00000001;
	numeric::xyzVector< Real > atmxyz = pose.residue( rsd_id ).xyz( atm_id );
	numeric::xyzVector< Real > coord_index = xyz2index_in_cell( atmxyz );
	numeric::xyzVector< Real > coord_index_dx = coord_index;
	numeric::xyzVector< Real > coord_index_dy = coord_index;
	numeric::xyzVector< Real > coord_index_dz = coord_index;
	coord_index_dx[0] += incre;
	coord_index_dy[1] += incre;
	coord_index_dz[2] += incre;
	numeric::xyzVector< Real > grad;
	Real score_atom = spline_interpolation( unweighted_score_coeff, coord_index );
	grad[0] = ( spline_interpolation( unweighted_score_coeff, coord_index_dx ) -
		score_atom ) / incre;
	grad[1] = ( spline_interpolation( unweighted_score_coeff, coord_index_dy ) -
		score_atom ) / incre;
	grad[2] = ( spline_interpolation( unweighted_score_coeff, coord_index_dz ) -
		score_atom ) / incre;
	numeric::xyzMatrix<core::Real> i2c_gradient = c2i;
	i2c_gradient.transpose();
	grad = weight * i2c_gradient * grad;
	return grad;
}

///////////////////////////////////////////////////////////////////////////
ElectronDensityAtomwise & get_density_map() {
	static ElectronDensityAtomwise theDensityMap;

	if ( !theDensityMap.isMapLoaded() ) {
		// Initialize ElectronDensity object
		theDensityMap.readMRCandResize();
	}

	return theDensityMap;
}
//////////////////////////////////////////////////////////////////////////

} // electron_density_atomwise
} // scoring
} // core

