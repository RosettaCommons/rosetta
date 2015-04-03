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

#include <core/scoring/fiber_diffraction/xray_scattering.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>


#include <iostream>
#include <map>

#include <utility/vector1.hh>


// C++ Headers

using basic::T;
using basic::Tracer;
static basic::Tracer TR("core.scoring.fiber_diffraction.xray_scattering");


namespace core {
namespace scoring {
namespace fiber_diffraction {

// x mod y, returns z in [0,y-1]
inline int pos_mod(int x,int y) {
	int r=x%y; if (r<0) r+=y;
	return r;
}
inline float pos_mod(float x,float y) {
	float r=std::fmod(x,y); if (r<0) r+=y;
	return r;
}
inline double pos_mod(double x,double y) {
	double r=std::fmod(x,y); if (r<0) r+=y;
	return r;
}
inline float  square(float  x) { return (x*x); }
inline double square(double x) { return (x*x); }

// weight from scattering factors
// single-gaussian approx
OneGaussianScattering get_A( std::string elt ) {
	static std::map< std::string, OneGaussianScattering > elt_db;

	//fpd these parameters were tuned by fitting a single gaussian
	//    to the 4-gaussian parameters in R^-1 space from 2-20A res
	if (elt_db.size() == 0) {
		elt_db["C"]  = OneGaussianScattering( 6, 4.88284);
		elt_db["N"]  = OneGaussianScattering( 7, 5.08287);
		elt_db["O"]  = OneGaussianScattering( 8, 4.92989);
		elt_db["F"]  = OneGaussianScattering( 9, 4.58837);
		elt_db["NA"] = OneGaussianScattering(11, 3.53389);
		elt_db["MG"] = OneGaussianScattering(12, 3.11868);
		elt_db["AL"] = OneGaussianScattering(13, 2.87655);
		elt_db["P"]  = OneGaussianScattering(15, 2.89121);
		elt_db["S"]  = OneGaussianScattering(16, 3.03242);
		elt_db["K"]  = OneGaussianScattering(19, 3.32675);
		elt_db["FE"] = OneGaussianScattering(26, 3.10519);

		elt_db["X"] = OneGaussianScattering(
			(int)(6*basic::options::option[ basic::options::OptionKeys::score::fiber_diffraction::centroid_density_mass ]()) , 4.88284);  // centroid
	}

	if (elt_db.find( elt ) == elt_db.end()) {
		// default to C
		TR.Warning << "[ WARNING ] Unknown atom " << elt << std::endl;
		return elt_db["C"];
	} else {
		return elt_db[ elt ];
	}
}

// weight from scattering factors
KromerMann get_km( std::string elt ) {
	static std::map< std::string, KromerMann > elt_db;

	if (elt_db.size() == 0) {
		//elt_db["C"]  = KromerMann(  0.215600,  2.310000, 1.020000, 1.588600, 0.865000, 20.843899, 10.207500,  0.568700,  51.651199);
		elt_db["C"]  = KromerMann(  0.286977, 2.26069, 1.05075, 1.56165, 0.839259, 22.6907, 9.75618, 0.656665, 55.5949);
		elt_db["N"]  = KromerMann(-11.528999, 12.212600, 3.132200, 2.012500, 1.166300,  0.005700,  9.893300, 28.997499,   0.582600);
		elt_db["O"]  = KromerMann(  0.250800,  3.048500, 2.286800, 1.546300, 0.867000, 13.277100,  5.701100,  0.323900,  32.908897);
		elt_db["F"]  = KromerMann(  0.277600,  3.539200, 2.641200, 1.517000, 1.024300, 10.282499,  4.294400,  0.261500,  26.147600);
		elt_db["NA"] = KromerMann(  0.676000,  4.762600, 3.173600, 1.267400, 1.112800,  3.285000,  8.842199,  0.313600, 129.423996);
		elt_db["MG"] = KromerMann(  0.858400,  5.420400, 2.173500, 1.226900, 2.307300,  2.827500, 79.261101,  0.380800,   7.193700);
		elt_db["AL"] = KromerMann(  1.115100, 6.420200, 1.900200, 1.593600, 1.964600, 3.038700,  0.742600, 31.547199, 85.088600);
		elt_db["P"]  = KromerMann(  1.114900, 6.434500, 4.179100, 1.780000, 1.490800, 1.906700, 27.157000,  0.526000, 68.164497);
		elt_db["S"]  = KromerMann(  0.866900, 6.905300, 5.203400, 1.437900, 1.586300, 1.467900, 22.215099,  0.253600, 56.172001);
		elt_db["K"]  = KromerMann(  1.422800, 8.218599 ,7.439800, 1.051900, 0.865900,12.794900, 0.774800 ,213.186996, 41.684097);
		elt_db["FE"] = KromerMann(  1.036900, 11.769500, 7.357300, 3.522200, 2.304500, 4.761100, 0.307200, 15.353500, 76.880501 );
		elt_db["X"] = KromerMann();  // centroid
	}

	/* XPLOR parameters
	if (elt_db.size() == 0) {
                elt_db["C"]  = KromerMann(  2.26069, 22.6907, 1.56165, 0.656665, 1.05075, 9.75618, 0.839259, 55.5949, 0.286977);
                elt_db["N"]  = KromerMann(  12.2126, 0.005700, 3.13220, 9.89330, 2.01250, 28.9975, 1.16630, 0.582600, -11.529);
                elt_db["O"]  = KromerMann(  3.04850, 13.2771, 2.28680, 5.70110, 1.54630, 0.323900, 0.867000, 32.9089, 0.250800);
                elt_db["F"]  = KromerMann(  0.277600,  3.539200, 2.641200, 1.517000, 1.024300, 10.282499,  4.294400,  0.261500,  26.147600);
                elt_db["NA"] = KromerMann(  0.676000,  4.762600, 3.173600, 1.267400, 1.112800,  3.285000,  8.842199,  0.313600, 129.423996);
                elt_db["MG"] = KromerMann(  0.858400,  5.420400, 2.173500, 1.226900, 2.307300,  2.827500, 79.261101,  0.380800,   7.193700);
                elt_db["AL"] = KromerMann(  1.115100, 6.420200, 1.900200, 1.593600, 1.964600, 3.038700,  0.742600, 31.547199, 85.088600);
                elt_db["P"]  = KromerMann(  1.114900, 6.434500, 4.179100, 1.780000, 1.490800, 1.906700, 27.157000,  0.526000, 68.164497);
                elt_db["S"]  = KromerMann(  6.90530, 1.46790, 5.20340, 22.2151, 1.43790, 0.253600, 1.58630, 56.1720, 0.866900);
                elt_db["K"]  = KromerMann(  1.422800, 8.218599 ,7.439800, 1.051900, 0.865900,12.794900, 0.774800 ,213.186996, 41.684097);
                elt_db["FE"] = KromerMann(  1.036900, 11.769500, 7.357300, 3.522200, 2.304500, 4.761100, 0.307200, 15.353500, 76.880501 );
                elt_db["X"] = KromerMann();  // centroid
        }*/


	if (elt_db.find( elt ) == elt_db.end()) {
		// default to C
		TR.Warning << "[ WARNING ] Unknown atom " << elt << std::endl;
		return elt_db["C"];
	} else {
		return elt_db[ elt ];
	}
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////


// AtomScattering::AtomScattering() { }
//
//
// // mask scattering
// AtomScattering::AtomScattering( core::Real mask, core::Real reso ) {
// 	init( mask, reso );
// 	numeric::xyzVector< core::Real > atm_j;
//
// 	// ??? should there be a resolution-based sigma on the mask?
//
// 	for (int z=1; z<=data.u3(); ++z) {
// 		atm_j[2] = (z-1)*i2c[2];
// 		for (int y=1; y<=data.u2(); ++y) {
// 			atm_j[1] = (y-1)*i2c[1];
// 			for (int x=1; x<=data.u1(); ++x) {
// 				atm_j[0] = (x-1)*i2c[0];
// 				core::Real d2 = atm_j.length_squared();
// 				core::Real inv_msk = 1/(1+exp( d2 - (mask*mask) ));
// 				data(x,y,z) = inv_msk;
// 			}
// 		}
// 	}
// }
//
// // atom scattering
// AtomScattering::AtomScattering( core::Real a, core::Real B, core::Real mask, core::Real reso ) {
// 	init( mask, reso );
//
//   // 1-gaussian approx
// 	core::Real k = (B<=0) ? square(M_PI/reso) : std::min( square(M_PI/reso) , 0.5/square(0.6+0.006*B)  );
// 	core::Real C = pow(k/M_PI,1.5);
// 	numeric::xyzVector< core::Real > atm_j;
//
// 	for (int z=1; z<=data.u3(); ++z) {
// 		atm_j[2] = (z-1)*i2c[2];
// 		for (int y=1; y<=data.u2(); ++y) {
// 			atm_j[1] = (y-1)*i2c[1];
// 			for (int x=1; x<=data.u1(); ++x) {
// 				atm_j[0] = (x-1)*i2c[0];
// 				core::Real d2 = atm_j.length_squared();
// 				core::Real atm = C*a*exp(-k*d2);
// 				data(x,y,z) = atm;
// 			}
// 		}
// 	}
// }
//
//
// // common init
// void AtomScattering::init( core::Real mask, core::Real reso ) {
// 	this->reso = reso;
//
// 	// extent is max(mask+1,3xreso) <<< should this be a function of atom mask
// 	core::Real extent = std::max( mask+1.0, 3*reso );
//
// 	// grid sampling (A/pt) 1/8 reso
// 	core::Real sampling = reso/8;
// 	core::Size ngrid = (core::Size) std::floor(extent/sampling+0.5);
// 	sampling = extent / (double) ngrid;
//
// 	grid[0] = grid[1] = grid[2] = ngrid;
// 	c2i[0] = c2i[1] = c2i[2] = 1/sampling;
// 	i2c[0] = i2c[1] = i2c[2] = sampling;
// 	data.dimension( (int)grid[0] , (int)grid[1] , (int)grid[2] );
// }
//
//
// // interp at a Cartesian point X
// core::Real AtomScattering::interp_linear( numeric::xyzVector< core::Real > const & X ) const {
// 	numeric::xyzVector< core::Real > idxX(X[0]*c2i[0],X[1]*c2i[1],X[2]*c2i[2]);  // cart -> idx
//
// 	// only compute density over one octant
// 	if (idxX[0] < 0) idxX[0] *= -1;
// 	if (idxX[1] < 0) idxX[1] *= -1;
// 	if (idxX[2] < 0) idxX[2] *= -1;
//
// 	int pt000[3], pt111[3];
// 	core::Real fpart[3],neg_fpart[3];
//
// 	// check for out of bounds
// 	if (idxX[0] >= grid[0] || idxX[1] >= grid[1] || idxX[2] >= grid[2]) return 0.0;
//
// 	// find bounding grid points
// 	pt000[0] = (int)(floor(idxX[0]))+1; pt111[0] = (pt000[0]+1);
// 	pt000[1] = (int)(floor(idxX[1]))+1; pt111[1] = (pt000[1]+1);
// 	pt000[2] = (int)(floor(idxX[2]))+1; pt111[2] = (pt000[2]+1);
//
// 	// interpolation coeffs
// 	fpart[0] = idxX[0]-floor(idxX[0]); neg_fpart[0] = 1-fpart[0];
// 	fpart[1] = idxX[1]-floor(idxX[1]); neg_fpart[1] = 1-fpart[1];
// 	fpart[2] = idxX[2]-floor(idxX[2]); neg_fpart[2] = 1-fpart[2];
//
// debug_assert( pt000[0] >= 1 && pt000[0] <= data.u1() );
// debug_assert( pt000[1] >= 1 && pt000[1] <= data.u2() );
// debug_assert( pt000[2] >= 1 && pt000[2] <= data.u3() );
// debug_assert( pt111[0] >= 1 && pt111[0] <= data.u1() );
// debug_assert( pt111[1] >= 1 && pt111[1] <= data.u2() );
// debug_assert( pt111[2] >= 1 && pt111[2] <= data.u3() );
//
//
// 	// interpolate
// 	core::Real retval = 0.0;
// 	retval += neg_fpart[0]*neg_fpart[1]*neg_fpart[2] * data(pt000[0],pt000[1],pt000[2]);
// 	retval += neg_fpart[0]*neg_fpart[1]*    fpart[2] * data(pt000[0],pt000[1],pt111[2]);
// 	retval += neg_fpart[0]*    fpart[1]*neg_fpart[2] * data(pt000[0],pt111[1],pt000[2]);
// 	retval += neg_fpart[0]*    fpart[1]*    fpart[2] * data(pt000[0],pt111[1],pt111[2]);
// 	retval += fpart[0]*neg_fpart[1]*neg_fpart[2] * data(pt111[0],pt000[1],pt000[2]);
// 	retval += fpart[0]*neg_fpart[1]*    fpart[2] * data(pt111[0],pt000[1],pt111[2]);
// 	retval += fpart[0]*    fpart[1]*neg_fpart[2] * data(pt111[0],pt111[1],pt000[2]);
// 	retval += fpart[0]*    fpart[1]*    fpart[2] * data(pt111[0],pt111[1],pt111[2]);
//
// 	return retval;
// }


/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////


bool factorsLTE5(int X) {
	while (X != 1 && X%2 == 0) X /= 2;
	while (X != 1 && X%3 == 0) X /= 3;
	while (X != 1 && X%5 == 0) X /= 5;

	return (X == 1);
}


bool factorsLTE19(int X) {
	while (X != 1 && X%2 == 0) X /= 2;
	while (X != 1 && X%3 == 0) X /= 3;
	while (X != 1 && X%5 == 0) X /= 5;
	while (X != 1 && X%7 == 0) X /= 7;
	while (X != 1 && X%11 == 0)	X /= 11;
	while (X != 1 && X%13 == 0)	X /= 13;
	while (X != 1 && X%17 == 0)	X /= 17;
	while (X != 1 && X%19 == 0)	X /= 19;

	return (X == 1);
}


int findSampling5(double MINSMP, int NMUL) {
	if (MINSMP <= 0) return NMUL;

	// multiple of nmul nearest minsmp
	int N = (int) floor( MINSMP/NMUL + 0.5 ) * NMUL;

	// increment until no factors >= 5
	while (!factorsLTE5(N))
		N += NMUL;

	return N;
}


int findSampling(double MINSMP, int NMUL) {
	if (MINSMP <= 0) return NMUL;

	// multiple of nmul nearest minsmp
	int N = (int) floor( MINSMP/NMUL + 0.5 ) * NMUL;

	// increment until no factors >= 19
	while (!factorsLTE19(N))
		N += NMUL;

	return N;
}


} // namespace constraints
} // namespace scoring
} // namespace core
