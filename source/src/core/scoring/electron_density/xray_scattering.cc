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
/// @author

#include <core/scoring/electron_density/xray_scattering.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/edensity.OptionKeys.gen.hh>

#include <iostream>
#include <map>

#include <utility/vector1.hh>


// C++ Headers

static basic::Tracer TR( "core.scoring.electron_density.xray_scattering" );


namespace core {
namespace scoring {
namespace electron_density {

// x mod y, returns z in [0,y-1]
inline int pos_mod(int x,int y) {
	int r=x%y; if ( r<0 ) r+=y;
	return r;
}
inline float pos_mod(float x,float y) {
	float r=std::fmod(x,y); if ( r<0 ) r+=y;
	return r;
}
inline double pos_mod(double x,double y) {
	double r=std::fmod(x,y); if ( r<0 ) r+=y;
	return r;
}
inline float  square(float  x) { return (x*x); }
inline double square(double x) { return (x*x); }

// weight from scattering factors
// single-gaussian approx
OneGaussianScattering get_A( std::string elt ) {
	static std::map< std::string, OneGaussianScattering > elt_db;

	// fit to parameters from Doyle and Turner, Acta Cryst A, 1968
	//   --> single gaussians fit in reciprocal space fom 20A-2A resolution
	//   --> TODO .. move to database, make parameters selectable
	if ( elt_db.empty() ) { //size() == 0) {
		if ( basic::options::option[ basic::options::OptionKeys::edensity::cryoem_scatterers ]() ) {
			// [1] electron scattering
			elt_db["C"]  = OneGaussianScattering(  6.00000, 7.10668);
			elt_db["N"]  = OneGaussianScattering(  5.28737, 6.03448);
			elt_db["O"]  = OneGaussianScattering(  4.74213, 5.17616);
			elt_db["NA"] = OneGaussianScattering( 11.42607, 6.58734);
			elt_db["MG"] = OneGaussianScattering( 12.45197, 7.34364);
			elt_db["P"]  = OneGaussianScattering( 13.12395, 7.48955);
			elt_db["S"]  = OneGaussianScattering( 12.34197, 7.05366);
			elt_db["K"]  = OneGaussianScattering( 21.48425, 7.09360);
			elt_db["CA"] = OneGaussianScattering( 23.70586, 7.47775);
			elt_db["FE"] = OneGaussianScattering( 17.13431, 5.96932);
			elt_db["CO"] = OneGaussianScattering( 15.70905, 5.56662); // copied from NI
			elt_db["NI"] = OneGaussianScattering( 15.70905, 5.56662);
			elt_db["ZN"] = OneGaussianScattering( 15.70905, 5.56662); // copied from NI
			elt_db["X"]  = OneGaussianScattering(
				static_cast<int>(
				6.0*basic::options::option[ basic::options::OptionKeys::edensity::centroid_density_mass ]()),
				4.88284);  // centroid
		} else {
			// [2] xray scattering
			elt_db["C"]  = OneGaussianScattering(  6.0, 4.88398);
			elt_db["N"]  = OneGaussianScattering(  7.0, 5.08400);
			elt_db["O"]  = OneGaussianScattering(  8.0, 4.92866);
			elt_db["NA"] = OneGaussianScattering( 11.0, 3.53009);
			elt_db["MG"] = OneGaussianScattering( 12.0, 3.11369);
			elt_db["P"]  = OneGaussianScattering( 15.0, 2.89121);
			elt_db["S"]  = OneGaussianScattering( 16.0, 3.03431);
			elt_db["K"]  = OneGaussianScattering( 19.0, 3.32509);
			elt_db["CA"] = OneGaussianScattering( 20.0, 3.25370);
			elt_db["FE"] = OneGaussianScattering( 26.0, 3.10464);
			elt_db["CO"] = OneGaussianScattering( 28.0, 2.95223); // copied from NI
			elt_db["NI"] = OneGaussianScattering( 28.0, 2.95223);
			elt_db["ZN"] = OneGaussianScattering( 28.0, 2.95223); // copied from NI
			elt_db["X"]  = OneGaussianScattering(
				static_cast<int>(
				6.0*basic::options::option[ basic::options::OptionKeys::edensity::centroid_density_mass ]()),
				4.88284);  // centroid
		}
	}


	if ( elt_db.find( elt ) == elt_db.end() ) {
		// default to C
		TR.Warning << "Unknown atom " << elt << std::endl;
		return elt_db["C"];
	} else {
		return elt_db[ elt ];
	}
}

// weight from scattering factors
KromerMann get_km( std::string elt ) {
	static std::map< std::string, KromerMann > elt_db;

	if ( elt_db.empty() ) { //size() == 0) {
		elt_db["C"]  = KromerMann(  0.215600,  2.310000, 1.020000, 1.588600, 0.865000, 20.843899, 10.207500,  0.568700,  51.651199);
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

	if ( elt_db.find( elt ) == elt_db.end() ) {
		// default to C
		TR.Warning << "Unknown atom " << elt << std::endl;
		return elt_db["C"];
	} else {
		return elt_db[ elt ];
	}
}

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////


bool factorsLTE5(int X) {
	while ( X != 1 && X%2 == 0 ) X /= 2;
	while ( X != 1 && X%3 == 0 ) X /= 3;
	while ( X != 1 && X%5 == 0 ) X /= 5;

	return (X == 1);
}


bool factorsLTE19(int X) {
	while ( X != 1 && X%2 == 0 ) X /= 2;
	while ( X != 1 && X%3 == 0 ) X /= 3;
	while ( X != 1 && X%5 == 0 ) X /= 5;
	while ( X != 1 && X%7 == 0 ) X /= 7;
	while ( X != 1 && X%11 == 0 ) X /= 11;
	while ( X != 1 && X%13 == 0 ) X /= 13;
	while ( X != 1 && X%17 == 0 ) X /= 17;
	while ( X != 1 && X%19 == 0 ) X /= 19;

	return (X == 1);
}


int findSampling5(double MINSMP, int NMUL) {
	if ( MINSMP <= 0 ) return NMUL;

	// multiple of nmul nearest minsmp
	int N = (int) floor( MINSMP/NMUL + 0.5 ) * NMUL;

	// increment until no factors >= 5
	while ( !factorsLTE5(N) )
			N += NMUL;

	return N;
}


int findSampling(double MINSMP, int NMUL) {
	if ( MINSMP <= 0 ) return NMUL;

	// multiple of nmul nearest minsmp
	int N = (int) floor( MINSMP/NMUL + 0.5 ) * NMUL;

	// increment until no factors >= 19
	while ( !factorsLTE19(N) )
			N += NMUL;

	return N;
}


} // namespace constraints
} // namespace scoring
} // namespace core
