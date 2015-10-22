// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file numeric/model_quality/maxsub.cc
/// @brief ab-initio fragment assembly protocol for proteins
/// @details Routines for calculating maxsub-based structural quality scores. Based on code originally
/// written by Charlie Strauss for rosetta++, ported over by James Thompson.
///
/// @author James Thompson

#include <numeric/model_quality/rms.hh>
#include <numeric/model_quality/maxsub.hh>
#include <numeric/model_quality/RmsData.hh>

// ObjexxFCL Headers
#include <ObjexxFCL/FArray1.hh>              // for FArray1
#include <ObjexxFCL/FArray1D.fwd.hh>         // for FArray1D_double, FArray1...
#include <ObjexxFCL/FArray1D.hh>             // for FArray1D
#include <ObjexxFCL/FArray2.hh>              // for FArray2
#include <ObjexxFCL/Star.hh>                 // for star
#include <numeric/numeric.functions.hh>      // for square
#include <ObjexxFCL/FArray1A.hh>
#include <ObjexxFCL/FArray2A.hh>
//#include <ObjexxFCL/format.hh>

// C++ Headers
#include <cmath>
#include <cstdlib>

namespace numeric {
namespace model_quality {

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
/// cems 2001.
/// this is the main rosetta entry point for this function.
/// it is a wrapper for max sub, converting the input arrays to
///  double precision and reducing them
/// to just calphas
/// it does max sub comparing the claphas passed into thosw of the
///  native and returns
/// the number of aligned residues, the rms of these and the log eval c
/// of the comparison
///
/// @param  x - [in/out]? -
/// @param  nali - [in/out]? -
/// @param  rms - [in/out]? -
/// @param  logeval - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
// void
// maxsub_native(
//  FArray3_float const & x,
//  int & nali,
//  float & rms,
//  float & logeval
// )
// {
//  using namespace misc;
//  using namespace native;
//
// // local
//  FArray2D_double xp( 3, total_residue );
//  FArray2D_double xe( 3, total_residue );
//  double mxrms,mxpsi,mxzscore,mxscore,mxeval;
//
//  if ( !get_native_exists() || ( files_paths::multi_chain && !design::dna_interface ) ) {
//   rms = 0.0;
//   logeval = 0.0;
//   nali = 0;
//   return;
//  }
//
//  int n_points = 0;
//  for ( int i = 1; i <= total_residue; ++i ) {
//   if ( native_occupancy( 2, i ) <= 0.0 ) continue;
//   n_points++;
//   for ( int k = 1; k <= 3; ++k ) {
//    xe(k,n_points) = native_ca(k,i);
//    xp(k,n_points) = x(k,2,i); // calphas
//   }
//  }
//  maxsub(n_points,xe,xp,mxrms,mxpsi,nali,mxzscore,mxeval,mxscore);
//
//  rms = mxrms; // double to float conversion
//  logeval = std::log(mxeval);
//
// }

////////////////////////////////////////////////////////////////////////////////
///
/// @brief identify the largest subset of CA atoms of a model that superimposes
/// "well" (under certain rms cutoff) over the experimental structure
///
/// @details
///     this function was adapted and improved by cem strauss from an
///     original template provided by angel ortiz.
///
///     Here applies a modification of Dani Fischer's heuristic algorithm
///     for finding the largest subset of residues for superimposition within
///     a threshold. The part that restraints the secondary structure
///     matching needs to be changed.
///
///     At this point, the algorithm works as follows: first, the residue
///     assignment is done on the basis of the global secondary structure
///     similarity. Then, with this assignment for residue pairs, the
///     heuristic procedure of Fisher is used.
///
///     A seed segment of 7 CA atoms was first superimposed onto the reference structure,
///     then the neighbor radius is gradually increased (up to "distance_tolerance"):to
///     include more CA atoms within the
///     range "t" to do a new superimposition, if the aligned rms is below the rms cutoff "rmstol"
///     the new alignment is accepted. This proceduare is repeated for every 7-residue seed segment
///     along the sequence until the largest subset of atoms which can yield an aligned rms below is
///     found."nsup" is the total number of CA atoms in the model, xe and xp are coordinates of native
///     and model respectively. "rms" is the final aligned rms based on the identified subset of atoms,
///     "nali" is the number of aligned atoms in the subset and "psi" is the ratio of nali/nsup; The rest
///     scores are metrics for quality of the alignment. Note that if no suitable alignment can be found
///     given the input rms cutoff, the rms of overall structure alignment is returned and "nali" is set to
///     zero. -- chu 2009/10
///
/// @param  nsup - [in/out]? -
/// @param  xe - [in/out]? -
/// @param  xp - [in/out]? -
/// @param  rms - [in/out]? -
/// @param  psi - [in/out]? -
/// @param  nali - [in/out]? -
/// @param  zscore - [in/out]? -
/// @param  evalue - [in/out]? -
/// @param  score - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
maxsub(
	int & nsup, // total number of residues to superimpose
	ObjexxFCL::FArray1A_double xe,
	ObjexxFCL::FArray1A_double xp,
	double & rms,
	double & psi,
	int & nali, // number of aligned residues
	double & zscore,
	double & evalue,
	double & score,
	double rmstol, //default = 4.0
	double distance_tolerance //default = 7.0
)
{
	using namespace ObjexxFCL;
	
	xe.dimension( 3*nsup );
	xp.dimension( 3*nsup );


	//     --- Vector alignment

	//int const maxres = { 3000 };
	//int const maxlen = { 2 * maxres };
	//int const maxfrag = { 100 };
	//double const angmax = { 60.0 };

	if ( distance_tolerance == -1.0 ) {
		distance_tolerance = ( rmstol * 7.0 ) / 4.0; // first guess, extrapolating from 4 -> 7
	}

	double distance_increment = distance_tolerance / 7.0; // maintain same number of cycles

	double znew;
	double am,as;

	FArray1D_double xp0( 3*nsup );
	FArray1D_double xe0( 3*nsup );
	FArray1D_double wmax( nsup );

	FArray1D_bool logical_w( nsup );

	numeric::model_quality::RmsData* rmsdata = RmsData::instance(); // get a pointer to the singleton class

	//------------------------------------------------------------------------------
	//     First selects the residues allowed to be superimposed
	//------------------------------------------------------------------------------
	// need to pass in xp and xe, nsup = total_res

	int nr = nsup;
	int l = 7; // length of peptide for Fisher's maxsub algorithm.

	for ( int i = 1; i <= nsup*3; ++i ) {
		xe0(i) = xe(i);
		xp0(i) = xp(i);
	}
	//------------------------------------------------------------------------------
	//     Now apply Fishers's maxsub algorithm. An heptapeptide is used.
	//     The algorithm is modified in that only pairs of residues with
	//     similar local secondary structures are allowed to be used as
	//     seed residues to superimpose.
	//------------------------------------------------------------------------------
	int smax = 0;
	for ( int i = 1; i <= nsup; ++i ) { // should this not be <=3*nsup?
		//w(i) = 0.0;
		wmax(i) = 0.0;
	}
	double rmsmax = 1000.0;


	for ( int i = 1; i <= nsup-l+1; ++i ) {

		rmsdata->clear_rms();
		// if ( matched(i)} ) {  // this line is Angel's variation on Danni Fischer's method.
		// only seed at points with matched SS.

		// build up a seed segment of length l
		int lmax = 0;
		for ( int j = 1; j <= nsup; ++j ) {
			if ( j >= i && j <= i+l-1 ) { // could do this without if statement
				rmsdata->add_rms(j,xp0,xe0);
				logical_w(j) = true; // w(j) = 1.0
				++lmax;
			} else {
				logical_w(j) = false; // w(j) = 0.0
			}
		}

		// find initial alignment using seed
		// rmsfitca3 rotates all of the residues but the rotation only aligns the
		// residues pushed into add_rms.
		rmsfitca3(nsup,xp0,xp,xe0,xe,rms);
		// chu 2009/10
		// accept the inital seed alignment in case the first distance incremental
		// step fails to find a better alignment
		if ( (lmax > smax) && (rms <= rmstol) ) {
			smax = lmax;
			for ( int n = 1; n <= nsup; ++n ) {
				if ( logical_w(n) ) {
					wmax(n) = 1.0;
				} else {
					wmax(n) = 0.0;
				}
				// wmax(n) = w(n);
			}
			rmsmax = rms;
		} else if ( (lmax == smax) && (rms < rmsmax) ) {
			smax = lmax;
			for ( int n = 1; n <= nsup; ++n ) {
				if ( logical_w(n) ) {
					wmax(n) = 1.0;
				} else {
					wmax(n) = 0.0;
				}
				//wmax(n) = w(n);
			}
			rmsmax = rms;
		}
		// next we iterate the following algorithm
		//  1) using current alignment find all atoms within a threshold t of being superimposed
		//  2) add these close ones to set to be aligned
		//  3) aling using the current set, then reorient all atoms.
		//  4) increment t by a small amount.
		//  5) repeat this until theshold = 7 angstroms.
		double t = 0.0;
		int const last = lmax;
		double min_d2 = square ( distance_tolerance + 1.0 );
		//  std::cout << min_d2 << std::endl;
		while ( t < distance_tolerance ) {  // for ( m = 1; m <= 7; ++m ) {
			t += distance_increment;
			//   std::cout << "maxsub: t = "<< t << " distance_tolerance = " << distance_tolerance << "\n";
			// increment threshold by one angstrom. (note this is ties to int(min_d))
			double t2 = t*t; // t squared


			// t = float(m);  // *7.0/float(7);  // huh? must be a relic???
			for ( int n = 1; n <= nsup; ++n ) {
				if ( !logical_w(n) ) { // if ( w(n) == 0.0} ) {
					int k = 3*(n-1);
					double const xpek1 = xp(k+1) - xe(k+1);
					double const xpek2 = xp(k+2) - xe(k+2);
					double const xpek3 = xp(k+3) - xe(k+3);
					double const d2 =
						( xpek1 * xpek1 ) + ( xpek2 * xpek2 ) + ( xpek3 * xpek3 );
					// squared distancne
					if ( d2 <= t2 ) { // is this atom within threshold?
						rmsdata->add_rms(n,xp0,xe0); // if so, add to list
						logical_w(n) = true; //w(n) = 1.0  // set membership flag
						++lmax; // keep a count of members
					} else {
						// if not below threshold, then find the closest atom
						if ( d2 <= min_d2 ) min_d2 = d2; //min_d = min(d,min_d)
					}
				}
			}
			// check if we added any residues on this iteration.
			// if not then 1) we dont need to refit the calphas 2) we can advance threshold level
			if ( lmax != last ) {
				rmsfitca3(nsup,xp0,xp,xe0,xe,rms);
			} else {
				t = std::sqrt( min_d2 ); //advance the threshold
				//    t = static_cast< int >(std::sqrt(min_d2)); // advance the threshold

				//   std::cout << i << " skipping " << t << ' ' <<
				// static_cast< int >(min_d2) << ' ' << lmax << ' ' << min_d2 << " new t = " << t<< std::endl;

			}


			// huh? logic here is confusing.
			// chu 2009/10
			// after each distance incremental step to include more atoms into the superimpostion,
			// check whether the new rms is still under the rmstol cutoff.
			// If more atoms can be aligned under the rmstol cutoff, accept this new alignment;
			// If in a different set (of the same size) of atoms can be aligned yielding lower rms, aceept this new alignment
			// Otherwise, ignore this new alignment and keep searching until done.
			if ( (lmax > smax) && (rms <= rmstol) ) {
				smax = lmax;
				for ( int n = 1; n <= nsup; ++n ) {
					if ( logical_w(n) ) {
						wmax(n) = 1.0;
					} else {
						wmax(n) = 0.0;
					}
					// wmax(n) = w(n);
				}
				rmsmax = rms;
			} else if ( (lmax == smax) && (rms < rmsmax) ) {
				smax = lmax;
				for ( int n = 1; n <= nsup; ++n ) {
					if ( logical_w(n) ) {
						wmax(n) = 1.0;
					} else {
						wmax(n) = 0.0;
					}
					//wmax(n) = w(n);
				}
				rmsmax = rms;
			}
		} // while( t< distance_tolerance) -- chu 2009/10 bugfix
	}

	//------------------------------------------------------------------------------
	//     --- Confirm final superimposition.
	//     --- first, compile regions without indels. Report rms
	//------------------------------------------------------------------------------
	// std::cout << "maxsub: next step " << std::endl;
	if ( smax > 1 ) {
		rmsfitca2(nr,xp,xe,wmax,smax,rms);
		// side effect sets xpc,zpc etc... via common
	} else {
		// if smax is less than 2 then basically we failed to find an alignement
		//  to make the best of a bad situation we simply revert to aligning all of the residue
		for ( int i = 1; i <= nr; ++i ) {
			wmax(i) = 1.0;
		}
		rmsfitca2(nr,xp,xe,wmax,nr,rms);
		// side effect sets xpc,zpc etc... via common
	}

	// std::cout << std::endl << "RMS = " << F( 7, 3, rms ) <<
	//  " SMAX = " << I( 5, smax ) << std::endl << std::endl;
	psi = ( static_cast< double >( smax ) / nsup ) * 100.0;

	//------------------------------------------------------------------------------
	//     --- Transform PSI to Levitt & Gerstein Sstr
	//     --- NEED TO BE REMOVED
	//------------------------------------------------------------------------------

	//     compute score without gaps

	score = 0.0;
	nali  = 0;

	for ( int i = 1; i <= nsup; ++i ) {
		if ( wmax(i) == 1.0 ) {
			double d = 0.0;
			int k = 3*(i-1);
			for ( int j = 1; j <= 3; ++j ) {
				double const xpekj = xp(k+j) - xe(k+j);
				d += xpekj * xpekj;
			}
			d = std::sqrt(d) / rmstol;
			score += 1.0 / ( 1.0 + ( d * d ) );
			++nali;
		}
	}
	//chu 2009/10 if no alignment is found, return nali as zero
	if ( smax <=1 ) nali=smax;
	//------------------------------------------------------------------------------
	//     --- All done. Report match statistics
	//------------------------------------------------------------------------------

	//------------------------------------------------------------------------------
	//     --- Report probabilities for random matches.
	//     --- These are preliminary values for the average and standard deviation
	//     --- of PSI as a function of NORM. These values are:
	//     --- m(L) = 759.31 * L**(-0.7545)
	//     --- s(L) = 393.32 * L**(-0.9009)
	//------------------------------------------------------------------------------

	//     am = 759.31 * std::pow( norm, -0.7545 );
	//     as = 393.32 * std::pow( norm, -0.9009 );
	//     am = 695.41 * std::pow( norm, -0.7278 );
	//     as = 340.00 * std::pow( norm, -0.9045 );
	//     EV fitting, using N > 70
	am = 747.29 * std::pow( static_cast< double >( nr ), -0.7971 );
	as = 124.99 * std::pow( static_cast< double >( nr ), -0.6882 );
	zscore = (psi-am)/as;
	//     this is the gaussian approach. Actually, extreme-value is more adequate
	evalue = 0.5 * erfcc(zscore/std::sqrt(2.0));
	//     here it is the EV approach
	znew   = 0.730*((1.2825755*zscore)+0.5772);
	evalue = 1.0-std::exp(-std::exp(-znew));
	//     due to numerical errors, the e-value is cutoff at 2.650E-14
	if ( evalue < 2.650E-14 ) evalue = 2.650E-14;

}

//////////////////////////////////////////////////////////////////////////////
///
/// @brief
///
/// @details
///
/// @param  x - [in/out]? -
///
/// @return
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///  (C) Copr. 1986-92 Numerical Recipes Software
///
/// @author
///
////////////////////////////////////////////////////////////////////////////////
double
erfcc( double x )
{

	double t,z;
	z = std::abs(x);
	t = 1.0/(1.0+0.5*z);
	double erfcc = t*std::exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+
		t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+t*(1.48851587+
		t*(-.82215223+t*.17087277)))))))));
	if ( x < 0.0 ) erfcc = 2.0 - erfcc;
	return erfcc;
}

////////////////////////////////////////////////////////////////////////////////
///
/// @brief
//    Calculate the center of geometry for the selected atoms ---
///
/// @details
///
/// @param  C - [in/out]? -
/// @param  WT - [in/out]? -
/// @param  NAT - [in/out]? -
/// @param  XC - [in/out]? -
/// @param  YC - [in/out]? -
/// @param  ZC - [in/out]? -
///
/// @global_read
///
/// @global_write
///
/// @remarks
///
/// @references
///
/// @author
///
/////////////////////////////////////////////////////////////////////////////////
void
COMAS(
	ObjexxFCL::FArray1A_double C,  // coordinates for each atom (must be of length NAT)
	ObjexxFCL::FArray1A_double WT, // weights on each atom (must be of length NAT)
	int NAT, // number of atoms provided
	double & XC, // x-coordinate for center of mass
	double & YC, // y-coordinate for center of mass
	double & ZC  // z-coordinate for center of mass
)
{
	using namespace ObjexxFCL;
	
	C.dimension( star );
	WT.dimension( star );

	//     local
	static double const ZERO = { 0.0 };

	double SUMX = ZERO;
	double SUMY = ZERO;
	double SUMZ = ZERO;
	double SUM  = ZERO;
	int i3 = 0;

	for ( int i = 1; i <= NAT; ++i ) {
		double const WT_i = WT(i);

		SUMX += C(i3+1) * WT_i;
		SUMY += C(i3+2) * WT_i;
		SUMZ += C(i3+3) * WT_i;
		SUM += WT_i;
		i3 += 3;
	}

	SUM = 1.0/SUM;

	XC = SUMX*SUM;

	YC = SUMY*SUM;

	ZC = SUMZ*SUM;


	i3 = 0;

	for ( int i = 1; i <= NAT; ++i ) {
		C(i3+1) -= XC;
		C(i3+2) -= YC;
		C(i3+3) -= ZC;
		i3 += 3;
	}

	// std::cout << "CENTER OF MASS:" << space( 19 ) <<
	//  F( 10, 4, XC ) << F( 10, 4, YC ) << F( 10, 4, ZC ) << std::endl;
} // COMAS

} // namespace model_quality
} // namespace numeric
