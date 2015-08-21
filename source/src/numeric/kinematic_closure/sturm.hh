// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   sturm.hh
/// @brief  implements sturm chain solver
/// @author Daniel J. Mandell

#ifndef INCLUDED_numeric_kinematic_closure_sturm_hh
#define INCLUDED_numeric_kinematic_closure_sturm_hh


// Rosetta Headers
#include <numeric/types.hh>

// C++ headers
#include <cmath>

// Utility headers

// option key includes


#include <utility/vector1_bool.hh>


// Constants
#define MAX_ORDER 16
#define MAXPOW 32
#define SMALL_ENOUGH 1.0e-18 // a coefficient smaller than SMALL_ENOUGH is considered to be zero (0.0)
#define MAX_SOLN 16 // maximum number of roots

namespace numeric {
namespace kinematic_closure {

/*
* structure type for representing a polynomial
*/
typedef   struct p {
	int ord;
	double coef[MAX_ORDER+1];
} poly;

extern double RELERROR;
extern int MAXIT, MAX_ITER_SECANT;

/* set termination criteria for polynomial solver */
inline
void initialize_sturm(double *tol_secant, int *max_iter_sturm, int *max_iter_secant)
{
	RELERROR = *tol_secant;
	MAXIT = *max_iter_sturm;
	MAX_ITER_SECANT = *max_iter_secant;
}

inline
double hyper_tan(double a, double x)
{
	double exp_x1, exp_x2, ax;

	ax = a*x;
	if ( ax > 100.0 ) {
		return(1.0);
	} else if ( ax < -100.0 ) {
		return(-1.0);
	} else {
		exp_x1 = exp(ax);
		exp_x2 = exp(-ax);
		return (exp_x1 - exp_x2)/(exp_x1 + exp_x2);
	}
}

/*
* modp
*
* calculates the modulus of u(x) / v(x) leaving it in r, it
*  returns 0 if r(x) is a constant.
*  note: this function assumes the leading coefficient of v
* is 1 or -1
*/

inline
int modp(poly* u, poly* v, poly* r)
{
	int  k, j;
	double *nr, *end, *uc;

	nr = r->coef;
	end = &u->coef[u->ord];

	uc = u->coef;
	while ( uc <= end ) {
		*nr++ = *uc++;
	}

	if ( v->coef[v->ord] < 0.0 ) {
		for ( k = u->ord - v->ord - 1; k >= 0; k -= 2 ) {
			r->coef[k] = -r->coef[k];
		}
		for ( k = u->ord - v->ord; k >= 0; k-- ) {
			for ( j = v->ord + k - 1; j >= k; j-- ) {
				r->coef[j] = -r->coef[j] - r->coef[v->ord + k] * v->coef[j - k];
			}
		}
	} else {
		for ( k = u->ord - v->ord; k >= 0; k-- ) {
			for ( j = v->ord + k - 1; j >= k; j-- ) {
				r->coef[j] -= r->coef[v->ord + k] * v->coef[j - k];
			}
		}
	}

	k = v->ord - 1;
	while ( k >= 0 && fabs(r->coef[k]) < SMALL_ENOUGH ) {
		r->coef[k] = 0.0;
		k--;
	}

	r->ord = (k < 0) ? 0 : k;
	return(r->ord);
}

/*
* buildsturm
*
* build up a sturm sequence for a polynomial in smat, returning
* the number of polynomials in the sequence
*/

inline
int buildsturm(int ord, poly* sseq)
{
	int  i;
	double f, *fp, *fc;
	poly *sp;

	sseq[0].ord = ord;
	sseq[1].ord = ord - 1;

	/*
	* calculate the derivative and normalise the leading
	* coefficient.
	*/
	f = fabs(sseq[0].coef[ord]*ord);

	fp = sseq[1].coef;
	fc = sseq[0].coef + 1;

	for ( i = 1; i <= ord; i++ ) {
		*fp++ = *fc++ * i / f;
	}


	/*
	* construct the rest of the Sturm sequence
	*/
	for ( sp = sseq + 2; modp(sp - 2, sp - 1, sp); sp++ ) {

		/*
		* reverse the sign and normalise
		*/

		f = -fabs(sp->coef[sp->ord]);
		for ( fp = &sp->coef[sp->ord]; fp >= sp->coef; fp-- ) {
			*fp /= f;
		}
	}

	sp->coef[0] = -sp->coef[0]; /* reverse the sign */

	return(sp - sseq);
}

/*
* numroots
*
* return the number of distinct real roots of the polynomial
* described in sseq.
*/

inline
int numroots(int np, poly* sseq, int* atneg, int* atpos)
{
	int  atposinf, atneginf;
	poly *s;
	double f, lf;

	atposinf = atneginf = 0;


	/*
	* changes at positive infinity
	*/
	lf = sseq[0].coef[sseq[0].ord];

	for ( s = sseq + 1; s <= sseq + np; s++ ) {
		f = s->coef[s->ord];
		if ( lf == 0.0 || lf * f < 0 ) {
			atposinf++;
		}
		lf = f;
	}

	//std::cout << "sturm.cc::numroots lf: " << lf << std::endl;

	/*
	* changes at negative infinity
	*/
	if ( sseq[0].ord & 1 ) {
		lf = -sseq[0].coef[sseq[0].ord];
	} else {
		lf = sseq[0].coef[sseq[0].ord];
	}

	for ( s = sseq + 1; s <= sseq + np; s++ ) {
		if ( s->ord & 1 ) {
			f = -s->coef[s->ord];
		} else {
			f = s->coef[s->ord];
		}
		if ( lf == 0.0 || lf * f < 0 ) {
			atneginf++;
		}
		lf = f;
	}

	//std::cout << "sturm.cc::numroots at neg inf: " << atneginf << std::endl;
	//std::cout << "sturm.cc::numroots at pos inf: " << atposinf << std::endl;

	*atneg = atneginf;
	*atpos = atposinf;
	return(atneginf - atposinf);
}

/*
* evalpoly
*
* evaluate polynomial defined in coef returning its value.
*/

inline
double evalpoly (int ord, double* coef, double x)
{
	double *fp, f;

	fp = &coef[ord];
	f = *fp;

	for ( fp--; fp >= coef; fp-- ) {
		f = x * f + *fp;
	}

	return(f);
}

/*
* numchanges
*
* return the number of sign changes in the Sturm sequence in
* sseq at the value a.
*/

inline
int numchanges(int np, poly* sseq, double a)
{
	int  changes;
	double f, lf;
	poly *s;

	changes = 0;

	lf = evalpoly(sseq[0].ord, sseq[0].coef, a);

	for ( s = sseq + 1; s <= sseq + np; s++ ) {
		f = evalpoly(s->ord, s->coef, a);
		if ( lf == 0.0 || lf * f < 0 ) {
			changes++;
		}
		lf = f;
	}

	return(changes);
}

/*
* modrf
*
* uses the modified regula-falsi method to evaluate the root
* in interval [a,b] of the polynomial described in coef. The
* root is returned is returned in *val. The routine returns zero
* if it can't converge.
*/

inline
int modrf(int ord, double *coef, double a, double b, double *val)
{
	int its;
	double fa, fb, x, fx, lfx;
	double *fp, *scoef, *ecoef;
	fx=0; // avoid uninitialized variable warning

	scoef = coef;
	ecoef = &coef[ord];

	fb = fa = *ecoef;
	for ( fp = ecoef - 1; fp >= scoef; fp-- ) {
		fa = a * fa + *fp;
		fb = b * fb + *fp;
	}

	/*
	* if there is no sign difference the method won't work
	*/
	if ( fa * fb > 0.0 ) {
		return(0);
	}
	/*  commented out to avoid duplicate solutions when the bounds are close to the roots

	if (fabs(fa) < RELERROR)
	{
	*val = a;
	return(1);
	}

	if (fabs(fb) < RELERROR)
	{
	*val = b;
	return(1);
	}
	*/

	lfx = fa;

	for ( its = 0; its < MAX_ITER_SECANT; its++ ) {
		x = (fb * a - fa * b) / (fb - fa);

		// constrain that x stays in the bounds
		if ( x < a || x > b ) {
			x = 0.5 * (a+b);
		}
		fx = *ecoef;
		for ( fp = ecoef - 1; fp >= scoef; fp-- ) {
			fx = x * fx + *fp;
		}

		if ( fabs(x) > RELERROR ) {
			if ( fabs(fx / x) < RELERROR ) {
				*val = x;
				return(1);
			}
		} else if ( fabs(fx) < RELERROR ) {
			*val = x;
			return(1);
		}

		if ( (fa * fx) < 0 ) {
			b = x;
			fb = fx;
			if ( (lfx * fx) > 0 ) {
				fa /= 2;
			}
		} else {
			a = x;
			fa = fx;
			if ( (lfx * fx) > 0 ) {
				fb /= 2;
			}
		}
		lfx = fx;
	}

	//std::cout << "sturm.cc::modrf overflow " << a << " " << b << " " << fx << std::endl;
	return(0);
}


/*
* sbisect
*
* uses a bisection based on the sturm sequence for the polynomial
* described in sseq to isolate intervals in which roots occur,
* the roots are returned in the roots array in order of magnitude.
*
* //DJM 6-20-07: added index to recurse through utility::vector1. Always
* //call sbisect with index=1. Will be increased during recursion.
*/
//void sbisect(int np, poly* sseq, double min, double max, int atmin, int atmax, utility::vector1<double>& roots, int index)

inline
void sbisect(int np, poly* sseq, double min, double max, int atmin, int atmax, double* roots)
{
	double mid=0.0;
	int  n1 = 0, n2 = 0, its, atmid, nroot;

	if ( (nroot = atmin - atmax) == 1 ) {

		/*
		* first try a less expensive technique.
		*/
		//if (modrf(sseq->ord, sseq->coef, min, max, &roots[index])) {
		if ( modrf(sseq->ord, sseq->coef, min, max, &roots[0]) ) {
			return;
		}

		/*
		* if we get here we have to evaluate the root the hard
		* way by using the Sturm sequence.
		*/
		for ( its = 0; its < MAXIT; its++ ) {
			mid = (min + max) / 2;
			atmid = numchanges(np, sseq, mid);

			if ( fabs(mid) > RELERROR ) {
				if ( fabs((max - min) / mid) < RELERROR ) {
					//roots[index] =  static_cast<Real> ( mid );
					roots[0] = mid;
					//std::cout << "sturm.cc::sbisect mid, min, max: " << mid << ", " << min << ", " << max << std::endl;
					return;
				}
			} else if ( fabs(max - min) < RELERROR ) {
				//roots[index] = static_cast<Real> ( mid );
				roots[0] = mid;
				//std::cout << "sturm.cc::sbisect mid, min, max: " << mid << ", " << min << ", " << max << std::endl;
				return;
			}

			if ( (atmin - atmid) == 0 ) {
				min = mid;
			} else {
				max = mid;
			}
		}
		if ( its == MAXIT ) {
			//roots[index] = static_cast<Real> ( mid );
			roots[0] = mid;
		}
		return;
	}

	/*
	* more than one root in the interval, we have to bisect...
	*/

	for ( its = 0; its < MAXIT; its++ ) {
		mid = (min + max) / 2;
		atmid = numchanges(np, sseq, mid);
		n1 = atmin - atmid;
		n2 = atmid - atmax;

		// DJM: this rarely occuring "artifact" causes crashing bugs so we take what we've got rather than proceed
		if ( n1 < 0 ) {
			//for (n1 = atmax; n1 < atmin; n1++) {
			// roots[index+(n1 - atmax)] = static_cast<Real> ( mid );
			//}
			//std::cout << "atmin was less than atmid (" << atmin << " vs " << atmid << "). Exiting sturm prematurely." <<
			//std::cout << "atmin was less than atmid (" << atmin << " vs " << atmid << "). But not exiting prematurely." <<
			// std::endl;
			//break;
		}

		if ( n1 != 0 && n2 != 0 ) {
			//sbisect(np, sseq, min, mid, atmin, atmid, roots, index);
			//sbisect(np, sseq, mid, max, atmid, atmax, roots, index+n1);
			sbisect(np, sseq, min, mid, atmin, atmid, roots);
			sbisect(np, sseq, mid, max, atmid, atmax, &roots[n1]);
			break;
		}

		if ( n1 == 0 ) {
			min = mid;
		} else {
			max = mid;
		}
	}

	if ( its == MAXIT ) {
		//std::cout << "Maximum bisection iterations reached in sturm. Returning mid value for remaining solutions." << std::endl;
		for ( n1 = atmax; n1 < atmin; n1++ ) {
			//roots[index+(n1 - atmax)] = static_cast<Real> ( mid );
			roots[n1 - atmax] = mid;
		}
	}
}

inline
void solve_sturm(const int& p_order, int& n_root, const utility::vector1<double>& poly_coeffs, utility::vector1<double>& roots)
{
	poly sseq[MAX_ORDER*2];
	double min, max;
	double roots_array[MAX_ORDER];
	//double* roots_array;
	//roots_array = new double[MAX_ORDER];
	int order, i, nroots, nchanges, np, atmin, atmax;

	order = p_order;
	// DJM: !!!here we are going from base1 vector to base0 array!!!
	for ( i = order; i >= 0; i-- ) {
		sseq[0].coef[i] = poly_coeffs[i+1];
	}

	/*
	* build the Sturm sequence
	*/
	np = buildsturm(order, sseq);

	/* // DJM: Debug
	if (basic::options::option[basic::options::OptionKeys::run::run_level] >= basic::options::verbose ) {
	//std::cout << "sturm.cc::solve_sturm Sturm sequence for: " << std::endl;
	for (int i = order; i >= 0; i--) {
	//std::cout << sseq[0].coef[i] << " ";
	}
	//std::cout << std::endl << std::endl;
	for (int i = 0; i <= np; i++) {
	for (int j = sseq[i].ord; j >= 0; j--) {
	//std::cout << sseq[i].coef[j] << " ";
	}
	//std::cout << std::endl;
	}
	}
	*/

	/*
	* get the number of real roots
	*/

	nroots = numroots(np, sseq, &atmin, &atmax);

	if ( (nroots == 0) || !( (0 <= nroots) && (nroots <= MAX_ORDER)) ) { // make sure we have some sensible roots
		n_root = 0;
		return;
	}

	roots.resize(nroots);

	//std::cout << "sturm.cc::solve_sturm: Number of real roots: " << nroots << std::endl;

	/*
	* calculate the bracket that the roots live in
	*/
	min = -1.0;
	nchanges = numchanges(np, sseq, min);

	for ( i = 0; nchanges != atmin && i != MAXPOW; i++ ) {
		min *= 10.0;
		nchanges = numchanges(np, sseq, min);
	}

	if ( nchanges != atmin ) {
		//std::cout << "sturm.cc::solve: unable to bracket all negative roots" << std::endl;
		atmin = nchanges;
	}

	max = 1.0;
	nchanges = numchanges(np, sseq, max);
	for ( i = 0; nchanges != atmax && i != MAXPOW; i++ ) {
		max *= 10.0;
		nchanges = numchanges(np, sseq, max);
	}

	if ( nchanges != atmax ) {
		atmax = nchanges;
	}

	nroots = atmin - atmax;

	//std::cout << "sturm.cc::solve_sturm atmin: " << atmin << std::endl;
	//std::cout << "sturm.cc::solve_sturm atmax: " << atmax << std::endl;
	//std::cout << "sturm.cc::solve_sturm nroots: " << nroots << std::endl;

	/*
	* perform the bisection.
	*/

	//sbisect(np, sseq, min, max, atmin, atmax, roots, 1); // DJM: Vector1 is base 1
	sbisect(np, sseq, min, max, atmin, atmax, roots_array);

	// DJM: copy the roots array into the roots vector
	for ( numeric::Size i=1; i<= roots.size(); i++ ) {
		roots[i]=roots_array[i-1];
	}
	//// DJM: now we can delete the roots array -- not if statically allocated
	//delete [] roots_array;
	//roots_array = NULL;

	n_root = nroots;

	/*
	* write out the roots...
	*/
	/* // DJM: debug
	if (basic::options::option[basic::options::OptionKeys::run::run_level] >= basic::options::verbose ) {
	if (nroots == 1) {
	//std::cout << std::endl << "sturm.cc::solve_sturm 1 distinct real root at x = " << roots[1] << std::endl;
	}
	else {
	//std::cout << nroots << " distinct real roots for x: " << std::endl;
	for (i = 1; i <= nroots; i++) {
	//std::cout << roots[i] << std::endl;
	}
	}
	}
	*/
	return;
}

} // end namespace kinematic_closure
} // end namespace numeric

#endif //INCLUDED_numeric_kinematic_closure_sturm_hh
