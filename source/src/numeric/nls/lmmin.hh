//
// Project:  LevenbergMarquardtLeastSquaresFitting
//
// File:     lmmin.h
//
// Contents: Public interface to the Levenberg-Marquardt core implementation.
//
// Author:   Joachim Wuttke 2004-2010
//
// Homepage: www.messen-und-deuten.de/lmfit


#ifndef INCLUDED_numeric_nls_LMMIN_HH
#define INCLUDED_numeric_nls_LMMIN_HH

namespace numeric {
namespace nls {

//#ifdef __cplusplus
//extern "C" {
//#endif


// Compact high-level interface.
// Collection of control (input) parameters.
struct lm_control_struct
{
	double ftol;      // relative error desired in the sum of squares.
	double xtol;      // relative error between last two approximations.
	double gtol;      // orthogonality desired between fvec and its derivs.
	double epsilon;   // step used to calculate the jacobian.
	double stepbound; // initial bound to steps in the outer loop.
	int maxcall;      // maximum number of iterations.
	int scale_diag;   // UNDOCUMENTED, TESTWISE automatical diag rescaling?
	int printflags;   // OR'ed to produce more noise
};

// Collection of status (output) parameters.
struct lm_status_struct {
	double fnorm;     // norm of the residue vector fvec.
	int nfev;       // actual number of iterations.
	int info;       // status of minimization.
};

// Recommended control parameter settings.
//const lm_control_struct lm_control_double;
//const lm_control_struct lm_control_float;

// Standard monitoring routine.
void lm_printout_std( int n_par, const double *par, int m_dat,
	const void *data, const double *fvec,
	int printflags, int iflag, int iter, int nfev);

// Refined calculation of Eucledian norm, typically used in printout routine.
double lm_enorm( int, const double * );

// The actual minimization.
void lmmin( int n_par, double *par, int m_dat, const void *data,
	void (*evaluate) (const double *par, int m_dat, const void *data,
	double *fvec, int *info),
	lm_status_struct *status,
	void (*printout) (int n_par, const double *par, int m_dat,
	const void *data, const double *fvec,
	int printflags, int iflag, int iter, int nfev) );


// Legacy low-level interface.
// Alternative to lm_minimize, allowing full control, and read-out
// of auxiliary arrays. For usage, see implementation of lm_minimize.
void lm_lmdif( int m, int n, double *x, double *fvec, double ftol,
	double xtol, double gtol, int maxfev, double epsfcn,
	double *diag, int mode, double factor, int *info, int *nfev,
	double *fjac, int *ipvt, double *qtf, double *wa1,
	double *wa2, double *wa3, double *wa4,
	void (*evaluate) (const double *par, int m_dat, const void *data,
	double *fvec, int *info),
	void (*printout) (int n_par, const double *par, int m_dat,
	const void *data, const double *fvec,
	int printflags, int iflag, int iter, int nfev),
	int printflags, const void *data );

//const char *lm_infmsg[];
//const char *lm_shortmsg[];

//#ifdef __cplusplus
//}
//#endif

}
}

#endif // INCLUDED_numeric_nls_LMMIN_HH
