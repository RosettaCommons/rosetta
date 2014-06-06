/*************************************************************************
Copyright (c) Sergey Bochkanov (ALGLIB project).

>>> SOURCE LICENSE >>>
This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation (www.fsf.org); either version 2 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License is available at
http://www.fsf.org/licensing/licenses
>>> END OF LICENSE >>>
*************************************************************************/
#include "stdafx.h"
#include "alglibinternal.h"

// disable some irrelevant warnings
#if (AE_COMPILER==AE_MSVC)
#pragma warning(disable:4100)
#pragma warning(disable:4127)
#pragma warning(disable:4702)
#pragma warning(disable:4996)
#endif
using namespace std;

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS IMPLEMENTATION OF C++ INTERFACE
//
/////////////////////////////////////////////////////////////////////////
namespace alglib
{


}

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS IMPLEMENTATION OF COMPUTATIONAL CORE
//
/////////////////////////////////////////////////////////////////////////
namespace alglib_impl
{




static void tsort_tagsortfastirec(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     /* Real    */ ae_vector* bufa,
     /* Integer */ ae_vector* bufb,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state);
static void tsort_tagsortfastrrec(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* bufa,
     /* Real    */ ae_vector* bufb,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state);
static void tsort_tagsortfastrec(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* bufa,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state);


















static void hsschur_internalauxschur(ae_bool wantt,
     ae_bool wantz,
     ae_int_t n,
     ae_int_t ilo,
     ae_int_t ihi,
     /* Real    */ ae_matrix* h,
     /* Real    */ ae_vector* wr,
     /* Real    */ ae_vector* wi,
     ae_int_t iloz,
     ae_int_t ihiz,
     /* Real    */ ae_matrix* z,
     /* Real    */ ae_vector* work,
     /* Real    */ ae_vector* workv3,
     /* Real    */ ae_vector* workc1,
     /* Real    */ ae_vector* works1,
     ae_int_t* info,
     ae_state *_state);
static void hsschur_aux2x2schur(double* a,
     double* b,
     double* c,
     double* d,
     double* rt1r,
     double* rt1i,
     double* rt2r,
     double* rt2i,
     double* cs,
     double* sn,
     ae_state *_state);
static double hsschur_extschursign(double a, double b, ae_state *_state);
static ae_int_t hsschur_extschursigntoone(double b, ae_state *_state);




static ae_bool safesolve_cbasicsolveandupdate(ae_complex alpha,
     ae_complex beta,
     double lnmax,
     double bnorm,
     double maxgrowth,
     double* xnorm,
     ae_complex* x,
     ae_state *_state);


static void xblas_xsum(/* Real    */ ae_vector* w,
     double mx,
     ae_int_t n,
     double* r,
     double* rerr,
     ae_state *_state);
static double xblas_xfastpow(double r, ae_int_t n, ae_state *_state);


static double linmin_ftol = 0.001;
static double linmin_xtol = 100*ae_machineepsilon;
static ae_int_t linmin_maxfev = 20;
static double linmin_stpmin = 1.0E-50;
static double linmin_defstpmax = 1.0E+50;
static double linmin_armijofactor = 1.3;
static void linmin_mcstep(double* stx,
     double* fx,
     double* dx,
     double* sty,
     double* fy,
     double* dy,
     double* stp,
     double fp,
     double dp,
     ae_bool* brackt,
     double stmin,
     double stmax,
     ae_int_t* info,
     ae_state *_state);


static ae_int_t ftbase_ftbaseplanentrysize = 8;
static ae_int_t ftbase_ftbasecffttask = 0;
static ae_int_t ftbase_ftbaserfhttask = 1;
static ae_int_t ftbase_ftbaserffttask = 2;
static ae_int_t ftbase_fftcooleytukeyplan = 0;
static ae_int_t ftbase_fftbluesteinplan = 1;
static ae_int_t ftbase_fftcodeletplan = 2;
static ae_int_t ftbase_fhtcooleytukeyplan = 3;
static ae_int_t ftbase_fhtcodeletplan = 4;
static ae_int_t ftbase_fftrealcooleytukeyplan = 5;
static ae_int_t ftbase_fftemptyplan = 6;
static ae_int_t ftbase_fhtn2plan = 999;
static ae_int_t ftbase_ftbaseupdatetw = 4;
static ae_int_t ftbase_ftbasecodeletrecommended = 5;
static double ftbase_ftbaseinefficiencyfactor = 1.3;
static ae_int_t ftbase_ftbasemaxsmoothfactor = 5;
static void ftbase_ftbasegenerateplanrec(ae_int_t n,
     ae_int_t tasktype,
     ftplan* plan,
     ae_int_t* plansize,
     ae_int_t* precomputedsize,
     ae_int_t* planarraysize,
     ae_int_t* tmpmemsize,
     ae_int_t* stackmemsize,
     ae_int_t stackptr,
     ae_state *_state);
static void ftbase_ftbaseprecomputeplanrec(ftplan* plan,
     ae_int_t entryoffset,
     ae_int_t stackptr,
     ae_state *_state);
static void ftbase_ffttwcalc(/* Real    */ ae_vector* a,
     ae_int_t aoffset,
     ae_int_t n1,
     ae_int_t n2,
     ae_state *_state);
static void ftbase_internalcomplexlintranspose(/* Real    */ ae_vector* a,
     ae_int_t m,
     ae_int_t n,
     ae_int_t astart,
     /* Real    */ ae_vector* buf,
     ae_state *_state);
static void ftbase_internalreallintranspose(/* Real    */ ae_vector* a,
     ae_int_t m,
     ae_int_t n,
     ae_int_t astart,
     /* Real    */ ae_vector* buf,
     ae_state *_state);
static void ftbase_ffticltrec(/* Real    */ ae_vector* a,
     ae_int_t astart,
     ae_int_t astride,
     /* Real    */ ae_vector* b,
     ae_int_t bstart,
     ae_int_t bstride,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state);
static void ftbase_fftirltrec(/* Real    */ ae_vector* a,
     ae_int_t astart,
     ae_int_t astride,
     /* Real    */ ae_vector* b,
     ae_int_t bstart,
     ae_int_t bstride,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state);
static void ftbase_ftbasefindsmoothrec(ae_int_t n,
     ae_int_t seed,
     ae_int_t leastfactor,
     ae_int_t* best,
     ae_state *_state);
static void ftbase_fftarrayresize(/* Integer */ ae_vector* a,
     ae_int_t* asize,
     ae_int_t newasize,
     ae_state *_state);
static void ftbase_reffht(/* Real    */ ae_vector* a,
     ae_int_t n,
     ae_int_t offs,
     ae_state *_state);









ae_int_t getrdfserializationcode(ae_state* ) //*_state)
{
    ae_int_t result;


    result = 1;
    return result;
}


ae_int_t getkdtreeserializationcode(ae_state* ) //*_state)
{
    ae_int_t result;


    result = 2;
    return result;
}


ae_int_t getmlpserializationcode(ae_state* ) //*_state)
{
    ae_int_t result;


    result = 3;
    return result;
}


ae_int_t getmlpeserializationcode(ae_state* ) //*_state)
{
    ae_int_t result;


    result = 4;
    return result;
}


ae_int_t getrbfserializationcode(ae_state* ) //*_state)
{
    ae_int_t result;


    result = 5;
    return result;
}




/*************************************************************************
This function compares two numbers for approximate equality, with tolerance
to errors as large as max(|a|,|b|)*tol.


  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool approxequalrel(double a, double b, double tol, ae_state *_state)
{
    ae_bool result;


    result = ae_fp_less_eq(ae_fabs(a-b, _state),ae_maxreal(ae_fabs(a, _state), ae_fabs(b, _state), _state)*tol);
    return result;
}


/*************************************************************************
This  function  generates  1-dimensional  general  interpolation task with
moderate Lipshitz constant (close to 1.0)

If N=1 then suborutine generates only one point at the middle of [A,B]

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void taskgenint1d(double a,
     double b,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t i;
    double h;

    ae_vector_clear(x);
    ae_vector_clear(y);

    ae_assert(n>=1, "TaskGenInterpolationEqdist1D: N<1!", _state);
    ae_vector_set_length(x, n, _state);
    ae_vector_set_length(y, n, _state);
    if( n>1 )
    {
        x->ptr.p_double[0] = a;
        y->ptr.p_double[0] = 2*ae_randomreal(_state)-1;
        h = (b-a)/(n-1);
        for(i=1; i<=n-1; i++)
        {
            if( i!=n-1 )
            {
                x->ptr.p_double[i] = a+(i+0.2*(2*ae_randomreal(_state)-1))*h;
            }
            else
            {
                x->ptr.p_double[i] = b;
            }
            y->ptr.p_double[i] = y->ptr.p_double[i-1]+(2*ae_randomreal(_state)-1)*(x->ptr.p_double[i]-x->ptr.p_double[i-1]);
        }
    }
    else
    {
        x->ptr.p_double[0] = 0.5*(a+b);
        y->ptr.p_double[0] = 2*ae_randomreal(_state)-1;
    }
}


/*************************************************************************
This function generates  1-dimensional equidistant interpolation task with
moderate Lipshitz constant (close to 1.0)

If N=1 then suborutine generates only one point at the middle of [A,B]

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void taskgenint1dequidist(double a,
     double b,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t i;
    double h;

    ae_vector_clear(x);
    ae_vector_clear(y);

    ae_assert(n>=1, "TaskGenInterpolationEqdist1D: N<1!", _state);
    ae_vector_set_length(x, n, _state);
    ae_vector_set_length(y, n, _state);
    if( n>1 )
    {
        x->ptr.p_double[0] = a;
        y->ptr.p_double[0] = 2*ae_randomreal(_state)-1;
        h = (b-a)/(n-1);
        for(i=1; i<=n-1; i++)
        {
            x->ptr.p_double[i] = a+i*h;
            y->ptr.p_double[i] = y->ptr.p_double[i-1]+(2*ae_randomreal(_state)-1)*h;
        }
    }
    else
    {
        x->ptr.p_double[0] = 0.5*(a+b);
        y->ptr.p_double[0] = 2*ae_randomreal(_state)-1;
    }
}


/*************************************************************************
This function generates  1-dimensional Chebyshev-1 interpolation task with
moderate Lipshitz constant (close to 1.0)

If N=1 then suborutine generates only one point at the middle of [A,B]

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void taskgenint1dcheb1(double a,
     double b,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t i;

    ae_vector_clear(x);
    ae_vector_clear(y);

    ae_assert(n>=1, "TaskGenInterpolation1DCheb1: N<1!", _state);
    ae_vector_set_length(x, n, _state);
    ae_vector_set_length(y, n, _state);
    if( n>1 )
    {
        for(i=0; i<=n-1; i++)
        {
            x->ptr.p_double[i] = 0.5*(b+a)+0.5*(b-a)*ae_cos(ae_pi*(2*i+1)/(2*n), _state);
            if( i==0 )
            {
                y->ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            }
            else
            {
                y->ptr.p_double[i] = y->ptr.p_double[i-1]+(2*ae_randomreal(_state)-1)*(x->ptr.p_double[i]-x->ptr.p_double[i-1]);
            }
        }
    }
    else
    {
        x->ptr.p_double[0] = 0.5*(a+b);
        y->ptr.p_double[0] = 2*ae_randomreal(_state)-1;
    }
}


/*************************************************************************
This function generates  1-dimensional Chebyshev-2 interpolation task with
moderate Lipshitz constant (close to 1.0)

If N=1 then suborutine generates only one point at the middle of [A,B]

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
void taskgenint1dcheb2(double a,
     double b,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t i;

    ae_vector_clear(x);
    ae_vector_clear(y);

    ae_assert(n>=1, "TaskGenInterpolation1DCheb2: N<1!", _state);
    ae_vector_set_length(x, n, _state);
    ae_vector_set_length(y, n, _state);
    if( n>1 )
    {
        for(i=0; i<=n-1; i++)
        {
            x->ptr.p_double[i] = 0.5*(b+a)+0.5*(b-a)*ae_cos(ae_pi*i/(n-1), _state);
            if( i==0 )
            {
                y->ptr.p_double[i] = 2*ae_randomreal(_state)-1;
            }
            else
            {
                y->ptr.p_double[i] = y->ptr.p_double[i-1]+(2*ae_randomreal(_state)-1)*(x->ptr.p_double[i]-x->ptr.p_double[i-1]);
            }
        }
    }
    else
    {
        x->ptr.p_double[0] = 0.5*(a+b);
        y->ptr.p_double[0] = 2*ae_randomreal(_state)-1;
    }
}


/*************************************************************************
This function checks that all values from X[] are distinct. It does more
than just usual floating point comparison:
* first, it calculates max(X) and min(X)
* second, it maps X[] from [min,max] to [1,2]
* only at this stage actual comparison is done

The meaning of such check is to ensure that all values are "distinct enough"
and will not cause interpolation subroutine to fail.

NOTE:
    X[] must be sorted by ascending (subroutine ASSERT's it)

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool aredistinct(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_state *_state)
{
    double a;
    double b;
    ae_int_t i;
    ae_bool nonsorted;
    ae_bool result;


    ae_assert(n>=1, "APSERVAreDistinct: internal error (N<1)", _state);
    if( n==1 )
    {
        
        /*
         * everything is alright, it is up to caller to decide whether it
         * can interpolate something with just one point
         */
        result = ae_true;
        return result;
    }
    a = x->ptr.p_double[0];
    b = x->ptr.p_double[0];
    nonsorted = ae_false;
    for(i=1; i<=n-1; i++)
    {
        a = ae_minreal(a, x->ptr.p_double[i], _state);
        b = ae_maxreal(b, x->ptr.p_double[i], _state);
        nonsorted = nonsorted||ae_fp_greater_eq(x->ptr.p_double[i-1],x->ptr.p_double[i]);
    }
    ae_assert(!nonsorted, "APSERVAreDistinct: internal error (not sorted)", _state);
    for(i=1; i<=n-1; i++)
    {
        if( ae_fp_eq((x->ptr.p_double[i]-a)/(b-a)+1,(x->ptr.p_double[i-1]-a)/(b-a)+1) )
        {
            result = ae_false;
            return result;
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
This function checks that two boolean values are the same (both  are  True 
or both are False).

  -- ALGLIB --
     Copyright 02.12.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool aresameboolean(ae_bool v1, ae_bool v2, ae_state *) //_state)
{
    ae_bool result;


    result = (v1&&v2)||(!v1&&!v2);
    return result;
}


/*************************************************************************
If Length(X)<N, resizes X

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void bvectorsetlengthatleast(/* Boolean */ ae_vector* x,
     ae_int_t n,
     ae_state *_state)
{


    if( x->cnt<n )
    {
        ae_vector_set_length(x, n, _state);
    }
}


/*************************************************************************
If Length(X)<N, resizes X

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void ivectorsetlengthatleast(/* Integer */ ae_vector* x,
     ae_int_t n,
     ae_state *_state)
{


    if( x->cnt<n )
    {
        ae_vector_set_length(x, n, _state);
    }
}


/*************************************************************************
If Length(X)<N, resizes X

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void rvectorsetlengthatleast(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_state *_state)
{


    if( x->cnt<n )
    {
        ae_vector_set_length(x, n, _state);
    }
}


/*************************************************************************
If Cols(X)<N or Rows(X)<M, resizes X

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void rmatrixsetlengthatleast(/* Real    */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state)
{


    if( m>0&&n>0 )
    {
        if( x->rows<m||x->cols<n )
        {
            ae_matrix_set_length(x, m, n, _state);
        }
    }
}


/*************************************************************************
Resizes X and:
* preserves old contents of X
* fills new elements by zeros

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void rmatrixresize(/* Real    */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix oldx;
    ae_int_t i;
    ae_int_t j;
    ae_int_t m2;
    ae_int_t n2;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init(&oldx, 0, 0, DT_REAL, _state, ae_true);

    m2 = x->rows;
    n2 = x->cols;
    ae_swap_matrices(x, &oldx);
    ae_matrix_set_length(x, m, n, _state);
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            if( i<m2&&j<n2 )
            {
                x->ptr.pp_double[i][j] = oldx.ptr.pp_double[i][j];
            }
            else
            {
                x->ptr.pp_double[i][j] = 0.0;
            }
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This function checks that length(X) is at least N and first N values  from
X[] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool isfinitevector(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_bool result;


    ae_assert(n>=0, "APSERVIsFiniteVector: internal error (N<0)", _state);
    if( n==0 )
    {
        result = ae_true;
        return result;
    }
    if( x->cnt<n )
    {
        result = ae_false;
        return result;
    }
    for(i=0; i<=n-1; i++)
    {
        if( !ae_isfinite(x->ptr.p_double[i], _state) )
        {
            result = ae_false;
            return result;
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
This function checks that first N values from X[] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool isfinitecvector(/* Complex */ ae_vector* z,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_bool result;


    ae_assert(n>=0, "APSERVIsFiniteCVector: internal error (N<0)", _state);
    for(i=0; i<=n-1; i++)
    {
        if( !ae_isfinite(z->ptr.p_complex[i].x, _state)||!ae_isfinite(z->ptr.p_complex[i].y, _state) )
        {
            result = ae_false;
            return result;
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
This function checks that size of X is at least MxN and values from
X[0..M-1,0..N-1] are finite.

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool apservisfinitematrix(/* Real    */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_bool result;


    ae_assert(n>=0, "APSERVIsFiniteMatrix: internal error (N<0)", _state);
    ae_assert(m>=0, "APSERVIsFiniteMatrix: internal error (M<0)", _state);
    if( m==0||n==0 )
    {
        result = ae_true;
        return result;
    }
    if( x->rows<m||x->cols<n )
    {
        result = ae_false;
        return result;
    }
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            if( !ae_isfinite(x->ptr.pp_double[i][j], _state) )
            {
                result = ae_false;
                return result;
            }
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
This function checks that all values from X[0..M-1,0..N-1] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool apservisfinitecmatrix(/* Complex */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_bool result;


    ae_assert(n>=0, "APSERVIsFiniteCMatrix: internal error (N<0)", _state);
    ae_assert(m>=0, "APSERVIsFiniteCMatrix: internal error (M<0)", _state);
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            if( !ae_isfinite(x->ptr.pp_complex[i][j].x, _state)||!ae_isfinite(x->ptr.pp_complex[i][j].y, _state) )
            {
                result = ae_false;
                return result;
            }
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
This function checks that size of X is at least NxN and all values from
upper/lower triangle of X[0..N-1,0..N-1] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool isfinitertrmatrix(/* Real    */ ae_matrix* x,
     ae_int_t n,
     ae_bool isupper,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j1;
    ae_int_t j2;
    ae_int_t j;
    ae_bool result;


    ae_assert(n>=0, "APSERVIsFiniteRTRMatrix: internal error (N<0)", _state);
    if( n==0 )
    {
        result = ae_true;
        return result;
    }
    if( x->rows<n||x->cols<n )
    {
        result = ae_false;
        return result;
    }
    for(i=0; i<=n-1; i++)
    {
        if( isupper )
        {
            j1 = i;
            j2 = n-1;
        }
        else
        {
            j1 = 0;
            j2 = i;
        }
        for(j=j1; j<=j2; j++)
        {
            if( !ae_isfinite(x->ptr.pp_double[i][j], _state) )
            {
                result = ae_false;
                return result;
            }
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
This function checks that all values from upper/lower triangle of
X[0..N-1,0..N-1] are finite

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool apservisfinitectrmatrix(/* Complex */ ae_matrix* x,
     ae_int_t n,
     ae_bool isupper,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j1;
    ae_int_t j2;
    ae_int_t j;
    ae_bool result;


    ae_assert(n>=0, "APSERVIsFiniteCTRMatrix: internal error (N<0)", _state);
    for(i=0; i<=n-1; i++)
    {
        if( isupper )
        {
            j1 = i;
            j2 = n-1;
        }
        else
        {
            j1 = 0;
            j2 = i;
        }
        for(j=j1; j<=j2; j++)
        {
            if( !ae_isfinite(x->ptr.pp_complex[i][j].x, _state)||!ae_isfinite(x->ptr.pp_complex[i][j].y, _state) )
            {
                result = ae_false;
                return result;
            }
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
This function checks that all values from X[0..M-1,0..N-1] are  finite  or
NaN's.

  -- ALGLIB --
     Copyright 18.06.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool apservisfiniteornanmatrix(/* Real    */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_bool result;


    ae_assert(n>=0, "APSERVIsFiniteOrNaNMatrix: internal error (N<0)", _state);
    ae_assert(m>=0, "APSERVIsFiniteOrNaNMatrix: internal error (M<0)", _state);
    for(i=0; i<=m-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            if( !(ae_isfinite(x->ptr.pp_double[i][j], _state)||ae_isnan(x->ptr.pp_double[i][j], _state)) )
            {
                result = ae_false;
                return result;
            }
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
Safe sqrt(x^2+y^2)

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************/
double safepythag2(double x, double y, ae_state *_state)
{
    double w;
    double xabs;
    double yabs;
    double z;
    double result;


    xabs = ae_fabs(x, _state);
    yabs = ae_fabs(y, _state);
    w = ae_maxreal(xabs, yabs, _state);
    z = ae_minreal(xabs, yabs, _state);
    if( ae_fp_eq(z,0) )
    {
        result = w;
    }
    else
    {
        result = w*ae_sqrt(1+ae_sqr(z/w, _state), _state);
    }
    return result;
}


/*************************************************************************
Safe sqrt(x^2+y^2)

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************/
double safepythag3(double x, double y, double z, ae_state *_state)
{
    double w;
    double result;


    w = ae_maxreal(ae_fabs(x, _state), ae_maxreal(ae_fabs(y, _state), ae_fabs(z, _state), _state), _state);
    if( ae_fp_eq(w,0) )
    {
        result = 0;
        return result;
    }
    x = x/w;
    y = y/w;
    z = z/w;
    result = w*ae_sqrt(ae_sqr(x, _state)+ae_sqr(y, _state)+ae_sqr(z, _state), _state);
    return result;
}


/*************************************************************************
Safe division.

This function attempts to calculate R=X/Y without overflow.

It returns:
* +1, if abs(X/Y)>=MaxRealNumber or undefined - overflow-like situation
      (no overlfow is generated, R is either NAN, PosINF, NegINF)
*  0, if MinRealNumber<abs(X/Y)<MaxRealNumber or X=0, Y<>0
      (R contains result, may be zero)
* -1, if 0<abs(X/Y)<MinRealNumber - underflow-like situation
      (R contains zero; it corresponds to underflow)

No overflow is generated in any case.

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************/
ae_int_t saferdiv(double x, double y, double* r, ae_state *_state)
{
    ae_int_t result;

    *r = 0;

    
    /*
     * Two special cases:
     * * Y=0
     * * X=0 and Y<>0
     */
    if( ae_fp_eq(y,0) )
    {
        result = 1;
        if( ae_fp_eq(x,0) )
        {
            *r = _state->v_nan;
        }
        if( ae_fp_greater(x,0) )
        {
            *r = _state->v_posinf;
        }
        if( ae_fp_less(x,0) )
        {
            *r = _state->v_neginf;
        }
        return result;
    }
    if( ae_fp_eq(x,0) )
    {
        *r = 0;
        result = 0;
        return result;
    }
    
    /*
     * make Y>0
     */
    if( ae_fp_less(y,0) )
    {
        x = -x;
        y = -y;
    }
    
    /*
     *
     */
    if( ae_fp_greater_eq(y,1) )
    {
        *r = x/y;
        if( ae_fp_less_eq(ae_fabs(*r, _state),ae_minrealnumber) )
        {
            result = -1;
            *r = 0;
        }
        else
        {
            result = 0;
        }
    }
    else
    {
        if( ae_fp_greater_eq(ae_fabs(x, _state),ae_maxrealnumber*y) )
        {
            if( ae_fp_greater(x,0) )
            {
                *r = _state->v_posinf;
            }
            else
            {
                *r = _state->v_neginf;
            }
            result = 1;
        }
        else
        {
            *r = x/y;
            result = 0;
        }
    }
    return result;
}


/*************************************************************************
This function calculates "safe" min(X/Y,V) for positive finite X, Y, V.
No overflow is generated in any case.

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************/
double safeminposrv(double x, double y, double v, ae_state *) //_state)
{
    double r;
    double result;


    if( ae_fp_greater_eq(y,1) )
    {
        
        /*
         * Y>=1, we can safely divide by Y
         */
        r = x/y;
        result = v;
        if( ae_fp_greater(v,r) )
        {
            result = r;
        }
        else
        {
            result = v;
        }
    }
    else
    {
        
        /*
         * Y<1, we can safely multiply by Y
         */
        if( ae_fp_less(x,v*y) )
        {
            result = x/y;
        }
        else
        {
            result = v;
        }
    }
    return result;
}


/*************************************************************************
This function makes periodic mapping of X to [A,B].

It accepts X, A, B (A>B). It returns T which lies in  [A,B] and integer K,
such that X = T + K*(B-A).

NOTES:
* K is represented as real value, although actually it is integer
* T is guaranteed to be in [A,B]
* T replaces X

  -- ALGLIB --
     Copyright by Bochkanov Sergey
*************************************************************************/
void apperiodicmap(double* x,
     double a,
     double b,
     double* k,
     ae_state *_state)
{

    *k = 0;

    ae_assert(ae_fp_less(a,b), "APPeriodicMap: internal error!", _state);
    *k = ae_ifloor((*x-a)/(b-a), _state);
    *x = *x-*k*(b-a);
    while(ae_fp_less(*x,a))
    {
        *x = *x+(b-a);
        *k = *k-1;
    }
    while(ae_fp_greater(*x,b))
    {
        *x = *x-(b-a);
        *k = *k+1;
    }
    *x = ae_maxreal(*x, a, _state);
    *x = ae_minreal(*x, b, _state);
}


/*************************************************************************
'bounds' value: maps X to [B1,B2]

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
double boundval(double x, double b1, double b2, ae_state *) //_state)
{
    double result;


    if( ae_fp_less_eq(x,b1) )
    {
        result = b1;
        return result;
    }
    if( ae_fp_greater_eq(x,b2) )
    {
        result = b2;
        return result;
    }
    result = x;
    return result;
}


/*************************************************************************
Allocation of serializer: complex value
*************************************************************************/
void alloccomplex(ae_serializer* s, ae_complex /*v*/, ae_state *) //_state)
{


    ae_serializer_alloc_entry(s);
    ae_serializer_alloc_entry(s);
}


/*************************************************************************
Serialization: complex value
*************************************************************************/
void serializecomplex(ae_serializer* s, ae_complex v, ae_state *_state)
{


    ae_serializer_serialize_double(s, v.x, _state);
    ae_serializer_serialize_double(s, v.y, _state);
}


/*************************************************************************
Unserialization: complex value
*************************************************************************/
ae_complex unserializecomplex(ae_serializer* s, ae_state *_state)
{
    ae_complex result;


    ae_serializer_unserialize_double(s, &result.x, _state);
    ae_serializer_unserialize_double(s, &result.y, _state);
    return result;
}


/*************************************************************************
Allocation of serializer: real array
*************************************************************************/
void allocrealarray(ae_serializer* s,
     /* Real    */ ae_vector* v,
     ae_int_t n,
     ae_state *) //_state)
{
    ae_int_t i;


    if( n<0 )
    {
        n = v->cnt;
    }
    ae_serializer_alloc_entry(s);
    for(i=0; i<=n-1; i++)
    {
        ae_serializer_alloc_entry(s);
    }
}


/*************************************************************************
Serialization: complex value
*************************************************************************/
void serializerealarray(ae_serializer* s,
     /* Real    */ ae_vector* v,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;


    if( n<0 )
    {
        n = v->cnt;
    }
    ae_serializer_serialize_int(s, n, _state);
    for(i=0; i<=n-1; i++)
    {
        ae_serializer_serialize_double(s, v->ptr.p_double[i], _state);
    }
}


/*************************************************************************
Unserialization: complex value
*************************************************************************/
void unserializerealarray(ae_serializer* s,
     /* Real    */ ae_vector* v,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t i;
    double t;

    ae_vector_clear(v);

    ae_serializer_unserialize_int(s, &n, _state);
    if( n==0 )
    {
        return;
    }
    ae_vector_set_length(v, n, _state);
    for(i=0; i<=n-1; i++)
    {
        ae_serializer_unserialize_double(s, &t, _state);
        v->ptr.p_double[i] = t;
    }
}


/*************************************************************************
Allocation of serializer: Integer array
*************************************************************************/
void allocintegerarray(ae_serializer* s,
     /* Integer */ ae_vector* v,
     ae_int_t n,
     ae_state *) //_state)
{
    ae_int_t i;


    if( n<0 )
    {
        n = v->cnt;
    }
    ae_serializer_alloc_entry(s);
    for(i=0; i<=n-1; i++)
    {
        ae_serializer_alloc_entry(s);
    }
}


/*************************************************************************
Serialization: Integer array
*************************************************************************/
void serializeintegerarray(ae_serializer* s,
     /* Integer */ ae_vector* v,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;


    if( n<0 )
    {
        n = v->cnt;
    }
    ae_serializer_serialize_int(s, n, _state);
    for(i=0; i<=n-1; i++)
    {
        ae_serializer_serialize_int(s, v->ptr.p_int[i], _state);
    }
}


/*************************************************************************
Unserialization: complex value
*************************************************************************/
void unserializeintegerarray(ae_serializer* s,
     /* Integer */ ae_vector* v,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t i;
    ae_int_t t;

    ae_vector_clear(v);

    ae_serializer_unserialize_int(s, &n, _state);
    if( n==0 )
    {
        return;
    }
    ae_vector_set_length(v, n, _state);
    for(i=0; i<=n-1; i++)
    {
        ae_serializer_unserialize_int(s, &t, _state);
        v->ptr.p_int[i] = t;
    }
}


/*************************************************************************
Allocation of serializer: real matrix
*************************************************************************/
void allocrealmatrix(ae_serializer* s,
     /* Real    */ ae_matrix* v,
     ae_int_t n0,
     ae_int_t n1,
     ae_state *) //_state)
{
    ae_int_t i;
    ae_int_t j;


    if( n0<0 )
    {
        n0 = v->rows;
    }
    if( n1<0 )
    {
        n1 = v->cols;
    }
    ae_serializer_alloc_entry(s);
    ae_serializer_alloc_entry(s);
    for(i=0; i<=n0-1; i++)
    {
        for(j=0; j<=n1-1; j++)
        {
            ae_serializer_alloc_entry(s);
        }
    }
}


/*************************************************************************
Serialization: complex value
*************************************************************************/
void serializerealmatrix(ae_serializer* s,
     /* Real    */ ae_matrix* v,
     ae_int_t n0,
     ae_int_t n1,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;


    if( n0<0 )
    {
        n0 = v->rows;
    }
    if( n1<0 )
    {
        n1 = v->cols;
    }
    ae_serializer_serialize_int(s, n0, _state);
    ae_serializer_serialize_int(s, n1, _state);
    for(i=0; i<=n0-1; i++)
    {
        for(j=0; j<=n1-1; j++)
        {
            ae_serializer_serialize_double(s, v->ptr.pp_double[i][j], _state);
        }
    }
}


/*************************************************************************
Unserialization: complex value
*************************************************************************/
void unserializerealmatrix(ae_serializer* s,
     /* Real    */ ae_matrix* v,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t n0;
    ae_int_t n1;
    double t;

    ae_matrix_clear(v);

    ae_serializer_unserialize_int(s, &n0, _state);
    ae_serializer_unserialize_int(s, &n1, _state);
    if( n0==0||n1==0 )
    {
        return;
    }
    ae_matrix_set_length(v, n0, n1, _state);
    for(i=0; i<=n0-1; i++)
    {
        for(j=0; j<=n1-1; j++)
        {
            ae_serializer_unserialize_double(s, &t, _state);
            v->ptr.pp_double[i][j] = t;
        }
    }
}


/*************************************************************************
Copy integer array
*************************************************************************/
void copyintegerarray(/* Integer */ ae_vector* src,
     /* Integer */ ae_vector* dst,
     ae_state *_state)
{
    ae_int_t i;

    ae_vector_clear(dst);

    if( src->cnt>0 )
    {
        ae_vector_set_length(dst, src->cnt, _state);
        for(i=0; i<=src->cnt-1; i++)
        {
            dst->ptr.p_int[i] = src->ptr.p_int[i];
        }
    }
}


/*************************************************************************
Copy real array
*************************************************************************/
void copyrealarray(/* Real    */ ae_vector* src,
     /* Real    */ ae_vector* dst,
     ae_state *_state)
{
    ae_int_t i;

    ae_vector_clear(dst);

    if( src->cnt>0 )
    {
        ae_vector_set_length(dst, src->cnt, _state);
        for(i=0; i<=src->cnt-1; i++)
        {
            dst->ptr.p_double[i] = src->ptr.p_double[i];
        }
    }
}


/*************************************************************************
Copy real matrix
*************************************************************************/
void copyrealmatrix(/* Real    */ ae_matrix* src,
     /* Real    */ ae_matrix* dst,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;

    ae_matrix_clear(dst);

    if( src->rows>0&&src->cols>0 )
    {
        ae_matrix_set_length(dst, src->rows, src->cols, _state);
        for(i=0; i<=src->rows-1; i++)
        {
            for(j=0; j<=src->cols-1; j++)
            {
                dst->ptr.pp_double[i][j] = src->ptr.pp_double[i][j];
            }
        }
    }
}


/*************************************************************************
This function searches integer array. Elements in this array are actually
records, each NRec elements wide. Each record has unique header - NHeader
integer values, which identify it. Records are lexicographically sorted by
header.

Records are identified by their index, not offset (offset = NRec*index).

This function searches A (records with indices [I0,I1)) for a record with
header B. It returns index of this record (not offset!), or -1 on failure.

  -- ALGLIB --
     Copyright 28.03.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t recsearch(/* Integer */ ae_vector* a,
     ae_int_t nrec,
     ae_int_t nheader,
     ae_int_t i0,
     ae_int_t i1,
     /* Integer */ ae_vector* b,
     ae_state *) //_state)
{
    ae_int_t mididx;
    ae_int_t cflag;
    ae_int_t k;
    ae_int_t offs;
    ae_int_t result;


    result = -1;
    for(;;)
    {
        if( i0>=i1 )
        {
            break;
        }
        mididx = (i0+i1)/2;
        offs = nrec*mididx;
        cflag = 0;
        for(k=0; k<=nheader-1; k++)
        {
            if( a->ptr.p_int[offs+k]<b->ptr.p_int[k] )
            {
                cflag = -1;
                break;
            }
            if( a->ptr.p_int[offs+k]>b->ptr.p_int[k] )
            {
                cflag = 1;
                break;
            }
        }
        if( cflag==0 )
        {
            result = mididx;
            return result;
        }
        if( cflag<0 )
        {
            i0 = mididx+1;
        }
        else
        {
            i1 = mididx;
        }
    }
    return result;
}


ae_bool _apbuffers_init(apbuffers* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->ia0, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->ia1, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->ia2, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->ia3, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->ra0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->ra1, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->ra2, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->ra3, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _apbuffers_init_copy(apbuffers* dst, apbuffers* src, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init_copy(&dst->ia0, &src->ia0, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->ia1, &src->ia1, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->ia2, &src->ia2, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->ia3, &src->ia3, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->ra0, &src->ra0, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->ra1, &src->ra1, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->ra2, &src->ra2, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->ra3, &src->ra3, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _apbuffers_clear(apbuffers* p)
{
    ae_vector_clear(&p->ia0);
    ae_vector_clear(&p->ia1);
    ae_vector_clear(&p->ia2);
    ae_vector_clear(&p->ia3);
    ae_vector_clear(&p->ra0);
    ae_vector_clear(&p->ra1);
    ae_vector_clear(&p->ra2);
    ae_vector_clear(&p->ra3);
}




/*************************************************************************
This function sorts array of real keys by ascending.

Its results are:
* sorted array A
* permutation tables P1, P2

Algorithm outputs permutation tables using two formats:
* as usual permutation of [0..N-1]. If P1[i]=j, then sorted A[i]  contains
  value which was moved there from J-th position.
* as a sequence of pairwise permutations. Sorted A[] may  be  obtained  by
  swaping A[i] and A[P2[i]] for all i from 0 to N-1.
  
INPUT PARAMETERS:
    A       -   unsorted array
    N       -   array size

OUPUT PARAMETERS:
    A       -   sorted array
    P1, P2  -   permutation tables, array[N]
    
NOTES:
    this function assumes that A[] is finite; it doesn't checks that
    condition. All other conditions (size of input arrays, etc.) are not
    checked too.

  -- ALGLIB --
     Copyright 14.05.2008 by Bochkanov Sergey
*************************************************************************/
void tagsort(/* Real    */ ae_vector* a,
     ae_int_t n,
     /* Integer */ ae_vector* p1,
     /* Integer */ ae_vector* p2,
     ae_state *_state)
{
    ae_frame _frame_block;
    apbuffers buf;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(p1);
    ae_vector_clear(p2);
    _apbuffers_init(&buf, _state, ae_true);

    tagsortbuf(a, n, p1, p2, &buf, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Buffered variant of TagSort, which accepts preallocated output arrays as
well as special structure for buffered allocations. If arrays are too
short, they are reallocated. If they are large enough, no memory
allocation is done.

It is intended to be used in the performance-critical parts of code, where
additional allocations can lead to severe performance degradation

  -- ALGLIB --
     Copyright 14.05.2008 by Bochkanov Sergey
*************************************************************************/
void tagsortbuf(/* Real    */ ae_vector* a,
     ae_int_t n,
     /* Integer */ ae_vector* p1,
     /* Integer */ ae_vector* p2,
     apbuffers* buf,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t lv;
    ae_int_t lp;
    ae_int_t rv;
    ae_int_t rp;


    
    /*
     * Special cases
     */
    if( n<=0 )
    {
        return;
    }
    if( n==1 )
    {
        ivectorsetlengthatleast(p1, 1, _state);
        ivectorsetlengthatleast(p2, 1, _state);
        p1->ptr.p_int[0] = 0;
        p2->ptr.p_int[0] = 0;
        return;
    }
    
    /*
     * General case, N>1: prepare permutations table P1
     */
    ivectorsetlengthatleast(p1, n, _state);
    for(i=0; i<=n-1; i++)
    {
        p1->ptr.p_int[i] = i;
    }
    
    /*
     * General case, N>1: sort, update P1
     */
    rvectorsetlengthatleast(&buf->ra0, n, _state);
    ivectorsetlengthatleast(&buf->ia0, n, _state);
    tagsortfasti(a, p1, &buf->ra0, &buf->ia0, n, _state);
    
    /*
     * General case, N>1: fill permutations table P2
     *
     * To fill P2 we maintain two arrays:
     * * PV (Buf.IA0), Position(Value). PV[i] contains position of I-th key at the moment
     * * VP (Buf.IA1), Value(Position). VP[i] contains key which has position I at the moment
     *
     * At each step we making permutation of two items:
     *   Left, which is given by position/value pair LP/LV
     *   and Right, which is given by RP/RV
     * and updating PV[] and VP[] correspondingly.
     */
    ivectorsetlengthatleast(&buf->ia0, n, _state);
    ivectorsetlengthatleast(&buf->ia1, n, _state);
    ivectorsetlengthatleast(p2, n, _state);
    for(i=0; i<=n-1; i++)
    {
        buf->ia0.ptr.p_int[i] = i;
        buf->ia1.ptr.p_int[i] = i;
    }
    for(i=0; i<=n-1; i++)
    {
        
        /*
         * calculate LP, LV, RP, RV
         */
        lp = i;
        lv = buf->ia1.ptr.p_int[lp];
        rv = p1->ptr.p_int[i];
        rp = buf->ia0.ptr.p_int[rv];
        
        /*
         * Fill P2
         */
        p2->ptr.p_int[i] = rp;
        
        /*
         * update PV and VP
         */
        buf->ia1.ptr.p_int[lp] = rv;
        buf->ia1.ptr.p_int[rp] = lv;
        buf->ia0.ptr.p_int[lv] = rp;
        buf->ia0.ptr.p_int[rv] = lp;
    }
}


/*************************************************************************
Same as TagSort, but optimized for real keys and integer labels.

A is sorted, and same permutations are applied to B.

NOTES:
1.  this function assumes that A[] is finite; it doesn't checks that
    condition. All other conditions (size of input arrays, etc.) are not
    checked too.
2.  this function uses two buffers, BufA and BufB, each is N elements large.
    They may be preallocated (which will save some time) or not, in which
    case function will automatically allocate memory.

  -- ALGLIB --
     Copyright 11.12.2008 by Bochkanov Sergey
*************************************************************************/
void tagsortfasti(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     /* Real    */ ae_vector* bufa,
     /* Integer */ ae_vector* bufb,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_bool isascending;
    ae_bool isdescending;
    double tmpr;
    ae_int_t tmpi;


    
    /*
     * Special case
     */
    if( n<=1 )
    {
        return;
    }
    
    /*
     * Test for already sorted set
     */
    isascending = ae_true;
    isdescending = ae_true;
    for(i=1; i<=n-1; i++)
    {
        isascending = isascending&&a->ptr.p_double[i]>=a->ptr.p_double[i-1];
        isdescending = isdescending&&a->ptr.p_double[i]<=a->ptr.p_double[i-1];
    }
    if( isascending )
    {
        return;
    }
    if( isdescending )
    {
        for(i=0; i<=n-1; i++)
        {
            j = n-1-i;
            if( j<=i )
            {
                break;
            }
            tmpr = a->ptr.p_double[i];
            a->ptr.p_double[i] = a->ptr.p_double[j];
            a->ptr.p_double[j] = tmpr;
            tmpi = b->ptr.p_int[i];
            b->ptr.p_int[i] = b->ptr.p_int[j];
            b->ptr.p_int[j] = tmpi;
        }
        return;
    }
    
    /*
     * General case
     */
    if( bufa->cnt<n )
    {
        ae_vector_set_length(bufa, n, _state);
    }
    if( bufb->cnt<n )
    {
        ae_vector_set_length(bufb, n, _state);
    }
    tsort_tagsortfastirec(a, b, bufa, bufb, 0, n-1, _state);
}


/*************************************************************************
Same as TagSort, but optimized for real keys and real labels.

A is sorted, and same permutations are applied to B.

NOTES:
1.  this function assumes that A[] is finite; it doesn't checks that
    condition. All other conditions (size of input arrays, etc.) are not
    checked too.
2.  this function uses two buffers, BufA and BufB, each is N elements large.
    They may be preallocated (which will save some time) or not, in which
    case function will automatically allocate memory.

  -- ALGLIB --
     Copyright 11.12.2008 by Bochkanov Sergey
*************************************************************************/
void tagsortfastr(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* bufa,
     /* Real    */ ae_vector* bufb,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_bool isascending;
    ae_bool isdescending;
    double tmpr;


    
    /*
     * Special case
     */
    if( n<=1 )
    {
        return;
    }
    
    /*
     * Test for already sorted set
     */
    isascending = ae_true;
    isdescending = ae_true;
    for(i=1; i<=n-1; i++)
    {
        isascending = isascending&&a->ptr.p_double[i]>=a->ptr.p_double[i-1];
        isdescending = isdescending&&a->ptr.p_double[i]<=a->ptr.p_double[i-1];
    }
    if( isascending )
    {
        return;
    }
    if( isdescending )
    {
        for(i=0; i<=n-1; i++)
        {
            j = n-1-i;
            if( j<=i )
            {
                break;
            }
            tmpr = a->ptr.p_double[i];
            a->ptr.p_double[i] = a->ptr.p_double[j];
            a->ptr.p_double[j] = tmpr;
            tmpr = b->ptr.p_double[i];
            b->ptr.p_double[i] = b->ptr.p_double[j];
            b->ptr.p_double[j] = tmpr;
        }
        return;
    }
    
    /*
     * General case
     */
    if( bufa->cnt<n )
    {
        ae_vector_set_length(bufa, n, _state);
    }
    if( bufb->cnt<n )
    {
        ae_vector_set_length(bufb, n, _state);
    }
    tsort_tagsortfastrrec(a, b, bufa, bufb, 0, n-1, _state);
}


/*************************************************************************
Same as TagSort, but optimized for real keys without labels.

A is sorted, and that's all.

NOTES:
1.  this function assumes that A[] is finite; it doesn't checks that
    condition. All other conditions (size of input arrays, etc.) are not
    checked too.
2.  this function uses buffer, BufA, which is N elements large. It may be
    preallocated (which will save some time) or not, in which case
    function will automatically allocate memory.

  -- ALGLIB --
     Copyright 11.12.2008 by Bochkanov Sergey
*************************************************************************/
void tagsortfast(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* bufa,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_bool isascending;
    ae_bool isdescending;
    double tmpr;


    
    /*
     * Special case
     */
    if( n<=1 )
    {
        return;
    }
    
    /*
     * Test for already sorted set
     */
    isascending = ae_true;
    isdescending = ae_true;
    for(i=1; i<=n-1; i++)
    {
        isascending = isascending&&a->ptr.p_double[i]>=a->ptr.p_double[i-1];
        isdescending = isdescending&&a->ptr.p_double[i]<=a->ptr.p_double[i-1];
    }
    if( isascending )
    {
        return;
    }
    if( isdescending )
    {
        for(i=0; i<=n-1; i++)
        {
            j = n-1-i;
            if( j<=i )
            {
                break;
            }
            tmpr = a->ptr.p_double[i];
            a->ptr.p_double[i] = a->ptr.p_double[j];
            a->ptr.p_double[j] = tmpr;
        }
        return;
    }
    
    /*
     * General case
     */
    if( bufa->cnt<n )
    {
        ae_vector_set_length(bufa, n, _state);
    }
    tsort_tagsortfastrec(a, bufa, 0, n-1, _state);
}


/*************************************************************************
Sorting function optimized for integer keys and real labels, can be used
to sort middle of the array

A is sorted, and same permutations are applied to B.

NOTES:
    this function assumes that A[] is finite; it doesn't checks that
    condition. All other conditions (size of input arrays, etc.) are not
    checked too.

  -- ALGLIB --
     Copyright 11.12.2008 by Bochkanov Sergey
*************************************************************************/
void tagsortmiddleir(/* Integer */ ae_vector* a,
     /* Real    */ ae_vector* b,
     ae_int_t offset,
     ae_int_t n,
     ae_state *) //_state)
{
    ae_int_t i;
    ae_int_t k;
    ae_int_t t;
    ae_int_t tmp;
    double tmpr;


    
    /*
     * Special cases
     */
    if( n<=1 )
    {
        return;
    }
    
    /*
     * General case, N>1: sort, update B
     */
    i = 2;
    do
    {
        t = i;
        while(t!=1)
        {
            k = t/2;
            if( a->ptr.p_int[offset+k-1]>=a->ptr.p_int[offset+t-1] )
            {
                t = 1;
            }
            else
            {
                tmp = a->ptr.p_int[offset+k-1];
                a->ptr.p_int[offset+k-1] = a->ptr.p_int[offset+t-1];
                a->ptr.p_int[offset+t-1] = tmp;
                tmpr = b->ptr.p_double[offset+k-1];
                b->ptr.p_double[offset+k-1] = b->ptr.p_double[offset+t-1];
                b->ptr.p_double[offset+t-1] = tmpr;
                t = k;
            }
        }
        i = i+1;
    }
    while(i<=n);
    i = n-1;
    do
    {
        tmp = a->ptr.p_int[offset+i];
        a->ptr.p_int[offset+i] = a->ptr.p_int[offset+0];
        a->ptr.p_int[offset+0] = tmp;
        tmpr = b->ptr.p_double[offset+i];
        b->ptr.p_double[offset+i] = b->ptr.p_double[offset+0];
        b->ptr.p_double[offset+0] = tmpr;
        t = 1;
        while(t!=0)
        {
            k = 2*t;
            if( k>i )
            {
                t = 0;
            }
            else
            {
                if( k<i )
                {
                    if( a->ptr.p_int[offset+k]>a->ptr.p_int[offset+k-1] )
                    {
                        k = k+1;
                    }
                }
                if( a->ptr.p_int[offset+t-1]>=a->ptr.p_int[offset+k-1] )
                {
                    t = 0;
                }
                else
                {
                    tmp = a->ptr.p_int[offset+k-1];
                    a->ptr.p_int[offset+k-1] = a->ptr.p_int[offset+t-1];
                    a->ptr.p_int[offset+t-1] = tmp;
                    tmpr = b->ptr.p_double[offset+k-1];
                    b->ptr.p_double[offset+k-1] = b->ptr.p_double[offset+t-1];
                    b->ptr.p_double[offset+t-1] = tmpr;
                    t = k;
                }
            }
        }
        i = i-1;
    }
    while(i>=1);
}


/*************************************************************************
Heap operations: adds element to the heap

PARAMETERS:
    A       -   heap itself, must be at least array[0..N]
    B       -   array of integer tags, which are updated according to
                permutations in the heap
    N       -   size of the heap (without new element).
                updated on output
    VA      -   value of the element being added
    VB      -   value of the tag

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void tagheappushi(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     ae_int_t* n,
     double va,
     ae_int_t vb,
     ae_state *) //_state)
{
    ae_int_t j;
    ae_int_t k;
    double v;


    if( *n<0 )
    {
        return;
    }
    
    /*
     * N=0 is a special case
     */
    if( *n==0 )
    {
        a->ptr.p_double[0] = va;
        b->ptr.p_int[0] = vb;
        *n = *n+1;
        return;
    }
    
    /*
     * add current point to the heap
     * (add to the bottom, then move up)
     *
     * we don't write point to the heap
     * until its final position is determined
     * (it allow us to reduce number of array access operations)
     */
    j = *n;
    *n = *n+1;
    while(j>0)
    {
        k = (j-1)/2;
        v = a->ptr.p_double[k];
        if( ae_fp_less(v,va) )
        {
            
            /*
             * swap with higher element
             */
            a->ptr.p_double[j] = v;
            b->ptr.p_int[j] = b->ptr.p_int[k];
            j = k;
        }
        else
        {
            
            /*
             * element in its place. terminate.
             */
            break;
        }
    }
    a->ptr.p_double[j] = va;
    b->ptr.p_int[j] = vb;
}


/*************************************************************************
Heap operations: replaces top element with new element
(which is moved down)

PARAMETERS:
    A       -   heap itself, must be at least array[0..N-1]
    B       -   array of integer tags, which are updated according to
                permutations in the heap
    N       -   size of the heap
    VA      -   value of the element which replaces top element
    VB      -   value of the tag

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void tagheapreplacetopi(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     ae_int_t n,
     double va,
     ae_int_t vb,
     ae_state *) //_state)
{
    ae_int_t j;
    ae_int_t k1;
    ae_int_t k2;
    double v;
    double v1;
    double v2;


    if( n<1 )
    {
        return;
    }
    
    /*
     * N=1 is a special case
     */
    if( n==1 )
    {
        a->ptr.p_double[0] = va;
        b->ptr.p_int[0] = vb;
        return;
    }
    
    /*
     * move down through heap:
     * * J  -   current element
     * * K1 -   first child (always exists)
     * * K2 -   second child (may not exists)
     *
     * we don't write point to the heap
     * until its final position is determined
     * (it allow us to reduce number of array access operations)
     */
    j = 0;
    k1 = 1;
    k2 = 2;
    while(k1<n)
    {
        if( k2>=n )
        {
            
            /*
             * only one child.
             *
             * swap and terminate (because this child
             * have no siblings due to heap structure)
             */
            v = a->ptr.p_double[k1];
            if( ae_fp_greater(v,va) )
            {
                a->ptr.p_double[j] = v;
                b->ptr.p_int[j] = b->ptr.p_int[k1];
                j = k1;
            }
            break;
        }
        else
        {
            
            /*
             * two childs
             */
            v1 = a->ptr.p_double[k1];
            v2 = a->ptr.p_double[k2];
            if( ae_fp_greater(v1,v2) )
            {
                if( ae_fp_less(va,v1) )
                {
                    a->ptr.p_double[j] = v1;
                    b->ptr.p_int[j] = b->ptr.p_int[k1];
                    j = k1;
                }
                else
                {
                    break;
                }
            }
            else
            {
                if( ae_fp_less(va,v2) )
                {
                    a->ptr.p_double[j] = v2;
                    b->ptr.p_int[j] = b->ptr.p_int[k2];
                    j = k2;
                }
                else
                {
                    break;
                }
            }
            k1 = 2*j+1;
            k2 = 2*j+2;
        }
    }
    a->ptr.p_double[j] = va;
    b->ptr.p_int[j] = vb;
}


/*************************************************************************
Heap operations: pops top element from the heap

PARAMETERS:
    A       -   heap itself, must be at least array[0..N-1]
    B       -   array of integer tags, which are updated according to
                permutations in the heap
    N       -   size of the heap, N>=1

On output top element is moved to A[N-1], B[N-1], heap is reordered, N is
decreased by 1.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void tagheappopi(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     ae_int_t* n,
     ae_state *_state)
{
    double va;
    ae_int_t vb;


    if( *n<1 )
    {
        return;
    }
    
    /*
     * N=1 is a special case
     */
    if( *n==1 )
    {
        *n = 0;
        return;
    }
    
    /*
     * swap top element and last element,
     * then reorder heap
     */
    va = a->ptr.p_double[*n-1];
    vb = b->ptr.p_int[*n-1];
    a->ptr.p_double[*n-1] = a->ptr.p_double[0];
    b->ptr.p_int[*n-1] = b->ptr.p_int[0];
    *n = *n-1;
    tagheapreplacetopi(a, b, *n, va, vb, _state);
}


/*************************************************************************
Search first element less than T in sorted array.

PARAMETERS:
    A - sorted array by ascending from 0 to N-1
    N - number of elements in array
    T - the desired element

RESULT:
    The very first element's index, which isn't less than T.
In the case when there aren't such elements, returns N.
*************************************************************************/
ae_int_t lowerbound(/* Real    */ ae_vector* a,
     ae_int_t n,
     double t,
     ae_state *) //_state)
{
    ae_int_t l;
    ae_int_t half;
    ae_int_t first;
    ae_int_t middle;
    ae_int_t result;


    l = n;
    first = 0;
    while(l>0)
    {
        half = l/2;
        middle = first+half;
        if( ae_fp_less(a->ptr.p_double[middle],t) )
        {
            first = middle+1;
            l = l-half-1;
        }
        else
        {
            l = half;
        }
    }
    result = first;
    return result;
}


/*************************************************************************
Search first element more than T in sorted array.

PARAMETERS:
    A - sorted array by ascending from 0 to N-1
    N - number of elements in array
    T - the desired element

    RESULT:
    The very first element's index, which more than T.
In the case when there aren't such elements, returns N.
*************************************************************************/
ae_int_t upperbound(/* Real    */ ae_vector* a,
     ae_int_t n,
     double t,
     ae_state *) //_state)
{
    ae_int_t l;
    ae_int_t half;
    ae_int_t first;
    ae_int_t middle;
    ae_int_t result;


    l = n;
    first = 0;
    while(l>0)
    {
        half = l/2;
        middle = first+half;
        if( ae_fp_less(t,a->ptr.p_double[middle]) )
        {
            l = half;
        }
        else
        {
            first = middle+1;
            l = l-half-1;
        }
    }
    result = first;
    return result;
}


/*************************************************************************
Internal TagSortFastI: sorts A[I1...I2] (both bounds are included),
applies same permutations to B.

  -- ALGLIB --
     Copyright 06.09.2010 by Bochkanov Sergey
*************************************************************************/
static void tsort_tagsortfastirec(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     /* Real    */ ae_vector* bufa,
     /* Integer */ ae_vector* bufb,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t cntless;
    ae_int_t cnteq;
    ae_int_t cntgreater;
    double tmpr;
    ae_int_t tmpi;
    double v0;
    double v1;
    double v2;
    double vp;


    
    /*
     * Fast exit
     */
    if( i2<=i1 )
    {
        return;
    }
    
    /*
     * Non-recursive sort for small arrays
     */
    if( i2-i1<=16 )
    {
        for(j=i1+1; j<=i2; j++)
        {
            
            /*
             * Search elements [I1..J-1] for place to insert Jth element.
             *
             * This code stops immediately if we can leave A[J] at J-th position
             * (all elements have same value of A[J] larger than any of them)
             */
            tmpr = a->ptr.p_double[j];
            tmpi = j;
            for(k=j-1; k>=i1; k--)
            {
                if( a->ptr.p_double[k]<=tmpr )
                {
                    break;
                }
                tmpi = k;
            }
            k = tmpi;
            
            /*
             * Insert Jth element into Kth position
             */
            if( k!=j )
            {
                tmpr = a->ptr.p_double[j];
                tmpi = b->ptr.p_int[j];
                for(i=j-1; i>=k; i--)
                {
                    a->ptr.p_double[i+1] = a->ptr.p_double[i];
                    b->ptr.p_int[i+1] = b->ptr.p_int[i];
                }
                a->ptr.p_double[k] = tmpr;
                b->ptr.p_int[k] = tmpi;
            }
        }
        return;
    }
    
    /*
     * Quicksort: choose pivot
     * Here we assume that I2-I1>=2
     */
    v0 = a->ptr.p_double[i1];
    v1 = a->ptr.p_double[i1+(i2-i1)/2];
    v2 = a->ptr.p_double[i2];
    if( v0>v1 )
    {
        tmpr = v1;
        v1 = v0;
        v0 = tmpr;
    }
    if( v1>v2 )
    {
        tmpr = v2;
        v2 = v1;
        v1 = tmpr;
    }
    if( v0>v1 )
    {
        tmpr = v1;
        v1 = v0;
        v0 = tmpr;
    }
    vp = v1;
    
    /*
     * now pass through A/B and:
     * * move elements that are LESS than VP to the left of A/B
     * * move elements that are EQUAL to VP to the right of BufA/BufB (in the reverse order)
     * * move elements that are GREATER than VP to the left of BufA/BufB (in the normal order
     * * move elements from the tail of BufA/BufB to the middle of A/B (restoring normal order)
     * * move elements from the left of BufA/BufB to the end of A/B
     */
    cntless = 0;
    cnteq = 0;
    cntgreater = 0;
    for(i=i1; i<=i2; i++)
    {
        v0 = a->ptr.p_double[i];
        if( v0<vp )
        {
            
            /*
             * LESS
             */
            k = i1+cntless;
            if( i!=k )
            {
                a->ptr.p_double[k] = v0;
                b->ptr.p_int[k] = b->ptr.p_int[i];
            }
            cntless = cntless+1;
            continue;
        }
        if( v0==vp )
        {
            
            /*
             * EQUAL
             */
            k = i2-cnteq;
            bufa->ptr.p_double[k] = v0;
            bufb->ptr.p_int[k] = b->ptr.p_int[i];
            cnteq = cnteq+1;
            continue;
        }
        
        /*
         * GREATER
         */
        k = i1+cntgreater;
        bufa->ptr.p_double[k] = v0;
        bufb->ptr.p_int[k] = b->ptr.p_int[i];
        cntgreater = cntgreater+1;
    }
    for(i=0; i<=cnteq-1; i++)
    {
        j = i1+cntless+cnteq-1-i;
        k = i2+i-(cnteq-1);
        a->ptr.p_double[j] = bufa->ptr.p_double[k];
        b->ptr.p_int[j] = bufb->ptr.p_int[k];
    }
    for(i=0; i<=cntgreater-1; i++)
    {
        j = i1+cntless+cnteq+i;
        k = i1+i;
        a->ptr.p_double[j] = bufa->ptr.p_double[k];
        b->ptr.p_int[j] = bufb->ptr.p_int[k];
    }
    
    /*
     * Sort left and right parts of the array (ignoring middle part)
     */
    tsort_tagsortfastirec(a, b, bufa, bufb, i1, i1+cntless-1, _state);
    tsort_tagsortfastirec(a, b, bufa, bufb, i1+cntless+cnteq, i2, _state);
}


/*************************************************************************
Internal TagSortFastR: sorts A[I1...I2] (both bounds are included),
applies same permutations to B.

  -- ALGLIB --
     Copyright 06.09.2010 by Bochkanov Sergey
*************************************************************************/
static void tsort_tagsortfastrrec(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* bufa,
     /* Real    */ ae_vector* bufb,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    double tmpr;
    double tmpr2;
    ae_int_t tmpi;
    ae_int_t cntless;
    ae_int_t cnteq;
    ae_int_t cntgreater;
    double v0;
    double v1;
    double v2;
    double vp;


    
    /*
     * Fast exit
     */
    if( i2<=i1 )
    {
        return;
    }
    
    /*
     * Non-recursive sort for small arrays
     */
    if( i2-i1<=16 )
    {
        for(j=i1+1; j<=i2; j++)
        {
            
            /*
             * Search elements [I1..J-1] for place to insert Jth element.
             *
             * This code stops immediatly if we can leave A[J] at J-th position
             * (all elements have same value of A[J] larger than any of them)
             */
            tmpr = a->ptr.p_double[j];
            tmpi = j;
            for(k=j-1; k>=i1; k--)
            {
                if( a->ptr.p_double[k]<=tmpr )
                {
                    break;
                }
                tmpi = k;
            }
            k = tmpi;
            
            /*
             * Insert Jth element into Kth position
             */
            if( k!=j )
            {
                tmpr = a->ptr.p_double[j];
                tmpr2 = b->ptr.p_double[j];
                for(i=j-1; i>=k; i--)
                {
                    a->ptr.p_double[i+1] = a->ptr.p_double[i];
                    b->ptr.p_double[i+1] = b->ptr.p_double[i];
                }
                a->ptr.p_double[k] = tmpr;
                b->ptr.p_double[k] = tmpr2;
            }
        }
        return;
    }
    
    /*
     * Quicksort: choose pivot
     * Here we assume that I2-I1>=16
     */
    v0 = a->ptr.p_double[i1];
    v1 = a->ptr.p_double[i1+(i2-i1)/2];
    v2 = a->ptr.p_double[i2];
    if( v0>v1 )
    {
        tmpr = v1;
        v1 = v0;
        v0 = tmpr;
    }
    if( v1>v2 )
    {
        tmpr = v2;
        v2 = v1;
        v1 = tmpr;
    }
    if( v0>v1 )
    {
        tmpr = v1;
        v1 = v0;
        v0 = tmpr;
    }
    vp = v1;
    
    /*
     * now pass through A/B and:
     * * move elements that are LESS than VP to the left of A/B
     * * move elements that are EQUAL to VP to the right of BufA/BufB (in the reverse order)
     * * move elements that are GREATER than VP to the left of BufA/BufB (in the normal order
     * * move elements from the tail of BufA/BufB to the middle of A/B (restoring normal order)
     * * move elements from the left of BufA/BufB to the end of A/B
     */
    cntless = 0;
    cnteq = 0;
    cntgreater = 0;
    for(i=i1; i<=i2; i++)
    {
        v0 = a->ptr.p_double[i];
        if( v0<vp )
        {
            
            /*
             * LESS
             */
            k = i1+cntless;
            if( i!=k )
            {
                a->ptr.p_double[k] = v0;
                b->ptr.p_double[k] = b->ptr.p_double[i];
            }
            cntless = cntless+1;
            continue;
        }
        if( v0==vp )
        {
            
            /*
             * EQUAL
             */
            k = i2-cnteq;
            bufa->ptr.p_double[k] = v0;
            bufb->ptr.p_double[k] = b->ptr.p_double[i];
            cnteq = cnteq+1;
            continue;
        }
        
        /*
         * GREATER
         */
        k = i1+cntgreater;
        bufa->ptr.p_double[k] = v0;
        bufb->ptr.p_double[k] = b->ptr.p_double[i];
        cntgreater = cntgreater+1;
    }
    for(i=0; i<=cnteq-1; i++)
    {
        j = i1+cntless+cnteq-1-i;
        k = i2+i-(cnteq-1);
        a->ptr.p_double[j] = bufa->ptr.p_double[k];
        b->ptr.p_double[j] = bufb->ptr.p_double[k];
    }
    for(i=0; i<=cntgreater-1; i++)
    {
        j = i1+cntless+cnteq+i;
        k = i1+i;
        a->ptr.p_double[j] = bufa->ptr.p_double[k];
        b->ptr.p_double[j] = bufb->ptr.p_double[k];
    }
    
    /*
     * Sort left and right parts of the array (ignoring middle part)
     */
    tsort_tagsortfastrrec(a, b, bufa, bufb, i1, i1+cntless-1, _state);
    tsort_tagsortfastrrec(a, b, bufa, bufb, i1+cntless+cnteq, i2, _state);
}


/*************************************************************************
Internal TagSortFastI: sorts A[I1...I2] (both bounds are included),
applies same permutations to B.

  -- ALGLIB --
     Copyright 06.09.2010 by Bochkanov Sergey
*************************************************************************/
static void tsort_tagsortfastrec(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* bufa,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state)
{
    ae_int_t cntless;
    ae_int_t cnteq;
    ae_int_t cntgreater;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    double tmpr;
    ae_int_t tmpi;
    double v0;
    double v1;
    double v2;
    double vp;


    
    /*
     * Fast exit
     */
    if( i2<=i1 )
    {
        return;
    }
    
    /*
     * Non-recursive sort for small arrays
     */
    if( i2-i1<=16 )
    {
        for(j=i1+1; j<=i2; j++)
        {
            
            /*
             * Search elements [I1..J-1] for place to insert Jth element.
             *
             * This code stops immediatly if we can leave A[J] at J-th position
             * (all elements have same value of A[J] larger than any of them)
             */
            tmpr = a->ptr.p_double[j];
            tmpi = j;
            for(k=j-1; k>=i1; k--)
            {
                if( a->ptr.p_double[k]<=tmpr )
                {
                    break;
                }
                tmpi = k;
            }
            k = tmpi;
            
            /*
             * Insert Jth element into Kth position
             */
            if( k!=j )
            {
                tmpr = a->ptr.p_double[j];
                for(i=j-1; i>=k; i--)
                {
                    a->ptr.p_double[i+1] = a->ptr.p_double[i];
                }
                a->ptr.p_double[k] = tmpr;
            }
        }
        return;
    }
    
    /*
     * Quicksort: choose pivot
     * Here we assume that I2-I1>=16
     */
    v0 = a->ptr.p_double[i1];
    v1 = a->ptr.p_double[i1+(i2-i1)/2];
    v2 = a->ptr.p_double[i2];
    if( v0>v1 )
    {
        tmpr = v1;
        v1 = v0;
        v0 = tmpr;
    }
    if( v1>v2 )
    {
        tmpr = v2;
        v2 = v1;
        v1 = tmpr;
    }
    if( v0>v1 )
    {
        tmpr = v1;
        v1 = v0;
        v0 = tmpr;
    }
    vp = v1;
    
    /*
     * now pass through A/B and:
     * * move elements that are LESS than VP to the left of A/B
     * * move elements that are EQUAL to VP to the right of BufA/BufB (in the reverse order)
     * * move elements that are GREATER than VP to the left of BufA/BufB (in the normal order
     * * move elements from the tail of BufA/BufB to the middle of A/B (restoring normal order)
     * * move elements from the left of BufA/BufB to the end of A/B
     */
    cntless = 0;
    cnteq = 0;
    cntgreater = 0;
    for(i=i1; i<=i2; i++)
    {
        v0 = a->ptr.p_double[i];
        if( v0<vp )
        {
            
            /*
             * LESS
             */
            k = i1+cntless;
            if( i!=k )
            {
                a->ptr.p_double[k] = v0;
            }
            cntless = cntless+1;
            continue;
        }
        if( v0==vp )
        {
            
            /*
             * EQUAL
             */
            k = i2-cnteq;
            bufa->ptr.p_double[k] = v0;
            cnteq = cnteq+1;
            continue;
        }
        
        /*
         * GREATER
         */
        k = i1+cntgreater;
        bufa->ptr.p_double[k] = v0;
        cntgreater = cntgreater+1;
    }
    for(i=0; i<=cnteq-1; i++)
    {
        j = i1+cntless+cnteq-1-i;
        k = i2+i-(cnteq-1);
        a->ptr.p_double[j] = bufa->ptr.p_double[k];
    }
    for(i=0; i<=cntgreater-1; i++)
    {
        j = i1+cntless+cnteq+i;
        k = i1+i;
        a->ptr.p_double[j] = bufa->ptr.p_double[k];
    }
    
    /*
     * Sort left and right parts of the array (ignoring middle part)
     */
    tsort_tagsortfastrec(a, bufa, i1, i1+cntless-1, _state);
    tsort_tagsortfastrec(a, bufa, i1+cntless+cnteq, i2, _state);
}




/*************************************************************************
Internal ranking subroutine
*************************************************************************/
void rankx(/* Real    */ ae_vector* x,
     ae_int_t n,
     apbuffers* buf,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t t;
    double tmp;
    ae_int_t tmpi;


    
    /*
     * Prepare
     */
    if( n<1 )
    {
        return;
    }
    if( n==1 )
    {
        x->ptr.p_double[0] = 1;
        return;
    }
    if( buf->ra1.cnt<n )
    {
        ae_vector_set_length(&buf->ra1, n, _state);
    }
    if( buf->ia1.cnt<n )
    {
        ae_vector_set_length(&buf->ia1, n, _state);
    }
    for(i=0; i<=n-1; i++)
    {
        buf->ra1.ptr.p_double[i] = x->ptr.p_double[i];
        buf->ia1.ptr.p_int[i] = i;
    }
    
    /*
     * sort {R, C}
     */
    if( n!=1 )
    {
        i = 2;
        do
        {
            t = i;
            while(t!=1)
            {
                k = t/2;
                if( ae_fp_greater_eq(buf->ra1.ptr.p_double[k-1],buf->ra1.ptr.p_double[t-1]) )
                {
                    t = 1;
                }
                else
                {
                    tmp = buf->ra1.ptr.p_double[k-1];
                    buf->ra1.ptr.p_double[k-1] = buf->ra1.ptr.p_double[t-1];
                    buf->ra1.ptr.p_double[t-1] = tmp;
                    tmpi = buf->ia1.ptr.p_int[k-1];
                    buf->ia1.ptr.p_int[k-1] = buf->ia1.ptr.p_int[t-1];
                    buf->ia1.ptr.p_int[t-1] = tmpi;
                    t = k;
                }
            }
            i = i+1;
        }
        while(i<=n);
        i = n-1;
        do
        {
            tmp = buf->ra1.ptr.p_double[i];
            buf->ra1.ptr.p_double[i] = buf->ra1.ptr.p_double[0];
            buf->ra1.ptr.p_double[0] = tmp;
            tmpi = buf->ia1.ptr.p_int[i];
            buf->ia1.ptr.p_int[i] = buf->ia1.ptr.p_int[0];
            buf->ia1.ptr.p_int[0] = tmpi;
            t = 1;
            while(t!=0)
            {
                k = 2*t;
                if( k>i )
                {
                    t = 0;
                }
                else
                {
                    if( k<i )
                    {
                        if( ae_fp_greater(buf->ra1.ptr.p_double[k],buf->ra1.ptr.p_double[k-1]) )
                        {
                            k = k+1;
                        }
                    }
                    if( ae_fp_greater_eq(buf->ra1.ptr.p_double[t-1],buf->ra1.ptr.p_double[k-1]) )
                    {
                        t = 0;
                    }
                    else
                    {
                        tmp = buf->ra1.ptr.p_double[k-1];
                        buf->ra1.ptr.p_double[k-1] = buf->ra1.ptr.p_double[t-1];
                        buf->ra1.ptr.p_double[t-1] = tmp;
                        tmpi = buf->ia1.ptr.p_int[k-1];
                        buf->ia1.ptr.p_int[k-1] = buf->ia1.ptr.p_int[t-1];
                        buf->ia1.ptr.p_int[t-1] = tmpi;
                        t = k;
                    }
                }
            }
            i = i-1;
        }
        while(i>=1);
    }
    
    /*
     * compute tied ranks
     */
    i = 0;
    while(i<=n-1)
    {
        j = i+1;
        while(j<=n-1)
        {
            if( ae_fp_neq(buf->ra1.ptr.p_double[j],buf->ra1.ptr.p_double[i]) )
            {
                break;
            }
            j = j+1;
        }
        for(k=i; k<=j-1; k++)
        {
            buf->ra1.ptr.p_double[k] = 1+(double)(i+j-1)/(double)2;
        }
        i = j;
    }
    
    /*
     * back to x
     */
    for(i=0; i<=n-1; i++)
    {
        x->ptr.p_double[buf->ia1.ptr.p_int[i]] = buf->ra1.ptr.p_double[i];
    }
}




/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool cmatrixrank1f(ae_int_t m,
     ae_int_t n,
     /* Complex */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     /* Complex */ ae_vector* u,
     ae_int_t iu,
     /* Complex */ ae_vector* v,
     ae_int_t iv,
     ae_state *) //_state)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_cmatrixrank1f(m, n, a, ia, ja, u, iu, v, iv);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixrank1f(ae_int_t m,
     ae_int_t n,
     /* Real    */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     /* Real    */ ae_vector* u,
     ae_int_t iu,
     /* Real    */ ae_vector* v,
     ae_int_t iv,
     ae_state *) //_state)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixrank1f(m, n, a, ia, ja, u, iu, v, iv);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool cmatrixmvf(ae_int_t ,//m,
     ae_int_t ,//n,
     /* Complex */ ae_matrix* ,//a,
     ae_int_t ,//ia,
     ae_int_t ,//ja,
     ae_int_t ,//opa,
     /* Complex */ ae_vector* ,// x,
     ae_int_t ,//ix,
     /* Complex */ ae_vector* ,//y,
     ae_int_t ,//iy,
     ae_state *) //_state)
{
    ae_bool result;


    result = ae_false;
    return result;
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixmvf(ae_int_t , //m,
     ae_int_t , //n,
     /* Real    */ ae_matrix* , //a,
     ae_int_t , //ia,
     ae_int_t , //ja,
     ae_int_t , //opa,
     /* Real    */ ae_vector* , //x,
     ae_int_t , //ix,
     /* Real    */ ae_vector* , //y,
     ae_int_t , //iy,
     ae_state *) //_state)
{
    ae_bool result;


    result = ae_false;
    return result;
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool cmatrixrighttrsmf(ae_int_t m,
     ae_int_t n,
     /* Complex */ ae_matrix* a,
     ae_int_t i1,
     ae_int_t j1,
     ae_bool isupper,
     ae_bool isunit,
     ae_int_t optype,
     /* Complex */ ae_matrix* x,
     ae_int_t i2,
     ae_int_t j2,
     ae_state *) //_state)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_cmatrixrighttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool cmatrixlefttrsmf(ae_int_t m,
     ae_int_t n,
     /* Complex */ ae_matrix* a,
     ae_int_t i1,
     ae_int_t j1,
     ae_bool isupper,
     ae_bool isunit,
     ae_int_t optype,
     /* Complex */ ae_matrix* x,
     ae_int_t i2,
     ae_int_t j2,
     ae_state *) //_state)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_cmatrixlefttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixrighttrsmf(ae_int_t m,
     ae_int_t n,
     /* Real    */ ae_matrix* a,
     ae_int_t i1,
     ae_int_t j1,
     ae_bool isupper,
     ae_bool isunit,
     ae_int_t optype,
     /* Real    */ ae_matrix* x,
     ae_int_t i2,
     ae_int_t j2,
     ae_state *) //_state)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixrighttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixlefttrsmf(ae_int_t m,
     ae_int_t n,
     /* Real    */ ae_matrix* a,
     ae_int_t i1,
     ae_int_t j1,
     ae_bool isupper,
     ae_bool isunit,
     ae_int_t optype,
     /* Real    */ ae_matrix* x,
     ae_int_t i2,
     ae_int_t j2,
     ae_state *) //_state)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixlefttrsmf(m, n, a, i1, j1, isupper, isunit, optype, x, i2, j2);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool cmatrixsyrkf(ae_int_t n,
     ae_int_t k,
     double alpha,
     /* Complex */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     ae_int_t optypea,
     double beta,
     /* Complex */ ae_matrix* c,
     ae_int_t ic,
     ae_int_t jc,
     ae_bool isupper,
     ae_state *) //_state)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_cmatrixsyrkf(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixsyrkf(ae_int_t n,
     ae_int_t k,
     double alpha,
     /* Real    */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     ae_int_t optypea,
     double beta,
     /* Real    */ ae_matrix* c,
     ae_int_t ic,
     ae_int_t jc,
     ae_bool isupper,
     ae_state *) //_state)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixsyrkf(n, k, alpha, a, ia, ja, optypea, beta, c, ic, jc, isupper);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixgemmf(ae_int_t m,
     ae_int_t n,
     ae_int_t k,
     double alpha,
     /* Real    */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     ae_int_t optypea,
     /* Real    */ ae_matrix* b,
     ae_int_t ib,
     ae_int_t jb,
     ae_int_t optypeb,
     double beta,
     /* Real    */ ae_matrix* c,
     ae_int_t ic,
     ae_int_t jc,
     ae_state* ) //*_state)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_rmatrixgemmf(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
#endif
}


/*************************************************************************
Fast kernel

  -- ALGLIB routine --
     19.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool cmatrixgemmf(ae_int_t m,
     ae_int_t n,
     ae_int_t k,
     ae_complex alpha,
     /* Complex */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     ae_int_t optypea,
     /* Complex */ ae_matrix* b,
     ae_int_t ib,
     ae_int_t jb,
     ae_int_t optypeb,
     ae_complex beta,
     /* Complex */ ae_matrix* c,
     ae_int_t ic,
     ae_int_t jc,
     ae_state *) //_state)
{
#ifndef ALGLIB_INTERCEPTS_ABLAS
    ae_bool result;


    result = ae_false;
    return result;
#else
    return _ialglib_i_cmatrixgemmf(m, n, k, alpha, a, ia, ja, optypea, b, ib, jb, optypeb, beta, c, ic, jc);
#endif
}




void hermitianmatrixvectormultiply(/* Complex */ ae_matrix* a,
     ae_bool isupper,
     ae_int_t i1,
     ae_int_t i2,
     /* Complex */ ae_vector* x,
     ae_complex alpha,
     /* Complex */ ae_vector* y,
     ae_state *) //_state)
{
    ae_int_t i;
    ae_int_t ba1;
    ae_int_t ba2;
    ae_int_t by1;
    ae_int_t by2;
    ae_int_t bx1;
    ae_int_t bx2;
    ae_int_t n;
    ae_complex v;


    n = i2-i1+1;
    if( n<=0 )
    {
        return;
    }
    
    /*
     * Let A = L + D + U, where
     *  L is strictly lower triangular (main diagonal is zero)
     *  D is diagonal
     *  U is strictly upper triangular (main diagonal is zero)
     *
     * A*x = L*x + D*x + U*x
     *
     * Calculate D*x first
     */
    for(i=i1; i<=i2; i++)
    {
        y->ptr.p_complex[i-i1+1] = ae_c_mul(a->ptr.pp_complex[i][i],x->ptr.p_complex[i-i1+1]);
    }
    
    /*
     * Add L*x + U*x
     */
    if( isupper )
    {
        for(i=i1; i<=i2-1; i++)
        {
            
            /*
             * Add L*x to the result
             */
            v = x->ptr.p_complex[i-i1+1];
            by1 = i-i1+2;
            by2 = n;
            ba1 = i+1;
            ba2 = i2;
            ae_v_caddc(&y->ptr.p_complex[by1], 1, &a->ptr.pp_complex[i][ba1], 1, "Conj", ae_v_len(by1,by2), v);
            
            /*
             * Add U*x to the result
             */
            bx1 = i-i1+2;
            bx2 = n;
            ba1 = i+1;
            ba2 = i2;
            v = ae_v_cdotproduct(&x->ptr.p_complex[bx1], 1, "N", &a->ptr.pp_complex[i][ba1], 1, "N", ae_v_len(bx1,bx2));
            y->ptr.p_complex[i-i1+1] = ae_c_add(y->ptr.p_complex[i-i1+1],v);
        }
    }
    else
    {
        for(i=i1+1; i<=i2; i++)
        {
            
            /*
             * Add L*x to the result
             */
            bx1 = 1;
            bx2 = i-i1;
            ba1 = i1;
            ba2 = i-1;
            v = ae_v_cdotproduct(&x->ptr.p_complex[bx1], 1, "N", &a->ptr.pp_complex[i][ba1], 1, "N", ae_v_len(bx1,bx2));
            y->ptr.p_complex[i-i1+1] = ae_c_add(y->ptr.p_complex[i-i1+1],v);
            
            /*
             * Add U*x to the result
             */
            v = x->ptr.p_complex[i-i1+1];
            by1 = 1;
            by2 = i-i1;
            ba1 = i1;
            ba2 = i-1;
            ae_v_caddc(&y->ptr.p_complex[by1], 1, &a->ptr.pp_complex[i][ba1], 1, "Conj", ae_v_len(by1,by2), v);
        }
    }
    ae_v_cmulc(&y->ptr.p_complex[1], 1, ae_v_len(1,n), alpha);
}


void hermitianrank2update(/* Complex */ ae_matrix* a,
     ae_bool isupper,
     ae_int_t i1,
     ae_int_t i2,
     /* Complex */ ae_vector* x,
     /* Complex */ ae_vector* y,
     /* Complex */ ae_vector* t,
     ae_complex alpha,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t tp1;
    ae_int_t tp2;
    ae_complex v;


    if( isupper )
    {
        for(i=i1; i<=i2; i++)
        {
            tp1 = i+1-i1;
            tp2 = i2-i1+1;
            v = ae_c_mul(alpha,x->ptr.p_complex[i+1-i1]);
            ae_v_cmovec(&t->ptr.p_complex[tp1], 1, &y->ptr.p_complex[tp1], 1, "Conj", ae_v_len(tp1,tp2), v);
            v = ae_c_mul(ae_c_conj(alpha, _state),y->ptr.p_complex[i+1-i1]);
            ae_v_caddc(&t->ptr.p_complex[tp1], 1, &x->ptr.p_complex[tp1], 1, "Conj", ae_v_len(tp1,tp2), v);
            ae_v_cadd(&a->ptr.pp_complex[i][i], 1, &t->ptr.p_complex[tp1], 1, "N", ae_v_len(i,i2));
        }
    }
    else
    {
        for(i=i1; i<=i2; i++)
        {
            tp1 = 1;
            tp2 = i+1-i1;
            v = ae_c_mul(alpha,x->ptr.p_complex[i+1-i1]);
            ae_v_cmovec(&t->ptr.p_complex[tp1], 1, &y->ptr.p_complex[tp1], 1, "Conj", ae_v_len(tp1,tp2), v);
            v = ae_c_mul(ae_c_conj(alpha, _state),y->ptr.p_complex[i+1-i1]);
            ae_v_caddc(&t->ptr.p_complex[tp1], 1, &x->ptr.p_complex[tp1], 1, "Conj", ae_v_len(tp1,tp2), v);
            ae_v_cadd(&a->ptr.pp_complex[i][i1], 1, &t->ptr.p_complex[tp1], 1, "N", ae_v_len(i1,i));
        }
    }
}




/*************************************************************************
Generation of an elementary reflection transformation

The subroutine generates elementary reflection H of order N, so that, for
a given X, the following equality holds true:

    ( X(1) )   ( Beta )
H * (  ..  ) = (  0   )
    ( X(n) )   (  0   )

where
              ( V(1) )
H = 1 - Tau * (  ..  ) * ( V(1), ..., V(n) )
              ( V(n) )

where the first component of vector V equals 1.

Input parameters:
    X   -   vector. Array whose index ranges within [1..N].
    N   -   reflection order.

Output parameters:
    X   -   components from 2 to N are replaced with vector V.
            The first component is replaced with parameter Beta.
    Tau -   scalar value Tau. If X is a null vector, Tau equals 0,
            otherwise 1 <= Tau <= 2.

This subroutine is the modification of the DLARFG subroutines from
the LAPACK library.

MODIFICATIONS:
    24.12.2005 sign(Alpha) was replaced with an analogous to the Fortran SIGN code.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************/
void generatereflection(/* Real    */ ae_vector* x,
     ae_int_t n,
     double* tau,
     ae_state *_state)
{
    ae_int_t j;
    double alpha;
    double xnorm;
    double v;
    double beta;
    double mx;
    double s;

    *tau = 0;

    if( n<=1 )
    {
        *tau = 0;
        return;
    }
    
    /*
     * Scale if needed (to avoid overflow/underflow during intermediate
     * calculations).
     */
    mx = 0;
    for(j=1; j<=n; j++)
    {
        mx = ae_maxreal(ae_fabs(x->ptr.p_double[j], _state), mx, _state);
    }
    s = 1;
    if( ae_fp_neq(mx,0) )
    {
        if( ae_fp_less_eq(mx,ae_minrealnumber/ae_machineepsilon) )
        {
            s = ae_minrealnumber/ae_machineepsilon;
            v = 1/s;
            ae_v_muld(&x->ptr.p_double[1], 1, ae_v_len(1,n), v);
            mx = mx*v;
        }
        else
        {
            if( ae_fp_greater_eq(mx,ae_maxrealnumber*ae_machineepsilon) )
            {
                s = ae_maxrealnumber*ae_machineepsilon;
                v = 1/s;
                ae_v_muld(&x->ptr.p_double[1], 1, ae_v_len(1,n), v);
                mx = mx*v;
            }
        }
    }
    
    /*
     * XNORM = DNRM2( N-1, X, INCX )
     */
    alpha = x->ptr.p_double[1];
    xnorm = 0;
    if( ae_fp_neq(mx,0) )
    {
        for(j=2; j<=n; j++)
        {
            xnorm = xnorm+ae_sqr(x->ptr.p_double[j]/mx, _state);
        }
        xnorm = ae_sqrt(xnorm, _state)*mx;
    }
    if( ae_fp_eq(xnorm,0) )
    {
        
        /*
         * H  =  I
         */
        *tau = 0;
        x->ptr.p_double[1] = x->ptr.p_double[1]*s;
        return;
    }
    
    /*
     * general case
     */
    mx = ae_maxreal(ae_fabs(alpha, _state), ae_fabs(xnorm, _state), _state);
    beta = -mx*ae_sqrt(ae_sqr(alpha/mx, _state)+ae_sqr(xnorm/mx, _state), _state);
    if( ae_fp_less(alpha,0) )
    {
        beta = -beta;
    }
    *tau = (beta-alpha)/beta;
    v = 1/(alpha-beta);
    ae_v_muld(&x->ptr.p_double[2], 1, ae_v_len(2,n), v);
    x->ptr.p_double[1] = beta;
    
    /*
     * Scale back outputs
     */
    x->ptr.p_double[1] = x->ptr.p_double[1]*s;
}


/*************************************************************************
Application of an elementary reflection to a rectangular matrix of size MxN

The algorithm pre-multiplies the matrix by an elementary reflection transformation
which is given by column V and scalar Tau (see the description of the
GenerateReflection procedure). Not the whole matrix but only a part of it
is transformed (rows from M1 to M2, columns from N1 to N2). Only the elements
of this submatrix are changed.

Input parameters:
    C       -   matrix to be transformed.
    Tau     -   scalar defining the transformation.
    V       -   column defining the transformation.
                Array whose index ranges within [1..M2-M1+1].
    M1, M2  -   range of rows to be transformed.
    N1, N2  -   range of columns to be transformed.
    WORK    -   working array whose indexes goes from N1 to N2.

Output parameters:
    C       -   the result of multiplying the input matrix C by the
                transformation matrix which is given by Tau and V.
                If N1>N2 or M1>M2, C is not modified.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************/
void applyreflectionfromtheleft(/* Real    */ ae_matrix* c,
     double tau,
     /* Real    */ ae_vector* v,
     ae_int_t m1,
     ae_int_t m2,
     ae_int_t n1,
     ae_int_t n2,
     /* Real    */ ae_vector* work,
     ae_state *) //_state)
{
    double t;
    ae_int_t i;
    ae_int_t vm;


    if( (ae_fp_eq(tau,0)||n1>n2)||m1>m2 )
    {
        return;
    }
    
    /*
     * w := C' * v
     */
    vm = m2-m1+1;
    for(i=n1; i<=n2; i++)
    {
        work->ptr.p_double[i] = 0;
    }
    for(i=m1; i<=m2; i++)
    {
        t = v->ptr.p_double[i+1-m1];
        ae_v_addd(&work->ptr.p_double[n1], 1, &c->ptr.pp_double[i][n1], 1, ae_v_len(n1,n2), t);
    }
    
    /*
     * C := C - tau * v * w'
     */
    for(i=m1; i<=m2; i++)
    {
        t = v->ptr.p_double[i-m1+1]*tau;
        ae_v_subd(&c->ptr.pp_double[i][n1], 1, &work->ptr.p_double[n1], 1, ae_v_len(n1,n2), t);
    }
}


/*************************************************************************
Application of an elementary reflection to a rectangular matrix of size MxN

The algorithm post-multiplies the matrix by an elementary reflection transformation
which is given by column V and scalar Tau (see the description of the
GenerateReflection procedure). Not the whole matrix but only a part of it
is transformed (rows from M1 to M2, columns from N1 to N2). Only the
elements of this submatrix are changed.

Input parameters:
    C       -   matrix to be transformed.
    Tau     -   scalar defining the transformation.
    V       -   column defining the transformation.
                Array whose index ranges within [1..N2-N1+1].
    M1, M2  -   range of rows to be transformed.
    N1, N2  -   range of columns to be transformed.
    WORK    -   working array whose indexes goes from M1 to M2.

Output parameters:
    C       -   the result of multiplying the input matrix C by the
                transformation matrix which is given by Tau and V.
                If N1>N2 or M1>M2, C is not modified.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************/
void applyreflectionfromtheright(/* Real    */ ae_matrix* c,
     double tau,
     /* Real    */ ae_vector* v,
     ae_int_t m1,
     ae_int_t m2,
     ae_int_t n1,
     ae_int_t n2,
     /* Real    */ ae_vector* ,//work,
     ae_state *) //_state)
{
    double t;
    ae_int_t i;
    ae_int_t vm;


    if( (ae_fp_eq(tau,0)||n1>n2)||m1>m2 )
    {
        return;
    }
    vm = n2-n1+1;
    for(i=m1; i<=m2; i++)
    {
        t = ae_v_dotproduct(&c->ptr.pp_double[i][n1], 1, &v->ptr.p_double[1], 1, ae_v_len(n1,n2));
        t = t*tau;
        ae_v_subd(&c->ptr.pp_double[i][n1], 1, &v->ptr.p_double[1], 1, ae_v_len(n1,n2), t);
    }
}




/*************************************************************************
Generation of an elementary complex reflection transformation

The subroutine generates elementary complex reflection H of  order  N,  so
that, for a given X, the following equality holds true:

     ( X(1) )   ( Beta )
H' * (  ..  ) = (  0   ),   H'*H = I,   Beta is a real number
     ( X(n) )   (  0   )

where

              ( V(1) )
H = 1 - Tau * (  ..  ) * ( conj(V(1)), ..., conj(V(n)) )
              ( V(n) )

where the first component of vector V equals 1.

Input parameters:
    X   -   vector. Array with elements [1..N].
    N   -   reflection order.

Output parameters:
    X   -   components from 2 to N are replaced by vector V.
            The first component is replaced with parameter Beta.
    Tau -   scalar value Tau.

This subroutine is the modification of CLARFG subroutines  from the LAPACK
library. It has similar functionality except for the fact that it  doesn’t
handle errors when intermediate results cause an overflow.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************/
void complexgeneratereflection(/* Complex */ ae_vector* x,
     ae_int_t n,
     ae_complex* tau,
     ae_state *_state)
{
    ae_int_t j;
    ae_complex alpha;
    double alphi;
    double alphr;
    double beta;
    double xnorm;
    double mx;
    ae_complex t;
    double s;
    ae_complex v;

    tau->x = 0;
    tau->y = 0;

    if( n<=0 )
    {
        *tau = ae_complex_from_d(0);
        return;
    }
    
    /*
     * Scale if needed (to avoid overflow/underflow during intermediate
     * calculations).
     */
    mx = 0;
    for(j=1; j<=n; j++)
    {
        mx = ae_maxreal(ae_c_abs(x->ptr.p_complex[j], _state), mx, _state);
    }
    s = 1;
    if( ae_fp_neq(mx,0) )
    {
        if( ae_fp_less(mx,1) )
        {
            s = ae_sqrt(ae_minrealnumber, _state);
            v = ae_complex_from_d(1/s);
            ae_v_cmulc(&x->ptr.p_complex[1], 1, ae_v_len(1,n), v);
        }
        else
        {
            s = ae_sqrt(ae_maxrealnumber, _state);
            v = ae_complex_from_d(1/s);
            ae_v_cmulc(&x->ptr.p_complex[1], 1, ae_v_len(1,n), v);
        }
    }
    
    /*
     * calculate
     */
    alpha = x->ptr.p_complex[1];
    mx = 0;
    for(j=2; j<=n; j++)
    {
        mx = ae_maxreal(ae_c_abs(x->ptr.p_complex[j], _state), mx, _state);
    }
    xnorm = 0;
    if( ae_fp_neq(mx,0) )
    {
        for(j=2; j<=n; j++)
        {
            t = ae_c_div_d(x->ptr.p_complex[j],mx);
            xnorm = xnorm+ae_c_mul(t,ae_c_conj(t, _state)).x;
        }
        xnorm = ae_sqrt(xnorm, _state)*mx;
    }
    alphr = alpha.x;
    alphi = alpha.y;
    if( ae_fp_eq(xnorm,0)&&ae_fp_eq(alphi,0) )
    {
        *tau = ae_complex_from_d(0);
        x->ptr.p_complex[1] = ae_c_mul_d(x->ptr.p_complex[1],s);
        return;
    }
    mx = ae_maxreal(ae_fabs(alphr, _state), ae_fabs(alphi, _state), _state);
    mx = ae_maxreal(mx, ae_fabs(xnorm, _state), _state);
    beta = -mx*ae_sqrt(ae_sqr(alphr/mx, _state)+ae_sqr(alphi/mx, _state)+ae_sqr(xnorm/mx, _state), _state);
    if( ae_fp_less(alphr,0) )
    {
        beta = -beta;
    }
    tau->x = (beta-alphr)/beta;
    tau->y = -alphi/beta;
    alpha = ae_c_d_div(1,ae_c_sub_d(alpha,beta));
    if( n>1 )
    {
        ae_v_cmulc(&x->ptr.p_complex[2], 1, ae_v_len(2,n), alpha);
    }
    alpha = ae_complex_from_d(beta);
    x->ptr.p_complex[1] = alpha;
    
    /*
     * Scale back
     */
    x->ptr.p_complex[1] = ae_c_mul_d(x->ptr.p_complex[1],s);
}


/*************************************************************************
Application of an elementary reflection to a rectangular matrix of size MxN

The  algorithm  pre-multiplies  the  matrix  by  an  elementary reflection
transformation  which  is  given  by  column  V  and  scalar  Tau (see the
description of the GenerateReflection). Not the whole matrix  but  only  a
part of it is transformed (rows from M1 to M2, columns from N1 to N2). Only
the elements of this submatrix are changed.

Note: the matrix is multiplied by H, not by H'.   If  it  is  required  to
multiply the matrix by H', it is necessary to pass Conj(Tau) instead of Tau.

Input parameters:
    C       -   matrix to be transformed.
    Tau     -   scalar defining transformation.
    V       -   column defining transformation.
                Array whose index ranges within [1..M2-M1+1]
    M1, M2  -   range of rows to be transformed.
    N1, N2  -   range of columns to be transformed.
    WORK    -   working array whose index goes from N1 to N2.

Output parameters:
    C       -   the result of multiplying the input matrix C by the
                transformation matrix which is given by Tau and V.
                If N1>N2 or M1>M2, C is not modified.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************/
void complexapplyreflectionfromtheleft(/* Complex */ ae_matrix* c,
     ae_complex tau,
     /* Complex */ ae_vector* v,
     ae_int_t m1,
     ae_int_t m2,
     ae_int_t n1,
     ae_int_t n2,
     /* Complex */ ae_vector* work,
     ae_state *_state)
{
    ae_complex t;
    ae_int_t i;
    ae_int_t vm;


    if( (ae_c_eq_d(tau,0)||n1>n2)||m1>m2 )
    {
        return;
    }
    
    /*
     * w := C^T * conj(v)
     */
    vm = m2-m1+1;
    for(i=n1; i<=n2; i++)
    {
        work->ptr.p_complex[i] = ae_complex_from_d(0);
    }
    for(i=m1; i<=m2; i++)
    {
        t = ae_c_conj(v->ptr.p_complex[i+1-m1], _state);
        ae_v_caddc(&work->ptr.p_complex[n1], 1, &c->ptr.pp_complex[i][n1], 1, "N", ae_v_len(n1,n2), t);
    }
    
    /*
     * C := C - tau * v * w^T
     */
    for(i=m1; i<=m2; i++)
    {
        t = ae_c_mul(v->ptr.p_complex[i-m1+1],tau);
        ae_v_csubc(&c->ptr.pp_complex[i][n1], 1, &work->ptr.p_complex[n1], 1, "N", ae_v_len(n1,n2), t);
    }
}


/*************************************************************************
Application of an elementary reflection to a rectangular matrix of size MxN

The  algorithm  post-multiplies  the  matrix  by  an elementary reflection
transformation  which  is  given  by  column  V  and  scalar  Tau (see the
description  of  the  GenerateReflection). Not the whole matrix but only a
part  of  it  is  transformed (rows from M1 to M2, columns from N1 to N2).
Only the elements of this submatrix are changed.

Input parameters:
    C       -   matrix to be transformed.
    Tau     -   scalar defining transformation.
    V       -   column defining transformation.
                Array whose index ranges within [1..N2-N1+1]
    M1, M2  -   range of rows to be transformed.
    N1, N2  -   range of columns to be transformed.
    WORK    -   working array whose index goes from M1 to M2.

Output parameters:
    C       -   the result of multiplying the input matrix C by the
                transformation matrix which is given by Tau and V.
                If N1>N2 or M1>M2, C is not modified.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     September 30, 1994
*************************************************************************/
void complexapplyreflectionfromtheright(/* Complex */ ae_matrix* c,
     ae_complex tau,
     /* Complex */ ae_vector* v,
     ae_int_t m1,
     ae_int_t m2,
     ae_int_t n1,
     ae_int_t n2,
     /* Complex */ ae_vector* work,
     ae_state *) //_state)
{
    ae_complex t;
    ae_int_t i;
    ae_int_t vm;


    if( (ae_c_eq_d(tau,0)||n1>n2)||m1>m2 )
    {
        return;
    }
    
    /*
     * w := C * v
     */
    vm = n2-n1+1;
    for(i=m1; i<=m2; i++)
    {
        t = ae_v_cdotproduct(&c->ptr.pp_complex[i][n1], 1, "N", &v->ptr.p_complex[1], 1, "N", ae_v_len(n1,n2));
        work->ptr.p_complex[i] = t;
    }
    
    /*
     * C := C - w * conj(v^T)
     */
    ae_v_cmove(&v->ptr.p_complex[1], 1, &v->ptr.p_complex[1], 1, "Conj", ae_v_len(1,vm));
    for(i=m1; i<=m2; i++)
    {
        t = ae_c_mul(work->ptr.p_complex[i],tau);
        ae_v_csubc(&c->ptr.pp_complex[i][n1], 1, &v->ptr.p_complex[1], 1, "N", ae_v_len(n1,n2), t);
    }
    ae_v_cmove(&v->ptr.p_complex[1], 1, &v->ptr.p_complex[1], 1, "Conj", ae_v_len(1,vm));
}




void symmetricmatrixvectormultiply(/* Real    */ ae_matrix* a,
     ae_bool isupper,
     ae_int_t i1,
     ae_int_t i2,
     /* Real    */ ae_vector* x,
     double alpha,
     /* Real    */ ae_vector* y,
     ae_state *) //_state)
{
    ae_int_t i;
    ae_int_t ba1;
    ae_int_t ba2;
    ae_int_t by1;
    ae_int_t by2;
    ae_int_t bx1;
    ae_int_t bx2;
    ae_int_t n;
    double v;


    n = i2-i1+1;
    if( n<=0 )
    {
        return;
    }
    
    /*
     * Let A = L + D + U, where
     *  L is strictly lower triangular (main diagonal is zero)
     *  D is diagonal
     *  U is strictly upper triangular (main diagonal is zero)
     *
     * A*x = L*x + D*x + U*x
     *
     * Calculate D*x first
     */
    for(i=i1; i<=i2; i++)
    {
        y->ptr.p_double[i-i1+1] = a->ptr.pp_double[i][i]*x->ptr.p_double[i-i1+1];
    }
    
    /*
     * Add L*x + U*x
     */
    if( isupper )
    {
        for(i=i1; i<=i2-1; i++)
        {
            
            /*
             * Add L*x to the result
             */
            v = x->ptr.p_double[i-i1+1];
            by1 = i-i1+2;
            by2 = n;
            ba1 = i+1;
            ba2 = i2;
            ae_v_addd(&y->ptr.p_double[by1], 1, &a->ptr.pp_double[i][ba1], 1, ae_v_len(by1,by2), v);
            
            /*
             * Add U*x to the result
             */
            bx1 = i-i1+2;
            bx2 = n;
            ba1 = i+1;
            ba2 = i2;
            v = ae_v_dotproduct(&x->ptr.p_double[bx1], 1, &a->ptr.pp_double[i][ba1], 1, ae_v_len(bx1,bx2));
            y->ptr.p_double[i-i1+1] = y->ptr.p_double[i-i1+1]+v;
        }
    }
    else
    {
        for(i=i1+1; i<=i2; i++)
        {
            
            /*
             * Add L*x to the result
             */
            bx1 = 1;
            bx2 = i-i1;
            ba1 = i1;
            ba2 = i-1;
            v = ae_v_dotproduct(&x->ptr.p_double[bx1], 1, &a->ptr.pp_double[i][ba1], 1, ae_v_len(bx1,bx2));
            y->ptr.p_double[i-i1+1] = y->ptr.p_double[i-i1+1]+v;
            
            /*
             * Add U*x to the result
             */
            v = x->ptr.p_double[i-i1+1];
            by1 = 1;
            by2 = i-i1;
            ba1 = i1;
            ba2 = i-1;
            ae_v_addd(&y->ptr.p_double[by1], 1, &a->ptr.pp_double[i][ba1], 1, ae_v_len(by1,by2), v);
        }
    }
    ae_v_muld(&y->ptr.p_double[1], 1, ae_v_len(1,n), alpha);
}


void symmetricrank2update(/* Real    */ ae_matrix* a,
     ae_bool isupper,
     ae_int_t i1,
     ae_int_t i2,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* t,
     double alpha,
     ae_state * ) //_state)
{
    ae_int_t i;
    ae_int_t tp1;
    ae_int_t tp2;
    double v;


    if( isupper )
    {
        for(i=i1; i<=i2; i++)
        {
            tp1 = i+1-i1;
            tp2 = i2-i1+1;
            v = x->ptr.p_double[i+1-i1];
            ae_v_moved(&t->ptr.p_double[tp1], 1, &y->ptr.p_double[tp1], 1, ae_v_len(tp1,tp2), v);
            v = y->ptr.p_double[i+1-i1];
            ae_v_addd(&t->ptr.p_double[tp1], 1, &x->ptr.p_double[tp1], 1, ae_v_len(tp1,tp2), v);
            ae_v_muld(&t->ptr.p_double[tp1], 1, ae_v_len(tp1,tp2), alpha);
            ae_v_add(&a->ptr.pp_double[i][i], 1, &t->ptr.p_double[tp1], 1, ae_v_len(i,i2));
        }
    }
    else
    {
        for(i=i1; i<=i2; i++)
        {
            tp1 = 1;
            tp2 = i+1-i1;
            v = x->ptr.p_double[i+1-i1];
            ae_v_moved(&t->ptr.p_double[tp1], 1, &y->ptr.p_double[tp1], 1, ae_v_len(tp1,tp2), v);
            v = y->ptr.p_double[i+1-i1];
            ae_v_addd(&t->ptr.p_double[tp1], 1, &x->ptr.p_double[tp1], 1, ae_v_len(tp1,tp2), v);
            ae_v_muld(&t->ptr.p_double[tp1], 1, ae_v_len(tp1,tp2), alpha);
            ae_v_add(&a->ptr.pp_double[i][i1], 1, &t->ptr.p_double[tp1], 1, ae_v_len(i1,i));
        }
    }
}




double vectornorm2(/* Real    */ ae_vector* x,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t ix;
    double absxi;
    double scl;
    double ssq;
    double result;


    n = i2-i1+1;
    if( n<1 )
    {
        result = 0;
        return result;
    }
    if( n==1 )
    {
        result = ae_fabs(x->ptr.p_double[i1], _state);
        return result;
    }
    scl = 0;
    ssq = 1;
    for(ix=i1; ix<=i2; ix++)
    {
        if( ae_fp_neq(x->ptr.p_double[ix],0) )
        {
            absxi = ae_fabs(x->ptr.p_double[ix], _state);
            if( ae_fp_less(scl,absxi) )
            {
                ssq = 1+ssq*ae_sqr(scl/absxi, _state);
                scl = absxi;
            }
            else
            {
                ssq = ssq+ae_sqr(absxi/scl, _state);
            }
        }
    }
    result = scl*ae_sqrt(ssq, _state);
    return result;
}


ae_int_t vectoridxabsmax(/* Real    */ ae_vector* x,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state)
{
    ae_int_t i;
    double a;
    ae_int_t result;


    result = i1;
    a = ae_fabs(x->ptr.p_double[result], _state);
    for(i=i1+1; i<=i2; i++)
    {
        if( ae_fp_greater(ae_fabs(x->ptr.p_double[i], _state),ae_fabs(x->ptr.p_double[result], _state)) )
        {
            result = i;
        }
    }
    return result;
}


ae_int_t columnidxabsmax(/* Real    */ ae_matrix* x,
     ae_int_t i1,
     ae_int_t i2,
     ae_int_t j,
     ae_state *_state)
{
    ae_int_t i;
    double a;
    ae_int_t result;


    result = i1;
    a = ae_fabs(x->ptr.pp_double[result][j], _state);
    for(i=i1+1; i<=i2; i++)
    {
        if( ae_fp_greater(ae_fabs(x->ptr.pp_double[i][j], _state),ae_fabs(x->ptr.pp_double[result][j], _state)) )
        {
            result = i;
        }
    }
    return result;
}


ae_int_t rowidxabsmax(/* Real    */ ae_matrix* x,
     ae_int_t j1,
     ae_int_t j2,
     ae_int_t i,
     ae_state *_state)
{
    ae_int_t j;
    double a;
    ae_int_t result;


    result = j1;
    a = ae_fabs(x->ptr.pp_double[i][result], _state);
    for(j=j1+1; j<=j2; j++)
    {
        if( ae_fp_greater(ae_fabs(x->ptr.pp_double[i][j], _state),ae_fabs(x->ptr.pp_double[i][result], _state)) )
        {
            result = j;
        }
    }
    return result;
}


double upperhessenberg1norm(/* Real    */ ae_matrix* a,
     ae_int_t i1,
     ae_int_t i2,
     ae_int_t j1,
     ae_int_t j2,
     /* Real    */ ae_vector* work,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    double result;


    ae_assert(i2-i1==j2-j1, "UpperHessenberg1Norm: I2-I1<>J2-J1!", _state);
    for(j=j1; j<=j2; j++)
    {
        work->ptr.p_double[j] = 0;
    }
    for(i=i1; i<=i2; i++)
    {
        for(j=ae_maxint(j1, j1+i-i1-1, _state); j<=j2; j++)
        {
            work->ptr.p_double[j] = work->ptr.p_double[j]+ae_fabs(a->ptr.pp_double[i][j], _state);
        }
    }
    result = 0;
    for(j=j1; j<=j2; j++)
    {
        result = ae_maxreal(result, work->ptr.p_double[j], _state);
    }
    return result;
}


void copymatrix(/* Real    */ ae_matrix* a,
     ae_int_t is1,
     ae_int_t is2,
     ae_int_t js1,
     ae_int_t js2,
     /* Real    */ ae_matrix* b,
     ae_int_t id1,
     ae_int_t id2,
     ae_int_t jd1,
     ae_int_t jd2,
     ae_state *_state)
{
    ae_int_t isrc;
    ae_int_t idst;


    //if( is1>is2||js1>js2 )
    //{
    //    return;
    //}
    ae_assert(is2-is1==id2-id1, "CopyMatrix: different sizes!", _state);
    ae_assert(js2-js1==jd2-jd1, "CopyMatrix: different sizes!", _state);
    for(isrc=is1; isrc<=is2; isrc++)
    {
        idst = isrc-is1+id1;
        ae_v_move(&b->ptr.pp_double[idst][jd1], 1, &a->ptr.pp_double[isrc][js1], 1, ae_v_len(jd1,jd2));
    }
}


void inplacetranspose(/* Real    */ ae_matrix* a,
     ae_int_t i1,
     ae_int_t i2,
     ae_int_t j1,
     ae_int_t j2,
     /* Real    */ ae_vector* work,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t ips;
    ae_int_t jps;
    ae_int_t l;


    if( i1>i2||j1>j2 )
    {
        return;
    }
    ae_assert(i1-i2==j1-j2, "InplaceTranspose error: incorrect array size!", _state);
    for(i=i1; i<=i2-1; i++)
    {
        j = j1+i-i1;
        ips = i+1;
        jps = j1+ips-i1;
        l = i2-i;
        ae_v_move(&work->ptr.p_double[1], 1, &a->ptr.pp_double[ips][j], a->stride, ae_v_len(1,l));
        ae_v_move(&a->ptr.pp_double[ips][j], a->stride, &a->ptr.pp_double[i][jps], 1, ae_v_len(ips,i2));
        ae_v_move(&a->ptr.pp_double[i][jps], 1, &work->ptr.p_double[1], 1, ae_v_len(jps,j2));
    }
}


void copyandtranspose(/* Real    */ ae_matrix* a,
     ae_int_t is1,
     ae_int_t is2,
     ae_int_t js1,
     ae_int_t js2,
     /* Real    */ ae_matrix* b,
     ae_int_t id1,
     ae_int_t id2,
     ae_int_t jd1,
     ae_int_t jd2,
     ae_state *_state)
{
    ae_int_t isrc;
    ae_int_t jdst;


    if( is1>is2||js1>js2 )
    {
        return;
    }
    ae_assert(is2-is1==jd2-jd1, "CopyAndTranspose: different sizes!", _state);
    ae_assert(js2-js1==id2-id1, "CopyAndTranspose: different sizes!", _state);
    for(isrc=is1; isrc<=is2; isrc++)
    {
        jdst = isrc-is1+jd1;
        ae_v_move(&b->ptr.pp_double[id1][jdst], b->stride, &a->ptr.pp_double[isrc][js1], 1, ae_v_len(id1,id2));
    }
}


void matrixvectormultiply(/* Real    */ ae_matrix* a,
     ae_int_t i1,
     ae_int_t i2,
     ae_int_t j1,
     ae_int_t j2,
     ae_bool trans,
     /* Real    */ ae_vector* x,
     ae_int_t ix1,
     ae_int_t ix2,
     double alpha,
     /* Real    */ ae_vector* y,
     ae_int_t iy1,
     ae_int_t iy2,
     double beta,
     ae_state *_state)
{
    ae_int_t i;
    double v;


    if( !trans )
    {
        
        /*
         * y := alpha*A*x + beta*y;
         */
        if( i1>i2||j1>j2 )
        {
            return;
        }
        ae_assert(j2-j1==ix2-ix1, "MatrixVectorMultiply: A and X dont match!", _state);
        ae_assert(i2-i1==iy2-iy1, "MatrixVectorMultiply: A and Y dont match!", _state);
        
        /*
         * beta*y
         */
        if( ae_fp_eq(beta,0) )
        {
            for(i=iy1; i<=iy2; i++)
            {
                y->ptr.p_double[i] = 0;
            }
        }
        else
        {
            ae_v_muld(&y->ptr.p_double[iy1], 1, ae_v_len(iy1,iy2), beta);
        }
        
        /*
         * alpha*A*x
         */
        for(i=i1; i<=i2; i++)
        {
            v = ae_v_dotproduct(&a->ptr.pp_double[i][j1], 1, &x->ptr.p_double[ix1], 1, ae_v_len(j1,j2));
            y->ptr.p_double[iy1+i-i1] = y->ptr.p_double[iy1+i-i1]+alpha*v;
        }
    }
    else
    {
        
        /*
         * y := alpha*A'*x + beta*y;
         */
        if( i1>i2||j1>j2 )
        {
            return;
        }
        ae_assert(i2-i1==ix2-ix1, "MatrixVectorMultiply: A and X dont match!", _state);
        ae_assert(j2-j1==iy2-iy1, "MatrixVectorMultiply: A and Y dont match!", _state);
        
        /*
         * beta*y
         */
        if( ae_fp_eq(beta,0) )
        {
            for(i=iy1; i<=iy2; i++)
            {
                y->ptr.p_double[i] = 0;
            }
        }
        else
        {
            ae_v_muld(&y->ptr.p_double[iy1], 1, ae_v_len(iy1,iy2), beta);
        }
        
        /*
         * alpha*A'*x
         */
        for(i=i1; i<=i2; i++)
        {
            v = alpha*x->ptr.p_double[ix1+i-i1];
            ae_v_addd(&y->ptr.p_double[iy1], 1, &a->ptr.pp_double[i][j1], 1, ae_v_len(iy1,iy2), v);
        }
    }
}


double pythag2(double x, double y, ae_state *_state)
{
    double w;
    double xabs;
    double yabs;
    double z;
    double result;


    xabs = ae_fabs(x, _state);
    yabs = ae_fabs(y, _state);
    w = ae_maxreal(xabs, yabs, _state);
    z = ae_minreal(xabs, yabs, _state);
    if( ae_fp_eq(z,0) )
    {
        result = w;
    }
    else
    {
        result = w*ae_sqrt(1+ae_sqr(z/w, _state), _state);
    }
    return result;
}


void matrixmatrixmultiply(/* Real    */ ae_matrix* a,
     ae_int_t ai1,
     ae_int_t ai2,
     ae_int_t aj1,
     ae_int_t aj2,
     ae_bool transa,
     /* Real    */ ae_matrix* b,
     ae_int_t bi1,
     ae_int_t bi2,
     ae_int_t bj1,
     ae_int_t bj2,
     ae_bool transb,
     double alpha,
     /* Real    */ ae_matrix* c,
     ae_int_t ci1,
     ae_int_t ci2,
     ae_int_t cj1,
     ae_int_t cj2,
     double beta,
     /* Real    */ ae_vector* work,
     ae_state *_state)
{
    ae_int_t arows;
    ae_int_t acols;
    ae_int_t brows;
    ae_int_t bcols;
    ae_int_t crows;
    ae_int_t ccols;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t l;
    ae_int_t r;
    double v;


    
    /*
     * Setup
     */
    if( !transa )
    {
        arows = ai2-ai1+1;
        acols = aj2-aj1+1;
    }
    else
    {
        arows = aj2-aj1+1;
        acols = ai2-ai1+1;
    }
    if( !transb )
    {
        brows = bi2-bi1+1;
        bcols = bj2-bj1+1;
    }
    else
    {
        brows = bj2-bj1+1;
        bcols = bi2-bi1+1;
    }
    ae_assert(acols==brows, "MatrixMatrixMultiply: incorrect matrix sizes!", _state);
    if( ((arows<=0||acols<=0)||brows<=0)||bcols<=0 )
    {
        return;
    }
    crows = arows;
    ccols = bcols;
    
    /*
     * Test WORK
     */
    i = ae_maxint(arows, acols, _state);
    i = ae_maxint(brows, i, _state);
    i = ae_maxint(i, bcols, _state);
    work->ptr.p_double[1] = 0;
    work->ptr.p_double[i] = 0;
    
    /*
     * Prepare C
     */
    if( ae_fp_eq(beta,0) )
    {
        for(i=ci1; i<=ci2; i++)
        {
            for(j=cj1; j<=cj2; j++)
            {
                c->ptr.pp_double[i][j] = 0;
            }
        }
    }
    else
    {
        for(i=ci1; i<=ci2; i++)
        {
            ae_v_muld(&c->ptr.pp_double[i][cj1], 1, ae_v_len(cj1,cj2), beta);
        }
    }
    
    /*
     * A*B
     */
    if( !transa&&!transb )
    {
        for(l=ai1; l<=ai2; l++)
        {
            for(r=bi1; r<=bi2; r++)
            {
                v = alpha*a->ptr.pp_double[l][aj1+r-bi1];
                k = ci1+l-ai1;
                ae_v_addd(&c->ptr.pp_double[k][cj1], 1, &b->ptr.pp_double[r][bj1], 1, ae_v_len(cj1,cj2), v);
            }
        }
        return;
    }
    
    /*
     * A*B'
     */
    if( !transa&&transb )
    {
        if( arows*acols<brows*bcols )
        {
            for(r=bi1; r<=bi2; r++)
            {
                for(l=ai1; l<=ai2; l++)
                {
                    v = ae_v_dotproduct(&a->ptr.pp_double[l][aj1], 1, &b->ptr.pp_double[r][bj1], 1, ae_v_len(aj1,aj2));
                    c->ptr.pp_double[ci1+l-ai1][cj1+r-bi1] = c->ptr.pp_double[ci1+l-ai1][cj1+r-bi1]+alpha*v;
                }
            }
            return;
        }
        else
        {
            for(l=ai1; l<=ai2; l++)
            {
                for(r=bi1; r<=bi2; r++)
                {
                    v = ae_v_dotproduct(&a->ptr.pp_double[l][aj1], 1, &b->ptr.pp_double[r][bj1], 1, ae_v_len(aj1,aj2));
                    c->ptr.pp_double[ci1+l-ai1][cj1+r-bi1] = c->ptr.pp_double[ci1+l-ai1][cj1+r-bi1]+alpha*v;
                }
            }
            return;
        }
    }
    
    /*
     * A'*B
     */
    if( transa&&!transb )
    {
        for(l=aj1; l<=aj2; l++)
        {
            for(r=bi1; r<=bi2; r++)
            {
                v = alpha*a->ptr.pp_double[ai1+r-bi1][l];
                k = ci1+l-aj1;
                ae_v_addd(&c->ptr.pp_double[k][cj1], 1, &b->ptr.pp_double[r][bj1], 1, ae_v_len(cj1,cj2), v);
            }
        }
        return;
    }
    
    /*
     * A'*B'
     */
    if( transa&&transb )
    {
        if( arows*acols<brows*bcols )
        {
            for(r=bi1; r<=bi2; r++)
            {
                k = cj1+r-bi1;
                for(i=1; i<=crows; i++)
                {
                    work->ptr.p_double[i] = 0.0;
                }
                for(l=ai1; l<=ai2; l++)
                {
                    v = alpha*b->ptr.pp_double[r][bj1+l-ai1];
                    ae_v_addd(&work->ptr.p_double[1], 1, &a->ptr.pp_double[l][aj1], 1, ae_v_len(1,crows), v);
                }
                ae_v_add(&c->ptr.pp_double[ci1][k], c->stride, &work->ptr.p_double[1], 1, ae_v_len(ci1,ci2));
            }
            return;
        }
        else
        {
            for(l=aj1; l<=aj2; l++)
            {
                k = ai2-ai1+1;
                ae_v_move(&work->ptr.p_double[1], 1, &a->ptr.pp_double[ai1][l], a->stride, ae_v_len(1,k));
                for(r=bi1; r<=bi2; r++)
                {
                    v = ae_v_dotproduct(&work->ptr.p_double[1], 1, &b->ptr.pp_double[r][bj1], 1, ae_v_len(1,k));
                    c->ptr.pp_double[ci1+l-aj1][cj1+r-bi1] = c->ptr.pp_double[ci1+l-aj1][cj1+r-bi1]+alpha*v;
                }
            }
            return;
        }
    }
}




/*************************************************************************
Application of a sequence of  elementary rotations to a matrix

The algorithm pre-multiplies the matrix by a sequence of rotation
transformations which is given by arrays C and S. Depending on the value
of the IsForward parameter either 1 and 2, 3 and 4 and so on (if IsForward=true)
rows are rotated, or the rows N and N-1, N-2 and N-3 and so on, are rotated.

Not the whole matrix but only a part of it is transformed (rows from M1 to
M2, columns from N1 to N2). Only the elements of this submatrix are changed.

Input parameters:
    IsForward   -   the sequence of the rotation application.
    M1,M2       -   the range of rows to be transformed.
    N1, N2      -   the range of columns to be transformed.
    C,S         -   transformation coefficients.
                    Array whose index ranges within [1..M2-M1].
    A           -   processed matrix.
    WORK        -   working array whose index ranges within [N1..N2].

Output parameters:
    A           -   transformed matrix.

Utility subroutine.
*************************************************************************/
void applyrotationsfromtheleft(ae_bool isforward,
     ae_int_t m1,
     ae_int_t m2,
     ae_int_t n1,
     ae_int_t n2,
     /* Real    */ ae_vector* c,
     /* Real    */ ae_vector* s,
     /* Real    */ ae_matrix* a,
     /* Real    */ ae_vector* work,
     ae_state *) //_state)
{
    ae_int_t j;
    ae_int_t jp1;
    double ctemp;
    double stemp;
    double temp;


    if( m1>m2||n1>n2 )
    {
        return;
    }
    
    /*
     * Form  P * A
     */
    if( isforward )
    {
        if( n1!=n2 )
        {
            
            /*
             * Common case: N1<>N2
             */
            for(j=m1; j<=m2-1; j++)
            {
                ctemp = c->ptr.p_double[j-m1+1];
                stemp = s->ptr.p_double[j-m1+1];
                if( ae_fp_neq(ctemp,1)||ae_fp_neq(stemp,0) )
                {
                    jp1 = j+1;
                    ae_v_moved(&work->ptr.p_double[n1], 1, &a->ptr.pp_double[jp1][n1], 1, ae_v_len(n1,n2), ctemp);
                    ae_v_subd(&work->ptr.p_double[n1], 1, &a->ptr.pp_double[j][n1], 1, ae_v_len(n1,n2), stemp);
                    ae_v_muld(&a->ptr.pp_double[j][n1], 1, ae_v_len(n1,n2), ctemp);
                    ae_v_addd(&a->ptr.pp_double[j][n1], 1, &a->ptr.pp_double[jp1][n1], 1, ae_v_len(n1,n2), stemp);
                    ae_v_move(&a->ptr.pp_double[jp1][n1], 1, &work->ptr.p_double[n1], 1, ae_v_len(n1,n2));
                }
            }
        }
        else
        {
            
            /*
             * Special case: N1=N2
             */
            for(j=m1; j<=m2-1; j++)
            {
                ctemp = c->ptr.p_double[j-m1+1];
                stemp = s->ptr.p_double[j-m1+1];
                if( ae_fp_neq(ctemp,1)||ae_fp_neq(stemp,0) )
                {
                    temp = a->ptr.pp_double[j+1][n1];
                    a->ptr.pp_double[j+1][n1] = ctemp*temp-stemp*a->ptr.pp_double[j][n1];
                    a->ptr.pp_double[j][n1] = stemp*temp+ctemp*a->ptr.pp_double[j][n1];
                }
            }
        }
    }
    else
    {
        if( n1!=n2 )
        {
            
            /*
             * Common case: N1<>N2
             */
            for(j=m2-1; j>=m1; j--)
            {
                ctemp = c->ptr.p_double[j-m1+1];
                stemp = s->ptr.p_double[j-m1+1];
                if( ae_fp_neq(ctemp,1)||ae_fp_neq(stemp,0) )
                {
                    jp1 = j+1;
                    ae_v_moved(&work->ptr.p_double[n1], 1, &a->ptr.pp_double[jp1][n1], 1, ae_v_len(n1,n2), ctemp);
                    ae_v_subd(&work->ptr.p_double[n1], 1, &a->ptr.pp_double[j][n1], 1, ae_v_len(n1,n2), stemp);
                    ae_v_muld(&a->ptr.pp_double[j][n1], 1, ae_v_len(n1,n2), ctemp);
                    ae_v_addd(&a->ptr.pp_double[j][n1], 1, &a->ptr.pp_double[jp1][n1], 1, ae_v_len(n1,n2), stemp);
                    ae_v_move(&a->ptr.pp_double[jp1][n1], 1, &work->ptr.p_double[n1], 1, ae_v_len(n1,n2));
                }
            }
        }
        else
        {
            
            /*
             * Special case: N1=N2
             */
            for(j=m2-1; j>=m1; j--)
            {
                ctemp = c->ptr.p_double[j-m1+1];
                stemp = s->ptr.p_double[j-m1+1];
                if( ae_fp_neq(ctemp,1)||ae_fp_neq(stemp,0) )
                {
                    temp = a->ptr.pp_double[j+1][n1];
                    a->ptr.pp_double[j+1][n1] = ctemp*temp-stemp*a->ptr.pp_double[j][n1];
                    a->ptr.pp_double[j][n1] = stemp*temp+ctemp*a->ptr.pp_double[j][n1];
                }
            }
        }
    }
}


/*************************************************************************
Application of a sequence of  elementary rotations to a matrix

The algorithm post-multiplies the matrix by a sequence of rotation
transformations which is given by arrays C and S. Depending on the value
of the IsForward parameter either 1 and 2, 3 and 4 and so on (if IsForward=true)
rows are rotated, or the rows N and N-1, N-2 and N-3 and so on are rotated.

Not the whole matrix but only a part of it is transformed (rows from M1
to M2, columns from N1 to N2). Only the elements of this submatrix are changed.

Input parameters:
    IsForward   -   the sequence of the rotation application.
    M1,M2       -   the range of rows to be transformed.
    N1, N2      -   the range of columns to be transformed.
    C,S         -   transformation coefficients.
                    Array whose index ranges within [1..N2-N1].
    A           -   processed matrix.
    WORK        -   working array whose index ranges within [M1..M2].

Output parameters:
    A           -   transformed matrix.

Utility subroutine.
*************************************************************************/
void applyrotationsfromtheright(ae_bool isforward,
     ae_int_t m1,
     ae_int_t m2,
     ae_int_t n1,
     ae_int_t n2,
     /* Real    */ ae_vector* c,
     /* Real    */ ae_vector* s,
     /* Real    */ ae_matrix* a,
     /* Real    */ ae_vector* work,
     ae_state *) //_state)
{
    ae_int_t j;
    ae_int_t jp1;
    double ctemp;
    double stemp;
    double temp;


    
    /*
     * Form A * P'
     */
    if( isforward )
    {
        if( m1!=m2 )
        {
            
            /*
             * Common case: M1<>M2
             */
            for(j=n1; j<=n2-1; j++)
            {
                ctemp = c->ptr.p_double[j-n1+1];
                stemp = s->ptr.p_double[j-n1+1];
                if( ae_fp_neq(ctemp,1)||ae_fp_neq(stemp,0) )
                {
                    jp1 = j+1;
                    ae_v_moved(&work->ptr.p_double[m1], 1, &a->ptr.pp_double[m1][jp1], a->stride, ae_v_len(m1,m2), ctemp);
                    ae_v_subd(&work->ptr.p_double[m1], 1, &a->ptr.pp_double[m1][j], a->stride, ae_v_len(m1,m2), stemp);
                    ae_v_muld(&a->ptr.pp_double[m1][j], a->stride, ae_v_len(m1,m2), ctemp);
                    ae_v_addd(&a->ptr.pp_double[m1][j], a->stride, &a->ptr.pp_double[m1][jp1], a->stride, ae_v_len(m1,m2), stemp);
                    ae_v_move(&a->ptr.pp_double[m1][jp1], a->stride, &work->ptr.p_double[m1], 1, ae_v_len(m1,m2));
                }
            }
        }
        else
        {
            
            /*
             * Special case: M1=M2
             */
            for(j=n1; j<=n2-1; j++)
            {
                ctemp = c->ptr.p_double[j-n1+1];
                stemp = s->ptr.p_double[j-n1+1];
                if( ae_fp_neq(ctemp,1)||ae_fp_neq(stemp,0) )
                {
                    temp = a->ptr.pp_double[m1][j+1];
                    a->ptr.pp_double[m1][j+1] = ctemp*temp-stemp*a->ptr.pp_double[m1][j];
                    a->ptr.pp_double[m1][j] = stemp*temp+ctemp*a->ptr.pp_double[m1][j];
                }
            }
        }
    }
    else
    {
        if( m1!=m2 )
        {
            
            /*
             * Common case: M1<>M2
             */
            for(j=n2-1; j>=n1; j--)
            {
                ctemp = c->ptr.p_double[j-n1+1];
                stemp = s->ptr.p_double[j-n1+1];
                if( ae_fp_neq(ctemp,1)||ae_fp_neq(stemp,0) )
                {
                    jp1 = j+1;
                    ae_v_moved(&work->ptr.p_double[m1], 1, &a->ptr.pp_double[m1][jp1], a->stride, ae_v_len(m1,m2), ctemp);
                    ae_v_subd(&work->ptr.p_double[m1], 1, &a->ptr.pp_double[m1][j], a->stride, ae_v_len(m1,m2), stemp);
                    ae_v_muld(&a->ptr.pp_double[m1][j], a->stride, ae_v_len(m1,m2), ctemp);
                    ae_v_addd(&a->ptr.pp_double[m1][j], a->stride, &a->ptr.pp_double[m1][jp1], a->stride, ae_v_len(m1,m2), stemp);
                    ae_v_move(&a->ptr.pp_double[m1][jp1], a->stride, &work->ptr.p_double[m1], 1, ae_v_len(m1,m2));
                }
            }
        }
        else
        {
            
            /*
             * Special case: M1=M2
             */
            for(j=n2-1; j>=n1; j--)
            {
                ctemp = c->ptr.p_double[j-n1+1];
                stemp = s->ptr.p_double[j-n1+1];
                if( ae_fp_neq(ctemp,1)||ae_fp_neq(stemp,0) )
                {
                    temp = a->ptr.pp_double[m1][j+1];
                    a->ptr.pp_double[m1][j+1] = ctemp*temp-stemp*a->ptr.pp_double[m1][j];
                    a->ptr.pp_double[m1][j] = stemp*temp+ctemp*a->ptr.pp_double[m1][j];
                }
            }
        }
    }
}


/*************************************************************************
The subroutine generates the elementary rotation, so that:

[  CS  SN  ]  .  [ F ]  =  [ R ]
[ -SN  CS  ]     [ G ]     [ 0 ]

CS**2 + SN**2 = 1
*************************************************************************/
void generaterotation(double f,
     double g,
     double* cs,
     double* sn,
     double* r,
     ae_state *_state)
{
    double f1;
    double g1;

    *cs = 0;
    *sn = 0;
    *r = 0;

    if( ae_fp_eq(g,0) )
    {
        *cs = 1;
        *sn = 0;
        *r = f;
    }
    else
    {
        if( ae_fp_eq(f,0) )
        {
            *cs = 0;
            *sn = 1;
            *r = g;
        }
        else
        {
            f1 = f;
            g1 = g;
            if( ae_fp_greater(ae_fabs(f1, _state),ae_fabs(g1, _state)) )
            {
                *r = ae_fabs(f1, _state)*ae_sqrt(1+ae_sqr(g1/f1, _state), _state);
            }
            else
            {
                *r = ae_fabs(g1, _state)*ae_sqrt(1+ae_sqr(f1/g1, _state), _state);
            }
            *cs = f1/(*r);
            *sn = g1/(*r);
            if( ae_fp_greater(ae_fabs(f, _state),ae_fabs(g, _state))&&ae_fp_less(*cs,0) )
            {
                *cs = -*cs;
                *sn = -*sn;
                *r = -*r;
            }
        }
    }
}




/*************************************************************************
Subroutine performing  the  Schur  decomposition  of  a  matrix  in  upper
Hessenberg form using the QR algorithm with multiple shifts.

The  source matrix  H  is  represented as  S'*H*S = T, where H - matrix in
upper Hessenberg form,  S - orthogonal matrix (Schur vectors),   T - upper
quasi-triangular matrix (with blocks of sizes  1x1  and  2x2  on  the main
diagonal).

Input parameters:
    H   -   matrix to be decomposed.
            Array whose indexes range within [1..N, 1..N].
    N   -   size of H, N>=0.


Output parameters:
    H   –   contains the matrix T.
            Array whose indexes range within [1..N, 1..N].
            All elements below the blocks on the main diagonal are equal
            to 0.
    S   -   contains Schur vectors.
            Array whose indexes range within [1..N, 1..N].

Note 1:
    The block structure of matrix T could be easily recognized: since  all
    the elements  below  the blocks are zeros, the elements a[i+1,i] which
    are equal to 0 show the block border.

Note 2:
    the algorithm  performance  depends  on  the  value  of  the  internal
    parameter NS of InternalSchurDecomposition  subroutine  which  defines
    the number of shifts in the QR algorithm (analog of  the  block  width
    in block matrix algorithms in linear algebra). If you require  maximum
    performance  on  your  machine,  it  is  recommended  to  adjust  this
    parameter manually.

Result:
    True, if the algorithm has converged and the parameters H and S contain
        the result.
    False, if the algorithm has not converged.

Algorithm implemented on the basis of subroutine DHSEQR (LAPACK 3.0 library).
*************************************************************************/
ae_bool upperhessenbergschurdecomposition(/* Real    */ ae_matrix* h,
     ae_int_t n,
     /* Real    */ ae_matrix* s,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector wi;
    ae_vector wr;
    ae_int_t info;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_clear(s);
    ae_vector_init(&wi, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&wr, 0, DT_REAL, _state, ae_true);

    internalschurdecomposition(h, n, 1, 2, &wr, &wi, s, &info, _state);
    result = info==0;
    ae_frame_leave(_state);
    return result;
}


void internalschurdecomposition(/* Real    */ ae_matrix* h,
     ae_int_t n,
     ae_int_t tneeded,
     ae_int_t zneeded,
     /* Real    */ ae_vector* wr,
     /* Real    */ ae_vector* wi,
     /* Real    */ ae_matrix* z,
     ae_int_t* info,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector work;
    ae_int_t i;
    ae_int_t i1;
    ae_int_t i2;
    ae_int_t ierr;
    ae_int_t ii;
    ae_int_t itemp;
    ae_int_t itn;
    ae_int_t its;
    ae_int_t j;
    ae_int_t k;
    ae_int_t l;
    ae_int_t maxb;
    ae_int_t nr;
    ae_int_t ns;
    ae_int_t nv;
    double absw;
    double ovfl;
    double smlnum;
    double tau;
    double temp;
    double tst1;
    double ulp;
    double unfl;
    ae_matrix s;
    ae_vector v;
    ae_vector vv;
    ae_vector workc1;
    ae_vector works1;
    ae_vector workv3;
    ae_vector tmpwr;
    ae_vector tmpwi;
    ae_bool initz;
    ae_bool wantt;
    ae_bool wantz;
    double cnst;
    ae_bool failflag;
    ae_int_t p1;
    ae_int_t p2;
    double vt;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(wr);
    ae_vector_clear(wi);
    *info = 0;
    ae_vector_init(&work, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&s, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&v, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&vv, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&workc1, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&works1, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&workv3, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmpwr, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmpwi, 0, DT_REAL, _state, ae_true);

    
    /*
     * Set the order of the multi-shift QR algorithm to be used.
     * If you want to tune algorithm, change this values
     */
    ns = 12;
    maxb = 50;
    
    /*
     * Now 2 < NS <= MAXB < NH.
     */
    maxb = ae_maxint(3, maxb, _state);
    ns = ae_minint(maxb, ns, _state);
    
    /*
     * Initialize
     */
    cnst = 1.5;
    ae_vector_set_length(&work, ae_maxint(n, 1, _state)+1, _state);
    ae_matrix_set_length(&s, ns+1, ns+1, _state);
    ae_vector_set_length(&v, ns+1+1, _state);
    ae_vector_set_length(&vv, ns+1+1, _state);
    ae_vector_set_length(wr, ae_maxint(n, 1, _state)+1, _state);
    ae_vector_set_length(wi, ae_maxint(n, 1, _state)+1, _state);
    ae_vector_set_length(&workc1, 1+1, _state);
    ae_vector_set_length(&works1, 1+1, _state);
    ae_vector_set_length(&workv3, 3+1, _state);
    ae_vector_set_length(&tmpwr, ae_maxint(n, 1, _state)+1, _state);
    ae_vector_set_length(&tmpwi, ae_maxint(n, 1, _state)+1, _state);
    ae_assert(n>=0, "InternalSchurDecomposition: incorrect N!", _state);
    ae_assert(tneeded==0||tneeded==1, "InternalSchurDecomposition: incorrect TNeeded!", _state);
    ae_assert((zneeded==0||zneeded==1)||zneeded==2, "InternalSchurDecomposition: incorrect ZNeeded!", _state);
    wantt = tneeded==1;
    initz = zneeded==2;
    wantz = zneeded!=0;
    *info = 0;
    
    /*
     * Initialize Z, if necessary
     */
    if( initz )
    {
        ae_matrix_set_length(z, n+1, n+1, _state);
        for(i=1; i<=n; i++)
        {
            for(j=1; j<=n; j++)
            {
                if( i==j )
                {
                    z->ptr.pp_double[i][j] = 1;
                }
                else
                {
                    z->ptr.pp_double[i][j] = 0;
                }
            }
        }
    }
    
    /*
     * Quick return if possible
     */
    if( n==0 )
    {
        ae_frame_leave(_state);
        return;
    }
    if( n==1 )
    {
        wr->ptr.p_double[1] = h->ptr.pp_double[1][1];
        wi->ptr.p_double[1] = 0;
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Set rows and columns 1 to N to zero below the first
     * subdiagonal.
     */
    for(j=1; j<=n-2; j++)
    {
        for(i=j+2; i<=n; i++)
        {
            h->ptr.pp_double[i][j] = 0;
        }
    }
    
    /*
     * Test if N is sufficiently small
     */
    if( (ns<=2||ns>n)||maxb>=n )
    {
        
        /*
         * Use the standard double-shift algorithm
         */
        hsschur_internalauxschur(wantt, wantz, n, 1, n, h, wr, wi, 1, n, z, &work, &workv3, &workc1, &works1, info, _state);
        
        /*
         * fill entries under diagonal blocks of T with zeros
         */
        if( wantt )
        {
            j = 1;
            while(j<=n)
            {
                if( ae_fp_eq(wi->ptr.p_double[j],0) )
                {
                    for(i=j+1; i<=n; i++)
                    {
                        h->ptr.pp_double[i][j] = 0;
                    }
                    j = j+1;
                }
                else
                {
                    for(i=j+2; i<=n; i++)
                    {
                        h->ptr.pp_double[i][j] = 0;
                        h->ptr.pp_double[i][j+1] = 0;
                    }
                    j = j+2;
                }
            }
        }
        ae_frame_leave(_state);
        return;
    }
    unfl = ae_minrealnumber;
    ovfl = 1/unfl;
    ulp = 2*ae_machineepsilon;
    smlnum = unfl*(n/ulp);
    
    /*
     * I1 and I2 are the indices of the first row and last column of H
     * to which transformations must be applied. If eigenvalues only are
     * being computed, I1 and I2 are set inside the main loop.
     */
    i1 = 1;
    i2 = n;
    
    /*
     * ITN is the total number of multiple-shift QR iterations allowed.
     */
    itn = 30*n;
    
    /*
     * The main loop begins here. I is the loop index and decreases from
     * IHI to ILO in steps of at most MAXB. Each iteration of the loop
     * works with the active submatrix in rows and columns L to I.
     * Eigenvalues I+1 to IHI have already converged. Either L = ILO or
     * H(L,L-1) is negligible so that the matrix splits.
     */
    i = n;
    for(;;)
    {
        l = 1;
        if( i<1 )
        {
            
            /*
             * fill entries under diagonal blocks of T with zeros
             */
            if( wantt )
            {
                j = 1;
                while(j<=n)
                {
                    if( ae_fp_eq(wi->ptr.p_double[j],0) )
                    {
                        for(i=j+1; i<=n; i++)
                        {
                            h->ptr.pp_double[i][j] = 0;
                        }
                        j = j+1;
                    }
                    else
                    {
                        for(i=j+2; i<=n; i++)
                        {
                            h->ptr.pp_double[i][j] = 0;
                            h->ptr.pp_double[i][j+1] = 0;
                        }
                        j = j+2;
                    }
                }
            }
            
            /*
             * Exit
             */
            ae_frame_leave(_state);
            return;
        }
        
        /*
         * Perform multiple-shift QR iterations on rows and columns ILO to I
         * until a submatrix of order at most MAXB splits off at the bottom
         * because a subdiagonal element has become negligible.
         */
        failflag = ae_true;
        for(its=0; its<=itn; its++)
        {
            
            /*
             * Look for a single small subdiagonal element.
             */
            for(k=i; k>=l+1; k--)
            {
                tst1 = ae_fabs(h->ptr.pp_double[k-1][k-1], _state)+ae_fabs(h->ptr.pp_double[k][k], _state);
                if( ae_fp_eq(tst1,0) )
                {
                    tst1 = upperhessenberg1norm(h, l, i, l, i, &work, _state);
                }
                if( ae_fp_less_eq(ae_fabs(h->ptr.pp_double[k][k-1], _state),ae_maxreal(ulp*tst1, smlnum, _state)) )
                {
                    break;
                }
            }
            l = k;
            if( l>1 )
            {
                
                /*
                 * H(L,L-1) is negligible.
                 */
                h->ptr.pp_double[l][l-1] = 0;
            }
            
            /*
             * Exit from loop if a submatrix of order <= MAXB has split off.
             */
            if( l>=i-maxb+1 )
            {
                failflag = ae_false;
                break;
            }
            
            /*
             * Now the active submatrix is in rows and columns L to I. If
             * eigenvalues only are being computed, only the active submatrix
             * need be transformed.
             */
            if( its==20||its==30 )
            {
                
                /*
                 * Exceptional shifts.
                 */
                for(ii=i-ns+1; ii<=i; ii++)
                {
                    wr->ptr.p_double[ii] = cnst*(ae_fabs(h->ptr.pp_double[ii][ii-1], _state)+ae_fabs(h->ptr.pp_double[ii][ii], _state));
                    wi->ptr.p_double[ii] = 0;
                }
            }
            else
            {
                
                /*
                 * Use eigenvalues of trailing submatrix of order NS as shifts.
                 */
                copymatrix(h, i-ns+1, i, i-ns+1, i, &s, 1, ns, 1, ns, _state);
                hsschur_internalauxschur(ae_false, ae_false, ns, 1, ns, &s, &tmpwr, &tmpwi, 1, ns, z, &work, &workv3, &workc1, &works1, &ierr, _state);
                for(p1=1; p1<=ns; p1++)
                {
                    wr->ptr.p_double[i-ns+p1] = tmpwr.ptr.p_double[p1];
                    wi->ptr.p_double[i-ns+p1] = tmpwi.ptr.p_double[p1];
                }
                if( ierr>0 )
                {
                    
                    /*
                     * If DLAHQR failed to compute all NS eigenvalues, use the
                     * unconverged diagonal elements as the remaining shifts.
                     */
                    for(ii=1; ii<=ierr; ii++)
                    {
                        wr->ptr.p_double[i-ns+ii] = s.ptr.pp_double[ii][ii];
                        wi->ptr.p_double[i-ns+ii] = 0;
                    }
                }
            }
            
            /*
             * Form the first column of (G-w(1)) (G-w(2)) . . . (G-w(ns))
             * where G is the Hessenberg submatrix H(L:I,L:I) and w is
             * the vector of shifts (stored in WR and WI). The result is
             * stored in the local array V.
             */
            v.ptr.p_double[1] = 1;
            for(ii=2; ii<=ns+1; ii++)
            {
                v.ptr.p_double[ii] = 0;
            }
            nv = 1;
            for(j=i-ns+1; j<=i; j++)
            {
                if( ae_fp_greater_eq(wi->ptr.p_double[j],0) )
                {
                    if( ae_fp_eq(wi->ptr.p_double[j],0) )
                    {
                        
                        /*
                         * real shift
                         */
                        p1 = nv+1;
                        ae_v_move(&vv.ptr.p_double[1], 1, &v.ptr.p_double[1], 1, ae_v_len(1,p1));
                        matrixvectormultiply(h, l, l+nv, l, l+nv-1, ae_false, &vv, 1, nv, 1.0, &v, 1, nv+1, -wr->ptr.p_double[j], _state);
                        nv = nv+1;
                    }
                    else
                    {
                        if( ae_fp_greater(wi->ptr.p_double[j],0) )
                        {
                            
                            /*
                             * complex conjugate pair of shifts
                             */
                            p1 = nv+1;
                            ae_v_move(&vv.ptr.p_double[1], 1, &v.ptr.p_double[1], 1, ae_v_len(1,p1));
                            matrixvectormultiply(h, l, l+nv, l, l+nv-1, ae_false, &v, 1, nv, 1.0, &vv, 1, nv+1, -2*wr->ptr.p_double[j], _state);
                            itemp = vectoridxabsmax(&vv, 1, nv+1, _state);
                            temp = 1/ae_maxreal(ae_fabs(vv.ptr.p_double[itemp], _state), smlnum, _state);
                            p1 = nv+1;
                            ae_v_muld(&vv.ptr.p_double[1], 1, ae_v_len(1,p1), temp);
                            absw = pythag2(wr->ptr.p_double[j], wi->ptr.p_double[j], _state);
                            temp = temp*absw*absw;
                            matrixvectormultiply(h, l, l+nv+1, l, l+nv, ae_false, &vv, 1, nv+1, 1.0, &v, 1, nv+2, temp, _state);
                            nv = nv+2;
                        }
                    }
                    
                    /*
                     * Scale V(1:NV) so that max(abs(V(i))) = 1. If V is zero,
                     * reset it to the unit vector.
                     */
                    itemp = vectoridxabsmax(&v, 1, nv, _state);
                    temp = ae_fabs(v.ptr.p_double[itemp], _state);
                    if( ae_fp_eq(temp,0) )
                    {
                        v.ptr.p_double[1] = 1;
                        for(ii=2; ii<=nv; ii++)
                        {
                            v.ptr.p_double[ii] = 0;
                        }
                    }
                    else
                    {
                        temp = ae_maxreal(temp, smlnum, _state);
                        vt = 1/temp;
                        ae_v_muld(&v.ptr.p_double[1], 1, ae_v_len(1,nv), vt);
                    }
                }
            }
            
            /*
             * Multiple-shift QR step
             */
            for(k=l; k<=i-1; k++)
            {
                
                /*
                 * The first iteration of this loop determines a reflection G
                 * from the vector V and applies it from left and right to H,
                 * thus creating a nonzero bulge below the subdiagonal.
                 *
                 * Each subsequent iteration determines a reflection G to
                 * restore the Hessenberg form in the (K-1)th column, and thus
                 * chases the bulge one step toward the bottom of the active
                 * submatrix. NR is the order of G.
                 */
                nr = ae_minint(ns+1, i-k+1, _state);
                if( k>l )
                {
                    p1 = k-1;
                    p2 = k+nr-1;
                    ae_v_move(&v.ptr.p_double[1], 1, &h->ptr.pp_double[k][p1], h->stride, ae_v_len(1,nr));
                }
                generatereflection(&v, nr, &tau, _state);
                if( k>l )
                {
                    h->ptr.pp_double[k][k-1] = v.ptr.p_double[1];
                    for(ii=k+1; ii<=i; ii++)
                    {
                        h->ptr.pp_double[ii][k-1] = 0;
                    }
                }
                v.ptr.p_double[1] = 1;
                
                /*
                 * Apply G from the left to transform the rows of the matrix in
                 * columns K to I2.
                 */
                applyreflectionfromtheleft(h, tau, &v, k, k+nr-1, k, i2, &work, _state);
                
                /*
                 * Apply G from the right to transform the columns of the
                 * matrix in rows I1 to min(K+NR,I).
                 */
                applyreflectionfromtheright(h, tau, &v, i1, ae_minint(k+nr, i, _state), k, k+nr-1, &work, _state);
                if( wantz )
                {
                    
                    /*
                     * Accumulate transformations in the matrix Z
                     */
                    applyreflectionfromtheright(z, tau, &v, 1, n, k, k+nr-1, &work, _state);
                }
            }
        }
        
        /*
         * Failure to converge in remaining number of iterations
         */
        if( failflag )
        {
            *info = i;
            ae_frame_leave(_state);
            return;
        }
        
        /*
         * A submatrix of order <= MAXB in rows and columns L to I has split
         * off. Use the double-shift QR algorithm to handle it.
         */
        hsschur_internalauxschur(wantt, wantz, n, l, i, h, wr, wi, 1, n, z, &work, &workv3, &workc1, &works1, info, _state);
        if( *info>0 )
        {
            ae_frame_leave(_state);
            return;
        }
        
        /*
         * Decrement number of remaining iterations, and return to start of
         * the main loop with a new value of I.
         */
        itn = itn-its;
        i = l-1;
    }
    ae_frame_leave(_state);
}


static void hsschur_internalauxschur(ae_bool wantt,
     ae_bool wantz,
     ae_int_t n,
     ae_int_t ilo,
     ae_int_t ihi,
     /* Real    */ ae_matrix* h,
     /* Real    */ ae_vector* wr,
     /* Real    */ ae_vector* wi,
     ae_int_t iloz,
     ae_int_t ihiz,
     /* Real    */ ae_matrix* z,
     /* Real    */ ae_vector* work,
     /* Real    */ ae_vector* workv3,
     /* Real    */ ae_vector* workc1,
     /* Real    */ ae_vector* works1,
     ae_int_t* info,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t i1;
    ae_int_t i2;
    ae_int_t itn;
    ae_int_t its;
    ae_int_t j;
    ae_int_t k;
    ae_int_t l;
    ae_int_t m;
    ae_int_t nh;
    ae_int_t nr;
    ae_int_t nz;
    double ave;
    double cs;
    double disc;
    double h00;
    double h10;
    double h11;
    double h12;
    double h21;
    double h22;
    double h33;
    double h33s;
    double h43h34;
    double h44;
    double h44s;
    double ovfl;
    double s;
    double smlnum;
    double sn;
    double sum;
    double t1;
    double t2;
    double t3;
    double tst1;
    double unfl;
    double v1;
    double v2;
    double v3;
    ae_bool failflag;
    double dat1;
    double dat2;
    ae_int_t p1;
    double him1im1;
    double him1i;
    double hiim1;
    double hii;
    double wrim1;
    double wri;
    double wiim1;
    double wii;
    double ulp;

    *info = 0;

    *info = 0;
    dat1 = 0.75;
    dat2 = -0.4375;
    ulp = ae_machineepsilon;
    
    /*
     * Quick return if possible
     */
    if( n==0 )
    {
        return;
    }
    if( ilo==ihi )
    {
        wr->ptr.p_double[ilo] = h->ptr.pp_double[ilo][ilo];
        wi->ptr.p_double[ilo] = 0;
        return;
    }
    nh = ihi-ilo+1;
    nz = ihiz-iloz+1;
    
    /*
     * Set machine-dependent constants for the stopping criterion.
     * If norm(H) <= sqrt(OVFL), overflow should not occur.
     */
    unfl = ae_minrealnumber;
    ovfl = 1/unfl;
    smlnum = unfl*(nh/ulp);
    
    /*
     * I1 and I2 are the indices of the first row and last column of H
     * to which transformations must be applied. If eigenvalues only are
     * being computed, I1 and I2 are set inside the main loop.
     */
    i1 = 1;
    i2 = n;
    
    /*
     * ITN is the total number of QR iterations allowed.
     */
    itn = 30*nh;
    
    /*
     * The main loop begins here. I is the loop index and decreases from
     * IHI to ILO in steps of 1 or 2. Each iteration of the loop works
     * with the active submatrix in rows and columns L to I.
     * Eigenvalues I+1 to IHI have already converged. Either L = ILO or
     * H(L,L-1) is negligible so that the matrix splits.
     */
    i = ihi;
    for(;;)
    {
        l = ilo;
        if( i<ilo )
        {
            return;
        }
        
        /*
         * Perform QR iterations on rows and columns ILO to I until a
         * submatrix of order 1 or 2 splits off at the bottom because a
         * subdiagonal element has become negligible.
         */
        failflag = ae_true;
        for(its=0; its<=itn; its++)
        {
            
            /*
             * Look for a single small subdiagonal element.
             */
            for(k=i; k>=l+1; k--)
            {
                tst1 = ae_fabs(h->ptr.pp_double[k-1][k-1], _state)+ae_fabs(h->ptr.pp_double[k][k], _state);
                if( ae_fp_eq(tst1,0) )
                {
                    tst1 = upperhessenberg1norm(h, l, i, l, i, work, _state);
                }
                if( ae_fp_less_eq(ae_fabs(h->ptr.pp_double[k][k-1], _state),ae_maxreal(ulp*tst1, smlnum, _state)) )
                {
                    break;
                }
            }
            l = k;
            if( l>ilo )
            {
                
                /*
                 * H(L,L-1) is negligible
                 */
                h->ptr.pp_double[l][l-1] = 0;
            }
            
            /*
             * Exit from loop if a submatrix of order 1 or 2 has split off.
             */
            if( l>=i-1 )
            {
                failflag = ae_false;
                break;
            }
            
            /*
             * Now the active submatrix is in rows and columns L to I. If
             * eigenvalues only are being computed, only the active submatrix
             * need be transformed.
             */
            if( its==10||its==20 )
            {
                
                /*
                 * Exceptional shift.
                 */
                s = ae_fabs(h->ptr.pp_double[i][i-1], _state)+ae_fabs(h->ptr.pp_double[i-1][i-2], _state);
                h44 = dat1*s+h->ptr.pp_double[i][i];
                h33 = h44;
                h43h34 = dat2*s*s;
            }
            else
            {
                
                /*
                 * Prepare to use Francis' double shift
                 * (i.e. 2nd degree generalized Rayleigh quotient)
                 */
                h44 = h->ptr.pp_double[i][i];
                h33 = h->ptr.pp_double[i-1][i-1];
                h43h34 = h->ptr.pp_double[i][i-1]*h->ptr.pp_double[i-1][i];
                s = h->ptr.pp_double[i-1][i-2]*h->ptr.pp_double[i-1][i-2];
                disc = (h33-h44)*0.5;
                disc = disc*disc+h43h34;
                if( ae_fp_greater(disc,0) )
                {
                    
                    /*
                     * Real roots: use Wilkinson's shift twice
                     */
                    disc = ae_sqrt(disc, _state);
                    ave = 0.5*(h33+h44);
                    if( ae_fp_greater(ae_fabs(h33, _state)-ae_fabs(h44, _state),0) )
                    {
                        h33 = h33*h44-h43h34;
                        h44 = h33/(hsschur_extschursign(disc, ave, _state)+ave);
                    }
                    else
                    {
                        h44 = hsschur_extschursign(disc, ave, _state)+ave;
                    }
                    h33 = h44;
                    h43h34 = 0;
                }
            }
            
            /*
             * Look for two consecutive small subdiagonal elements.
             */
            for(m=i-2; m>=l; m--)
            {
                
                /*
                 * Determine the effect of starting the double-shift QR
                 * iteration at row M, and see if this would make H(M,M-1)
                 * negligible.
                 */
                h11 = h->ptr.pp_double[m][m];
                h22 = h->ptr.pp_double[m+1][m+1];
                h21 = h->ptr.pp_double[m+1][m];
                h12 = h->ptr.pp_double[m][m+1];
                h44s = h44-h11;
                h33s = h33-h11;
                v1 = (h33s*h44s-h43h34)/h21+h12;
                v2 = h22-h11-h33s-h44s;
                v3 = h->ptr.pp_double[m+2][m+1];
                s = ae_fabs(v1, _state)+ae_fabs(v2, _state)+ae_fabs(v3, _state);
                v1 = v1/s;
                v2 = v2/s;
                v3 = v3/s;
                workv3->ptr.p_double[1] = v1;
                workv3->ptr.p_double[2] = v2;
                workv3->ptr.p_double[3] = v3;
                if( m==l )
                {
                    break;
                }
                h00 = h->ptr.pp_double[m-1][m-1];
                h10 = h->ptr.pp_double[m][m-1];
                tst1 = ae_fabs(v1, _state)*(ae_fabs(h00, _state)+ae_fabs(h11, _state)+ae_fabs(h22, _state));
                if( ae_fp_less_eq(ae_fabs(h10, _state)*(ae_fabs(v2, _state)+ae_fabs(v3, _state)),ulp*tst1) )
                {
                    break;
                }
            }
            
            /*
             * Double-shift QR step
             */
            for(k=m; k<=i-1; k++)
            {
                
                /*
                 * The first iteration of this loop determines a reflection G
                 * from the vector V and applies it from left and right to H,
                 * thus creating a nonzero bulge below the subdiagonal.
                 *
                 * Each subsequent iteration determines a reflection G to
                 * restore the Hessenberg form in the (K-1)th column, and thus
                 * chases the bulge one step toward the bottom of the active
                 * submatrix. NR is the order of G.
                 */
                nr = ae_minint(3, i-k+1, _state);
                if( k>m )
                {
                    for(p1=1; p1<=nr; p1++)
                    {
                        workv3->ptr.p_double[p1] = h->ptr.pp_double[k+p1-1][k-1];
                    }
                }
                generatereflection(workv3, nr, &t1, _state);
                if( k>m )
                {
                    h->ptr.pp_double[k][k-1] = workv3->ptr.p_double[1];
                    h->ptr.pp_double[k+1][k-1] = 0;
                    if( k<i-1 )
                    {
                        h->ptr.pp_double[k+2][k-1] = 0;
                    }
                }
                else
                {
                    if( m>l )
                    {
                        h->ptr.pp_double[k][k-1] = -h->ptr.pp_double[k][k-1];
                    }
                }
                v2 = workv3->ptr.p_double[2];
                t2 = t1*v2;
                if( nr==3 )
                {
                    v3 = workv3->ptr.p_double[3];
                    t3 = t1*v3;
                    
                    /*
                     * Apply G from the left to transform the rows of the matrix
                     * in columns K to I2.
                     */
                    for(j=k; j<=i2; j++)
                    {
                        sum = h->ptr.pp_double[k][j]+v2*h->ptr.pp_double[k+1][j]+v3*h->ptr.pp_double[k+2][j];
                        h->ptr.pp_double[k][j] = h->ptr.pp_double[k][j]-sum*t1;
                        h->ptr.pp_double[k+1][j] = h->ptr.pp_double[k+1][j]-sum*t2;
                        h->ptr.pp_double[k+2][j] = h->ptr.pp_double[k+2][j]-sum*t3;
                    }
                    
                    /*
                     * Apply G from the right to transform the columns of the
                     * matrix in rows I1 to min(K+3,I).
                     */
                    for(j=i1; j<=ae_minint(k+3, i, _state); j++)
                    {
                        sum = h->ptr.pp_double[j][k]+v2*h->ptr.pp_double[j][k+1]+v3*h->ptr.pp_double[j][k+2];
                        h->ptr.pp_double[j][k] = h->ptr.pp_double[j][k]-sum*t1;
                        h->ptr.pp_double[j][k+1] = h->ptr.pp_double[j][k+1]-sum*t2;
                        h->ptr.pp_double[j][k+2] = h->ptr.pp_double[j][k+2]-sum*t3;
                    }
                    if( wantz )
                    {
                        
                        /*
                         * Accumulate transformations in the matrix Z
                         */
                        for(j=iloz; j<=ihiz; j++)
                        {
                            sum = z->ptr.pp_double[j][k]+v2*z->ptr.pp_double[j][k+1]+v3*z->ptr.pp_double[j][k+2];
                            z->ptr.pp_double[j][k] = z->ptr.pp_double[j][k]-sum*t1;
                            z->ptr.pp_double[j][k+1] = z->ptr.pp_double[j][k+1]-sum*t2;
                            z->ptr.pp_double[j][k+2] = z->ptr.pp_double[j][k+2]-sum*t3;
                        }
                    }
                }
                else
                {
                    if( nr==2 )
                    {
                        
                        /*
                         * Apply G from the left to transform the rows of the matrix
                         * in columns K to I2.
                         */
                        for(j=k; j<=i2; j++)
                        {
                            sum = h->ptr.pp_double[k][j]+v2*h->ptr.pp_double[k+1][j];
                            h->ptr.pp_double[k][j] = h->ptr.pp_double[k][j]-sum*t1;
                            h->ptr.pp_double[k+1][j] = h->ptr.pp_double[k+1][j]-sum*t2;
                        }
                        
                        /*
                         * Apply G from the right to transform the columns of the
                         * matrix in rows I1 to min(K+3,I).
                         */
                        for(j=i1; j<=i; j++)
                        {
                            sum = h->ptr.pp_double[j][k]+v2*h->ptr.pp_double[j][k+1];
                            h->ptr.pp_double[j][k] = h->ptr.pp_double[j][k]-sum*t1;
                            h->ptr.pp_double[j][k+1] = h->ptr.pp_double[j][k+1]-sum*t2;
                        }
                        if( wantz )
                        {
                            
                            /*
                             * Accumulate transformations in the matrix Z
                             */
                            for(j=iloz; j<=ihiz; j++)
                            {
                                sum = z->ptr.pp_double[j][k]+v2*z->ptr.pp_double[j][k+1];
                                z->ptr.pp_double[j][k] = z->ptr.pp_double[j][k]-sum*t1;
                                z->ptr.pp_double[j][k+1] = z->ptr.pp_double[j][k+1]-sum*t2;
                            }
                        }
                    }
                }
            }
        }
        if( failflag )
        {
            
            /*
             * Failure to converge in remaining number of iterations
             */
            *info = i;
            return;
        }
        if( l==i )
        {
            
            /*
             * H(I,I-1) is negligible: one eigenvalue has converged.
             */
            wr->ptr.p_double[i] = h->ptr.pp_double[i][i];
            wi->ptr.p_double[i] = 0;
        }
        else
        {
            if( l==i-1 )
            {
                
                /*
                 * H(I-1,I-2) is negligible: a pair of eigenvalues have converged.
                 *
                 *        Transform the 2-by-2 submatrix to standard Schur form,
                 *        and compute and store the eigenvalues.
                 */
                him1im1 = h->ptr.pp_double[i-1][i-1];
                him1i = h->ptr.pp_double[i-1][i];
                hiim1 = h->ptr.pp_double[i][i-1];
                hii = h->ptr.pp_double[i][i];
                hsschur_aux2x2schur(&him1im1, &him1i, &hiim1, &hii, &wrim1, &wiim1, &wri, &wii, &cs, &sn, _state);
                wr->ptr.p_double[i-1] = wrim1;
                wi->ptr.p_double[i-1] = wiim1;
                wr->ptr.p_double[i] = wri;
                wi->ptr.p_double[i] = wii;
                h->ptr.pp_double[i-1][i-1] = him1im1;
                h->ptr.pp_double[i-1][i] = him1i;
                h->ptr.pp_double[i][i-1] = hiim1;
                h->ptr.pp_double[i][i] = hii;
                if( wantt )
                {
                    
                    /*
                     * Apply the transformation to the rest of H.
                     */
                    if( i2>i )
                    {
                        workc1->ptr.p_double[1] = cs;
                        works1->ptr.p_double[1] = sn;
                        applyrotationsfromtheleft(ae_true, i-1, i, i+1, i2, workc1, works1, h, work, _state);
                    }
                    workc1->ptr.p_double[1] = cs;
                    works1->ptr.p_double[1] = sn;
                    applyrotationsfromtheright(ae_true, i1, i-2, i-1, i, workc1, works1, h, work, _state);
                }
                if( wantz )
                {
                    
                    /*
                     * Apply the transformation to Z.
                     */
                    workc1->ptr.p_double[1] = cs;
                    works1->ptr.p_double[1] = sn;
                    applyrotationsfromtheright(ae_true, iloz, iloz+nz-1, i-1, i, workc1, works1, z, work, _state);
                }
            }
        }
        
        /*
         * Decrement number of remaining iterations, and return to start of
         * the main loop with new value of I.
         */
        itn = itn-its;
        i = l-1;
    }
}


static void hsschur_aux2x2schur(double* a,
     double* b,
     double* c,
     double* d,
     double* rt1r,
     double* rt1i,
     double* rt2r,
     double* rt2i,
     double* cs,
     double* sn,
     ae_state *_state)
{
    double multpl;
    double aa;
    double bb;
    double bcmax;
    double bcmis;
    double cc;
    double cs1;
    double dd;
    double eps;
    double p;
    double sab;
    double sac;
    double scl;
    double sigma;
    double sn1;
    double tau;
    double temp;
    double z;

    *rt1r = 0;
    *rt1i = 0;
    *rt2r = 0;
    *rt2i = 0;
    *cs = 0;
    *sn = 0;

    multpl = 4.0;
    eps = ae_machineepsilon;
    if( ae_fp_eq(*c,0) )
    {
        *cs = 1;
        *sn = 0;
    }
    else
    {
        if( ae_fp_eq(*b,0) )
        {
            
            /*
             * Swap rows and columns
             */
            *cs = 0;
            *sn = 1;
            temp = *d;
            *d = *a;
            *a = temp;
            *b = -*c;
            *c = 0;
        }
        else
        {
            if( ae_fp_eq(*a-(*d),0)&&hsschur_extschursigntoone(*b, _state)!=hsschur_extschursigntoone(*c, _state) )
            {
                *cs = 1;
                *sn = 0;
            }
            else
            {
                temp = *a-(*d);
                p = 0.5*temp;
                bcmax = ae_maxreal(ae_fabs(*b, _state), ae_fabs(*c, _state), _state);
                bcmis = ae_minreal(ae_fabs(*b, _state), ae_fabs(*c, _state), _state)*hsschur_extschursigntoone(*b, _state)*hsschur_extschursigntoone(*c, _state);
                scl = ae_maxreal(ae_fabs(p, _state), bcmax, _state);
                z = p/scl*p+bcmax/scl*bcmis;
                
                /*
                 * If Z is of the order of the machine accuracy, postpone the
                 * decision on the nature of eigenvalues
                 */
                if( ae_fp_greater_eq(z,multpl*eps) )
                {
                    
                    /*
                     * Real eigenvalues. Compute A and D.
                     */
                    z = p+hsschur_extschursign(ae_sqrt(scl, _state)*ae_sqrt(z, _state), p, _state);
                    *a = *d+z;
                    *d = *d-bcmax/z*bcmis;
                    
                    /*
                     * Compute B and the rotation matrix
                     */
                    tau = pythag2(*c, z, _state);
                    *cs = z/tau;
                    *sn = *c/tau;
                    *b = *b-(*c);
                    *c = 0;
                }
                else
                {
                    
                    /*
                     * Complex eigenvalues, or real (almost) equal eigenvalues.
                     * Make diagonal elements equal.
                     */
                    sigma = *b+(*c);
                    tau = pythag2(sigma, temp, _state);
                    *cs = ae_sqrt(0.5*(1+ae_fabs(sigma, _state)/tau), _state);
                    *sn = -p/(tau*(*cs))*hsschur_extschursign(1, sigma, _state);
                    
                    /*
                     * Compute [ AA  BB ] = [ A  B ] [ CS -SN ]
                     *         [ CC  DD ]   [ C  D ] [ SN  CS ]
                     */
                    aa = *a*(*cs)+*b*(*sn);
                    bb = -*a*(*sn)+*b*(*cs);
                    cc = *c*(*cs)+*d*(*sn);
                    dd = -*c*(*sn)+*d*(*cs);
                    
                    /*
                     * Compute [ A  B ] = [ CS  SN ] [ AA  BB ]
                     *         [ C  D ]   [-SN  CS ] [ CC  DD ]
                     */
                    *a = aa*(*cs)+cc*(*sn);
                    *b = bb*(*cs)+dd*(*sn);
                    *c = -aa*(*sn)+cc*(*cs);
                    *d = -bb*(*sn)+dd*(*cs);
                    temp = 0.5*(*a+(*d));
                    *a = temp;
                    *d = temp;
                    if( ae_fp_neq(*c,0) )
                    {
                        if( ae_fp_neq(*b,0) )
                        {
                            if( hsschur_extschursigntoone(*b, _state)==hsschur_extschursigntoone(*c, _state) )
                            {
                                
                                /*
                                 * Real eigenvalues: reduce to upper triangular form
                                 */
                                sab = ae_sqrt(ae_fabs(*b, _state), _state);
                                sac = ae_sqrt(ae_fabs(*c, _state), _state);
                                p = hsschur_extschursign(sab*sac, *c, _state);
                                tau = 1/ae_sqrt(ae_fabs(*b+(*c), _state), _state);
                                *a = temp+p;
                                *d = temp-p;
                                *b = *b-(*c);
                                *c = 0;
                                cs1 = sab*tau;
                                sn1 = sac*tau;
                                temp = *cs*cs1-*sn*sn1;
                                *sn = *cs*sn1+*sn*cs1;
                                *cs = temp;
                            }
                        }
                        else
                        {
                            *b = -*c;
                            *c = 0;
                            temp = *cs;
                            *cs = -*sn;
                            *sn = temp;
                        }
                    }
                }
            }
        }
    }
    
    /*
     * Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).
     */
    *rt1r = *a;
    *rt2r = *d;
    if( ae_fp_eq(*c,0) )
    {
        *rt1i = 0;
        *rt2i = 0;
    }
    else
    {
        *rt1i = ae_sqrt(ae_fabs(*b, _state), _state)*ae_sqrt(ae_fabs(*c, _state), _state);
        *rt2i = -*rt1i;
    }
}


static double hsschur_extschursign(double a, double b, ae_state *_state)
{
    double result;


    if( ae_fp_greater_eq(b,0) )
    {
        result = ae_fabs(a, _state);
    }
    else
    {
        result = -ae_fabs(a, _state);
    }
    return result;
}


static ae_int_t hsschur_extschursigntoone(double b, ae_state *) //_state)
{
    ae_int_t result;


    if( ae_fp_greater_eq(b,0) )
    {
        result = 1;
    }
    else
    {
        result = -1;
    }
    return result;
}




/*************************************************************************
Utility subroutine performing the "safe" solution of system of linear
equations with triangular coefficient matrices.

The subroutine uses scaling and solves the scaled system A*x=s*b (where  s
is  a  scalar  value)  instead  of  A*x=b,  choosing  s  so  that x can be
represented by a floating-point number. The closer the system  gets  to  a
singular, the less s is. If the system is singular, s=0 and x contains the
non-trivial solution of equation A*x=0.

The feature of an algorithm is that it could not cause an  overflow  or  a
division by zero regardless of the matrix used as the input.

The algorithm can solve systems of equations with  upper/lower  triangular
matrices,  with/without unit diagonal, and systems of type A*x=b or A'*x=b
(where A' is a transposed matrix A).

Input parameters:
    A       -   system matrix. Array whose indexes range within [0..N-1, 0..N-1].
    N       -   size of matrix A.
    X       -   right-hand member of a system.
                Array whose index ranges within [0..N-1].
    IsUpper -   matrix type. If it is True, the system matrix is the upper
                triangular and is located in  the  corresponding  part  of
                matrix A.
    Trans   -   problem type. If it is True, the problem to be  solved  is
                A'*x=b, otherwise it is A*x=b.
    Isunit  -   matrix type. If it is True, the system matrix has  a  unit
                diagonal (the elements on the main diagonal are  not  used
                in the calculation process), otherwise the matrix is considered
                to be a general triangular matrix.

Output parameters:
    X       -   solution. Array whose index ranges within [0..N-1].
    S       -   scaling factor.

  -- LAPACK auxiliary routine (version 3.0) --
     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
     Courant Institute, Argonne National Lab, and Rice University
     June 30, 1992
*************************************************************************/
void rmatrixtrsafesolve(/* Real    */ ae_matrix* a,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     double* s,
     ae_bool isupper,
     ae_bool istrans,
     ae_bool isunit,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_bool normin;
    ae_vector cnorm;
    ae_matrix a1;
    ae_vector x1;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    *s = 0;
    ae_vector_init(&cnorm, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&a1, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&x1, 0, DT_REAL, _state, ae_true);

    
    /*
     * From 0-based to 1-based
     */
    normin = ae_false;
    ae_matrix_set_length(&a1, n+1, n+1, _state);
    ae_vector_set_length(&x1, n+1, _state);
    for(i=1; i<=n; i++)
    {
        ae_v_move(&a1.ptr.pp_double[i][1], 1, &a->ptr.pp_double[i-1][0], 1, ae_v_len(1,n));
    }
    ae_v_move(&x1.ptr.p_double[1], 1, &x->ptr.p_double[0], 1, ae_v_len(1,n));
    
    /*
     * Solve 1-based
     */
    safesolvetriangular(&a1, n, &x1, s, isupper, istrans, isunit, normin, &cnorm, _state);
    
    /*
     * From 1-based to 0-based
     */
    ae_v_move(&x->ptr.p_double[0], 1, &x1.ptr.p_double[1], 1, ae_v_len(0,n-1));
    ae_frame_leave(_state);
}


/*************************************************************************
Obsolete 1-based subroutine.
See RMatrixTRSafeSolve for 0-based replacement.
*************************************************************************/
void safesolvetriangular(/* Real    */ ae_matrix* a,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     double* s,
     ae_bool isupper,
     ae_bool istrans,
     ae_bool isunit,
     ae_bool normin,
     /* Real    */ ae_vector* cnorm,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t imax;
    ae_int_t j;
    ae_int_t jfirst;
    ae_int_t jinc;
    ae_int_t jlast;
    ae_int_t jm1;
    ae_int_t jp1;
    ae_int_t ip1;
    ae_int_t im1;
    ae_int_t k;
    ae_int_t flg;
    double v;
    double vd;
    double bignum;
    double grow;
    double rec;
    double smlnum;
    double sumj;
    double tjj;
    double tjjs;
    double tmax;
    double tscal;
    double uscal;
    double xbnd;
    double xj;
    double xmax;
    ae_bool notran;
    ae_bool upper;
    ae_bool nounit;

    *s = 0;

    upper = isupper;
    notran = !istrans;
    nounit = !isunit;
    
    /*
     * these initializers are not really necessary,
     * but without them compiler complains about uninitialized locals
     */
    tjjs = 0;
    
    /*
     * Quick return if possible
     */
    if( n==0 )
    {
        return;
    }
    
    /*
     * Determine machine dependent parameters to control overflow.
     */
    smlnum = ae_minrealnumber/(ae_machineepsilon*2);
    bignum = 1/smlnum;
    *s = 1;
    if( !normin )
    {
        ae_vector_set_length(cnorm, n+1, _state);
        
        /*
         * Compute the 1-norm of each column, not including the diagonal.
         */
        if( upper )
        {
            
            /*
             * A is upper triangular.
             */
            for(j=1; j<=n; j++)
            {
                v = 0;
                for(k=1; k<=j-1; k++)
                {
                    v = v+ae_fabs(a->ptr.pp_double[k][j], _state);
                }
                cnorm->ptr.p_double[j] = v;
            }
        }
        else
        {
            
            /*
             * A is lower triangular.
             */
            for(j=1; j<=n-1; j++)
            {
                v = 0;
                for(k=j+1; k<=n; k++)
                {
                    v = v+ae_fabs(a->ptr.pp_double[k][j], _state);
                }
                cnorm->ptr.p_double[j] = v;
            }
            cnorm->ptr.p_double[n] = 0;
        }
    }
    
    /*
     * Scale the column norms by TSCAL if the maximum element in CNORM is
     * greater than BIGNUM.
     */
    imax = 1;
    for(k=2; k<=n; k++)
    {
        if( ae_fp_greater(cnorm->ptr.p_double[k],cnorm->ptr.p_double[imax]) )
        {
            imax = k;
        }
    }
    tmax = cnorm->ptr.p_double[imax];
    if( ae_fp_less_eq(tmax,bignum) )
    {
        tscal = 1;
    }
    else
    {
        tscal = 1/(smlnum*tmax);
        ae_v_muld(&cnorm->ptr.p_double[1], 1, ae_v_len(1,n), tscal);
    }
    
    /*
     * Compute a bound on the computed solution vector to see if the
     * Level 2 BLAS routine DTRSV can be used.
     */
    j = 1;
    for(k=2; k<=n; k++)
    {
        if( ae_fp_greater(ae_fabs(x->ptr.p_double[k], _state),ae_fabs(x->ptr.p_double[j], _state)) )
        {
            j = k;
        }
    }
    xmax = ae_fabs(x->ptr.p_double[j], _state);
    xbnd = xmax;
    if( notran )
    {
        
        /*
         * Compute the growth in A * x = b.
         */
        if( upper )
        {
            jfirst = n;
            jlast = 1;
            jinc = -1;
        }
        else
        {
            jfirst = 1;
            jlast = n;
            jinc = 1;
        }
        if( ae_fp_neq(tscal,1) )
        {
            grow = 0;
        }
        else
        {
            if( nounit )
            {
                
                /*
                 * A is non-unit triangular.
                 *
                 * Compute GROW = 1/G(j) and XBND = 1/M(j).
                 * Initially, G(0) = max{x(i), i=1,...,n}.
                 */
                grow = 1/ae_maxreal(xbnd, smlnum, _state);
                xbnd = grow;
                j = jfirst;
                while((jinc>0&&j<=jlast)||(jinc<0&&j>=jlast))
                {
                    
                    /*
                     * Exit the loop if the growth factor is too small.
                     */
                    if( ae_fp_less_eq(grow,smlnum) )
                    {
                        break;
                    }
                    
                    /*
                     * M(j) = G(j-1) / abs(A(j,j))
                     */
                    tjj = ae_fabs(a->ptr.pp_double[j][j], _state);
                    xbnd = ae_minreal(xbnd, ae_minreal(1, tjj, _state)*grow, _state);
                    if( ae_fp_greater_eq(tjj+cnorm->ptr.p_double[j],smlnum) )
                    {
                        
                        /*
                         * G(j) = G(j-1)*( 1 + CNORM(j) / abs(A(j,j)) )
                         */
                        grow = grow*(tjj/(tjj+cnorm->ptr.p_double[j]));
                    }
                    else
                    {
                        
                        /*
                         * G(j) could overflow, set GROW to 0.
                         */
                        grow = 0;
                    }
                    if( j==jlast )
                    {
                        grow = xbnd;
                    }
                    j = j+jinc;
                }
            }
            else
            {
                
                /*
                 * A is unit triangular.
                 *
                 * Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
                 */
                grow = ae_minreal(1, 1/ae_maxreal(xbnd, smlnum, _state), _state);
                j = jfirst;
                while((jinc>0&&j<=jlast)||(jinc<0&&j>=jlast))
                {
                    
                    /*
                     * Exit the loop if the growth factor is too small.
                     */
                    if( ae_fp_less_eq(grow,smlnum) )
                    {
                        break;
                    }
                    
                    /*
                     * G(j) = G(j-1)*( 1 + CNORM(j) )
                     */
                    grow = grow*(1/(1+cnorm->ptr.p_double[j]));
                    j = j+jinc;
                }
            }
        }
    }
    else
    {
        
        /*
         * Compute the growth in A' * x = b.
         */
        if( upper )
        {
            jfirst = 1;
            jlast = n;
            jinc = 1;
        }
        else
        {
            jfirst = n;
            jlast = 1;
            jinc = -1;
        }
        if( ae_fp_neq(tscal,1) )
        {
            grow = 0;
        }
        else
        {
            if( nounit )
            {
                
                /*
                 * A is non-unit triangular.
                 *
                 * Compute GROW = 1/G(j) and XBND = 1/M(j).
                 * Initially, M(0) = max{x(i), i=1,...,n}.
                 */
                grow = 1/ae_maxreal(xbnd, smlnum, _state);
                xbnd = grow;
                j = jfirst;
                while((jinc>0&&j<=jlast)||(jinc<0&&j>=jlast))
                {
                    
                    /*
                     * Exit the loop if the growth factor is too small.
                     */
                    if( ae_fp_less_eq(grow,smlnum) )
                    {
                        break;
                    }
                    
                    /*
                     * G(j) = max( G(j-1), M(j-1)*( 1 + CNORM(j) ) )
                     */
                    xj = 1+cnorm->ptr.p_double[j];
                    grow = ae_minreal(grow, xbnd/xj, _state);
                    
                    /*
                     * M(j) = M(j-1)*( 1 + CNORM(j) ) / abs(A(j,j))
                     */
                    tjj = ae_fabs(a->ptr.pp_double[j][j], _state);
                    if( ae_fp_greater(xj,tjj) )
                    {
                        xbnd = xbnd*(tjj/xj);
                    }
                    if( j==jlast )
                    {
                        grow = ae_minreal(grow, xbnd, _state);
                    }
                    j = j+jinc;
                }
            }
            else
            {
                
                /*
                 * A is unit triangular.
                 *
                 * Compute GROW = 1/G(j), where G(0) = max{x(i), i=1,...,n}.
                 */
                grow = ae_minreal(1, 1/ae_maxreal(xbnd, smlnum, _state), _state);
                j = jfirst;
                while((jinc>0&&j<=jlast)||(jinc<0&&j>=jlast))
                {
                    
                    /*
                     * Exit the loop if the growth factor is too small.
                     */
                    if( ae_fp_less_eq(grow,smlnum) )
                    {
                        break;
                    }
                    
                    /*
                     * G(j) = ( 1 + CNORM(j) )*G(j-1)
                     */
                    xj = 1+cnorm->ptr.p_double[j];
                    grow = grow/xj;
                    j = j+jinc;
                }
            }
        }
    }
    if( ae_fp_greater(grow*tscal,smlnum) )
    {
        
        /*
         * Use the Level 2 BLAS solve if the reciprocal of the bound on
         * elements of X is not too small.
         */
        if( (upper&&notran)||(!upper&&!notran) )
        {
            if( nounit )
            {
                vd = a->ptr.pp_double[n][n];
            }
            else
            {
                vd = 1;
            }
            x->ptr.p_double[n] = x->ptr.p_double[n]/vd;
            for(i=n-1; i>=1; i--)
            {
                ip1 = i+1;
                if( upper )
                {
                    v = ae_v_dotproduct(&a->ptr.pp_double[i][ip1], 1, &x->ptr.p_double[ip1], 1, ae_v_len(ip1,n));
                }
                else
                {
                    v = ae_v_dotproduct(&a->ptr.pp_double[ip1][i], a->stride, &x->ptr.p_double[ip1], 1, ae_v_len(ip1,n));
                }
                if( nounit )
                {
                    vd = a->ptr.pp_double[i][i];
                }
                else
                {
                    vd = 1;
                }
                x->ptr.p_double[i] = (x->ptr.p_double[i]-v)/vd;
            }
        }
        else
        {
            if( nounit )
            {
                vd = a->ptr.pp_double[1][1];
            }
            else
            {
                vd = 1;
            }
            x->ptr.p_double[1] = x->ptr.p_double[1]/vd;
            for(i=2; i<=n; i++)
            {
                im1 = i-1;
                if( upper )
                {
                    v = ae_v_dotproduct(&a->ptr.pp_double[1][i], a->stride, &x->ptr.p_double[1], 1, ae_v_len(1,im1));
                }
                else
                {
                    v = ae_v_dotproduct(&a->ptr.pp_double[i][1], 1, &x->ptr.p_double[1], 1, ae_v_len(1,im1));
                }
                if( nounit )
                {
                    vd = a->ptr.pp_double[i][i];
                }
                else
                {
                    vd = 1;
                }
                x->ptr.p_double[i] = (x->ptr.p_double[i]-v)/vd;
            }
        }
    }
    else
    {
        
        /*
         * Use a Level 1 BLAS solve, scaling intermediate results.
         */
        if( ae_fp_greater(xmax,bignum) )
        {
            
            /*
             * Scale X so that its components are less than or equal to
             * BIGNUM in absolute value.
             */
            *s = bignum/xmax;
            ae_v_muld(&x->ptr.p_double[1], 1, ae_v_len(1,n), *s);
            xmax = bignum;
        }
        if( notran )
        {
            
            /*
             * Solve A * x = b
             */
            j = jfirst;
            while((jinc>0&&j<=jlast)||(jinc<0&&j>=jlast))
            {
                
                /*
                 * Compute x(j) = b(j) / A(j,j), scaling x if necessary.
                 */
                xj = ae_fabs(x->ptr.p_double[j], _state);
                flg = 0;
                if( nounit )
                {
                    tjjs = a->ptr.pp_double[j][j]*tscal;
                }
                else
                {
                    tjjs = tscal;
                    if( ae_fp_eq(tscal,1) )
                    {
                        flg = 100;
                    }
                }
                if( flg!=100 )
                {
                    tjj = ae_fabs(tjjs, _state);
                    if( ae_fp_greater(tjj,smlnum) )
                    {
                        
                        /*
                         * abs(A(j,j)) > SMLNUM:
                         */
                        if( ae_fp_less(tjj,1) )
                        {
                            if( ae_fp_greater(xj,tjj*bignum) )
                            {
                                
                                /*
                                 * Scale x by 1/b(j).
                                 */
                                rec = 1/xj;
                                ae_v_muld(&x->ptr.p_double[1], 1, ae_v_len(1,n), rec);
                                *s = *s*rec;
                                xmax = xmax*rec;
                            }
                        }
                        x->ptr.p_double[j] = x->ptr.p_double[j]/tjjs;
                        xj = ae_fabs(x->ptr.p_double[j], _state);
                    }
                    else
                    {
                        if( ae_fp_greater(tjj,0) )
                        {
                            
                            /*
                             * 0 < abs(A(j,j)) <= SMLNUM:
                             */
                            if( ae_fp_greater(xj,tjj*bignum) )
                            {
                                
                                /*
                                 * Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM
                                 * to avoid overflow when dividing by A(j,j).
                                 */
                                rec = tjj*bignum/xj;
                                if( ae_fp_greater(cnorm->ptr.p_double[j],1) )
                                {
                                    
                                    /*
                                     * Scale by 1/CNORM(j) to avoid overflow when
                                     * multiplying x(j) times column j.
                                     */
                                    rec = rec/cnorm->ptr.p_double[j];
                                }
                                ae_v_muld(&x->ptr.p_double[1], 1, ae_v_len(1,n), rec);
                                *s = *s*rec;
                                xmax = xmax*rec;
                            }
                            x->ptr.p_double[j] = x->ptr.p_double[j]/tjjs;
                            xj = ae_fabs(x->ptr.p_double[j], _state);
                        }
                        else
                        {
                            
                            /*
                             * A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                             * scale = 0, and compute a solution to A*x = 0.
                             */
                            for(i=1; i<=n; i++)
                            {
                                x->ptr.p_double[i] = 0;
                            }
                            x->ptr.p_double[j] = 1;
                            xj = 1;
                            *s = 0;
                            xmax = 0;
                        }
                    }
                }
                
                /*
                 * Scale x if necessary to avoid overflow when adding a
                 * multiple of column j of A.
                 */
                if( ae_fp_greater(xj,1) )
                {
                    rec = 1/xj;
                    if( ae_fp_greater(cnorm->ptr.p_double[j],(bignum-xmax)*rec) )
                    {
                        
                        /*
                         * Scale x by 1/(2*abs(x(j))).
                         */
                        rec = rec*0.5;
                        ae_v_muld(&x->ptr.p_double[1], 1, ae_v_len(1,n), rec);
                        *s = *s*rec;
                    }
                }
                else
                {
                    if( ae_fp_greater(xj*cnorm->ptr.p_double[j],bignum-xmax) )
                    {
                        
                        /*
                         * Scale x by 1/2.
                         */
                        ae_v_muld(&x->ptr.p_double[1], 1, ae_v_len(1,n), 0.5);
                        *s = *s*0.5;
                    }
                }
                if( upper )
                {
                    if( j>1 )
                    {
                        
                        /*
                         * Compute the update
                         * x(1:j-1) := x(1:j-1) - x(j) * A(1:j-1,j)
                         */
                        v = x->ptr.p_double[j]*tscal;
                        jm1 = j-1;
                        ae_v_subd(&x->ptr.p_double[1], 1, &a->ptr.pp_double[1][j], a->stride, ae_v_len(1,jm1), v);
                        i = 1;
                        for(k=2; k<=j-1; k++)
                        {
                            if( ae_fp_greater(ae_fabs(x->ptr.p_double[k], _state),ae_fabs(x->ptr.p_double[i], _state)) )
                            {
                                i = k;
                            }
                        }
                        xmax = ae_fabs(x->ptr.p_double[i], _state);
                    }
                }
                else
                {
                    if( j<n )
                    {
                        
                        /*
                         * Compute the update
                         * x(j+1:n) := x(j+1:n) - x(j) * A(j+1:n,j)
                         */
                        jp1 = j+1;
                        v = x->ptr.p_double[j]*tscal;
                        ae_v_subd(&x->ptr.p_double[jp1], 1, &a->ptr.pp_double[jp1][j], a->stride, ae_v_len(jp1,n), v);
                        i = j+1;
                        for(k=j+2; k<=n; k++)
                        {
                            if( ae_fp_greater(ae_fabs(x->ptr.p_double[k], _state),ae_fabs(x->ptr.p_double[i], _state)) )
                            {
                                i = k;
                            }
                        }
                        xmax = ae_fabs(x->ptr.p_double[i], _state);
                    }
                }
                j = j+jinc;
            }
        }
        else
        {
            
            /*
             * Solve A' * x = b
             */
            j = jfirst;
            while((jinc>0&&j<=jlast)||(jinc<0&&j>=jlast))
            {
                
                /*
                 * Compute x(j) = b(j) - sum A(k,j)*x(k).
                 *   k<>j
                 */
                xj = ae_fabs(x->ptr.p_double[j], _state);
                uscal = tscal;
                rec = 1/ae_maxreal(xmax, 1, _state);
                if( ae_fp_greater(cnorm->ptr.p_double[j],(bignum-xj)*rec) )
                {
                    
                    /*
                     * If x(j) could overflow, scale x by 1/(2*XMAX).
                     */
                    rec = rec*0.5;
                    if( nounit )
                    {
                        tjjs = a->ptr.pp_double[j][j]*tscal;
                    }
                    else
                    {
                        tjjs = tscal;
                    }
                    tjj = ae_fabs(tjjs, _state);
                    if( ae_fp_greater(tjj,1) )
                    {
                        
                        /*
                         * Divide by A(j,j) when scaling x if A(j,j) > 1.
                         */
                        rec = ae_minreal(1, rec*tjj, _state);
                        uscal = uscal/tjjs;
                    }
                    if( ae_fp_less(rec,1) )
                    {
                        ae_v_muld(&x->ptr.p_double[1], 1, ae_v_len(1,n), rec);
                        *s = *s*rec;
                        xmax = xmax*rec;
                    }
                }
                sumj = 0;
                if( ae_fp_eq(uscal,1) )
                {
                    
                    /*
                     * If the scaling needed for A in the dot product is 1,
                     * call DDOT to perform the dot product.
                     */
                    if( upper )
                    {
                        if( j>1 )
                        {
                            jm1 = j-1;
                            sumj = ae_v_dotproduct(&a->ptr.pp_double[1][j], a->stride, &x->ptr.p_double[1], 1, ae_v_len(1,jm1));
                        }
                        else
                        {
                            sumj = 0;
                        }
                    }
                    else
                    {
                        if( j<n )
                        {
                            jp1 = j+1;
                            sumj = ae_v_dotproduct(&a->ptr.pp_double[jp1][j], a->stride, &x->ptr.p_double[jp1], 1, ae_v_len(jp1,n));
                        }
                    }
                }
                else
                {
                    
                    /*
                     * Otherwise, use in-line code for the dot product.
                     */
                    if( upper )
                    {
                        for(i=1; i<=j-1; i++)
                        {
                            v = a->ptr.pp_double[i][j]*uscal;
                            sumj = sumj+v*x->ptr.p_double[i];
                        }
                    }
                    else
                    {
                        if( j<n )
                        {
                            for(i=j+1; i<=n; i++)
                            {
                                v = a->ptr.pp_double[i][j]*uscal;
                                sumj = sumj+v*x->ptr.p_double[i];
                            }
                        }
                    }
                }
                if( ae_fp_eq(uscal,tscal) )
                {
                    
                    /*
                     * Compute x(j) := ( x(j) - sumj ) / A(j,j) if 1/A(j,j)
                     * was not used to scale the dotproduct.
                     */
                    x->ptr.p_double[j] = x->ptr.p_double[j]-sumj;
                    xj = ae_fabs(x->ptr.p_double[j], _state);
                    flg = 0;
                    if( nounit )
                    {
                        tjjs = a->ptr.pp_double[j][j]*tscal;
                    }
                    else
                    {
                        tjjs = tscal;
                        if( ae_fp_eq(tscal,1) )
                        {
                            flg = 150;
                        }
                    }
                    
                    /*
                     * Compute x(j) = x(j) / A(j,j), scaling if necessary.
                     */
                    if( flg!=150 )
                    {
                        tjj = ae_fabs(tjjs, _state);
                        if( ae_fp_greater(tjj,smlnum) )
                        {
                            
                            /*
                             * abs(A(j,j)) > SMLNUM:
                             */
                            if( ae_fp_less(tjj,1) )
                            {
                                if( ae_fp_greater(xj,tjj*bignum) )
                                {
                                    
                                    /*
                                     * Scale X by 1/abs(x(j)).
                                     */
                                    rec = 1/xj;
                                    ae_v_muld(&x->ptr.p_double[1], 1, ae_v_len(1,n), rec);
                                    *s = *s*rec;
                                    xmax = xmax*rec;
                                }
                            }
                            x->ptr.p_double[j] = x->ptr.p_double[j]/tjjs;
                        }
                        else
                        {
                            if( ae_fp_greater(tjj,0) )
                            {
                                
                                /*
                                 * 0 < abs(A(j,j)) <= SMLNUM:
                                 */
                                if( ae_fp_greater(xj,tjj*bignum) )
                                {
                                    
                                    /*
                                     * Scale x by (1/abs(x(j)))*abs(A(j,j))*BIGNUM.
                                     */
                                    rec = tjj*bignum/xj;
                                    ae_v_muld(&x->ptr.p_double[1], 1, ae_v_len(1,n), rec);
                                    *s = *s*rec;
                                    xmax = xmax*rec;
                                }
                                x->ptr.p_double[j] = x->ptr.p_double[j]/tjjs;
                            }
                            else
                            {
                                
                                /*
                                 * A(j,j) = 0:  Set x(1:n) = 0, x(j) = 1, and
                                 * scale = 0, and compute a solution to A'*x = 0.
                                 */
                                for(i=1; i<=n; i++)
                                {
                                    x->ptr.p_double[i] = 0;
                                }
                                x->ptr.p_double[j] = 1;
                                *s = 0;
                                xmax = 0;
                            }
                        }
                    }
                }
                else
                {
                    
                    /*
                     * Compute x(j) := x(j) / A(j,j)  - sumj if the dot
                     * product has already been divided by 1/A(j,j).
                     */
                    x->ptr.p_double[j] = x->ptr.p_double[j]/tjjs-sumj;
                }
                xmax = ae_maxreal(xmax, ae_fabs(x->ptr.p_double[j], _state), _state);
                j = j+jinc;
            }
        }
        *s = *s/tscal;
    }
    
    /*
     * Scale the column norms by 1/TSCAL for return.
     */
    if( ae_fp_neq(tscal,1) )
    {
        v = 1/tscal;
        ae_v_muld(&cnorm->ptr.p_double[1], 1, ae_v_len(1,n), v);
    }
}




/*************************************************************************
Real implementation of CMatrixScaledTRSafeSolve

  -- ALGLIB routine --
     21.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool rmatrixscaledtrsafesolve(/* Real    */ ae_matrix* a,
     double sa,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     ae_bool isupper,
     ae_int_t trans,
     ae_bool isunit,
     double maxgrowth,
     ae_state *_state)
{
    ae_frame _frame_block;
    double lnmax;
    double nrmb;
    double nrmx;
    ae_int_t i;
    ae_complex alpha;
    ae_complex beta;
    double vr;
    ae_complex cx;
    ae_vector tmp;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);

    ae_assert(n>0, "RMatrixTRSafeSolve: incorrect N!", _state);
    ae_assert(trans==0||trans==1, "RMatrixTRSafeSolve: incorrect Trans!", _state);
    result = ae_true;
    lnmax = ae_log(ae_maxrealnumber, _state);
    
    /*
     * Quick return if possible
     */
    if( n<=0 )
    {
        ae_frame_leave(_state);
        return result;
    }
    
    /*
     * Load norms: right part and X
     */
    nrmb = 0;
    for(i=0; i<=n-1; i++)
    {
        nrmb = ae_maxreal(nrmb, ae_fabs(x->ptr.p_double[i], _state), _state);
    }
    nrmx = 0;
    
    /*
     * Solve
     */
    ae_vector_set_length(&tmp, n, _state);
    result = ae_true;
    if( isupper&&trans==0 )
    {
        
        /*
         * U*x = b
         */
        for(i=n-1; i>=0; i--)
        {
            
            /*
             * Task is reduced to alpha*x[i] = beta
             */
            if( isunit )
            {
                alpha = ae_complex_from_d(sa);
            }
            else
            {
                alpha = ae_complex_from_d(a->ptr.pp_double[i][i]*sa);
            }
            if( i<n-1 )
            {
                ae_v_moved(&tmp.ptr.p_double[i+1], 1, &a->ptr.pp_double[i][i+1], 1, ae_v_len(i+1,n-1), sa);
                vr = ae_v_dotproduct(&tmp.ptr.p_double[i+1], 1, &x->ptr.p_double[i+1], 1, ae_v_len(i+1,n-1));
                beta = ae_complex_from_d(x->ptr.p_double[i]-vr);
            }
            else
            {
                beta = ae_complex_from_d(x->ptr.p_double[i]);
            }
            
            /*
             * solve alpha*x[i] = beta
             */
            result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &cx, _state);
            if( !result )
            {
                ae_frame_leave(_state);
                return result;
            }
            x->ptr.p_double[i] = cx.x;
        }
        ae_frame_leave(_state);
        return result;
    }
    if( !isupper&&trans==0 )
    {
        
        /*
         * L*x = b
         */
        for(i=0; i<=n-1; i++)
        {
            
            /*
             * Task is reduced to alpha*x[i] = beta
             */
            if( isunit )
            {
                alpha = ae_complex_from_d(sa);
            }
            else
            {
                alpha = ae_complex_from_d(a->ptr.pp_double[i][i]*sa);
            }
            if( i>0 )
            {
                ae_v_moved(&tmp.ptr.p_double[0], 1, &a->ptr.pp_double[i][0], 1, ae_v_len(0,i-1), sa);
                vr = ae_v_dotproduct(&tmp.ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,i-1));
                beta = ae_complex_from_d(x->ptr.p_double[i]-vr);
            }
            else
            {
                beta = ae_complex_from_d(x->ptr.p_double[i]);
            }
            
            /*
             * solve alpha*x[i] = beta
             */
            result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &cx, _state);
            if( !result )
            {
                ae_frame_leave(_state);
                return result;
            }
            x->ptr.p_double[i] = cx.x;
        }
        ae_frame_leave(_state);
        return result;
    }
    if( isupper&&trans==1 )
    {
        
        /*
         * U^T*x = b
         */
        for(i=0; i<=n-1; i++)
        {
            
            /*
             * Task is reduced to alpha*x[i] = beta
             */
            if( isunit )
            {
                alpha = ae_complex_from_d(sa);
            }
            else
            {
                alpha = ae_complex_from_d(a->ptr.pp_double[i][i]*sa);
            }
            beta = ae_complex_from_d(x->ptr.p_double[i]);
            
            /*
             * solve alpha*x[i] = beta
             */
            result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &cx, _state);
            if( !result )
            {
                ae_frame_leave(_state);
                return result;
            }
            x->ptr.p_double[i] = cx.x;
            
            /*
             * update the rest of right part
             */
            if( i<n-1 )
            {
                vr = cx.x;
                ae_v_moved(&tmp.ptr.p_double[i+1], 1, &a->ptr.pp_double[i][i+1], 1, ae_v_len(i+1,n-1), sa);
                ae_v_subd(&x->ptr.p_double[i+1], 1, &tmp.ptr.p_double[i+1], 1, ae_v_len(i+1,n-1), vr);
            }
        }
        ae_frame_leave(_state);
        return result;
    }
    if( !isupper&&trans==1 )
    {
        
        /*
         * L^T*x = b
         */
        for(i=n-1; i>=0; i--)
        {
            
            /*
             * Task is reduced to alpha*x[i] = beta
             */
            if( isunit )
            {
                alpha = ae_complex_from_d(sa);
            }
            else
            {
                alpha = ae_complex_from_d(a->ptr.pp_double[i][i]*sa);
            }
            beta = ae_complex_from_d(x->ptr.p_double[i]);
            
            /*
             * solve alpha*x[i] = beta
             */
            result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &cx, _state);
            if( !result )
            {
                ae_frame_leave(_state);
                return result;
            }
            x->ptr.p_double[i] = cx.x;
            
            /*
             * update the rest of right part
             */
            if( i>0 )
            {
                vr = cx.x;
                ae_v_moved(&tmp.ptr.p_double[0], 1, &a->ptr.pp_double[i][0], 1, ae_v_len(0,i-1), sa);
                ae_v_subd(&x->ptr.p_double[0], 1, &tmp.ptr.p_double[0], 1, ae_v_len(0,i-1), vr);
            }
        }
        ae_frame_leave(_state);
        return result;
    }
    result = ae_false;
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
Internal subroutine for safe solution of

    SA*op(A)=b
    
where  A  is  NxN  upper/lower  triangular/unitriangular  matrix, op(A) is
either identity transform, transposition or Hermitian transposition, SA is
a scaling factor such that max(|SA*A[i,j]|) is close to 1.0 in magnutude.

This subroutine  limits  relative  growth  of  solution  (in inf-norm)  by
MaxGrowth,  returning  False  if  growth  exceeds MaxGrowth. Degenerate or
near-degenerate matrices are handled correctly (False is returned) as long
as MaxGrowth is significantly less than MaxRealNumber/norm(b).

  -- ALGLIB routine --
     21.01.2010
     Bochkanov Sergey
*************************************************************************/
ae_bool cmatrixscaledtrsafesolve(/* Complex */ ae_matrix* a,
     double sa,
     ae_int_t n,
     /* Complex */ ae_vector* x,
     ae_bool isupper,
     ae_int_t trans,
     ae_bool isunit,
     double maxgrowth,
     ae_state *_state)
{
    ae_frame _frame_block;
    double lnmax;
    double nrmb;
    double nrmx;
    ae_int_t i;
    ae_complex alpha;
    ae_complex beta;
    ae_complex vc;
    ae_vector tmp;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&tmp, 0, DT_COMPLEX, _state, ae_true);

    ae_assert(n>0, "CMatrixTRSafeSolve: incorrect N!", _state);
    ae_assert((trans==0||trans==1)||trans==2, "CMatrixTRSafeSolve: incorrect Trans!", _state);
    result = ae_true;
    lnmax = ae_log(ae_maxrealnumber, _state);
    
    /*
     * Quick return if possible
     */
    if( n<=0 )
    {
        ae_frame_leave(_state);
        return result;
    }
    
    /*
     * Load norms: right part and X
     */
    nrmb = 0;
    for(i=0; i<=n-1; i++)
    {
        nrmb = ae_maxreal(nrmb, ae_c_abs(x->ptr.p_complex[i], _state), _state);
    }
    nrmx = 0;
    
    /*
     * Solve
     */
    ae_vector_set_length(&tmp, n, _state);
    result = ae_true;
    if( isupper&&trans==0 )
    {
        
        /*
         * U*x = b
         */
        for(i=n-1; i>=0; i--)
        {
            
            /*
             * Task is reduced to alpha*x[i] = beta
             */
            if( isunit )
            {
                alpha = ae_complex_from_d(sa);
            }
            else
            {
                alpha = ae_c_mul_d(a->ptr.pp_complex[i][i],sa);
            }
            if( i<n-1 )
            {
                ae_v_cmoved(&tmp.ptr.p_complex[i+1], 1, &a->ptr.pp_complex[i][i+1], 1, "N", ae_v_len(i+1,n-1), sa);
                vc = ae_v_cdotproduct(&tmp.ptr.p_complex[i+1], 1, "N", &x->ptr.p_complex[i+1], 1, "N", ae_v_len(i+1,n-1));
                beta = ae_c_sub(x->ptr.p_complex[i],vc);
            }
            else
            {
                beta = x->ptr.p_complex[i];
            }
            
            /*
             * solve alpha*x[i] = beta
             */
            result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc, _state);
            if( !result )
            {
                ae_frame_leave(_state);
                return result;
            }
            x->ptr.p_complex[i] = vc;
        }
        ae_frame_leave(_state);
        return result;
    }
    if( !isupper&&trans==0 )
    {
        
        /*
         * L*x = b
         */
        for(i=0; i<=n-1; i++)
        {
            
            /*
             * Task is reduced to alpha*x[i] = beta
             */
            if( isunit )
            {
                alpha = ae_complex_from_d(sa);
            }
            else
            {
                alpha = ae_c_mul_d(a->ptr.pp_complex[i][i],sa);
            }
            if( i>0 )
            {
                ae_v_cmoved(&tmp.ptr.p_complex[0], 1, &a->ptr.pp_complex[i][0], 1, "N", ae_v_len(0,i-1), sa);
                vc = ae_v_cdotproduct(&tmp.ptr.p_complex[0], 1, "N", &x->ptr.p_complex[0], 1, "N", ae_v_len(0,i-1));
                beta = ae_c_sub(x->ptr.p_complex[i],vc);
            }
            else
            {
                beta = x->ptr.p_complex[i];
            }
            
            /*
             * solve alpha*x[i] = beta
             */
            result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc, _state);
            if( !result )
            {
                ae_frame_leave(_state);
                return result;
            }
            x->ptr.p_complex[i] = vc;
        }
        ae_frame_leave(_state);
        return result;
    }
    if( isupper&&trans==1 )
    {
        
        /*
         * U^T*x = b
         */
        for(i=0; i<=n-1; i++)
        {
            
            /*
             * Task is reduced to alpha*x[i] = beta
             */
            if( isunit )
            {
                alpha = ae_complex_from_d(sa);
            }
            else
            {
                alpha = ae_c_mul_d(a->ptr.pp_complex[i][i],sa);
            }
            beta = x->ptr.p_complex[i];
            
            /*
             * solve alpha*x[i] = beta
             */
            result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc, _state);
            if( !result )
            {
                ae_frame_leave(_state);
                return result;
            }
            x->ptr.p_complex[i] = vc;
            
            /*
             * update the rest of right part
             */
            if( i<n-1 )
            {
                ae_v_cmoved(&tmp.ptr.p_complex[i+1], 1, &a->ptr.pp_complex[i][i+1], 1, "N", ae_v_len(i+1,n-1), sa);
                ae_v_csubc(&x->ptr.p_complex[i+1], 1, &tmp.ptr.p_complex[i+1], 1, "N", ae_v_len(i+1,n-1), vc);
            }
        }
        ae_frame_leave(_state);
        return result;
    }
    if( !isupper&&trans==1 )
    {
        
        /*
         * L^T*x = b
         */
        for(i=n-1; i>=0; i--)
        {
            
            /*
             * Task is reduced to alpha*x[i] = beta
             */
            if( isunit )
            {
                alpha = ae_complex_from_d(sa);
            }
            else
            {
                alpha = ae_c_mul_d(a->ptr.pp_complex[i][i],sa);
            }
            beta = x->ptr.p_complex[i];
            
            /*
             * solve alpha*x[i] = beta
             */
            result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc, _state);
            if( !result )
            {
                ae_frame_leave(_state);
                return result;
            }
            x->ptr.p_complex[i] = vc;
            
            /*
             * update the rest of right part
             */
            if( i>0 )
            {
                ae_v_cmoved(&tmp.ptr.p_complex[0], 1, &a->ptr.pp_complex[i][0], 1, "N", ae_v_len(0,i-1), sa);
                ae_v_csubc(&x->ptr.p_complex[0], 1, &tmp.ptr.p_complex[0], 1, "N", ae_v_len(0,i-1), vc);
            }
        }
        ae_frame_leave(_state);
        return result;
    }
    if( isupper&&trans==2 )
    {
        
        /*
         * U^H*x = b
         */
        for(i=0; i<=n-1; i++)
        {
            
            /*
             * Task is reduced to alpha*x[i] = beta
             */
            if( isunit )
            {
                alpha = ae_complex_from_d(sa);
            }
            else
            {
                alpha = ae_c_mul_d(ae_c_conj(a->ptr.pp_complex[i][i], _state),sa);
            }
            beta = x->ptr.p_complex[i];
            
            /*
             * solve alpha*x[i] = beta
             */
            result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc, _state);
            if( !result )
            {
                ae_frame_leave(_state);
                return result;
            }
            x->ptr.p_complex[i] = vc;
            
            /*
             * update the rest of right part
             */
            if( i<n-1 )
            {
                ae_v_cmoved(&tmp.ptr.p_complex[i+1], 1, &a->ptr.pp_complex[i][i+1], 1, "Conj", ae_v_len(i+1,n-1), sa);
                ae_v_csubc(&x->ptr.p_complex[i+1], 1, &tmp.ptr.p_complex[i+1], 1, "N", ae_v_len(i+1,n-1), vc);
            }
        }
        ae_frame_leave(_state);
        return result;
    }
    if( !isupper&&trans==2 )
    {
        
        /*
         * L^T*x = b
         */
        for(i=n-1; i>=0; i--)
        {
            
            /*
             * Task is reduced to alpha*x[i] = beta
             */
            if( isunit )
            {
                alpha = ae_complex_from_d(sa);
            }
            else
            {
                alpha = ae_c_mul_d(ae_c_conj(a->ptr.pp_complex[i][i], _state),sa);
            }
            beta = x->ptr.p_complex[i];
            
            /*
             * solve alpha*x[i] = beta
             */
            result = safesolve_cbasicsolveandupdate(alpha, beta, lnmax, nrmb, maxgrowth, &nrmx, &vc, _state);
            if( !result )
            {
                ae_frame_leave(_state);
                return result;
            }
            x->ptr.p_complex[i] = vc;
            
            /*
             * update the rest of right part
             */
            if( i>0 )
            {
                ae_v_cmoved(&tmp.ptr.p_complex[0], 1, &a->ptr.pp_complex[i][0], 1, "Conj", ae_v_len(0,i-1), sa);
                ae_v_csubc(&x->ptr.p_complex[0], 1, &tmp.ptr.p_complex[0], 1, "N", ae_v_len(0,i-1), vc);
            }
        }
        ae_frame_leave(_state);
        return result;
    }
    result = ae_false;
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
complex basic solver-updater for reduced linear system

    alpha*x[i] = beta

solves this equation and updates it in overlfow-safe manner (keeping track
of relative growth of solution).

Parameters:
    Alpha   -   alpha
    Beta    -   beta
    LnMax   -   precomputed Ln(MaxRealNumber)
    BNorm   -   inf-norm of b (right part of original system)
    MaxGrowth-  maximum growth of norm(x) relative to norm(b)
    XNorm   -   inf-norm of other components of X (which are already processed)
                it is updated by CBasicSolveAndUpdate.
    X       -   solution

  -- ALGLIB routine --
     26.01.2009
     Bochkanov Sergey
*************************************************************************/
static ae_bool safesolve_cbasicsolveandupdate(ae_complex alpha,
     ae_complex beta,
     double lnmax,
     double bnorm,
     double maxgrowth,
     double* xnorm,
     ae_complex* x,
     ae_state *_state)
{
    double v;
    ae_bool result;

    x->x = 0;
    x->y = 0;

    result = ae_false;
    if( ae_c_eq_d(alpha,0) )
    {
        return result;
    }
    if( ae_c_neq_d(beta,0) )
    {
        
        /*
         * alpha*x[i]=beta
         */
        v = ae_log(ae_c_abs(beta, _state), _state)-ae_log(ae_c_abs(alpha, _state), _state);
        if( ae_fp_greater(v,lnmax) )
        {
            return result;
        }
        *x = ae_c_div(beta,alpha);
    }
    else
    {
        
        /*
         * alpha*x[i]=0
         */
        *x = ae_complex_from_d(0);
    }
    
    /*
     * update NrmX, test growth limit
     */
    *xnorm = ae_maxreal(*xnorm, ae_c_abs(*x, _state), _state);
    if( ae_fp_greater(*xnorm,maxgrowth*bnorm) )
    {
        return result;
    }
    result = ae_true;
    return result;
}




/*************************************************************************
More precise dot-product. Absolute error of  subroutine  result  is  about
1 ulp of max(MX,V), where:
    MX = max( |a[i]*b[i]| )
    V  = |(a,b)|

INPUT PARAMETERS
    A       -   array[0..N-1], vector 1
    B       -   array[0..N-1], vector 2
    N       -   vectors length, N<2^29.
    Temp    -   array[0..N-1], pre-allocated temporary storage

OUTPUT PARAMETERS
    R       -   (A,B)
    RErr    -   estimate of error. This estimate accounts for both  errors
                during  calculation  of  (A,B)  and  errors  introduced by
                rounding of A and B to fit in double (about 1 ulp).

  -- ALGLIB --
     Copyright 24.08.2009 by Bochkanov Sergey
*************************************************************************/
void xdot(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* b,
     ae_int_t n,
     /* Real    */ ae_vector* temp,
     double* r,
     double* rerr,
     ae_state *_state)
{
    ae_int_t i;
    double mx;
    double v;

    *r = 0;
    *rerr = 0;

    
    /*
     * special cases:
     * * N=0
     */
    if( n==0 )
    {
        *r = 0;
        *rerr = 0;
        return;
    }
    mx = 0;
    for(i=0; i<=n-1; i++)
    {
        v = a->ptr.p_double[i]*b->ptr.p_double[i];
        temp->ptr.p_double[i] = v;
        mx = ae_maxreal(mx, ae_fabs(v, _state), _state);
    }
    if( ae_fp_eq(mx,0) )
    {
        *r = 0;
        *rerr = 0;
        return;
    }
    xblas_xsum(temp, mx, n, r, rerr, _state);
}


/*************************************************************************
More precise complex dot-product. Absolute error of  subroutine  result is
about 1 ulp of max(MX,V), where:
    MX = max( |a[i]*b[i]| )
    V  = |(a,b)|

INPUT PARAMETERS
    A       -   array[0..N-1], vector 1
    B       -   array[0..N-1], vector 2
    N       -   vectors length, N<2^29.
    Temp    -   array[0..2*N-1], pre-allocated temporary storage

OUTPUT PARAMETERS
    R       -   (A,B)
    RErr    -   estimate of error. This estimate accounts for both  errors
                during  calculation  of  (A,B)  and  errors  introduced by
                rounding of A and B to fit in double (about 1 ulp).

  -- ALGLIB --
     Copyright 27.01.2010 by Bochkanov Sergey
*************************************************************************/
void xcdot(/* Complex */ ae_vector* a,
     /* Complex */ ae_vector* b,
     ae_int_t n,
     /* Real    */ ae_vector* temp,
     ae_complex* r,
     double* rerr,
     ae_state *_state)
{
    ae_int_t i;
    double mx;
    double v;
    double rerrx;
    double rerry;

    r->x = 0;
    r->y = 0;
    *rerr = 0;

    
    /*
     * special cases:
     * * N=0
     */
    if( n==0 )
    {
        *r = ae_complex_from_d(0);
        *rerr = 0;
        return;
    }
    
    /*
     * calculate real part
     */
    mx = 0;
    for(i=0; i<=n-1; i++)
    {
        v = a->ptr.p_complex[i].x*b->ptr.p_complex[i].x;
        temp->ptr.p_double[2*i+0] = v;
        mx = ae_maxreal(mx, ae_fabs(v, _state), _state);
        v = -a->ptr.p_complex[i].y*b->ptr.p_complex[i].y;
        temp->ptr.p_double[2*i+1] = v;
        mx = ae_maxreal(mx, ae_fabs(v, _state), _state);
    }
    if( ae_fp_eq(mx,0) )
    {
        r->x = 0;
        rerrx = 0;
    }
    else
    {
        xblas_xsum(temp, mx, 2*n, &r->x, &rerrx, _state);
    }
    
    /*
     * calculate imaginary part
     */
    mx = 0;
    for(i=0; i<=n-1; i++)
    {
        v = a->ptr.p_complex[i].x*b->ptr.p_complex[i].y;
        temp->ptr.p_double[2*i+0] = v;
        mx = ae_maxreal(mx, ae_fabs(v, _state), _state);
        v = a->ptr.p_complex[i].y*b->ptr.p_complex[i].x;
        temp->ptr.p_double[2*i+1] = v;
        mx = ae_maxreal(mx, ae_fabs(v, _state), _state);
    }
    if( ae_fp_eq(mx,0) )
    {
        r->y = 0;
        rerry = 0;
    }
    else
    {
        xblas_xsum(temp, mx, 2*n, &r->y, &rerry, _state);
    }
    
    /*
     * total error
     */
    if( ae_fp_eq(rerrx,0)&&ae_fp_eq(rerry,0) )
    {
        *rerr = 0;
    }
    else
    {
        *rerr = ae_maxreal(rerrx, rerry, _state)*ae_sqrt(1+ae_sqr(ae_minreal(rerrx, rerry, _state)/ae_maxreal(rerrx, rerry, _state), _state), _state);
    }
}


/*************************************************************************
Internal subroutine for extra-precise calculation of SUM(w[i]).

INPUT PARAMETERS:
    W   -   array[0..N-1], values to be added
            W is modified during calculations.
    MX  -   max(W[i])
    N   -   array size
    
OUTPUT PARAMETERS:
    R   -   SUM(w[i])
    RErr-   error estimate for R

  -- ALGLIB --
     Copyright 24.08.2009 by Bochkanov Sergey
*************************************************************************/
static void xblas_xsum(/* Real    */ ae_vector* w,
     double mx,
     ae_int_t n,
     double* r,
     double* rerr,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t k;
    ae_int_t ks;
    double v;
    double s;
    double ln2;
    double chunk;
    double invchunk;
    ae_bool allzeros;

    *r = 0;
    *rerr = 0;

    
    /*
     * special cases:
     * * N=0
     * * N is too large to use integer arithmetics
     */
    if( n==0 )
    {
        *r = 0;
        *rerr = 0;
        return;
    }
    if( ae_fp_eq(mx,0) )
    {
        *r = 0;
        *rerr = 0;
        return;
    }
    ae_assert(n<536870912, "XDot: N is too large!", _state);
    
    /*
     * Prepare
     */
    ln2 = ae_log(2, _state);
    *rerr = mx*ae_machineepsilon;
    
    /*
     * 1. find S such that 0.5<=S*MX<1
     * 2. multiply W by S, so task is normalized in some sense
     * 3. S:=1/S so we can obtain original vector multiplying by S
     */
    k = ae_round(ae_log(mx, _state)/ln2, _state);
    s = xblas_xfastpow(2, -k, _state);
    while(ae_fp_greater_eq(s*mx,1))
    {
        s = 0.5*s;
    }
    while(ae_fp_less(s*mx,0.5))
    {
        s = 2*s;
    }
    ae_v_muld(&w->ptr.p_double[0], 1, ae_v_len(0,n-1), s);
    s = 1/s;
    
    /*
     * find Chunk=2^M such that N*Chunk<2^29
     *
     * we have chosen upper limit (2^29) with enough space left
     * to tolerate possible problems with rounding and N's close
     * to the limit, so we don't want to be very strict here.
     */
    k = ae_trunc(ae_log((double)536870912/(double)n, _state)/ln2, _state);
    chunk = xblas_xfastpow(2, k, _state);
    if( ae_fp_less(chunk,2) )
    {
        chunk = 2;
    }
    invchunk = 1/chunk;
    
    /*
     * calculate result
     */
    *r = 0;
    ae_v_muld(&w->ptr.p_double[0], 1, ae_v_len(0,n-1), chunk);
    for(;;)
    {
        s = s*invchunk;
        allzeros = ae_true;
        ks = 0;
        for(i=0; i<=n-1; i++)
        {
            v = w->ptr.p_double[i];
            k = ae_trunc(v, _state);
            if( ae_fp_neq(v,k) )
            {
                allzeros = ae_false;
            }
            w->ptr.p_double[i] = chunk*(v-k);
            ks = ks+k;
        }
        *r = *r+s*ks;
        v = ae_fabs(*r, _state);
        if( allzeros||ae_fp_eq(s*n+mx,mx) )
        {
            break;
        }
    }
    
    /*
     * correct error
     */
    *rerr = ae_maxreal(*rerr, ae_fabs(*r, _state)*ae_machineepsilon, _state);
}


/*************************************************************************
Fast Pow

  -- ALGLIB --
     Copyright 24.08.2009 by Bochkanov Sergey
*************************************************************************/
static double xblas_xfastpow(double r, ae_int_t n, ae_state *_state)
{
    double result;


    result = 0;
    if( n>0 )
    {
        if( n%2==0 )
        {
            result = ae_sqr(xblas_xfastpow(r, n/2, _state), _state);
        }
        else
        {
            result = r*xblas_xfastpow(r, n-1, _state);
        }
        return result;
    }
    if( n==0 )
    {
        result = 1;
    }
    if( n<0 )
    {
        result = xblas_xfastpow(1/r, -n, _state);
    }
    return result;
}




/*************************************************************************
Normalizes direction/step pair: makes |D|=1, scales Stp.
If |D|=0, it returns, leavind D/Stp unchanged.

  -- ALGLIB --
     Copyright 01.04.2010 by Bochkanov Sergey
*************************************************************************/
void linminnormalized(/* Real    */ ae_vector* d,
     double* stp,
     ae_int_t n,
     ae_state *_state)
{
    double mx;
    double s;
    ae_int_t i;


    
    /*
     * first, scale D to avoid underflow/overflow durng squaring
     */
    mx = 0;
    for(i=0; i<=n-1; i++)
    {
        mx = ae_maxreal(mx, ae_fabs(d->ptr.p_double[i], _state), _state);
    }
    if( ae_fp_eq(mx,0) )
    {
        return;
    }
    s = 1/mx;
    ae_v_muld(&d->ptr.p_double[0], 1, ae_v_len(0,n-1), s);
    *stp = *stp/s;
    
    /*
     * normalize D
     */
    s = ae_v_dotproduct(&d->ptr.p_double[0], 1, &d->ptr.p_double[0], 1, ae_v_len(0,n-1));
    s = 1/ae_sqrt(s, _state);
    ae_v_muld(&d->ptr.p_double[0], 1, ae_v_len(0,n-1), s);
    *stp = *stp/s;
}


/*************************************************************************
THE  PURPOSE  OF  MCSRCH  IS  TO  FIND A STEP WHICH SATISFIES A SUFFICIENT
DECREASE CONDITION AND A CURVATURE CONDITION.

AT EACH STAGE THE SUBROUTINE  UPDATES  AN  INTERVAL  OF  UNCERTAINTY  WITH
ENDPOINTS  STX  AND  STY.  THE INTERVAL OF UNCERTAINTY IS INITIALLY CHOSEN
SO THAT IT CONTAINS A MINIMIZER OF THE MODIFIED FUNCTION

    F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).

IF  A STEP  IS OBTAINED FOR  WHICH THE MODIFIED FUNCTION HAS A NONPOSITIVE
FUNCTION  VALUE  AND  NONNEGATIVE  DERIVATIVE,   THEN   THE   INTERVAL  OF
UNCERTAINTY IS CHOSEN SO THAT IT CONTAINS A MINIMIZER OF F(X+STP*S).

THE  ALGORITHM  IS  DESIGNED TO FIND A STEP WHICH SATISFIES THE SUFFICIENT
DECREASE CONDITION

    F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S),

AND THE CURVATURE CONDITION

    ABS(GRADF(X+STP*S)'S)) .LE. GTOL*ABS(GRADF(X)'S).

IF  FTOL  IS  LESS  THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION IS BOUNDED
BELOW,  THEN  THERE  IS  ALWAYS  A  STEP  WHICH SATISFIES BOTH CONDITIONS.
IF  NO  STEP  CAN BE FOUND  WHICH  SATISFIES  BOTH  CONDITIONS,  THEN  THE
ALGORITHM  USUALLY STOPS  WHEN  ROUNDING ERRORS  PREVENT FURTHER PROGRESS.
IN THIS CASE STP ONLY SATISFIES THE SUFFICIENT DECREASE CONDITION.


:::::::::::::IMPORTANT NOTES:::::::::::::

NOTE 1:

This routine  guarantees that it will stop at the last point where function
value was calculated. It won't make several additional function evaluations
after finding good point. So if you store function evaluations requested by
this routine, you can be sure that last one is the point where we've stopped.

NOTE 2:

when 0<StpMax<StpMin, algorithm will terminate with INFO=5 and Stp=0.0
:::::::::::::::::::::::::::::::::::::::::


PARAMETERS DESCRIPRION

STAGE IS ZERO ON FIRST CALL, ZERO ON FINAL EXIT

N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF VARIABLES.

X IS  AN  ARRAY  OF  LENGTH N. ON INPUT IT MUST CONTAIN THE BASE POINT FOR
THE LINE SEARCH. ON OUTPUT IT CONTAINS X+STP*S.

F IS  A  VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F AT X. ON OUTPUT
IT CONTAINS THE VALUE OF F AT X + STP*S.

G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE GRADIENT OF F AT X.
ON OUTPUT IT CONTAINS THE GRADIENT OF F AT X + STP*S.

S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE SEARCH DIRECTION.

STP  IS  A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN INITIAL ESTIMATE
OF A SATISFACTORY STEP. ON OUTPUT STP CONTAINS THE FINAL ESTIMATE.

FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. TERMINATION OCCURS WHEN THE
SUFFICIENT DECREASE CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE
SATISFIED.

XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS WHEN THE RELATIVE
WIDTH OF THE INTERVAL OF UNCERTAINTY IS AT MOST XTOL.

STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH SPECIFY LOWER  AND
UPPER BOUNDS FOR THE STEP.

MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION OCCURS WHEN THE
NUMBER OF CALLS TO FCN IS AT LEAST MAXFEV BY THE END OF AN ITERATION.

INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
    INFO = 0  IMPROPER INPUT PARAMETERS.

    INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
              DIRECTIONAL DERIVATIVE CONDITION HOLD.

    INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
              IS AT MOST XTOL.

    INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.

    INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.

    INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.

    INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
              THERE MAY NOT BE A STEP WHICH SATISFIES THE
              SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
              TOLERANCES MAY BE TOO SMALL.

NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF CALLS TO FCN.

WA IS A WORK ARRAY OF LENGTH N.

ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
JORGE J. MORE', DAVID J. THUENTE
*************************************************************************/
void mcsrch(ae_int_t n,
     /* Real    */ ae_vector* x,
     double* f,
     /* Real    */ ae_vector* g,
     /* Real    */ ae_vector* s,
     double* stp,
     double stpmax,
     double gtol,
     ae_int_t* info,
     ae_int_t* nfev,
     /* Real    */ ae_vector* wa,
     linminstate* state,
     ae_int_t* stage,
     ae_state *_state)
{
    double v;
    double p5;
    double p66;
    double zero;


    
    /*
     * init
     */
    p5 = 0.5;
    p66 = 0.66;
    state->xtrapf = 4.0;
    zero = 0;
    if( ae_fp_eq(stpmax,0) )
    {
        stpmax = linmin_defstpmax;
    }
    if( ae_fp_less(*stp,linmin_stpmin) )
    {
        *stp = linmin_stpmin;
    }
    if( ae_fp_greater(*stp,stpmax) )
    {
        *stp = stpmax;
    }
    
    /*
     * Main cycle
     */
    for(;;)
    {
        if( *stage==0 )
        {
            
            /*
             * NEXT
             */
            *stage = 2;
            continue;
        }
        if( *stage==2 )
        {
            state->infoc = 1;
            *info = 0;
            
            /*
             *     CHECK THE INPUT PARAMETERS FOR ERRORS.
             */
            if( ae_fp_less(stpmax,linmin_stpmin)&&ae_fp_greater(stpmax,0) )
            {
                *info = 5;
                *stp = 0.0;
                return;
            }
            if( ((((((n<=0||ae_fp_less_eq(*stp,0))||ae_fp_less(linmin_ftol,0))||ae_fp_less(gtol,zero))||ae_fp_less(linmin_xtol,zero))||ae_fp_less(linmin_stpmin,zero))||ae_fp_less(stpmax,linmin_stpmin))||linmin_maxfev<=0 )
            {
                *stage = 0;
                return;
            }
            
            /*
             *     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
             *     AND CHECK THAT S IS A DESCENT DIRECTION.
             */
            v = ae_v_dotproduct(&g->ptr.p_double[0], 1, &s->ptr.p_double[0], 1, ae_v_len(0,n-1));
            state->dginit = v;
            if( ae_fp_greater_eq(state->dginit,0) )
            {
                *stage = 0;
                return;
            }
            
            /*
             *     INITIALIZE LOCAL VARIABLES.
             */
            state->brackt = ae_false;
            state->stage1 = ae_true;
            *nfev = 0;
            state->finit = *f;
            state->dgtest = linmin_ftol*state->dginit;
            state->width = stpmax-linmin_stpmin;
            state->width1 = state->width/p5;
            ae_v_move(&wa->ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,n-1));
            
            /*
             *     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
             *     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
             *     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
             *     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
             *     THE INTERVAL OF UNCERTAINTY.
             *     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
             *     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
             */
            state->stx = 0;
            state->fx = state->finit;
            state->dgx = state->dginit;
            state->sty = 0;
            state->fy = state->finit;
            state->dgy = state->dginit;
            
            /*
             * NEXT
             */
            *stage = 3;
            continue;
        }
        if( *stage==3 )
        {
            
            /*
             *     START OF ITERATION.
             *
             *     SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
             *     TO THE PRESENT INTERVAL OF UNCERTAINTY.
             */
            if( state->brackt )
            {
                if( ae_fp_less(state->stx,state->sty) )
                {
                    state->stmin = state->stx;
                    state->stmax = state->sty;
                }
                else
                {
                    state->stmin = state->sty;
                    state->stmax = state->stx;
                }
            }
            else
            {
                state->stmin = state->stx;
                state->stmax = *stp+state->xtrapf*(*stp-state->stx);
            }
            
            /*
             *        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
             */
            if( ae_fp_greater(*stp,stpmax) )
            {
                *stp = stpmax;
            }
            if( ae_fp_less(*stp,linmin_stpmin) )
            {
                *stp = linmin_stpmin;
            }
            
            /*
             *        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
             *        STP BE THE LOWEST POINT OBTAINED SO FAR.
             */
            if( (((state->brackt&&(ae_fp_less_eq(*stp,state->stmin)||ae_fp_greater_eq(*stp,state->stmax)))||*nfev>=linmin_maxfev-1)||state->infoc==0)||(state->brackt&&ae_fp_less_eq(state->stmax-state->stmin,linmin_xtol*state->stmax)) )
            {
                *stp = state->stx;
            }
            
            /*
             *        EVALUATE THE FUNCTION AND GRADIENT AT STP
             *        AND COMPUTE THE DIRECTIONAL DERIVATIVE.
             */
            ae_v_move(&x->ptr.p_double[0], 1, &wa->ptr.p_double[0], 1, ae_v_len(0,n-1));
            ae_v_addd(&x->ptr.p_double[0], 1, &s->ptr.p_double[0], 1, ae_v_len(0,n-1), *stp);
            
            /*
             * NEXT
             */
            *stage = 4;
            return;
        }
        if( *stage==4 )
        {
            *info = 0;
            *nfev = *nfev+1;
            v = ae_v_dotproduct(&g->ptr.p_double[0], 1, &s->ptr.p_double[0], 1, ae_v_len(0,n-1));
            state->dg = v;
            state->ftest1 = state->finit+*stp*state->dgtest;
            
            /*
             *        TEST FOR CONVERGENCE.
             */
            if( (state->brackt&&(ae_fp_less_eq(*stp,state->stmin)||ae_fp_greater_eq(*stp,state->stmax)))||state->infoc==0 )
            {
                *info = 6;
            }
            if( (ae_fp_eq(*stp,stpmax)&&ae_fp_less_eq(*f,state->ftest1))&&ae_fp_less_eq(state->dg,state->dgtest) )
            {
                *info = 5;
            }
            if( ae_fp_eq(*stp,linmin_stpmin)&&(ae_fp_greater(*f,state->ftest1)||ae_fp_greater_eq(state->dg,state->dgtest)) )
            {
                *info = 4;
            }
            if( *nfev>=linmin_maxfev )
            {
                *info = 3;
            }
            if( state->brackt&&ae_fp_less_eq(state->stmax-state->stmin,linmin_xtol*state->stmax) )
            {
                *info = 2;
            }
            if( ae_fp_less_eq(*f,state->ftest1)&&ae_fp_less_eq(ae_fabs(state->dg, _state),-gtol*state->dginit) )
            {
                *info = 1;
            }
            
            /*
             *        CHECK FOR TERMINATION.
             */
            if( *info!=0 )
            {
                *stage = 0;
                return;
            }
            
            /*
             *        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
             *        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
             */
            if( (state->stage1&&ae_fp_less_eq(*f,state->ftest1))&&ae_fp_greater_eq(state->dg,ae_minreal(linmin_ftol, gtol, _state)*state->dginit) )
            {
                state->stage1 = ae_false;
            }
            
            /*
             *        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
             *        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
             *        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
             *        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
             *        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
             */
            if( (state->stage1&&ae_fp_less_eq(*f,state->fx))&&ae_fp_greater(*f,state->ftest1) )
            {
                
                /*
                 *           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
                 */
                state->fm = *f-*stp*state->dgtest;
                state->fxm = state->fx-state->stx*state->dgtest;
                state->fym = state->fy-state->sty*state->dgtest;
                state->dgm = state->dg-state->dgtest;
                state->dgxm = state->dgx-state->dgtest;
                state->dgym = state->dgy-state->dgtest;
                
                /*
                 *           CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
                 *           AND TO COMPUTE THE NEW STEP.
                 */
                linmin_mcstep(&state->stx, &state->fxm, &state->dgxm, &state->sty, &state->fym, &state->dgym, stp, state->fm, state->dgm, &state->brackt, state->stmin, state->stmax, &state->infoc, _state);
                
                /*
                 *           RESET THE FUNCTION AND GRADIENT VALUES FOR F.
                 */
                state->fx = state->fxm+state->stx*state->dgtest;
                state->fy = state->fym+state->sty*state->dgtest;
                state->dgx = state->dgxm+state->dgtest;
                state->dgy = state->dgym+state->dgtest;
            }
            else
            {
                
                /*
                 *           CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
                 *           AND TO COMPUTE THE NEW STEP.
                 */
                linmin_mcstep(&state->stx, &state->fx, &state->dgx, &state->sty, &state->fy, &state->dgy, stp, *f, state->dg, &state->brackt, state->stmin, state->stmax, &state->infoc, _state);
            }
            
            /*
             *        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
             *        INTERVAL OF UNCERTAINTY.
             */
            if( state->brackt )
            {
                if( ae_fp_greater_eq(ae_fabs(state->sty-state->stx, _state),p66*state->width1) )
                {
                    *stp = state->stx+p5*(state->sty-state->stx);
                }
                state->width1 = state->width;
                state->width = ae_fabs(state->sty-state->stx, _state);
            }
            
            /*
             *  NEXT.
             */
            *stage = 3;
            continue;
        }
    }
}


/*************************************************************************
These functions perform Armijo line search using  at  most  FMAX  function
evaluations.  It  doesn't  enforce  some  kind  of  " sufficient decrease"
criterion - it just tries different Armijo steps and returns optimum found
so far.

Optimization is done using F-rcomm interface:
* ArmijoCreate initializes State structure
  (reusing previously allocated buffers)
* ArmijoIteration is subsequently called
* ArmijoResults returns results

INPUT PARAMETERS:
    N       -   problem size
    X       -   array[N], starting point
    F       -   F(X+S*STP)
    S       -   step direction, S>0
    STP     -   step length
    STPMAX  -   maximum value for STP or zero (if no limit is imposed)
    FMAX    -   maximum number of function evaluations
    State   -   optimization state

  -- ALGLIB --
     Copyright 05.10.2010 by Bochkanov Sergey
*************************************************************************/
void armijocreate(ae_int_t n,
     /* Real    */ ae_vector* x,
     double f,
     /* Real    */ ae_vector* s,
     double stp,
     double stpmax,
     ae_int_t fmax,
     armijostate* state,
     ae_state *_state)
{


    if( state->x.cnt<n )
    {
        ae_vector_set_length(&state->x, n, _state);
    }
    if( state->xbase.cnt<n )
    {
        ae_vector_set_length(&state->xbase, n, _state);
    }
    if( state->s.cnt<n )
    {
        ae_vector_set_length(&state->s, n, _state);
    }
    state->stpmax = stpmax;
    state->fmax = fmax;
    state->stplen = stp;
    state->fcur = f;
    state->n = n;
    ae_v_move(&state->xbase.ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_move(&state->s.ptr.p_double[0], 1, &s->ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_vector_set_length(&state->rstate.ia, 0+1, _state);
    ae_vector_set_length(&state->rstate.ra, 0+1, _state);
    state->rstate.stage = -1;
}


/*************************************************************************
This is rcomm-based search function

  -- ALGLIB --
     Copyright 05.10.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool armijoiteration(armijostate* state, ae_state *) //_state)
{
    double v;
    ae_int_t n;
    ae_bool result;


    
    /*
     * Reverse communication preparations
     * I know it looks ugly, but it works the same way
     * anywhere from C++ to Python.
     *
     * This code initializes locals by:
     * * random values determined during code
     *   generation - on first subroutine call
     * * values from previous call - on subsequent calls
     */
    if( state->rstate.stage>=0 )
    {
        n = state->rstate.ia.ptr.p_int[0];
        v = state->rstate.ra.ptr.p_double[0];
    }
    else
    {
        n = -983;
        v = -989;
    }
    if( state->rstate.stage==0 )
    {
        goto lbl_0;
    }
    if( state->rstate.stage==1 )
    {
        goto lbl_1;
    }
    if( state->rstate.stage==2 )
    {
        goto lbl_2;
    }
    if( state->rstate.stage==3 )
    {
        goto lbl_3;
    }
    
    /*
     * Routine body
     */
    if( (ae_fp_less_eq(state->stplen,0)||ae_fp_less(state->stpmax,0))||state->fmax<2 )
    {
        state->info = 0;
        result = ae_false;
        return result;
    }
    if( ae_fp_less_eq(state->stplen,linmin_stpmin) )
    {
        state->info = 4;
        result = ae_false;
        return result;
    }
    n = state->n;
    state->nfev = 0;
    
    /*
     * We always need F
     */
    state->needf = ae_true;
    
    /*
     * Bound StpLen
     */
    if( ae_fp_greater(state->stplen,state->stpmax)&&ae_fp_neq(state->stpmax,0) )
    {
        state->stplen = state->stpmax;
    }
    
    /*
     * Increase length
     */
    v = state->stplen*linmin_armijofactor;
    if( ae_fp_greater(v,state->stpmax)&&ae_fp_neq(state->stpmax,0) )
    {
        v = state->stpmax;
    }
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_addd(&state->x.ptr.p_double[0], 1, &state->s.ptr.p_double[0], 1, ae_v_len(0,n-1), v);
    state->rstate.stage = 0;
    goto lbl_rcomm;
lbl_0:
    state->nfev = state->nfev+1;
    if( ae_fp_greater_eq(state->f,state->fcur) )
    {
        goto lbl_4;
    }
    state->stplen = v;
    state->fcur = state->f;
lbl_6:
    if( ae_false )
    {
        goto lbl_7;
    }
    
    /*
     * test stopping conditions
     */
    if( state->nfev>=state->fmax )
    {
        state->info = 3;
        result = ae_false;
        return result;
    }
    if( ae_fp_greater_eq(state->stplen,state->stpmax) )
    {
        state->info = 5;
        result = ae_false;
        return result;
    }
    
    /*
     * evaluate F
     */
    v = state->stplen*linmin_armijofactor;
    if( ae_fp_greater(v,state->stpmax)&&ae_fp_neq(state->stpmax,0) )
    {
        v = state->stpmax;
    }
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_addd(&state->x.ptr.p_double[0], 1, &state->s.ptr.p_double[0], 1, ae_v_len(0,n-1), v);
    state->rstate.stage = 1;
    goto lbl_rcomm;
lbl_1:
    state->nfev = state->nfev+1;
    
    /*
     * make decision
     */
    if( ae_fp_less(state->f,state->fcur) )
    {
        state->stplen = v;
        state->fcur = state->f;
    }
    else
    {
        state->info = 1;
        result = ae_false;
        return result;
    }
    goto lbl_6;
lbl_7:
lbl_4:
    
    /*
     * Decrease length
     */
    v = state->stplen/linmin_armijofactor;
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_addd(&state->x.ptr.p_double[0], 1, &state->s.ptr.p_double[0], 1, ae_v_len(0,n-1), v);
    state->rstate.stage = 2;
    goto lbl_rcomm;
lbl_2:
    state->nfev = state->nfev+1;
    if( ae_fp_greater_eq(state->f,state->fcur) )
    {
        goto lbl_8;
    }
    state->stplen = state->stplen/linmin_armijofactor;
    state->fcur = state->f;
lbl_10:
    if( ae_false )
    {
        goto lbl_11;
    }
    
    /*
     * test stopping conditions
     */
    if( state->nfev>=state->fmax )
    {
        state->info = 3;
        result = ae_false;
        return result;
    }
    if( ae_fp_less_eq(state->stplen,linmin_stpmin) )
    {
        state->info = 4;
        result = ae_false;
        return result;
    }
    
    /*
     * evaluate F
     */
    v = state->stplen/linmin_armijofactor;
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_addd(&state->x.ptr.p_double[0], 1, &state->s.ptr.p_double[0], 1, ae_v_len(0,n-1), v);
    state->rstate.stage = 3;
    goto lbl_rcomm;
lbl_3:
    state->nfev = state->nfev+1;
    
    /*
     * make decision
     */
    if( ae_fp_less(state->f,state->fcur) )
    {
        state->stplen = state->stplen/linmin_armijofactor;
        state->fcur = state->f;
    }
    else
    {
        state->info = 1;
        result = ae_false;
        return result;
    }
    goto lbl_10;
lbl_11:
lbl_8:
    
    /*
     * Nothing to be done
     */
    state->info = 1;
    result = ae_false;
    return result;
    
    /*
     * Saving state
     */
lbl_rcomm:
    result = ae_true;
    state->rstate.ia.ptr.p_int[0] = n;
    state->rstate.ra.ptr.p_double[0] = v;
    return result;
}


/*************************************************************************
Results of Armijo search

OUTPUT PARAMETERS:
    INFO    -   on output it is set to one of the return codes:
                * 0     improper input params
                * 1     optimum step is found with at most FMAX evaluations
                * 3     FMAX evaluations were used,
                        X contains optimum found so far
                * 4     step is at lower bound STPMIN
                * 5     step is at upper bound
    STP     -   step length (in case of failure it is still returned)
    F       -   function value (in case of failure it is still returned)

  -- ALGLIB --
     Copyright 05.10.2010 by Bochkanov Sergey
*************************************************************************/
void armijoresults(armijostate* state,
     ae_int_t* info,
     double* stp,
     double* f,
     ae_state *) //_state)
{


    *info = state->info;
    *stp = state->stplen;
    *f = state->fcur;
}


static void linmin_mcstep(double* stx,
     double* fx,
     double* dx,
     double* sty,
     double* fy,
     double* dy,
     double* stp,
     double fp,
     double dp,
     ae_bool* brackt,
     double stmin,
     double stmax,
     ae_int_t* info,
     ae_state *_state)
{
    ae_bool bound;
    double gamma;
    double p;
    double q;
    double r;
    double s;
    double sgnd;
    double stpc;
    double stpf;
    double stpq;
    double theta;


    *info = 0;
    
    /*
     *     CHECK THE INPUT PARAMETERS FOR ERRORS.
     */
    if( ((*brackt&&(ae_fp_less_eq(*stp,ae_minreal(*stx, *sty, _state))||ae_fp_greater_eq(*stp,ae_maxreal(*stx, *sty, _state))))||ae_fp_greater_eq(*dx*(*stp-(*stx)),0))||ae_fp_less(stmax,stmin) )
    {
        return;
    }
    
    /*
     *     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.
     */
    sgnd = dp*(*dx/ae_fabs(*dx, _state));
    
    /*
     *     FIRST CASE. A HIGHER FUNCTION VALUE.
     *     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
     *     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
     *     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.
     */
    if( ae_fp_greater(fp,*fx) )
    {
        *info = 1;
        bound = ae_true;
        theta = 3*(*fx-fp)/(*stp-(*stx))+(*dx)+dp;
        s = ae_maxreal(ae_fabs(theta, _state), ae_maxreal(ae_fabs(*dx, _state), ae_fabs(dp, _state), _state), _state);
        gamma = s*ae_sqrt(ae_sqr(theta/s, _state)-*dx/s*(dp/s), _state);
        if( ae_fp_less(*stp,*stx) )
        {
            gamma = -gamma;
        }
        p = gamma-(*dx)+theta;
        q = gamma-(*dx)+gamma+dp;
        r = p/q;
        stpc = *stx+r*(*stp-(*stx));
        stpq = *stx+*dx/((*fx-fp)/(*stp-(*stx))+(*dx))/2*(*stp-(*stx));
        if( ae_fp_less(ae_fabs(stpc-(*stx), _state),ae_fabs(stpq-(*stx), _state)) )
        {
            stpf = stpc;
        }
        else
        {
            stpf = stpc+(stpq-stpc)/2;
        }
        *brackt = ae_true;
    }
    else
    {
        if( ae_fp_less(sgnd,0) )
        {
            
            /*
             *     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
             *     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
             *     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
             *     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
             */
            *info = 2;
            bound = ae_false;
            theta = 3*(*fx-fp)/(*stp-(*stx))+(*dx)+dp;
            s = ae_maxreal(ae_fabs(theta, _state), ae_maxreal(ae_fabs(*dx, _state), ae_fabs(dp, _state), _state), _state);
            gamma = s*ae_sqrt(ae_sqr(theta/s, _state)-*dx/s*(dp/s), _state);
            if( ae_fp_greater(*stp,*stx) )
            {
                gamma = -gamma;
            }
            p = gamma-dp+theta;
            q = gamma-dp+gamma+(*dx);
            r = p/q;
            stpc = *stp+r*(*stx-(*stp));
            stpq = *stp+dp/(dp-(*dx))*(*stx-(*stp));
            if( ae_fp_greater(ae_fabs(stpc-(*stp), _state),ae_fabs(stpq-(*stp), _state)) )
            {
                stpf = stpc;
            }
            else
            {
                stpf = stpq;
            }
            *brackt = ae_true;
        }
        else
        {
            if( ae_fp_less(ae_fabs(dp, _state),ae_fabs(*dx, _state)) )
            {
                
                /*
                 *     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
                 *     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
                 *     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
                 *     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
                 *     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
                 *     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
                 *     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
                 *     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.
                 */
                *info = 3;
                bound = ae_true;
                theta = 3*(*fx-fp)/(*stp-(*stx))+(*dx)+dp;
                s = ae_maxreal(ae_fabs(theta, _state), ae_maxreal(ae_fabs(*dx, _state), ae_fabs(dp, _state), _state), _state);
                
                /*
                 *        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
                 *        TO INFINITY IN THE DIRECTION OF THE STEP.
                 */
                gamma = s*ae_sqrt(ae_maxreal(0, ae_sqr(theta/s, _state)-*dx/s*(dp/s), _state), _state);
                if( ae_fp_greater(*stp,*stx) )
                {
                    gamma = -gamma;
                }
                p = gamma-dp+theta;
                q = gamma+(*dx-dp)+gamma;
                r = p/q;
                if( ae_fp_less(r,0)&&ae_fp_neq(gamma,0) )
                {
                    stpc = *stp+r*(*stx-(*stp));
                }
                else
                {
                    if( ae_fp_greater(*stp,*stx) )
                    {
                        stpc = stmax;
                    }
                    else
                    {
                        stpc = stmin;
                    }
                }
                stpq = *stp+dp/(dp-(*dx))*(*stx-(*stp));
                if( *brackt )
                {
                    if( ae_fp_less(ae_fabs(*stp-stpc, _state),ae_fabs(*stp-stpq, _state)) )
                    {
                        stpf = stpc;
                    }
                    else
                    {
                        stpf = stpq;
                    }
                }
                else
                {
                    if( ae_fp_greater(ae_fabs(*stp-stpc, _state),ae_fabs(*stp-stpq, _state)) )
                    {
                        stpf = stpc;
                    }
                    else
                    {
                        stpf = stpq;
                    }
                }
            }
            else
            {
                
                /*
                 *     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
                 *     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
                 *     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
                 *     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
                 */
                *info = 4;
                bound = ae_false;
                if( *brackt )
                {
                    theta = 3*(fp-(*fy))/(*sty-(*stp))+(*dy)+dp;
                    s = ae_maxreal(ae_fabs(theta, _state), ae_maxreal(ae_fabs(*dy, _state), ae_fabs(dp, _state), _state), _state);
                    gamma = s*ae_sqrt(ae_sqr(theta/s, _state)-*dy/s*(dp/s), _state);
                    if( ae_fp_greater(*stp,*sty) )
                    {
                        gamma = -gamma;
                    }
                    p = gamma-dp+theta;
                    q = gamma-dp+gamma+(*dy);
                    r = p/q;
                    stpc = *stp+r*(*sty-(*stp));
                    stpf = stpc;
                }
                else
                {
                    if( ae_fp_greater(*stp,*stx) )
                    {
                        stpf = stmax;
                    }
                    else
                    {
                        stpf = stmin;
                    }
                }
            }
        }
    }
    
    /*
     *     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
     *     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.
     */
    if( ae_fp_greater(fp,*fx) )
    {
        *sty = *stp;
        *fy = fp;
        *dy = dp;
    }
    else
    {
        if( ae_fp_less(sgnd,0.0) )
        {
            *sty = *stx;
            *fy = *fx;
            *dy = *dx;
        }
        *stx = *stp;
        *fx = fp;
        *dx = dp;
    }
    
    /*
     *     COMPUTE THE NEW STEP AND SAFEGUARD IT.
     */
    stpf = ae_minreal(stmax, stpf, _state);
    stpf = ae_maxreal(stmin, stpf, _state);
    *stp = stpf;
    if( *brackt&&bound )
    {
        if( ae_fp_greater(*sty,*stx) )
        {
            *stp = ae_minreal(*stx+0.66*(*sty-(*stx)), *stp, _state);
        }
        else
        {
            *stp = ae_maxreal(*stx+0.66*(*sty-(*stx)), *stp, _state);
        }
    }
}


ae_bool _linminstate_init(linminstate* /*p*/, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    return ae_true;
}


ae_bool _linminstate_init_copy(linminstate* dst, linminstate* src, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    dst->brackt = src->brackt;
    dst->stage1 = src->stage1;
    dst->infoc = src->infoc;
    dst->dg = src->dg;
    dst->dgm = src->dgm;
    dst->dginit = src->dginit;
    dst->dgtest = src->dgtest;
    dst->dgx = src->dgx;
    dst->dgxm = src->dgxm;
    dst->dgy = src->dgy;
    dst->dgym = src->dgym;
    dst->finit = src->finit;
    dst->ftest1 = src->ftest1;
    dst->fm = src->fm;
    dst->fx = src->fx;
    dst->fxm = src->fxm;
    dst->fy = src->fy;
    dst->fym = src->fym;
    dst->stx = src->stx;
    dst->sty = src->sty;
    dst->stmin = src->stmin;
    dst->stmax = src->stmax;
    dst->width = src->width;
    dst->width1 = src->width1;
    dst->xtrapf = src->xtrapf;
    return ae_true;
}


void _linminstate_clear(linminstate* ) //p)
{
}


ae_bool _armijostate_init(armijostate* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->x, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xbase, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->s, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !_rcommstate_init(&p->rstate, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _armijostate_init_copy(armijostate* dst, armijostate* src, ae_state *_state, ae_bool make_automatic)
{
    dst->needf = src->needf;
    if( !ae_vector_init_copy(&dst->x, &src->x, _state, make_automatic) )
        return ae_false;
    dst->f = src->f;
    dst->n = src->n;
    if( !ae_vector_init_copy(&dst->xbase, &src->xbase, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->s, &src->s, _state, make_automatic) )
        return ae_false;
    dst->stplen = src->stplen;
    dst->fcur = src->fcur;
    dst->stpmax = src->stpmax;
    dst->fmax = src->fmax;
    dst->nfev = src->nfev;
    dst->info = src->info;
    if( !_rcommstate_init_copy(&dst->rstate, &src->rstate, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _armijostate_clear(armijostate* p)
{
    ae_vector_clear(&p->x);
    ae_vector_clear(&p->xbase);
    ae_vector_clear(&p->s);
    _rcommstate_clear(&p->rstate);
}




/*************************************************************************
This subroutine generates FFT plan - a decomposition of a N-length FFT to
the more simpler operations. Plan consists of the root entry and the child
entries.

Subroutine parameters:
    N               task size
    
Output parameters:
    Plan            plan

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
void ftbasegeneratecomplexfftplan(ae_int_t n,
     ftplan* plan,
     ae_state *_state)
{
    ae_int_t planarraysize;
    ae_int_t plansize;
    ae_int_t precomputedsize;
    ae_int_t tmpmemsize;
    ae_int_t stackmemsize;
    ae_int_t stackptr;

    _ftplan_clear(plan);

    planarraysize = 1;
    plansize = 0;
    precomputedsize = 0;
    stackmemsize = 0;
    stackptr = 0;
    tmpmemsize = 2*n;
    ae_vector_set_length(&plan->plan, planarraysize, _state);
    ftbase_ftbasegenerateplanrec(n, ftbase_ftbasecffttask, plan, &plansize, &precomputedsize, &planarraysize, &tmpmemsize, &stackmemsize, stackptr, _state);
    ae_assert(stackptr==0, "Internal error in FTBaseGenerateComplexFFTPlan: stack ptr!", _state);
    ae_vector_set_length(&plan->stackbuf, ae_maxint(stackmemsize, 1, _state), _state);
    ae_vector_set_length(&plan->tmpbuf, ae_maxint(tmpmemsize, 1, _state), _state);
    ae_vector_set_length(&plan->precomputed, ae_maxint(precomputedsize, 1, _state), _state);
    stackptr = 0;
    ftbase_ftbaseprecomputeplanrec(plan, 0, stackptr, _state);
    ae_assert(stackptr==0, "Internal error in FTBaseGenerateComplexFFTPlan: stack ptr!", _state);
}


/*************************************************************************
Generates real FFT plan
*************************************************************************/
void ftbasegeneraterealfftplan(ae_int_t n, ftplan* plan, ae_state *_state)
{
    ae_int_t planarraysize;
    ae_int_t plansize;
    ae_int_t precomputedsize;
    ae_int_t tmpmemsize;
    ae_int_t stackmemsize;
    ae_int_t stackptr;

    _ftplan_clear(plan);

    planarraysize = 1;
    plansize = 0;
    precomputedsize = 0;
    stackmemsize = 0;
    stackptr = 0;
    tmpmemsize = 2*n;
    ae_vector_set_length(&plan->plan, planarraysize, _state);
    ftbase_ftbasegenerateplanrec(n, ftbase_ftbaserffttask, plan, &plansize, &precomputedsize, &planarraysize, &tmpmemsize, &stackmemsize, stackptr, _state);
    ae_assert(stackptr==0, "Internal error in FTBaseGenerateRealFFTPlan: stack ptr!", _state);
    ae_vector_set_length(&plan->stackbuf, ae_maxint(stackmemsize, 1, _state), _state);
    ae_vector_set_length(&plan->tmpbuf, ae_maxint(tmpmemsize, 1, _state), _state);
    ae_vector_set_length(&plan->precomputed, ae_maxint(precomputedsize, 1, _state), _state);
    stackptr = 0;
    ftbase_ftbaseprecomputeplanrec(plan, 0, stackptr, _state);
    ae_assert(stackptr==0, "Internal error in FTBaseGenerateRealFFTPlan: stack ptr!", _state);
}


/*************************************************************************
Generates real FHT plan
*************************************************************************/
void ftbasegeneraterealfhtplan(ae_int_t n, ftplan* plan, ae_state *_state)
{
    ae_int_t planarraysize;
    ae_int_t plansize;
    ae_int_t precomputedsize;
    ae_int_t tmpmemsize;
    ae_int_t stackmemsize;
    ae_int_t stackptr;

    _ftplan_clear(plan);

    planarraysize = 1;
    plansize = 0;
    precomputedsize = 0;
    stackmemsize = 0;
    stackptr = 0;
    tmpmemsize = n;
    ae_vector_set_length(&plan->plan, planarraysize, _state);
    ftbase_ftbasegenerateplanrec(n, ftbase_ftbaserfhttask, plan, &plansize, &precomputedsize, &planarraysize, &tmpmemsize, &stackmemsize, stackptr, _state);
    ae_assert(stackptr==0, "Internal error in FTBaseGenerateRealFHTPlan: stack ptr!", _state);
    ae_vector_set_length(&plan->stackbuf, ae_maxint(stackmemsize, 1, _state), _state);
    ae_vector_set_length(&plan->tmpbuf, ae_maxint(tmpmemsize, 1, _state), _state);
    ae_vector_set_length(&plan->precomputed, ae_maxint(precomputedsize, 1, _state), _state);
    stackptr = 0;
    ftbase_ftbaseprecomputeplanrec(plan, 0, stackptr, _state);
    ae_assert(stackptr==0, "Internal error in FTBaseGenerateRealFHTPlan: stack ptr!", _state);
}


/*************************************************************************
This subroutine executes FFT/FHT plan.

If Plan is a:
* complex FFT plan  -   sizeof(A)=2*N,
                        A contains interleaved real/imaginary values
* real FFT plan     -   sizeof(A)=2*N,
                        A contains real values interleaved with zeros
* real FHT plan     -   sizeof(A)=2*N,
                        A contains real values interleaved with zeros

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
void ftbaseexecuteplan(/* Real    */ ae_vector* a,
     ae_int_t aoffset,
     ae_int_t , //n,
     ftplan* plan,
     ae_state *_state)
{
    ae_int_t stackptr;


    stackptr = 0;
    ftbaseexecuteplanrec(a, aoffset, plan, 0, stackptr, _state);
}


/*************************************************************************
Recurrent subroutine for the FTBaseExecutePlan

Parameters:
    A           FFT'ed array
    AOffset     offset of the FFT'ed part (distance is measured in doubles)

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
void ftbaseexecuteplanrec(/* Real    */ ae_vector* a,
     ae_int_t aoffset,
     ftplan* plan,
     ae_int_t entryoffset,
     ae_int_t stackptr,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t n1;
    ae_int_t n2;
    ae_int_t n;
    ae_int_t m;
    ae_int_t offs;
    ae_int_t offs1;
    ae_int_t offs2;
    ae_int_t offsa;
    ae_int_t offsb;
    ae_int_t offsp;
    double hk;
    double hnk;
    double x;
    double y;
    double bx;
    double by;
    ae_vector emptyarray;
    double a0x;
    double a0y;
    double a1x;
    double a1y;
    double a2x;
    double a2y;
    double a3x;
    double a3y;
    double v0;
    double v1;
    double v2;
    double v3;
    double t1x;
    double t1y;
    double t2x;
    double t2y;
    double t3x;
    double t3y;
    double t4x;
    double t4y;
    double t5x;
    double t5y;
    double m1x;
    double m1y;
    double m2x;
    double m2y;
    double m3x;
    double m3y;
    double m4x;
    double m4y;
    double m5x;
    double m5y;
    double s1x;
    double s1y;
    double s2x;
    double s2y;
    double s3x;
    double s3y;
    double s4x;
    double s4y;
    double s5x;
    double s5y;
    double c1;
    double c2;
    double c3;
    double c4;
    double c5;
    ae_vector tmp;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&emptyarray, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);

    if( plan->plan.ptr.p_int[entryoffset+3]==ftbase_fftemptyplan )
    {
        ae_frame_leave(_state);
        return;
    }
    if( plan->plan.ptr.p_int[entryoffset+3]==ftbase_fftcooleytukeyplan )
    {
        
        /*
         * Cooley-Tukey plan
         * * transposition
         * * row-wise FFT
         * * twiddle factors:
         *   - TwBase is a basis twiddle factor for I=1, J=1
         *   - TwRow is a twiddle factor for a second element in a row (J=1)
         *   - Tw is a twiddle factor for a current element
         * * transposition again
         * * row-wise FFT again
         */
        n1 = plan->plan.ptr.p_int[entryoffset+1];
        n2 = plan->plan.ptr.p_int[entryoffset+2];
        ftbase_internalcomplexlintranspose(a, n1, n2, aoffset, &plan->tmpbuf, _state);
        for(i=0; i<=n2-1; i++)
        {
            ftbaseexecuteplanrec(a, aoffset+i*n1*2, plan, plan->plan.ptr.p_int[entryoffset+5], stackptr, _state);
        }
        ftbase_ffttwcalc(a, aoffset, n1, n2, _state);
        ftbase_internalcomplexlintranspose(a, n2, n1, aoffset, &plan->tmpbuf, _state);
        for(i=0; i<=n1-1; i++)
        {
            ftbaseexecuteplanrec(a, aoffset+i*n2*2, plan, plan->plan.ptr.p_int[entryoffset+6], stackptr, _state);
        }
        ftbase_internalcomplexlintranspose(a, n1, n2, aoffset, &plan->tmpbuf, _state);
        ae_frame_leave(_state);
        return;
    }
    if( plan->plan.ptr.p_int[entryoffset+3]==ftbase_fftrealcooleytukeyplan )
    {
        
        /*
         * Cooley-Tukey plan
         * * transposition
         * * row-wise FFT
         * * twiddle factors:
         *   - TwBase is a basis twiddle factor for I=1, J=1
         *   - TwRow is a twiddle factor for a second element in a row (J=1)
         *   - Tw is a twiddle factor for a current element
         * * transposition again
         * * row-wise FFT again
         */
        n1 = plan->plan.ptr.p_int[entryoffset+1];
        n2 = plan->plan.ptr.p_int[entryoffset+2];
        ftbase_internalcomplexlintranspose(a, n2, n1, aoffset, &plan->tmpbuf, _state);
        for(i=0; i<=n1/2-1; i++)
        {
            
            /*
             * pack two adjacent smaller real FFT's together,
             * make one complex FFT,
             * unpack result
             */
            offs = aoffset+2*i*n2*2;
            for(k=0; k<=n2-1; k++)
            {
                a->ptr.p_double[offs+2*k+1] = a->ptr.p_double[offs+2*n2+2*k+0];
            }
            ftbaseexecuteplanrec(a, offs, plan, plan->plan.ptr.p_int[entryoffset+6], stackptr, _state);
            plan->tmpbuf.ptr.p_double[0] = a->ptr.p_double[offs+0];
            plan->tmpbuf.ptr.p_double[1] = 0;
            plan->tmpbuf.ptr.p_double[2*n2+0] = a->ptr.p_double[offs+1];
            plan->tmpbuf.ptr.p_double[2*n2+1] = 0;
            for(k=1; k<=n2-1; k++)
            {
                offs1 = 2*k;
                offs2 = 2*n2+2*k;
                hk = a->ptr.p_double[offs+2*k+0];
                hnk = a->ptr.p_double[offs+2*(n2-k)+0];
                plan->tmpbuf.ptr.p_double[offs1+0] = 0.5*(hk+hnk);
                plan->tmpbuf.ptr.p_double[offs2+1] = -0.5*(hk-hnk);
                hk = a->ptr.p_double[offs+2*k+1];
                hnk = a->ptr.p_double[offs+2*(n2-k)+1];
                plan->tmpbuf.ptr.p_double[offs2+0] = 0.5*(hk+hnk);
                plan->tmpbuf.ptr.p_double[offs1+1] = 0.5*(hk-hnk);
            }
            ae_v_move(&a->ptr.p_double[offs], 1, &plan->tmpbuf.ptr.p_double[0], 1, ae_v_len(offs,offs+2*n2*2-1));
        }
        if( n1%2!=0 )
        {
            ftbaseexecuteplanrec(a, aoffset+(n1-1)*n2*2, plan, plan->plan.ptr.p_int[entryoffset+6], stackptr, _state);
        }
        ftbase_ffttwcalc(a, aoffset, n2, n1, _state);
        ftbase_internalcomplexlintranspose(a, n1, n2, aoffset, &plan->tmpbuf, _state);
        for(i=0; i<=n2-1; i++)
        {
            ftbaseexecuteplanrec(a, aoffset+i*n1*2, plan, plan->plan.ptr.p_int[entryoffset+5], stackptr, _state);
        }
        ftbase_internalcomplexlintranspose(a, n2, n1, aoffset, &plan->tmpbuf, _state);
        ae_frame_leave(_state);
        return;
    }
    if( plan->plan.ptr.p_int[entryoffset+3]==ftbase_fhtcooleytukeyplan )
    {
        
        /*
         * Cooley-Tukey FHT plan:
         * * transpose                    \
         * * smaller FHT's                |
         * * pre-process                  |
         * * multiply by twiddle factors  | corresponds to multiplication by H1
         * * post-process                 |
         * * transpose again              /
         * * multiply by H2 (smaller FHT's)
         * * final transposition
         *
         * For more details see Vitezslav Vesely, "Fast algorithms
         * of Fourier and Hartley transform and their implementation in MATLAB",
         * page 31.
         */
        n1 = plan->plan.ptr.p_int[entryoffset+1];
        n2 = plan->plan.ptr.p_int[entryoffset+2];
        n = n1*n2;
        ftbase_internalreallintranspose(a, n1, n2, aoffset, &plan->tmpbuf, _state);
        for(i=0; i<=n2-1; i++)
        {
            ftbaseexecuteplanrec(a, aoffset+i*n1, plan, plan->plan.ptr.p_int[entryoffset+5], stackptr, _state);
        }
        for(i=0; i<=n2-1; i++)
        {
            for(j=0; j<=n1-1; j++)
            {
                offsa = aoffset+i*n1;
                hk = a->ptr.p_double[offsa+j];
                hnk = a->ptr.p_double[offsa+(n1-j)%n1];
                offs = 2*(i*n1+j);
                plan->tmpbuf.ptr.p_double[offs+0] = -0.5*(hnk-hk);
                plan->tmpbuf.ptr.p_double[offs+1] = 0.5*(hk+hnk);
            }
        }
        ftbase_ffttwcalc(&plan->tmpbuf, 0, n1, n2, _state);
        for(j=0; j<=n1-1; j++)
        {
            a->ptr.p_double[aoffset+j] = plan->tmpbuf.ptr.p_double[2*j+0]+plan->tmpbuf.ptr.p_double[2*j+1];
        }
        if( n2%2==0 )
        {
            offs = 2*(n2/2)*n1;
            offsa = aoffset+n2/2*n1;
            for(j=0; j<=n1-1; j++)
            {
                a->ptr.p_double[offsa+j] = plan->tmpbuf.ptr.p_double[offs+2*j+0]+plan->tmpbuf.ptr.p_double[offs+2*j+1];
            }
        }
        for(i=1; i<=(n2+1)/2-1; i++)
        {
            offs = 2*i*n1;
            offs2 = 2*(n2-i)*n1;
            offsa = aoffset+i*n1;
            for(j=0; j<=n1-1; j++)
            {
                a->ptr.p_double[offsa+j] = plan->tmpbuf.ptr.p_double[offs+2*j+1]+plan->tmpbuf.ptr.p_double[offs2+2*j+0];
            }
            offsa = aoffset+(n2-i)*n1;
            for(j=0; j<=n1-1; j++)
            {
                a->ptr.p_double[offsa+j] = plan->tmpbuf.ptr.p_double[offs+2*j+0]+plan->tmpbuf.ptr.p_double[offs2+2*j+1];
            }
        }
        ftbase_internalreallintranspose(a, n2, n1, aoffset, &plan->tmpbuf, _state);
        for(i=0; i<=n1-1; i++)
        {
            ftbaseexecuteplanrec(a, aoffset+i*n2, plan, plan->plan.ptr.p_int[entryoffset+6], stackptr, _state);
        }
        ftbase_internalreallintranspose(a, n1, n2, aoffset, &plan->tmpbuf, _state);
        ae_frame_leave(_state);
        return;
    }
    if( plan->plan.ptr.p_int[entryoffset+3]==ftbase_fhtn2plan )
    {
        
        /*
         * Cooley-Tukey FHT plan
         */
        n1 = plan->plan.ptr.p_int[entryoffset+1];
        n2 = plan->plan.ptr.p_int[entryoffset+2];
        n = n1*n2;
        ftbase_reffht(a, n, aoffset, _state);
        ae_frame_leave(_state);
        return;
    }
    if( plan->plan.ptr.p_int[entryoffset+3]==ftbase_fftcodeletplan )
    {
        n1 = plan->plan.ptr.p_int[entryoffset+1];
        n2 = plan->plan.ptr.p_int[entryoffset+2];
        n = n1*n2;
        if( n==2 )
        {
            a0x = a->ptr.p_double[aoffset+0];
            a0y = a->ptr.p_double[aoffset+1];
            a1x = a->ptr.p_double[aoffset+2];
            a1y = a->ptr.p_double[aoffset+3];
            v0 = a0x+a1x;
            v1 = a0y+a1y;
            v2 = a0x-a1x;
            v3 = a0y-a1y;
            a->ptr.p_double[aoffset+0] = v0;
            a->ptr.p_double[aoffset+1] = v1;
            a->ptr.p_double[aoffset+2] = v2;
            a->ptr.p_double[aoffset+3] = v3;
            ae_frame_leave(_state);
            return;
        }
        if( n==3 )
        {
            offs = plan->plan.ptr.p_int[entryoffset+7];
            c1 = plan->precomputed.ptr.p_double[offs+0];
            c2 = plan->precomputed.ptr.p_double[offs+1];
            a0x = a->ptr.p_double[aoffset+0];
            a0y = a->ptr.p_double[aoffset+1];
            a1x = a->ptr.p_double[aoffset+2];
            a1y = a->ptr.p_double[aoffset+3];
            a2x = a->ptr.p_double[aoffset+4];
            a2y = a->ptr.p_double[aoffset+5];
            t1x = a1x+a2x;
            t1y = a1y+a2y;
            a0x = a0x+t1x;
            a0y = a0y+t1y;
            m1x = c1*t1x;
            m1y = c1*t1y;
            m2x = c2*(a1y-a2y);
            m2y = c2*(a2x-a1x);
            s1x = a0x+m1x;
            s1y = a0y+m1y;
            a1x = s1x+m2x;
            a1y = s1y+m2y;
            a2x = s1x-m2x;
            a2y = s1y-m2y;
            a->ptr.p_double[aoffset+0] = a0x;
            a->ptr.p_double[aoffset+1] = a0y;
            a->ptr.p_double[aoffset+2] = a1x;
            a->ptr.p_double[aoffset+3] = a1y;
            a->ptr.p_double[aoffset+4] = a2x;
            a->ptr.p_double[aoffset+5] = a2y;
            ae_frame_leave(_state);
            return;
        }
        if( n==4 )
        {
            a0x = a->ptr.p_double[aoffset+0];
            a0y = a->ptr.p_double[aoffset+1];
            a1x = a->ptr.p_double[aoffset+2];
            a1y = a->ptr.p_double[aoffset+3];
            a2x = a->ptr.p_double[aoffset+4];
            a2y = a->ptr.p_double[aoffset+5];
            a3x = a->ptr.p_double[aoffset+6];
            a3y = a->ptr.p_double[aoffset+7];
            t1x = a0x+a2x;
            t1y = a0y+a2y;
            t2x = a1x+a3x;
            t2y = a1y+a3y;
            m2x = a0x-a2x;
            m2y = a0y-a2y;
            m3x = a1y-a3y;
            m3y = a3x-a1x;
            a->ptr.p_double[aoffset+0] = t1x+t2x;
            a->ptr.p_double[aoffset+1] = t1y+t2y;
            a->ptr.p_double[aoffset+4] = t1x-t2x;
            a->ptr.p_double[aoffset+5] = t1y-t2y;
            a->ptr.p_double[aoffset+2] = m2x+m3x;
            a->ptr.p_double[aoffset+3] = m2y+m3y;
            a->ptr.p_double[aoffset+6] = m2x-m3x;
            a->ptr.p_double[aoffset+7] = m2y-m3y;
            ae_frame_leave(_state);
            return;
        }
        if( n==5 )
        {
            offs = plan->plan.ptr.p_int[entryoffset+7];
            c1 = plan->precomputed.ptr.p_double[offs+0];
            c2 = plan->precomputed.ptr.p_double[offs+1];
            c3 = plan->precomputed.ptr.p_double[offs+2];
            c4 = plan->precomputed.ptr.p_double[offs+3];
            c5 = plan->precomputed.ptr.p_double[offs+4];
            t1x = a->ptr.p_double[aoffset+2]+a->ptr.p_double[aoffset+8];
            t1y = a->ptr.p_double[aoffset+3]+a->ptr.p_double[aoffset+9];
            t2x = a->ptr.p_double[aoffset+4]+a->ptr.p_double[aoffset+6];
            t2y = a->ptr.p_double[aoffset+5]+a->ptr.p_double[aoffset+7];
            t3x = a->ptr.p_double[aoffset+2]-a->ptr.p_double[aoffset+8];
            t3y = a->ptr.p_double[aoffset+3]-a->ptr.p_double[aoffset+9];
            t4x = a->ptr.p_double[aoffset+6]-a->ptr.p_double[aoffset+4];
            t4y = a->ptr.p_double[aoffset+7]-a->ptr.p_double[aoffset+5];
            t5x = t1x+t2x;
            t5y = t1y+t2y;
            a->ptr.p_double[aoffset+0] = a->ptr.p_double[aoffset+0]+t5x;
            a->ptr.p_double[aoffset+1] = a->ptr.p_double[aoffset+1]+t5y;
            m1x = c1*t5x;
            m1y = c1*t5y;
            m2x = c2*(t1x-t2x);
            m2y = c2*(t1y-t2y);
            m3x = -c3*(t3y+t4y);
            m3y = c3*(t3x+t4x);
            m4x = -c4*t4y;
            m4y = c4*t4x;
            m5x = -c5*t3y;
            m5y = c5*t3x;
            s3x = m3x-m4x;
            s3y = m3y-m4y;
            s5x = m3x+m5x;
            s5y = m3y+m5y;
            s1x = a->ptr.p_double[aoffset+0]+m1x;
            s1y = a->ptr.p_double[aoffset+1]+m1y;
            s2x = s1x+m2x;
            s2y = s1y+m2y;
            s4x = s1x-m2x;
            s4y = s1y-m2y;
            a->ptr.p_double[aoffset+2] = s2x+s3x;
            a->ptr.p_double[aoffset+3] = s2y+s3y;
            a->ptr.p_double[aoffset+4] = s4x+s5x;
            a->ptr.p_double[aoffset+5] = s4y+s5y;
            a->ptr.p_double[aoffset+6] = s4x-s5x;
            a->ptr.p_double[aoffset+7] = s4y-s5y;
            a->ptr.p_double[aoffset+8] = s2x-s3x;
            a->ptr.p_double[aoffset+9] = s2y-s3y;
            ae_frame_leave(_state);
            return;
        }
    }
    if( plan->plan.ptr.p_int[entryoffset+3]==ftbase_fhtcodeletplan )
    {
        n1 = plan->plan.ptr.p_int[entryoffset+1];
        n2 = plan->plan.ptr.p_int[entryoffset+2];
        n = n1*n2;
        if( n==2 )
        {
            a0x = a->ptr.p_double[aoffset+0];
            a1x = a->ptr.p_double[aoffset+1];
            a->ptr.p_double[aoffset+0] = a0x+a1x;
            a->ptr.p_double[aoffset+1] = a0x-a1x;
            ae_frame_leave(_state);
            return;
        }
        if( n==3 )
        {
            offs = plan->plan.ptr.p_int[entryoffset+7];
            c1 = plan->precomputed.ptr.p_double[offs+0];
            c2 = plan->precomputed.ptr.p_double[offs+1];
            a0x = a->ptr.p_double[aoffset+0];
            a1x = a->ptr.p_double[aoffset+1];
            a2x = a->ptr.p_double[aoffset+2];
            t1x = a1x+a2x;
            a0x = a0x+t1x;
            m1x = c1*t1x;
            m2y = c2*(a2x-a1x);
            s1x = a0x+m1x;
            a->ptr.p_double[aoffset+0] = a0x;
            a->ptr.p_double[aoffset+1] = s1x-m2y;
            a->ptr.p_double[aoffset+2] = s1x+m2y;
            ae_frame_leave(_state);
            return;
        }
        if( n==4 )
        {
            a0x = a->ptr.p_double[aoffset+0];
            a1x = a->ptr.p_double[aoffset+1];
            a2x = a->ptr.p_double[aoffset+2];
            a3x = a->ptr.p_double[aoffset+3];
            t1x = a0x+a2x;
            t2x = a1x+a3x;
            m2x = a0x-a2x;
            m3y = a3x-a1x;
            a->ptr.p_double[aoffset+0] = t1x+t2x;
            a->ptr.p_double[aoffset+1] = m2x-m3y;
            a->ptr.p_double[aoffset+2] = t1x-t2x;
            a->ptr.p_double[aoffset+3] = m2x+m3y;
            ae_frame_leave(_state);
            return;
        }
        if( n==5 )
        {
            offs = plan->plan.ptr.p_int[entryoffset+7];
            c1 = plan->precomputed.ptr.p_double[offs+0];
            c2 = plan->precomputed.ptr.p_double[offs+1];
            c3 = plan->precomputed.ptr.p_double[offs+2];
            c4 = plan->precomputed.ptr.p_double[offs+3];
            c5 = plan->precomputed.ptr.p_double[offs+4];
            t1x = a->ptr.p_double[aoffset+1]+a->ptr.p_double[aoffset+4];
            t2x = a->ptr.p_double[aoffset+2]+a->ptr.p_double[aoffset+3];
            t3x = a->ptr.p_double[aoffset+1]-a->ptr.p_double[aoffset+4];
            t4x = a->ptr.p_double[aoffset+3]-a->ptr.p_double[aoffset+2];
            t5x = t1x+t2x;
            v0 = a->ptr.p_double[aoffset+0]+t5x;
            a->ptr.p_double[aoffset+0] = v0;
            m2x = c2*(t1x-t2x);
            m3y = c3*(t3x+t4x);
            s3y = m3y-c4*t4x;
            s5y = m3y+c5*t3x;
            s1x = v0+c1*t5x;
            s2x = s1x+m2x;
            s4x = s1x-m2x;
            a->ptr.p_double[aoffset+1] = s2x-s3y;
            a->ptr.p_double[aoffset+2] = s4x-s5y;
            a->ptr.p_double[aoffset+3] = s4x+s5y;
            a->ptr.p_double[aoffset+4] = s2x+s3y;
            ae_frame_leave(_state);
            return;
        }
    }
    if( plan->plan.ptr.p_int[entryoffset+3]==ftbase_fftbluesteinplan )
    {
        
        /*
         * Bluestein plan:
         * 1. multiply by precomputed coefficients
         * 2. make convolution: forward FFT, multiplication by precomputed FFT
         *    and backward FFT. backward FFT is represented as
         *
         *        invfft(x) = fft(x')'/M
         *
         *    for performance reasons reduction of inverse FFT to
         *    forward FFT is merged with multiplication of FFT components
         *    and last stage of Bluestein's transformation.
         * 3. post-multiplication by Bluestein factors
         */
        n = plan->plan.ptr.p_int[entryoffset+1];
        m = plan->plan.ptr.p_int[entryoffset+4];
        offs = plan->plan.ptr.p_int[entryoffset+7];
        for(i=stackptr+2*n; i<=stackptr+2*m-1; i++)
        {
            plan->stackbuf.ptr.p_double[i] = 0;
        }
        offsp = offs+2*m;
        offsa = aoffset;
        offsb = stackptr;
        for(i=0; i<=n-1; i++)
        {
            bx = plan->precomputed.ptr.p_double[offsp+0];
            by = plan->precomputed.ptr.p_double[offsp+1];
            x = a->ptr.p_double[offsa+0];
            y = a->ptr.p_double[offsa+1];
            plan->stackbuf.ptr.p_double[offsb+0] = x*bx-y*(-by);
            plan->stackbuf.ptr.p_double[offsb+1] = x*(-by)+y*bx;
            offsp = offsp+2;
            offsa = offsa+2;
            offsb = offsb+2;
        }
        ftbaseexecuteplanrec(&plan->stackbuf, stackptr, plan, plan->plan.ptr.p_int[entryoffset+5], stackptr+2*2*m, _state);
        offsb = stackptr;
        offsp = offs;
        for(i=0; i<=m-1; i++)
        {
            x = plan->stackbuf.ptr.p_double[offsb+0];
            y = plan->stackbuf.ptr.p_double[offsb+1];
            bx = plan->precomputed.ptr.p_double[offsp+0];
            by = plan->precomputed.ptr.p_double[offsp+1];
            plan->stackbuf.ptr.p_double[offsb+0] = x*bx-y*by;
            plan->stackbuf.ptr.p_double[offsb+1] = -(x*by+y*bx);
            offsb = offsb+2;
            offsp = offsp+2;
        }
        ftbaseexecuteplanrec(&plan->stackbuf, stackptr, plan, plan->plan.ptr.p_int[entryoffset+5], stackptr+2*2*m, _state);
        offsb = stackptr;
        offsp = offs+2*m;
        offsa = aoffset;
        for(i=0; i<=n-1; i++)
        {
            x = plan->stackbuf.ptr.p_double[offsb+0]/m;
            y = -plan->stackbuf.ptr.p_double[offsb+1]/m;
            bx = plan->precomputed.ptr.p_double[offsp+0];
            by = plan->precomputed.ptr.p_double[offsp+1];
            a->ptr.p_double[offsa+0] = x*bx-y*(-by);
            a->ptr.p_double[offsa+1] = x*(-by)+y*bx;
            offsp = offsp+2;
            offsa = offsa+2;
            offsb = offsb+2;
        }
        ae_frame_leave(_state);
        return;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Returns good factorization N=N1*N2.

Usually N1<=N2 (but not always - small N's may be exception).
if N1<>1 then N2<>1.

Factorization is chosen depending on task type and codelets we have.

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
void ftbasefactorize(ae_int_t n,
     ae_int_t , //tasktype,
     ae_int_t* n1,
     ae_int_t* n2,
     ae_state *) //_state)
{
    ae_int_t j;

    *n1 = 0;
    *n2 = 0;

    *n1 = 0;
    *n2 = 0;
    
    /*
     * try to find good codelet
     */
    if( *n1*(*n2)!=n )
    {
        for(j=ftbase_ftbasecodeletrecommended; j>=2; j--)
        {
            if( n%j==0 )
            {
                *n1 = j;
                *n2 = n/j;
                break;
            }
        }
    }
    
    /*
     * try to factorize N
     */
    if( *n1*(*n2)!=n )
    {
        for(j=ftbase_ftbasecodeletrecommended+1; j<=n-1; j++)
        {
            if( n%j==0 )
            {
                *n1 = j;
                *n2 = n/j;
                break;
            }
        }
    }
    
    /*
     * looks like N is prime :(
     */
    if( *n1*(*n2)!=n )
    {
        *n1 = 1;
        *n2 = n;
    }
    
    /*
     * normalize
     */
    if( *n2==1&&*n1!=1 )
    {
        *n2 = *n1;
        *n1 = 1;
    }
}


/*************************************************************************
Is number smooth?

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool ftbaseissmooth(ae_int_t n, ae_state *) //_state)
{
    ae_int_t i;
    ae_bool result;


    for(i=2; i<=ftbase_ftbasemaxsmoothfactor; i++)
    {
        while(n%i==0)
        {
            n = n/i;
        }
    }
    result = n==1;
    return result;
}


/*************************************************************************
Returns smallest smooth (divisible only by 2, 3, 5) number that is greater
than or equal to max(N,2)

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
ae_int_t ftbasefindsmooth(ae_int_t n, ae_state *_state)
{
    ae_int_t best;
    ae_int_t result;


    best = 2;
    while(best<n)
    {
        best = 2*best;
    }
    ftbase_ftbasefindsmoothrec(n, 1, 2, &best, _state);
    result = best;
    return result;
}


/*************************************************************************
Returns  smallest  smooth  (divisible only by 2, 3, 5) even number that is
greater than or equal to max(N,2)

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
ae_int_t ftbasefindsmootheven(ae_int_t n, ae_state *_state)
{
    ae_int_t best;
    ae_int_t result;


    best = 2;
    while(best<n)
    {
        best = 2*best;
    }
    ftbase_ftbasefindsmoothrec(n, 2, 2, &best, _state);
    result = best;
    return result;
}


/*************************************************************************
Returns estimate of FLOP count for the FFT.

It is only an estimate based on operations count for the PERFECT FFT
and relative inefficiency of the algorithm actually used.

N should be power of 2, estimates are badly wrong for non-power-of-2 N's.

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
double ftbasegetflopestimate(ae_int_t n, ae_state *_state)
{
    double result;


    result = ftbase_ftbaseinefficiencyfactor*(4*n*ae_log(n, _state)/ae_log(2, _state)-6*n+8);
    return result;
}


/*************************************************************************
Recurrent subroutine for the FFTGeneratePlan:

PARAMETERS:
    N                   plan size
    IsReal              whether input is real or not.
                        subroutine MUST NOT ignore this flag because real
                        inputs comes with non-initialized imaginary parts,
                        so ignoring this flag will result in corrupted output
    HalfOut             whether full output or only half of it from 0 to
                        floor(N/2) is needed. This flag may be ignored if
                        doing so will simplify calculations
    Plan                plan array
    PlanSize            size of used part (in integers)
    PrecomputedSize     size of precomputed array allocated yet
    PlanArraySize       plan array size (actual)
    TmpMemSize          temporary memory required size
    BluesteinMemSize    temporary memory required size

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ftbasegenerateplanrec(ae_int_t n,
     ae_int_t tasktype,
     ftplan* plan,
     ae_int_t* plansize,
     ae_int_t* precomputedsize,
     ae_int_t* planarraysize,
     ae_int_t* tmpmemsize,
     ae_int_t* stackmemsize,
     ae_int_t stackptr,
     ae_state *_state)
{
    ae_int_t k;
    ae_int_t m;
    ae_int_t n1;
    ae_int_t n2;
    ae_int_t esize;
    ae_int_t entryoffset;


    
    /*
     * prepare
     */
    if( *plansize+ftbase_ftbaseplanentrysize>(*planarraysize) )
    {
        ftbase_fftarrayresize(&plan->plan, planarraysize, 8*(*planarraysize), _state);
    }
    entryoffset = *plansize;
    esize = ftbase_ftbaseplanentrysize;
    *plansize = *plansize+esize;
    
    /*
     * if N=1, generate empty plan and exit
     */
    if( n==1 )
    {
        plan->plan.ptr.p_int[entryoffset+0] = esize;
        plan->plan.ptr.p_int[entryoffset+1] = -1;
        plan->plan.ptr.p_int[entryoffset+2] = -1;
        plan->plan.ptr.p_int[entryoffset+3] = ftbase_fftemptyplan;
        plan->plan.ptr.p_int[entryoffset+4] = -1;
        plan->plan.ptr.p_int[entryoffset+5] = -1;
        plan->plan.ptr.p_int[entryoffset+6] = -1;
        plan->plan.ptr.p_int[entryoffset+7] = -1;
        return;
    }
    
    /*
     * generate plans
     */
    ftbasefactorize(n, tasktype, &n1, &n2, _state);
    if( tasktype==ftbase_ftbasecffttask||tasktype==ftbase_ftbaserffttask )
    {
        
        /*
         * complex FFT plans
         */
        if( n1!=1 )
        {
            
            /*
             * Cooley-Tukey plan (real or complex)
             *
             * Note that child plans are COMPLEX
             * (whether plan itself is complex or not).
             */
            *tmpmemsize = ae_maxint(*tmpmemsize, 2*n1*n2, _state);
            plan->plan.ptr.p_int[entryoffset+0] = esize;
            plan->plan.ptr.p_int[entryoffset+1] = n1;
            plan->plan.ptr.p_int[entryoffset+2] = n2;
            if( tasktype==ftbase_ftbasecffttask )
            {
                plan->plan.ptr.p_int[entryoffset+3] = ftbase_fftcooleytukeyplan;
            }
            else
            {
                plan->plan.ptr.p_int[entryoffset+3] = ftbase_fftrealcooleytukeyplan;
            }
            plan->plan.ptr.p_int[entryoffset+4] = 0;
            plan->plan.ptr.p_int[entryoffset+5] = *plansize;
            ftbase_ftbasegenerateplanrec(n1, ftbase_ftbasecffttask, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr, _state);
            plan->plan.ptr.p_int[entryoffset+6] = *plansize;
            ftbase_ftbasegenerateplanrec(n2, ftbase_ftbasecffttask, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr, _state);
            plan->plan.ptr.p_int[entryoffset+7] = -1;
            return;
        }
        else
        {
            if( ((n==2||n==3)||n==4)||n==5 )
            {
                
                /*
                 * hard-coded plan
                 */
                plan->plan.ptr.p_int[entryoffset+0] = esize;
                plan->plan.ptr.p_int[entryoffset+1] = n1;
                plan->plan.ptr.p_int[entryoffset+2] = n2;
                plan->plan.ptr.p_int[entryoffset+3] = ftbase_fftcodeletplan;
                plan->plan.ptr.p_int[entryoffset+4] = 0;
                plan->plan.ptr.p_int[entryoffset+5] = -1;
                plan->plan.ptr.p_int[entryoffset+6] = -1;
                plan->plan.ptr.p_int[entryoffset+7] = *precomputedsize;
                if( n==3 )
                {
                    *precomputedsize = *precomputedsize+2;
                }
                if( n==5 )
                {
                    *precomputedsize = *precomputedsize+5;
                }
                return;
            }
            else
            {
                
                /*
                 * Bluestein's plan
                 *
                 * Select such M that M>=2*N-1, M is composite, and M's
                 * factors are 2, 3, 5
                 */
                k = 2*n2-1;
                m = ftbasefindsmooth(k, _state);
                *tmpmemsize = ae_maxint(*tmpmemsize, 2*m, _state);
                plan->plan.ptr.p_int[entryoffset+0] = esize;
                plan->plan.ptr.p_int[entryoffset+1] = n2;
                plan->plan.ptr.p_int[entryoffset+2] = -1;
                plan->plan.ptr.p_int[entryoffset+3] = ftbase_fftbluesteinplan;
                plan->plan.ptr.p_int[entryoffset+4] = m;
                plan->plan.ptr.p_int[entryoffset+5] = *plansize;
                stackptr = stackptr+2*2*m;
                *stackmemsize = ae_maxint(*stackmemsize, stackptr, _state);
                ftbase_ftbasegenerateplanrec(m, ftbase_ftbasecffttask, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr, _state);
                stackptr = stackptr-2*2*m;
                plan->plan.ptr.p_int[entryoffset+6] = -1;
                plan->plan.ptr.p_int[entryoffset+7] = *precomputedsize;
                *precomputedsize = *precomputedsize+2*m+2*n;
                return;
            }
        }
    }
    if( tasktype==ftbase_ftbaserfhttask )
    {
        
        /*
         * real FHT plans
         */
        if( n1!=1 )
        {
            
            /*
             * Cooley-Tukey plan
             *
             */
            *tmpmemsize = ae_maxint(*tmpmemsize, 2*n1*n2, _state);
            plan->plan.ptr.p_int[entryoffset+0] = esize;
            plan->plan.ptr.p_int[entryoffset+1] = n1;
            plan->plan.ptr.p_int[entryoffset+2] = n2;
            plan->plan.ptr.p_int[entryoffset+3] = ftbase_fhtcooleytukeyplan;
            plan->plan.ptr.p_int[entryoffset+4] = 0;
            plan->plan.ptr.p_int[entryoffset+5] = *plansize;
            ftbase_ftbasegenerateplanrec(n1, tasktype, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr, _state);
            plan->plan.ptr.p_int[entryoffset+6] = *plansize;
            ftbase_ftbasegenerateplanrec(n2, tasktype, plan, plansize, precomputedsize, planarraysize, tmpmemsize, stackmemsize, stackptr, _state);
            plan->plan.ptr.p_int[entryoffset+7] = -1;
            return;
        }
        else
        {
            
            /*
             * N2 plan
             */
            plan->plan.ptr.p_int[entryoffset+0] = esize;
            plan->plan.ptr.p_int[entryoffset+1] = n1;
            plan->plan.ptr.p_int[entryoffset+2] = n2;
            plan->plan.ptr.p_int[entryoffset+3] = ftbase_fhtn2plan;
            plan->plan.ptr.p_int[entryoffset+4] = 0;
            plan->plan.ptr.p_int[entryoffset+5] = -1;
            plan->plan.ptr.p_int[entryoffset+6] = -1;
            plan->plan.ptr.p_int[entryoffset+7] = -1;
            if( ((n==2||n==3)||n==4)||n==5 )
            {
                
                /*
                 * hard-coded plan
                 */
                plan->plan.ptr.p_int[entryoffset+0] = esize;
                plan->plan.ptr.p_int[entryoffset+1] = n1;
                plan->plan.ptr.p_int[entryoffset+2] = n2;
                plan->plan.ptr.p_int[entryoffset+3] = ftbase_fhtcodeletplan;
                plan->plan.ptr.p_int[entryoffset+4] = 0;
                plan->plan.ptr.p_int[entryoffset+5] = -1;
                plan->plan.ptr.p_int[entryoffset+6] = -1;
                plan->plan.ptr.p_int[entryoffset+7] = *precomputedsize;
                if( n==3 )
                {
                    *precomputedsize = *precomputedsize+2;
                }
                if( n==5 )
                {
                    *precomputedsize = *precomputedsize+5;
                }
                return;
            }
            return;
        }
    }
}


/*************************************************************************
Recurrent subroutine for precomputing FFT plans

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ftbaseprecomputeplanrec(ftplan* plan,
     ae_int_t entryoffset,
     ae_int_t stackptr,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t n1;
    ae_int_t n2;
    ae_int_t n;
    ae_int_t m;
    ae_int_t offs;
    double v;
    ae_vector emptyarray;
    double bx;
    double by;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&emptyarray, 0, DT_REAL, _state, ae_true);

    if( (plan->plan.ptr.p_int[entryoffset+3]==ftbase_fftcooleytukeyplan||plan->plan.ptr.p_int[entryoffset+3]==ftbase_fftrealcooleytukeyplan)||plan->plan.ptr.p_int[entryoffset+3]==ftbase_fhtcooleytukeyplan )
    {
        ftbase_ftbaseprecomputeplanrec(plan, plan->plan.ptr.p_int[entryoffset+5], stackptr, _state);
        ftbase_ftbaseprecomputeplanrec(plan, plan->plan.ptr.p_int[entryoffset+6], stackptr, _state);
        ae_frame_leave(_state);
        return;
    }
    if( plan->plan.ptr.p_int[entryoffset+3]==ftbase_fftcodeletplan||plan->plan.ptr.p_int[entryoffset+3]==ftbase_fhtcodeletplan )
    {
        n1 = plan->plan.ptr.p_int[entryoffset+1];
        n2 = plan->plan.ptr.p_int[entryoffset+2];
        n = n1*n2;
        if( n==3 )
        {
            offs = plan->plan.ptr.p_int[entryoffset+7];
            plan->precomputed.ptr.p_double[offs+0] = ae_cos(2*ae_pi/3, _state)-1;
            plan->precomputed.ptr.p_double[offs+1] = ae_sin(2*ae_pi/3, _state);
            ae_frame_leave(_state);
            return;
        }
        if( n==5 )
        {
            offs = plan->plan.ptr.p_int[entryoffset+7];
            v = 2*ae_pi/5;
            plan->precomputed.ptr.p_double[offs+0] = (ae_cos(v, _state)+ae_cos(2*v, _state))/2-1;
            plan->precomputed.ptr.p_double[offs+1] = (ae_cos(v, _state)-ae_cos(2*v, _state))/2;
            plan->precomputed.ptr.p_double[offs+2] = -ae_sin(v, _state);
            plan->precomputed.ptr.p_double[offs+3] = -(ae_sin(v, _state)+ae_sin(2*v, _state));
            plan->precomputed.ptr.p_double[offs+4] = ae_sin(v, _state)-ae_sin(2*v, _state);
            ae_frame_leave(_state);
            return;
        }
    }
    if( plan->plan.ptr.p_int[entryoffset+3]==ftbase_fftbluesteinplan )
    {
        ftbase_ftbaseprecomputeplanrec(plan, plan->plan.ptr.p_int[entryoffset+5], stackptr, _state);
        n = plan->plan.ptr.p_int[entryoffset+1];
        m = plan->plan.ptr.p_int[entryoffset+4];
        offs = plan->plan.ptr.p_int[entryoffset+7];
        for(i=0; i<=2*m-1; i++)
        {
            plan->precomputed.ptr.p_double[offs+i] = 0;
        }
        for(i=0; i<=n-1; i++)
        {
            bx = ae_cos(ae_pi*ae_sqr(i, _state)/n, _state);
            by = ae_sin(ae_pi*ae_sqr(i, _state)/n, _state);
            plan->precomputed.ptr.p_double[offs+2*i+0] = bx;
            plan->precomputed.ptr.p_double[offs+2*i+1] = by;
            plan->precomputed.ptr.p_double[offs+2*m+2*i+0] = bx;
            plan->precomputed.ptr.p_double[offs+2*m+2*i+1] = by;
            if( i>0 )
            {
                plan->precomputed.ptr.p_double[offs+2*(m-i)+0] = bx;
                plan->precomputed.ptr.p_double[offs+2*(m-i)+1] = by;
            }
        }
        ftbaseexecuteplanrec(&plan->precomputed, offs, plan, plan->plan.ptr.p_int[entryoffset+5], stackptr, _state);
        ae_frame_leave(_state);
        return;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Twiddle factors calculation

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ffttwcalc(/* Real    */ ae_vector* a,
     ae_int_t aoffset,
     ae_int_t n1,
     ae_int_t n2,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t n;
    ae_int_t idx;
    ae_int_t offs;
    double x;
    double y;
    double twxm1;
    double twy;
    double twbasexm1;
    double twbasey;
    double twrowxm1;
    double twrowy;
    double tmpx;
    double tmpy;
    double v;


    n = n1*n2;
    v = -2*ae_pi/n;
    twbasexm1 = -2*ae_sqr(ae_sin(0.5*v, _state), _state);
    twbasey = ae_sin(v, _state);
    twrowxm1 = 0;
    twrowy = 0;
    for(i=0; i<=n2-1; i++)
    {
        twxm1 = 0;
        twy = 0;
        for(j=0; j<=n1-1; j++)
        {
            idx = i*n1+j;
            offs = aoffset+2*idx;
            x = a->ptr.p_double[offs+0];
            y = a->ptr.p_double[offs+1];
            tmpx = x*twxm1-y*twy;
            tmpy = x*twy+y*twxm1;
            a->ptr.p_double[offs+0] = x+tmpx;
            a->ptr.p_double[offs+1] = y+tmpy;
            
            /*
             * update Tw: Tw(new) = Tw(old)*TwRow
             */
            if( j<n1-1 )
            {
                if( j%ftbase_ftbaseupdatetw==0 )
                {
                    v = -2*ae_pi*i*(j+1)/n;
                    twxm1 = -2*ae_sqr(ae_sin(0.5*v, _state), _state);
                    twy = ae_sin(v, _state);
                }
                else
                {
                    tmpx = twrowxm1+twxm1*twrowxm1-twy*twrowy;
                    tmpy = twrowy+twxm1*twrowy+twy*twrowxm1;
                    twxm1 = twxm1+tmpx;
                    twy = twy+tmpy;
                }
            }
        }
        
        /*
         * update TwRow: TwRow(new) = TwRow(old)*TwBase
         */
        if( i<n2-1 )
        {
            if( j%ftbase_ftbaseupdatetw==0 )
            {
                v = -2*ae_pi*(i+1)/n;
                twrowxm1 = -2*ae_sqr(ae_sin(0.5*v, _state), _state);
                twrowy = ae_sin(v, _state);
            }
            else
            {
                tmpx = twbasexm1+twrowxm1*twbasexm1-twrowy*twbasey;
                tmpy = twbasey+twrowxm1*twbasey+twrowy*twbasexm1;
                twrowxm1 = twrowxm1+tmpx;
                twrowy = twrowy+tmpy;
            }
        }
    }
}


/*************************************************************************
Linear transpose: transpose complex matrix stored in 1-dimensional array

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_internalcomplexlintranspose(/* Real    */ ae_vector* a,
     ae_int_t m,
     ae_int_t n,
     ae_int_t astart,
     /* Real    */ ae_vector* buf,
     ae_state *_state)
{


    ftbase_ffticltrec(a, astart, n, buf, 0, m, m, n, _state);
    ae_v_move(&a->ptr.p_double[astart], 1, &buf->ptr.p_double[0], 1, ae_v_len(astart,astart+2*m*n-1));
}


/*************************************************************************
Linear transpose: transpose real matrix stored in 1-dimensional array

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_internalreallintranspose(/* Real    */ ae_vector* a,
     ae_int_t m,
     ae_int_t n,
     ae_int_t astart,
     /* Real    */ ae_vector* buf,
     ae_state *_state)
{


    ftbase_fftirltrec(a, astart, n, buf, 0, m, m, n, _state);
    ae_v_move(&a->ptr.p_double[astart], 1, &buf->ptr.p_double[0], 1, ae_v_len(astart,astart+m*n-1));
}


/*************************************************************************
Recurrent subroutine for a InternalComplexLinTranspose

Write A^T to B, where:
* A is m*n complex matrix stored in array A as pairs of real/image values,
  beginning from AStart position, with AStride stride
* B is n*m complex matrix stored in array B as pairs of real/image values,
  beginning from BStart position, with BStride stride
stride is measured in complex numbers, i.e. in real/image pairs.

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ffticltrec(/* Real    */ ae_vector* a,
     ae_int_t astart,
     ae_int_t astride,
     /* Real    */ ae_vector* b,
     ae_int_t bstart,
     ae_int_t bstride,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t idx1;
    ae_int_t idx2;
    ae_int_t m2;
    ae_int_t m1;
    ae_int_t n1;


    if( m==0||n==0 )
    {
        return;
    }
    if( ae_maxint(m, n, _state)<=8 )
    {
        m2 = 2*bstride;
        for(i=0; i<=m-1; i++)
        {
            idx1 = bstart+2*i;
            idx2 = astart+2*i*astride;
            for(j=0; j<=n-1; j++)
            {
                b->ptr.p_double[idx1+0] = a->ptr.p_double[idx2+0];
                b->ptr.p_double[idx1+1] = a->ptr.p_double[idx2+1];
                idx1 = idx1+m2;
                idx2 = idx2+2;
            }
        }
        return;
    }
    if( n>m )
    {
        
        /*
         * New partition:
         *
         * "A^T -> B" becomes "(A1 A2)^T -> ( B1 )
         *                                  ( B2 )
         */
        n1 = n/2;
        if( n-n1>=8&&n1%8!=0 )
        {
            n1 = n1+(8-n1%8);
        }
        ae_assert(n-n1>0, "Assertion failed", _state);
        ftbase_ffticltrec(a, astart, astride, b, bstart, bstride, m, n1, _state);
        ftbase_ffticltrec(a, astart+2*n1, astride, b, bstart+2*n1*bstride, bstride, m, n-n1, _state);
    }
    else
    {
        
        /*
         * New partition:
         *
         * "A^T -> B" becomes "( A1 )^T -> ( B1 B2 )
         *                     ( A2 )
         */
        m1 = m/2;
        if( m-m1>=8&&m1%8!=0 )
        {
            m1 = m1+(8-m1%8);
        }
        ae_assert(m-m1>0, "Assertion failed", _state);
        ftbase_ffticltrec(a, astart, astride, b, bstart, bstride, m1, n, _state);
        ftbase_ffticltrec(a, astart+2*m1*astride, astride, b, bstart+2*m1, bstride, m-m1, n, _state);
    }
}


/*************************************************************************
Recurrent subroutine for a InternalRealLinTranspose


  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_fftirltrec(/* Real    */ ae_vector* a,
     ae_int_t astart,
     ae_int_t astride,
     /* Real    */ ae_vector* b,
     ae_int_t bstart,
     ae_int_t bstride,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t idx1;
    ae_int_t idx2;
    ae_int_t m1;
    ae_int_t n1;


    if( m==0||n==0 )
    {
        return;
    }
    if( ae_maxint(m, n, _state)<=8 )
    {
        for(i=0; i<=m-1; i++)
        {
            idx1 = bstart+i;
            idx2 = astart+i*astride;
            for(j=0; j<=n-1; j++)
            {
                b->ptr.p_double[idx1] = a->ptr.p_double[idx2];
                idx1 = idx1+bstride;
                idx2 = idx2+1;
            }
        }
        return;
    }
    if( n>m )
    {
        
        /*
         * New partition:
         *
         * "A^T -> B" becomes "(A1 A2)^T -> ( B1 )
         *                                  ( B2 )
         */
        n1 = n/2;
        if( n-n1>=8&&n1%8!=0 )
        {
            n1 = n1+(8-n1%8);
        }
        ae_assert(n-n1>0, "Assertion failed", _state);
        ftbase_fftirltrec(a, astart, astride, b, bstart, bstride, m, n1, _state);
        ftbase_fftirltrec(a, astart+n1, astride, b, bstart+n1*bstride, bstride, m, n-n1, _state);
    }
    else
    {
        
        /*
         * New partition:
         *
         * "A^T -> B" becomes "( A1 )^T -> ( B1 B2 )
         *                     ( A2 )
         */
        m1 = m/2;
        if( m-m1>=8&&m1%8!=0 )
        {
            m1 = m1+(8-m1%8);
        }
        ae_assert(m-m1>0, "Assertion failed", _state);
        ftbase_fftirltrec(a, astart, astride, b, bstart, bstride, m1, n, _state);
        ftbase_fftirltrec(a, astart+m1*astride, astride, b, bstart+m1, bstride, m-m1, n, _state);
    }
}


/*************************************************************************
recurrent subroutine for FFTFindSmoothRec

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_ftbasefindsmoothrec(ae_int_t n,
     ae_int_t seed,
     ae_int_t leastfactor,
     ae_int_t* best,
     ae_state *_state)
{


    ae_assert(ftbase_ftbasemaxsmoothfactor<=5, "FTBaseFindSmoothRec: internal error!", _state);
    if( seed>=n )
    {
        *best = ae_minint(*best, seed, _state);
        return;
    }
    if( leastfactor<=2 )
    {
        ftbase_ftbasefindsmoothrec(n, seed*2, 2, best, _state);
    }
    if( leastfactor<=3 )
    {
        ftbase_ftbasefindsmoothrec(n, seed*3, 3, best, _state);
    }
    if( leastfactor<=5 )
    {
        ftbase_ftbasefindsmoothrec(n, seed*5, 5, best, _state);
    }
}


/*************************************************************************
Internal subroutine: array resize

  -- ALGLIB --
     Copyright 01.05.2009 by Bochkanov Sergey
*************************************************************************/
static void ftbase_fftarrayresize(/* Integer */ ae_vector* a,
     ae_int_t* asize,
     ae_int_t newasize,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector tmp;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&tmp, 0, DT_INT, _state, ae_true);

    ae_vector_set_length(&tmp, *asize, _state);
    for(i=0; i<=*asize-1; i++)
    {
        tmp.ptr.p_int[i] = a->ptr.p_int[i];
    }
    ae_vector_set_length(a, newasize, _state);
    for(i=0; i<=*asize-1; i++)
    {
        a->ptr.p_int[i] = tmp.ptr.p_int[i];
    }
    *asize = newasize;
    ae_frame_leave(_state);
}


/*************************************************************************
Reference FHT stub
*************************************************************************/
static void ftbase_reffht(/* Real    */ ae_vector* a,
     ae_int_t n,
     ae_int_t offs,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector buf;
    ae_int_t i;
    ae_int_t j;
    double v;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&buf, 0, DT_REAL, _state, ae_true);

    ae_assert(n>0, "RefFHTR1D: incorrect N!", _state);
    ae_vector_set_length(&buf, n, _state);
    for(i=0; i<=n-1; i++)
    {
        v = 0;
        for(j=0; j<=n-1; j++)
        {
            v = v+a->ptr.p_double[offs+j]*(ae_cos(2*ae_pi*i*j/n, _state)+ae_sin(2*ae_pi*i*j/n, _state));
        }
        buf.ptr.p_double[i] = v;
    }
    for(i=0; i<=n-1; i++)
    {
        a->ptr.p_double[offs+i] = buf.ptr.p_double[i];
    }
    ae_frame_leave(_state);
}


ae_bool _ftplan_init(ftplan* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->plan, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->precomputed, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tmpbuf, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->stackbuf, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _ftplan_init_copy(ftplan* dst, ftplan* src, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init_copy(&dst->plan, &src->plan, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->precomputed, &src->precomputed, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tmpbuf, &src->tmpbuf, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->stackbuf, &src->stackbuf, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _ftplan_clear(ftplan* p)
{
    ae_vector_clear(&p->plan);
    ae_vector_clear(&p->precomputed);
    ae_vector_clear(&p->tmpbuf);
    ae_vector_clear(&p->stackbuf);
}




double nulog1p(double x, ae_state *_state)
{
    double z;
    double lp;
    double lq;
    double result;


    z = 1.0+x;
    if( ae_fp_less(z,0.70710678118654752440)||ae_fp_greater(z,1.41421356237309504880) )
    {
        result = ae_log(z, _state);
        return result;
    }
    z = x*x;
    lp = 4.5270000862445199635215E-5;
    lp = lp*x+4.9854102823193375972212E-1;
    lp = lp*x+6.5787325942061044846969E0;
    lp = lp*x+2.9911919328553073277375E1;
    lp = lp*x+6.0949667980987787057556E1;
    lp = lp*x+5.7112963590585538103336E1;
    lp = lp*x+2.0039553499201281259648E1;
    lq = 1.0000000000000000000000E0;
    lq = lq*x+1.5062909083469192043167E1;
    lq = lq*x+8.3047565967967209469434E1;
    lq = lq*x+2.2176239823732856465394E2;
    lq = lq*x+3.0909872225312059774938E2;
    lq = lq*x+2.1642788614495947685003E2;
    lq = lq*x+6.0118660497603843919306E1;
    z = -0.5*z+x*(z*lp/lq);
    result = x+z;
    return result;
}


double nuexpm1(double x, ae_state *_state)
{
    double r;
    double xx;
    double ep;
    double eq;
    double result;


    if( ae_fp_less(x,-0.5)||ae_fp_greater(x,0.5) )
    {
        result = ae_exp(x, _state)-1.0;
        return result;
    }
    xx = x*x;
    ep = 1.2617719307481059087798E-4;
    ep = ep*xx+3.0299440770744196129956E-2;
    ep = ep*xx+9.9999999999999999991025E-1;
    eq = 3.0019850513866445504159E-6;
    eq = eq*xx+2.5244834034968410419224E-3;
    eq = eq*xx+2.2726554820815502876593E-1;
    eq = eq*xx+2.0000000000000000000897E0;
    r = x*ep;
    r = r/(eq-r);
    result = r+r;
    return result;
}


double nucosm1(double x, ae_state *_state)
{
    double xx;
    double c;
    double result;


    if( ae_fp_less(x,-0.25*ae_pi)||ae_fp_greater(x,0.25*ae_pi) )
    {
        result = ae_cos(x, _state)-1;
        return result;
    }
    xx = x*x;
    c = 4.7377507964246204691685E-14;
    c = c*xx-1.1470284843425359765671E-11;
    c = c*xx+2.0876754287081521758361E-9;
    c = c*xx-2.7557319214999787979814E-7;
    c = c*xx+2.4801587301570552304991E-5;
    c = c*xx-1.3888888888888872993737E-3;
    c = c*xx+4.1666666666666666609054E-2;
    result = -0.5*xx+xx*xx*c;
    return result;
}





}

