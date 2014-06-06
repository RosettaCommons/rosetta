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
#include "optimization.h"

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




/*************************************************************************
This object stores state of the nonlinear CG optimizer.

You should use ALGLIB functions to work with this object.
*************************************************************************/
_mincgstate_owner::_mincgstate_owner()
{
    p_struct = (alglib_impl::mincgstate*)alglib_impl::ae_malloc(sizeof(alglib_impl::mincgstate), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_mincgstate_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_mincgstate_owner::_mincgstate_owner(const _mincgstate_owner &rhs)
{
    p_struct = (alglib_impl::mincgstate*)alglib_impl::ae_malloc(sizeof(alglib_impl::mincgstate), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_mincgstate_init_copy(p_struct, const_cast<alglib_impl::mincgstate*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_mincgstate_owner& _mincgstate_owner::operator=(const _mincgstate_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_mincgstate_clear(p_struct);
    if( !alglib_impl::_mincgstate_init_copy(p_struct, const_cast<alglib_impl::mincgstate*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_mincgstate_owner::~_mincgstate_owner()
{
    alglib_impl::_mincgstate_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::mincgstate* _mincgstate_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::mincgstate* _mincgstate_owner::c_ptr() const
{
    return const_cast<alglib_impl::mincgstate*>(p_struct);
}
mincgstate::mincgstate() : _mincgstate_owner() ,needf(p_struct->needf),needfg(p_struct->needfg),xupdated(p_struct->xupdated),f(p_struct->f),g(&p_struct->g),x(&p_struct->x)
{
}

mincgstate::mincgstate(const mincgstate &rhs):_mincgstate_owner(rhs) ,needf(p_struct->needf),needfg(p_struct->needfg),xupdated(p_struct->xupdated),f(p_struct->f),g(&p_struct->g),x(&p_struct->x)
{
}

mincgstate& mincgstate::operator=(const mincgstate &rhs)
{
    if( this==&rhs )
        return *this;
    _mincgstate_owner::operator=(rhs);
    return *this;
}

mincgstate::~mincgstate()
{
}


/*************************************************************************

*************************************************************************/
_mincgreport_owner::_mincgreport_owner()
{
    p_struct = (alglib_impl::mincgreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::mincgreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_mincgreport_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_mincgreport_owner::_mincgreport_owner(const _mincgreport_owner &rhs)
{
    p_struct = (alglib_impl::mincgreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::mincgreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_mincgreport_init_copy(p_struct, const_cast<alglib_impl::mincgreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_mincgreport_owner& _mincgreport_owner::operator=(const _mincgreport_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_mincgreport_clear(p_struct);
    if( !alglib_impl::_mincgreport_init_copy(p_struct, const_cast<alglib_impl::mincgreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_mincgreport_owner::~_mincgreport_owner()
{
    alglib_impl::_mincgreport_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::mincgreport* _mincgreport_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::mincgreport* _mincgreport_owner::c_ptr() const
{
    return const_cast<alglib_impl::mincgreport*>(p_struct);
}
mincgreport::mincgreport() : _mincgreport_owner() ,iterationscount(p_struct->iterationscount),nfev(p_struct->nfev),varidx(p_struct->varidx),terminationtype(p_struct->terminationtype)
{
}

mincgreport::mincgreport(const mincgreport &rhs):_mincgreport_owner(rhs) ,iterationscount(p_struct->iterationscount),nfev(p_struct->nfev),varidx(p_struct->varidx),terminationtype(p_struct->terminationtype)
{
}

mincgreport& mincgreport::operator=(const mincgreport &rhs)
{
    if( this==&rhs )
        return *this;
    _mincgreport_owner::operator=(rhs);
    return *this;
}

mincgreport::~mincgreport()
{
}

/*************************************************************************
        NONLINEAR CONJUGATE GRADIENT METHOD

DESCRIPTION:
The subroutine minimizes function F(x) of N arguments by using one of  the
nonlinear conjugate gradient methods.

These CG methods are globally convergent (even on non-convex functions) as
long as grad(f) is Lipschitz continuous in  a  some  neighborhood  of  the
L = { x : f(x)<=f(x0) }.


REQUIREMENTS:
Algorithm will request following information during its operation:
* function value F and its gradient G (simultaneously) at given point X


USAGE:
1. User initializes algorithm state with MinCGCreate() call
2. User tunes solver parameters with MinCGSetCond(), MinCGSetStpMax() and
   other functions
3. User calls MinCGOptimize() function which takes algorithm  state   and
   pointer (delegate, etc.) to callback function which calculates F/G.
4. User calls MinCGResults() to get solution
5. Optionally, user may call MinCGRestartFrom() to solve another  problem
   with same N but another starting point and/or another function.
   MinCGRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   problem dimension, N>0:
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    X       -   starting point, array[0..N-1].

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 25.03.2010 by Bochkanov Sergey
*************************************************************************/
void mincgcreate(const ae_int_t n, const real_1d_array &x, mincgstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mincgcreate(n, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::mincgstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
        NONLINEAR CONJUGATE GRADIENT METHOD

DESCRIPTION:
The subroutine minimizes function F(x) of N arguments by using one of  the
nonlinear conjugate gradient methods.

These CG methods are globally convergent (even on non-convex functions) as
long as grad(f) is Lipschitz continuous in  a  some  neighborhood  of  the
L = { x : f(x)<=f(x0) }.


REQUIREMENTS:
Algorithm will request following information during its operation:
* function value F and its gradient G (simultaneously) at given point X


USAGE:
1. User initializes algorithm state with MinCGCreate() call
2. User tunes solver parameters with MinCGSetCond(), MinCGSetStpMax() and
   other functions
3. User calls MinCGOptimize() function which takes algorithm  state   and
   pointer (delegate, etc.) to callback function which calculates F/G.
4. User calls MinCGResults() to get solution
5. Optionally, user may call MinCGRestartFrom() to solve another  problem
   with same N but another starting point and/or another function.
   MinCGRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   problem dimension, N>0:
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    X       -   starting point, array[0..N-1].

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 25.03.2010 by Bochkanov Sergey
*************************************************************************/
void mincgcreate(const real_1d_array &x, mincgstate &state)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;

    n = x.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mincgcreate(n, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::mincgstate*>(state.c_ptr()), &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
The subroutine is finite difference variant of MinCGCreate(). It uses
finite differences in order to differentiate target function.

Description below contains information which is specific to this function
only. We recommend to read comments on MinCGCreate() in order to get more
information about creation of CG optimizer.

INPUT PARAMETERS:
    N       -   problem dimension, N>0:
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    X       -   starting point, array[0..N-1].
    DiffStep-   differentiation step, >0

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTES:
1. algorithm uses 4-point central formula for differentiation.
2. differentiation step along I-th axis is equal to DiffStep*S[I] where
   S[] is scaling vector which can be set by MinCGSetScale() call.
3. we recommend you to use moderate values of  differentiation  step.  Too
   large step will result in too large truncation  errors, while too small
   step will result in too large numerical  errors.  1.0E-6  can  be  good
   value to start with.
4. Numerical  differentiation  is   very   inefficient  -   one   gradient
   calculation needs 4*N function evaluations. This function will work for
   any N - either small (1...10), moderate (10...100) or  large  (100...).
   However, performance penalty will be too severe for any N's except  for
   small ones.
   We should also say that code which relies on numerical  differentiation
   is  less  robust  and  precise.  L-BFGS  needs  exact  gradient values.
   Imprecise  gradient may slow down  convergence,  especially  on  highly
   nonlinear problems.
   Thus  we  recommend to use this function for fast prototyping on small-
   dimensional problems only, and to implement analytical gradient as soon
   as possible.

  -- ALGLIB --
     Copyright 16.05.2011 by Bochkanov Sergey
*************************************************************************/
void mincgcreatef(const ae_int_t n, const real_1d_array &x, const double diffstep, mincgstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mincgcreatef(n, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), diffstep, const_cast<alglib_impl::mincgstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
The subroutine is finite difference variant of MinCGCreate(). It uses
finite differences in order to differentiate target function.

Description below contains information which is specific to this function
only. We recommend to read comments on MinCGCreate() in order to get more
information about creation of CG optimizer.

INPUT PARAMETERS:
    N       -   problem dimension, N>0:
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    X       -   starting point, array[0..N-1].
    DiffStep-   differentiation step, >0

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTES:
1. algorithm uses 4-point central formula for differentiation.
2. differentiation step along I-th axis is equal to DiffStep*S[I] where
   S[] is scaling vector which can be set by MinCGSetScale() call.
3. we recommend you to use moderate values of  differentiation  step.  Too
   large step will result in too large truncation  errors, while too small
   step will result in too large numerical  errors.  1.0E-6  can  be  good
   value to start with.
4. Numerical  differentiation  is   very   inefficient  -   one   gradient
   calculation needs 4*N function evaluations. This function will work for
   any N - either small (1...10), moderate (10...100) or  large  (100...).
   However, performance penalty will be too severe for any N's except  for
   small ones.
   We should also say that code which relies on numerical  differentiation
   is  less  robust  and  precise.  L-BFGS  needs  exact  gradient values.
   Imprecise  gradient may slow down  convergence,  especially  on  highly
   nonlinear problems.
   Thus  we  recommend to use this function for fast prototyping on small-
   dimensional problems only, and to implement analytical gradient as soon
   as possible.

  -- ALGLIB --
     Copyright 16.05.2011 by Bochkanov Sergey
*************************************************************************/
void mincgcreatef(const real_1d_array &x, const double diffstep, mincgstate &state)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;

    n = x.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mincgcreatef(n, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), diffstep, const_cast<alglib_impl::mincgstate*>(state.c_ptr()), &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets stopping conditions for CG optimization algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    EpsG    -   >=0
                The  subroutine  finishes  its  work   if   the  condition
                |v|<EpsG is satisfied, where:
                * |.| means Euclidian norm
                * v - scaled gradient vector, v[i]=g[i]*s[i]
                * g - gradient
                * s - scaling coefficients set by MinCGSetScale()
    EpsF    -   >=0
                The  subroutine  finishes  its work if on k+1-th iteration
                the  condition  |F(k+1)-F(k)|<=EpsF*max{|F(k)|,|F(k+1)|,1}
                is satisfied.
    EpsX    -   >=0
                The subroutine finishes its work if  on  k+1-th  iteration
                the condition |v|<=EpsX is fulfilled, where:
                * |.| means Euclidian norm
                * v - scaled step vector, v[i]=dx[i]/s[i]
                * dx - ste pvector, dx=X(k+1)-X(k)
                * s - scaling coefficients set by MinCGSetScale()
    MaxIts  -   maximum number of iterations. If MaxIts=0, the  number  of
                iterations is unlimited.

Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to
automatic stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsetcond(const mincgstate &state, const double epsg, const double epsf, const double epsx, const ae_int_t maxits)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mincgsetcond(const_cast<alglib_impl::mincgstate*>(state.c_ptr()), epsg, epsf, epsx, maxits, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets scaling coefficients for CG optimizer.

ALGLIB optimizers use scaling matrices to test stopping  conditions  (step
size and gradient are scaled before comparison with tolerances).  Scale of
the I-th variable is a translation invariant measure of:
a) "how large" the variable is
b) how large the step should be to make significant changes in the function

Scaling is also used by finite difference variant of CG optimizer  -  step
along I-th axis is equal to DiffStep*S[I].

In   most   optimizers  (and  in  the  CG  too)  scaling is NOT a form  of
preconditioning. It just  affects  stopping  conditions.  You  should  set
preconditioner by separate call to one of the MinCGSetPrec...() functions.

There  is  special  preconditioning  mode, however,  which  uses   scaling
coefficients to form diagonal preconditioning matrix. You  can  turn  this
mode on, if you want.   But  you should understand that scaling is not the
same thing as preconditioning - these are two different, although  related
forms of tuning solver.

INPUT PARAMETERS:
    State   -   structure stores algorithm state
    S       -   array[N], non-zero scaling coefficients
                S[i] may be negative, sign doesn't matter.

  -- ALGLIB --
     Copyright 14.01.2011 by Bochkanov Sergey
*************************************************************************/
void mincgsetscale(const mincgstate &state, const real_1d_array &s)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mincgsetscale(const_cast<alglib_impl::mincgstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(s.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    NeedXRep-   whether iteration reports are needed or not

If NeedXRep is True, algorithm will call rep() callback function if  it is
provided to MinCGOptimize().

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsetxrep(const mincgstate &state, const bool needxrep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mincgsetxrep(const_cast<alglib_impl::mincgstate*>(state.c_ptr()), needxrep, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets CG algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    CGType  -   algorithm type:
                * -1    automatic selection of the best algorithm
                * 0     DY (Dai and Yuan) algorithm
                * 1     Hybrid DY-HS algorithm

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsetcgtype(const mincgstate &state, const ae_int_t cgtype)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mincgsetcgtype(const_cast<alglib_impl::mincgstate*>(state.c_ptr()), cgtype, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets maximum step length

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    StpMax  -   maximum step length, >=0. Set StpMax to 0.0,  if you don't
                want to limit step length.

Use this subroutine when you optimize target function which contains exp()
or  other  fast  growing  functions,  and optimization algorithm makes too
large  steps  which  leads  to overflow. This function allows us to reject
steps  that  are  too  large  (and  therefore  expose  us  to the possible
overflow) without actually calculating function value at the x+stp*d.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsetstpmax(const mincgstate &state, const double stpmax)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mincgsetstpmax(const_cast<alglib_impl::mincgstate*>(state.c_ptr()), stpmax, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function allows to suggest initial step length to the CG algorithm.

Suggested  step  length  is used as starting point for the line search. It
can be useful when you have  badly  scaled  problem,  i.e.  when  ||grad||
(which is used as initial estimate for the first step) is many  orders  of
magnitude different from the desired step.

Line search  may  fail  on  such problems without good estimate of initial
step length. Imagine, for example, problem with ||grad||=10^50 and desired
step equal to 0.1 Line  search function will use 10^50  as  initial  step,
then  it  will  decrease step length by 2 (up to 20 attempts) and will get
10^44, which is still too large.

This function allows us to tell than line search should  be  started  from
some moderate step length, like 1.0, so algorithm will be able  to  detect
desired step length in a several searches.

Default behavior (when no step is suggested) is to use preconditioner,  if
it is available, to generate initial estimate of step length.

This function influences only first iteration of algorithm. It  should  be
called between MinCGCreate/MinCGRestartFrom() call and MinCGOptimize call.
Suggested step is ignored if you have preconditioner.

INPUT PARAMETERS:
    State   -   structure used to store algorithm state.
    Stp     -   initial estimate of the step length.
                Can be zero (no estimate).

  -- ALGLIB --
     Copyright 30.07.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsuggeststep(const mincgstate &state, const double stp)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mincgsuggeststep(const_cast<alglib_impl::mincgstate*>(state.c_ptr()), stp, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Modification of the preconditioner: preconditioning is turned off.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTE:  you  can  change  preconditioner  "on  the  fly",  during algorithm
iterations.

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsetprecdefault(const mincgstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mincgsetprecdefault(const_cast<alglib_impl::mincgstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Modification  of  the  preconditioner:  diagonal of approximate Hessian is
used.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    D       -   diagonal of the approximate Hessian, array[0..N-1],
                (if larger, only leading N elements are used).

NOTE:  you  can  change  preconditioner  "on  the  fly",  during algorithm
iterations.

NOTE 2: D[i] should be positive. Exception will be thrown otherwise.

NOTE 3: you should pass diagonal of approximate Hessian - NOT ITS INVERSE.

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsetprecdiag(const mincgstate &state, const real_1d_array &d)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mincgsetprecdiag(const_cast<alglib_impl::mincgstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(d.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Modification of the preconditioner: scale-based diagonal preconditioning.

This preconditioning mode can be useful when you  don't  have  approximate
diagonal of Hessian, but you know that your  variables  are  badly  scaled
(for  example,  one  variable is in [1,10], and another in [1000,100000]),
and most part of the ill-conditioning comes from different scales of vars.

In this case simple  scale-based  preconditioner,  with H[i] = 1/(s[i]^2),
can greatly improve convergence.

IMPRTANT: you should set scale of your variables with MinCGSetScale() call
(before or after MinCGSetPrecScale() call). Without knowledge of the scale
of your variables scale-based preconditioner will be just unit matrix.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTE:  you  can  change  preconditioner  "on  the  fly",  during algorithm
iterations.

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsetprecscale(const mincgstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mincgsetprecscale(const_cast<alglib_impl::mincgstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function provides reverse communication interface
Reverse communication interface is not documented or recommended to use.
See below for functions which provide better documented API
*************************************************************************/
bool mincgiteration(const mincgstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        ae_bool result = alglib_impl::mincgiteration(const_cast<alglib_impl::mincgstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<bool*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}


void mincgoptimize(mincgstate &state,
    void (*func)(const real_1d_array &x, double &func, void *ptr),
    void  (*rep)(const real_1d_array &x, double func, void *ptr), 
    void *ptr)
{
    alglib_impl::ae_state _alglib_env_state;
    if( func==NULL )
        throw ap_error("ALGLIB: error in 'mincgoptimize()' (func is NULL)");
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        while( alglib_impl::mincgiteration(state.c_ptr(), &_alglib_env_state) )
        {
            if( state.needf )
            {
                func(state.x, state.f, ptr);
                continue;
            }
            if( state.xupdated )
            {
                if( rep!=NULL )
                    rep(state.x, state.f, ptr);
                continue;
            }
            throw ap_error("ALGLIB: error in 'mincgoptimize' (some derivatives were not provided?)");
        }
        alglib_impl::ae_state_clear(&_alglib_env_state);
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}


void mincgoptimize(mincgstate &state,
    void (*grad)(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr),
    void  (*rep)(const real_1d_array &x, double func, void *ptr), 
    void *ptr)
{
    alglib_impl::ae_state _alglib_env_state;
    if( grad==NULL )
        throw ap_error("ALGLIB: error in 'mincgoptimize()' (grad is NULL)");
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        while( alglib_impl::mincgiteration(state.c_ptr(), &_alglib_env_state) )
        {
            if( state.needfg )
            {
                grad(state.x, state.f, state.g, ptr);
                continue;
            }
            if( state.xupdated )
            {
                if( rep!=NULL )
                    rep(state.x, state.f, ptr);
                continue;
            }
            throw ap_error("ALGLIB: error in 'mincgoptimize' (some derivatives were not provided?)");
        }
        alglib_impl::ae_state_clear(&_alglib_env_state);
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}



/*************************************************************************
Conjugate gradient results

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    X       -   array[0..N-1], solution
    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    * -7    gradient verification failed.
                            See MinCGSetGradientCheck() for more information.
                    *  1    relative function improvement is no more than
                            EpsF.
                    *  2    relative step is no more than EpsX.
                    *  4    gradient norm is no more than EpsG
                    *  5    MaxIts steps was taken
                    *  7    stopping conditions are too stringent,
                            further improvement is impossible,
                            we return best X found so far
                    *  8    terminated by user
                * Rep.IterationsCount contains iterations count
                * NFEV countains number of function calculations

  -- ALGLIB --
     Copyright 20.04.2009 by Bochkanov Sergey
*************************************************************************/
void mincgresults(const mincgstate &state, real_1d_array &x, mincgreport &rep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mincgresults(const_cast<alglib_impl::mincgstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::mincgreport*>(rep.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Conjugate gradient results

Buffered implementation of MinCGResults(), which uses pre-allocated buffer
to store X[]. If buffer size is  too  small,  it  resizes  buffer.  It  is
intended to be used in the inner cycles of performance critical algorithms
where array reallocation penalty is too large to be ignored.

  -- ALGLIB --
     Copyright 20.04.2009 by Bochkanov Sergey
*************************************************************************/
void mincgresultsbuf(const mincgstate &state, real_1d_array &x, mincgreport &rep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mincgresultsbuf(const_cast<alglib_impl::mincgstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::mincgreport*>(rep.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This  subroutine  restarts  CG  algorithm from new point. All optimization
parameters are left unchanged.

This  function  allows  to  solve multiple  optimization  problems  (which
must have same number of dimensions) without object reallocation penalty.

INPUT PARAMETERS:
    State   -   structure used to store algorithm state.
    X       -   new starting point.

  -- ALGLIB --
     Copyright 30.07.2010 by Bochkanov Sergey
*************************************************************************/
void mincgrestartfrom(const mincgstate &state, const real_1d_array &x)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mincgrestartfrom(const_cast<alglib_impl::mincgstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************

This  subroutine  turns  on  verification  of  the  user-supplied analytic
gradient:
* user calls this subroutine before optimization begins
* MinCGOptimize() is called
* prior to  actual  optimization, for each component  of  parameters being
  optimized X[i] algorithm performs following steps:
  * two trial steps are made to X[i]-TestStep*S[i] and X[i]+TestStep*S[i],
    where X[i] is i-th component of the initial point and S[i] is a  scale
    of i-th parameter
  * F(X) is evaluated at these trial points
  * we perform one more evaluation in the middle point of the interval
  * we  build  cubic  model using function values and derivatives at trial
    points and we compare its prediction with actual value in  the  middle
    point
  * in case difference between prediction and actual value is higher  than
    some predetermined threshold, algorithm stops with completion code -7;
    Rep.VarIdx is set to index of the parameter with incorrect derivative.
* after verification is over, algorithm proceeds to the actual optimization.

NOTE 1: verification  needs  N (parameters count) gradient evaluations. It
        is very costly and you should use  it  only  for  low  dimensional
        problems,  when  you  want  to  be  sure  that  you've   correctly
        calculated  analytic  derivatives.  You  should  not use it in the
        production code (unless you want to check derivatives provided  by
        some third party).

NOTE 2: you  should  carefully  choose  TestStep. Value which is too large
        (so large that function behaviour is significantly non-cubic) will
        lead to false alarms. You may use  different  step  for  different
        parameters by means of setting scale with MinCGSetScale().

NOTE 3: this function may lead to false positives. In case it reports that
        I-th  derivative was calculated incorrectly, you may decrease test
        step  and  try  one  more  time  - maybe your function changes too
        sharply  and  your  step  is  too  large for such rapidly chanding
        function.

INPUT PARAMETERS:
    State       -   structure used to store algorithm state
    TestStep    -   verification step:
                    * TestStep=0 turns verification off
                    * TestStep>0 activates verification

  -- ALGLIB --
     Copyright 31.05.2012 by Bochkanov Sergey
*************************************************************************/
void mincgsetgradientcheck(const mincgstate &state, const double teststep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mincgsetgradientcheck(const_cast<alglib_impl::mincgstate*>(state.c_ptr()), teststep, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This object stores nonlinear optimizer state.
You should use functions provided by MinBLEIC subpackage to work with this
object
*************************************************************************/
_minbleicstate_owner::_minbleicstate_owner()
{
    p_struct = (alglib_impl::minbleicstate*)alglib_impl::ae_malloc(sizeof(alglib_impl::minbleicstate), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_minbleicstate_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_minbleicstate_owner::_minbleicstate_owner(const _minbleicstate_owner &rhs)
{
    p_struct = (alglib_impl::minbleicstate*)alglib_impl::ae_malloc(sizeof(alglib_impl::minbleicstate), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_minbleicstate_init_copy(p_struct, const_cast<alglib_impl::minbleicstate*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_minbleicstate_owner& _minbleicstate_owner::operator=(const _minbleicstate_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_minbleicstate_clear(p_struct);
    if( !alglib_impl::_minbleicstate_init_copy(p_struct, const_cast<alglib_impl::minbleicstate*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_minbleicstate_owner::~_minbleicstate_owner()
{
    alglib_impl::_minbleicstate_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::minbleicstate* _minbleicstate_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::minbleicstate* _minbleicstate_owner::c_ptr() const
{
    return const_cast<alglib_impl::minbleicstate*>(p_struct);
}
minbleicstate::minbleicstate() : _minbleicstate_owner() ,needf(p_struct->needf),needfg(p_struct->needfg),xupdated(p_struct->xupdated),f(p_struct->f),g(&p_struct->g),x(&p_struct->x)
{
}

minbleicstate::minbleicstate(const minbleicstate &rhs):_minbleicstate_owner(rhs) ,needf(p_struct->needf),needfg(p_struct->needfg),xupdated(p_struct->xupdated),f(p_struct->f),g(&p_struct->g),x(&p_struct->x)
{
}

minbleicstate& minbleicstate::operator=(const minbleicstate &rhs)
{
    if( this==&rhs )
        return *this;
    _minbleicstate_owner::operator=(rhs);
    return *this;
}

minbleicstate::~minbleicstate()
{
}


/*************************************************************************
This structure stores optimization report:
* InnerIterationsCount      number of inner iterations
* OuterIterationsCount      number of outer iterations
* NFEV                      number of gradient evaluations
* TerminationType           termination type (see below)

TERMINATION CODES

TerminationType field contains completion code, which can be:
  -10   unsupported combination of algorithm settings:
        1) StpMax is set to non-zero value,
        AND 2) non-default preconditioner is used.
        You can't use both features at the same moment,
        so you have to choose one of them (and to turn
        off another one).
  -7    gradient verification failed.
        See MinBLEICSetGradientCheck() for more information.
  -3    inconsistent constraints. Feasible point is
        either nonexistent or too hard to find. Try to
        restart optimizer with better initial
        approximation
   4    conditions on constraints are fulfilled
        with error less than or equal to EpsC
   5    MaxIts steps was taken
   7    stopping conditions are too stringent,
        further improvement is impossible,
        X contains best point found so far.

ADDITIONAL FIELDS

There are additional fields which can be used for debugging:
* DebugEqErr                error in the equality constraints (2-norm)
* DebugFS                   f, calculated at projection of initial point
                            to the feasible set
* DebugFF                   f, calculated at the final point
* DebugDX                   |X_start-X_final|
*************************************************************************/
_minbleicreport_owner::_minbleicreport_owner()
{
    p_struct = (alglib_impl::minbleicreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::minbleicreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_minbleicreport_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_minbleicreport_owner::_minbleicreport_owner(const _minbleicreport_owner &rhs)
{
    p_struct = (alglib_impl::minbleicreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::minbleicreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_minbleicreport_init_copy(p_struct, const_cast<alglib_impl::minbleicreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_minbleicreport_owner& _minbleicreport_owner::operator=(const _minbleicreport_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_minbleicreport_clear(p_struct);
    if( !alglib_impl::_minbleicreport_init_copy(p_struct, const_cast<alglib_impl::minbleicreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_minbleicreport_owner::~_minbleicreport_owner()
{
    alglib_impl::_minbleicreport_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::minbleicreport* _minbleicreport_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::minbleicreport* _minbleicreport_owner::c_ptr() const
{
    return const_cast<alglib_impl::minbleicreport*>(p_struct);
}
minbleicreport::minbleicreport() : _minbleicreport_owner() ,inneriterationscount(p_struct->inneriterationscount),outeriterationscount(p_struct->outeriterationscount),nfev(p_struct->nfev),varidx(p_struct->varidx),terminationtype(p_struct->terminationtype),debugeqerr(p_struct->debugeqerr),debugfs(p_struct->debugfs),debugff(p_struct->debugff),debugdx(p_struct->debugdx),debugfeasqpits(p_struct->debugfeasqpits),debugfeasgpaits(p_struct->debugfeasgpaits)
{
}

minbleicreport::minbleicreport(const minbleicreport &rhs):_minbleicreport_owner(rhs) ,inneriterationscount(p_struct->inneriterationscount),outeriterationscount(p_struct->outeriterationscount),nfev(p_struct->nfev),varidx(p_struct->varidx),terminationtype(p_struct->terminationtype),debugeqerr(p_struct->debugeqerr),debugfs(p_struct->debugfs),debugff(p_struct->debugff),debugdx(p_struct->debugdx),debugfeasqpits(p_struct->debugfeasqpits),debugfeasgpaits(p_struct->debugfeasgpaits)
{
}

minbleicreport& minbleicreport::operator=(const minbleicreport &rhs)
{
    if( this==&rhs )
        return *this;
    _minbleicreport_owner::operator=(rhs);
    return *this;
}

minbleicreport::~minbleicreport()
{
}

/*************************************************************************
                     BOUND CONSTRAINED OPTIMIZATION
       WITH ADDITIONAL LINEAR EQUALITY AND INEQUALITY CONSTRAINTS

DESCRIPTION:
The  subroutine  minimizes  function   F(x)  of N arguments subject to any
combination of:
* bound constraints
* linear inequality constraints
* linear equality constraints

REQUIREMENTS:
* user must provide function value and gradient
* starting point X0 must be feasible or
  not too far away from the feasible set
* grad(f) must be Lipschitz continuous on a level set:
  L = { x : f(x)<=f(x0) }
* function must be defined everywhere on the feasible set F

USAGE:

Constrained optimization if far more complex than the unconstrained one.
Here we give very brief outline of the BLEIC optimizer. We strongly recommend
you to read examples in the ALGLIB Reference Manual and to read ALGLIB User Guide
on optimization, which is available at http://www.alglib.net/optimization/

1. User initializes algorithm state with MinBLEICCreate() call

2. USer adds boundary and/or linear constraints by calling
   MinBLEICSetBC() and MinBLEICSetLC() functions.

3. User sets stopping conditions for underlying unconstrained solver
   with MinBLEICSetInnerCond() call.
   This function controls accuracy of underlying optimization algorithm.

4. User sets stopping conditions for outer iteration by calling
   MinBLEICSetOuterCond() function.
   This function controls handling of boundary and inequality constraints.

5. Additionally, user may set limit on number of internal iterations
   by MinBLEICSetMaxIts() call.
   This function allows to prevent algorithm from looping forever.

6. User calls MinBLEICOptimize() function which takes algorithm  state and
   pointer (delegate, etc.) to callback function which calculates F/G.

7. User calls MinBLEICResults() to get solution

8. Optionally user may call MinBLEICRestartFrom() to solve another problem
   with same N but another starting point.
   MinBLEICRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   problem dimension, N>0:
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size ofX
    X       -   starting point, array[N]:
                * it is better to set X to a feasible point
                * but X can be infeasible, in which case algorithm will try
                  to find feasible point first, using X as initial
                  approximation.

OUTPUT PARAMETERS:
    State   -   structure stores algorithm state

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleiccreate(const ae_int_t n, const real_1d_array &x, minbleicstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleiccreate(n, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
                     BOUND CONSTRAINED OPTIMIZATION
       WITH ADDITIONAL LINEAR EQUALITY AND INEQUALITY CONSTRAINTS

DESCRIPTION:
The  subroutine  minimizes  function   F(x)  of N arguments subject to any
combination of:
* bound constraints
* linear inequality constraints
* linear equality constraints

REQUIREMENTS:
* user must provide function value and gradient
* starting point X0 must be feasible or
  not too far away from the feasible set
* grad(f) must be Lipschitz continuous on a level set:
  L = { x : f(x)<=f(x0) }
* function must be defined everywhere on the feasible set F

USAGE:

Constrained optimization if far more complex than the unconstrained one.
Here we give very brief outline of the BLEIC optimizer. We strongly recommend
you to read examples in the ALGLIB Reference Manual and to read ALGLIB User Guide
on optimization, which is available at http://www.alglib.net/optimization/

1. User initializes algorithm state with MinBLEICCreate() call

2. USer adds boundary and/or linear constraints by calling
   MinBLEICSetBC() and MinBLEICSetLC() functions.

3. User sets stopping conditions for underlying unconstrained solver
   with MinBLEICSetInnerCond() call.
   This function controls accuracy of underlying optimization algorithm.

4. User sets stopping conditions for outer iteration by calling
   MinBLEICSetOuterCond() function.
   This function controls handling of boundary and inequality constraints.

5. Additionally, user may set limit on number of internal iterations
   by MinBLEICSetMaxIts() call.
   This function allows to prevent algorithm from looping forever.

6. User calls MinBLEICOptimize() function which takes algorithm  state and
   pointer (delegate, etc.) to callback function which calculates F/G.

7. User calls MinBLEICResults() to get solution

8. Optionally user may call MinBLEICRestartFrom() to solve another problem
   with same N but another starting point.
   MinBLEICRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   problem dimension, N>0:
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size ofX
    X       -   starting point, array[N]:
                * it is better to set X to a feasible point
                * but X can be infeasible, in which case algorithm will try
                  to find feasible point first, using X as initial
                  approximation.

OUTPUT PARAMETERS:
    State   -   structure stores algorithm state

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleiccreate(const real_1d_array &x, minbleicstate &state)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;

    n = x.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleiccreate(n, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
The subroutine is finite difference variant of MinBLEICCreate().  It  uses
finite differences in order to differentiate target function.

Description below contains information which is specific to  this function
only. We recommend to read comments on MinBLEICCreate() in  order  to  get
more information about creation of BLEIC optimizer.

INPUT PARAMETERS:
    N       -   problem dimension, N>0:
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    X       -   starting point, array[0..N-1].
    DiffStep-   differentiation step, >0

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTES:
1. algorithm uses 4-point central formula for differentiation.
2. differentiation step along I-th axis is equal to DiffStep*S[I] where
   S[] is scaling vector which can be set by MinBLEICSetScale() call.
3. we recommend you to use moderate values of  differentiation  step.  Too
   large step will result in too large truncation  errors, while too small
   step will result in too large numerical  errors.  1.0E-6  can  be  good
   value to start with.
4. Numerical  differentiation  is   very   inefficient  -   one   gradient
   calculation needs 4*N function evaluations. This function will work for
   any N - either small (1...10), moderate (10...100) or  large  (100...).
   However, performance penalty will be too severe for any N's except  for
   small ones.
   We should also say that code which relies on numerical  differentiation
   is  less  robust and precise. CG needs exact gradient values. Imprecise
   gradient may slow  down  convergence, especially  on  highly  nonlinear
   problems.
   Thus  we  recommend to use this function for fast prototyping on small-
   dimensional problems only, and to implement analytical gradient as soon
   as possible.

  -- ALGLIB --
     Copyright 16.05.2011 by Bochkanov Sergey
*************************************************************************/
void minbleiccreatef(const ae_int_t n, const real_1d_array &x, const double diffstep, minbleicstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleiccreatef(n, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), diffstep, const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
The subroutine is finite difference variant of MinBLEICCreate().  It  uses
finite differences in order to differentiate target function.

Description below contains information which is specific to  this function
only. We recommend to read comments on MinBLEICCreate() in  order  to  get
more information about creation of BLEIC optimizer.

INPUT PARAMETERS:
    N       -   problem dimension, N>0:
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    X       -   starting point, array[0..N-1].
    DiffStep-   differentiation step, >0

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTES:
1. algorithm uses 4-point central formula for differentiation.
2. differentiation step along I-th axis is equal to DiffStep*S[I] where
   S[] is scaling vector which can be set by MinBLEICSetScale() call.
3. we recommend you to use moderate values of  differentiation  step.  Too
   large step will result in too large truncation  errors, while too small
   step will result in too large numerical  errors.  1.0E-6  can  be  good
   value to start with.
4. Numerical  differentiation  is   very   inefficient  -   one   gradient
   calculation needs 4*N function evaluations. This function will work for
   any N - either small (1...10), moderate (10...100) or  large  (100...).
   However, performance penalty will be too severe for any N's except  for
   small ones.
   We should also say that code which relies on numerical  differentiation
   is  less  robust and precise. CG needs exact gradient values. Imprecise
   gradient may slow  down  convergence, especially  on  highly  nonlinear
   problems.
   Thus  we  recommend to use this function for fast prototyping on small-
   dimensional problems only, and to implement analytical gradient as soon
   as possible.

  -- ALGLIB --
     Copyright 16.05.2011 by Bochkanov Sergey
*************************************************************************/
void minbleiccreatef(const real_1d_array &x, const double diffstep, minbleicstate &state)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;

    n = x.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleiccreatef(n, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), diffstep, const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets boundary constraints for BLEIC optimizer.

Boundary constraints are inactive by default (after initial creation).
They are preserved after algorithm restart with MinBLEICRestartFrom().

INPUT PARAMETERS:
    State   -   structure stores algorithm state
    BndL    -   lower bounds, array[N].
                If some (all) variables are unbounded, you may specify
                very small number or -INF.
    BndU    -   upper bounds, array[N].
                If some (all) variables are unbounded, you may specify
                very large number or +INF.

NOTE 1: it is possible to specify BndL[i]=BndU[i]. In this case I-th
variable will be "frozen" at X[i]=BndL[i]=BndU[i].

NOTE 2: this solver has following useful properties:
* bound constraints are always satisfied exactly
* function is evaluated only INSIDE area specified by  bound  constraints,
  even  when  numerical  differentiation is used (algorithm adjusts  nodes
  according to boundary constraints)

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetbc(const minbleicstate &state, const real_1d_array &bndl, const real_1d_array &bndu)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleicsetbc(const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(bndl.c_ptr()), const_cast<alglib_impl::ae_vector*>(bndu.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets linear constraints for BLEIC optimizer.

Linear constraints are inactive by default (after initial creation).
They are preserved after algorithm restart with MinBLEICRestartFrom().

INPUT PARAMETERS:
    State   -   structure previously allocated with MinBLEICCreate call.
    C       -   linear constraints, array[K,N+1].
                Each row of C represents one constraint, either equality
                or inequality (see below):
                * first N elements correspond to coefficients,
                * last element corresponds to the right part.
                All elements of C (including right part) must be finite.
    CT      -   type of constraints, array[K]:
                * if CT[i]>0, then I-th constraint is C[i,*]*x >= C[i,n+1]
                * if CT[i]=0, then I-th constraint is C[i,*]*x  = C[i,n+1]
                * if CT[i]<0, then I-th constraint is C[i,*]*x <= C[i,n+1]
    K       -   number of equality/inequality constraints, K>=0:
                * if given, only leading K elements of C/CT are used
                * if not given, automatically determined from sizes of C/CT

NOTE 1: linear (non-bound) constraints are satisfied only approximately:
* there always exists some minor violation (about Epsilon in magnitude)
  due to rounding errors
* numerical differentiation, if used, may  lead  to  function  evaluations
  outside  of the feasible  area,   because   algorithm  does  NOT  change
  numerical differentiation formula according to linear constraints.
If you want constraints to be  satisfied  exactly, try to reformulate your
problem  in  such  manner  that  all constraints will become boundary ones
(this kind of constraints is always satisfied exactly, both in  the  final
solution and in all intermediate points).

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetlc(const minbleicstate &state, const real_2d_array &c, const integer_1d_array &ct, const ae_int_t k)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleicsetlc(const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), const_cast<alglib_impl::ae_matrix*>(c.c_ptr()), const_cast<alglib_impl::ae_vector*>(ct.c_ptr()), k, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets linear constraints for BLEIC optimizer.

Linear constraints are inactive by default (after initial creation).
They are preserved after algorithm restart with MinBLEICRestartFrom().

INPUT PARAMETERS:
    State   -   structure previously allocated with MinBLEICCreate call.
    C       -   linear constraints, array[K,N+1].
                Each row of C represents one constraint, either equality
                or inequality (see below):
                * first N elements correspond to coefficients,
                * last element corresponds to the right part.
                All elements of C (including right part) must be finite.
    CT      -   type of constraints, array[K]:
                * if CT[i]>0, then I-th constraint is C[i,*]*x >= C[i,n+1]
                * if CT[i]=0, then I-th constraint is C[i,*]*x  = C[i,n+1]
                * if CT[i]<0, then I-th constraint is C[i,*]*x <= C[i,n+1]
    K       -   number of equality/inequality constraints, K>=0:
                * if given, only leading K elements of C/CT are used
                * if not given, automatically determined from sizes of C/CT

NOTE 1: linear (non-bound) constraints are satisfied only approximately:
* there always exists some minor violation (about Epsilon in magnitude)
  due to rounding errors
* numerical differentiation, if used, may  lead  to  function  evaluations
  outside  of the feasible  area,   because   algorithm  does  NOT  change
  numerical differentiation formula according to linear constraints.
If you want constraints to be  satisfied  exactly, try to reformulate your
problem  in  such  manner  that  all constraints will become boundary ones
(this kind of constraints is always satisfied exactly, both in  the  final
solution and in all intermediate points).

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetlc(const minbleicstate &state, const real_2d_array &c, const integer_1d_array &ct)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t k;
    if( (c.rows()!=ct.length()))
        throw ap_error("Error while calling 'minbleicsetlc': looks like one of arguments has wrong size");
    k = c.rows();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleicsetlc(const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), const_cast<alglib_impl::ae_matrix*>(c.c_ptr()), const_cast<alglib_impl::ae_vector*>(ct.c_ptr()), k, &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets stopping conditions for the underlying nonlinear CG
optimizer. It controls overall accuracy of solution. These conditions
should be strict enough in order for algorithm to converge.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    EpsG    -   >=0
                The  subroutine  finishes  its  work   if   the  condition
                |v|<EpsG is satisfied, where:
                * |.| means Euclidian norm
                * v - scaled gradient vector, v[i]=g[i]*s[i]
                * g - gradient
                * s - scaling coefficients set by MinBLEICSetScale()
    EpsF    -   >=0
                The  subroutine  finishes  its work if on k+1-th iteration
                the  condition  |F(k+1)-F(k)|<=EpsF*max{|F(k)|,|F(k+1)|,1}
                is satisfied.
    EpsX    -   >=0
                The subroutine finishes its work if  on  k+1-th  iteration
                the condition |v|<=EpsX is fulfilled, where:
                * |.| means Euclidian norm
                * v - scaled step vector, v[i]=dx[i]/s[i]
                * dx - ste pvector, dx=X(k+1)-X(k)
                * s - scaling coefficients set by MinBLEICSetScale()

Passing EpsG=0, EpsF=0 and EpsX=0 (simultaneously) will lead to
automatic stopping criterion selection.

These conditions are used to terminate inner iterations. However, you
need to tune termination conditions for outer iterations too.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetinnercond(const minbleicstate &state, const double epsg, const double epsf, const double epsx)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleicsetinnercond(const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), epsg, epsf, epsx, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets stopping conditions for outer iteration of BLEIC algo.

These conditions control accuracy of constraint handling and amount of
infeasibility allowed in the solution.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    EpsX    -   >0, stopping condition on outer iteration step length
    EpsI    -   >0, stopping condition on infeasibility

Both EpsX and EpsI must be non-zero.

MEANING OF EpsX

EpsX  is  a  stopping  condition for outer iterations. Algorithm will stop
when  solution  of  the  current  modified  subproblem will be within EpsX
(using 2-norm) of the previous solution.

MEANING OF EpsI

EpsI controls feasibility properties -  algorithm  won't  stop  until  all
inequality constraints will be satisfied with error (distance from current
point to the feasible area) at most EpsI.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetoutercond(const minbleicstate &state, const double epsx, const double epsi)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleicsetoutercond(const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), epsx, epsi, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets scaling coefficients for BLEIC optimizer.

ALGLIB optimizers use scaling matrices to test stopping  conditions  (step
size and gradient are scaled before comparison with tolerances).  Scale of
the I-th variable is a translation invariant measure of:
a) "how large" the variable is
b) how large the step should be to make significant changes in the function

Scaling is also used by finite difference variant of the optimizer  - step
along I-th axis is equal to DiffStep*S[I].

In  most  optimizers  (and  in  the  BLEIC  too)  scaling is NOT a form of
preconditioning. It just  affects  stopping  conditions.  You  should  set
preconditioner  by  separate  call  to  one  of  the  MinBLEICSetPrec...()
functions.

There is a special  preconditioning  mode, however,  which  uses   scaling
coefficients to form diagonal preconditioning matrix. You  can  turn  this
mode on, if you want.   But  you should understand that scaling is not the
same thing as preconditioning - these are two different, although  related
forms of tuning solver.

INPUT PARAMETERS:
    State   -   structure stores algorithm state
    S       -   array[N], non-zero scaling coefficients
                S[i] may be negative, sign doesn't matter.

  -- ALGLIB --
     Copyright 14.01.2011 by Bochkanov Sergey
*************************************************************************/
void minbleicsetscale(const minbleicstate &state, const real_1d_array &s)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleicsetscale(const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(s.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Modification of the preconditioner: preconditioning is turned off.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetprecdefault(const minbleicstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleicsetprecdefault(const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Modification  of  the  preconditioner:  diagonal of approximate Hessian is
used.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    D       -   diagonal of the approximate Hessian, array[0..N-1],
                (if larger, only leading N elements are used).

NOTE 1: D[i] should be positive. Exception will be thrown otherwise.

NOTE 2: you should pass diagonal of approximate Hessian - NOT ITS INVERSE.

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetprecdiag(const minbleicstate &state, const real_1d_array &d)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleicsetprecdiag(const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(d.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Modification of the preconditioner: scale-based diagonal preconditioning.

This preconditioning mode can be useful when you  don't  have  approximate
diagonal of Hessian, but you know that your  variables  are  badly  scaled
(for  example,  one  variable is in [1,10], and another in [1000,100000]),
and most part of the ill-conditioning comes from different scales of vars.

In this case simple  scale-based  preconditioner,  with H[i] = 1/(s[i]^2),
can greatly improve convergence.

IMPRTANT: you should set scale of your variables  with  MinBLEICSetScale()
call  (before  or after MinBLEICSetPrecScale() call). Without knowledge of
the scale of your variables scale-based preconditioner will be  just  unit
matrix.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetprecscale(const minbleicstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleicsetprecscale(const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function allows to stop algorithm after specified number of inner
iterations.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    MaxIts  -   maximum number of inner iterations.
                If MaxIts=0, the number of iterations is unlimited.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetmaxits(const minbleicstate &state, const ae_int_t maxits)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleicsetmaxits(const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), maxits, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    NeedXRep-   whether iteration reports are needed or not

If NeedXRep is True, algorithm will call rep() callback function if  it is
provided to MinBLEICOptimize().

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetxrep(const minbleicstate &state, const bool needxrep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleicsetxrep(const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), needxrep, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets maximum step length

IMPORTANT: this feature is hard to combine with preconditioning. You can't
set upper limit on step length, when you solve optimization  problem  with
linear (non-boundary) constraints AND preconditioner turned on.

When  non-boundary  constraints  are  present,  you  have to either a) use
preconditioner, or b) use upper limit on step length.  YOU CAN'T USE BOTH!
In this case algorithm will terminate with appropriate error code.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    StpMax  -   maximum step length, >=0. Set StpMax to 0.0,  if you don't
                want to limit step length.

Use this subroutine when you optimize target function which contains exp()
or  other  fast  growing  functions,  and optimization algorithm makes too
large  steps  which  lead   to overflow. This function allows us to reject
steps  that  are  too  large  (and  therefore  expose  us  to the possible
overflow) without actually calculating function value at the x+stp*d.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetstpmax(const minbleicstate &state, const double stpmax)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleicsetstpmax(const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), stpmax, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function provides reverse communication interface
Reverse communication interface is not documented or recommended to use.
See below for functions which provide better documented API
*************************************************************************/
bool minbleiciteration(const minbleicstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        ae_bool result = alglib_impl::minbleiciteration(const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<bool*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}


void minbleicoptimize(minbleicstate &state,
    void (*func)(const real_1d_array &x, double &func, void *ptr),
    void  (*rep)(const real_1d_array &x, double func, void *ptr), 
    void *ptr)
{
    alglib_impl::ae_state _alglib_env_state;
    if( func==NULL )
        throw ap_error("ALGLIB: error in 'minbleicoptimize()' (func is NULL)");
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        while( alglib_impl::minbleiciteration(state.c_ptr(), &_alglib_env_state) )
        {
            if( state.needf )
            {
                func(state.x, state.f, ptr);
                continue;
            }
            if( state.xupdated )
            {
                if( rep!=NULL )
                    rep(state.x, state.f, ptr);
                continue;
            }
            throw ap_error("ALGLIB: error in 'minbleicoptimize' (some derivatives were not provided?)");
        }
        alglib_impl::ae_state_clear(&_alglib_env_state);
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}


void minbleicoptimize(minbleicstate &state,
    void (*grad)(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr),
    void  (*rep)(const real_1d_array &x, double func, void *ptr), 
    void *ptr)
{
    alglib_impl::ae_state _alglib_env_state;
    if( grad==NULL )
        throw ap_error("ALGLIB: error in 'minbleicoptimize()' (grad is NULL)");
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        while( alglib_impl::minbleiciteration(state.c_ptr(), &_alglib_env_state) )
        {
            if( state.needfg )
            {
                grad(state.x, state.f, state.g, ptr);
                continue;
            }
            if( state.xupdated )
            {
                if( rep!=NULL )
                    rep(state.x, state.f, ptr);
                continue;
            }
            throw ap_error("ALGLIB: error in 'minbleicoptimize' (some derivatives were not provided?)");
        }
        alglib_impl::ae_state_clear(&_alglib_env_state);
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}



/*************************************************************************
BLEIC results

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    X       -   array[0..N-1], solution
    Rep     -   optimization report. You should check Rep.TerminationType
                in  order  to  distinguish  successful  termination  from
                unsuccessful one.
                More information about fields of this  structure  can  be
                found in the comments on MinBLEICReport datatype.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicresults(const minbleicstate &state, real_1d_array &x, minbleicreport &rep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleicresults(const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minbleicreport*>(rep.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
BLEIC results

Buffered implementation of MinBLEICResults() which uses pre-allocated buffer
to store X[]. If buffer size is  too  small,  it  resizes  buffer.  It  is
intended to be used in the inner cycles of performance critical algorithms
where array reallocation penalty is too large to be ignored.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicresultsbuf(const minbleicstate &state, real_1d_array &x, minbleicreport &rep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleicresultsbuf(const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minbleicreport*>(rep.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This subroutine restarts algorithm from new point.
All optimization parameters (including constraints) are left unchanged.

This  function  allows  to  solve multiple  optimization  problems  (which
must have  same number of dimensions) without object reallocation penalty.

INPUT PARAMETERS:
    State   -   structure previously allocated with MinBLEICCreate call.
    X       -   new starting point.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicrestartfrom(const minbleicstate &state, const real_1d_array &x)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleicrestartfrom(const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This  subroutine  turns  on  verification  of  the  user-supplied analytic
gradient:
* user calls this subroutine before optimization begins
* MinBLEICOptimize() is called
* prior to  actual  optimization, for each component  of  parameters being
  optimized X[i] algorithm performs following steps:
  * two trial steps are made to X[i]-TestStep*S[i] and X[i]+TestStep*S[i],
    where X[i] is i-th component of the initial point and S[i] is a  scale
    of i-th parameter
  * if needed, steps are bounded with respect to constraints on X[]
  * F(X) is evaluated at these trial points
  * we perform one more evaluation in the middle point of the interval
  * we  build  cubic  model using function values and derivatives at trial
    points and we compare its prediction with actual value in  the  middle
    point
  * in case difference between prediction and actual value is higher  than
    some predetermined threshold, algorithm stops with completion code -7;
    Rep.VarIdx is set to index of the parameter with incorrect derivative.
* after verification is over, algorithm proceeds to the actual optimization.

NOTE 1: verification  needs  N (parameters count) gradient evaluations. It
        is very costly and you should use  it  only  for  low  dimensional
        problems,  when  you  want  to  be  sure  that  you've   correctly
        calculated  analytic  derivatives.  You  should  not use it in the
        production code (unless you want to check derivatives provided  by
        some third party).

NOTE 2: you  should  carefully  choose  TestStep. Value which is too large
        (so large that function behaviour is significantly non-cubic) will
        lead to false alarms. You may use  different  step  for  different
        parameters by means of setting scale with MinBLEICSetScale().

NOTE 3: this function may lead to false positives. In case it reports that
        I-th  derivative was calculated incorrectly, you may decrease test
        step  and  try  one  more  time  - maybe your function changes too
        sharply  and  your  step  is  too  large for such rapidly chanding
        function.

INPUT PARAMETERS:
    State       -   structure used to store algorithm state
    TestStep    -   verification step:
                    * TestStep=0 turns verification off
                    * TestStep>0 activates verification

  -- ALGLIB --
     Copyright 15.06.2012 by Bochkanov Sergey
*************************************************************************/
void minbleicsetgradientcheck(const minbleicstate &state, const double teststep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleicsetgradientcheck(const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), teststep, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************

*************************************************************************/
_minlbfgsstate_owner::_minlbfgsstate_owner()
{
    p_struct = (alglib_impl::minlbfgsstate*)alglib_impl::ae_malloc(sizeof(alglib_impl::minlbfgsstate), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_minlbfgsstate_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_minlbfgsstate_owner::_minlbfgsstate_owner(const _minlbfgsstate_owner &rhs)
{
    p_struct = (alglib_impl::minlbfgsstate*)alglib_impl::ae_malloc(sizeof(alglib_impl::minlbfgsstate), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_minlbfgsstate_init_copy(p_struct, const_cast<alglib_impl::minlbfgsstate*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_minlbfgsstate_owner& _minlbfgsstate_owner::operator=(const _minlbfgsstate_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_minlbfgsstate_clear(p_struct);
    if( !alglib_impl::_minlbfgsstate_init_copy(p_struct, const_cast<alglib_impl::minlbfgsstate*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_minlbfgsstate_owner::~_minlbfgsstate_owner()
{
    alglib_impl::_minlbfgsstate_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::minlbfgsstate* _minlbfgsstate_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::minlbfgsstate* _minlbfgsstate_owner::c_ptr() const
{
    return const_cast<alglib_impl::minlbfgsstate*>(p_struct);
}
minlbfgsstate::minlbfgsstate() : _minlbfgsstate_owner() ,needf(p_struct->needf),needfg(p_struct->needfg),xupdated(p_struct->xupdated),f(p_struct->f),g(&p_struct->g),x(&p_struct->x)
{
}

minlbfgsstate::minlbfgsstate(const minlbfgsstate &rhs):_minlbfgsstate_owner(rhs) ,needf(p_struct->needf),needfg(p_struct->needfg),xupdated(p_struct->xupdated),f(p_struct->f),g(&p_struct->g),x(&p_struct->x)
{
}

minlbfgsstate& minlbfgsstate::operator=(const minlbfgsstate &rhs)
{
    if( this==&rhs )
        return *this;
    _minlbfgsstate_owner::operator=(rhs);
    return *this;
}

minlbfgsstate::~minlbfgsstate()
{
}


/*************************************************************************

*************************************************************************/
_minlbfgsreport_owner::_minlbfgsreport_owner()
{
    p_struct = (alglib_impl::minlbfgsreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::minlbfgsreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_minlbfgsreport_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_minlbfgsreport_owner::_minlbfgsreport_owner(const _minlbfgsreport_owner &rhs)
{
    p_struct = (alglib_impl::minlbfgsreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::minlbfgsreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_minlbfgsreport_init_copy(p_struct, const_cast<alglib_impl::minlbfgsreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_minlbfgsreport_owner& _minlbfgsreport_owner::operator=(const _minlbfgsreport_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_minlbfgsreport_clear(p_struct);
    if( !alglib_impl::_minlbfgsreport_init_copy(p_struct, const_cast<alglib_impl::minlbfgsreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_minlbfgsreport_owner::~_minlbfgsreport_owner()
{
    alglib_impl::_minlbfgsreport_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::minlbfgsreport* _minlbfgsreport_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::minlbfgsreport* _minlbfgsreport_owner::c_ptr() const
{
    return const_cast<alglib_impl::minlbfgsreport*>(p_struct);
}
minlbfgsreport::minlbfgsreport() : _minlbfgsreport_owner() ,iterationscount(p_struct->iterationscount),nfev(p_struct->nfev),varidx(p_struct->varidx),terminationtype(p_struct->terminationtype)
{
}

minlbfgsreport::minlbfgsreport(const minlbfgsreport &rhs):_minlbfgsreport_owner(rhs) ,iterationscount(p_struct->iterationscount),nfev(p_struct->nfev),varidx(p_struct->varidx),terminationtype(p_struct->terminationtype)
{
}

minlbfgsreport& minlbfgsreport::operator=(const minlbfgsreport &rhs)
{
    if( this==&rhs )
        return *this;
    _minlbfgsreport_owner::operator=(rhs);
    return *this;
}

minlbfgsreport::~minlbfgsreport()
{
}

/*************************************************************************
        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION

DESCRIPTION:
The subroutine minimizes function F(x) of N arguments by  using  a  quasi-
Newton method (LBFGS scheme) which is optimized to use  a  minimum  amount
of memory.
The subroutine generates the approximation of an inverse Hessian matrix by
using information about the last M steps of the algorithm  (instead of N).
It lessens a required amount of memory from a value  of  order  N^2  to  a
value of order 2*N*M.


REQUIREMENTS:
Algorithm will request following information during its operation:
* function value F and its gradient G (simultaneously) at given point X


USAGE:
1. User initializes algorithm state with MinLBFGSCreate() call
2. User tunes solver parameters with MinLBFGSSetCond() MinLBFGSSetStpMax()
   and other functions
3. User calls MinLBFGSOptimize() function which takes algorithm  state and
   pointer (delegate, etc.) to callback function which calculates F/G.
4. User calls MinLBFGSResults() to get solution
5. Optionally user may call MinLBFGSRestartFrom() to solve another problem
   with same N/M but another starting point and/or another function.
   MinLBFGSRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   problem dimension. N>0
    M       -   number of corrections in the BFGS scheme of Hessian
                approximation update. Recommended value:  3<=M<=7. The smaller
                value causes worse convergence, the bigger will  not  cause  a
                considerably better convergence, but will cause a fall in  the
                performance. M<=N.
    X       -   initial solution approximation, array[0..N-1].


OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state


NOTES:
1. you may tune stopping conditions with MinLBFGSSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLBFGSSetStpMax() function to bound algorithm's  steps.  However,
   L-BFGS rarely needs such a tuning.


  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgscreate(const ae_int_t n, const ae_int_t m, const real_1d_array &x, minlbfgsstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlbfgscreate(n, m, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minlbfgsstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION

DESCRIPTION:
The subroutine minimizes function F(x) of N arguments by  using  a  quasi-
Newton method (LBFGS scheme) which is optimized to use  a  minimum  amount
of memory.
The subroutine generates the approximation of an inverse Hessian matrix by
using information about the last M steps of the algorithm  (instead of N).
It lessens a required amount of memory from a value  of  order  N^2  to  a
value of order 2*N*M.


REQUIREMENTS:
Algorithm will request following information during its operation:
* function value F and its gradient G (simultaneously) at given point X


USAGE:
1. User initializes algorithm state with MinLBFGSCreate() call
2. User tunes solver parameters with MinLBFGSSetCond() MinLBFGSSetStpMax()
   and other functions
3. User calls MinLBFGSOptimize() function which takes algorithm  state and
   pointer (delegate, etc.) to callback function which calculates F/G.
4. User calls MinLBFGSResults() to get solution
5. Optionally user may call MinLBFGSRestartFrom() to solve another problem
   with same N/M but another starting point and/or another function.
   MinLBFGSRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   problem dimension. N>0
    M       -   number of corrections in the BFGS scheme of Hessian
                approximation update. Recommended value:  3<=M<=7. The smaller
                value causes worse convergence, the bigger will  not  cause  a
                considerably better convergence, but will cause a fall in  the
                performance. M<=N.
    X       -   initial solution approximation, array[0..N-1].


OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state


NOTES:
1. you may tune stopping conditions with MinLBFGSSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLBFGSSetStpMax() function to bound algorithm's  steps.  However,
   L-BFGS rarely needs such a tuning.


  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgscreate(const ae_int_t m, const real_1d_array &x, minlbfgsstate &state)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;

    n = x.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlbfgscreate(n, m, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minlbfgsstate*>(state.c_ptr()), &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
The subroutine is finite difference variant of MinLBFGSCreate().  It  uses
finite differences in order to differentiate target function.

Description below contains information which is specific to  this function
only. We recommend to read comments on MinLBFGSCreate() in  order  to  get
more information about creation of LBFGS optimizer.

INPUT PARAMETERS:
    N       -   problem dimension, N>0:
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    M       -   number of corrections in the BFGS scheme of Hessian
                approximation update. Recommended value:  3<=M<=7. The smaller
                value causes worse convergence, the bigger will  not  cause  a
                considerably better convergence, but will cause a fall in  the
                performance. M<=N.
    X       -   starting point, array[0..N-1].
    DiffStep-   differentiation step, >0

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTES:
1. algorithm uses 4-point central formula for differentiation.
2. differentiation step along I-th axis is equal to DiffStep*S[I] where
   S[] is scaling vector which can be set by MinLBFGSSetScale() call.
3. we recommend you to use moderate values of  differentiation  step.  Too
   large step will result in too large truncation  errors, while too small
   step will result in too large numerical  errors.  1.0E-6  can  be  good
   value to start with.
4. Numerical  differentiation  is   very   inefficient  -   one   gradient
   calculation needs 4*N function evaluations. This function will work for
   any N - either small (1...10), moderate (10...100) or  large  (100...).
   However, performance penalty will be too severe for any N's except  for
   small ones.
   We should also say that code which relies on numerical  differentiation
   is   less  robust  and  precise.  LBFGS  needs  exact  gradient values.
   Imprecise gradient may slow  down  convergence,  especially  on  highly
   nonlinear problems.
   Thus  we  recommend to use this function for fast prototyping on small-
   dimensional problems only, and to implement analytical gradient as soon
   as possible.

  -- ALGLIB --
     Copyright 16.05.2011 by Bochkanov Sergey
*************************************************************************/
void minlbfgscreatef(const ae_int_t n, const ae_int_t m, const real_1d_array &x, const double diffstep, minlbfgsstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlbfgscreatef(n, m, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), diffstep, const_cast<alglib_impl::minlbfgsstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
The subroutine is finite difference variant of MinLBFGSCreate().  It  uses
finite differences in order to differentiate target function.

Description below contains information which is specific to  this function
only. We recommend to read comments on MinLBFGSCreate() in  order  to  get
more information about creation of LBFGS optimizer.

INPUT PARAMETERS:
    N       -   problem dimension, N>0:
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    M       -   number of corrections in the BFGS scheme of Hessian
                approximation update. Recommended value:  3<=M<=7. The smaller
                value causes worse convergence, the bigger will  not  cause  a
                considerably better convergence, but will cause a fall in  the
                performance. M<=N.
    X       -   starting point, array[0..N-1].
    DiffStep-   differentiation step, >0

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTES:
1. algorithm uses 4-point central formula for differentiation.
2. differentiation step along I-th axis is equal to DiffStep*S[I] where
   S[] is scaling vector which can be set by MinLBFGSSetScale() call.
3. we recommend you to use moderate values of  differentiation  step.  Too
   large step will result in too large truncation  errors, while too small
   step will result in too large numerical  errors.  1.0E-6  can  be  good
   value to start with.
4. Numerical  differentiation  is   very   inefficient  -   one   gradient
   calculation needs 4*N function evaluations. This function will work for
   any N - either small (1...10), moderate (10...100) or  large  (100...).
   However, performance penalty will be too severe for any N's except  for
   small ones.
   We should also say that code which relies on numerical  differentiation
   is   less  robust  and  precise.  LBFGS  needs  exact  gradient values.
   Imprecise gradient may slow  down  convergence,  especially  on  highly
   nonlinear problems.
   Thus  we  recommend to use this function for fast prototyping on small-
   dimensional problems only, and to implement analytical gradient as soon
   as possible.

  -- ALGLIB --
     Copyright 16.05.2011 by Bochkanov Sergey
*************************************************************************/
void minlbfgscreatef(const ae_int_t m, const real_1d_array &x, const double diffstep, minlbfgsstate &state)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;

    n = x.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlbfgscreatef(n, m, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), diffstep, const_cast<alglib_impl::minlbfgsstate*>(state.c_ptr()), &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets stopping conditions for L-BFGS optimization algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    EpsG    -   >=0
                The  subroutine  finishes  its  work   if   the  condition
                |v|<EpsG is satisfied, where:
                * |.| means Euclidian norm
                * v - scaled gradient vector, v[i]=g[i]*s[i]
                * g - gradient
                * s - scaling coefficients set by MinLBFGSSetScale()
    EpsF    -   >=0
                The  subroutine  finishes  its work if on k+1-th iteration
                the  condition  |F(k+1)-F(k)|<=EpsF*max{|F(k)|,|F(k+1)|,1}
                is satisfied.
    EpsX    -   >=0
                The subroutine finishes its work if  on  k+1-th  iteration
                the condition |v|<=EpsX is fulfilled, where:
                * |.| means Euclidian norm
                * v - scaled step vector, v[i]=dx[i]/s[i]
                * dx - ste pvector, dx=X(k+1)-X(k)
                * s - scaling coefficients set by MinLBFGSSetScale()
    MaxIts  -   maximum number of iterations. If MaxIts=0, the  number  of
                iterations is unlimited.

Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to
automatic stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetcond(const minlbfgsstate &state, const double epsg, const double epsf, const double epsx, const ae_int_t maxits)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlbfgssetcond(const_cast<alglib_impl::minlbfgsstate*>(state.c_ptr()), epsg, epsf, epsx, maxits, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    NeedXRep-   whether iteration reports are needed or not

If NeedXRep is True, algorithm will call rep() callback function if  it is
provided to MinLBFGSOptimize().


  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetxrep(const minlbfgsstate &state, const bool needxrep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlbfgssetxrep(const_cast<alglib_impl::minlbfgsstate*>(state.c_ptr()), needxrep, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets maximum step length

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    StpMax  -   maximum step length, >=0. Set StpMax to 0.0 (default),  if
                you don't want to limit step length.

Use this subroutine when you optimize target function which contains exp()
or  other  fast  growing  functions,  and optimization algorithm makes too
large  steps  which  leads  to overflow. This function allows us to reject
steps  that  are  too  large  (and  therefore  expose  us  to the possible
overflow) without actually calculating function value at the x+stp*d.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetstpmax(const minlbfgsstate &state, const double stpmax)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlbfgssetstpmax(const_cast<alglib_impl::minlbfgsstate*>(state.c_ptr()), stpmax, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets scaling coefficients for LBFGS optimizer.

ALGLIB optimizers use scaling matrices to test stopping  conditions  (step
size and gradient are scaled before comparison with tolerances).  Scale of
the I-th variable is a translation invariant measure of:
a) "how large" the variable is
b) how large the step should be to make significant changes in the function

Scaling is also used by finite difference variant of the optimizer  - step
along I-th axis is equal to DiffStep*S[I].

In  most  optimizers  (and  in  the  LBFGS  too)  scaling is NOT a form of
preconditioning. It just  affects  stopping  conditions.  You  should  set
preconditioner  by  separate  call  to  one  of  the  MinLBFGSSetPrec...()
functions.

There  is  special  preconditioning  mode, however,  which  uses   scaling
coefficients to form diagonal preconditioning matrix. You  can  turn  this
mode on, if you want.   But  you should understand that scaling is not the
same thing as preconditioning - these are two different, although  related
forms of tuning solver.

INPUT PARAMETERS:
    State   -   structure stores algorithm state
    S       -   array[N], non-zero scaling coefficients
                S[i] may be negative, sign doesn't matter.

  -- ALGLIB --
     Copyright 14.01.2011 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetscale(const minlbfgsstate &state, const real_1d_array &s)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlbfgssetscale(const_cast<alglib_impl::minlbfgsstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(s.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Modification  of  the  preconditioner:  default  preconditioner    (simple
scaling, same for all elements of X) is used.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTE:  you  can  change  preconditioner  "on  the  fly",  during algorithm
iterations.

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetprecdefault(const minlbfgsstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlbfgssetprecdefault(const_cast<alglib_impl::minlbfgsstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Modification of the preconditioner: Cholesky factorization of  approximate
Hessian is used.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    P       -   triangular preconditioner, Cholesky factorization of
                the approximate Hessian. array[0..N-1,0..N-1],
                (if larger, only leading N elements are used).
    IsUpper -   whether upper or lower triangle of P is given
                (other triangle is not referenced)

After call to this function preconditioner is changed to P  (P  is  copied
into the internal buffer).

NOTE:  you  can  change  preconditioner  "on  the  fly",  during algorithm
iterations.

NOTE 2:  P  should  be nonsingular. Exception will be thrown otherwise.

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetpreccholesky(const minlbfgsstate &state, const real_2d_array &p, const bool isupper)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlbfgssetpreccholesky(const_cast<alglib_impl::minlbfgsstate*>(state.c_ptr()), const_cast<alglib_impl::ae_matrix*>(p.c_ptr()), isupper, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Modification  of  the  preconditioner:  diagonal of approximate Hessian is
used.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    D       -   diagonal of the approximate Hessian, array[0..N-1],
                (if larger, only leading N elements are used).

NOTE:  you  can  change  preconditioner  "on  the  fly",  during algorithm
iterations.

NOTE 2: D[i] should be positive. Exception will be thrown otherwise.

NOTE 3: you should pass diagonal of approximate Hessian - NOT ITS INVERSE.

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetprecdiag(const minlbfgsstate &state, const real_1d_array &d)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlbfgssetprecdiag(const_cast<alglib_impl::minlbfgsstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(d.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Modification of the preconditioner: scale-based diagonal preconditioning.

This preconditioning mode can be useful when you  don't  have  approximate
diagonal of Hessian, but you know that your  variables  are  badly  scaled
(for  example,  one  variable is in [1,10], and another in [1000,100000]),
and most part of the ill-conditioning comes from different scales of vars.

In this case simple  scale-based  preconditioner,  with H[i] = 1/(s[i]^2),
can greatly improve convergence.

IMPRTANT: you should set scale of your variables  with  MinLBFGSSetScale()
call  (before  or after MinLBFGSSetPrecScale() call). Without knowledge of
the scale of your variables scale-based preconditioner will be  just  unit
matrix.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetprecscale(const minlbfgsstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlbfgssetprecscale(const_cast<alglib_impl::minlbfgsstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function provides reverse communication interface
Reverse communication interface is not documented or recommended to use.
See below for functions which provide better documented API
*************************************************************************/
bool minlbfgsiteration(const minlbfgsstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        ae_bool result = alglib_impl::minlbfgsiteration(const_cast<alglib_impl::minlbfgsstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<bool*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}


void minlbfgsoptimize(minlbfgsstate &state,
    void (*func)(const real_1d_array &x, double &func, void *ptr),
    void  (*rep)(const real_1d_array &x, double func, void *ptr), 
    void *ptr)
{
    alglib_impl::ae_state _alglib_env_state;
    if( func==NULL )
        throw ap_error("ALGLIB: error in 'minlbfgsoptimize()' (func is NULL)");
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        while( alglib_impl::minlbfgsiteration(state.c_ptr(), &_alglib_env_state) )
        {
            if( state.needf )
            {
                func(state.x, state.f, ptr);
                continue;
            }
            if( state.xupdated )
            {
                if( rep!=NULL )
                    rep(state.x, state.f, ptr);
                continue;
            }
            throw ap_error("ALGLIB: error in 'minlbfgsoptimize' (some derivatives were not provided?)");
        }
        alglib_impl::ae_state_clear(&_alglib_env_state);
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}


void minlbfgsoptimize(minlbfgsstate &state,
    void (*grad)(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr),
    void  (*rep)(const real_1d_array &x, double func, void *ptr), 
    void *ptr)
{
    alglib_impl::ae_state _alglib_env_state;
    if( grad==NULL )
        throw ap_error("ALGLIB: error in 'minlbfgsoptimize()' (grad is NULL)");
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        while( alglib_impl::minlbfgsiteration(state.c_ptr(), &_alglib_env_state) )
        {
            if( state.needfg )
            {
                grad(state.x, state.f, state.g, ptr);
                continue;
            }
            if( state.xupdated )
            {
                if( rep!=NULL )
                    rep(state.x, state.f, ptr);
                continue;
            }
            throw ap_error("ALGLIB: error in 'minlbfgsoptimize' (some derivatives were not provided?)");
        }
        alglib_impl::ae_state_clear(&_alglib_env_state);
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}



/*************************************************************************
L-BFGS algorithm results

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    X       -   array[0..N-1], solution
    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    * -7    gradient verification failed.
                            See MinLBFGSSetGradientCheck() for more information.
                    * -2    rounding errors prevent further improvement.
                            X contains best point found.
                    * -1    incorrect parameters were specified
                    *  1    relative function improvement is no more than
                            EpsF.
                    *  2    relative step is no more than EpsX.
                    *  4    gradient norm is no more than EpsG
                    *  5    MaxIts steps was taken
                    *  7    stopping conditions are too stringent,
                            further improvement is impossible
                * Rep.IterationsCount contains iterations count
                * NFEV countains number of function calculations

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgsresults(const minlbfgsstate &state, real_1d_array &x, minlbfgsreport &rep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlbfgsresults(const_cast<alglib_impl::minlbfgsstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minlbfgsreport*>(rep.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
L-BFGS algorithm results

Buffered implementation of MinLBFGSResults which uses pre-allocated buffer
to store X[]. If buffer size is  too  small,  it  resizes  buffer.  It  is
intended to be used in the inner cycles of performance critical algorithms
where array reallocation penalty is too large to be ignored.

  -- ALGLIB --
     Copyright 20.08.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgsresultsbuf(const minlbfgsstate &state, real_1d_array &x, minlbfgsreport &rep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlbfgsresultsbuf(const_cast<alglib_impl::minlbfgsstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minlbfgsreport*>(rep.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This  subroutine restarts LBFGS algorithm from new point. All optimization
parameters are left unchanged.

This  function  allows  to  solve multiple  optimization  problems  (which
must have same number of dimensions) without object reallocation penalty.

INPUT PARAMETERS:
    State   -   structure used to store algorithm state
    X       -   new starting point.

  -- ALGLIB --
     Copyright 30.07.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgsrestartfrom(const minlbfgsstate &state, const real_1d_array &x)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlbfgsrestartfrom(const_cast<alglib_impl::minlbfgsstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This  subroutine  turns  on  verification  of  the  user-supplied analytic
gradient:
* user calls this subroutine before optimization begins
* MinLBFGSOptimize() is called
* prior to  actual  optimization, for each component  of  parameters being
  optimized X[i] algorithm performs following steps:
  * two trial steps are made to X[i]-TestStep*S[i] and X[i]+TestStep*S[i],
    where X[i] is i-th component of the initial point and S[i] is a  scale
    of i-th parameter
  * if needed, steps are bounded with respect to constraints on X[]
  * F(X) is evaluated at these trial points
  * we perform one more evaluation in the middle point of the interval
  * we  build  cubic  model using function values and derivatives at trial
    points and we compare its prediction with actual value in  the  middle
    point
  * in case difference between prediction and actual value is higher  than
    some predetermined threshold, algorithm stops with completion code -7;
    Rep.VarIdx is set to index of the parameter with incorrect derivative.
* after verification is over, algorithm proceeds to the actual optimization.

NOTE 1: verification  needs  N (parameters count) gradient evaluations. It
        is very costly and you should use  it  only  for  low  dimensional
        problems,  when  you  want  to  be  sure  that  you've   correctly
        calculated  analytic  derivatives.  You  should  not use it in the
        production code (unless you want to check derivatives provided  by
        some third party).

NOTE 2: you  should  carefully  choose  TestStep. Value which is too large
        (so large that function behaviour is significantly non-cubic) will
        lead to false alarms. You may use  different  step  for  different
        parameters by means of setting scale with MinLBFGSSetScale().

NOTE 3: this function may lead to false positives. In case it reports that
        I-th  derivative was calculated incorrectly, you may decrease test
        step  and  try  one  more  time  - maybe your function changes too
        sharply  and  your  step  is  too  large for such rapidly chanding
        function.

INPUT PARAMETERS:
    State       -   structure used to store algorithm state
    TestStep    -   verification step:
                    * TestStep=0 turns verification off
                    * TestStep>0 activates verification

  -- ALGLIB --
     Copyright 24.05.2012 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetgradientcheck(const minlbfgsstate &state, const double teststep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlbfgssetgradientcheck(const_cast<alglib_impl::minlbfgsstate*>(state.c_ptr()), teststep, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}



/*************************************************************************
This object stores nonlinear optimizer state.
You should use functions provided by MinQP subpackage to work with this
object
*************************************************************************/
_minqpstate_owner::_minqpstate_owner()
{
    p_struct = (alglib_impl::minqpstate*)alglib_impl::ae_malloc(sizeof(alglib_impl::minqpstate), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_minqpstate_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_minqpstate_owner::_minqpstate_owner(const _minqpstate_owner &rhs)
{
    p_struct = (alglib_impl::minqpstate*)alglib_impl::ae_malloc(sizeof(alglib_impl::minqpstate), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_minqpstate_init_copy(p_struct, const_cast<alglib_impl::minqpstate*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_minqpstate_owner& _minqpstate_owner::operator=(const _minqpstate_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_minqpstate_clear(p_struct);
    if( !alglib_impl::_minqpstate_init_copy(p_struct, const_cast<alglib_impl::minqpstate*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_minqpstate_owner::~_minqpstate_owner()
{
    alglib_impl::_minqpstate_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::minqpstate* _minqpstate_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::minqpstate* _minqpstate_owner::c_ptr() const
{
    return const_cast<alglib_impl::minqpstate*>(p_struct);
}
minqpstate::minqpstate() : _minqpstate_owner() 
{
}

minqpstate::minqpstate(const minqpstate &rhs):_minqpstate_owner(rhs) 
{
}

minqpstate& minqpstate::operator=(const minqpstate &rhs)
{
    if( this==&rhs )
        return *this;
    _minqpstate_owner::operator=(rhs);
    return *this;
}

minqpstate::~minqpstate()
{
}


/*************************************************************************
This structure stores optimization report:
* InnerIterationsCount      number of inner iterations
* OuterIterationsCount      number of outer iterations
* NCholesky                 number of Cholesky decomposition
* NMV                       number of matrix-vector products
                            (only products calculated as part of iterative
                            process are counted)
* TerminationType           completion code (see below)

Completion codes:
* -5    inappropriate solver was used:
        * Cholesky solver for semidefinite or indefinite problems
        * Cholesky solver for problems with non-boundary constraints
* -3    inconsistent constraints (or, maybe, feasible point is
        too hard to find). If you are sure that constraints are feasible,
        try to restart optimizer with better initial approximation.
* -1    solver error
*  4    successful completion
*  5    MaxIts steps was taken
*  7    stopping conditions are too stringent,
        further improvement is impossible,
        X contains best point found so far.
*************************************************************************/
_minqpreport_owner::_minqpreport_owner()
{
    p_struct = (alglib_impl::minqpreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::minqpreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_minqpreport_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_minqpreport_owner::_minqpreport_owner(const _minqpreport_owner &rhs)
{
    p_struct = (alglib_impl::minqpreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::minqpreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_minqpreport_init_copy(p_struct, const_cast<alglib_impl::minqpreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_minqpreport_owner& _minqpreport_owner::operator=(const _minqpreport_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_minqpreport_clear(p_struct);
    if( !alglib_impl::_minqpreport_init_copy(p_struct, const_cast<alglib_impl::minqpreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_minqpreport_owner::~_minqpreport_owner()
{
    alglib_impl::_minqpreport_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::minqpreport* _minqpreport_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::minqpreport* _minqpreport_owner::c_ptr() const
{
    return const_cast<alglib_impl::minqpreport*>(p_struct);
}
minqpreport::minqpreport() : _minqpreport_owner() ,inneriterationscount(p_struct->inneriterationscount),outeriterationscount(p_struct->outeriterationscount),nmv(p_struct->nmv),ncholesky(p_struct->ncholesky),terminationtype(p_struct->terminationtype)
{
}

minqpreport::minqpreport(const minqpreport &rhs):_minqpreport_owner(rhs) ,inneriterationscount(p_struct->inneriterationscount),outeriterationscount(p_struct->outeriterationscount),nmv(p_struct->nmv),ncholesky(p_struct->ncholesky),terminationtype(p_struct->terminationtype)
{
}

minqpreport& minqpreport::operator=(const minqpreport &rhs)
{
    if( this==&rhs )
        return *this;
    _minqpreport_owner::operator=(rhs);
    return *this;
}

minqpreport::~minqpreport()
{
}

/*************************************************************************
                    CONSTRAINED QUADRATIC PROGRAMMING

The subroutine creates QP optimizer. After initial creation,  it  contains
default optimization problem with zero quadratic and linear terms  and  no
constraints. You should set quadratic/linear terms with calls to functions
provided by MinQP subpackage.

INPUT PARAMETERS:
    N       -   problem size

OUTPUT PARAMETERS:
    State   -   optimizer with zero quadratic/linear terms
                and no constraints

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpcreate(const ae_int_t n, minqpstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minqpcreate(n, const_cast<alglib_impl::minqpstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets linear term for QP solver.

By default, linear term is zero.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    B       -   linear term, array[N].

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpsetlinearterm(const minqpstate &state, const real_1d_array &b)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minqpsetlinearterm(const_cast<alglib_impl::minqpstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(b.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets quadratic term for QP solver.

By default quadratic term is zero.

IMPORTANT: this solver minimizes following  function:
    f(x) = 0.5*x'*A*x + b'*x.
Note that quadratic term has 0.5 before it. So if  you  want  to  minimize
    f(x) = x^2 + x
you should rewrite your problem as follows:
    f(x) = 0.5*(2*x^2) + x
and your matrix A will be equal to [[2.0]], not to [[1.0]]

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    A       -   matrix, array[N,N]
    IsUpper -   (optional) storage type:
                * if True, symmetric matrix  A  is  given  by  its  upper
                  triangle, and the lower triangle isn’t used
                * if False, symmetric matrix  A  is  given  by  its lower
                  triangle, and the upper triangle isn’t used
                * if not given, both lower and upper  triangles  must  be
                  filled.

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpsetquadraticterm(const minqpstate &state, const real_2d_array &a, const bool isupper)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minqpsetquadraticterm(const_cast<alglib_impl::minqpstate*>(state.c_ptr()), const_cast<alglib_impl::ae_matrix*>(a.c_ptr()), isupper, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets quadratic term for QP solver.

By default quadratic term is zero.

IMPORTANT: this solver minimizes following  function:
    f(x) = 0.5*x'*A*x + b'*x.
Note that quadratic term has 0.5 before it. So if  you  want  to  minimize
    f(x) = x^2 + x
you should rewrite your problem as follows:
    f(x) = 0.5*(2*x^2) + x
and your matrix A will be equal to [[2.0]], not to [[1.0]]

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    A       -   matrix, array[N,N]
    IsUpper -   (optional) storage type:
                * if True, symmetric matrix  A  is  given  by  its  upper
                  triangle, and the lower triangle isn’t used
                * if False, symmetric matrix  A  is  given  by  its lower
                  triangle, and the upper triangle isn’t used
                * if not given, both lower and upper  triangles  must  be
                  filled.

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpsetquadraticterm(const minqpstate &state, const real_2d_array &a)
{
    alglib_impl::ae_state _alglib_env_state;    
    bool isupper;
    if( !alglib_impl::ae_is_symmetric(const_cast<alglib_impl::ae_matrix*>(a.c_ptr())) )
        throw ap_error("'a' parameter is not symmetric matrix");
    isupper = false;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minqpsetquadraticterm(const_cast<alglib_impl::minqpstate*>(state.c_ptr()), const_cast<alglib_impl::ae_matrix*>(a.c_ptr()), isupper, &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets starting point for QP solver. It is useful to have
good initial approximation to the solution, because it will increase
speed of convergence and identification of active constraints.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    X       -   starting point, array[N].

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpsetstartingpoint(const minqpstate &state, const real_1d_array &x)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minqpsetstartingpoint(const_cast<alglib_impl::minqpstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This  function sets origin for QP solver. By default, following QP program
is solved:

    min(0.5*x'*A*x+b'*x)

This function allows to solve different problem:

    min(0.5*(x-x_origin)'*A*(x-x_origin)+b'*(x-x_origin))

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    XOrigin -   origin, array[N].

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpsetorigin(const minqpstate &state, const real_1d_array &xorigin)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minqpsetorigin(const_cast<alglib_impl::minqpstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(xorigin.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function tells solver to use Cholesky-based algorithm.

Cholesky-based algorithm can be used when:
* problem is convex
* there is no constraints or only boundary constraints are present

This algorithm has O(N^3) complexity for unconstrained problem and  is  up
to several times slower on bound constrained  problems  (these  additional
iterations are needed to identify active constraints).

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpsetalgocholesky(const minqpstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minqpsetalgocholesky(const_cast<alglib_impl::minqpstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets boundary constraints for QP solver

Boundary constraints are inactive by default (after initial creation).
After  being  set,  they  are  preserved  until explicitly turned off with
another SetBC() call.

INPUT PARAMETERS:
    State   -   structure stores algorithm state
    BndL    -   lower bounds, array[N].
                If some (all) variables are unbounded, you may specify
                very small number or -INF (latter is recommended because
                it will allow solver to use better algorithm).
    BndU    -   upper bounds, array[N].
                If some (all) variables are unbounded, you may specify
                very large number or +INF (latter is recommended because
                it will allow solver to use better algorithm).

NOTE: it is possible to specify BndL[i]=BndU[i]. In this case I-th
variable will be "frozen" at X[i]=BndL[i]=BndU[i].

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpsetbc(const minqpstate &state, const real_1d_array &bndl, const real_1d_array &bndu)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minqpsetbc(const_cast<alglib_impl::minqpstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(bndl.c_ptr()), const_cast<alglib_impl::ae_vector*>(bndu.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets linear constraints for QP optimizer.

Linear constraints are inactive by default (after initial creation).

INPUT PARAMETERS:
    State   -   structure previously allocated with MinQPCreate call.
    C       -   linear constraints, array[K,N+1].
                Each row of C represents one constraint, either equality
                or inequality (see below):
                * first N elements correspond to coefficients,
                * last element corresponds to the right part.
                All elements of C (including right part) must be finite.
    CT      -   type of constraints, array[K]:
                * if CT[i]>0, then I-th constraint is C[i,*]*x >= C[i,n+1]
                * if CT[i]=0, then I-th constraint is C[i,*]*x  = C[i,n+1]
                * if CT[i]<0, then I-th constraint is C[i,*]*x <= C[i,n+1]
    K       -   number of equality/inequality constraints, K>=0:
                * if given, only leading K elements of C/CT are used
                * if not given, automatically determined from sizes of C/CT

NOTE 1: linear (non-bound) constraints are satisfied only approximately  -
        there always exists some minor violation (about 10^10...10^13) due
        to rounding errors.

  -- ALGLIB --
     Copyright 19.06.2012 by Bochkanov Sergey
*************************************************************************/
void minqpsetlc(const minqpstate &state, const real_2d_array &c, const integer_1d_array &ct, const ae_int_t k)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minqpsetlc(const_cast<alglib_impl::minqpstate*>(state.c_ptr()), const_cast<alglib_impl::ae_matrix*>(c.c_ptr()), const_cast<alglib_impl::ae_vector*>(ct.c_ptr()), k, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets linear constraints for QP optimizer.

Linear constraints are inactive by default (after initial creation).

INPUT PARAMETERS:
    State   -   structure previously allocated with MinQPCreate call.
    C       -   linear constraints, array[K,N+1].
                Each row of C represents one constraint, either equality
                or inequality (see below):
                * first N elements correspond to coefficients,
                * last element corresponds to the right part.
                All elements of C (including right part) must be finite.
    CT      -   type of constraints, array[K]:
                * if CT[i]>0, then I-th constraint is C[i,*]*x >= C[i,n+1]
                * if CT[i]=0, then I-th constraint is C[i,*]*x  = C[i,n+1]
                * if CT[i]<0, then I-th constraint is C[i,*]*x <= C[i,n+1]
    K       -   number of equality/inequality constraints, K>=0:
                * if given, only leading K elements of C/CT are used
                * if not given, automatically determined from sizes of C/CT

NOTE 1: linear (non-bound) constraints are satisfied only approximately  -
        there always exists some minor violation (about 10^10...10^13) due
        to rounding errors.

  -- ALGLIB --
     Copyright 19.06.2012 by Bochkanov Sergey
*************************************************************************/
void minqpsetlc(const minqpstate &state, const real_2d_array &c, const integer_1d_array &ct)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t k;
    if( (c.rows()!=ct.length()))
        throw ap_error("Error while calling 'minqpsetlc': looks like one of arguments has wrong size");
    k = c.rows();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minqpsetlc(const_cast<alglib_impl::minqpstate*>(state.c_ptr()), const_cast<alglib_impl::ae_matrix*>(c.c_ptr()), const_cast<alglib_impl::ae_vector*>(ct.c_ptr()), k, &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function solves quadratic programming problem.
You should call it after setting solver options with MinQPSet...() calls.

INPUT PARAMETERS:
    State   -   algorithm state

You should use MinQPResults() function to access results after calls
to this function.

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpoptimize(const minqpstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minqpoptimize(const_cast<alglib_impl::minqpstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
QP solver results

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    X       -   array[0..N-1], solution
    Rep     -   optimization report. You should check Rep.TerminationType,
                which contains completion code, and you may check  another
                fields which contain another information  about  algorithm
                functioning.

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpresults(const minqpstate &state, real_1d_array &x, minqpreport &rep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minqpresults(const_cast<alglib_impl::minqpstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minqpreport*>(rep.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
QP results

Buffered implementation of MinQPResults() which uses pre-allocated  buffer
to store X[]. If buffer size is  too  small,  it  resizes  buffer.  It  is
intended to be used in the inner cycles of performance critical algorithms
where array reallocation penalty is too large to be ignored.

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpresultsbuf(const minqpstate &state, real_1d_array &x, minqpreport &rep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minqpresultsbuf(const_cast<alglib_impl::minqpstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minqpreport*>(rep.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Levenberg-Marquardt optimizer.

This structure should be created using one of the MinLMCreate???()
functions. You should not access its fields directly; use ALGLIB functions
to work with it.
*************************************************************************/
_minlmstate_owner::_minlmstate_owner()
{
    p_struct = (alglib_impl::minlmstate*)alglib_impl::ae_malloc(sizeof(alglib_impl::minlmstate), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_minlmstate_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_minlmstate_owner::_minlmstate_owner(const _minlmstate_owner &rhs)
{
    p_struct = (alglib_impl::minlmstate*)alglib_impl::ae_malloc(sizeof(alglib_impl::minlmstate), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_minlmstate_init_copy(p_struct, const_cast<alglib_impl::minlmstate*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_minlmstate_owner& _minlmstate_owner::operator=(const _minlmstate_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_minlmstate_clear(p_struct);
    if( !alglib_impl::_minlmstate_init_copy(p_struct, const_cast<alglib_impl::minlmstate*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_minlmstate_owner::~_minlmstate_owner()
{
    alglib_impl::_minlmstate_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::minlmstate* _minlmstate_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::minlmstate* _minlmstate_owner::c_ptr() const
{
    return const_cast<alglib_impl::minlmstate*>(p_struct);
}
minlmstate::minlmstate() : _minlmstate_owner() ,needf(p_struct->needf),needfg(p_struct->needfg),needfgh(p_struct->needfgh),needfi(p_struct->needfi),needfij(p_struct->needfij),xupdated(p_struct->xupdated),f(p_struct->f),fi(&p_struct->fi),g(&p_struct->g),h(&p_struct->h),j(&p_struct->j),x(&p_struct->x)
{
}

minlmstate::minlmstate(const minlmstate &rhs):_minlmstate_owner(rhs) ,needf(p_struct->needf),needfg(p_struct->needfg),needfgh(p_struct->needfgh),needfi(p_struct->needfi),needfij(p_struct->needfij),xupdated(p_struct->xupdated),f(p_struct->f),fi(&p_struct->fi),g(&p_struct->g),h(&p_struct->h),j(&p_struct->j),x(&p_struct->x)
{
}

minlmstate& minlmstate::operator=(const minlmstate &rhs)
{
    if( this==&rhs )
        return *this;
    _minlmstate_owner::operator=(rhs);
    return *this;
}

minlmstate::~minlmstate()
{
}


/*************************************************************************
Optimization report, filled by MinLMResults() function

FIELDS:
* TerminationType, completetion code:
    * -7    derivative correctness check failed;
            see Rep.WrongNum, Rep.WrongI, Rep.WrongJ for
            more information.
    *  1    relative function improvement is no more than
            EpsF.
    *  2    relative step is no more than EpsX.
    *  4    gradient is no more than EpsG.
    *  5    MaxIts steps was taken
    *  7    stopping conditions are too stringent,
            further improvement is impossible
* IterationsCount, contains iterations count
* NFunc, number of function calculations
* NJac, number of Jacobi matrix calculations
* NGrad, number of gradient calculations
* NHess, number of Hessian calculations
* NCholesky, number of Cholesky decomposition calculations
*************************************************************************/
_minlmreport_owner::_minlmreport_owner()
{
    p_struct = (alglib_impl::minlmreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::minlmreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_minlmreport_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_minlmreport_owner::_minlmreport_owner(const _minlmreport_owner &rhs)
{
    p_struct = (alglib_impl::minlmreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::minlmreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_minlmreport_init_copy(p_struct, const_cast<alglib_impl::minlmreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_minlmreport_owner& _minlmreport_owner::operator=(const _minlmreport_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_minlmreport_clear(p_struct);
    if( !alglib_impl::_minlmreport_init_copy(p_struct, const_cast<alglib_impl::minlmreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_minlmreport_owner::~_minlmreport_owner()
{
    alglib_impl::_minlmreport_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::minlmreport* _minlmreport_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::minlmreport* _minlmreport_owner::c_ptr() const
{
    return const_cast<alglib_impl::minlmreport*>(p_struct);
}
minlmreport::minlmreport() : _minlmreport_owner() ,iterationscount(p_struct->iterationscount),terminationtype(p_struct->terminationtype),funcidx(p_struct->funcidx),varidx(p_struct->varidx),nfunc(p_struct->nfunc),njac(p_struct->njac),ngrad(p_struct->ngrad),nhess(p_struct->nhess),ncholesky(p_struct->ncholesky)
{
}

minlmreport::minlmreport(const minlmreport &rhs):_minlmreport_owner(rhs) ,iterationscount(p_struct->iterationscount),terminationtype(p_struct->terminationtype),funcidx(p_struct->funcidx),varidx(p_struct->varidx),nfunc(p_struct->nfunc),njac(p_struct->njac),ngrad(p_struct->ngrad),nhess(p_struct->nhess),ncholesky(p_struct->ncholesky)
{
}

minlmreport& minlmreport::operator=(const minlmreport &rhs)
{
    if( this==&rhs )
        return *this;
    _minlmreport_owner::operator=(rhs);
    return *this;
}

minlmreport::~minlmreport()
{
}

/*************************************************************************
                IMPROVED LEVENBERG-MARQUARDT METHOD FOR
                 NON-LINEAR LEAST SQUARES OPTIMIZATION

DESCRIPTION:
This function is used to find minimum of function which is represented  as
sum of squares:
    F(x) = f[0]^2(x[0],...,x[n-1]) + ... + f[m-1]^2(x[0],...,x[n-1])
using value of function vector f[] and Jacobian of f[].


REQUIREMENTS:
This algorithm will request following information during its operation:

* function vector f[] at given point X
* function vector f[] and Jacobian of f[] (simultaneously) at given point

There are several overloaded versions of  MinLMOptimize()  function  which
correspond  to  different LM-like optimization algorithms provided by this
unit. You should choose version which accepts fvec()  and jac() callbacks.
First  one  is used to calculate f[] at given point, second one calculates
f[] and Jacobian df[i]/dx[j].

You can try to initialize MinLMState structure with VJ  function and  then
use incorrect version  of  MinLMOptimize()  (for  example,  version  which
works  with  general  form function and does not provide Jacobian), but it
will  lead  to  exception  being  thrown  after first attempt to calculate
Jacobian.


USAGE:
1. User initializes algorithm state with MinLMCreateVJ() call
2. User tunes solver parameters with MinLMSetCond(),  MinLMSetStpMax() and
   other functions
3. User calls MinLMOptimize() function which  takes algorithm  state   and
   callback functions.
4. User calls MinLMResults() to get solution
5. Optionally, user may call MinLMRestartFrom() to solve  another  problem
   with same N/M but another starting point and/or another function.
   MinLMRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   dimension, N>1
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    M       -   number of functions f[i]
    X       -   initial solution, array[0..N-1]

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTES:
1. you may tune stopping conditions with MinLMSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLMSetStpMax() function to bound algorithm's steps.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatevj(const ae_int_t n, const ae_int_t m, const real_1d_array &x, minlmstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmcreatevj(n, m, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minlmstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
                IMPROVED LEVENBERG-MARQUARDT METHOD FOR
                 NON-LINEAR LEAST SQUARES OPTIMIZATION

DESCRIPTION:
This function is used to find minimum of function which is represented  as
sum of squares:
    F(x) = f[0]^2(x[0],...,x[n-1]) + ... + f[m-1]^2(x[0],...,x[n-1])
using value of function vector f[] and Jacobian of f[].


REQUIREMENTS:
This algorithm will request following information during its operation:

* function vector f[] at given point X
* function vector f[] and Jacobian of f[] (simultaneously) at given point

There are several overloaded versions of  MinLMOptimize()  function  which
correspond  to  different LM-like optimization algorithms provided by this
unit. You should choose version which accepts fvec()  and jac() callbacks.
First  one  is used to calculate f[] at given point, second one calculates
f[] and Jacobian df[i]/dx[j].

You can try to initialize MinLMState structure with VJ  function and  then
use incorrect version  of  MinLMOptimize()  (for  example,  version  which
works  with  general  form function and does not provide Jacobian), but it
will  lead  to  exception  being  thrown  after first attempt to calculate
Jacobian.


USAGE:
1. User initializes algorithm state with MinLMCreateVJ() call
2. User tunes solver parameters with MinLMSetCond(),  MinLMSetStpMax() and
   other functions
3. User calls MinLMOptimize() function which  takes algorithm  state   and
   callback functions.
4. User calls MinLMResults() to get solution
5. Optionally, user may call MinLMRestartFrom() to solve  another  problem
   with same N/M but another starting point and/or another function.
   MinLMRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   dimension, N>1
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    M       -   number of functions f[i]
    X       -   initial solution, array[0..N-1]

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTES:
1. you may tune stopping conditions with MinLMSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLMSetStpMax() function to bound algorithm's steps.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatevj(const ae_int_t m, const real_1d_array &x, minlmstate &state)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;

    n = x.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmcreatevj(n, m, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minlmstate*>(state.c_ptr()), &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
                IMPROVED LEVENBERG-MARQUARDT METHOD FOR
                 NON-LINEAR LEAST SQUARES OPTIMIZATION

DESCRIPTION:
This function is used to find minimum of function which is represented  as
sum of squares:
    F(x) = f[0]^2(x[0],...,x[n-1]) + ... + f[m-1]^2(x[0],...,x[n-1])
using value of function vector f[] only. Finite differences  are  used  to
calculate Jacobian.


REQUIREMENTS:
This algorithm will request following information during its operation:
* function vector f[] at given point X

There are several overloaded versions of  MinLMOptimize()  function  which
correspond  to  different LM-like optimization algorithms provided by this
unit. You should choose version which accepts fvec() callback.

You can try to initialize MinLMState structure with VJ  function and  then
use incorrect version  of  MinLMOptimize()  (for  example,  version  which
works with general form function and does not accept function vector), but
it will  lead  to  exception being thrown after first attempt to calculate
Jacobian.


USAGE:
1. User initializes algorithm state with MinLMCreateV() call
2. User tunes solver parameters with MinLMSetCond(),  MinLMSetStpMax() and
   other functions
3. User calls MinLMOptimize() function which  takes algorithm  state   and
   callback functions.
4. User calls MinLMResults() to get solution
5. Optionally, user may call MinLMRestartFrom() to solve  another  problem
   with same N/M but another starting point and/or another function.
   MinLMRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   dimension, N>1
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    M       -   number of functions f[i]
    X       -   initial solution, array[0..N-1]
    DiffStep-   differentiation step, >0

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

See also MinLMIteration, MinLMResults.

NOTES:
1. you may tune stopping conditions with MinLMSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLMSetStpMax() function to bound algorithm's steps.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatev(const ae_int_t n, const ae_int_t m, const real_1d_array &x, const double diffstep, minlmstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmcreatev(n, m, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), diffstep, const_cast<alglib_impl::minlmstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
                IMPROVED LEVENBERG-MARQUARDT METHOD FOR
                 NON-LINEAR LEAST SQUARES OPTIMIZATION

DESCRIPTION:
This function is used to find minimum of function which is represented  as
sum of squares:
    F(x) = f[0]^2(x[0],...,x[n-1]) + ... + f[m-1]^2(x[0],...,x[n-1])
using value of function vector f[] only. Finite differences  are  used  to
calculate Jacobian.


REQUIREMENTS:
This algorithm will request following information during its operation:
* function vector f[] at given point X

There are several overloaded versions of  MinLMOptimize()  function  which
correspond  to  different LM-like optimization algorithms provided by this
unit. You should choose version which accepts fvec() callback.

You can try to initialize MinLMState structure with VJ  function and  then
use incorrect version  of  MinLMOptimize()  (for  example,  version  which
works with general form function and does not accept function vector), but
it will  lead  to  exception being thrown after first attempt to calculate
Jacobian.


USAGE:
1. User initializes algorithm state with MinLMCreateV() call
2. User tunes solver parameters with MinLMSetCond(),  MinLMSetStpMax() and
   other functions
3. User calls MinLMOptimize() function which  takes algorithm  state   and
   callback functions.
4. User calls MinLMResults() to get solution
5. Optionally, user may call MinLMRestartFrom() to solve  another  problem
   with same N/M but another starting point and/or another function.
   MinLMRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   dimension, N>1
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    M       -   number of functions f[i]
    X       -   initial solution, array[0..N-1]
    DiffStep-   differentiation step, >0

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

See also MinLMIteration, MinLMResults.

NOTES:
1. you may tune stopping conditions with MinLMSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLMSetStpMax() function to bound algorithm's steps.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatev(const ae_int_t m, const real_1d_array &x, const double diffstep, minlmstate &state)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;

    n = x.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmcreatev(n, m, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), diffstep, const_cast<alglib_impl::minlmstate*>(state.c_ptr()), &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
    LEVENBERG-MARQUARDT-LIKE METHOD FOR NON-LINEAR OPTIMIZATION

DESCRIPTION:
This  function  is  used  to  find  minimum  of general form (not "sum-of-
-squares") function
    F = F(x[0], ..., x[n-1])
using  its  gradient  and  Hessian.  Levenberg-Marquardt modification with
L-BFGS pre-optimization and internal pre-conditioned  L-BFGS  optimization
after each Levenberg-Marquardt step is used.


REQUIREMENTS:
This algorithm will request following information during its operation:

* function value F at given point X
* F and gradient G (simultaneously) at given point X
* F, G and Hessian H (simultaneously) at given point X

There are several overloaded versions of  MinLMOptimize()  function  which
correspond  to  different LM-like optimization algorithms provided by this
unit. You should choose version which accepts func(),  grad()  and  hess()
function pointers. First pointer is used to calculate F  at  given  point,
second  one  calculates  F(x)  and  grad F(x),  third one calculates F(x),
grad F(x), hess F(x).

You can try to initialize MinLMState structure with FGH-function and  then
use incorrect version of MinLMOptimize() (for example, version which  does
not provide Hessian matrix), but it will lead to  exception  being  thrown
after first attempt to calculate Hessian.


USAGE:
1. User initializes algorithm state with MinLMCreateFGH() call
2. User tunes solver parameters with MinLMSetCond(),  MinLMSetStpMax() and
   other functions
3. User calls MinLMOptimize() function which  takes algorithm  state   and
   pointers (delegates, etc.) to callback functions.
4. User calls MinLMResults() to get solution
5. Optionally, user may call MinLMRestartFrom() to solve  another  problem
   with same N but another starting point and/or another function.
   MinLMRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   dimension, N>1
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    X       -   initial solution, array[0..N-1]

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTES:
1. you may tune stopping conditions with MinLMSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLMSetStpMax() function to bound algorithm's steps.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatefgh(const ae_int_t n, const real_1d_array &x, minlmstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmcreatefgh(n, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minlmstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
    LEVENBERG-MARQUARDT-LIKE METHOD FOR NON-LINEAR OPTIMIZATION

DESCRIPTION:
This  function  is  used  to  find  minimum  of general form (not "sum-of-
-squares") function
    F = F(x[0], ..., x[n-1])
using  its  gradient  and  Hessian.  Levenberg-Marquardt modification with
L-BFGS pre-optimization and internal pre-conditioned  L-BFGS  optimization
after each Levenberg-Marquardt step is used.


REQUIREMENTS:
This algorithm will request following information during its operation:

* function value F at given point X
* F and gradient G (simultaneously) at given point X
* F, G and Hessian H (simultaneously) at given point X

There are several overloaded versions of  MinLMOptimize()  function  which
correspond  to  different LM-like optimization algorithms provided by this
unit. You should choose version which accepts func(),  grad()  and  hess()
function pointers. First pointer is used to calculate F  at  given  point,
second  one  calculates  F(x)  and  grad F(x),  third one calculates F(x),
grad F(x), hess F(x).

You can try to initialize MinLMState structure with FGH-function and  then
use incorrect version of MinLMOptimize() (for example, version which  does
not provide Hessian matrix), but it will lead to  exception  being  thrown
after first attempt to calculate Hessian.


USAGE:
1. User initializes algorithm state with MinLMCreateFGH() call
2. User tunes solver parameters with MinLMSetCond(),  MinLMSetStpMax() and
   other functions
3. User calls MinLMOptimize() function which  takes algorithm  state   and
   pointers (delegates, etc.) to callback functions.
4. User calls MinLMResults() to get solution
5. Optionally, user may call MinLMRestartFrom() to solve  another  problem
   with same N but another starting point and/or another function.
   MinLMRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   dimension, N>1
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    X       -   initial solution, array[0..N-1]

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTES:
1. you may tune stopping conditions with MinLMSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLMSetStpMax() function to bound algorithm's steps.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatefgh(const real_1d_array &x, minlmstate &state)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;

    n = x.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmcreatefgh(n, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minlmstate*>(state.c_ptr()), &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets stopping conditions for Levenberg-Marquardt optimization
algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    EpsG    -   >=0
                The  subroutine  finishes  its  work   if   the  condition
                |v|<EpsG is satisfied, where:
                * |.| means Euclidian norm
                * v - scaled gradient vector, v[i]=g[i]*s[i]
                * g - gradient
                * s - scaling coefficients set by MinLMSetScale()
    EpsF    -   >=0
                The  subroutine  finishes  its work if on k+1-th iteration
                the  condition  |F(k+1)-F(k)|<=EpsF*max{|F(k)|,|F(k+1)|,1}
                is satisfied.
    EpsX    -   >=0
                The subroutine finishes its work if  on  k+1-th  iteration
                the condition |v|<=EpsX is fulfilled, where:
                * |.| means Euclidian norm
                * v - scaled step vector, v[i]=dx[i]/s[i]
                * dx - ste pvector, dx=X(k+1)-X(k)
                * s - scaling coefficients set by MinLMSetScale()
    MaxIts  -   maximum number of iterations. If MaxIts=0, the  number  of
                iterations   is    unlimited.   Only   Levenberg-Marquardt
                iterations  are  counted  (L-BFGS/CG  iterations  are  NOT
                counted because their cost is very low compared to that of
                LM).

Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to
automatic stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlmsetcond(const minlmstate &state, const double epsg, const double epsf, const double epsx, const ae_int_t maxits)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmsetcond(const_cast<alglib_impl::minlmstate*>(state.c_ptr()), epsg, epsf, epsx, maxits, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    NeedXRep-   whether iteration reports are needed or not

If NeedXRep is True, algorithm will call rep() callback function if  it is
provided to MinLMOptimize(). Both Levenberg-Marquardt and internal  L-BFGS
iterations are reported.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlmsetxrep(const minlmstate &state, const bool needxrep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmsetxrep(const_cast<alglib_impl::minlmstate*>(state.c_ptr()), needxrep, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets maximum step length

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    StpMax  -   maximum step length, >=0. Set StpMax to 0.0,  if you don't
                want to limit step length.

Use this subroutine when you optimize target function which contains exp()
or  other  fast  growing  functions,  and optimization algorithm makes too
large  steps  which  leads  to overflow. This function allows us to reject
steps  that  are  too  large  (and  therefore  expose  us  to the possible
overflow) without actually calculating function value at the x+stp*d.

NOTE: non-zero StpMax leads to moderate  performance  degradation  because
intermediate  step  of  preconditioned L-BFGS optimization is incompatible
with limits on step size.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlmsetstpmax(const minlmstate &state, const double stpmax)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmsetstpmax(const_cast<alglib_impl::minlmstate*>(state.c_ptr()), stpmax, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets scaling coefficients for LM optimizer.

ALGLIB optimizers use scaling matrices to test stopping  conditions  (step
size and gradient are scaled before comparison with tolerances).  Scale of
the I-th variable is a translation invariant measure of:
a) "how large" the variable is
b) how large the step should be to make significant changes in the function

Generally, scale is NOT considered to be a form of preconditioner.  But LM
optimizer is unique in that it uses scaling matrix both  in  the  stopping
condition tests and as Marquardt damping factor.

Proper scaling is very important for the algorithm performance. It is less
important for the quality of results, but still has some influence (it  is
easier  to  converge  when  variables  are  properly  scaled, so premature
stopping is possible when very badly scalled variables are  combined  with
relaxed stopping conditions).

INPUT PARAMETERS:
    State   -   structure stores algorithm state
    S       -   array[N], non-zero scaling coefficients
                S[i] may be negative, sign doesn't matter.

  -- ALGLIB --
     Copyright 14.01.2011 by Bochkanov Sergey
*************************************************************************/
void minlmsetscale(const minlmstate &state, const real_1d_array &s)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmsetscale(const_cast<alglib_impl::minlmstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(s.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function sets boundary constraints for LM optimizer

Boundary constraints are inactive by default (after initial creation).
They are preserved until explicitly turned off with another SetBC() call.

INPUT PARAMETERS:
    State   -   structure stores algorithm state
    BndL    -   lower bounds, array[N].
                If some (all) variables are unbounded, you may specify
                very small number or -INF (latter is recommended because
                it will allow solver to use better algorithm).
    BndU    -   upper bounds, array[N].
                If some (all) variables are unbounded, you may specify
                very large number or +INF (latter is recommended because
                it will allow solver to use better algorithm).

NOTE 1: it is possible to specify BndL[i]=BndU[i]. In this case I-th
variable will be "frozen" at X[i]=BndL[i]=BndU[i].

NOTE 2: this solver has following useful properties:
* bound constraints are always satisfied exactly
* function is evaluated only INSIDE area specified by bound constraints
  or at its boundary

  -- ALGLIB --
     Copyright 14.01.2011 by Bochkanov Sergey
*************************************************************************/
void minlmsetbc(const minlmstate &state, const real_1d_array &bndl, const real_1d_array &bndu)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmsetbc(const_cast<alglib_impl::minlmstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(bndl.c_ptr()), const_cast<alglib_impl::ae_vector*>(bndu.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function is used to change acceleration settings

You can choose between three acceleration strategies:
* AccType=0, no acceleration.
* AccType=1, secant updates are used to update quadratic model after  each
  iteration. After fixed number of iterations (or after  model  breakdown)
  we  recalculate  quadratic  model  using  analytic  Jacobian  or  finite
  differences. Number of secant-based iterations depends  on  optimization
  settings: about 3 iterations - when we have analytic Jacobian, up to 2*N
  iterations - when we use finite differences to calculate Jacobian.

AccType=1 is recommended when Jacobian  calculation  cost  is  prohibitive
high (several Mx1 function vector calculations  followed  by  several  NxN
Cholesky factorizations are faster than calculation of one M*N  Jacobian).
It should also be used when we have no Jacobian, because finite difference
approximation takes too much time to compute.

Table below list  optimization  protocols  (XYZ  protocol  corresponds  to
MinLMCreateXYZ) and acceleration types they support (and use by  default).

ACCELERATION TYPES SUPPORTED BY OPTIMIZATION PROTOCOLS:

protocol    0   1   comment
V           +   +
VJ          +   +
FGH         +

DAFAULT VALUES:

protocol    0   1   comment
V               x   without acceleration it is so slooooooooow
VJ          x
FGH         x

NOTE: this  function should be called before optimization. Attempt to call
it during algorithm iterations may result in unexpected behavior.

NOTE: attempt to call this function with unsupported protocol/acceleration
combination will result in exception being thrown.

  -- ALGLIB --
     Copyright 14.10.2010 by Bochkanov Sergey
*************************************************************************/
void minlmsetacctype(const minlmstate &state, const ae_int_t acctype)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmsetacctype(const_cast<alglib_impl::minlmstate*>(state.c_ptr()), acctype, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function provides reverse communication interface
Reverse communication interface is not documented or recommended to use.
See below for functions which provide better documented API
*************************************************************************/
bool minlmiteration(const minlmstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        ae_bool result = alglib_impl::minlmiteration(const_cast<alglib_impl::minlmstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<bool*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}


void minlmoptimize(minlmstate &state,
    void (*fvec)(const real_1d_array &x, real_1d_array &fi, void *ptr),
    void  (*rep)(const real_1d_array &x, double func, void *ptr), 
    void *ptr)
{
    alglib_impl::ae_state _alglib_env_state;
    if( fvec==NULL )
        throw ap_error("ALGLIB: error in 'minlmoptimize()' (fvec is NULL)");
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        while( alglib_impl::minlmiteration(state.c_ptr(), &_alglib_env_state) )
        {
            if( state.needfi )
            {
                fvec(state.x, state.fi, ptr);
                continue;
            }
            if( state.xupdated )
            {
                if( rep!=NULL )
                    rep(state.x, state.f, ptr);
                continue;
            }
            throw ap_error("ALGLIB: error in 'minlmoptimize' (some derivatives were not provided?)");
        }
        alglib_impl::ae_state_clear(&_alglib_env_state);
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}


void minlmoptimize(minlmstate &state,
    void (*fvec)(const real_1d_array &x, real_1d_array &fi, void *ptr),
    void  (*jac)(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr),
    void  (*rep)(const real_1d_array &x, double func, void *ptr), 
    void *ptr)
{
    alglib_impl::ae_state _alglib_env_state;
    if( fvec==NULL )
        throw ap_error("ALGLIB: error in 'minlmoptimize()' (fvec is NULL)");
    if( jac==NULL )
        throw ap_error("ALGLIB: error in 'minlmoptimize()' (jac is NULL)");
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        while( alglib_impl::minlmiteration(state.c_ptr(), &_alglib_env_state) )
        {
            if( state.needfi )
            {
                fvec(state.x, state.fi, ptr);
                continue;
            }
            if( state.needfij )
            {
                jac(state.x, state.fi, state.j, ptr);
                continue;
            }
            if( state.xupdated )
            {
                if( rep!=NULL )
                    rep(state.x, state.f, ptr);
                continue;
            }
            throw ap_error("ALGLIB: error in 'minlmoptimize' (some derivatives were not provided?)");
        }
        alglib_impl::ae_state_clear(&_alglib_env_state);
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}


void minlmoptimize(minlmstate &state,
    void (*func)(const real_1d_array &x, double &func, void *ptr),
    void (*grad)(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr),
    void (*hess)(const real_1d_array &x, double &func, real_1d_array &grad, real_2d_array &hess, void *ptr),
    void  (*rep)(const real_1d_array &x, double func, void *ptr), 
    void *ptr)
{
    alglib_impl::ae_state _alglib_env_state;
    if( func==NULL )
        throw ap_error("ALGLIB: error in 'minlmoptimize()' (func is NULL)");
    if( grad==NULL )
        throw ap_error("ALGLIB: error in 'minlmoptimize()' (grad is NULL)");
    if( hess==NULL )
        throw ap_error("ALGLIB: error in 'minlmoptimize()' (hess is NULL)");
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        while( alglib_impl::minlmiteration(state.c_ptr(), &_alglib_env_state) )
        {
            if( state.needf )
            {
                func(state.x, state.f, ptr);
                continue;
            }
            if( state.needfg )
            {
                grad(state.x, state.f, state.g, ptr);
                continue;
            }
            if( state.needfgh )
            {
                hess(state.x, state.f, state.g, state.h, ptr);
                continue;
            }
            if( state.xupdated )
            {
                if( rep!=NULL )
                    rep(state.x, state.f, ptr);
                continue;
            }
            throw ap_error("ALGLIB: error in 'minlmoptimize' (some derivatives were not provided?)");
        }
        alglib_impl::ae_state_clear(&_alglib_env_state);
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}


void minlmoptimize(minlmstate &state,
    void (*func)(const real_1d_array &x, double &func, void *ptr),
    void  (*jac)(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr),
    void  (*rep)(const real_1d_array &x, double func, void *ptr), 
    void *ptr)
{
    alglib_impl::ae_state _alglib_env_state;
    if( func==NULL )
        throw ap_error("ALGLIB: error in 'minlmoptimize()' (func is NULL)");
    if( jac==NULL )
        throw ap_error("ALGLIB: error in 'minlmoptimize()' (jac is NULL)");
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        while( alglib_impl::minlmiteration(state.c_ptr(), &_alglib_env_state) )
        {
            if( state.needf )
            {
                func(state.x, state.f, ptr);
                continue;
            }
            if( state.needfij )
            {
                jac(state.x, state.fi, state.j, ptr);
                continue;
            }
            if( state.xupdated )
            {
                if( rep!=NULL )
                    rep(state.x, state.f, ptr);
                continue;
            }
            throw ap_error("ALGLIB: error in 'minlmoptimize' (some derivatives were not provided?)");
        }
        alglib_impl::ae_state_clear(&_alglib_env_state);
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}


void minlmoptimize(minlmstate &state,
    void (*func)(const real_1d_array &x, double &func, void *ptr),
    void (*grad)(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr),
    void  (*jac)(const real_1d_array &x, real_1d_array &fi, real_2d_array &jac, void *ptr),
    void  (*rep)(const real_1d_array &x, double func, void *ptr), 
    void *ptr)
{
    alglib_impl::ae_state _alglib_env_state;
    if( func==NULL )
        throw ap_error("ALGLIB: error in 'minlmoptimize()' (func is NULL)");
    if( grad==NULL )
        throw ap_error("ALGLIB: error in 'minlmoptimize()' (grad is NULL)");
    if( jac==NULL )
        throw ap_error("ALGLIB: error in 'minlmoptimize()' (jac is NULL)");
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        while( alglib_impl::minlmiteration(state.c_ptr(), &_alglib_env_state) )
        {
            if( state.needf )
            {
                func(state.x, state.f, ptr);
                continue;
            }
            if( state.needfg )
            {
                grad(state.x, state.f, state.g, ptr);
                continue;
            }
            if( state.needfij )
            {
                jac(state.x, state.fi, state.j, ptr);
                continue;
            }
            if( state.xupdated )
            {
                if( rep!=NULL )
                    rep(state.x, state.f, ptr);
                continue;
            }
            throw ap_error("ALGLIB: error in 'minlmoptimize' (some derivatives were not provided?)");
        }
        alglib_impl::ae_state_clear(&_alglib_env_state);
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}



/*************************************************************************
Levenberg-Marquardt algorithm results

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    X       -   array[0..N-1], solution
    Rep     -   optimization report;
                see comments for this structure for more info.

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmresults(const minlmstate &state, real_1d_array &x, minlmreport &rep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmresults(const_cast<alglib_impl::minlmstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minlmreport*>(rep.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Levenberg-Marquardt algorithm results

Buffered implementation of MinLMResults(), which uses pre-allocated buffer
to store X[]. If buffer size is  too  small,  it  resizes  buffer.  It  is
intended to be used in the inner cycles of performance critical algorithms
where array reallocation penalty is too large to be ignored.

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmresultsbuf(const minlmstate &state, real_1d_array &x, minlmreport &rep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmresultsbuf(const_cast<alglib_impl::minlmstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minlmreport*>(rep.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This  subroutine  restarts  LM  algorithm from new point. All optimization
parameters are left unchanged.

This  function  allows  to  solve multiple  optimization  problems  (which
must have same number of dimensions) without object reallocation penalty.

INPUT PARAMETERS:
    State   -   structure used for reverse communication previously
                allocated with MinLMCreateXXX call.
    X       -   new starting point.

  -- ALGLIB --
     Copyright 30.07.2010 by Bochkanov Sergey
*************************************************************************/
void minlmrestartfrom(const minlmstate &state, const real_1d_array &x)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmrestartfrom(const_cast<alglib_impl::minlmstate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This is obsolete function.

Since ALGLIB 3.3 it is equivalent to MinLMCreateVJ().

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatevgj(const ae_int_t n, const ae_int_t m, const real_1d_array &x, minlmstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmcreatevgj(n, m, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minlmstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This is obsolete function.

Since ALGLIB 3.3 it is equivalent to MinLMCreateVJ().

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatevgj(const ae_int_t m, const real_1d_array &x, minlmstate &state)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;

    n = x.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmcreatevgj(n, m, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minlmstate*>(state.c_ptr()), &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This is obsolete function.

Since ALGLIB 3.3 it is equivalent to MinLMCreateFJ().

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatefgj(const ae_int_t n, const ae_int_t m, const real_1d_array &x, minlmstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmcreatefgj(n, m, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minlmstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This is obsolete function.

Since ALGLIB 3.3 it is equivalent to MinLMCreateFJ().

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatefgj(const ae_int_t m, const real_1d_array &x, minlmstate &state)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;

    n = x.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmcreatefgj(n, m, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minlmstate*>(state.c_ptr()), &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function is considered obsolete since ALGLIB 3.1.0 and is present for
backward  compatibility  only.  We  recommend  to use MinLMCreateVJ, which
provides similar, but more consistent and feature-rich interface.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatefj(const ae_int_t n, const ae_int_t m, const real_1d_array &x, minlmstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmcreatefj(n, m, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minlmstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function is considered obsolete since ALGLIB 3.1.0 and is present for
backward  compatibility  only.  We  recommend  to use MinLMCreateVJ, which
provides similar, but more consistent and feature-rich interface.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatefj(const ae_int_t m, const real_1d_array &x, minlmstate &state)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;

    n = x.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmcreatefj(n, m, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minlmstate*>(state.c_ptr()), &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This  subroutine  turns  on  verification  of  the  user-supplied analytic
gradient:
* user calls this subroutine before optimization begins
* MinLMOptimize() is called
* prior to actual optimization, for  each  function Fi and each  component
  of parameters  being  optimized X[j] algorithm performs following steps:
  * two trial steps are made to X[j]-TestStep*S[j] and X[j]+TestStep*S[j],
    where X[j] is j-th parameter and S[j] is a scale of j-th parameter
  * if needed, steps are bounded with respect to constraints on X[]
  * Fi(X) is evaluated at these trial points
  * we perform one more evaluation in the middle point of the interval
  * we  build  cubic  model using function values and derivatives at trial
    points and we compare its prediction with actual value in  the  middle
    point
  * in case difference between prediction and actual value is higher  than
    some predetermined threshold, algorithm stops with completion code -7;
    Rep.VarIdx is set to index of the parameter with incorrect derivative,
    Rep.FuncIdx is set to index of the function.
* after verification is over, algorithm proceeds to the actual optimization.

NOTE 1: verification  needs  N (parameters count) Jacobian evaluations. It
        is  very  costly  and  you  should use it only for low dimensional
        problems,  when  you  want  to  be  sure  that  you've   correctly
        calculated  analytic  derivatives.  You should not  use  it in the
        production code  (unless  you  want  to check derivatives provided
        by some third party).

NOTE 2: you  should  carefully  choose  TestStep. Value which is too large
        (so large that function behaviour is significantly non-cubic) will
        lead to false alarms. You may use  different  step  for  different
        parameters by means of setting scale with MinLMSetScale().

NOTE 3: this function may lead to false positives. In case it reports that
        I-th  derivative was calculated incorrectly, you may decrease test
        step  and  try  one  more  time  - maybe your function changes too
        sharply  and  your  step  is  too  large for such rapidly chanding
        function.

INPUT PARAMETERS:
    State       -   structure used to store algorithm state
    TestStep    -   verification step:
                    * TestStep=0 turns verification off
                    * TestStep>0 activates verification

  -- ALGLIB --
     Copyright 15.06.2012 by Bochkanov Sergey
*************************************************************************/
void minlmsetgradientcheck(const minlmstate &state, const double teststep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlmsetgradientcheck(const_cast<alglib_impl::minlmstate*>(state.c_ptr()), teststep, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************

*************************************************************************/
_minasastate_owner::_minasastate_owner()
{
    p_struct = (alglib_impl::minasastate*)alglib_impl::ae_malloc(sizeof(alglib_impl::minasastate), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_minasastate_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_minasastate_owner::_minasastate_owner(const _minasastate_owner &rhs)
{
    p_struct = (alglib_impl::minasastate*)alglib_impl::ae_malloc(sizeof(alglib_impl::minasastate), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_minasastate_init_copy(p_struct, const_cast<alglib_impl::minasastate*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_minasastate_owner& _minasastate_owner::operator=(const _minasastate_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_minasastate_clear(p_struct);
    if( !alglib_impl::_minasastate_init_copy(p_struct, const_cast<alglib_impl::minasastate*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_minasastate_owner::~_minasastate_owner()
{
    alglib_impl::_minasastate_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::minasastate* _minasastate_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::minasastate* _minasastate_owner::c_ptr() const
{
    return const_cast<alglib_impl::minasastate*>(p_struct);
}
minasastate::minasastate() : _minasastate_owner() ,needfg(p_struct->needfg),xupdated(p_struct->xupdated),f(p_struct->f),g(&p_struct->g),x(&p_struct->x)
{
}

minasastate::minasastate(const minasastate &rhs):_minasastate_owner(rhs) ,needfg(p_struct->needfg),xupdated(p_struct->xupdated),f(p_struct->f),g(&p_struct->g),x(&p_struct->x)
{
}

minasastate& minasastate::operator=(const minasastate &rhs)
{
    if( this==&rhs )
        return *this;
    _minasastate_owner::operator=(rhs);
    return *this;
}

minasastate::~minasastate()
{
}


/*************************************************************************

*************************************************************************/
_minasareport_owner::_minasareport_owner()
{
    p_struct = (alglib_impl::minasareport*)alglib_impl::ae_malloc(sizeof(alglib_impl::minasareport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_minasareport_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_minasareport_owner::_minasareport_owner(const _minasareport_owner &rhs)
{
    p_struct = (alglib_impl::minasareport*)alglib_impl::ae_malloc(sizeof(alglib_impl::minasareport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_minasareport_init_copy(p_struct, const_cast<alglib_impl::minasareport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_minasareport_owner& _minasareport_owner::operator=(const _minasareport_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_minasareport_clear(p_struct);
    if( !alglib_impl::_minasareport_init_copy(p_struct, const_cast<alglib_impl::minasareport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_minasareport_owner::~_minasareport_owner()
{
    alglib_impl::_minasareport_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::minasareport* _minasareport_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::minasareport* _minasareport_owner::c_ptr() const
{
    return const_cast<alglib_impl::minasareport*>(p_struct);
}
minasareport::minasareport() : _minasareport_owner() ,iterationscount(p_struct->iterationscount),nfev(p_struct->nfev),terminationtype(p_struct->terminationtype),activeconstraints(p_struct->activeconstraints)
{
}

minasareport::minasareport(const minasareport &rhs):_minasareport_owner(rhs) ,iterationscount(p_struct->iterationscount),nfev(p_struct->nfev),terminationtype(p_struct->terminationtype),activeconstraints(p_struct->activeconstraints)
{
}

minasareport& minasareport::operator=(const minasareport &rhs)
{
    if( this==&rhs )
        return *this;
    _minasareport_owner::operator=(rhs);
    return *this;
}

minasareport::~minasareport()
{
}

/*************************************************************************
Obsolete function, use MinLBFGSSetPrecDefault() instead.

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetdefaultpreconditioner(const minlbfgsstate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlbfgssetdefaultpreconditioner(const_cast<alglib_impl::minlbfgsstate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Obsolete function, use MinLBFGSSetCholeskyPreconditioner() instead.

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetcholeskypreconditioner(const minlbfgsstate &state, const real_2d_array &p, const bool isupper)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minlbfgssetcholeskypreconditioner(const_cast<alglib_impl::minlbfgsstate*>(state.c_ptr()), const_cast<alglib_impl::ae_matrix*>(p.c_ptr()), isupper, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This is obsolete function which was used by previous version of the  BLEIC
optimizer. It does nothing in the current version of BLEIC.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetbarrierwidth(const minbleicstate &state, const double mu)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleicsetbarrierwidth(const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), mu, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This is obsolete function which was used by previous version of the  BLEIC
optimizer. It does nothing in the current version of BLEIC.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetbarrierdecay(const minbleicstate &state, const double mudecay)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minbleicsetbarrierdecay(const_cast<alglib_impl::minbleicstate*>(state.c_ptr()), mudecay, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Obsolete optimization algorithm.
Was replaced by MinBLEIC subpackage.

  -- ALGLIB --
     Copyright 25.03.2010 by Bochkanov Sergey
*************************************************************************/
void minasacreate(const ae_int_t n, const real_1d_array &x, const real_1d_array &bndl, const real_1d_array &bndu, minasastate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minasacreate(n, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::ae_vector*>(bndl.c_ptr()), const_cast<alglib_impl::ae_vector*>(bndu.c_ptr()), const_cast<alglib_impl::minasastate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Obsolete optimization algorithm.
Was replaced by MinBLEIC subpackage.

  -- ALGLIB --
     Copyright 25.03.2010 by Bochkanov Sergey
*************************************************************************/
void minasacreate(const real_1d_array &x, const real_1d_array &bndl, const real_1d_array &bndu, minasastate &state)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;
    if( (x.length()!=bndl.length()) || (x.length()!=bndu.length()))
        throw ap_error("Error while calling 'minasacreate': looks like one of arguments has wrong size");
    n = x.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minasacreate(n, const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::ae_vector*>(bndl.c_ptr()), const_cast<alglib_impl::ae_vector*>(bndu.c_ptr()), const_cast<alglib_impl::minasastate*>(state.c_ptr()), &_alglib_env_state);

        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Obsolete optimization algorithm.
Was replaced by MinBLEIC subpackage.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minasasetcond(const minasastate &state, const double epsg, const double epsf, const double epsx, const ae_int_t maxits)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minasasetcond(const_cast<alglib_impl::minasastate*>(state.c_ptr()), epsg, epsf, epsx, maxits, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Obsolete optimization algorithm.
Was replaced by MinBLEIC subpackage.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minasasetxrep(const minasastate &state, const bool needxrep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minasasetxrep(const_cast<alglib_impl::minasastate*>(state.c_ptr()), needxrep, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Obsolete optimization algorithm.
Was replaced by MinBLEIC subpackage.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minasasetalgorithm(const minasastate &state, const ae_int_t algotype)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minasasetalgorithm(const_cast<alglib_impl::minasastate*>(state.c_ptr()), algotype, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Obsolete optimization algorithm.
Was replaced by MinBLEIC subpackage.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minasasetstpmax(const minasastate &state, const double stpmax)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minasasetstpmax(const_cast<alglib_impl::minasastate*>(state.c_ptr()), stpmax, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This function provides reverse communication interface
Reverse communication interface is not documented or recommended to use.
See below for functions which provide better documented API
*************************************************************************/
bool minasaiteration(const minasastate &state)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        ae_bool result = alglib_impl::minasaiteration(const_cast<alglib_impl::minasastate*>(state.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<bool*>(&result));
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}


void minasaoptimize(minasastate &state,
    void (*grad)(const real_1d_array &x, double &func, real_1d_array &grad, void *ptr),
    void  (*rep)(const real_1d_array &x, double func, void *ptr), 
    void *ptr)
{
    alglib_impl::ae_state _alglib_env_state;
    if( grad==NULL )
        throw ap_error("ALGLIB: error in 'minasaoptimize()' (grad is NULL)");
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        while( alglib_impl::minasaiteration(state.c_ptr(), &_alglib_env_state) )
        {
            if( state.needfg )
            {
                grad(state.x, state.f, state.g, ptr);
                continue;
            }
            if( state.xupdated )
            {
                if( rep!=NULL )
                    rep(state.x, state.f, ptr);
                continue;
            }
            throw ap_error("ALGLIB: error in 'minasaoptimize' (some derivatives were not provided?)");
        }
        alglib_impl::ae_state_clear(&_alglib_env_state);
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}



/*************************************************************************
Obsolete optimization algorithm.
Was replaced by MinBLEIC subpackage.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void minasaresults(const minasastate &state, real_1d_array &x, minasareport &rep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minasaresults(const_cast<alglib_impl::minasastate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minasareport*>(rep.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Obsolete optimization algorithm.
Was replaced by MinBLEIC subpackage.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void minasaresultsbuf(const minasastate &state, real_1d_array &x, minasareport &rep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minasaresultsbuf(const_cast<alglib_impl::minasastate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::minasareport*>(rep.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Obsolete optimization algorithm.
Was replaced by MinBLEIC subpackage.

  -- ALGLIB --
     Copyright 30.07.2010 by Bochkanov Sergey
*************************************************************************/
void minasarestartfrom(const minasastate &state, const real_1d_array &x, const real_1d_array &bndl, const real_1d_array &bndu)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::minasarestartfrom(const_cast<alglib_impl::minasastate*>(state.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::ae_vector*>(bndl.c_ptr()), const_cast<alglib_impl::ae_vector*>(bndu.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return;
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(_alglib_env_state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}
}

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS IMPLEMENTATION OF COMPUTATIONAL CORE
//
/////////////////////////////////////////////////////////////////////////
namespace alglib_impl
{


static ae_int_t mincg_rscountdownlen = 10;
static double mincg_gtol = 0.3;
static void mincg_clearrequestfields(mincgstate* state, ae_state *_state);
static void mincg_preconditionedmultiply(mincgstate* state,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* work0,
     /* Real    */ ae_vector* work1,
     ae_state *_state);
static double mincg_preconditionedmultiply2(mincgstate* state,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* work0,
     /* Real    */ ae_vector* work1,
     ae_state *_state);
static void mincg_mincginitinternal(ae_int_t n,
     double diffstep,
     mincgstate* state,
     ae_state *_state);


//static double minbleic_svdtol = 100;
static double minbleic_maxouterits = 20;
static double minbleic_maxnonmonotoniclen = 1.0E5;
static void minbleic_clearrequestfields(minbleicstate* state,
     ae_state *_state);
static void minbleic_unscalepoint(minbleicstate* state,
     /* Real    */ ae_vector* xscaled,
     /* Real    */ ae_vector* xunscaled,
     ae_state *_state);
static void minbleic_projectpointandunscale(minbleicstate* state,
     /* Real    */ ae_vector* xscaled,
     /* Real    */ ae_vector* xunscaled,
     /* Real    */ ae_vector* rscaled,
     double* rnorm2,
     ae_state *_state);
static void minbleic_scalegradientandexpand(minbleicstate* state,
     /* Real    */ ae_vector* gunscaled,
     /* Real    */ ae_vector* gscaled,
     ae_state *_state);
static void minbleic_modifytargetfunction(minbleicstate* state,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* r,
     double rnorm2,
     double* f,
     /* Real    */ ae_vector* g,
     double* gnorm,
     double* mpgnorm,
     ae_state *_state);
static ae_bool minbleic_additionalcheckforconstraints(minbleicstate* state,
     /* Real    */ ae_vector* x,
     ae_state *_state);
static void minbleic_rebuildcexe(minbleicstate* state, ae_state *_state);
static void minbleic_makegradientprojection(minbleicstate* state,
     /* Real    */ ae_vector* pg,
     ae_state *_state);
static ae_bool minbleic_prepareconstraintmatrix(minbleicstate* state,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* g,
     /* Real    */ ae_vector* px,
     /* Real    */ ae_vector* pg,
     ae_state *_state);
static void minbleic_minbleicinitinternal(ae_int_t n,
     /* Real    */ ae_vector* x,
     double diffstep,
     minbleicstate* state,
     ae_state *_state);


static double minlbfgs_gtol = 0.4;
static void minlbfgs_clearrequestfields(minlbfgsstate* state,
     ae_state *_state);


static ae_int_t cqmodels_newtonrefinementits = 3;
static ae_bool cqmodels_cqmrebuild(convexquadraticmodel* s,
     ae_state *_state);
static void cqmodels_cqmsolveea(convexquadraticmodel* s,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* tmp,
     ae_state *_state);


static ae_int_t minqp_badqpitsbeforetermination = 5;
static ae_int_t minqp_maxlagrangeits = 10;
static ae_int_t minqp_minqpdropconstraints(minqpstate* state,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* buf,
     ae_state *_state);
static ae_int_t minqp_minqpboundedstepandactivation(minqpstate* state,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* xn,
     /* Real    */ ae_vector* buf,
     ae_state *_state);
static double minqp_minqpmodelvalue(convexquadraticmodel* a,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* xc,
     ae_int_t n,
     /* Real    */ ae_vector* tmp,
     ae_state *_state);
static double minqp_minqpgradientandconstraineddescent(minqpstate* state,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* gc,
     /* Real    */ ae_vector* d,
     ae_state *_state);
static ae_bool minqp_minqpconstrainedoptimum(convexquadraticmodel* a,
     double anorm,
     /* Real    */ ae_vector* b,
     /* Boolean */ ae_vector* activeb,
     /* Real    */ ae_vector* xc,
     ae_int_t n,
     /* Real    */ ae_matrix* c,
     /* Real    */ ae_vector* r,
     ae_int_t k,
     /* Real    */ ae_vector* xn,
     /* Real    */ ae_vector* tmp,
     /* Real    */ ae_vector* lagrangec,
     ae_state *_state);
static void minqp_orthogonalizeactiveconstraints(/* Real    */ ae_vector* xc,
     /* Boolean */ ae_vector* activeb,
     ae_int_t n,
     /* Real    */ ae_matrix* cleic,
     /* Boolean */ ae_vector* activelin,
     ae_int_t ctotal,
     /* Real    */ ae_matrix* activecm,
     /* Real    */ ae_vector* activecr,
     ae_int_t* nactive,
     ae_state *_state);


//static ae_int_t minlm_lmmodefj = 0;
//static ae_int_t minlm_lmmodefgj = 1;
//static ae_int_t minlm_lmmodefgh = 2;
//static ae_int_t minlm_lmflagnoprelbfgs = 1;
//static ae_int_t minlm_lmflagnointlbfgs = 2;
//static ae_int_t minlm_lmprelbfgsm = 5;
//static ae_int_t minlm_lmintlbfgsits = 5;
//static ae_int_t minlm_lbfgsnorealloc = 1;
static double minlm_lambdaup = 2.0;
static double minlm_lambdadown = 0.33;
static double minlm_suspiciousnu = 16;
static ae_int_t minlm_smallmodelage = 3;
static ae_int_t minlm_additers = 5;
static void minlm_lmprepare(ae_int_t n,
     ae_int_t m,
     ae_bool havegrad,
     minlmstate* state,
     ae_state *_state);
static void minlm_clearrequestfields(minlmstate* state, ae_state *_state);
static ae_bool minlm_increaselambda(double* lambdav,
     double* nu,
     ae_state *_state);
static void minlm_decreaselambda(double* lambdav,
     double* nu,
     ae_state *_state);
static double minlm_boundedscaledantigradnorm(minlmstate* state,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* g,
     ae_state *_state);


static ae_int_t mincomp_n1 = 2;
static ae_int_t mincomp_n2 = 2;
static double mincomp_stpmin = 1.0E-300;
static double mincomp_gtol = 0.3;
static double mincomp_gpaftol = 0.0001;
static double mincomp_gpadecay = 0.5;
static double mincomp_asarho = 0.5;
static double mincomp_asaboundedantigradnorm(minasastate* state,
     ae_state *_state);
static double mincomp_asaginorm(minasastate* state, ae_state *_state);
static double mincomp_asad1norm(minasastate* state, ae_state *_state);
static ae_bool mincomp_asauisempty(minasastate* state, ae_state *_state);
static void mincomp_clearrequestfields(minasastate* state,
     ae_state *_state);





/*************************************************************************
This subroutine is used to prepare threshold value which will be used for
trimming of the target function (see comments on TrimFunction() for more
information).

This function accepts only one parameter: function value at the starting
point. It returns threshold which will be used for trimming.

  -- ALGLIB --
     Copyright 10.05.2011 by Bochkanov Sergey
*************************************************************************/
void trimprepare(double f, double* threshold, ae_state *_state)
{

    *threshold = 0;

    *threshold = 10*(ae_fabs(f, _state)+1);
}


/*************************************************************************
This subroutine is used to "trim" target function, i.e. to do following
transformation:

                   { {F,G}          if F<Threshold
    {F_tr, G_tr} = {
                   { {Threshold, 0} if F>=Threshold
                   
Such transformation allows us to  solve  problems  with  singularities  by
redefining function in such way that it becomes bounded from above.

  -- ALGLIB --
     Copyright 10.05.2011 by Bochkanov Sergey
*************************************************************************/
void trimfunction(double* f,
     /* Real    */ ae_vector* g,
     ae_int_t n,
     double threshold,
     ae_state *) //_state)
{
    ae_int_t i;


    if( ae_fp_greater_eq(*f,threshold) )
    {
        *f = threshold;
        for(i=0; i<=n-1; i++)
        {
            g->ptr.p_double[i] = 0.0;
        }
    }
}


/*************************************************************************
This function enforces boundary constraints in the X.

This function correctly (although a bit inefficient) handles BL[i] which
are -INF and BU[i] which are +INF.

We have NMain+NSlack  dimensional  X,  with first NMain components bounded
by BL/BU, and next NSlack ones bounded by non-negativity constraints.

INPUT PARAMETERS
    X       -   array[NMain+NSlack], point
    BL      -   array[NMain], lower bounds
                (may contain -INF, when bound is not present)
    HaveBL  -   array[NMain], if HaveBL[i] is False,
                then i-th bound is not present
    BU      -   array[NMain], upper bounds
                (may contain +INF, when bound is not present)
    HaveBU  -   array[NMain], if HaveBU[i] is False,
                then i-th bound is not present

OUTPUT PARAMETERS
    X       -   X with all constraints being enforced

It returns True when constraints are consistent,
False - when constraints are inconsistent.

  -- ALGLIB --
     Copyright 10.01.2012 by Bochkanov Sergey
*************************************************************************/
ae_bool enforceboundaryconstraints(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* bl,
     /* Boolean */ ae_vector* havebl,
     /* Real    */ ae_vector* bu,
     /* Boolean */ ae_vector* havebu,
     ae_int_t nmain,
     ae_int_t nslack,
     ae_state *) //_state)
{
    ae_int_t i;
    ae_bool result;


    result = ae_false;
    for(i=0; i<=nmain-1; i++)
    {
        if( (havebl->ptr.p_bool[i]&&havebu->ptr.p_bool[i])&&ae_fp_greater(bl->ptr.p_double[i],bu->ptr.p_double[i]) )
        {
            return result;
        }
        if( havebl->ptr.p_bool[i]&&ae_fp_less(x->ptr.p_double[i],bl->ptr.p_double[i]) )
        {
            x->ptr.p_double[i] = bl->ptr.p_double[i];
        }
        if( havebu->ptr.p_bool[i]&&ae_fp_greater(x->ptr.p_double[i],bu->ptr.p_double[i]) )
        {
            x->ptr.p_double[i] = bu->ptr.p_double[i];
        }
    }
    for(i=0; i<=nslack-1; i++)
    {
        if( ae_fp_less(x->ptr.p_double[nmain+i],0) )
        {
            x->ptr.p_double[nmain+i] = 0;
        }
    }
    result = ae_true;
    return result;
}


/*************************************************************************
This function projects gradient into feasible area of boundary constrained
optimization  problem.  X  can  be  infeasible  with  respect  to boundary
constraints.  We  have  NMain+NSlack  dimensional  X,   with  first  NMain 
components bounded by BL/BU, and next NSlack ones bounded by non-negativity
constraints.

INPUT PARAMETERS
    X       -   array[NMain+NSlack], point
    G       -   array[NMain+NSlack], gradient
    BL      -   lower bounds (may contain -INF, when bound is not present)
    HaveBL  -   if HaveBL[i] is False, then i-th bound is not present
    BU      -   upper bounds (may contain +INF, when bound is not present)
    HaveBU  -   if HaveBU[i] is False, then i-th bound is not present

OUTPUT PARAMETERS
    G       -   projection of G. Components of G which satisfy one of the
                following
                    (1) (X[I]<=BndL[I]) and (G[I]>0), OR
                    (2) (X[I]>=BndU[I]) and (G[I]<0)
                are replaced by zeros.

NOTE 1: this function assumes that constraints are feasible. It throws
exception otherwise.

NOTE 2: in fact, projection of ANTI-gradient is calculated,  because  this
function trims components of -G which points outside of the feasible area.
However, working with -G is considered confusing, because all optimization
source work with G.

  -- ALGLIB --
     Copyright 10.01.2012 by Bochkanov Sergey
*************************************************************************/
void projectgradientintobc(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* g,
     /* Real    */ ae_vector* bl,
     /* Boolean */ ae_vector* havebl,
     /* Real    */ ae_vector* bu,
     /* Boolean */ ae_vector* havebu,
     ae_int_t nmain,
     ae_int_t nslack,
     ae_state *_state)
{
    ae_int_t i;


    for(i=0; i<=nmain-1; i++)
    {
        ae_assert((!havebl->ptr.p_bool[i]||!havebu->ptr.p_bool[i])||ae_fp_less_eq(bl->ptr.p_double[i],bu->ptr.p_double[i]), "ProjectGradientIntoBC: internal error (infeasible constraints)", _state);
        if( (havebl->ptr.p_bool[i]&&ae_fp_less_eq(x->ptr.p_double[i],bl->ptr.p_double[i]))&&ae_fp_greater(g->ptr.p_double[i],0) )
        {
            g->ptr.p_double[i] = 0;
        }
        if( (havebu->ptr.p_bool[i]&&ae_fp_greater_eq(x->ptr.p_double[i],bu->ptr.p_double[i]))&&ae_fp_less(g->ptr.p_double[i],0) )
        {
            g->ptr.p_double[i] = 0;
        }
    }
    for(i=0; i<=nslack-1; i++)
    {
        if( ae_fp_less_eq(x->ptr.p_double[nmain+i],0)&&ae_fp_greater(g->ptr.p_double[nmain+i],0) )
        {
            g->ptr.p_double[nmain+i] = 0;
        }
    }
}


/*************************************************************************
Given
    a) initial point X0[NMain+NSlack]
       (feasible with respect to bound constraints)
    b) step vector alpha*D[NMain+NSlack]
    c) boundary constraints BndL[NMain], BndU[NMain]
    d) implicit non-negativity constraints for slack variables
this  function  calculates  bound  on  the step length subject to boundary
constraints.

It returns:
    *  MaxStepLen - such step length that X0+MaxStepLen*alpha*D is exactly
       at the boundary given by constraints
    *  VariableToFreeze - index of the constraint to be activated,
       0 <= VariableToFreeze < NMain+NSlack
    *  ValueToFreeze - value of the corresponding constraint.

Notes:
    * it is possible that several constraints can be activated by the step
      at once. In such cases only one constraint is returned. It is caller
      responsibility to check other constraints. This function makes  sure
      that we activate at least one constraint, and everything else is the
      responsibility of the caller.
    * steps smaller than MaxStepLen still can activate constraints due  to
      numerical errors. Thus purpose of this  function  is  not  to  guard 
      against accidental activation of the constraints - quite the reverse, 
      its purpose is to activate at least constraint upon performing  step
      which is too long.
    * in case there is no constraints to activate, we return negative
      VariableToFreeze and zero MaxStepLen and ValueToFreeze.
    * this function assumes that constraints are consistent; it throws
      exception otherwise.

INPUT PARAMETERS
    X           -   array[NMain+NSlack], point. Must be feasible with respect 
                    to bound constraints (exception will be thrown otherwise)
    D           -   array[NMain+NSlack], step direction
    alpha       -   scalar multiplier before D, alpha<>0
    BndL        -   lower bounds, array[NMain]
                    (may contain -INF, when bound is not present)
    HaveBndL    -   array[NMain], if HaveBndL[i] is False,
                    then i-th bound is not present
    BndU        -   array[NMain], upper bounds
                    (may contain +INF, when bound is not present)
    HaveBndU    -   array[NMain], if HaveBndU[i] is False,
                    then i-th bound is not present
    NMain       -   number of main variables
    NSlack      -   number of slack variables
    
OUTPUT PARAMETERS
    VariableToFreeze:
                    * negative value     = step is unbounded, ValueToFreeze=0,
                                           MaxStepLen=0.
                    * non-negative value = at least one constraint, given by
                                           this parameter, will  be  activated
                                           upon performing maximum step.
    ValueToFreeze-  value of the variable which will be constrained
    MaxStepLen  -   maximum length of the step. Can be zero when step vector
                    looks outside of the feasible area.

  -- ALGLIB --
     Copyright 10.01.2012 by Bochkanov Sergey
*************************************************************************/
void calculatestepbound(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* d,
     double alpha,
     /* Real    */ ae_vector* bndl,
     /* Boolean */ ae_vector* havebndl,
     /* Real    */ ae_vector* bndu,
     /* Boolean */ ae_vector* havebndu,
     ae_int_t nmain,
     ae_int_t nslack,
     ae_int_t* variabletofreeze,
     double* valuetofreeze,
     double* maxsteplen,
     ae_state *_state)
{
    ae_int_t i;
    double prevmax;
    double initval;

    *variabletofreeze = 0;
    *valuetofreeze = 0;
    *maxsteplen = 0;

    ae_assert(ae_fp_neq(alpha,0), "CalculateStepBound: zero alpha", _state);
    *variabletofreeze = -1;
    initval = ae_maxrealnumber;
    *maxsteplen = initval;
    for(i=0; i<=nmain-1; i++)
    {
        if( havebndl->ptr.p_bool[i]&&ae_fp_less(alpha*d->ptr.p_double[i],0) )
        {
            ae_assert(ae_fp_greater_eq(x->ptr.p_double[i],bndl->ptr.p_double[i]), "CalculateStepBound: infeasible X", _state);
            prevmax = *maxsteplen;
            *maxsteplen = safeminposrv(x->ptr.p_double[i]-bndl->ptr.p_double[i], -alpha*d->ptr.p_double[i], *maxsteplen, _state);
            if( ae_fp_less(*maxsteplen,prevmax) )
            {
                *variabletofreeze = i;
                *valuetofreeze = bndl->ptr.p_double[i];
            }
        }
        if( havebndu->ptr.p_bool[i]&&ae_fp_greater(alpha*d->ptr.p_double[i],0) )
        {
            ae_assert(ae_fp_less_eq(x->ptr.p_double[i],bndu->ptr.p_double[i]), "CalculateStepBound: infeasible X", _state);
            prevmax = *maxsteplen;
            *maxsteplen = safeminposrv(bndu->ptr.p_double[i]-x->ptr.p_double[i], alpha*d->ptr.p_double[i], *maxsteplen, _state);
            if( ae_fp_less(*maxsteplen,prevmax) )
            {
                *variabletofreeze = i;
                *valuetofreeze = bndu->ptr.p_double[i];
            }
        }
    }
    for(i=0; i<=nslack-1; i++)
    {
        if( ae_fp_less(alpha*d->ptr.p_double[nmain+i],0) )
        {
            ae_assert(ae_fp_greater_eq(x->ptr.p_double[nmain+i],0), "CalculateStepBound: infeasible X", _state);
            prevmax = *maxsteplen;
            *maxsteplen = safeminposrv(x->ptr.p_double[nmain+i], -alpha*d->ptr.p_double[nmain+i], *maxsteplen, _state);
            if( ae_fp_less(*maxsteplen,prevmax) )
            {
                *variabletofreeze = nmain+i;
                *valuetofreeze = 0;
            }
        }
    }
    if( ae_fp_eq(*maxsteplen,initval) )
    {
        *valuetofreeze = 0;
        *maxsteplen = 0;
    }
}


/*************************************************************************
This function postprocesses bounded step by:
* analysing step length (whether it is equal to MaxStepLen) and activating 
  constraint given by VariableToFreeze if needed
* checking for additional bound constraints to activate

This function uses final point of the step, quantities calculated  by  the
CalculateStepBound()  function.  As  result,  it  returns  point  which is 
exactly feasible with respect to boundary constraints.

NOTE 1: this function does NOT handle and check linear equality constraints
NOTE 2: when StepTaken=MaxStepLen we always activate at least one constraint

INPUT PARAMETERS
    X           -   array[NMain+NSlack], final point to postprocess
    XPrev       -   array[NMain+NSlack], initial point
    BndL        -   lower bounds, array[NMain]
                    (may contain -INF, when bound is not present)
    HaveBndL    -   array[NMain], if HaveBndL[i] is False,
                    then i-th bound is not present
    BndU        -   array[NMain], upper bounds
                    (may contain +INF, when bound is not present)
    HaveBndU    -   array[NMain], if HaveBndU[i] is False,
                    then i-th bound is not present
    NMain       -   number of main variables
    NSlack      -   number of slack variables
    VariableToFreeze-result of CalculateStepBound()
    ValueToFreeze-  result of CalculateStepBound()
    StepTaken   -   actual step length (actual step is equal to the possibly 
                    non-unit step direction vector times this parameter).
                    StepTaken<=MaxStepLen.
    MaxStepLen  -   result of CalculateStepBound()
    
OUTPUT PARAMETERS
    X           -   point bounded with respect to constraints.
                    components corresponding to active constraints are exactly
                    equal to the boundary values.
                    
RESULT:
    number of constraints activated in addition to previously active ones.
    Constraints which were DEACTIVATED are ignored (do not influence
    function value).

  -- ALGLIB --
     Copyright 10.01.2012 by Bochkanov Sergey
*************************************************************************/
ae_int_t postprocessboundedstep(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* xprev,
     /* Real    */ ae_vector* bndl,
     /* Boolean */ ae_vector* havebndl,
     /* Real    */ ae_vector* bndu,
     /* Boolean */ ae_vector* havebndu,
     ae_int_t nmain,
     ae_int_t nslack,
     ae_int_t variabletofreeze,
     double valuetofreeze,
     double steptaken,
     double maxsteplen,
     ae_state *_state)
{
    ae_int_t i;
    ae_bool wasactivated;
    ae_int_t result;


    ae_assert(variabletofreeze<0||ae_fp_less_eq(steptaken,maxsteplen), "Assertion failed", _state);
    
    /*
     * Activate constraints
     */
    if( variabletofreeze>=0&&ae_fp_eq(steptaken,maxsteplen) )
    {
        x->ptr.p_double[variabletofreeze] = valuetofreeze;
    }
    for(i=0; i<=nmain-1; i++)
    {
        if( havebndl->ptr.p_bool[i]&&ae_fp_less(x->ptr.p_double[i],bndl->ptr.p_double[i]) )
        {
            x->ptr.p_double[i] = bndl->ptr.p_double[i];
        }
        if( havebndu->ptr.p_bool[i]&&ae_fp_greater(x->ptr.p_double[i],bndu->ptr.p_double[i]) )
        {
            x->ptr.p_double[i] = bndu->ptr.p_double[i];
        }
    }
    for(i=0; i<=nslack-1; i++)
    {
        if( ae_fp_less_eq(x->ptr.p_double[nmain+i],0) )
        {
            x->ptr.p_double[nmain+i] = 0;
        }
    }
    
    /*
     * Calculate number of constraints being activated
     */
    result = 0;
    for(i=0; i<=nmain-1; i++)
    {
        wasactivated = ae_fp_neq(x->ptr.p_double[i],xprev->ptr.p_double[i])&&((havebndl->ptr.p_bool[i]&&ae_fp_eq(x->ptr.p_double[i],bndl->ptr.p_double[i]))||(havebndu->ptr.p_bool[i]&&ae_fp_eq(x->ptr.p_double[i],bndu->ptr.p_double[i])));
        wasactivated = wasactivated||variabletofreeze==i;
        if( wasactivated )
        {
            result = result+1;
        }
    }
    for(i=0; i<=nslack-1; i++)
    {
        wasactivated = ae_fp_neq(x->ptr.p_double[nmain+i],xprev->ptr.p_double[nmain+i])&&ae_fp_eq(x->ptr.p_double[nmain+i],0.0);
        wasactivated = wasactivated||variabletofreeze==nmain+i;
        if( wasactivated )
        {
            result = result+1;
        }
    }
    return result;
}


/*************************************************************************
The  purpose  of  this  function is to prevent algorithm from "unsticking" 
from  the  active  bound  constraints  because  of  numerical noise in the
gradient or Hessian.

It is done by zeroing some components of the search direction D.  D[i]  is
zeroed when both (a) and (b) are true:
a) corresponding X[i] is exactly at the boundary
b) |D[i]*S[i]| <= DropTol*Sqrt(SUM(D[i]^2*S[I]^2))

D  can  be  step  direction , antigradient, gradient, or anything similar. 
Sign of D does not matter, nor matters step length.

NOTE 1: boundary constraints are expected to be consistent, as well as X
        is expected to be feasible. Exception will be thrown otherwise.

INPUT PARAMETERS
    D           -   array[NMain+NSlack], direction
    X           -   array[NMain+NSlack], current point
    BndL        -   lower bounds, array[NMain]
                    (may contain -INF, when bound is not present)
    HaveBndL    -   array[NMain], if HaveBndL[i] is False,
                    then i-th bound is not present
    BndU        -   array[NMain], upper bounds
                    (may contain +INF, when bound is not present)
    HaveBndU    -   array[NMain], if HaveBndU[i] is False,
                    then i-th bound is not present
    S           -   array[NMain+NSlack], scaling of the variables
    NMain       -   number of main variables
    NSlack      -   number of slack variables
    DropTol     -   drop tolerance, >=0
    
OUTPUT PARAMETERS
    X           -   point bounded with respect to constraints.
                    components corresponding to active constraints are exactly
                    equal to the boundary values.

  -- ALGLIB --
     Copyright 10.01.2012 by Bochkanov Sergey
*************************************************************************/
void filterdirection(/* Real    */ ae_vector* d,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* bndl,
     /* Boolean */ ae_vector* havebndl,
     /* Real    */ ae_vector* bndu,
     /* Boolean */ ae_vector* havebndu,
     /* Real    */ ae_vector* s,
     ae_int_t nmain,
     ae_int_t nslack,
     double droptol,
     ae_state *_state)
{
    ae_int_t i;
    double scalednorm;
    ae_bool isactive;


    scalednorm = 0.0;
    for(i=0; i<=nmain+nslack-1; i++)
    {
        scalednorm = scalednorm+ae_sqr(d->ptr.p_double[i]*s->ptr.p_double[i], _state);
    }
    scalednorm = ae_sqrt(scalednorm, _state);
    for(i=0; i<=nmain-1; i++)
    {
        ae_assert(!havebndl->ptr.p_bool[i]||ae_fp_greater_eq(x->ptr.p_double[i],bndl->ptr.p_double[i]), "FilterDirection: infeasible point", _state);
        ae_assert(!havebndu->ptr.p_bool[i]||ae_fp_less_eq(x->ptr.p_double[i],bndu->ptr.p_double[i]), "FilterDirection: infeasible point", _state);
        isactive = (havebndl->ptr.p_bool[i]&&ae_fp_eq(x->ptr.p_double[i],bndl->ptr.p_double[i]))||(havebndu->ptr.p_bool[i]&&ae_fp_eq(x->ptr.p_double[i],bndu->ptr.p_double[i]));
        if( isactive&&ae_fp_less_eq(ae_fabs(d->ptr.p_double[i]*s->ptr.p_double[i], _state),droptol*scalednorm) )
        {
            d->ptr.p_double[i] = 0.0;
        }
    }
    for(i=0; i<=nslack-1; i++)
    {
        ae_assert(ae_fp_greater_eq(x->ptr.p_double[nmain+i],0), "FilterDirection: infeasible point", _state);
        if( ae_fp_eq(x->ptr.p_double[nmain+i],0)&&ae_fp_less_eq(ae_fabs(d->ptr.p_double[nmain+i]*s->ptr.p_double[nmain+i], _state),droptol*scalednorm) )
        {
            d->ptr.p_double[nmain+i] = 0.0;
        }
    }
}


/*************************************************************************
This function returns number of bound constraints whose state was  changed
(either activated or deactivated) when making step from XPrev to X.

Constraints are considered:
* active - when we are exactly at the boundary
* inactive - when we are not at the boundary

You should note that antigradient direction is NOT taken into account when
we make decions on the constraint status.

INPUT PARAMETERS
    X           -   array[NMain+NSlack], final point.
                    Must be feasible with respect to bound constraints.
    XPrev       -   array[NMain+NSlack], initial point.
                    Must be feasible with respect to bound constraints.
    BndL        -   lower bounds, array[NMain]
                    (may contain -INF, when bound is not present)
    HaveBndL    -   array[NMain], if HaveBndL[i] is False,
                    then i-th bound is not present
    BndU        -   array[NMain], upper bounds
                    (may contain +INF, when bound is not present)
    HaveBndU    -   array[NMain], if HaveBndU[i] is False,
                    then i-th bound is not present
    NMain       -   number of main variables
    NSlack      -   number of slack variables
    
RESULT:
    number of constraints whose state was changed.

  -- ALGLIB --
     Copyright 10.01.2012 by Bochkanov Sergey
*************************************************************************/
ae_int_t numberofchangedconstraints(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* xprev,
     /* Real    */ ae_vector* bndl,
     /* Boolean */ ae_vector* havebndl,
     /* Real    */ ae_vector* bndu,
     /* Boolean */ ae_vector* havebndu,
     ae_int_t nmain,
     ae_int_t nslack,
     ae_state *) //_state)
{
    ae_int_t i;
    ae_bool statuschanged;
    ae_int_t result;


    result = 0;
    for(i=0; i<=nmain-1; i++)
    {
        if( ae_fp_neq(x->ptr.p_double[i],xprev->ptr.p_double[i]) )
        {
            statuschanged = ae_false;
            if( havebndl->ptr.p_bool[i]&&(ae_fp_eq(x->ptr.p_double[i],bndl->ptr.p_double[i])||ae_fp_eq(xprev->ptr.p_double[i],bndl->ptr.p_double[i])) )
            {
                statuschanged = ae_true;
            }
            if( havebndu->ptr.p_bool[i]&&(ae_fp_eq(x->ptr.p_double[i],bndu->ptr.p_double[i])||ae_fp_eq(xprev->ptr.p_double[i],bndu->ptr.p_double[i])) )
            {
                statuschanged = ae_true;
            }
            if( statuschanged )
            {
                result = result+1;
            }
        }
    }
    for(i=0; i<=nslack-1; i++)
    {
        if( ae_fp_neq(x->ptr.p_double[nmain+i],xprev->ptr.p_double[nmain+i])&&(ae_fp_eq(x->ptr.p_double[nmain+i],0)||ae_fp_eq(xprev->ptr.p_double[nmain+i],0)) )
        {
            result = result+1;
        }
    }
    return result;
}


/*************************************************************************
This function finds feasible point of  (NMain+NSlack)-dimensional  problem
subject to NMain explicit boundary constraints (some  constraints  can  be
omitted), NSlack implicit non-negativity constraints,  K  linear  equality
constraints.

INPUT PARAMETERS
    X           -   array[NMain+NSlack], initial point.
    BndL        -   lower bounds, array[NMain]
                    (may contain -INF, when bound is not present)
    HaveBndL    -   array[NMain], if HaveBndL[i] is False,
                    then i-th bound is not present
    BndU        -   array[NMain], upper bounds
                    (may contain +INF, when bound is not present)
    HaveBndU    -   array[NMain], if HaveBndU[i] is False,
                    then i-th bound is not present
    NMain       -   number of main variables
    NSlack      -   number of slack variables
    CE          -   array[K,NMain+NSlack+1], equality  constraints CE*x=b.
                    Rows contain constraints, first  NMain+NSlack  columns
                    contain coefficients before X[], last  column  contain
                    right part.
    K           -   number of linear constraints
    EpsI        -   infeasibility (error in the right part) allowed in the
                    solution

OUTPUT PARAMETERS:
    X           -   feasible point or best infeasible point found before
                    algorithm termination
    QPIts       -   number of QP iterations (for debug purposes)
    GPAIts      -   number of GPA iterations (for debug purposes)
    
RESULT:
    True in case X is feasible, False - if it is infeasible.

  -- ALGLIB --
     Copyright 20.01.2012 by Bochkanov Sergey
*************************************************************************/
ae_bool findfeasiblepoint(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* bndl,
     /* Boolean */ ae_vector* havebndl,
     /* Real    */ ae_vector* bndu,
     /* Boolean */ ae_vector* havebndu,
     ae_int_t nmain,
     ae_int_t nslack,
     /* Real    */ ae_matrix* ce,
     ae_int_t k,
     double epsi,
     ae_int_t* qpits,
     ae_int_t* gpaits,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix _ce;
    ae_int_t i;
    ae_int_t j;
    ae_int_t idx0;
    ae_int_t idx1;
    ae_vector permx;
    ae_vector xn;
    ae_vector xa;
    ae_vector newtonstep;
    ae_vector g;
    ae_vector pg;
    ae_matrix a;
    double armijostep;
    double armijobeststep;
    double armijobestfeas;
    double v;
    double mx;
    double feaserr;
    double feasold;
    double feasnew;
    double pgnorm;
    double vn;
    double vd;
    double stp;
    ae_int_t vartofreeze;
    double valtofreeze;
    double maxsteplen;
    ae_bool werechangesinconstraints;
    ae_bool stage1isover;
    ae_bool converged;
    ae_vector activeconstraints;
    ae_vector tmpk;
    ae_vector colnorms;
    ae_int_t nactive;
    ae_int_t nfree;
    ae_int_t nsvd;
    ae_vector p1;
    ae_vector p2;
    apbuffers buf;
    ae_vector w;
    ae_vector s;
    ae_matrix u;
    ae_matrix vt;
    ae_int_t itscount;
    ae_int_t itswithintolerance;
    ae_int_t maxitswithintolerance;
    ae_int_t gparuns;
    ae_int_t maxgparuns;
    ae_int_t maxarmijoruns;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init_copy(&_ce, ce, _state, ae_true);
    ce = &_ce;
    *qpits = 0;
    *gpaits = 0;
    ae_vector_init(&permx, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xn, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xa, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&newtonstep, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&g, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&pg, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&a, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&activeconstraints, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmpk, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&colnorms, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&p1, 0, DT_INT, _state, ae_true);
    ae_vector_init(&p2, 0, DT_INT, _state, ae_true);
    _apbuffers_init(&buf, _state, ae_true);
    ae_vector_init(&w, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&s, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&u, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&vt, 0, 0, DT_REAL, _state, ae_true);

    maxitswithintolerance = 3;
    maxgparuns = 3;
    maxarmijoruns = 5;
    *qpits = 0;
    *gpaits = 0;
    
    /*
     * Initial enforcement of the feasibility with respect to boundary constraints
     * NOTE: after this block we assume that boundary constraints are consistent.
     */
    if( !enforceboundaryconstraints(x, bndl, havebndl, bndu, havebndu, nmain, nslack, _state) )
    {
        result = ae_false;
        ae_frame_leave(_state);
        return result;
    }
    if( k==0 )
    {
        
        /*
         * No linear constraints, we can exit right now
         */
        result = ae_true;
        ae_frame_leave(_state);
        return result;
    }
    
    /*
     * Scale rows of CE in such way that max(CE[i,0..nmain+nslack-1])=1 for any i=0..k-1
     */
    for(i=0; i<=k-1; i++)
    {
        v = 0.0;
        for(j=0; j<=nmain+nslack-1; j++)
        {
            v = ae_maxreal(v, ae_fabs(ce->ptr.pp_double[i][j], _state), _state);
        }
        if( ae_fp_neq(v,0) )
        {
            v = 1/v;
            ae_v_muld(&ce->ptr.pp_double[i][0], 1, ae_v_len(0,nmain+nslack), v);
        }
    }
    
    /*
     * Allocate temporaries
     */
    ae_vector_set_length(&xn, nmain+nslack, _state);
    ae_vector_set_length(&xa, nmain+nslack, _state);
    ae_vector_set_length(&permx, nmain+nslack, _state);
    ae_vector_set_length(&g, nmain+nslack, _state);
    ae_vector_set_length(&pg, nmain+nslack, _state);
    ae_vector_set_length(&tmpk, k, _state);
    ae_matrix_set_length(&a, k, nmain+nslack, _state);
    ae_vector_set_length(&activeconstraints, nmain+nslack, _state);
    ae_vector_set_length(&newtonstep, nmain+nslack, _state);
    ae_vector_set_length(&s, nmain+nslack, _state);
    ae_vector_set_length(&colnorms, nmain+nslack, _state);
    for(i=0; i<=nmain+nslack-1; i++)
    {
        s.ptr.p_double[i] = 1.0;
        colnorms.ptr.p_double[i] = 0.0;
        for(j=0; j<=k-1; j++)
        {
            colnorms.ptr.p_double[i] = colnorms.ptr.p_double[i]+ae_sqr(ce->ptr.pp_double[j][i], _state);
        }
    }
    
    /*
     * K>0, we have linear equality constraints combined with bound constraints.
     *
     * Try to find feasible point as minimizer of the quadratic function
     *     F(x) = 0.5*||CE*x-b||^2 = 0.5*x'*(CE'*CE)*x - (b'*CE)*x + 0.5*b'*b
     * subject to boundary constraints given by BL, BU and non-negativity of
     * the slack variables. BTW, we drop constant term because it does not
     * actually influences on the solution.
     *
     * Below we will assume that K>0.
     */
    itswithintolerance = 0;
    itscount = 0;
    for(;;)
    {
        
        /*
         * Stage 0: check for exact convergence
         */
        converged = ae_true;
        feaserr = 0;
        for(i=0; i<=k-1; i++)
        {
            
            /*
             * Calculate:
             * * V - error in the right part
             * * MX - maximum term in the left part
             *
             * Terminate if error in the right part is not greater than 100*Eps*MX.
             *
             * IMPORTANT: we must perform check for non-strict inequality, i.e. to use <= instead of <.
             *            it will allow us to easily handle situations with zero rows of CE.
             */
            mx = 0;
            v = -ce->ptr.pp_double[i][nmain+nslack];
            for(j=0; j<=nmain+nslack-1; j++)
            {
                mx = ae_maxreal(mx, ae_fabs(ce->ptr.pp_double[i][j]*x->ptr.p_double[j], _state), _state);
                v = v+ce->ptr.pp_double[i][j]*x->ptr.p_double[j];
            }
            feaserr = feaserr+ae_sqr(v, _state);
            converged = converged&&ae_fp_less_eq(ae_fabs(v, _state),100*ae_machineepsilon*mx);
        }
        feaserr = ae_sqrt(feaserr, _state);
        if( converged )
        {
            result = ae_fp_less_eq(feaserr,epsi);
            ae_frame_leave(_state);
            return result;
        }
        
        /*
         * Stage 1: equality constrained quadratic programming
         *
         * * treat active bound constraints as equality ones (constraint is considered 
         *   active when we are at the boundary, independently of the antigradient direction)
         * * calculate unrestricted Newton step to point XM (which may be infeasible)
         *   calculate MaxStepLen = largest step in direction of XM which retains feasibility.
         * * perform bounded step from X to XN:
         *   a) XN=XM                  (if XM is feasible)
         *   b) XN=X-MaxStepLen*(XM-X) (otherwise)
         * * X := XN
         * * if XM (Newton step subject to currently active constraints) was feasible, goto Stage 2
         * * repeat Stage 1
         *
         * NOTE 1: in order to solve constrained qudratic subproblem we will have to reorder
         *         variables in such way that ones corresponding to inactive constraints will
         *         be first, and active ones will be last in the list. CE and X are now
         *                                                       [ xi ]
         *         separated into two parts: CE = [CEi CEa], x = [    ], where CEi/Xi correspond
         *                                                       [ xa ]
         *         to INACTIVE constraints, and CEa/Xa correspond to the ACTIVE ones.
         *
         *         Now, instead of F=0.5*x'*(CE'*CE)*x - (b'*CE)*x + 0.5*b'*b, we have
         *         F(xi) = 0.5*(CEi*xi,CEi*xi) + (CEa*xa-b,CEi*xi) + (0.5*CEa*xa-b,CEa*xa).
         *         Here xa is considered constant, i.e. we optimize with respect to xi, leaving xa fixed.
         *
         *         We can solve it by performing SVD of CEi and calculating pseudoinverse of the
         *         Hessian matrix. Of course, we do NOT calculate pseudoinverse explicitly - we
         *         just use singular vectors to perform implicit multiplication by it.
         *
         */
        for(;;)
        {
            
            /*
             * Calculate G - gradient subject to equality constraints,
             * multiply it by inverse of the Hessian diagonal to obtain initial
             * step vector.
             *
             * Bound step subject to constraints which can be activated,
             * run Armijo search with increasing step size.
             * Search is terminated when feasibility error stops to decrease.
             *
             * NOTE: it is important to test for "stops to decrease" instead
             * of "starts to increase" in order to correctly handle cases with
             * zero CE.
             */
            armijobeststep = 0.0;
            armijobestfeas = 0.0;
            for(i=0; i<=nmain+nslack-1; i++)
            {
                g.ptr.p_double[i] = 0;
            }
            for(i=0; i<=k-1; i++)
            {
                v = ae_v_dotproduct(&ce->ptr.pp_double[i][0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
                v = v-ce->ptr.pp_double[i][nmain+nslack];
                armijobestfeas = armijobestfeas+ae_sqr(v, _state);
                ae_v_addd(&g.ptr.p_double[0], 1, &ce->ptr.pp_double[i][0], 1, ae_v_len(0,nmain+nslack-1), v);
            }
            armijobestfeas = ae_sqrt(armijobestfeas, _state);
            for(i=0; i<=nmain-1; i++)
            {
                if( havebndl->ptr.p_bool[i]&&ae_fp_eq(x->ptr.p_double[i],bndl->ptr.p_double[i]) )
                {
                    g.ptr.p_double[i] = 0.0;
                }
                if( havebndu->ptr.p_bool[i]&&ae_fp_eq(x->ptr.p_double[i],bndu->ptr.p_double[i]) )
                {
                    g.ptr.p_double[i] = 0.0;
                }
            }
            for(i=0; i<=nslack-1; i++)
            {
                if( ae_fp_eq(x->ptr.p_double[nmain+i],0.0) )
                {
                    g.ptr.p_double[nmain+i] = 0.0;
                }
            }
            v = 0.0;
            for(i=0; i<=nmain+nslack-1; i++)
            {
                if( ae_fp_neq(ae_sqr(colnorms.ptr.p_double[i], _state),0) )
                {
                    newtonstep.ptr.p_double[i] = -g.ptr.p_double[i]/ae_sqr(colnorms.ptr.p_double[i], _state);
                }
                else
                {
                    newtonstep.ptr.p_double[i] = 0.0;
                }
                v = v+ae_sqr(newtonstep.ptr.p_double[i], _state);
            }
            if( ae_fp_eq(v,0) )
            {
                
                /*
                 * Constrained gradient is zero, QP iterations are over
                 */
                break;
            }
            calculatestepbound(x, &newtonstep, 1.0, bndl, havebndl, bndu, havebndu, nmain, nslack, &vartofreeze, &valtofreeze, &maxsteplen, _state);
            if( vartofreeze>=0&&ae_fp_eq(maxsteplen,0) )
            {
                
                /*
                 * Can not perform step, QP iterations are over
                 */
                break;
            }
            if( vartofreeze>=0 )
            {
                armijostep = ae_minreal(1.0, maxsteplen, _state);
            }
            else
            {
                armijostep = 1;
            }
            for(;;)
            {
                ae_v_move(&xa.ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
                ae_v_addd(&xa.ptr.p_double[0], 1, &newtonstep.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1), armijostep);
                enforceboundaryconstraints(&xa, bndl, havebndl, bndu, havebndu, nmain, nslack, _state);
                feaserr = 0.0;
                for(i=0; i<=k-1; i++)
                {
                    v = ae_v_dotproduct(&ce->ptr.pp_double[i][0], 1, &xa.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
                    v = v-ce->ptr.pp_double[i][nmain+nslack];
                    feaserr = feaserr+ae_sqr(v, _state);
                }
                feaserr = ae_sqrt(feaserr, _state);
                if( ae_fp_greater_eq(feaserr,armijobestfeas) )
                {
                    break;
                }
                armijobestfeas = feaserr;
                armijobeststep = armijostep;
                armijostep = 2.0*armijostep;
            }
            ae_v_addd(&x->ptr.p_double[0], 1, &newtonstep.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1), armijobeststep);
            enforceboundaryconstraints(x, bndl, havebndl, bndu, havebndu, nmain, nslack, _state);
            
            /*
             * Determine number of active and free constraints
             */
            nactive = 0;
            for(i=0; i<=nmain-1; i++)
            {
                activeconstraints.ptr.p_double[i] = 0;
                if( havebndl->ptr.p_bool[i]&&ae_fp_eq(x->ptr.p_double[i],bndl->ptr.p_double[i]) )
                {
                    activeconstraints.ptr.p_double[i] = 1;
                }
                if( havebndu->ptr.p_bool[i]&&ae_fp_eq(x->ptr.p_double[i],bndu->ptr.p_double[i]) )
                {
                    activeconstraints.ptr.p_double[i] = 1;
                }
                if( ae_fp_greater(activeconstraints.ptr.p_double[i],0) )
                {
                    nactive = nactive+1;
                }
            }
            for(i=0; i<=nslack-1; i++)
            {
                activeconstraints.ptr.p_double[nmain+i] = 0;
                if( ae_fp_eq(x->ptr.p_double[nmain+i],0.0) )
                {
                    activeconstraints.ptr.p_double[nmain+i] = 1;
                }
                if( ae_fp_greater(activeconstraints.ptr.p_double[nmain+i],0) )
                {
                    nactive = nactive+1;
                }
            }
            nfree = nmain+nslack-nactive;
            if( nfree==0 )
            {
                break;
            }
            *qpits = *qpits+1;
            
            /*
             * Reorder variables
             */
            tagsortbuf(&activeconstraints, nmain+nslack, &p1, &p2, &buf, _state);
            for(i=0; i<=k-1; i++)
            {
                for(j=0; j<=nmain+nslack-1; j++)
                {
                    a.ptr.pp_double[i][j] = ce->ptr.pp_double[i][j];
                }
            }
            for(j=0; j<=nmain+nslack-1; j++)
            {
                permx.ptr.p_double[j] = x->ptr.p_double[j];
            }
            for(j=0; j<=nmain+nslack-1; j++)
            {
                if( p2.ptr.p_int[j]!=j )
                {
                    idx0 = p2.ptr.p_int[j];
                    idx1 = j;
                    for(i=0; i<=k-1; i++)
                    {
                        v = a.ptr.pp_double[i][idx0];
                        a.ptr.pp_double[i][idx0] = a.ptr.pp_double[i][idx1];
                        a.ptr.pp_double[i][idx1] = v;
                    }
                    v = permx.ptr.p_double[idx0];
                    permx.ptr.p_double[idx0] = permx.ptr.p_double[idx1];
                    permx.ptr.p_double[idx1] = v;
                }
            }
            
            /*
             * Calculate (unprojected) gradient:
             * G(xi) = CEi'*(CEi*xi + CEa*xa - b)
             */
            for(i=0; i<=nfree-1; i++)
            {
                g.ptr.p_double[i] = 0;
            }
            for(i=0; i<=k-1; i++)
            {
                v = ae_v_dotproduct(&a.ptr.pp_double[i][0], 1, &permx.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
                tmpk.ptr.p_double[i] = v-ce->ptr.pp_double[i][nmain+nslack];
            }
            for(i=0; i<=k-1; i++)
            {
                v = tmpk.ptr.p_double[i];
                ae_v_addd(&g.ptr.p_double[0], 1, &a.ptr.pp_double[i][0], 1, ae_v_len(0,nfree-1), v);
            }
            
            /*
             * Calculate Newton step using SVD of CEi:
             *     F(xi)  = 0.5*xi'*H*xi + g'*xi    (Taylor decomposition)
             *     XN     = -H^(-1)*g               (new point, solution of the QP subproblem)
             *     H      = CEi'*CEi                
             *     CEi    = U*W*V'                  (SVD of CEi)
             *     H      = V*W^2*V'                 
             *     H^(-1) = V*W^(-2)*V'
             *     step     = -V*W^(-2)*V'*g          (it is better to perform multiplication from right to left)
             *
             * NOTE 1: we do NOT need left singular vectors to perform Newton step.
             */
            nsvd = ae_minint(k, nfree, _state);
            if( !rmatrixsvd(&a, k, nfree, 0, 1, 2, &w, &u, &vt, _state) )
            {
                result = ae_false;
                ae_frame_leave(_state);
                return result;
            }
            for(i=0; i<=nsvd-1; i++)
            {
                v = ae_v_dotproduct(&vt.ptr.pp_double[i][0], 1, &g.ptr.p_double[0], 1, ae_v_len(0,nfree-1));
                tmpk.ptr.p_double[i] = v;
            }
            for(i=0; i<=nsvd-1; i++)
            {
                
                /*
                 * It is important to have strict ">" in order to correctly 
                 * handle zero singular values.
                 */
                if( ae_fp_greater(ae_sqr(w.ptr.p_double[i], _state),ae_sqr(w.ptr.p_double[0], _state)*(nmain+nslack)*ae_machineepsilon) )
                {
                    tmpk.ptr.p_double[i] = tmpk.ptr.p_double[i]/ae_sqr(w.ptr.p_double[i], _state);
                }
                else
                {
                    tmpk.ptr.p_double[i] = 0;
                }
            }
            for(i=0; i<=nmain+nslack-1; i++)
            {
                newtonstep.ptr.p_double[i] = 0;
            }
            for(i=0; i<=nsvd-1; i++)
            {
                v = tmpk.ptr.p_double[i];
                ae_v_subd(&newtonstep.ptr.p_double[0], 1, &vt.ptr.pp_double[i][0], 1, ae_v_len(0,nfree-1), v);
            }
            for(j=nmain+nslack-1; j>=0; j--)
            {
                if( p2.ptr.p_int[j]!=j )
                {
                    idx0 = p2.ptr.p_int[j];
                    idx1 = j;
                    v = newtonstep.ptr.p_double[idx0];
                    newtonstep.ptr.p_double[idx0] = newtonstep.ptr.p_double[idx1];
                    newtonstep.ptr.p_double[idx1] = v;
                }
            }
            
            /*
             * NewtonStep contains Newton step subject to active bound constraints.
             *
             * Such step leads us to the minimizer of the equality constrained F,
             * but such minimizer may be infeasible because some constraints which
             * are inactive at the initial point can be violated at the solution.
             *
             * Thus, we perform optimization in two stages:
             * a) perform bounded Newton step, i.e. step in the Newton direction
             *    until activation of the first constraint
             * b) in case (MaxStepLen>0)and(MaxStepLen<1), perform additional iteration
             *    of the Armijo line search in the rest of the Newton direction.
             */
            calculatestepbound(x, &newtonstep, 1.0, bndl, havebndl, bndu, havebndu, nmain, nslack, &vartofreeze, &valtofreeze, &maxsteplen, _state);
            if( vartofreeze>=0&&ae_fp_eq(maxsteplen,0) )
            {
                
                /*
                 * Activation of the constraints prevent us from performing step,
                 * QP iterations are over
                 */
                break;
            }
            if( vartofreeze>=0 )
            {
                v = ae_minreal(1.0, maxsteplen, _state);
            }
            else
            {
                v = 1.0;
            }
            ae_v_moved(&xn.ptr.p_double[0], 1, &newtonstep.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1), v);
            ae_v_add(&xn.ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
            postprocessboundedstep(&xn, x, bndl, havebndl, bndu, havebndu, nmain, nslack, vartofreeze, valtofreeze, v, maxsteplen, _state);
            if( ae_fp_greater(maxsteplen,0)&&ae_fp_less(maxsteplen,1) )
            {
                
                /*
                 * Newton step was restricted by activation of the constraints,
                 * perform Armijo iteration.
                 *
                 * Initial estimate for best step is zero step. We try different
                 * step sizes, from the 1-MaxStepLen (residual of the full Newton
                 * step) to progressively smaller and smaller steps.
                 */
                armijobeststep = 0.0;
                armijobestfeas = 0.0;
                for(i=0; i<=k-1; i++)
                {
                    v = ae_v_dotproduct(&ce->ptr.pp_double[i][0], 1, &xn.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
                    v = v-ce->ptr.pp_double[i][nmain+nslack];
                    armijobestfeas = armijobestfeas+ae_sqr(v, _state);
                }
                armijobestfeas = ae_sqrt(armijobestfeas, _state);
                armijostep = 1-maxsteplen;
                for(j=0; j<=maxarmijoruns-1; j++)
                {
                    ae_v_move(&xa.ptr.p_double[0], 1, &xn.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
                    ae_v_addd(&xa.ptr.p_double[0], 1, &newtonstep.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1), armijostep);
                    enforceboundaryconstraints(&xa, bndl, havebndl, bndu, havebndu, nmain, nslack, _state);
                    feaserr = 0.0;
                    for(i=0; i<=k-1; i++)
                    {
                        v = ae_v_dotproduct(&ce->ptr.pp_double[i][0], 1, &xa.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
                        v = v-ce->ptr.pp_double[i][nmain+nslack];
                        feaserr = feaserr+ae_sqr(v, _state);
                    }
                    feaserr = ae_sqrt(feaserr, _state);
                    if( ae_fp_less(feaserr,armijobestfeas) )
                    {
                        armijobestfeas = feaserr;
                        armijobeststep = armijostep;
                    }
                    armijostep = 0.5*armijostep;
                }
                ae_v_move(&xa.ptr.p_double[0], 1, &xn.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
                ae_v_addd(&xa.ptr.p_double[0], 1, &newtonstep.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1), armijobeststep);
                enforceboundaryconstraints(&xa, bndl, havebndl, bndu, havebndu, nmain, nslack, _state);
            }
            else
            {
                
                /*
                 * Armijo iteration is not performed
                 */
                ae_v_move(&xa.ptr.p_double[0], 1, &xn.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
            }
            stage1isover = ae_fp_greater_eq(maxsteplen,1)||ae_fp_eq(maxsteplen,0);
            
            /*
             * Calculate feasibility errors for old and new X.
             * These quantinies are used for debugging purposes only.
             * However, we can leave them in release code because performance impact is insignificant.
             *
             * Update X. Exit if needed.
             */
            feasold = 0;
            feasnew = 0;
            for(i=0; i<=k-1; i++)
            {
                v = ae_v_dotproduct(&ce->ptr.pp_double[i][0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
                feasold = feasold+ae_sqr(v-ce->ptr.pp_double[i][nmain+nslack], _state);
                v = ae_v_dotproduct(&ce->ptr.pp_double[i][0], 1, &xa.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
                feasnew = feasnew+ae_sqr(v-ce->ptr.pp_double[i][nmain+nslack], _state);
            }
            feasold = ae_sqrt(feasold, _state);
            feasnew = ae_sqrt(feasnew, _state);
            if( ae_fp_greater_eq(feasnew,feasold) )
            {
                break;
            }
            ae_v_move(&x->ptr.p_double[0], 1, &xa.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
            if( stage1isover )
            {
                break;
            }
        }
        
        /*
         * Stage 2: gradient projection algorithm (GPA)
         *
         * * calculate feasibility error (with respect to linear equality constraints)
         * * calculate gradient G of F, project it into feasible area (G => PG)
         * * exit if norm(PG) is exactly zero or feasibility error is smaller than EpsC
         * * let XM be exact minimum of F along -PG (XM may be infeasible).
         *   calculate MaxStepLen = largest step in direction of -PG which retains feasibility.
         * * perform bounded step from X to XN:
         *   a) XN=XM              (if XM is feasible)
         *   b) XN=X-MaxStepLen*PG (otherwise)
         * * X := XN
         * * stop after specified number of iterations or when no new constraints was activated
         *
         * NOTES:
         * * grad(F) = (CE'*CE)*x - (b'*CE)^T
         * * CE[i] denotes I-th row of CE
         * * XM = X+stp*(-PG) where stp=(grad(F(X)),PG)/(CE*PG,CE*PG).
         *   Here PG is a projected gradient, but in fact it can be arbitrary non-zero 
         *   direction vector - formula for minimum of F along PG still will be correct.
         */
        werechangesinconstraints = ae_false;
        for(gparuns=1; gparuns<=k; gparuns++)
        {
            
            /*
             * calculate feasibility error and G
             */
            feaserr = 0;
            for(i=0; i<=nmain+nslack-1; i++)
            {
                g.ptr.p_double[i] = 0;
            }
            for(i=0; i<=k-1; i++)
            {
                
                /*
                 * G += CE[i]^T * (CE[i]*x-b[i])
                 */
                v = ae_v_dotproduct(&ce->ptr.pp_double[i][0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
                v = v-ce->ptr.pp_double[i][nmain+nslack];
                feaserr = feaserr+ae_sqr(v, _state);
                ae_v_addd(&g.ptr.p_double[0], 1, &ce->ptr.pp_double[i][0], 1, ae_v_len(0,nmain+nslack-1), v);
            }
            
            /*
             * project G, filter it (strip numerical noise)
             */
            ae_v_move(&pg.ptr.p_double[0], 1, &g.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
            projectgradientintobc(x, &pg, bndl, havebndl, bndu, havebndu, nmain, nslack, _state);
            filterdirection(&pg, x, bndl, havebndl, bndu, havebndu, &s, nmain, nslack, 1.0E-9, _state);
            for(i=0; i<=nmain+nslack-1; i++)
            {
                if( ae_fp_neq(ae_sqr(colnorms.ptr.p_double[i], _state),0) )
                {
                    pg.ptr.p_double[i] = pg.ptr.p_double[i]/ae_sqr(colnorms.ptr.p_double[i], _state);
                }
                else
                {
                    pg.ptr.p_double[i] = 0.0;
                }
            }
            
            /*
             * Check GNorm and feasibility.
             * Exit when GNorm is exactly zero.
             */
            pgnorm = ae_v_dotproduct(&pg.ptr.p_double[0], 1, &pg.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
            feaserr = ae_sqrt(feaserr, _state);
            pgnorm = ae_sqrt(pgnorm, _state);
            if( ae_fp_eq(pgnorm,0) )
            {
                result = ae_fp_less_eq(feaserr,epsi);
                ae_frame_leave(_state);
                return result;
            }
            
            /*
             * calculate planned step length
             */
            vn = ae_v_dotproduct(&g.ptr.p_double[0], 1, &pg.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
            vd = 0;
            for(i=0; i<=k-1; i++)
            {
                v = ae_v_dotproduct(&ce->ptr.pp_double[i][0], 1, &pg.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
                vd = vd+ae_sqr(v, _state);
            }
            stp = vn/vd;
            
            /*
             * Calculate step bound.
             * Perform bounded step and post-process it
             */
            calculatestepbound(x, &pg, -1.0, bndl, havebndl, bndu, havebndu, nmain, nslack, &vartofreeze, &valtofreeze, &maxsteplen, _state);
            if( vartofreeze>=0&&ae_fp_eq(maxsteplen,0) )
            {
                result = ae_false;
                ae_frame_leave(_state);
                return result;
            }
            if( vartofreeze>=0 )
            {
                v = ae_minreal(stp, maxsteplen, _state);
            }
            else
            {
                v = stp;
            }
            ae_v_move(&xn.ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
            ae_v_subd(&xn.ptr.p_double[0], 1, &pg.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1), v);
            postprocessboundedstep(&xn, x, bndl, havebndl, bndu, havebndu, nmain, nslack, vartofreeze, valtofreeze, v, maxsteplen, _state);
            
            /*
             * update X
             * check stopping criteria
             */
            werechangesinconstraints = werechangesinconstraints||numberofchangedconstraints(&xn, x, bndl, havebndl, bndu, havebndu, nmain, nslack, _state)>0;
            ae_v_move(&x->ptr.p_double[0], 1, &xn.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
            *gpaits = *gpaits+1;
            if( !werechangesinconstraints )
            {
                break;
            }
        }
        
        /*
         * Stage 3: decide to stop algorithm or not to stop
         *
         * 1. we can stop when last GPA run did NOT changed constraints status.
         *    It means that we've found final set of the active constraints even
         *    before GPA made its run. And it means that Newton step moved us to
         *    the minimum subject to the present constraints.
         *    Depending on feasibility error, True or False is returned.
         */
        feaserr = 0;
        for(i=0; i<=k-1; i++)
        {
            v = ae_v_dotproduct(&ce->ptr.pp_double[i][0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
            v = v-ce->ptr.pp_double[i][nmain+nslack];
            feaserr = feaserr+ae_sqr(v, _state);
        }
        feaserr = ae_sqrt(feaserr, _state);
        if( ae_fp_less_eq(feaserr,epsi) )
        {
            itswithintolerance = itswithintolerance+1;
        }
        else
        {
            itswithintolerance = 0;
        }
        if( !werechangesinconstraints||itswithintolerance>=maxitswithintolerance )
        {
            result = ae_fp_less_eq(feaserr,epsi);
            ae_frame_leave(_state);
            return result;
        }
        itscount = itscount+1;
    }
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
    This function check, that input derivatives are right. First it scale
parameters DF0 and DF1 from segment [A;B] to [0;1]. Than it build Hermite
spline and derivative of it in 0,5. Search scale as Max(DF0,DF1, |F0-F1|).
Right derivative has to satisfy condition:
    |H-F|/S<=0,01, |H'-F'|/S<=0,01.
    
INPUT PARAMETERS:
    F0  -   function's value in X-TestStep point;
    DF0 -   derivative's value in X-TestStep point;
    F1  -   function's value in X+TestStep point;
    DF1 -   derivative's value in X+TestStep point;
    F   -   testing function's value;
    DF  -   testing derivative's value;
   Width-   width of verification segment.

RESULT:
    If input derivatives is right then function returns true, else 
    function returns false.
    
  -- ALGLIB --
     Copyright 29.05.2012 by Bochkanov Sergey
*************************************************************************/
ae_bool derivativecheck(double f0,
     double df0,
     double f1,
     double df1,
     double f,
     double df,
     double width,
     ae_state *_state)
{
    double s;
    double h;
    double dh;
    ae_bool result;


    df = width*df;
    df0 = width*df0;
    df1 = width*df1;
    s = ae_maxreal(ae_maxreal(ae_fabs(df0, _state), ae_fabs(df1, _state), _state), ae_fabs(f1-f0, _state), _state);
    h = 0.5*f0+0.125*df0+0.5*f1-0.125*df1;
    dh = -1.5*f0-0.25*df0+1.5*f1-0.25*df1;
    if( ae_fp_neq(s,0) )
    {
        if( ae_fp_greater(ae_fabs(h-f, _state)/s,0.001)||ae_fp_greater(ae_fabs(dh-df, _state)/s,0.001) )
        {
            result = ae_false;
            return result;
        }
    }
    else
    {
        if( ae_fp_neq(h-f,0.0)||ae_fp_neq(dh-df,0.0) )
        {
            result = ae_false;
            return result;
        }
    }
    result = ae_true;
    return result;
}




/*************************************************************************
        NONLINEAR CONJUGATE GRADIENT METHOD

DESCRIPTION:
The subroutine minimizes function F(x) of N arguments by using one of  the
nonlinear conjugate gradient methods.

These CG methods are globally convergent (even on non-convex functions) as
long as grad(f) is Lipschitz continuous in  a  some  neighborhood  of  the
L = { x : f(x)<=f(x0) }.


REQUIREMENTS:
Algorithm will request following information during its operation:
* function value F and its gradient G (simultaneously) at given point X


USAGE:
1. User initializes algorithm state with MinCGCreate() call
2. User tunes solver parameters with MinCGSetCond(), MinCGSetStpMax() and
   other functions
3. User calls MinCGOptimize() function which takes algorithm  state   and
   pointer (delegate, etc.) to callback function which calculates F/G.
4. User calls MinCGResults() to get solution
5. Optionally, user may call MinCGRestartFrom() to solve another  problem
   with same N but another starting point and/or another function.
   MinCGRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   problem dimension, N>0:
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    X       -   starting point, array[0..N-1].

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 25.03.2010 by Bochkanov Sergey
*************************************************************************/
void mincgcreate(ae_int_t n,
     /* Real    */ ae_vector* x,
     mincgstate* state,
     ae_state *_state)
{

    _mincgstate_clear(state);

    ae_assert(n>=1, "MinCGCreate: N too small!", _state);
    ae_assert(x->cnt>=n, "MinCGCreate: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "MinCGCreate: X contains infinite or NaN values!", _state);
    mincg_mincginitinternal(n, 0.0, state, _state);
    mincgrestartfrom(state, x, _state);
}


/*************************************************************************
The subroutine is finite difference variant of MinCGCreate(). It uses
finite differences in order to differentiate target function.

Description below contains information which is specific to this function
only. We recommend to read comments on MinCGCreate() in order to get more
information about creation of CG optimizer.

INPUT PARAMETERS:
    N       -   problem dimension, N>0:
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    X       -   starting point, array[0..N-1].
    DiffStep-   differentiation step, >0

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTES:
1. algorithm uses 4-point central formula for differentiation.
2. differentiation step along I-th axis is equal to DiffStep*S[I] where
   S[] is scaling vector which can be set by MinCGSetScale() call.
3. we recommend you to use moderate values of  differentiation  step.  Too
   large step will result in too large truncation  errors, while too small
   step will result in too large numerical  errors.  1.0E-6  can  be  good
   value to start with.
4. Numerical  differentiation  is   very   inefficient  -   one   gradient
   calculation needs 4*N function evaluations. This function will work for
   any N - either small (1...10), moderate (10...100) or  large  (100...).
   However, performance penalty will be too severe for any N's except  for
   small ones.
   We should also say that code which relies on numerical  differentiation
   is  less  robust  and  precise.  L-BFGS  needs  exact  gradient values.
   Imprecise  gradient may slow down  convergence,  especially  on  highly
   nonlinear problems.
   Thus  we  recommend to use this function for fast prototyping on small-
   dimensional problems only, and to implement analytical gradient as soon
   as possible.

  -- ALGLIB --
     Copyright 16.05.2011 by Bochkanov Sergey
*************************************************************************/
void mincgcreatef(ae_int_t n,
     /* Real    */ ae_vector* x,
     double diffstep,
     mincgstate* state,
     ae_state *_state)
{

    _mincgstate_clear(state);

    ae_assert(n>=1, "MinCGCreateF: N too small!", _state);
    ae_assert(x->cnt>=n, "MinCGCreateF: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "MinCGCreateF: X contains infinite or NaN values!", _state);
    ae_assert(ae_isfinite(diffstep, _state), "MinCGCreateF: DiffStep is infinite or NaN!", _state);
    ae_assert(ae_fp_greater(diffstep,0), "MinCGCreateF: DiffStep is non-positive!", _state);
    mincg_mincginitinternal(n, diffstep, state, _state);
    mincgrestartfrom(state, x, _state);
}


/*************************************************************************
This function sets stopping conditions for CG optimization algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    EpsG    -   >=0
                The  subroutine  finishes  its  work   if   the  condition
                |v|<EpsG is satisfied, where:
                * |.| means Euclidian norm
                * v - scaled gradient vector, v[i]=g[i]*s[i]
                * g - gradient
                * s - scaling coefficients set by MinCGSetScale()
    EpsF    -   >=0
                The  subroutine  finishes  its work if on k+1-th iteration
                the  condition  |F(k+1)-F(k)|<=EpsF*max{|F(k)|,|F(k+1)|,1}
                is satisfied.
    EpsX    -   >=0
                The subroutine finishes its work if  on  k+1-th  iteration
                the condition |v|<=EpsX is fulfilled, where:
                * |.| means Euclidian norm
                * v - scaled step vector, v[i]=dx[i]/s[i]
                * dx - ste pvector, dx=X(k+1)-X(k)
                * s - scaling coefficients set by MinCGSetScale()
    MaxIts  -   maximum number of iterations. If MaxIts=0, the  number  of
                iterations is unlimited.

Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to
automatic stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsetcond(mincgstate* state,
     double epsg,
     double epsf,
     double epsx,
     ae_int_t maxits,
     ae_state *_state)
{


    ae_assert(ae_isfinite(epsg, _state), "MinCGSetCond: EpsG is not finite number!", _state);
    ae_assert(ae_fp_greater_eq(epsg,0), "MinCGSetCond: negative EpsG!", _state);
    ae_assert(ae_isfinite(epsf, _state), "MinCGSetCond: EpsF is not finite number!", _state);
    ae_assert(ae_fp_greater_eq(epsf,0), "MinCGSetCond: negative EpsF!", _state);
    ae_assert(ae_isfinite(epsx, _state), "MinCGSetCond: EpsX is not finite number!", _state);
    ae_assert(ae_fp_greater_eq(epsx,0), "MinCGSetCond: negative EpsX!", _state);
    ae_assert(maxits>=0, "MinCGSetCond: negative MaxIts!", _state);
    if( ((ae_fp_eq(epsg,0)&&ae_fp_eq(epsf,0))&&ae_fp_eq(epsx,0))&&maxits==0 )
    {
        epsx = 1.0E-6;
    }
    state->epsg = epsg;
    state->epsf = epsf;
    state->epsx = epsx;
    state->maxits = maxits;
}


/*************************************************************************
This function sets scaling coefficients for CG optimizer.

ALGLIB optimizers use scaling matrices to test stopping  conditions  (step
size and gradient are scaled before comparison with tolerances).  Scale of
the I-th variable is a translation invariant measure of:
a) "how large" the variable is
b) how large the step should be to make significant changes in the function

Scaling is also used by finite difference variant of CG optimizer  -  step
along I-th axis is equal to DiffStep*S[I].

In   most   optimizers  (and  in  the  CG  too)  scaling is NOT a form  of
preconditioning. It just  affects  stopping  conditions.  You  should  set
preconditioner by separate call to one of the MinCGSetPrec...() functions.

There  is  special  preconditioning  mode, however,  which  uses   scaling
coefficients to form diagonal preconditioning matrix. You  can  turn  this
mode on, if you want.   But  you should understand that scaling is not the
same thing as preconditioning - these are two different, although  related
forms of tuning solver.

INPUT PARAMETERS:
    State   -   structure stores algorithm state
    S       -   array[N], non-zero scaling coefficients
                S[i] may be negative, sign doesn't matter.

  -- ALGLIB --
     Copyright 14.01.2011 by Bochkanov Sergey
*************************************************************************/
void mincgsetscale(mincgstate* state,
     /* Real    */ ae_vector* s,
     ae_state *_state)
{
    ae_int_t i;


    ae_assert(s->cnt>=state->n, "MinCGSetScale: Length(S)<N", _state);
    for(i=0; i<=state->n-1; i++)
    {
        ae_assert(ae_isfinite(s->ptr.p_double[i], _state), "MinCGSetScale: S contains infinite or NAN elements", _state);
        ae_assert(ae_fp_neq(s->ptr.p_double[i],0), "MinCGSetScale: S contains zero elements", _state);
        state->s.ptr.p_double[i] = ae_fabs(s->ptr.p_double[i], _state);
    }
}


/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    NeedXRep-   whether iteration reports are needed or not

If NeedXRep is True, algorithm will call rep() callback function if  it is
provided to MinCGOptimize().

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsetxrep(mincgstate* state, ae_bool needxrep, ae_state *) //_state)
{


    state->xrep = needxrep;
}


/*************************************************************************
This function turns on/off line search reports.
These reports are described in more details in developer-only  comments on
MinCGState object.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    NeedDRep-   whether line search reports are needed or not

This function is intended for private use only. Turning it on artificially
may cause program failure.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsetdrep(mincgstate* state, ae_bool needdrep, ae_state *) //_state)
{


    state->drep = needdrep;
}


/*************************************************************************
This function sets CG algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    CGType  -   algorithm type:
                * -1    automatic selection of the best algorithm
                * 0     DY (Dai and Yuan) algorithm
                * 1     Hybrid DY-HS algorithm

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsetcgtype(mincgstate* state, ae_int_t cgtype, ae_state *_state)
{


    ae_assert(cgtype>=-1&&cgtype<=1, "MinCGSetCGType: incorrect CGType!", _state);
    if( cgtype==-1 )
    {
        cgtype = 1;
    }
    state->cgtype = cgtype;
}


/*************************************************************************
This function sets maximum step length

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    StpMax  -   maximum step length, >=0. Set StpMax to 0.0,  if you don't
                want to limit step length.

Use this subroutine when you optimize target function which contains exp()
or  other  fast  growing  functions,  and optimization algorithm makes too
large  steps  which  leads  to overflow. This function allows us to reject
steps  that  are  too  large  (and  therefore  expose  us  to the possible
overflow) without actually calculating function value at the x+stp*d.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsetstpmax(mincgstate* state, double stpmax, ae_state *_state)
{


    ae_assert(ae_isfinite(stpmax, _state), "MinCGSetStpMax: StpMax is not finite!", _state);
    ae_assert(ae_fp_greater_eq(stpmax,0), "MinCGSetStpMax: StpMax<0!", _state);
    state->stpmax = stpmax;
}


/*************************************************************************
This function allows to suggest initial step length to the CG algorithm.

Suggested  step  length  is used as starting point for the line search. It
can be useful when you have  badly  scaled  problem,  i.e.  when  ||grad||
(which is used as initial estimate for the first step) is many  orders  of
magnitude different from the desired step.

Line search  may  fail  on  such problems without good estimate of initial
step length. Imagine, for example, problem with ||grad||=10^50 and desired
step equal to 0.1 Line  search function will use 10^50  as  initial  step,
then  it  will  decrease step length by 2 (up to 20 attempts) and will get
10^44, which is still too large.

This function allows us to tell than line search should  be  started  from
some moderate step length, like 1.0, so algorithm will be able  to  detect
desired step length in a several searches.

Default behavior (when no step is suggested) is to use preconditioner,  if
it is available, to generate initial estimate of step length.

This function influences only first iteration of algorithm. It  should  be
called between MinCGCreate/MinCGRestartFrom() call and MinCGOptimize call.
Suggested step is ignored if you have preconditioner.

INPUT PARAMETERS:
    State   -   structure used to store algorithm state.
    Stp     -   initial estimate of the step length.
                Can be zero (no estimate).

  -- ALGLIB --
     Copyright 30.07.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsuggeststep(mincgstate* state, double stp, ae_state *_state)
{


    ae_assert(ae_isfinite(stp, _state), "MinCGSuggestStep: Stp is infinite or NAN", _state);
    ae_assert(ae_fp_greater_eq(stp,0), "MinCGSuggestStep: Stp<0", _state);
    state->suggestedstep = stp;
}


/*************************************************************************
Modification of the preconditioner: preconditioning is turned off.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTE:  you  can  change  preconditioner  "on  the  fly",  during algorithm
iterations.

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsetprecdefault(mincgstate* state, ae_state *) //_state)
{


    state->prectype = 0;
    state->innerresetneeded = ae_true;
}


/*************************************************************************
Modification  of  the  preconditioner:  diagonal of approximate Hessian is
used.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    D       -   diagonal of the approximate Hessian, array[0..N-1],
                (if larger, only leading N elements are used).

NOTE:  you  can  change  preconditioner  "on  the  fly",  during algorithm
iterations.

NOTE 2: D[i] should be positive. Exception will be thrown otherwise.

NOTE 3: you should pass diagonal of approximate Hessian - NOT ITS INVERSE.

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsetprecdiag(mincgstate* state,
     /* Real    */ ae_vector* d,
     ae_state *_state)
{
    ae_int_t i;


    ae_assert(d->cnt>=state->n, "MinCGSetPrecDiag: D is too short", _state);
    for(i=0; i<=state->n-1; i++)
    {
        ae_assert(ae_isfinite(d->ptr.p_double[i], _state), "MinCGSetPrecDiag: D contains infinite or NAN elements", _state);
        ae_assert(ae_fp_greater(d->ptr.p_double[i],0), "MinCGSetPrecDiag: D contains non-positive elements", _state);
    }
    mincgsetprecdiagfast(state, d, _state);
}


/*************************************************************************
Modification of the preconditioner: scale-based diagonal preconditioning.

This preconditioning mode can be useful when you  don't  have  approximate
diagonal of Hessian, but you know that your  variables  are  badly  scaled
(for  example,  one  variable is in [1,10], and another in [1000,100000]),
and most part of the ill-conditioning comes from different scales of vars.

In this case simple  scale-based  preconditioner,  with H[i] = 1/(s[i]^2),
can greatly improve convergence.

IMPRTANT: you should set scale of your variables with MinCGSetScale() call
(before or after MinCGSetPrecScale() call). Without knowledge of the scale
of your variables scale-based preconditioner will be just unit matrix.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTE:  you  can  change  preconditioner  "on  the  fly",  during algorithm
iterations.

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsetprecscale(mincgstate* state, ae_state *) //_state)
{


    state->prectype = 3;
    state->innerresetneeded = ae_true;
}


/*************************************************************************
NOTES:

1. This function has two different implementations: one which  uses  exact
   (analytical) user-supplied  gradient, and one which uses function value
   only  and  numerically  differentiates  function  in  order  to  obtain
   gradient.
   
   Depending  on  the  specific  function  used to create optimizer object
   (either MinCGCreate()  for analytical gradient  or  MinCGCreateF()  for
   numerical differentiation) you should  choose  appropriate  variant  of
   MinCGOptimize() - one which accepts function AND gradient or one  which
   accepts function ONLY.

   Be careful to choose variant of MinCGOptimize()  which  corresponds  to
   your optimization scheme! Table below lists different  combinations  of
   callback (function/gradient) passed  to  MinCGOptimize()  and  specific
   function used to create optimizer.
   

                  |         USER PASSED TO MinCGOptimize()
   CREATED WITH   |  function only   |  function and gradient
   ------------------------------------------------------------
   MinCGCreateF() |     work                FAIL
   MinCGCreate()  |     FAIL                work

   Here "FAIL" denotes inappropriate combinations  of  optimizer  creation
   function and MinCGOptimize() version. Attemps to use  such  combination
   (for  example,  to create optimizer with  MinCGCreateF()  and  to  pass
   gradient information to MinCGOptimize()) will lead to  exception  being
   thrown. Either  you  did  not  pass  gradient when it WAS needed or you
   passed gradient when it was NOT needed.

  -- ALGLIB --
     Copyright 20.04.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool mincgiteration(mincgstate* state, ae_state *_state)
{
    ae_int_t n;
    ae_int_t i;
    double betak;
    double v;
    double vv;
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
        i = state->rstate.ia.ptr.p_int[1];
        betak = state->rstate.ra.ptr.p_double[0];
        v = state->rstate.ra.ptr.p_double[1];
        vv = state->rstate.ra.ptr.p_double[2];
    }
    else
    {
        n = -983;
        i = -989;
        betak = -834;
        v = 900;
        vv = -287;
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
    if( state->rstate.stage==4 )
    {
        goto lbl_4;
    }
    if( state->rstate.stage==5 )
    {
        goto lbl_5;
    }
    if( state->rstate.stage==6 )
    {
        goto lbl_6;
    }
    if( state->rstate.stage==7 )
    {
        goto lbl_7;
    }
    if( state->rstate.stage==8 )
    {
        goto lbl_8;
    }
    if( state->rstate.stage==9 )
    {
        goto lbl_9;
    }
    if( state->rstate.stage==10 )
    {
        goto lbl_10;
    }
    if( state->rstate.stage==11 )
    {
        goto lbl_11;
    }
    if( state->rstate.stage==12 )
    {
        goto lbl_12;
    }
    if( state->rstate.stage==13 )
    {
        goto lbl_13;
    }
    if( state->rstate.stage==14 )
    {
        goto lbl_14;
    }
    if( state->rstate.stage==15 )
    {
        goto lbl_15;
    }
    if( state->rstate.stage==16 )
    {
        goto lbl_16;
    }
    if( state->rstate.stage==17 )
    {
        goto lbl_17;
    }
    if( state->rstate.stage==18 )
    {
        goto lbl_18;
    }
    if( state->rstate.stage==19 )
    {
        goto lbl_19;
    }
    
    /*
     * Routine body
     */
    
    /*
     * Prepare
     */
    n = state->n;
    state->repterminationtype = 0;
    state->repiterationscount = 0;
    state->repvaridx = -1;
    state->repnfev = 0;
    state->debugrestartscount = 0;
    
    /*
     *  Check, that transferred derivative value is right
     */
    mincg_clearrequestfields(state, _state);
    if( !(ae_fp_eq(state->diffstep,0)&&ae_fp_greater(state->teststep,0)) )
    {
        goto lbl_20;
    }
    state->needfg = ae_true;
    i = 0;
lbl_22:
    if( i>n-1 )
    {
        goto lbl_24;
    }
    v = state->x.ptr.p_double[i];
    state->x.ptr.p_double[i] = v-state->teststep*state->s.ptr.p_double[i];
    state->rstate.stage = 0;
    goto lbl_rcomm;
lbl_0:
    state->fm1 = state->f;
    state->fp1 = state->g.ptr.p_double[i];
    state->x.ptr.p_double[i] = v+state->teststep*state->s.ptr.p_double[i];
    state->rstate.stage = 1;
    goto lbl_rcomm;
lbl_1:
    state->fm2 = state->f;
    state->fp2 = state->g.ptr.p_double[i];
    state->x.ptr.p_double[i] = v;
    state->rstate.stage = 2;
    goto lbl_rcomm;
lbl_2:
    
    /*
     * 2*State.TestStep   -   scale parameter
     * width of segment [Xi-TestStep;Xi+TestStep]
     */
    if( !derivativecheck(state->fm1, state->fp1, state->fm2, state->fp2, state->f, state->g.ptr.p_double[i], 2*state->teststep, _state) )
    {
        state->repvaridx = i;
        state->repterminationtype = -7;
        result = ae_false;
        return result;
    }
    i = i+1;
    goto lbl_22;
lbl_24:
    state->needfg = ae_false;
lbl_20:
    
    /*
     * Preparations continue:
     * * set XK
     * * calculate F/G
     * * set DK to -G
     * * powerup algo (it may change preconditioner)
     * * apply preconditioner to DK
     * * report update of X
     * * check stopping conditions for G
     */
    ae_v_move(&state->xk.ptr.p_double[0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->terminationneeded = ae_false;
    mincg_clearrequestfields(state, _state);
    if( ae_fp_neq(state->diffstep,0) )
    {
        goto lbl_25;
    }
    state->needfg = ae_true;
    state->rstate.stage = 3;
    goto lbl_rcomm;
lbl_3:
    state->needfg = ae_false;
    goto lbl_26;
lbl_25:
    state->needf = ae_true;
    state->rstate.stage = 4;
    goto lbl_rcomm;
lbl_4:
    state->fbase = state->f;
    i = 0;
lbl_27:
    if( i>n-1 )
    {
        goto lbl_29;
    }
    v = state->x.ptr.p_double[i];
    state->x.ptr.p_double[i] = v-state->diffstep*state->s.ptr.p_double[i];
    state->rstate.stage = 5;
    goto lbl_rcomm;
lbl_5:
    state->fm2 = state->f;
    state->x.ptr.p_double[i] = v-0.5*state->diffstep*state->s.ptr.p_double[i];
    state->rstate.stage = 6;
    goto lbl_rcomm;
lbl_6:
    state->fm1 = state->f;
    state->x.ptr.p_double[i] = v+0.5*state->diffstep*state->s.ptr.p_double[i];
    state->rstate.stage = 7;
    goto lbl_rcomm;
lbl_7:
    state->fp1 = state->f;
    state->x.ptr.p_double[i] = v+state->diffstep*state->s.ptr.p_double[i];
    state->rstate.stage = 8;
    goto lbl_rcomm;
lbl_8:
    state->fp2 = state->f;
    state->x.ptr.p_double[i] = v;
    state->g.ptr.p_double[i] = (8*(state->fp1-state->fm1)-(state->fp2-state->fm2))/(6*state->diffstep*state->s.ptr.p_double[i]);
    i = i+1;
    goto lbl_27;
lbl_29:
    state->f = state->fbase;
    state->needf = ae_false;
lbl_26:
    if( !state->drep )
    {
        goto lbl_30;
    }
    
    /*
     * Report algorithm powerup (if needed)
     */
    mincg_clearrequestfields(state, _state);
    state->algpowerup = ae_true;
    state->rstate.stage = 9;
    goto lbl_rcomm;
lbl_9:
    state->algpowerup = ae_false;
lbl_30:
    trimprepare(state->f, &state->trimthreshold, _state);
    ae_v_moveneg(&state->dk.ptr.p_double[0], 1, &state->g.ptr.p_double[0], 1, ae_v_len(0,n-1));
    mincg_preconditionedmultiply(state, &state->dk, &state->work0, &state->work1, _state);
    if( !state->xrep )
    {
        goto lbl_32;
    }
    mincg_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 10;
    goto lbl_rcomm;
lbl_10:
    state->xupdated = ae_false;
lbl_32:
    if( state->terminationneeded )
    {
        ae_v_move(&state->xn.ptr.p_double[0], 1, &state->xk.ptr.p_double[0], 1, ae_v_len(0,n-1));
        state->repterminationtype = 8;
        result = ae_false;
        return result;
    }
    v = 0;
    for(i=0; i<=n-1; i++)
    {
        v = v+ae_sqr(state->g.ptr.p_double[i]*state->s.ptr.p_double[i], _state);
    }
    if( ae_fp_less_eq(ae_sqrt(v, _state),state->epsg) )
    {
        ae_v_move(&state->xn.ptr.p_double[0], 1, &state->xk.ptr.p_double[0], 1, ae_v_len(0,n-1));
        state->repterminationtype = 4;
        result = ae_false;
        return result;
    }
    state->repnfev = 1;
    state->k = 0;
    state->fold = state->f;
    
    /*
     * Choose initial step.
     * Apply preconditioner, if we have something other than default.
     */
    if( state->prectype==2||state->prectype==3 )
    {
        
        /*
         * because we use preconditioner, step length must be equal
         * to the norm of DK
         */
        v = ae_v_dotproduct(&state->dk.ptr.p_double[0], 1, &state->dk.ptr.p_double[0], 1, ae_v_len(0,n-1));
        state->lastgoodstep = ae_sqrt(v, _state);
    }
    else
    {
        
        /*
         * No preconditioner is used, we try to use suggested step
         */
        if( ae_fp_greater(state->suggestedstep,0) )
        {
            state->lastgoodstep = state->suggestedstep;
        }
        else
        {
            state->lastgoodstep = 1.0;
        }
    }
    
    /*
     * Main cycle
     */
    state->rstimer = mincg_rscountdownlen;
lbl_34:
    if( ae_false )
    {
        goto lbl_35;
    }
    
    /*
     * * clear reset flag
     * * clear termination flag
     * * store G[k] for later calculation of Y[k]
     * * prepare starting point and direction and step length for line search
     */
    state->innerresetneeded = ae_false;
    state->terminationneeded = ae_false;
    ae_v_moveneg(&state->yk.ptr.p_double[0], 1, &state->g.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_move(&state->d.ptr.p_double[0], 1, &state->dk.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xk.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->mcstage = 0;
    state->stp = 1.0;
    linminnormalized(&state->d, &state->stp, n, _state);
    if( ae_fp_neq(state->lastgoodstep,0) )
    {
        state->stp = state->lastgoodstep;
    }
    state->curstpmax = state->stpmax;
    
    /*
     * Report beginning of line search (if needed)
     * Terminate algorithm, if user request was detected
     */
    if( !state->drep )
    {
        goto lbl_36;
    }
    mincg_clearrequestfields(state, _state);
    state->lsstart = ae_true;
    state->rstate.stage = 11;
    goto lbl_rcomm;
lbl_11:
    state->lsstart = ae_false;
lbl_36:
    if( state->terminationneeded )
    {
        ae_v_move(&state->xn.ptr.p_double[0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,n-1));
        state->repterminationtype = 8;
        result = ae_false;
        return result;
    }
    
    /*
     * Minimization along D
     */
    mcsrch(n, &state->x, &state->f, &state->g, &state->d, &state->stp, state->curstpmax, mincg_gtol, &state->mcinfo, &state->nfev, &state->work0, &state->lstate, &state->mcstage, _state);
lbl_38:
    if( state->mcstage==0 )
    {
        goto lbl_39;
    }
    
    /*
     * Calculate function/gradient using either
     * analytical gradient supplied by user
     * or finite difference approximation.
     *
     * "Trim" function in order to handle near-singularity points.
     */
    mincg_clearrequestfields(state, _state);
    if( ae_fp_neq(state->diffstep,0) )
    {
        goto lbl_40;
    }
    state->needfg = ae_true;
    state->rstate.stage = 12;
    goto lbl_rcomm;
lbl_12:
    state->needfg = ae_false;
    goto lbl_41;
lbl_40:
    state->needf = ae_true;
    state->rstate.stage = 13;
    goto lbl_rcomm;
lbl_13:
    state->fbase = state->f;
    i = 0;
lbl_42:
    if( i>n-1 )
    {
        goto lbl_44;
    }
    v = state->x.ptr.p_double[i];
    state->x.ptr.p_double[i] = v-state->diffstep*state->s.ptr.p_double[i];
    state->rstate.stage = 14;
    goto lbl_rcomm;
lbl_14:
    state->fm2 = state->f;
    state->x.ptr.p_double[i] = v-0.5*state->diffstep*state->s.ptr.p_double[i];
    state->rstate.stage = 15;
    goto lbl_rcomm;
lbl_15:
    state->fm1 = state->f;
    state->x.ptr.p_double[i] = v+0.5*state->diffstep*state->s.ptr.p_double[i];
    state->rstate.stage = 16;
    goto lbl_rcomm;
lbl_16:
    state->fp1 = state->f;
    state->x.ptr.p_double[i] = v+state->diffstep*state->s.ptr.p_double[i];
    state->rstate.stage = 17;
    goto lbl_rcomm;
lbl_17:
    state->fp2 = state->f;
    state->x.ptr.p_double[i] = v;
    state->g.ptr.p_double[i] = (8*(state->fp1-state->fm1)-(state->fp2-state->fm2))/(6*state->diffstep*state->s.ptr.p_double[i]);
    i = i+1;
    goto lbl_42;
lbl_44:
    state->f = state->fbase;
    state->needf = ae_false;
lbl_41:
    trimfunction(&state->f, &state->g, n, state->trimthreshold, _state);
    
    /*
     * Call MCSRCH again
     */
    mcsrch(n, &state->x, &state->f, &state->g, &state->d, &state->stp, state->curstpmax, mincg_gtol, &state->mcinfo, &state->nfev, &state->work0, &state->lstate, &state->mcstage, _state);
    goto lbl_38;
lbl_39:
    
    /*
     * * report end of line search
     * * store current point to XN
     * * report iteration
     * * terminate algorithm if user request was detected
     */
    if( !state->drep )
    {
        goto lbl_45;
    }
    
    /*
     * Report end of line search (if needed)
     */
    mincg_clearrequestfields(state, _state);
    state->lsend = ae_true;
    state->rstate.stage = 18;
    goto lbl_rcomm;
lbl_18:
    state->lsend = ae_false;
lbl_45:
    ae_v_move(&state->xn.ptr.p_double[0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,n-1));
    if( !state->xrep )
    {
        goto lbl_47;
    }
    mincg_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 19;
    goto lbl_rcomm;
lbl_19:
    state->xupdated = ae_false;
lbl_47:
    if( state->terminationneeded )
    {
        ae_v_move(&state->xn.ptr.p_double[0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,n-1));
        state->repterminationtype = 8;
        result = ae_false;
        return result;
    }
    
    /*
     * Line search is finished.
     * * calculate BetaK
     * * calculate DN
     * * update timers
     * * calculate step length:
     *   * LastScaledStep is ALWAYS calculated because it is used in the stopping criteria
     *   * LastGoodStep is updated only when MCINFO is equal to 1 (Wolfe conditions hold).
     *     See below for more explanation.
     */
    if( state->mcinfo==1&&!state->innerresetneeded )
    {
        
        /*
         * Standard Wolfe conditions hold
         * Calculate Y[K] and D[K]'*Y[K]
         */
        ae_v_add(&state->yk.ptr.p_double[0], 1, &state->g.ptr.p_double[0], 1, ae_v_len(0,n-1));
        vv = ae_v_dotproduct(&state->yk.ptr.p_double[0], 1, &state->dk.ptr.p_double[0], 1, ae_v_len(0,n-1));
        
        /*
         * Calculate BetaK according to DY formula
         */
        v = mincg_preconditionedmultiply2(state, &state->g, &state->g, &state->work0, &state->work1, _state);
        state->betady = v/vv;
        
        /*
         * Calculate BetaK according to HS formula
         */
        v = mincg_preconditionedmultiply2(state, &state->g, &state->yk, &state->work0, &state->work1, _state);
        state->betahs = v/vv;
        
        /*
         * Choose BetaK
         */
        if( state->cgtype==0 )
        {
            betak = state->betady;
        }
        if( state->cgtype==1 )
        {
            betak = ae_maxreal(0, ae_minreal(state->betady, state->betahs, _state), _state);
        }
    }
    else
    {
        
        /*
         * Something is wrong (may be function is too wild or too flat)
         * or we just have to restart algo.
         *
         * We'll set BetaK=0, which will restart CG algorithm.
         * We can stop later (during normal checks) if stopping conditions are met.
         */
        betak = 0;
        state->debugrestartscount = state->debugrestartscount+1;
    }
    if( state->repiterationscount>0&&state->repiterationscount%(3+n)==0 )
    {
        
        /*
         * clear Beta every N iterations
         */
        betak = 0;
    }
    if( state->mcinfo==1||state->mcinfo==5 )
    {
        state->rstimer = mincg_rscountdownlen;
    }
    else
    {
        state->rstimer = state->rstimer-1;
    }
    ae_v_moveneg(&state->dn.ptr.p_double[0], 1, &state->g.ptr.p_double[0], 1, ae_v_len(0,n-1));
    mincg_preconditionedmultiply(state, &state->dn, &state->work0, &state->work1, _state);
    ae_v_addd(&state->dn.ptr.p_double[0], 1, &state->dk.ptr.p_double[0], 1, ae_v_len(0,n-1), betak);
    state->lastscaledstep = 0.0;
    for(i=0; i<=n-1; i++)
    {
        state->lastscaledstep = state->lastscaledstep+ae_sqr(state->d.ptr.p_double[i]/state->s.ptr.p_double[i], _state);
    }
    state->lastscaledstep = state->stp*ae_sqrt(state->lastscaledstep, _state);
    if( state->mcinfo==1 )
    {
        
        /*
         * Step is good (Wolfe conditions hold), update LastGoodStep.
         *
         * This check for MCINFO=1 is essential because sometimes in the
         * constrained optimization setting we may take very short steps
         * (like 1E-15) because we were very close to boundary of the
         * feasible area. Such short step does not mean that we've converged
         * to the solution - it was so short because we were close to the
         * boundary and there was a limit on step length.
         *
         * So having such short step is quite normal situation. However, we
         * should NOT start next iteration from step whose initial length is
         * estimated as 1E-15 because it may lead to the failure of the
         * linear minimizer (step is too short, function does not changes,
         * line search stagnates).
         */
        state->lastgoodstep = 0;
        for(i=0; i<=n-1; i++)
        {
            state->lastgoodstep = state->lastgoodstep+ae_sqr(state->d.ptr.p_double[i], _state);
        }
        state->lastgoodstep = state->stp*ae_sqrt(state->lastgoodstep, _state);
    }
    
    /*
     * Update information.
     * Check stopping conditions.
     */
    state->repnfev = state->repnfev+state->nfev;
    state->repiterationscount = state->repiterationscount+1;
    if( state->repiterationscount>=state->maxits&&state->maxits>0 )
    {
        
        /*
         * Too many iterations
         */
        state->repterminationtype = 5;
        result = ae_false;
        return result;
    }
    v = 0;
    for(i=0; i<=n-1; i++)
    {
        v = v+ae_sqr(state->g.ptr.p_double[i]*state->s.ptr.p_double[i], _state);
    }
    if( ae_fp_less_eq(ae_sqrt(v, _state),state->epsg) )
    {
        
        /*
         * Gradient is small enough
         */
        state->repterminationtype = 4;
        result = ae_false;
        return result;
    }
    if( !state->innerresetneeded )
    {
        
        /*
         * These conditions are checked only when no inner reset was requested by user
         */
        if( ae_fp_less_eq(state->fold-state->f,state->epsf*ae_maxreal(ae_fabs(state->fold, _state), ae_maxreal(ae_fabs(state->f, _state), 1.0, _state), _state)) )
        {
            
            /*
             * F(k+1)-F(k) is small enough
             */
            state->repterminationtype = 1;
            result = ae_false;
            return result;
        }
        if( ae_fp_less_eq(state->lastscaledstep,state->epsx) )
        {
            
            /*
             * X(k+1)-X(k) is small enough
             */
            state->repterminationtype = 2;
            result = ae_false;
            return result;
        }
    }
    if( state->rstimer<=0 )
    {
        
        /*
         * Too many subsequent restarts
         */
        state->repterminationtype = 7;
        result = ae_false;
        return result;
    }
    
    /*
     * Shift Xk/Dk, update other information
     */
    ae_v_move(&state->xk.ptr.p_double[0], 1, &state->xn.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_move(&state->dk.ptr.p_double[0], 1, &state->dn.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->fold = state->f;
    state->k = state->k+1;
    goto lbl_34;
lbl_35:
    result = ae_false;
    return result;
    
    /*
     * Saving state
     */
lbl_rcomm:
    result = ae_true;
    state->rstate.ia.ptr.p_int[0] = n;
    state->rstate.ia.ptr.p_int[1] = i;
    state->rstate.ra.ptr.p_double[0] = betak;
    state->rstate.ra.ptr.p_double[1] = v;
    state->rstate.ra.ptr.p_double[2] = vv;
    return result;
}


/*************************************************************************
Conjugate gradient results

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    X       -   array[0..N-1], solution
    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    * -7    gradient verification failed.
                            See MinCGSetGradientCheck() for more information.
                    *  1    relative function improvement is no more than
                            EpsF.
                    *  2    relative step is no more than EpsX.
                    *  4    gradient norm is no more than EpsG
                    *  5    MaxIts steps was taken
                    *  7    stopping conditions are too stringent,
                            further improvement is impossible,
                            we return best X found so far
                    *  8    terminated by user
                * Rep.IterationsCount contains iterations count
                * NFEV countains number of function calculations

  -- ALGLIB --
     Copyright 20.04.2009 by Bochkanov Sergey
*************************************************************************/
void mincgresults(mincgstate* state,
     /* Real    */ ae_vector* x,
     mincgreport* rep,
     ae_state *_state)
{

    ae_vector_clear(x);
    _mincgreport_clear(rep);

    mincgresultsbuf(state, x, rep, _state);
}


/*************************************************************************
Conjugate gradient results

Buffered implementation of MinCGResults(), which uses pre-allocated buffer
to store X[]. If buffer size is  too  small,  it  resizes  buffer.  It  is
intended to be used in the inner cycles of performance critical algorithms
where array reallocation penalty is too large to be ignored.

  -- ALGLIB --
     Copyright 20.04.2009 by Bochkanov Sergey
*************************************************************************/
void mincgresultsbuf(mincgstate* state,
     /* Real    */ ae_vector* x,
     mincgreport* rep,
     ae_state *_state)
{


    if( x->cnt<state->n )
    {
        ae_vector_set_length(x, state->n, _state);
    }
    ae_v_move(&x->ptr.p_double[0], 1, &state->xn.ptr.p_double[0], 1, ae_v_len(0,state->n-1));
    rep->iterationscount = state->repiterationscount;
    rep->nfev = state->repnfev;
    rep->varidx = state->repvaridx;
    rep->terminationtype = state->repterminationtype;
}


/*************************************************************************
This  subroutine  restarts  CG  algorithm from new point. All optimization
parameters are left unchanged.

This  function  allows  to  solve multiple  optimization  problems  (which
must have same number of dimensions) without object reallocation penalty.

INPUT PARAMETERS:
    State   -   structure used to store algorithm state.
    X       -   new starting point.

  -- ALGLIB --
     Copyright 30.07.2010 by Bochkanov Sergey
*************************************************************************/
void mincgrestartfrom(mincgstate* state,
     /* Real    */ ae_vector* x,
     ae_state *_state)
{


    ae_assert(x->cnt>=state->n, "MinCGRestartFrom: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, state->n, _state), "MinCGCreate: X contains infinite or NaN values!", _state);
    ae_v_move(&state->x.ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,state->n-1));
    mincgsuggeststep(state, 0.0, _state);
    ae_vector_set_length(&state->rstate.ia, 1+1, _state);
    ae_vector_set_length(&state->rstate.ra, 2+1, _state);
    state->rstate.stage = -1;
    mincg_clearrequestfields(state, _state);
}


/*************************************************************************
Faster version of MinCGSetPrecDiag(), for time-critical parts of code,
without safety checks.

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsetprecdiagfast(mincgstate* state,
     /* Real    */ ae_vector* d,
     ae_state *_state)
{
    ae_int_t i;


    rvectorsetlengthatleast(&state->diagh, state->n, _state);
    rvectorsetlengthatleast(&state->diaghl2, state->n, _state);
    state->prectype = 2;
    state->vcnt = 0;
    state->innerresetneeded = ae_true;
    for(i=0; i<=state->n-1; i++)
    {
        state->diagh.ptr.p_double[i] = d->ptr.p_double[i];
        state->diaghl2.ptr.p_double[i] = 0.0;
    }
}


/*************************************************************************
This function sets low-rank preconditioner for Hessian matrix  H=D+V'*C*V,
where:
* H is a Hessian matrix, which is approximated by D/V/C
* D=D1+D2 is a diagonal matrix, which includes two positive definite terms:
  * constant term D1 (is not updated or infrequently updated)
  * variable term D2 (can be cheaply updated from iteration to iteration)
* V is a low-rank correction
* C is a diagonal factor of low-rank correction

Preconditioner P is calculated using approximate Woodburry formula:
    P  = D^(-1) - D^(-1)*V'*(C^(-1)+V*D1^(-1)*V')^(-1)*V*D^(-1)
       = D^(-1) - D^(-1)*VC'*VC*D^(-1),
where
    VC = sqrt(B)*V
    B  = (C^(-1)+V*D1^(-1)*V')^(-1)
    
Note that B is calculated using constant term (D1) only,  which  allows us
to update D2 without recalculation of B or   VC.  Such  preconditioner  is
exact when D2 is zero. When D2 is non-zero, it is only approximation,  but
very good and cheap one.

This function accepts D1, V, C.
D2 is set to zero by default.

Cost of this update is O(N*VCnt*VCnt), but D2 can be updated in just O(N)
by MinCGSetPrecVarPart.

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsetpreclowrankfast(mincgstate* state,
     /* Real    */ ae_vector* d1,
     /* Real    */ ae_vector* c,
     /* Real    */ ae_matrix* v,
     ae_int_t vcnt,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t n;
    double t;
    ae_matrix b;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init(&b, 0, 0, DT_REAL, _state, ae_true);

    if( vcnt==0 )
    {
        mincgsetprecdiagfast(state, d1, _state);
        ae_frame_leave(_state);
        return;
    }
    n = state->n;
    ae_matrix_set_length(&b, vcnt, vcnt, _state);
    rvectorsetlengthatleast(&state->diagh, n, _state);
    rvectorsetlengthatleast(&state->diaghl2, n, _state);
    rmatrixsetlengthatleast(&state->vcorr, vcnt, n, _state);
    state->prectype = 2;
    state->vcnt = vcnt;
    state->innerresetneeded = ae_true;
    for(i=0; i<=n-1; i++)
    {
        state->diagh.ptr.p_double[i] = d1->ptr.p_double[i];
        state->diaghl2.ptr.p_double[i] = 0.0;
    }
    for(i=0; i<=vcnt-1; i++)
    {
        for(j=i; j<=vcnt-1; j++)
        {
            t = 0;
            for(k=0; k<=n-1; k++)
            {
                t = t+v->ptr.pp_double[i][k]*v->ptr.pp_double[j][k]/d1->ptr.p_double[k];
            }
            b.ptr.pp_double[i][j] = t;
        }
        b.ptr.pp_double[i][i] = b.ptr.pp_double[i][i]+1.0/c->ptr.p_double[i];
    }
    if( !spdmatrixcholeskyrec(&b, 0, vcnt, ae_true, &state->work0, _state) )
    {
        state->vcnt = 0;
        ae_frame_leave(_state);
        return;
    }
    for(i=0; i<=vcnt-1; i++)
    {
        ae_v_move(&state->vcorr.ptr.pp_double[i][0], 1, &v->ptr.pp_double[i][0], 1, ae_v_len(0,n-1));
        for(j=0; j<=i-1; j++)
        {
            t = b.ptr.pp_double[j][i];
            ae_v_subd(&state->vcorr.ptr.pp_double[i][0], 1, &state->vcorr.ptr.pp_double[j][0], 1, ae_v_len(0,n-1), t);
        }
        t = 1/b.ptr.pp_double[i][i];
        ae_v_muld(&state->vcorr.ptr.pp_double[i][0], 1, ae_v_len(0,n-1), t);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This function updates variable part (diagonal matrix D2)
of low-rank preconditioner.

This update is very cheap and takes just O(N) time.

It has no effect with default preconditioner.

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void mincgsetprecvarpart(mincgstate* state,
     /* Real    */ ae_vector* d2,
     ae_state *) //_state)
{
    ae_int_t i;
    ae_int_t n;


    n = state->n;
    for(i=0; i<=n-1; i++)
    {
        state->diaghl2.ptr.p_double[i] = d2->ptr.p_double[i];
    }
}


/*************************************************************************

This  subroutine  turns  on  verification  of  the  user-supplied analytic
gradient:
* user calls this subroutine before optimization begins
* MinCGOptimize() is called
* prior to  actual  optimization, for each component  of  parameters being
  optimized X[i] algorithm performs following steps:
  * two trial steps are made to X[i]-TestStep*S[i] and X[i]+TestStep*S[i],
    where X[i] is i-th component of the initial point and S[i] is a  scale
    of i-th parameter
  * F(X) is evaluated at these trial points
  * we perform one more evaluation in the middle point of the interval
  * we  build  cubic  model using function values and derivatives at trial
    points and we compare its prediction with actual value in  the  middle
    point
  * in case difference between prediction and actual value is higher  than
    some predetermined threshold, algorithm stops with completion code -7;
    Rep.VarIdx is set to index of the parameter with incorrect derivative.
* after verification is over, algorithm proceeds to the actual optimization.

NOTE 1: verification  needs  N (parameters count) gradient evaluations. It
        is very costly and you should use  it  only  for  low  dimensional
        problems,  when  you  want  to  be  sure  that  you've   correctly
        calculated  analytic  derivatives.  You  should  not use it in the
        production code (unless you want to check derivatives provided  by
        some third party).

NOTE 2: you  should  carefully  choose  TestStep. Value which is too large
        (so large that function behaviour is significantly non-cubic) will
        lead to false alarms. You may use  different  step  for  different
        parameters by means of setting scale with MinCGSetScale().

NOTE 3: this function may lead to false positives. In case it reports that
        I-th  derivative was calculated incorrectly, you may decrease test
        step  and  try  one  more  time  - maybe your function changes too
        sharply  and  your  step  is  too  large for such rapidly chanding
        function.

INPUT PARAMETERS:
    State       -   structure used to store algorithm state
    TestStep    -   verification step:
                    * TestStep=0 turns verification off
                    * TestStep>0 activates verification

  -- ALGLIB --
     Copyright 31.05.2012 by Bochkanov Sergey
*************************************************************************/
void mincgsetgradientcheck(mincgstate* state,
     double teststep,
     ae_state *_state)
{


    ae_assert(ae_isfinite(teststep, _state), "MinCGSetGradientCheck: TestStep contains NaN or Infinite", _state);
    ae_assert(ae_fp_greater_eq(teststep,0), "MinCGSetGradientCheck: invalid argument TestStep(TestStep<0)", _state);
    state->teststep = teststep;
}


/*************************************************************************
Clears request fileds (to be sure that we don't forgot to clear something)
*************************************************************************/
static void mincg_clearrequestfields(mincgstate* state, ae_state *) //_state)
{


    state->needf = ae_false;
    state->needfg = ae_false;
    state->xupdated = ae_false;
    state->lsstart = ae_false;
    state->lsend = ae_false;
    state->algpowerup = ae_false;
}


/*************************************************************************
This function calculates preconditioned product H^(-1)*x and stores result
back into X. Work0[] and Work1[] are used as temporaries (size must be at
least N; this function doesn't allocate arrays).

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
static void mincg_preconditionedmultiply(mincgstate* state,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* work0,
     /* Real    */ ae_vector* work1,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t n;
    ae_int_t vcnt;
    double v;


    n = state->n;
    vcnt = state->vcnt;
    if( state->prectype==0 )
    {
        return;
    }
    if( state->prectype==3 )
    {
        for(i=0; i<=n-1; i++)
        {
            x->ptr.p_double[i] = x->ptr.p_double[i]*state->s.ptr.p_double[i]*state->s.ptr.p_double[i];
        }
        return;
    }
    ae_assert(state->prectype==2, "MinCG: internal error (unexpected PrecType)", _state);
    
    /*
     * handle part common for VCnt=0 and VCnt<>0
     */
    for(i=0; i<=n-1; i++)
    {
        x->ptr.p_double[i] = x->ptr.p_double[i]/(state->diagh.ptr.p_double[i]+state->diaghl2.ptr.p_double[i]);
    }
    
    /*
     * if VCnt>0
     */
    if( vcnt>0 )
    {
        for(i=0; i<=vcnt-1; i++)
        {
            v = ae_v_dotproduct(&state->vcorr.ptr.pp_double[i][0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,n-1));
            work0->ptr.p_double[i] = v;
        }
        for(i=0; i<=n-1; i++)
        {
            work1->ptr.p_double[i] = 0;
        }
        for(i=0; i<=vcnt-1; i++)
        {
            v = work0->ptr.p_double[i];
            ae_v_addd(&state->work1.ptr.p_double[0], 1, &state->vcorr.ptr.pp_double[i][0], 1, ae_v_len(0,n-1), v);
        }
        for(i=0; i<=n-1; i++)
        {
            x->ptr.p_double[i] = x->ptr.p_double[i]-state->work1.ptr.p_double[i]/(state->diagh.ptr.p_double[i]+state->diaghl2.ptr.p_double[i]);
        }
    }
}


/*************************************************************************
This function calculates preconditioned product x'*H^(-1)*y. Work0[] and
Work1[] are used as temporaries (size must be at least N; this function
doesn't allocate arrays).

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
static double mincg_preconditionedmultiply2(mincgstate* state,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* work0,
     /* Real    */ ae_vector* work1,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t n;
    ae_int_t vcnt;
    double v0;
    double v1;
    double result;


    n = state->n;
    vcnt = state->vcnt;
    
    /*
     * no preconditioning
     */
    if( state->prectype==0 )
    {
        v0 = ae_v_dotproduct(&x->ptr.p_double[0], 1, &y->ptr.p_double[0], 1, ae_v_len(0,n-1));
        result = v0;
        return result;
    }
    if( state->prectype==3 )
    {
        result = 0;
        for(i=0; i<=n-1; i++)
        {
            result = result+x->ptr.p_double[i]*state->s.ptr.p_double[i]*state->s.ptr.p_double[i]*y->ptr.p_double[i];
        }
        return result;
    }
    ae_assert(state->prectype==2, "MinCG: internal error (unexpected PrecType)", _state);
    
    /*
     * low rank preconditioning
     */
    result = 0.0;
    for(i=0; i<=n-1; i++)
    {
        result = result+x->ptr.p_double[i]*y->ptr.p_double[i]/(state->diagh.ptr.p_double[i]+state->diaghl2.ptr.p_double[i]);
    }
    if( vcnt>0 )
    {
        for(i=0; i<=n-1; i++)
        {
            work0->ptr.p_double[i] = x->ptr.p_double[i]/(state->diagh.ptr.p_double[i]+state->diaghl2.ptr.p_double[i]);
            work1->ptr.p_double[i] = y->ptr.p_double[i]/(state->diagh.ptr.p_double[i]+state->diaghl2.ptr.p_double[i]);
        }
        for(i=0; i<=vcnt-1; i++)
        {
            v0 = ae_v_dotproduct(&work0->ptr.p_double[0], 1, &state->vcorr.ptr.pp_double[i][0], 1, ae_v_len(0,n-1));
            v1 = ae_v_dotproduct(&work1->ptr.p_double[0], 1, &state->vcorr.ptr.pp_double[i][0], 1, ae_v_len(0,n-1));
            result = result-v0*v1;
        }
    }
    return result;
}


/*************************************************************************
Internal initialization subroutine

  -- ALGLIB --
     Copyright 16.05.2011 by Bochkanov Sergey
*************************************************************************/
static void mincg_mincginitinternal(ae_int_t n,
     double diffstep,
     mincgstate* state,
     ae_state *_state)
{
    ae_int_t i;


    
    /*
     * Initialize
     */
    state->teststep = 0;
    state->n = n;
    state->diffstep = diffstep;
    mincgsetcond(state, 0, 0, 0, 0, _state);
    mincgsetxrep(state, ae_false, _state);
    mincgsetdrep(state, ae_false, _state);
    mincgsetstpmax(state, 0, _state);
    mincgsetcgtype(state, -1, _state);
    mincgsetprecdefault(state, _state);
    ae_vector_set_length(&state->xk, n, _state);
    ae_vector_set_length(&state->dk, n, _state);
    ae_vector_set_length(&state->xn, n, _state);
    ae_vector_set_length(&state->dn, n, _state);
    ae_vector_set_length(&state->x, n, _state);
    ae_vector_set_length(&state->d, n, _state);
    ae_vector_set_length(&state->g, n, _state);
    ae_vector_set_length(&state->work0, n, _state);
    ae_vector_set_length(&state->work1, n, _state);
    ae_vector_set_length(&state->yk, n, _state);
    ae_vector_set_length(&state->s, n, _state);
    for(i=0; i<=n-1; i++)
    {
        state->s.ptr.p_double[i] = 1.0;
    }
}


ae_bool _mincgstate_init(mincgstate* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->diagh, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->diaghl2, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->vcorr, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->s, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xk, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->dk, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xn, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->dn, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->d, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->yk, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->x, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->g, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !_rcommstate_init(&p->rstate, _state, make_automatic) )
        return ae_false;
    if( !_linminstate_init(&p->lstate, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->work0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->work1, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _mincgstate_init_copy(mincgstate* dst, mincgstate* src, ae_state *_state, ae_bool make_automatic)
{
    dst->n = src->n;
    dst->epsg = src->epsg;
    dst->epsf = src->epsf;
    dst->epsx = src->epsx;
    dst->maxits = src->maxits;
    dst->stpmax = src->stpmax;
    dst->suggestedstep = src->suggestedstep;
    dst->xrep = src->xrep;
    dst->drep = src->drep;
    dst->cgtype = src->cgtype;
    dst->prectype = src->prectype;
    if( !ae_vector_init_copy(&dst->diagh, &src->diagh, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->diaghl2, &src->diaghl2, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->vcorr, &src->vcorr, _state, make_automatic) )
        return ae_false;
    dst->vcnt = src->vcnt;
    if( !ae_vector_init_copy(&dst->s, &src->s, _state, make_automatic) )
        return ae_false;
    dst->diffstep = src->diffstep;
    dst->nfev = src->nfev;
    dst->mcstage = src->mcstage;
    dst->k = src->k;
    if( !ae_vector_init_copy(&dst->xk, &src->xk, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->dk, &src->dk, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->xn, &src->xn, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->dn, &src->dn, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->d, &src->d, _state, make_automatic) )
        return ae_false;
    dst->fold = src->fold;
    dst->stp = src->stp;
    dst->curstpmax = src->curstpmax;
    if( !ae_vector_init_copy(&dst->yk, &src->yk, _state, make_automatic) )
        return ae_false;
    dst->lastgoodstep = src->lastgoodstep;
    dst->lastscaledstep = src->lastscaledstep;
    dst->mcinfo = src->mcinfo;
    dst->innerresetneeded = src->innerresetneeded;
    dst->terminationneeded = src->terminationneeded;
    dst->trimthreshold = src->trimthreshold;
    dst->rstimer = src->rstimer;
    if( !ae_vector_init_copy(&dst->x, &src->x, _state, make_automatic) )
        return ae_false;
    dst->f = src->f;
    if( !ae_vector_init_copy(&dst->g, &src->g, _state, make_automatic) )
        return ae_false;
    dst->needf = src->needf;
    dst->needfg = src->needfg;
    dst->xupdated = src->xupdated;
    dst->algpowerup = src->algpowerup;
    dst->lsstart = src->lsstart;
    dst->lsend = src->lsend;
    dst->teststep = src->teststep;
    if( !_rcommstate_init_copy(&dst->rstate, &src->rstate, _state, make_automatic) )
        return ae_false;
    dst->repiterationscount = src->repiterationscount;
    dst->repnfev = src->repnfev;
    dst->repvaridx = src->repvaridx;
    dst->repterminationtype = src->repterminationtype;
    dst->debugrestartscount = src->debugrestartscount;
    if( !_linminstate_init_copy(&dst->lstate, &src->lstate, _state, make_automatic) )
        return ae_false;
    dst->fbase = src->fbase;
    dst->fm2 = src->fm2;
    dst->fm1 = src->fm1;
    dst->fp1 = src->fp1;
    dst->fp2 = src->fp2;
    dst->betahs = src->betahs;
    dst->betady = src->betady;
    if( !ae_vector_init_copy(&dst->work0, &src->work0, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->work1, &src->work1, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _mincgstate_clear(mincgstate* p)
{
    ae_vector_clear(&p->diagh);
    ae_vector_clear(&p->diaghl2);
    ae_matrix_clear(&p->vcorr);
    ae_vector_clear(&p->s);
    ae_vector_clear(&p->xk);
    ae_vector_clear(&p->dk);
    ae_vector_clear(&p->xn);
    ae_vector_clear(&p->dn);
    ae_vector_clear(&p->d);
    ae_vector_clear(&p->yk);
    ae_vector_clear(&p->x);
    ae_vector_clear(&p->g);
    _rcommstate_clear(&p->rstate);
    _linminstate_clear(&p->lstate);
    ae_vector_clear(&p->work0);
    ae_vector_clear(&p->work1);
}


ae_bool _mincgreport_init(mincgreport* /*p*/, ae_state */*_state*/, ae_bool) // make_automatic)
{
    return ae_true;
}


ae_bool _mincgreport_init_copy(mincgreport* dst, mincgreport* src, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    dst->iterationscount = src->iterationscount;
    dst->nfev = src->nfev;
    dst->varidx = src->varidx;
    dst->terminationtype = src->terminationtype;
    return ae_true;
}


void _mincgreport_clear(mincgreport*) // p)
{
}




/*************************************************************************
                     BOUND CONSTRAINED OPTIMIZATION
       WITH ADDITIONAL LINEAR EQUALITY AND INEQUALITY CONSTRAINTS

DESCRIPTION:
The  subroutine  minimizes  function   F(x)  of N arguments subject to any
combination of:
* bound constraints
* linear inequality constraints
* linear equality constraints

REQUIREMENTS:
* user must provide function value and gradient
* starting point X0 must be feasible or
  not too far away from the feasible set
* grad(f) must be Lipschitz continuous on a level set:
  L = { x : f(x)<=f(x0) }
* function must be defined everywhere on the feasible set F

USAGE:

Constrained optimization if far more complex than the unconstrained one.
Here we give very brief outline of the BLEIC optimizer. We strongly recommend
you to read examples in the ALGLIB Reference Manual and to read ALGLIB User Guide
on optimization, which is available at http://www.alglib.net/optimization/

1. User initializes algorithm state with MinBLEICCreate() call

2. USer adds boundary and/or linear constraints by calling
   MinBLEICSetBC() and MinBLEICSetLC() functions.

3. User sets stopping conditions for underlying unconstrained solver
   with MinBLEICSetInnerCond() call.
   This function controls accuracy of underlying optimization algorithm.

4. User sets stopping conditions for outer iteration by calling
   MinBLEICSetOuterCond() function.
   This function controls handling of boundary and inequality constraints.

5. Additionally, user may set limit on number of internal iterations
   by MinBLEICSetMaxIts() call.
   This function allows to prevent algorithm from looping forever.

6. User calls MinBLEICOptimize() function which takes algorithm  state and
   pointer (delegate, etc.) to callback function which calculates F/G.

7. User calls MinBLEICResults() to get solution

8. Optionally user may call MinBLEICRestartFrom() to solve another problem
   with same N but another starting point.
   MinBLEICRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   problem dimension, N>0:
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size ofX
    X       -   starting point, array[N]:
                * it is better to set X to a feasible point
                * but X can be infeasible, in which case algorithm will try
                  to find feasible point first, using X as initial
                  approximation.

OUTPUT PARAMETERS:
    State   -   structure stores algorithm state

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleiccreate(ae_int_t n,
     /* Real    */ ae_vector* x,
     minbleicstate* state,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix c;
    ae_vector ct;

    ae_frame_make(_state, &_frame_block);
    _minbleicstate_clear(state);
    ae_matrix_init(&c, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&ct, 0, DT_INT, _state, ae_true);

    ae_assert(n>=1, "MinBLEICCreate: N<1", _state);
    ae_assert(x->cnt>=n, "MinBLEICCreate: Length(X)<N", _state);
    ae_assert(isfinitevector(x, n, _state), "MinBLEICCreate: X contains infinite or NaN values!", _state);
    minbleic_minbleicinitinternal(n, x, 0.0, state, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
The subroutine is finite difference variant of MinBLEICCreate().  It  uses
finite differences in order to differentiate target function.

Description below contains information which is specific to  this function
only. We recommend to read comments on MinBLEICCreate() in  order  to  get
more information about creation of BLEIC optimizer.

INPUT PARAMETERS:
    N       -   problem dimension, N>0:
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    X       -   starting point, array[0..N-1].
    DiffStep-   differentiation step, >0

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTES:
1. algorithm uses 4-point central formula for differentiation.
2. differentiation step along I-th axis is equal to DiffStep*S[I] where
   S[] is scaling vector which can be set by MinBLEICSetScale() call.
3. we recommend you to use moderate values of  differentiation  step.  Too
   large step will result in too large truncation  errors, while too small
   step will result in too large numerical  errors.  1.0E-6  can  be  good
   value to start with.
4. Numerical  differentiation  is   very   inefficient  -   one   gradient
   calculation needs 4*N function evaluations. This function will work for
   any N - either small (1...10), moderate (10...100) or  large  (100...).
   However, performance penalty will be too severe for any N's except  for
   small ones.
   We should also say that code which relies on numerical  differentiation
   is  less  robust and precise. CG needs exact gradient values. Imprecise
   gradient may slow  down  convergence, especially  on  highly  nonlinear
   problems.
   Thus  we  recommend to use this function for fast prototyping on small-
   dimensional problems only, and to implement analytical gradient as soon
   as possible.

  -- ALGLIB --
     Copyright 16.05.2011 by Bochkanov Sergey
*************************************************************************/
void minbleiccreatef(ae_int_t n,
     /* Real    */ ae_vector* x,
     double diffstep,
     minbleicstate* state,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix c;
    ae_vector ct;

    ae_frame_make(_state, &_frame_block);
    _minbleicstate_clear(state);
    ae_matrix_init(&c, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&ct, 0, DT_INT, _state, ae_true);

    ae_assert(n>=1, "MinBLEICCreateF: N<1", _state);
    ae_assert(x->cnt>=n, "MinBLEICCreateF: Length(X)<N", _state);
    ae_assert(isfinitevector(x, n, _state), "MinBLEICCreateF: X contains infinite or NaN values!", _state);
    ae_assert(ae_isfinite(diffstep, _state), "MinBLEICCreateF: DiffStep is infinite or NaN!", _state);
    ae_assert(ae_fp_greater(diffstep,0), "MinBLEICCreateF: DiffStep is non-positive!", _state);
    minbleic_minbleicinitinternal(n, x, diffstep, state, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
This function sets boundary constraints for BLEIC optimizer.

Boundary constraints are inactive by default (after initial creation).
They are preserved after algorithm restart with MinBLEICRestartFrom().

INPUT PARAMETERS:
    State   -   structure stores algorithm state
    BndL    -   lower bounds, array[N].
                If some (all) variables are unbounded, you may specify
                very small number or -INF.
    BndU    -   upper bounds, array[N].
                If some (all) variables are unbounded, you may specify
                very large number or +INF.

NOTE 1: it is possible to specify BndL[i]=BndU[i]. In this case I-th
variable will be "frozen" at X[i]=BndL[i]=BndU[i].

NOTE 2: this solver has following useful properties:
* bound constraints are always satisfied exactly
* function is evaluated only INSIDE area specified by  bound  constraints,
  even  when  numerical  differentiation is used (algorithm adjusts  nodes
  according to boundary constraints)

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetbc(minbleicstate* state,
     /* Real    */ ae_vector* bndl,
     /* Real    */ ae_vector* bndu,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t n;


    n = state->nmain;
    ae_assert(bndl->cnt>=n, "MinBLEICSetBC: Length(BndL)<N", _state);
    ae_assert(bndu->cnt>=n, "MinBLEICSetBC: Length(BndU)<N", _state);
    for(i=0; i<=n-1; i++)
    {
        ae_assert(ae_isfinite(bndl->ptr.p_double[i], _state)||ae_isneginf(bndl->ptr.p_double[i], _state), "MinBLEICSetBC: BndL contains NAN or +INF", _state);
        ae_assert(ae_isfinite(bndu->ptr.p_double[i], _state)||ae_isposinf(bndu->ptr.p_double[i], _state), "MinBLEICSetBC: BndL contains NAN or -INF", _state);
        state->bndloriginal.ptr.p_double[i] = bndl->ptr.p_double[i];
        state->hasbndl.ptr.p_bool[i] = ae_isfinite(bndl->ptr.p_double[i], _state);
        state->bnduoriginal.ptr.p_double[i] = bndu->ptr.p_double[i];
        state->hasbndu.ptr.p_bool[i] = ae_isfinite(bndu->ptr.p_double[i], _state);
    }
}


/*************************************************************************
This function sets linear constraints for BLEIC optimizer.

Linear constraints are inactive by default (after initial creation).
They are preserved after algorithm restart with MinBLEICRestartFrom().

INPUT PARAMETERS:
    State   -   structure previously allocated with MinBLEICCreate call.
    C       -   linear constraints, array[K,N+1].
                Each row of C represents one constraint, either equality
                or inequality (see below):
                * first N elements correspond to coefficients,
                * last element corresponds to the right part.
                All elements of C (including right part) must be finite.
    CT      -   type of constraints, array[K]:
                * if CT[i]>0, then I-th constraint is C[i,*]*x >= C[i,n+1]
                * if CT[i]=0, then I-th constraint is C[i,*]*x  = C[i,n+1]
                * if CT[i]<0, then I-th constraint is C[i,*]*x <= C[i,n+1]
    K       -   number of equality/inequality constraints, K>=0:
                * if given, only leading K elements of C/CT are used
                * if not given, automatically determined from sizes of C/CT

NOTE 1: linear (non-bound) constraints are satisfied only approximately:
* there always exists some minor violation (about Epsilon in magnitude)
  due to rounding errors
* numerical differentiation, if used, may  lead  to  function  evaluations
  outside  of the feasible  area,   because   algorithm  does  NOT  change
  numerical differentiation formula according to linear constraints.
If you want constraints to be  satisfied  exactly, try to reformulate your
problem  in  such  manner  that  all constraints will become boundary ones
(this kind of constraints is always satisfied exactly, both in  the  final
solution and in all intermediate points).

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetlc(minbleicstate* state,
     /* Real    */ ae_matrix* c,
     /* Integer */ ae_vector* ct,
     ae_int_t k,
     ae_state *_state)
{
    ae_int_t nmain;
    ae_int_t i;


    nmain = state->nmain;
    
    /*
     * First, check for errors in the inputs
     */
    ae_assert(k>=0, "MinBLEICSetLC: K<0", _state);
    ae_assert(c->cols>=nmain+1||k==0, "MinBLEICSetLC: Cols(C)<N+1", _state);
    ae_assert(c->rows>=k, "MinBLEICSetLC: Rows(C)<K", _state);
    ae_assert(ct->cnt>=k, "MinBLEICSetLC: Length(CT)<K", _state);
    ae_assert(apservisfinitematrix(c, k, nmain+1, _state), "MinBLEICSetLC: C contains infinite or NaN values!", _state);
    
    /*
     * Determine number of constraints,
     * allocate space and copy
     */
    state->cecnt = k;
    rmatrixsetlengthatleast(&state->ceoriginal, state->cecnt, nmain+1, _state);
    ivectorsetlengthatleast(&state->ct, state->cecnt, _state);
    for(i=0; i<=k-1; i++)
    {
        state->ct.ptr.p_int[i] = ct->ptr.p_int[i];
        ae_v_move(&state->ceoriginal.ptr.pp_double[i][0], 1, &c->ptr.pp_double[i][0], 1, ae_v_len(0,nmain));
    }
}


/*************************************************************************
This function sets stopping conditions for the underlying nonlinear CG
optimizer. It controls overall accuracy of solution. These conditions
should be strict enough in order for algorithm to converge.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    EpsG    -   >=0
                The  subroutine  finishes  its  work   if   the  condition
                |v|<EpsG is satisfied, where:
                * |.| means Euclidian norm
                * v - scaled gradient vector, v[i]=g[i]*s[i]
                * g - gradient
                * s - scaling coefficients set by MinBLEICSetScale()
    EpsF    -   >=0
                The  subroutine  finishes  its work if on k+1-th iteration
                the  condition  |F(k+1)-F(k)|<=EpsF*max{|F(k)|,|F(k+1)|,1}
                is satisfied.
    EpsX    -   >=0
                The subroutine finishes its work if  on  k+1-th  iteration
                the condition |v|<=EpsX is fulfilled, where:
                * |.| means Euclidian norm
                * v - scaled step vector, v[i]=dx[i]/s[i]
                * dx - ste pvector, dx=X(k+1)-X(k)
                * s - scaling coefficients set by MinBLEICSetScale()

Passing EpsG=0, EpsF=0 and EpsX=0 (simultaneously) will lead to
automatic stopping criterion selection.

These conditions are used to terminate inner iterations. However, you
need to tune termination conditions for outer iterations too.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetinnercond(minbleicstate* state,
     double epsg,
     double epsf,
     double epsx,
     ae_state *_state)
{


    ae_assert(ae_isfinite(epsg, _state), "MinBLEICSetInnerCond: EpsG is not finite number", _state);
    ae_assert(ae_fp_greater_eq(epsg,0), "MinBLEICSetInnerCond: negative EpsG", _state);
    ae_assert(ae_isfinite(epsf, _state), "MinBLEICSetInnerCond: EpsF is not finite number", _state);
    ae_assert(ae_fp_greater_eq(epsf,0), "MinBLEICSetInnerCond: negative EpsF", _state);
    ae_assert(ae_isfinite(epsx, _state), "MinBLEICSetInnerCond: EpsX is not finite number", _state);
    ae_assert(ae_fp_greater_eq(epsx,0), "MinBLEICSetInnerCond: negative EpsX", _state);
    state->innerepsg = epsg;
    state->innerepsf = epsf;
    state->innerepsx = epsx;
}


/*************************************************************************
This function sets stopping conditions for outer iteration of BLEIC algo.

These conditions control accuracy of constraint handling and amount of
infeasibility allowed in the solution.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    EpsX    -   >0, stopping condition on outer iteration step length
    EpsI    -   >0, stopping condition on infeasibility
    
Both EpsX and EpsI must be non-zero.

MEANING OF EpsX

EpsX  is  a  stopping  condition for outer iterations. Algorithm will stop
when  solution  of  the  current  modified  subproblem will be within EpsX
(using 2-norm) of the previous solution.

MEANING OF EpsI

EpsI controls feasibility properties -  algorithm  won't  stop  until  all
inequality constraints will be satisfied with error (distance from current
point to the feasible area) at most EpsI.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetoutercond(minbleicstate* state,
     double epsx,
     double epsi,
     ae_state *_state)
{


    ae_assert(ae_isfinite(epsx, _state), "MinBLEICSetOuterCond: EpsX is not finite number", _state);
    ae_assert(ae_fp_greater(epsx,0), "MinBLEICSetOuterCond: non-positive EpsX", _state);
    ae_assert(ae_isfinite(epsi, _state), "MinBLEICSetOuterCond: EpsI is not finite number", _state);
    ae_assert(ae_fp_greater(epsi,0), "MinBLEICSetOuterCond: non-positive EpsI", _state);
    state->outerepsx = epsx;
    state->outerepsi = epsi;
}


/*************************************************************************
This function sets scaling coefficients for BLEIC optimizer.

ALGLIB optimizers use scaling matrices to test stopping  conditions  (step
size and gradient are scaled before comparison with tolerances).  Scale of
the I-th variable is a translation invariant measure of:
a) "how large" the variable is
b) how large the step should be to make significant changes in the function

Scaling is also used by finite difference variant of the optimizer  - step
along I-th axis is equal to DiffStep*S[I].

In  most  optimizers  (and  in  the  BLEIC  too)  scaling is NOT a form of
preconditioning. It just  affects  stopping  conditions.  You  should  set
preconditioner  by  separate  call  to  one  of  the  MinBLEICSetPrec...()
functions.

There is a special  preconditioning  mode, however,  which  uses   scaling
coefficients to form diagonal preconditioning matrix. You  can  turn  this
mode on, if you want.   But  you should understand that scaling is not the
same thing as preconditioning - these are two different, although  related
forms of tuning solver.

INPUT PARAMETERS:
    State   -   structure stores algorithm state
    S       -   array[N], non-zero scaling coefficients
                S[i] may be negative, sign doesn't matter.

  -- ALGLIB --
     Copyright 14.01.2011 by Bochkanov Sergey
*************************************************************************/
void minbleicsetscale(minbleicstate* state,
     /* Real    */ ae_vector* s,
     ae_state *_state)
{
    ae_int_t i;


    ae_assert(s->cnt>=state->nmain, "MinBLEICSetScale: Length(S)<N", _state);
    for(i=0; i<=state->nmain-1; i++)
    {
        ae_assert(ae_isfinite(s->ptr.p_double[i], _state), "MinBLEICSetScale: S contains infinite or NAN elements", _state);
        ae_assert(ae_fp_neq(s->ptr.p_double[i],0), "MinBLEICSetScale: S contains zero elements", _state);
        state->soriginal.ptr.p_double[i] = ae_fabs(s->ptr.p_double[i], _state);
    }
}


/*************************************************************************
Modification of the preconditioner: preconditioning is turned off.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetprecdefault(minbleicstate* state, ae_state *) //_state)
{


    state->prectype = 0;
}


/*************************************************************************
Modification  of  the  preconditioner:  diagonal of approximate Hessian is
used.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    D       -   diagonal of the approximate Hessian, array[0..N-1],
                (if larger, only leading N elements are used).

NOTE 1: D[i] should be positive. Exception will be thrown otherwise.

NOTE 2: you should pass diagonal of approximate Hessian - NOT ITS INVERSE.

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetprecdiag(minbleicstate* state,
     /* Real    */ ae_vector* d,
     ae_state *_state)
{
    ae_int_t i;


    ae_assert(d->cnt>=state->nmain, "MinBLEICSetPrecDiag: D is too short", _state);
    for(i=0; i<=state->nmain-1; i++)
    {
        ae_assert(ae_isfinite(d->ptr.p_double[i], _state), "MinBLEICSetPrecDiag: D contains infinite or NAN elements", _state);
        ae_assert(ae_fp_greater(d->ptr.p_double[i],0), "MinBLEICSetPrecDiag: D contains non-positive elements", _state);
    }
    rvectorsetlengthatleast(&state->diaghoriginal, state->nmain, _state);
    state->prectype = 2;
    for(i=0; i<=state->nmain-1; i++)
    {
        state->diaghoriginal.ptr.p_double[i] = d->ptr.p_double[i];
    }
}


/*************************************************************************
Modification of the preconditioner: scale-based diagonal preconditioning.

This preconditioning mode can be useful when you  don't  have  approximate
diagonal of Hessian, but you know that your  variables  are  badly  scaled
(for  example,  one  variable is in [1,10], and another in [1000,100000]),
and most part of the ill-conditioning comes from different scales of vars.

In this case simple  scale-based  preconditioner,  with H[i] = 1/(s[i]^2),
can greatly improve convergence.

IMPRTANT: you should set scale of your variables  with  MinBLEICSetScale()
call  (before  or after MinBLEICSetPrecScale() call). Without knowledge of
the scale of your variables scale-based preconditioner will be  just  unit
matrix.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetprecscale(minbleicstate* state, ae_state *) //_state)
{


    state->prectype = 3;
}


/*************************************************************************
This function allows to stop algorithm after specified number of inner
iterations.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    MaxIts  -   maximum number of inner iterations.
                If MaxIts=0, the number of iterations is unlimited.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetmaxits(minbleicstate* state,
     ae_int_t maxits,
     ae_state *_state)
{


    ae_assert(maxits>=0, "MinBLEICSetMaxIts: negative MaxIts!", _state);
    state->maxits = maxits;
}


/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    NeedXRep-   whether iteration reports are needed or not

If NeedXRep is True, algorithm will call rep() callback function if  it is
provided to MinBLEICOptimize().

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetxrep(minbleicstate* state,
     ae_bool needxrep,
     ae_state *) //_state)
{


    state->xrep = needxrep;
}


/*************************************************************************
This function sets maximum step length

IMPORTANT: this feature is hard to combine with preconditioning. You can't
set upper limit on step length, when you solve optimization  problem  with
linear (non-boundary) constraints AND preconditioner turned on.

When  non-boundary  constraints  are  present,  you  have to either a) use
preconditioner, or b) use upper limit on step length.  YOU CAN'T USE BOTH!
In this case algorithm will terminate with appropriate error code.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    StpMax  -   maximum step length, >=0. Set StpMax to 0.0,  if you don't
                want to limit step length.

Use this subroutine when you optimize target function which contains exp()
or  other  fast  growing  functions,  and optimization algorithm makes too
large  steps  which  lead   to overflow. This function allows us to reject
steps  that  are  too  large  (and  therefore  expose  us  to the possible
overflow) without actually calculating function value at the x+stp*d.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetstpmax(minbleicstate* state,
     double stpmax,
     ae_state *_state)
{


    ae_assert(ae_isfinite(stpmax, _state), "MinBLEICSetStpMax: StpMax is not finite!", _state);
    ae_assert(ae_fp_greater_eq(stpmax,0), "MinBLEICSetStpMax: StpMax<0!", _state);
    state->stpmax = stpmax;
}


/*************************************************************************
NOTES:

1. This function has two different implementations: one which  uses  exact
   (analytical) user-supplied gradient,  and one which uses function value
   only  and  numerically  differentiates  function  in  order  to  obtain
   gradient.

   Depending  on  the  specific  function  used to create optimizer object
   (either  MinBLEICCreate() for analytical gradient or  MinBLEICCreateF()
   for numerical differentiation) you should choose appropriate variant of
   MinBLEICOptimize() - one  which  accepts  function  AND gradient or one
   which accepts function ONLY.

   Be careful to choose variant of MinBLEICOptimize() which corresponds to
   your optimization scheme! Table below lists different  combinations  of
   callback (function/gradient) passed to MinBLEICOptimize()  and specific
   function used to create optimizer.


                     |         USER PASSED TO MinBLEICOptimize()
   CREATED WITH      |  function only   |  function and gradient
   ------------------------------------------------------------
   MinBLEICCreateF() |     work                FAIL
   MinBLEICCreate()  |     FAIL                work

   Here "FAIL" denotes inappropriate combinations  of  optimizer  creation
   function  and  MinBLEICOptimize()  version.   Attemps   to   use   such
   combination (for  example,  to  create optimizer with MinBLEICCreateF()
   and  to  pass  gradient  information  to  MinCGOptimize()) will lead to
   exception being thrown. Either  you  did  not pass gradient when it WAS
   needed or you passed gradient when it was NOT needed.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
ae_bool minbleiciteration(minbleicstate* state, ae_state *_state)
{
    ae_int_t nmain;
    ae_int_t nslack;
    ae_int_t m;
    ae_int_t i;
    ae_int_t j;
    double v;
    double vv;
    ae_bool b;
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
        nmain = state->rstate.ia.ptr.p_int[0];
        nslack = state->rstate.ia.ptr.p_int[1];
        m = state->rstate.ia.ptr.p_int[2];
        i = state->rstate.ia.ptr.p_int[3];
        j = state->rstate.ia.ptr.p_int[4];
        b = state->rstate.ba.ptr.p_bool[0];
        v = state->rstate.ra.ptr.p_double[0];
        vv = state->rstate.ra.ptr.p_double[1];
    }
    else
    {
        nmain = -983;
        nslack = -989;
        m = -834;
        i = 900;
        j = -287;
        b = ae_false;
        v = 214;
        vv = -338;
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
    if( state->rstate.stage==4 )
    {
        goto lbl_4;
    }
    if( state->rstate.stage==5 )
    {
        goto lbl_5;
    }
    if( state->rstate.stage==6 )
    {
        goto lbl_6;
    }
    if( state->rstate.stage==7 )
    {
        goto lbl_7;
    }
    if( state->rstate.stage==8 )
    {
        goto lbl_8;
    }
    if( state->rstate.stage==9 )
    {
        goto lbl_9;
    }
    if( state->rstate.stage==10 )
    {
        goto lbl_10;
    }
    if( state->rstate.stage==11 )
    {
        goto lbl_11;
    }
    if( state->rstate.stage==12 )
    {
        goto lbl_12;
    }
    if( state->rstate.stage==13 )
    {
        goto lbl_13;
    }
    if( state->rstate.stage==14 )
    {
        goto lbl_14;
    }
    if( state->rstate.stage==15 )
    {
        goto lbl_15;
    }
    
    /*
     * Routine body
     */
    
    /*
     * Prepare:
     * * calculate number of slack variables
     * * initialize locals
     * * initialize debug fields
     * * make quick check
     */
    nmain = state->nmain;
    nslack = 0;
    for(i=0; i<=state->cecnt-1; i++)
    {
        if( state->ct.ptr.p_int[i]!=0 )
        {
            nslack = nslack+1;
        }
    }
    state->nslack = nslack;
    state->repterminationtype = 0;
    state->repinneriterationscount = 0;
    state->repouteriterationscount = 0;
    state->repnfev = 0;
    state->repvaridx = -1;
    state->repdebugeqerr = 0.0;
    state->repdebugfs = _state->v_nan;
    state->repdebugff = _state->v_nan;
    state->repdebugdx = _state->v_nan;
    if( ae_fp_neq(state->stpmax,0)&&state->prectype!=0 )
    {
        state->repterminationtype = -10;
        result = ae_false;
        return result;
    }
    
    /*
     * allocate
     */
    rvectorsetlengthatleast(&state->r, nmain+nslack, _state);
    rvectorsetlengthatleast(&state->diagh, nmain+nslack, _state);
    rvectorsetlengthatleast(&state->tmp0, nmain+nslack, _state);
    rvectorsetlengthatleast(&state->tmp1, nmain+nslack, _state);
    rvectorsetlengthatleast(&state->tmp2, nmain+nslack, _state);
    rmatrixsetlengthatleast(&state->cecurrent, state->cecnt, nmain+nslack+1, _state);
    bvectorsetlengthatleast(&state->activeconstraints, nmain+nslack, _state);
    rvectorsetlengthatleast(&state->constrainedvalues, nmain+nslack, _state);
    rvectorsetlengthatleast(&state->lastg, nmain+nslack, _state);
    rvectorsetlengthatleast(&state->xe, nmain+nslack, _state);
    rvectorsetlengthatleast(&state->xcur, nmain+nslack, _state);
    rvectorsetlengthatleast(&state->xprev, nmain+nslack, _state);
    rvectorsetlengthatleast(&state->xend, nmain, _state);
    
    /*
     * Create/restart optimizer.
     *
     * State.OptDim is used to determine current state of optimizer.
     */
    if( state->optdim!=nmain+nslack )
    {
        for(i=0; i<=nmain+nslack-1; i++)
        {
            state->tmp1.ptr.p_double[i] = 0.0;
        }
        mincgcreate(nmain+nslack, &state->tmp1, &state->cgstate, _state);
        state->optdim = nmain+nslack;
    }
    
    /*
     * Prepare transformation.
     *
     * MinBLEIC's handling of preconditioner matrix is somewhat unusual -
     * instead of incorporating it into algorithm and making implicit
     * scaling (as most optimizers do) BLEIC optimizer uses explicit
     * scaling - it solves problem in the scaled parameters space S,
     * making transition between scaled (S) and unscaled (X) variables
     * every time we ask for function value.
     *
     * Following fields are calculated here:
     * * TransformS         X[i] = TransformS[i]*S[i], array[NMain]
     * * SEffective         "effective" scale of the variables after
     *                      transformation, array[NMain+NSlack]
     */
    rvectorsetlengthatleast(&state->transforms, nmain, _state);
    for(i=0; i<=nmain-1; i++)
    {
        if( state->prectype==2 )
        {
            state->transforms.ptr.p_double[i] = 1/ae_sqrt(state->diaghoriginal.ptr.p_double[i], _state);
            continue;
        }
        if( state->prectype==3 )
        {
            state->transforms.ptr.p_double[i] = state->soriginal.ptr.p_double[i];
            continue;
        }
        state->transforms.ptr.p_double[i] = 1;
    }
    rvectorsetlengthatleast(&state->seffective, nmain+nslack, _state);
    for(i=0; i<=nmain-1; i++)
    {
        state->seffective.ptr.p_double[i] = state->soriginal.ptr.p_double[i]/state->transforms.ptr.p_double[i];
    }
    for(i=0; i<=nslack-1; i++)
    {
        state->seffective.ptr.p_double[nmain+i] = 1;
    }
    mincgsetscale(&state->cgstate, &state->seffective, _state);
    
    /*
     * Pre-process constraints
     * * check consistency of bound constraints
     * * add slack vars, convert problem to the bound/equality
     *   constrained one
     *
     * We calculate here:
     * * BndLEffective - lower bounds after transformation of variables (see above)
     * * BndUEffective - upper bounds after transformation of variables (see above)
     * * CEEffective - matrix of equality constraints for transformed variables
     */
    for(i=0; i<=nmain-1; i++)
    {
        if( state->hasbndl.ptr.p_bool[i] )
        {
            state->bndleffective.ptr.p_double[i] = state->bndloriginal.ptr.p_double[i]/state->transforms.ptr.p_double[i];
        }
        if( state->hasbndu.ptr.p_bool[i] )
        {
            state->bndueffective.ptr.p_double[i] = state->bnduoriginal.ptr.p_double[i]/state->transforms.ptr.p_double[i];
        }
    }
    for(i=0; i<=nmain-1; i++)
    {
        if( state->hasbndl.ptr.p_bool[i]&&state->hasbndu.ptr.p_bool[i] )
        {
            if( ae_fp_greater(state->bndleffective.ptr.p_double[i],state->bndueffective.ptr.p_double[i]) )
            {
                state->repterminationtype = -3;
                result = ae_false;
                return result;
            }
        }
    }
    rmatrixsetlengthatleast(&state->ceeffective, state->cecnt, nmain+nslack+1, _state);
    m = 0;
    for(i=0; i<=state->cecnt-1; i++)
    {
        
        /*
         * NOTE: when we add slack variable, we use V = max(abs(CE[i,...])) as
         * coefficient before it in order to make linear equations better
         * conditioned.
         */
        v = 0;
        for(j=0; j<=nmain-1; j++)
        {
            state->ceeffective.ptr.pp_double[i][j] = state->ceoriginal.ptr.pp_double[i][j]*state->transforms.ptr.p_double[j];
            v = ae_maxreal(v, ae_fabs(state->ceeffective.ptr.pp_double[i][j], _state), _state);
        }
        if( ae_fp_eq(v,0) )
        {
            v = 1;
        }
        for(j=0; j<=nslack-1; j++)
        {
            state->ceeffective.ptr.pp_double[i][nmain+j] = 0.0;
        }
        state->ceeffective.ptr.pp_double[i][nmain+nslack] = state->ceoriginal.ptr.pp_double[i][nmain];
        if( state->ct.ptr.p_int[i]<0 )
        {
            state->ceeffective.ptr.pp_double[i][nmain+m] = v;
            m = m+1;
        }
        if( state->ct.ptr.p_int[i]>0 )
        {
            state->ceeffective.ptr.pp_double[i][nmain+m] = -v;
            m = m+1;
        }
    }
    
    /*
     * Find feasible point.
     *
     * 0. Convert from unscaled values (as stored in XStart) to scaled
     *    ones
     * 1. calculate values of slack variables such that starting
     *    point satisfies inequality constraints (after conversion to
     *    equality ones) as much as possible.
     * 2. use PrepareConstraintMatrix() function, which forces X
     *    to be strictly feasible.
     */
    for(i=0; i<=nmain-1; i++)
    {
        state->tmp0.ptr.p_double[i] = state->xstart.ptr.p_double[i]/state->transforms.ptr.p_double[i];
    }
    m = 0;
    for(i=0; i<=state->cecnt-1; i++)
    {
        v = ae_v_dotproduct(&state->ceeffective.ptr.pp_double[i][0], 1, &state->tmp0.ptr.p_double[0], 1, ae_v_len(0,nmain-1));
        if( state->ct.ptr.p_int[i]<0 )
        {
            state->tmp0.ptr.p_double[nmain+m] = state->ceeffective.ptr.pp_double[i][nmain+nslack]-v;
            m = m+1;
        }
        if( state->ct.ptr.p_int[i]>0 )
        {
            state->tmp0.ptr.p_double[nmain+m] = v-state->ceeffective.ptr.pp_double[i][nmain+nslack];
            m = m+1;
        }
    }
    for(i=0; i<=nmain+nslack-1; i++)
    {
        state->tmp1.ptr.p_double[i] = 0;
    }
    for(i=0; i<=nmain+nslack-1; i++)
    {
        state->activeconstraints.ptr.p_bool[i] = ae_false;
    }
    b = findfeasiblepoint(&state->tmp0, &state->bndleffective, &state->hasbndl, &state->bndueffective, &state->hasbndu, nmain, nslack, &state->ceeffective, state->cecnt, state->outerepsi, &state->repdebugfeasqpits, &state->repdebugfeasgpaits, _state);
    ae_v_move(&state->xcur.ptr.p_double[0], 1, &state->tmp0.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
    state->repdebugeqerr = 0.0;
    for(i=0; i<=state->cecnt-1; i++)
    {
        v = ae_v_dotproduct(&state->ceeffective.ptr.pp_double[i][0], 1, &state->xcur.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
        state->repdebugeqerr = state->repdebugeqerr+ae_sqr(v-state->ceeffective.ptr.pp_double[i][nmain+nslack], _state);
    }
    state->repdebugeqerr = ae_sqrt(state->repdebugeqerr, _state);
    if( !b )
    {
        state->repterminationtype = -3;
        result = ae_false;
        return result;
    }
    
    /*
     *  Check, that transferred derivative value is right
     */
    minbleic_unscalepoint(state, &state->xcur, &state->x, _state);
    minbleic_clearrequestfields(state, _state);
    if( !(ae_fp_eq(state->diffstep,0)&&ae_fp_greater(state->teststep,0)) )
    {
        goto lbl_16;
    }
    state->needfg = ae_true;
    i = 0;
lbl_18:
    if( i>nmain-1 )
    {
        goto lbl_20;
    }
    ae_assert((state->hasbndl.ptr.p_bool[i]&&ae_fp_less_eq(state->bndloriginal.ptr.p_double[i],state->x.ptr.p_double[i]))||!state->hasbndl.ptr.p_bool[i], "MinBLEICIteration: internal error(State.X is out of bounds)", _state);
    ae_assert((state->hasbndu.ptr.p_bool[i]&&ae_fp_less_eq(state->x.ptr.p_double[i],state->bnduoriginal.ptr.p_double[i]))||!state->hasbndu.ptr.p_bool[i], "MinBLEICIteration: internal error(State.X is out of bounds)", _state);
    v = state->x.ptr.p_double[i];
    state->x.ptr.p_double[i] = v-state->teststep*state->soriginal.ptr.p_double[i];
    if( state->hasbndl.ptr.p_bool[i] )
    {
        state->x.ptr.p_double[i] = ae_maxreal(state->x.ptr.p_double[i], state->bndloriginal.ptr.p_double[i], _state);
    }
    state->xm1 = state->x.ptr.p_double[i];
    state->rstate.stage = 0;
    goto lbl_rcomm;
lbl_0:
    state->fm1 = state->f;
    state->fp1 = state->g.ptr.p_double[i];
    state->x.ptr.p_double[i] = v+state->teststep*state->soriginal.ptr.p_double[i];
    if( state->hasbndu.ptr.p_bool[i] )
    {
        state->x.ptr.p_double[i] = ae_minreal(state->x.ptr.p_double[i], state->bnduoriginal.ptr.p_double[i], _state);
    }
    state->xp1 = state->x.ptr.p_double[i];
    state->rstate.stage = 1;
    goto lbl_rcomm;
lbl_1:
    state->fm2 = state->f;
    state->fp2 = state->g.ptr.p_double[i];
    state->x.ptr.p_double[i] = (state->xm1+state->xp1)/2;
    if( state->hasbndl.ptr.p_bool[i] )
    {
        state->x.ptr.p_double[i] = ae_maxreal(state->x.ptr.p_double[i], state->bndloriginal.ptr.p_double[i], _state);
    }
    if( state->hasbndu.ptr.p_bool[i] )
    {
        state->x.ptr.p_double[i] = ae_minreal(state->x.ptr.p_double[i], state->bnduoriginal.ptr.p_double[i], _state);
    }
    state->rstate.stage = 2;
    goto lbl_rcomm;
lbl_2:
    state->x.ptr.p_double[i] = v;
    if( !derivativecheck(state->fm1, state->fp1, state->fm2, state->fp2, state->f, state->g.ptr.p_double[i], state->xp1-state->xm1, _state) )
    {
        state->repvaridx = i;
        state->repterminationtype = -7;
        result = ae_false;
        return result;
    }
    i = i+1;
    goto lbl_18;
lbl_20:
    state->needfg = ae_false;
lbl_16:
    
    /*
     * Initialize RepDebugFS with function value at initial point
     */
    minbleic_unscalepoint(state, &state->xcur, &state->x, _state);
    minbleic_clearrequestfields(state, _state);
    if( ae_fp_neq(state->diffstep,0) )
    {
        goto lbl_21;
    }
    state->needfg = ae_true;
    state->rstate.stage = 3;
    goto lbl_rcomm;
lbl_3:
    state->needfg = ae_false;
    goto lbl_22;
lbl_21:
    state->needf = ae_true;
    state->rstate.stage = 4;
    goto lbl_rcomm;
lbl_4:
    state->needf = ae_false;
lbl_22:
    trimprepare(state->f, &state->trimthreshold, _state);
    state->repnfev = state->repnfev+1;
    state->repdebugfs = state->f;
    
    /*
     * Outer cycle
     */
    state->itsleft = state->maxits;
    ae_v_move(&state->xprev.ptr.p_double[0], 1, &state->xcur.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
lbl_23:
    if( ae_false )
    {
        goto lbl_24;
    }
    ae_assert(state->prectype==0||ae_fp_eq(state->stpmax,0), "MinBLEIC: internal error (-10)", _state);
    
    /*
     * Inner cycle: CG with projections and penalty functions
     */
    ae_v_move(&state->tmp0.ptr.p_double[0], 1, &state->xcur.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
    for(i=0; i<=nmain+nslack-1; i++)
    {
        state->tmp1.ptr.p_double[i] = 0;
        state->activeconstraints.ptr.p_bool[i] = ae_false;
    }
    if( !minbleic_prepareconstraintmatrix(state, &state->tmp0, &state->tmp1, &state->xcur, &state->tmp2, _state) )
    {
        state->repterminationtype = -3;
        result = ae_false;
        return result;
    }
    for(i=0; i<=nmain+nslack-1; i++)
    {
        state->activeconstraints.ptr.p_bool[i] = ae_false;
    }
    minbleic_rebuildcexe(state, _state);
    mincgrestartfrom(&state->cgstate, &state->xcur, _state);
    mincgsetcond(&state->cgstate, state->innerepsg, state->innerepsf, state->innerepsx, state->itsleft, _state);
    mincgsetxrep(&state->cgstate, state->xrep, _state);
    mincgsetdrep(&state->cgstate, ae_true, _state);
    mincgsetstpmax(&state->cgstate, state->stpmax, _state);
lbl_25:
    if( !mincgiteration(&state->cgstate, _state) )
    {
        goto lbl_26;
    }
    
    /*
     * process different requests/reports of inner optimizer
     */
    if( state->cgstate.algpowerup )
    {
        for(i=0; i<=nmain+nslack-1; i++)
        {
            state->activeconstraints.ptr.p_bool[i] = ae_false;
        }
        do
        {
            minbleic_rebuildcexe(state, _state);
            ae_v_move(&state->tmp1.ptr.p_double[0], 1, &state->cgstate.g.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
            minbleic_makegradientprojection(state, &state->tmp1, _state);
            b = ae_false;
            for(i=0; i<=nmain-1; i++)
            {
                if( !state->activeconstraints.ptr.p_bool[i] )
                {
                    if( state->hasbndl.ptr.p_bool[i] )
                    {
                        if( ae_fp_eq(state->cgstate.x.ptr.p_double[i],state->bndleffective.ptr.p_double[i])&&ae_fp_greater_eq(state->tmp1.ptr.p_double[i],0) )
                        {
                            state->activeconstraints.ptr.p_bool[i] = ae_true;
                            state->constrainedvalues.ptr.p_double[i] = state->bndleffective.ptr.p_double[i];
                            b = ae_true;
                        }
                    }
                    if( state->hasbndu.ptr.p_bool[i] )
                    {
                        if( ae_fp_eq(state->cgstate.x.ptr.p_double[i],state->bndueffective.ptr.p_double[i])&&ae_fp_less_eq(state->tmp1.ptr.p_double[i],0) )
                        {
                            state->activeconstraints.ptr.p_bool[i] = ae_true;
                            state->constrainedvalues.ptr.p_double[i] = state->bndueffective.ptr.p_double[i];
                            b = ae_true;
                        }
                    }
                }
            }
            for(i=0; i<=nslack-1; i++)
            {
                if( !state->activeconstraints.ptr.p_bool[nmain+i] )
                {
                    if( ae_fp_eq(state->cgstate.x.ptr.p_double[nmain+i],0)&&ae_fp_greater_eq(state->tmp1.ptr.p_double[nmain+i],0) )
                    {
                        state->activeconstraints.ptr.p_bool[nmain+i] = ae_true;
                        state->constrainedvalues.ptr.p_double[nmain+i] = 0;
                        b = ae_true;
                    }
                }
            }
        }
        while(b);
        ae_v_move(&state->cgstate.g.ptr.p_double[0], 1, &state->tmp1.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
        goto lbl_25;
    }
    if( state->cgstate.lsstart )
    {
        
        /*
         * Beginning of the line search: set upper limit on step size
         * to prevent algo from leaving feasible area.
         */
        if( ae_fp_eq(state->cgstate.curstpmax,0) )
        {
            state->cgstate.curstpmax = 1.0E50;
        }
        state->variabletofreeze = -1;
        for(i=0; i<=nmain-1; i++)
        {
            if( state->hasbndl.ptr.p_bool[i]&&ae_fp_less(state->cgstate.d.ptr.p_double[i],0) )
            {
                v = state->cgstate.curstpmax;
                vv = state->cgstate.x.ptr.p_double[i]-state->bndleffective.ptr.p_double[i];
                if( ae_fp_less(vv,0) )
                {
                    vv = 0;
                }
                state->cgstate.curstpmax = safeminposrv(vv, -state->cgstate.d.ptr.p_double[i], state->cgstate.curstpmax, _state);
                if( ae_fp_less(state->cgstate.curstpmax,v) )
                {
                    state->variabletofreeze = i;
                    state->valuetofreeze = state->bndleffective.ptr.p_double[i];
                }
            }
            if( state->hasbndu.ptr.p_bool[i]&&ae_fp_greater(state->cgstate.d.ptr.p_double[i],0) )
            {
                v = state->cgstate.curstpmax;
                vv = state->bndueffective.ptr.p_double[i]-state->cgstate.x.ptr.p_double[i];
                if( ae_fp_less(vv,0) )
                {
                    vv = 0;
                }
                state->cgstate.curstpmax = safeminposrv(vv, state->cgstate.d.ptr.p_double[i], state->cgstate.curstpmax, _state);
                if( ae_fp_less(state->cgstate.curstpmax,v) )
                {
                    state->variabletofreeze = i;
                    state->valuetofreeze = state->bndueffective.ptr.p_double[i];
                }
            }
        }
        for(i=0; i<=nslack-1; i++)
        {
            if( ae_fp_less(state->cgstate.d.ptr.p_double[nmain+i],0) )
            {
                v = state->cgstate.curstpmax;
                vv = state->cgstate.x.ptr.p_double[nmain+i];
                if( ae_fp_less(vv,0) )
                {
                    vv = 0;
                }
                state->cgstate.curstpmax = safeminposrv(vv, -state->cgstate.d.ptr.p_double[nmain+i], state->cgstate.curstpmax, _state);
                if( ae_fp_less(state->cgstate.curstpmax,v) )
                {
                    state->variabletofreeze = nmain+i;
                    state->valuetofreeze = 0;
                }
            }
        }
        if( ae_fp_eq(state->cgstate.curstpmax,0) )
        {
            state->activeconstraints.ptr.p_bool[state->variabletofreeze] = ae_true;
            state->constrainedvalues.ptr.p_double[state->variabletofreeze] = state->valuetofreeze;
            state->cgstate.x.ptr.p_double[state->variabletofreeze] = state->valuetofreeze;
            state->cgstate.terminationneeded = ae_true;
        }
        goto lbl_25;
    }
    if( state->cgstate.lsend )
    {
        
        /*
         * Line search just finished.
         * Maybe we should activate some constraints?
         */
        if( state->variabletofreeze>=0 )
        {
            
            /*
             * We have candidate constraint for activation, which can be
             * activated in two cases:
             *
             * 1. monotonic step with Stp=StpMax is performed (MCINFO=5),
             *    we just have to post-process this step to ensure that
             *    constraint is activated
             *
             * 2. MCINFO<>5 and MCINFO<>1, and BOTH of the following holds:
             *    2.a) ScaledNorm(State.CGState.G)>InnerEpsG (strict inequality to handle EpsG=0)
             *    2.b) ScaledNorm(StpMax*State.CGState.D)<=MaxNonMonotonicLen*MachineEpsilon;
             *
             *    What do these conditions mean?
             *    a) Having MCINFO which is not 5 or 1 means that it is
             *       possible that our function is actually increasing in
             *       the direction of D instead of decreasing.
             *    b) clause (2.a) means that gradient value is large enough
             *       to consider it non-zero (EpsG=0 means that any non-zero
             *       gradient value will be large enough).
             *    c) clause (2.b) means that step length is very small and
             *       function non-monotonicity can be just a consequence of
             *       numerical noise.
             *
             *    Given all of the listed above, it may be worth trying to
             *    perform step which increases function value.
             *
             * NOTE: some non-monotonic algorithms can be caught in an eternal
             *       loop. However, this algorithm is guaranteed from having
             *       such eternal loop because non-monotonic step always activates
             *       at least one constraint, and in order to have eternal loop
             *       we need to deactivate this constraint, which can be done
             *       only at the beginning of the outer iteration. And we have
             *       upper limit on the number of outer iterations, so we can
             *       be sure that algorithm will stop eventually.
             */
            if( state->cgstate.mcinfo!=5&&state->cgstate.mcinfo!=1 )
            {
                v = 0;
                vv = 0;
                for(i=0; i<=nmain+nslack-1; i++)
                {
                    v = v+ae_sqr(state->cgstate.g.ptr.p_double[i]*state->seffective.ptr.p_double[i], _state);
                    vv = vv+ae_sqr(state->cgstate.curstpmax*state->cgstate.d.ptr.p_double[i]/state->seffective.ptr.p_double[i], _state);
                }
                v = ae_sqrt(v, _state);
                vv = ae_sqrt(vv, _state);
                b = ae_fp_greater(v,state->innerepsg)&&ae_fp_less_eq(vv,minbleic_maxnonmonotoniclen*ae_machineepsilon);
            }
            else
            {
                b = state->cgstate.mcinfo==5;
            }
            if( b )
            {
                state->activeconstraints.ptr.p_bool[state->variabletofreeze] = ae_true;
                state->constrainedvalues.ptr.p_double[state->variabletofreeze] = state->valuetofreeze;
            }
        }
        else
        {
            b = ae_false;
        }
        
        /*
         * Additional activation of constraints
         */
        b = b||minbleic_additionalcheckforconstraints(state, &state->cgstate.x, _state);
        
        /*
         * If at least one constraint was activated we have to rebuild constraint matrices
         */
        if( b )
        {
            ae_v_move(&state->tmp0.ptr.p_double[0], 1, &state->cgstate.x.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
            ae_v_move(&state->tmp1.ptr.p_double[0], 1, &state->lastg.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
            if( !minbleic_prepareconstraintmatrix(state, &state->tmp0, &state->tmp1, &state->cgstate.x, &state->cgstate.g, _state) )
            {
                state->repterminationtype = -3;
                result = ae_false;
                return result;
            }
            state->cgstate.innerresetneeded = ae_true;
        }
        goto lbl_25;
    }
    if( !state->cgstate.needfg )
    {
        goto lbl_27;
    }
    ae_v_move(&state->tmp1.ptr.p_double[0], 1, &state->cgstate.x.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
    minbleic_projectpointandunscale(state, &state->tmp1, &state->x, &state->r, &vv, _state);
    minbleic_clearrequestfields(state, _state);
    if( ae_fp_neq(state->diffstep,0) )
    {
        goto lbl_29;
    }
    state->needfg = ae_true;
    state->rstate.stage = 5;
    goto lbl_rcomm;
lbl_5:
    state->needfg = ae_false;
    goto lbl_30;
lbl_29:
    state->needf = ae_true;
    state->rstate.stage = 6;
    goto lbl_rcomm;
lbl_6:
    state->fbase = state->f;
    i = 0;
lbl_31:
    if( i>nmain-1 )
    {
        goto lbl_33;
    }
    v = state->x.ptr.p_double[i];
    b = ae_false;
    if( state->hasbndl.ptr.p_bool[i] )
    {
        b = b||ae_fp_less(v-state->diffstep*state->soriginal.ptr.p_double[i],state->bndloriginal.ptr.p_double[i]);
    }
    if( state->hasbndu.ptr.p_bool[i] )
    {
        b = b||ae_fp_greater(v+state->diffstep*state->soriginal.ptr.p_double[i],state->bnduoriginal.ptr.p_double[i]);
    }
    if( b )
    {
        goto lbl_34;
    }
    state->x.ptr.p_double[i] = v-state->diffstep*state->soriginal.ptr.p_double[i];
    state->rstate.stage = 7;
    goto lbl_rcomm;
lbl_7:
    state->fm2 = state->f;
    state->x.ptr.p_double[i] = v-0.5*state->diffstep*state->soriginal.ptr.p_double[i];
    state->rstate.stage = 8;
    goto lbl_rcomm;
lbl_8:
    state->fm1 = state->f;
    state->x.ptr.p_double[i] = v+0.5*state->diffstep*state->soriginal.ptr.p_double[i];
    state->rstate.stage = 9;
    goto lbl_rcomm;
lbl_9:
    state->fp1 = state->f;
    state->x.ptr.p_double[i] = v+state->diffstep*state->soriginal.ptr.p_double[i];
    state->rstate.stage = 10;
    goto lbl_rcomm;
lbl_10:
    state->fp2 = state->f;
    state->g.ptr.p_double[i] = (8*(state->fp1-state->fm1)-(state->fp2-state->fm2))/(6*state->diffstep*state->soriginal.ptr.p_double[i]);
    goto lbl_35;
lbl_34:
    state->xm1 = ae_maxreal(v-state->diffstep*state->soriginal.ptr.p_double[i], state->bndloriginal.ptr.p_double[i], _state);
    state->x.ptr.p_double[i] = state->xm1;
    state->rstate.stage = 11;
    goto lbl_rcomm;
lbl_11:
    state->fm1 = state->f;
    state->xp1 = ae_minreal(v+state->diffstep*state->soriginal.ptr.p_double[i], state->bnduoriginal.ptr.p_double[i], _state);
    state->x.ptr.p_double[i] = state->xp1;
    state->rstate.stage = 12;
    goto lbl_rcomm;
lbl_12:
    state->fp1 = state->f;
    state->g.ptr.p_double[i] = (state->fp1-state->fm1)/(state->xp1-state->xm1);
lbl_35:
    state->x.ptr.p_double[i] = v;
    i = i+1;
    goto lbl_31;
lbl_33:
    state->f = state->fbase;
    state->needf = ae_false;
lbl_30:
    if( ae_fp_less(state->f,state->trimthreshold) )
    {
        
        /*
         * normal processing
         */
        state->cgstate.f = state->f;
        minbleic_scalegradientandexpand(state, &state->g, &state->cgstate.g, _state);
        ae_v_move(&state->lastg.ptr.p_double[0], 1, &state->cgstate.g.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
        minbleic_modifytargetfunction(state, &state->tmp1, &state->r, vv, &state->cgstate.f, &state->cgstate.g, &state->gnorm, &state->mpgnorm, _state);
    }
    else
    {
        
        /*
         * function value is too high, trim it
         */
        state->cgstate.f = state->trimthreshold;
        for(i=0; i<=nmain+nslack-1; i++)
        {
            state->cgstate.g.ptr.p_double[i] = 0.0;
        }
    }
    goto lbl_25;
lbl_27:
    if( !state->cgstate.xupdated )
    {
        goto lbl_36;
    }
    
    /*
     * Report
     */
    minbleic_unscalepoint(state, &state->cgstate.x, &state->x, _state);
    state->f = state->cgstate.f;
    minbleic_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 13;
    goto lbl_rcomm;
lbl_13:
    state->xupdated = ae_false;
    goto lbl_25;
lbl_36:
    goto lbl_25;
lbl_26:
    mincgresults(&state->cgstate, &state->xcur, &state->cgrep, _state);
    minbleic_unscalepoint(state, &state->xcur, &state->xend, _state);
    state->repinneriterationscount = state->repinneriterationscount+state->cgrep.iterationscount;
    state->repouteriterationscount = state->repouteriterationscount+1;
    state->repnfev = state->repnfev+state->cgrep.nfev;
    
    /*
     * Update RepDebugFF with function value at current point
     */
    minbleic_unscalepoint(state, &state->xcur, &state->x, _state);
    minbleic_clearrequestfields(state, _state);
    if( ae_fp_neq(state->diffstep,0) )
    {
        goto lbl_38;
    }
    state->needfg = ae_true;
    state->rstate.stage = 14;
    goto lbl_rcomm;
lbl_14:
    state->needfg = ae_false;
    goto lbl_39;
lbl_38:
    state->needf = ae_true;
    state->rstate.stage = 15;
    goto lbl_rcomm;
lbl_15:
    state->needf = ae_false;
lbl_39:
    state->repnfev = state->repnfev+1;
    state->repdebugff = state->f;
    
    /*
     * Check for stopping:
     * * "normal", outer step size is small enough, infeasibility is within bounds
     * * "inconsistent",  if Lagrange multipliers increased beyond threshold given by MaxLagrangeMul
     * * "too stringent", in other cases
     */
    v = 0;
    for(i=0; i<=nmain-1; i++)
    {
        v = v+ae_sqr((state->xcur.ptr.p_double[i]-state->xprev.ptr.p_double[i])/state->seffective.ptr.p_double[i], _state);
    }
    v = ae_sqrt(v, _state);
    if( ae_fp_less_eq(v,state->outerepsx) )
    {
        state->repterminationtype = 4;
        goto lbl_24;
    }
    if( state->maxits>0 )
    {
        state->itsleft = state->itsleft-state->cgrep.iterationscount;
        if( state->itsleft<=0 )
        {
            state->repterminationtype = 5;
            goto lbl_24;
        }
    }
    if( ae_fp_greater_eq(state->repouteriterationscount,minbleic_maxouterits) )
    {
        
        /*
         * Stopped because of reaching limit on the number of outer iterations.
         *
         * IMPORTANT: DO NOT REMOVE THIS STOPPING CRITERIA because it
         *            saves algorithm from eternal loop with non-monotonic
         *            steps being made (first step increases function value,
         *            next step moves us back). See comments on the code which
         *            handles State.CGState.LSEnd for more info on this subject.
         */
        state->repterminationtype = 5;
        goto lbl_24;
    }
    
    /*
     * Next iteration
     */
    ae_v_move(&state->xprev.ptr.p_double[0], 1, &state->xcur.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
    goto lbl_23;
lbl_24:
    
    /*
     * We've stopped, fill debug information
     */
    state->repdebugeqerr = 0.0;
    for(i=0; i<=state->cecnt-1; i++)
    {
        v = ae_v_dotproduct(&state->ceeffective.ptr.pp_double[i][0], 1, &state->xcur.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
        state->repdebugeqerr = state->repdebugeqerr+ae_sqr(v-state->ceeffective.ptr.pp_double[i][nmain+nslack], _state);
    }
    state->repdebugeqerr = ae_sqrt(state->repdebugeqerr, _state);
    state->repdebugdx = 0;
    for(i=0; i<=nmain-1; i++)
    {
        state->repdebugdx = state->repdebugdx+ae_sqr(state->xcur.ptr.p_double[i]-state->xstart.ptr.p_double[i], _state);
    }
    state->repdebugdx = ae_sqrt(state->repdebugdx, _state);
    result = ae_false;
    return result;
    
    /*
     * Saving state
     */
lbl_rcomm:
    result = ae_true;
    state->rstate.ia.ptr.p_int[0] = nmain;
    state->rstate.ia.ptr.p_int[1] = nslack;
    state->rstate.ia.ptr.p_int[2] = m;
    state->rstate.ia.ptr.p_int[3] = i;
    state->rstate.ia.ptr.p_int[4] = j;
    state->rstate.ba.ptr.p_bool[0] = b;
    state->rstate.ra.ptr.p_double[0] = v;
    state->rstate.ra.ptr.p_double[1] = vv;
    return result;
}


/*************************************************************************
BLEIC results

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    X       -   array[0..N-1], solution
    Rep     -   optimization report. You should check Rep.TerminationType
                in  order  to  distinguish  successful  termination  from
                unsuccessful one.
                More information about fields of this  structure  can  be
                found in the comments on MinBLEICReport datatype.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicresults(minbleicstate* state,
     /* Real    */ ae_vector* x,
     minbleicreport* rep,
     ae_state *_state)
{

    ae_vector_clear(x);
    _minbleicreport_clear(rep);

    minbleicresultsbuf(state, x, rep, _state);
}


/*************************************************************************
BLEIC results

Buffered implementation of MinBLEICResults() which uses pre-allocated buffer
to store X[]. If buffer size is  too  small,  it  resizes  buffer.  It  is
intended to be used in the inner cycles of performance critical algorithms
where array reallocation penalty is too large to be ignored.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicresultsbuf(minbleicstate* state,
     /* Real    */ ae_vector* x,
     minbleicreport* rep,
     ae_state *_state)
{
    ae_int_t i;


    if( x->cnt<state->nmain )
    {
        ae_vector_set_length(x, state->nmain, _state);
    }
    rep->inneriterationscount = state->repinneriterationscount;
    rep->outeriterationscount = state->repouteriterationscount;
    rep->nfev = state->repnfev;
    rep->varidx = state->repvaridx;
    rep->terminationtype = state->repterminationtype;
    if( state->repterminationtype>0 )
    {
        ae_v_move(&x->ptr.p_double[0], 1, &state->xend.ptr.p_double[0], 1, ae_v_len(0,state->nmain-1));
    }
    else
    {
        for(i=0; i<=state->nmain-1; i++)
        {
            x->ptr.p_double[i] = _state->v_nan;
        }
    }
    rep->debugeqerr = state->repdebugeqerr;
    rep->debugfs = state->repdebugfs;
    rep->debugff = state->repdebugff;
    rep->debugdx = state->repdebugdx;
    rep->debugfeasqpits = state->repdebugfeasqpits;
    rep->debugfeasgpaits = state->repdebugfeasgpaits;
}


/*************************************************************************
This subroutine restarts algorithm from new point.
All optimization parameters (including constraints) are left unchanged.

This  function  allows  to  solve multiple  optimization  problems  (which
must have  same number of dimensions) without object reallocation penalty.

INPUT PARAMETERS:
    State   -   structure previously allocated with MinBLEICCreate call.
    X       -   new starting point.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicrestartfrom(minbleicstate* state,
     /* Real    */ ae_vector* x,
     ae_state *_state)
{
    ae_int_t n;


    n = state->nmain;
    
    /*
     * First, check for errors in the inputs
     */
    ae_assert(x->cnt>=n, "MinBLEICRestartFrom: Length(X)<N", _state);
    ae_assert(isfinitevector(x, n, _state), "MinBLEICRestartFrom: X contains infinite or NaN values!", _state);
    
    /*
     * Set XC
     */
    ae_v_move(&state->xstart.ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,n-1));
    
    /*
     * prepare RComm facilities
     */
    ae_vector_set_length(&state->rstate.ia, 4+1, _state);
    ae_vector_set_length(&state->rstate.ba, 0+1, _state);
    ae_vector_set_length(&state->rstate.ra, 1+1, _state);
    state->rstate.stage = -1;
    minbleic_clearrequestfields(state, _state);
}


/*************************************************************************
This  subroutine  turns  on  verification  of  the  user-supplied analytic
gradient:
* user calls this subroutine before optimization begins
* MinBLEICOptimize() is called
* prior to  actual  optimization, for each component  of  parameters being
  optimized X[i] algorithm performs following steps:
  * two trial steps are made to X[i]-TestStep*S[i] and X[i]+TestStep*S[i],
    where X[i] is i-th component of the initial point and S[i] is a  scale
    of i-th parameter
  * if needed, steps are bounded with respect to constraints on X[]
  * F(X) is evaluated at these trial points
  * we perform one more evaluation in the middle point of the interval
  * we  build  cubic  model using function values and derivatives at trial
    points and we compare its prediction with actual value in  the  middle
    point
  * in case difference between prediction and actual value is higher  than
    some predetermined threshold, algorithm stops with completion code -7;
    Rep.VarIdx is set to index of the parameter with incorrect derivative.
* after verification is over, algorithm proceeds to the actual optimization.

NOTE 1: verification  needs  N (parameters count) gradient evaluations. It
        is very costly and you should use  it  only  for  low  dimensional
        problems,  when  you  want  to  be  sure  that  you've   correctly
        calculated  analytic  derivatives.  You  should  not use it in the
        production code (unless you want to check derivatives provided  by
        some third party).

NOTE 2: you  should  carefully  choose  TestStep. Value which is too large
        (so large that function behaviour is significantly non-cubic) will
        lead to false alarms. You may use  different  step  for  different
        parameters by means of setting scale with MinBLEICSetScale().

NOTE 3: this function may lead to false positives. In case it reports that
        I-th  derivative was calculated incorrectly, you may decrease test
        step  and  try  one  more  time  - maybe your function changes too
        sharply  and  your  step  is  too  large for such rapidly chanding
        function.

INPUT PARAMETERS:
    State       -   structure used to store algorithm state
    TestStep    -   verification step:
                    * TestStep=0 turns verification off
                    * TestStep>0 activates verification

  -- ALGLIB --
     Copyright 15.06.2012 by Bochkanov Sergey
*************************************************************************/
void minbleicsetgradientcheck(minbleicstate* state,
     double teststep,
     ae_state *_state)
{


    ae_assert(ae_isfinite(teststep, _state), "MinBLEICSetGradientCheck: TestStep contains NaN or Infinite", _state);
    ae_assert(ae_fp_greater_eq(teststep,0), "MinBLEICSetGradientCheck: invalid argument TestStep(TestStep<0)", _state);
    state->teststep = teststep;
}


/*************************************************************************
Clears request fileds (to be sure that we don't forget to clear something)
*************************************************************************/
static void minbleic_clearrequestfields(minbleicstate* state,
     ae_state *) //_state)
{


    state->needf = ae_false;
    state->needfg = ae_false;
    state->xupdated = ae_false;
}


/*************************************************************************
This functions "unscales" point, i.e. it makes transformation  from scaled
variables to unscaled ones. Only leading NMain variables are copied from
XUnscaled to XScaled.
*************************************************************************/
static void minbleic_unscalepoint(minbleicstate* state,
     /* Real    */ ae_vector* xscaled,
     /* Real    */ ae_vector* xunscaled,
     ae_state *) //_state)
{
    ae_int_t i;
    double v;


    for(i=0; i<=state->nmain-1; i++)
    {
        v = xscaled->ptr.p_double[i]*state->transforms.ptr.p_double[i];
        if( state->hasbndl.ptr.p_bool[i] )
        {
            if( ae_fp_less(v,state->bndloriginal.ptr.p_double[i]) )
            {
                v = state->bndloriginal.ptr.p_double[i];
            }
        }
        if( state->hasbndu.ptr.p_bool[i] )
        {
            if( ae_fp_greater(v,state->bnduoriginal.ptr.p_double[i]) )
            {
                v = state->bnduoriginal.ptr.p_double[i];
            }
        }
        xunscaled->ptr.p_double[i] = v;
    }
}


/*************************************************************************
This function:
1. makes projection of XScaled into equality constrained subspace
   (X is modified in-place)
2. stores residual from the projection into R
3. unscales projected XScaled and stores result into XUnscaled with
   additional enforcement
It calculates set of additional values which are used later for
modification of the target function F.

INPUT PARAMETERS:
    State   -   optimizer state (we use its fields to get information
                about constraints)
    X       -   vector being projected
    R       -   preallocated buffer, used to store residual from projection

OUTPUT PARAMETERS:
    X       -   projection of input X
    R       -   residual
    RNorm   -   residual norm squared, used later to modify target function
*************************************************************************/
static void minbleic_projectpointandunscale(minbleicstate* state,
     /* Real    */ ae_vector* xscaled,
     /* Real    */ ae_vector* xunscaled,
     /* Real    */ ae_vector* rscaled,
     double* rnorm2,
     ae_state *_state)
{
    double v;
    ae_int_t i;
    ae_int_t nmain;
    ae_int_t nslack;

    *rnorm2 = 0;

    nmain = state->nmain;
    nslack = state->nslack;
    
    /*
     * * subtract XE from XScaled
     * * project XScaled
     * * calculate norm of deviation from null space, store it in RNorm2
     * * calculate residual from projection, store it in R
     * * add XE to XScaled
     * * unscale variables
     */
    ae_v_sub(&xscaled->ptr.p_double[0], 1, &state->xe.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
    *rnorm2 = 0;
    for(i=0; i<=nmain+nslack-1; i++)
    {
        rscaled->ptr.p_double[i] = 0;
    }
    for(i=0; i<=nmain+nslack-1; i++)
    {
        if( state->activeconstraints.ptr.p_bool[i] )
        {
            v = xscaled->ptr.p_double[i];
            xscaled->ptr.p_double[i] = 0;
            rscaled->ptr.p_double[i] = rscaled->ptr.p_double[i]+v;
            *rnorm2 = *rnorm2+ae_sqr(v, _state);
        }
    }
    for(i=0; i<=state->cecnt-1; i++)
    {
        v = ae_v_dotproduct(&xscaled->ptr.p_double[0], 1, &state->cecurrent.ptr.pp_double[i][0], 1, ae_v_len(0,nmain+nslack-1));
        ae_v_subd(&xscaled->ptr.p_double[0], 1, &state->cecurrent.ptr.pp_double[i][0], 1, ae_v_len(0,nmain+nslack-1), v);
        ae_v_addd(&rscaled->ptr.p_double[0], 1, &state->cecurrent.ptr.pp_double[i][0], 1, ae_v_len(0,nmain+nslack-1), v);
        *rnorm2 = *rnorm2+ae_sqr(v, _state);
    }
    ae_v_add(&xscaled->ptr.p_double[0], 1, &state->xe.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
    minbleic_unscalepoint(state, xscaled, xunscaled, _state);
}


/*************************************************************************
This function scales and copies NMain elements of GUnscaled into GScaled.
Other NSlack components of GScaled are set to zero.
*************************************************************************/
static void minbleic_scalegradientandexpand(minbleicstate* state,
     /* Real    */ ae_vector* gunscaled,
     /* Real    */ ae_vector* gscaled,
     ae_state *) //_state)
{
    ae_int_t i;


    for(i=0; i<=state->nmain-1; i++)
    {
        gscaled->ptr.p_double[i] = gunscaled->ptr.p_double[i]*state->transforms.ptr.p_double[i];
    }
    for(i=0; i<=state->nslack-1; i++)
    {
        gscaled->ptr.p_double[state->nmain+i] = 0;
    }
}


/*************************************************************************
This subroutine applies modifications to the target function given by
its value F and gradient G at the projected point X which lies in the
equality constrained subspace.

Following modifications are applied:
* modified barrier functions to handle inequality constraints
  (both F and G are modified)
* projection of gradient into equality constrained subspace
  (only G is modified)
* quadratic penalty for deviations from equality constrained subspace
  (both F and G are modified)

It also calculates gradient norm (three different norms for three
different types of gradient), feasibility and complementary slackness
errors.

INPUT PARAMETERS:
    State   -   optimizer state (we use its fields to get information
                about constraints)
    X       -   point (projected into equality constrained subspace)
    R       -   residual from projection
    RNorm2  -   residual norm squared
    F       -   function value at X
    G       -   function gradient at X

OUTPUT PARAMETERS:
    F       -   modified function value at X
    G       -   modified function gradient at X
    GNorm   -   2-norm of unmodified G
    MPGNorm -   2-norm of modified G
    MBA     -   minimum argument of barrier functions.
                If X is strictly feasible, it is greater than zero.
                If X lies on a boundary, it is zero.
                It is negative for infeasible X.
    FIErr   -   2-norm of feasibility error with respect to
                inequality/bound constraints
    CSErr   -   2-norm of complementarity slackness error
*************************************************************************/
static void minbleic_modifytargetfunction(minbleicstate* state,
     /* Real    */ ae_vector* , //x,
     /* Real    */ ae_vector* r,
     double rnorm2,
     double* f,
     /* Real    */ ae_vector* g,
     double* gnorm,
     double* mpgnorm,
     ae_state *_state)
{
    double v;
    ae_int_t i;
    ae_int_t nmain;
    ae_int_t nslack;
    ae_bool hasconstraints;

    *gnorm = 0;
    *mpgnorm = 0;

    nmain = state->nmain;
    nslack = state->nslack;
    hasconstraints = ae_false;
    
    /*
     * GNorm
     */
    v = ae_v_dotproduct(&g->ptr.p_double[0], 1, &g->ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
    *gnorm = ae_sqrt(v, _state);
    
    /*
     * Process equality constraints:
     * * modify F to handle penalty term for equality constraints
     * * project gradient on null space of equality constraints
     * * add penalty term for equality constraints to gradient
     */
    *f = *f+rnorm2;
    for(i=0; i<=nmain+nslack-1; i++)
    {
        if( state->activeconstraints.ptr.p_bool[i] )
        {
            g->ptr.p_double[i] = 0;
        }
    }
    for(i=0; i<=state->cecnt-1; i++)
    {
        v = ae_v_dotproduct(&g->ptr.p_double[0], 1, &state->cecurrent.ptr.pp_double[i][0], 1, ae_v_len(0,nmain+nslack-1));
        ae_v_subd(&g->ptr.p_double[0], 1, &state->cecurrent.ptr.pp_double[i][0], 1, ae_v_len(0,nmain+nslack-1), v);
    }
    ae_v_addd(&g->ptr.p_double[0], 1, &r->ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1), 2);
    
    /*
     * MPGNorm
     */
    v = ae_v_dotproduct(&g->ptr.p_double[0], 1, &g->ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
    *mpgnorm = ae_sqrt(v, _state);
}


/*************************************************************************
This function makes additional check for constraints which can be activated.

We try activate constraints one by one, but it is possible that several
constraints should be activated during one iteration. It this case only
one of them (probably last) will be activated. This function will fix it -
it will pass through constraints and activate those which are at the boundary
or beyond it.

It will return True, if at least one constraint was activated by this function.
*************************************************************************/
static ae_bool minbleic_additionalcheckforconstraints(minbleicstate* state,
     /* Real    */ ae_vector* x,
     ae_state *) //_state)
{
    ae_int_t i;
    ae_int_t nmain;
    ae_int_t nslack;
    ae_bool result;


    result = ae_false;
    nmain = state->nmain;
    nslack = state->nslack;
    for(i=0; i<=nmain-1; i++)
    {
        if( !state->activeconstraints.ptr.p_bool[i] )
        {
            if( state->hasbndl.ptr.p_bool[i] )
            {
                if( ae_fp_less_eq(x->ptr.p_double[i],state->bndleffective.ptr.p_double[i]) )
                {
                    state->activeconstraints.ptr.p_bool[i] = ae_true;
                    state->constrainedvalues.ptr.p_double[i] = state->bndleffective.ptr.p_double[i];
                    result = ae_true;
                }
            }
            if( state->hasbndu.ptr.p_bool[i] )
            {
                if( ae_fp_greater_eq(x->ptr.p_double[i],state->bndueffective.ptr.p_double[i]) )
                {
                    state->activeconstraints.ptr.p_bool[i] = ae_true;
                    state->constrainedvalues.ptr.p_double[i] = state->bndueffective.ptr.p_double[i];
                    result = ae_true;
                }
            }
        }
    }
    for(i=0; i<=nslack-1; i++)
    {
        if( !state->activeconstraints.ptr.p_bool[nmain+i] )
        {
            if( ae_fp_less_eq(x->ptr.p_double[nmain+i],0) )
            {
                state->activeconstraints.ptr.p_bool[nmain+i] = ae_true;
                state->constrainedvalues.ptr.p_double[nmain+i] = 0;
                result = ae_true;
            }
        }
    }
    return result;
}


/*************************************************************************
This function rebuilds CECurrent and XE according to current set of
active bound constraints.
*************************************************************************/
static void minbleic_rebuildcexe(minbleicstate* state, ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t nmain;
    ae_int_t nslack;
    double v;


    nmain = state->nmain;
    nslack = state->nslack;
    rmatrixcopy(state->cecnt, nmain+nslack+1, &state->ceeffective, 0, 0, &state->cecurrent, 0, 0, _state);
    for(i=0; i<=state->cecnt-1; i++)
    {
        
        /*
         * "Subtract" active bound constraints from I-th linear constraint
         */
        for(j=0; j<=nmain+nslack-1; j++)
        {
            if( state->activeconstraints.ptr.p_bool[j] )
            {
                state->cecurrent.ptr.pp_double[i][nmain+nslack] = state->cecurrent.ptr.pp_double[i][nmain+nslack]-state->cecurrent.ptr.pp_double[i][j]*state->constrainedvalues.ptr.p_double[j];
                state->cecurrent.ptr.pp_double[i][j] = 0.0;
            }
        }
        
        /*
         * Reorthogonalize I-th constraint with respect to previous ones
         * NOTE: we also update right part, which is CECurrent[...,NMain+NSlack].
         */
        for(k=0; k<=i-1; k++)
        {
            v = ae_v_dotproduct(&state->cecurrent.ptr.pp_double[k][0], 1, &state->cecurrent.ptr.pp_double[i][0], 1, ae_v_len(0,nmain+nslack-1));
            ae_v_subd(&state->cecurrent.ptr.pp_double[i][0], 1, &state->cecurrent.ptr.pp_double[k][0], 1, ae_v_len(0,nmain+nslack), v);
        }
        
        /*
         * Calculate norm of I-th row of CECurrent. Fill by zeros, if it is
         * too small. Normalize otherwise.
         *
         * NOTE: we also scale last column of CECurrent (right part)
         */
        v = ae_v_dotproduct(&state->cecurrent.ptr.pp_double[i][0], 1, &state->cecurrent.ptr.pp_double[i][0], 1, ae_v_len(0,nmain+nslack-1));
        v = ae_sqrt(v, _state);
        if( ae_fp_greater(v,10000*ae_machineepsilon) )
        {
            v = 1/v;
            ae_v_muld(&state->cecurrent.ptr.pp_double[i][0], 1, ae_v_len(0,nmain+nslack), v);
        }
        else
        {
            for(j=0; j<=nmain+nslack; j++)
            {
                state->cecurrent.ptr.pp_double[i][j] = 0;
            }
        }
    }
    for(j=0; j<=nmain+nslack-1; j++)
    {
        state->xe.ptr.p_double[j] = 0;
    }
    for(i=0; i<=nmain+nslack-1; i++)
    {
        if( state->activeconstraints.ptr.p_bool[i] )
        {
            state->xe.ptr.p_double[i] = state->xe.ptr.p_double[i]+state->constrainedvalues.ptr.p_double[i];
        }
    }
    for(i=0; i<=state->cecnt-1; i++)
    {
        v = state->cecurrent.ptr.pp_double[i][nmain+nslack];
        ae_v_addd(&state->xe.ptr.p_double[0], 1, &state->cecurrent.ptr.pp_double[i][0], 1, ae_v_len(0,nmain+nslack-1), v);
    }
}


/*************************************************************************
This function projects gradient onto equality constrained subspace
*************************************************************************/
static void minbleic_makegradientprojection(minbleicstate* state,
     /* Real    */ ae_vector* pg,
     ae_state *) //_state)
{
    ae_int_t i;
    ae_int_t nmain;
    ae_int_t nslack;
    double v;


    nmain = state->nmain;
    nslack = state->nslack;
    for(i=0; i<=nmain+nslack-1; i++)
    {
        if( state->activeconstraints.ptr.p_bool[i] )
        {
            pg->ptr.p_double[i] = 0;
        }
    }
    for(i=0; i<=state->cecnt-1; i++)
    {
        v = ae_v_dotproduct(&pg->ptr.p_double[0], 1, &state->cecurrent.ptr.pp_double[i][0], 1, ae_v_len(0,nmain+nslack-1));
        ae_v_subd(&pg->ptr.p_double[0], 1, &state->cecurrent.ptr.pp_double[i][0], 1, ae_v_len(0,nmain+nslack-1), v);
    }
}


/*************************************************************************
This function prepares equality constrained subproblem:

1. X is used to activate constraints (if there are ones which are still
   inactive, but should be activated).
2. constraints matrix CEOrt is copied to CECurrent and modified  according
   to the list of active bound constraints (corresponding elements are
   filled by zeros and reorthogonalized).
3. XE - least squares solution of equality constraints - is recalculated
4. X is copied to PX and projected onto equality constrained subspace
5. inactive constraints are checked against PX - if there is at least one
   which should be activated, we activate it and move back to (2)
6. as result, PX is feasible with respect to bound constraints - step (5)
   guarantees it. But PX can be infeasible with respect to equality ones,
   because step (2) is done without checks for consistency. As the final
   step, we check that PX is feasible. If not, we return False. True is
   returned otherwise.

If this algorithm returned True, then:
* X is not changed
* PX contains projection of X onto constrained subspace
* G is not changed
* PG contains projection of G onto constrained subspace
* PX is feasible with respect to all constraints
* all constraints which are active at PX, are activated
*************************************************************************/
static ae_bool minbleic_prepareconstraintmatrix(minbleicstate* state,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* g,
     /* Real    */ ae_vector* px,
     /* Real    */ ae_vector* pg,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t nmain;
    ae_int_t nslack;
    double v;
    double ferr;
    ae_bool result;


    nmain = state->nmain;
    nslack = state->nslack;
    result = ae_true;
    
    /*
     * Step 1
     */
    minbleic_additionalcheckforconstraints(state, x, _state);
    
    /*
     * Steps 2-5
     */
    do
    {
        
        /*
         * Steps 2-3
         */
        minbleic_rebuildcexe(state, _state);
        
        /*
         * Step 4
         *
         * Calculate PX, PG
         */
        ae_v_move(&px->ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
        ae_v_sub(&px->ptr.p_double[0], 1, &state->xe.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
        ae_v_move(&pg->ptr.p_double[0], 1, &g->ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
        for(i=0; i<=nmain+nslack-1; i++)
        {
            if( state->activeconstraints.ptr.p_bool[i] )
            {
                px->ptr.p_double[i] = 0;
                pg->ptr.p_double[i] = 0;
            }
        }
        for(i=0; i<=state->cecnt-1; i++)
        {
            v = ae_v_dotproduct(&px->ptr.p_double[0], 1, &state->cecurrent.ptr.pp_double[i][0], 1, ae_v_len(0,nmain+nslack-1));
            ae_v_subd(&px->ptr.p_double[0], 1, &state->cecurrent.ptr.pp_double[i][0], 1, ae_v_len(0,nmain+nslack-1), v);
            v = ae_v_dotproduct(&pg->ptr.p_double[0], 1, &state->cecurrent.ptr.pp_double[i][0], 1, ae_v_len(0,nmain+nslack-1));
            ae_v_subd(&pg->ptr.p_double[0], 1, &state->cecurrent.ptr.pp_double[i][0], 1, ae_v_len(0,nmain+nslack-1), v);
        }
        ae_v_add(&px->ptr.p_double[0], 1, &state->xe.ptr.p_double[0], 1, ae_v_len(0,nmain+nslack-1));
        
        /*
         * Step 5 (loop condition below)
         */
    }
    while(minbleic_additionalcheckforconstraints(state, px, _state));
    
    /*
     * Step 6
     */
    ferr = 0;
    for(i=0; i<=state->cecnt-1; i++)
    {
        v = ae_v_dotproduct(&px->ptr.p_double[0], 1, &state->ceeffective.ptr.pp_double[i][0], 1, ae_v_len(0,nmain+nslack-1));
        v = v-state->ceeffective.ptr.pp_double[i][nmain+nslack];
        ferr = ae_maxreal(ferr, ae_fabs(v, _state), _state);
    }
    result = ae_fp_less_eq(ferr,state->outerepsi);
    return result;
}


/*************************************************************************
Internal initialization subroutine
*************************************************************************/
static void minbleic_minbleicinitinternal(ae_int_t n,
     /* Real    */ ae_vector* x,
     double diffstep,
     minbleicstate* state,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_matrix c;
    ae_vector ct;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init(&c, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&ct, 0, DT_INT, _state, ae_true);

    
    /*
     * Initialize
     */
    state->teststep = 0;
    state->nmain = n;
    state->optdim = 0;
    state->diffstep = diffstep;
    ae_vector_set_length(&state->bndloriginal, n, _state);
    ae_vector_set_length(&state->bndleffective, n, _state);
    ae_vector_set_length(&state->hasbndl, n, _state);
    ae_vector_set_length(&state->bnduoriginal, n, _state);
    ae_vector_set_length(&state->bndueffective, n, _state);
    ae_vector_set_length(&state->hasbndu, n, _state);
    ae_vector_set_length(&state->xstart, n, _state);
    ae_vector_set_length(&state->soriginal, n, _state);
    ae_vector_set_length(&state->x, n, _state);
    ae_vector_set_length(&state->g, n, _state);
    for(i=0; i<=n-1; i++)
    {
        state->bndloriginal.ptr.p_double[i] = _state->v_neginf;
        state->hasbndl.ptr.p_bool[i] = ae_false;
        state->bnduoriginal.ptr.p_double[i] = _state->v_posinf;
        state->hasbndu.ptr.p_bool[i] = ae_false;
        state->soriginal.ptr.p_double[i] = 1.0;
    }
    minbleicsetlc(state, &c, &ct, 0, _state);
    minbleicsetinnercond(state, 0.0, 0.0, 0.0, _state);
    minbleicsetoutercond(state, 1.0E-6, 1.0E-6, _state);
    minbleicsetmaxits(state, 0, _state);
    minbleicsetxrep(state, ae_false, _state);
    minbleicsetstpmax(state, 0.0, _state);
    minbleicsetprecdefault(state, _state);
    minbleicrestartfrom(state, x, _state);
    ae_frame_leave(_state);
}


ae_bool _minbleicstate_init(minbleicstate* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->diaghoriginal, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->diagh, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->x, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->g, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !_rcommstate_init(&p->rstate, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xcur, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xprev, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xstart, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xend, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->lastg, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->ceoriginal, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->ceeffective, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->cecurrent, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->ct, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xe, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->hasbndl, 0, DT_BOOL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->hasbndu, 0, DT_BOOL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->bndloriginal, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->bnduoriginal, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->bndleffective, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->bndueffective, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->activeconstraints, 0, DT_BOOL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->constrainedvalues, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->transforms, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->seffective, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->soriginal, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->w, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tmp0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tmp1, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tmp2, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->r, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->lmmatrix, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !_mincgstate_init(&p->cgstate, _state, make_automatic) )
        return ae_false;
    if( !_mincgreport_init(&p->cgrep, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _minbleicstate_init_copy(minbleicstate* dst, minbleicstate* src, ae_state *_state, ae_bool make_automatic)
{
    dst->nmain = src->nmain;
    dst->nslack = src->nslack;
    dst->innerepsg = src->innerepsg;
    dst->innerepsf = src->innerepsf;
    dst->innerepsx = src->innerepsx;
    dst->outerepsx = src->outerepsx;
    dst->outerepsi = src->outerepsi;
    dst->maxits = src->maxits;
    dst->xrep = src->xrep;
    dst->stpmax = src->stpmax;
    dst->diffstep = src->diffstep;
    dst->prectype = src->prectype;
    if( !ae_vector_init_copy(&dst->diaghoriginal, &src->diaghoriginal, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->diagh, &src->diagh, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->x, &src->x, _state, make_automatic) )
        return ae_false;
    dst->f = src->f;
    if( !ae_vector_init_copy(&dst->g, &src->g, _state, make_automatic) )
        return ae_false;
    dst->needf = src->needf;
    dst->needfg = src->needfg;
    dst->xupdated = src->xupdated;
    dst->teststep = src->teststep;
    if( !_rcommstate_init_copy(&dst->rstate, &src->rstate, _state, make_automatic) )
        return ae_false;
    dst->repinneriterationscount = src->repinneriterationscount;
    dst->repouteriterationscount = src->repouteriterationscount;
    dst->repnfev = src->repnfev;
    dst->repvaridx = src->repvaridx;
    dst->repterminationtype = src->repterminationtype;
    dst->repdebugeqerr = src->repdebugeqerr;
    dst->repdebugfs = src->repdebugfs;
    dst->repdebugff = src->repdebugff;
    dst->repdebugdx = src->repdebugdx;
    dst->repdebugfeasqpits = src->repdebugfeasqpits;
    dst->repdebugfeasgpaits = src->repdebugfeasgpaits;
    if( !ae_vector_init_copy(&dst->xcur, &src->xcur, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->xprev, &src->xprev, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->xstart, &src->xstart, _state, make_automatic) )
        return ae_false;
    dst->itsleft = src->itsleft;
    if( !ae_vector_init_copy(&dst->xend, &src->xend, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->lastg, &src->lastg, _state, make_automatic) )
        return ae_false;
    dst->trimthreshold = src->trimthreshold;
    if( !ae_matrix_init_copy(&dst->ceoriginal, &src->ceoriginal, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->ceeffective, &src->ceeffective, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->cecurrent, &src->cecurrent, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->ct, &src->ct, _state, make_automatic) )
        return ae_false;
    dst->cecnt = src->cecnt;
    dst->cedim = src->cedim;
    if( !ae_vector_init_copy(&dst->xe, &src->xe, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->hasbndl, &src->hasbndl, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->hasbndu, &src->hasbndu, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->bndloriginal, &src->bndloriginal, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->bnduoriginal, &src->bnduoriginal, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->bndleffective, &src->bndleffective, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->bndueffective, &src->bndueffective, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->activeconstraints, &src->activeconstraints, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->constrainedvalues, &src->constrainedvalues, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->transforms, &src->transforms, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->seffective, &src->seffective, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->soriginal, &src->soriginal, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->w, &src->w, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tmp0, &src->tmp0, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tmp1, &src->tmp1, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tmp2, &src->tmp2, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->r, &src->r, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->lmmatrix, &src->lmmatrix, _state, make_automatic) )
        return ae_false;
    dst->errfeas = src->errfeas;
    dst->gnorm = src->gnorm;
    dst->mpgnorm = src->mpgnorm;
    dst->mba = src->mba;
    dst->variabletofreeze = src->variabletofreeze;
    dst->valuetofreeze = src->valuetofreeze;
    dst->fbase = src->fbase;
    dst->fm2 = src->fm2;
    dst->fm1 = src->fm1;
    dst->fp1 = src->fp1;
    dst->fp2 = src->fp2;
    dst->xm1 = src->xm1;
    dst->xp1 = src->xp1;
    if( !_mincgstate_init_copy(&dst->cgstate, &src->cgstate, _state, make_automatic) )
        return ae_false;
    if( !_mincgreport_init_copy(&dst->cgrep, &src->cgrep, _state, make_automatic) )
        return ae_false;
    dst->optdim = src->optdim;
    return ae_true;
}


void _minbleicstate_clear(minbleicstate* p)
{
    ae_vector_clear(&p->diaghoriginal);
    ae_vector_clear(&p->diagh);
    ae_vector_clear(&p->x);
    ae_vector_clear(&p->g);
    _rcommstate_clear(&p->rstate);
    ae_vector_clear(&p->xcur);
    ae_vector_clear(&p->xprev);
    ae_vector_clear(&p->xstart);
    ae_vector_clear(&p->xend);
    ae_vector_clear(&p->lastg);
    ae_matrix_clear(&p->ceoriginal);
    ae_matrix_clear(&p->ceeffective);
    ae_matrix_clear(&p->cecurrent);
    ae_vector_clear(&p->ct);
    ae_vector_clear(&p->xe);
    ae_vector_clear(&p->hasbndl);
    ae_vector_clear(&p->hasbndu);
    ae_vector_clear(&p->bndloriginal);
    ae_vector_clear(&p->bnduoriginal);
    ae_vector_clear(&p->bndleffective);
    ae_vector_clear(&p->bndueffective);
    ae_vector_clear(&p->activeconstraints);
    ae_vector_clear(&p->constrainedvalues);
    ae_vector_clear(&p->transforms);
    ae_vector_clear(&p->seffective);
    ae_vector_clear(&p->soriginal);
    ae_vector_clear(&p->w);
    ae_vector_clear(&p->tmp0);
    ae_vector_clear(&p->tmp1);
    ae_vector_clear(&p->tmp2);
    ae_vector_clear(&p->r);
    ae_matrix_clear(&p->lmmatrix);
    _mincgstate_clear(&p->cgstate);
    _mincgreport_clear(&p->cgrep);
}


ae_bool _minbleicreport_init(minbleicreport* /*p*/, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    return ae_true;
}


ae_bool _minbleicreport_init_copy(minbleicreport* dst, minbleicreport* src, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    dst->inneriterationscount = src->inneriterationscount;
    dst->outeriterationscount = src->outeriterationscount;
    dst->nfev = src->nfev;
    dst->varidx = src->varidx;
    dst->terminationtype = src->terminationtype;
    dst->debugeqerr = src->debugeqerr;
    dst->debugfs = src->debugfs;
    dst->debugff = src->debugff;
    dst->debugdx = src->debugdx;
    dst->debugfeasqpits = src->debugfeasqpits;
    dst->debugfeasgpaits = src->debugfeasgpaits;
    return ae_true;
}


void _minbleicreport_clear(minbleicreport*) // p)
{
}




/*************************************************************************
        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION

DESCRIPTION:
The subroutine minimizes function F(x) of N arguments by  using  a  quasi-
Newton method (LBFGS scheme) which is optimized to use  a  minimum  amount
of memory.
The subroutine generates the approximation of an inverse Hessian matrix by
using information about the last M steps of the algorithm  (instead of N).
It lessens a required amount of memory from a value  of  order  N^2  to  a
value of order 2*N*M.


REQUIREMENTS:
Algorithm will request following information during its operation:
* function value F and its gradient G (simultaneously) at given point X


USAGE:
1. User initializes algorithm state with MinLBFGSCreate() call
2. User tunes solver parameters with MinLBFGSSetCond() MinLBFGSSetStpMax()
   and other functions
3. User calls MinLBFGSOptimize() function which takes algorithm  state and
   pointer (delegate, etc.) to callback function which calculates F/G.
4. User calls MinLBFGSResults() to get solution
5. Optionally user may call MinLBFGSRestartFrom() to solve another problem
   with same N/M but another starting point and/or another function.
   MinLBFGSRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   problem dimension. N>0
    M       -   number of corrections in the BFGS scheme of Hessian
                approximation update. Recommended value:  3<=M<=7. The smaller
                value causes worse convergence, the bigger will  not  cause  a
                considerably better convergence, but will cause a fall in  the
                performance. M<=N.
    X       -   initial solution approximation, array[0..N-1].


OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state
    

NOTES:
1. you may tune stopping conditions with MinLBFGSSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLBFGSSetStpMax() function to bound algorithm's  steps.  However,
   L-BFGS rarely needs such a tuning.


  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgscreate(ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_vector* x,
     minlbfgsstate* state,
     ae_state *_state)
{

    _minlbfgsstate_clear(state);

    ae_assert(n>=1, "MinLBFGSCreate: N<1!", _state);
    ae_assert(m>=1, "MinLBFGSCreate: M<1", _state);
    ae_assert(m<=n, "MinLBFGSCreate: M>N", _state);
    ae_assert(x->cnt>=n, "MinLBFGSCreate: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "MinLBFGSCreate: X contains infinite or NaN values!", _state);
    minlbfgscreatex(n, m, x, 0, 0.0, state, _state);
}


/*************************************************************************
The subroutine is finite difference variant of MinLBFGSCreate().  It  uses
finite differences in order to differentiate target function.

Description below contains information which is specific to  this function
only. We recommend to read comments on MinLBFGSCreate() in  order  to  get
more information about creation of LBFGS optimizer.

INPUT PARAMETERS:
    N       -   problem dimension, N>0:
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    M       -   number of corrections in the BFGS scheme of Hessian
                approximation update. Recommended value:  3<=M<=7. The smaller
                value causes worse convergence, the bigger will  not  cause  a
                considerably better convergence, but will cause a fall in  the
                performance. M<=N.
    X       -   starting point, array[0..N-1].
    DiffStep-   differentiation step, >0

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTES:
1. algorithm uses 4-point central formula for differentiation.
2. differentiation step along I-th axis is equal to DiffStep*S[I] where
   S[] is scaling vector which can be set by MinLBFGSSetScale() call.
3. we recommend you to use moderate values of  differentiation  step.  Too
   large step will result in too large truncation  errors, while too small
   step will result in too large numerical  errors.  1.0E-6  can  be  good
   value to start with.
4. Numerical  differentiation  is   very   inefficient  -   one   gradient
   calculation needs 4*N function evaluations. This function will work for
   any N - either small (1...10), moderate (10...100) or  large  (100...).
   However, performance penalty will be too severe for any N's except  for
   small ones.
   We should also say that code which relies on numerical  differentiation
   is   less  robust  and  precise.  LBFGS  needs  exact  gradient values.
   Imprecise gradient may slow  down  convergence,  especially  on  highly
   nonlinear problems.
   Thus  we  recommend to use this function for fast prototyping on small-
   dimensional problems only, and to implement analytical gradient as soon
   as possible.

  -- ALGLIB --
     Copyright 16.05.2011 by Bochkanov Sergey
*************************************************************************/
void minlbfgscreatef(ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_vector* x,
     double diffstep,
     minlbfgsstate* state,
     ae_state *_state)
{

    _minlbfgsstate_clear(state);

    ae_assert(n>=1, "MinLBFGSCreateF: N too small!", _state);
    ae_assert(m>=1, "MinLBFGSCreateF: M<1", _state);
    ae_assert(m<=n, "MinLBFGSCreateF: M>N", _state);
    ae_assert(x->cnt>=n, "MinLBFGSCreateF: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "MinLBFGSCreateF: X contains infinite or NaN values!", _state);
    ae_assert(ae_isfinite(diffstep, _state), "MinLBFGSCreateF: DiffStep is infinite or NaN!", _state);
    ae_assert(ae_fp_greater(diffstep,0), "MinLBFGSCreateF: DiffStep is non-positive!", _state);
    minlbfgscreatex(n, m, x, 0, diffstep, state, _state);
}


/*************************************************************************
This function sets stopping conditions for L-BFGS optimization algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    EpsG    -   >=0
                The  subroutine  finishes  its  work   if   the  condition
                |v|<EpsG is satisfied, where:
                * |.| means Euclidian norm
                * v - scaled gradient vector, v[i]=g[i]*s[i]
                * g - gradient
                * s - scaling coefficients set by MinLBFGSSetScale()
    EpsF    -   >=0
                The  subroutine  finishes  its work if on k+1-th iteration
                the  condition  |F(k+1)-F(k)|<=EpsF*max{|F(k)|,|F(k+1)|,1}
                is satisfied.
    EpsX    -   >=0
                The subroutine finishes its work if  on  k+1-th  iteration
                the condition |v|<=EpsX is fulfilled, where:
                * |.| means Euclidian norm
                * v - scaled step vector, v[i]=dx[i]/s[i]
                * dx - ste pvector, dx=X(k+1)-X(k)
                * s - scaling coefficients set by MinLBFGSSetScale()
    MaxIts  -   maximum number of iterations. If MaxIts=0, the  number  of
                iterations is unlimited.

Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to
automatic stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetcond(minlbfgsstate* state,
     double epsg,
     double epsf,
     double epsx,
     ae_int_t maxits,
     ae_state *_state)
{


    ae_assert(ae_isfinite(epsg, _state), "MinLBFGSSetCond: EpsG is not finite number!", _state);
    ae_assert(ae_fp_greater_eq(epsg,0), "MinLBFGSSetCond: negative EpsG!", _state);
    ae_assert(ae_isfinite(epsf, _state), "MinLBFGSSetCond: EpsF is not finite number!", _state);
    ae_assert(ae_fp_greater_eq(epsf,0), "MinLBFGSSetCond: negative EpsF!", _state);
    ae_assert(ae_isfinite(epsx, _state), "MinLBFGSSetCond: EpsX is not finite number!", _state);
    ae_assert(ae_fp_greater_eq(epsx,0), "MinLBFGSSetCond: negative EpsX!", _state);
    ae_assert(maxits>=0, "MinLBFGSSetCond: negative MaxIts!", _state);
    if( ((ae_fp_eq(epsg,0)&&ae_fp_eq(epsf,0))&&ae_fp_eq(epsx,0))&&maxits==0 )
    {
        epsx = 1.0E-6;
    }
    state->epsg = epsg;
    state->epsf = epsf;
    state->epsx = epsx;
    state->maxits = maxits;
}


/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    NeedXRep-   whether iteration reports are needed or not

If NeedXRep is True, algorithm will call rep() callback function if  it is
provided to MinLBFGSOptimize().


  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetxrep(minlbfgsstate* state,
     ae_bool needxrep,
     ae_state *) //_state)
{


    state->xrep = needxrep;
}


/*************************************************************************
This function sets maximum step length

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    StpMax  -   maximum step length, >=0. Set StpMax to 0.0 (default),  if
                you don't want to limit step length.

Use this subroutine when you optimize target function which contains exp()
or  other  fast  growing  functions,  and optimization algorithm makes too
large  steps  which  leads  to overflow. This function allows us to reject
steps  that  are  too  large  (and  therefore  expose  us  to the possible
overflow) without actually calculating function value at the x+stp*d.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetstpmax(minlbfgsstate* state,
     double stpmax,
     ae_state *_state)
{


    ae_assert(ae_isfinite(stpmax, _state), "MinLBFGSSetStpMax: StpMax is not finite!", _state);
    ae_assert(ae_fp_greater_eq(stpmax,0), "MinLBFGSSetStpMax: StpMax<0!", _state);
    state->stpmax = stpmax;
}


/*************************************************************************
This function sets scaling coefficients for LBFGS optimizer.

ALGLIB optimizers use scaling matrices to test stopping  conditions  (step
size and gradient are scaled before comparison with tolerances).  Scale of
the I-th variable is a translation invariant measure of:
a) "how large" the variable is
b) how large the step should be to make significant changes in the function

Scaling is also used by finite difference variant of the optimizer  - step
along I-th axis is equal to DiffStep*S[I].

In  most  optimizers  (and  in  the  LBFGS  too)  scaling is NOT a form of
preconditioning. It just  affects  stopping  conditions.  You  should  set
preconditioner  by  separate  call  to  one  of  the  MinLBFGSSetPrec...()
functions.

There  is  special  preconditioning  mode, however,  which  uses   scaling
coefficients to form diagonal preconditioning matrix. You  can  turn  this
mode on, if you want.   But  you should understand that scaling is not the
same thing as preconditioning - these are two different, although  related
forms of tuning solver.

INPUT PARAMETERS:
    State   -   structure stores algorithm state
    S       -   array[N], non-zero scaling coefficients
                S[i] may be negative, sign doesn't matter.

  -- ALGLIB --
     Copyright 14.01.2011 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetscale(minlbfgsstate* state,
     /* Real    */ ae_vector* s,
     ae_state *_state)
{
    ae_int_t i;


    ae_assert(s->cnt>=state->n, "MinLBFGSSetScale: Length(S)<N", _state);
    for(i=0; i<=state->n-1; i++)
    {
        ae_assert(ae_isfinite(s->ptr.p_double[i], _state), "MinLBFGSSetScale: S contains infinite or NAN elements", _state);
        ae_assert(ae_fp_neq(s->ptr.p_double[i],0), "MinLBFGSSetScale: S contains zero elements", _state);
        state->s.ptr.p_double[i] = ae_fabs(s->ptr.p_double[i], _state);
    }
}


/*************************************************************************
Extended subroutine for internal use only.

Accepts additional parameters:

    Flags - additional settings:
            * Flags = 0     means no additional settings
            * Flags = 1     "do not allocate memory". used when solving
                            a many subsequent tasks with  same N/M  values.
                            First  call MUST  be without this flag bit set,
                            subsequent  calls   of   MinLBFGS   with   same
                            MinLBFGSState structure can set Flags to 1.
    DiffStep - numerical differentiation step

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgscreatex(ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_vector* x,
     ae_int_t flags,
     double diffstep,
     minlbfgsstate* state,
     ae_state *_state)
{
    ae_bool allocatemem;
    ae_int_t i;


    ae_assert(n>=1, "MinLBFGS: N too small!", _state);
    ae_assert(m>=1, "MinLBFGS: M too small!", _state);
    ae_assert(m<=n, "MinLBFGS: M too large!", _state);
    
    /*
     * Initialize
     */
    state->teststep = 0;
    state->diffstep = diffstep;
    state->n = n;
    state->m = m;
    allocatemem = flags%2==0;
    flags = flags/2;
    if( allocatemem )
    {
        ae_vector_set_length(&state->rho, m, _state);
        ae_vector_set_length(&state->theta, m, _state);
        ae_matrix_set_length(&state->yk, m, n, _state);
        ae_matrix_set_length(&state->sk, m, n, _state);
        ae_vector_set_length(&state->d, n, _state);
        ae_vector_set_length(&state->x, n, _state);
        ae_vector_set_length(&state->s, n, _state);
        ae_vector_set_length(&state->g, n, _state);
        ae_vector_set_length(&state->work, n, _state);
    }
    minlbfgssetcond(state, 0, 0, 0, 0, _state);
    minlbfgssetxrep(state, ae_false, _state);
    minlbfgssetstpmax(state, 0, _state);
    minlbfgsrestartfrom(state, x, _state);
    for(i=0; i<=n-1; i++)
    {
        state->s.ptr.p_double[i] = 1.0;
    }
    state->prectype = 0;
}


/*************************************************************************
Modification  of  the  preconditioner:  default  preconditioner    (simple
scaling, same for all elements of X) is used.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTE:  you  can  change  preconditioner  "on  the  fly",  during algorithm
iterations.

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetprecdefault(minlbfgsstate* state, ae_state *) //_state)
{


    state->prectype = 0;
}


/*************************************************************************
Modification of the preconditioner: Cholesky factorization of  approximate
Hessian is used.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    P       -   triangular preconditioner, Cholesky factorization of
                the approximate Hessian. array[0..N-1,0..N-1],
                (if larger, only leading N elements are used).
    IsUpper -   whether upper or lower triangle of P is given
                (other triangle is not referenced)

After call to this function preconditioner is changed to P  (P  is  copied
into the internal buffer).

NOTE:  you  can  change  preconditioner  "on  the  fly",  during algorithm
iterations.

NOTE 2:  P  should  be nonsingular. Exception will be thrown otherwise.

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetpreccholesky(minlbfgsstate* state,
     /* Real    */ ae_matrix* p,
     ae_bool isupper,
     ae_state *_state)
{
    ae_int_t i;
    double mx;


    ae_assert(isfinitertrmatrix(p, state->n, isupper, _state), "MinLBFGSSetPrecCholesky: P contains infinite or NAN values!", _state);
    mx = 0;
    for(i=0; i<=state->n-1; i++)
    {
        mx = ae_maxreal(mx, ae_fabs(p->ptr.pp_double[i][i], _state), _state);
    }
    ae_assert(ae_fp_greater(mx,0), "MinLBFGSSetPrecCholesky: P is strictly singular!", _state);
    if( state->denseh.rows<state->n||state->denseh.cols<state->n )
    {
        ae_matrix_set_length(&state->denseh, state->n, state->n, _state);
    }
    state->prectype = 1;
    if( isupper )
    {
        rmatrixcopy(state->n, state->n, p, 0, 0, &state->denseh, 0, 0, _state);
    }
    else
    {
        rmatrixtranspose(state->n, state->n, p, 0, 0, &state->denseh, 0, 0, _state);
    }
}


/*************************************************************************
Modification  of  the  preconditioner:  diagonal of approximate Hessian is
used.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    D       -   diagonal of the approximate Hessian, array[0..N-1],
                (if larger, only leading N elements are used).

NOTE:  you  can  change  preconditioner  "on  the  fly",  during algorithm
iterations.

NOTE 2: D[i] should be positive. Exception will be thrown otherwise.

NOTE 3: you should pass diagonal of approximate Hessian - NOT ITS INVERSE.

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetprecdiag(minlbfgsstate* state,
     /* Real    */ ae_vector* d,
     ae_state *_state)
{
    ae_int_t i;


    ae_assert(d->cnt>=state->n, "MinLBFGSSetPrecDiag: D is too short", _state);
    for(i=0; i<=state->n-1; i++)
    {
        ae_assert(ae_isfinite(d->ptr.p_double[i], _state), "MinLBFGSSetPrecDiag: D contains infinite or NAN elements", _state);
        ae_assert(ae_fp_greater(d->ptr.p_double[i],0), "MinLBFGSSetPrecDiag: D contains non-positive elements", _state);
    }
    rvectorsetlengthatleast(&state->diagh, state->n, _state);
    state->prectype = 2;
    for(i=0; i<=state->n-1; i++)
    {
        state->diagh.ptr.p_double[i] = d->ptr.p_double[i];
    }
}


/*************************************************************************
Modification of the preconditioner: scale-based diagonal preconditioning.

This preconditioning mode can be useful when you  don't  have  approximate
diagonal of Hessian, but you know that your  variables  are  badly  scaled
(for  example,  one  variable is in [1,10], and another in [1000,100000]),
and most part of the ill-conditioning comes from different scales of vars.

In this case simple  scale-based  preconditioner,  with H[i] = 1/(s[i]^2),
can greatly improve convergence.

IMPRTANT: you should set scale of your variables  with  MinLBFGSSetScale()
call  (before  or after MinLBFGSSetPrecScale() call). Without knowledge of
the scale of your variables scale-based preconditioner will be  just  unit
matrix.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetprecscale(minlbfgsstate* state, ae_state *) //_state)
{


    state->prectype = 3;
}


/*************************************************************************
NOTES:

1. This function has two different implementations: one which  uses  exact
   (analytical) user-supplied gradient,  and one which uses function value
   only  and  numerically  differentiates  function  in  order  to  obtain
   gradient.

   Depending  on  the  specific  function  used to create optimizer object
   (either MinLBFGSCreate() for analytical gradient  or  MinLBFGSCreateF()
   for numerical differentiation) you should choose appropriate variant of
   MinLBFGSOptimize() - one  which  accepts  function  AND gradient or one
   which accepts function ONLY.

   Be careful to choose variant of MinLBFGSOptimize() which corresponds to
   your optimization scheme! Table below lists different  combinations  of
   callback (function/gradient) passed to MinLBFGSOptimize()  and specific
   function used to create optimizer.


                     |         USER PASSED TO MinLBFGSOptimize()
   CREATED WITH      |  function only   |  function and gradient
   ------------------------------------------------------------
   MinLBFGSCreateF() |     work                FAIL
   MinLBFGSCreate()  |     FAIL                work

   Here "FAIL" denotes inappropriate combinations  of  optimizer  creation
   function  and  MinLBFGSOptimize()  version.   Attemps   to   use   such
   combination (for example, to create optimizer with MinLBFGSCreateF() and
   to pass gradient information to MinCGOptimize()) will lead to exception
   being thrown. Either  you  did  not pass gradient when it WAS needed or
   you passed gradient when it was NOT needed.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool minlbfgsiteration(minlbfgsstate* state, ae_state *_state)
{
    ae_int_t n;
    ae_int_t m;
    ae_int_t i;
    ae_int_t j;
    ae_int_t ic;
    ae_int_t mcinfo;
    double v;
    double vv;
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
        m = state->rstate.ia.ptr.p_int[1];
        i = state->rstate.ia.ptr.p_int[2];
        j = state->rstate.ia.ptr.p_int[3];
        ic = state->rstate.ia.ptr.p_int[4];
        mcinfo = state->rstate.ia.ptr.p_int[5];
        v = state->rstate.ra.ptr.p_double[0];
        vv = state->rstate.ra.ptr.p_double[1];
    }
    else
    {
        n = -983;
        m = -989;
        i = -834;
        j = 900;
        ic = -287;
        mcinfo = 364;
        v = 214;
        vv = -338;
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
    if( state->rstate.stage==4 )
    {
        goto lbl_4;
    }
    if( state->rstate.stage==5 )
    {
        goto lbl_5;
    }
    if( state->rstate.stage==6 )
    {
        goto lbl_6;
    }
    if( state->rstate.stage==7 )
    {
        goto lbl_7;
    }
    if( state->rstate.stage==8 )
    {
        goto lbl_8;
    }
    if( state->rstate.stage==9 )
    {
        goto lbl_9;
    }
    if( state->rstate.stage==10 )
    {
        goto lbl_10;
    }
    if( state->rstate.stage==11 )
    {
        goto lbl_11;
    }
    if( state->rstate.stage==12 )
    {
        goto lbl_12;
    }
    if( state->rstate.stage==13 )
    {
        goto lbl_13;
    }
    if( state->rstate.stage==14 )
    {
        goto lbl_14;
    }
    if( state->rstate.stage==15 )
    {
        goto lbl_15;
    }
    if( state->rstate.stage==16 )
    {
        goto lbl_16;
    }
    
    /*
     * Routine body
     */
    
    /*
     * Unload frequently used variables from State structure
     * (just for typing convinience)
     */
    n = state->n;
    m = state->m;
    state->repterminationtype = 0;
    state->repiterationscount = 0;
    state->repvaridx = -1;
    state->repnfev = 0;
    
    /*
     *  Check, that transferred derivative value is right
     */
    minlbfgs_clearrequestfields(state, _state);
    if( !(ae_fp_eq(state->diffstep,0)&&ae_fp_greater(state->teststep,0)) )
    {
        goto lbl_17;
    }
    state->needfg = ae_true;
    i = 0;
lbl_19:
    if( i>n-1 )
    {
        goto lbl_21;
    }
    v = state->x.ptr.p_double[i];
    state->x.ptr.p_double[i] = v-state->teststep*state->s.ptr.p_double[i];
    state->rstate.stage = 0;
    goto lbl_rcomm;
lbl_0:
    state->fm1 = state->f;
    state->fp1 = state->g.ptr.p_double[i];
    state->x.ptr.p_double[i] = v+state->teststep*state->s.ptr.p_double[i];
    state->rstate.stage = 1;
    goto lbl_rcomm;
lbl_1:
    state->fm2 = state->f;
    state->fp2 = state->g.ptr.p_double[i];
    state->x.ptr.p_double[i] = v;
    state->rstate.stage = 2;
    goto lbl_rcomm;
lbl_2:
    
    /*
     * 2*State.TestStep   -   scale parameter
     * width of segment [Xi-TestStep;Xi+TestStep]
     */
    if( !derivativecheck(state->fm1, state->fp1, state->fm2, state->fp2, state->f, state->g.ptr.p_double[i], 2*state->teststep, _state) )
    {
        state->repvaridx = i;
        state->repterminationtype = -7;
        result = ae_false;
        return result;
    }
    i = i+1;
    goto lbl_19;
lbl_21:
    state->needfg = ae_false;
lbl_17:
    
    /*
     * Calculate F/G at the initial point
     */
    minlbfgs_clearrequestfields(state, _state);
    if( ae_fp_neq(state->diffstep,0) )
    {
        goto lbl_22;
    }
    state->needfg = ae_true;
    state->rstate.stage = 3;
    goto lbl_rcomm;
lbl_3:
    state->needfg = ae_false;
    goto lbl_23;
lbl_22:
    state->needf = ae_true;
    state->rstate.stage = 4;
    goto lbl_rcomm;
lbl_4:
    state->fbase = state->f;
    i = 0;
lbl_24:
    if( i>n-1 )
    {
        goto lbl_26;
    }
    v = state->x.ptr.p_double[i];
    state->x.ptr.p_double[i] = v-state->diffstep*state->s.ptr.p_double[i];
    state->rstate.stage = 5;
    goto lbl_rcomm;
lbl_5:
    state->fm2 = state->f;
    state->x.ptr.p_double[i] = v-0.5*state->diffstep*state->s.ptr.p_double[i];
    state->rstate.stage = 6;
    goto lbl_rcomm;
lbl_6:
    state->fm1 = state->f;
    state->x.ptr.p_double[i] = v+0.5*state->diffstep*state->s.ptr.p_double[i];
    state->rstate.stage = 7;
    goto lbl_rcomm;
lbl_7:
    state->fp1 = state->f;
    state->x.ptr.p_double[i] = v+state->diffstep*state->s.ptr.p_double[i];
    state->rstate.stage = 8;
    goto lbl_rcomm;
lbl_8:
    state->fp2 = state->f;
    state->x.ptr.p_double[i] = v;
    state->g.ptr.p_double[i] = (8*(state->fp1-state->fm1)-(state->fp2-state->fm2))/(6*state->diffstep*state->s.ptr.p_double[i]);
    i = i+1;
    goto lbl_24;
lbl_26:
    state->f = state->fbase;
    state->needf = ae_false;
lbl_23:
    trimprepare(state->f, &state->trimthreshold, _state);
    if( !state->xrep )
    {
        goto lbl_27;
    }
    minlbfgs_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 9;
    goto lbl_rcomm;
lbl_9:
    state->xupdated = ae_false;
lbl_27:
    state->repnfev = 1;
    state->fold = state->f;
    v = 0;
    for(i=0; i<=n-1; i++)
    {
        v = v+ae_sqr(state->g.ptr.p_double[i]*state->s.ptr.p_double[i], _state);
    }
    if( ae_fp_less_eq(ae_sqrt(v, _state),state->epsg) )
    {
        state->repterminationtype = 4;
        result = ae_false;
        return result;
    }
    
    /*
     * Choose initial step and direction.
     * Apply preconditioner, if we have something other than default.
     */
    ae_v_moveneg(&state->d.ptr.p_double[0], 1, &state->g.ptr.p_double[0], 1, ae_v_len(0,n-1));
    if( state->prectype==0 )
    {
        
        /*
         * Default preconditioner is used, but we can't use it before iterations will start
         */
        v = ae_v_dotproduct(&state->g.ptr.p_double[0], 1, &state->g.ptr.p_double[0], 1, ae_v_len(0,n-1));
        v = ae_sqrt(v, _state);
        if( ae_fp_eq(state->stpmax,0) )
        {
            state->stp = ae_minreal(1.0/v, 1, _state);
        }
        else
        {
            state->stp = ae_minreal(1.0/v, state->stpmax, _state);
        }
    }
    if( state->prectype==1 )
    {
        
        /*
         * Cholesky preconditioner is used
         */
        fblscholeskysolve(&state->denseh, 1.0, n, ae_true, &state->d, &state->autobuf, _state);
        state->stp = 1;
    }
    if( state->prectype==2 )
    {
        
        /*
         * diagonal approximation is used
         */
        for(i=0; i<=n-1; i++)
        {
            state->d.ptr.p_double[i] = state->d.ptr.p_double[i]/state->diagh.ptr.p_double[i];
        }
        state->stp = 1;
    }
    if( state->prectype==3 )
    {
        
        /*
         * scale-based preconditioner is used
         */
        for(i=0; i<=n-1; i++)
        {
            state->d.ptr.p_double[i] = state->d.ptr.p_double[i]*state->s.ptr.p_double[i]*state->s.ptr.p_double[i];
        }
        state->stp = 1;
    }
    
    /*
     * Main cycle
     */
    state->k = 0;
lbl_29:
    if( ae_false )
    {
        goto lbl_30;
    }
    
    /*
     * Main cycle: prepare to 1-D line search
     */
    state->p = state->k%m;
    state->q = ae_minint(state->k, m-1, _state);
    
    /*
     * Store X[k], G[k]
     */
    ae_v_moveneg(&state->sk.ptr.pp_double[state->p][0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_moveneg(&state->yk.ptr.pp_double[state->p][0], 1, &state->g.ptr.p_double[0], 1, ae_v_len(0,n-1));
    
    /*
     * Minimize F(x+alpha*d)
     * Calculate S[k], Y[k]
     */
    state->mcstage = 0;
    if( state->k!=0 )
    {
        state->stp = 1.0;
    }
    linminnormalized(&state->d, &state->stp, n, _state);
    mcsrch(n, &state->x, &state->f, &state->g, &state->d, &state->stp, state->stpmax, minlbfgs_gtol, &mcinfo, &state->nfev, &state->work, &state->lstate, &state->mcstage, _state);
lbl_31:
    if( state->mcstage==0 )
    {
        goto lbl_32;
    }
    minlbfgs_clearrequestfields(state, _state);
    if( ae_fp_neq(state->diffstep,0) )
    {
        goto lbl_33;
    }
    state->needfg = ae_true;
    state->rstate.stage = 10;
    goto lbl_rcomm;
lbl_10:
    state->needfg = ae_false;
    goto lbl_34;
lbl_33:
    state->needf = ae_true;
    state->rstate.stage = 11;
    goto lbl_rcomm;
lbl_11:
    state->fbase = state->f;
    i = 0;
lbl_35:
    if( i>n-1 )
    {
        goto lbl_37;
    }
    v = state->x.ptr.p_double[i];
    state->x.ptr.p_double[i] = v-state->diffstep*state->s.ptr.p_double[i];
    state->rstate.stage = 12;
    goto lbl_rcomm;
lbl_12:
    state->fm2 = state->f;
    state->x.ptr.p_double[i] = v-0.5*state->diffstep*state->s.ptr.p_double[i];
    state->rstate.stage = 13;
    goto lbl_rcomm;
lbl_13:
    state->fm1 = state->f;
    state->x.ptr.p_double[i] = v+0.5*state->diffstep*state->s.ptr.p_double[i];
    state->rstate.stage = 14;
    goto lbl_rcomm;
lbl_14:
    state->fp1 = state->f;
    state->x.ptr.p_double[i] = v+state->diffstep*state->s.ptr.p_double[i];
    state->rstate.stage = 15;
    goto lbl_rcomm;
lbl_15:
    state->fp2 = state->f;
    state->x.ptr.p_double[i] = v;
    state->g.ptr.p_double[i] = (8*(state->fp1-state->fm1)-(state->fp2-state->fm2))/(6*state->diffstep*state->s.ptr.p_double[i]);
    i = i+1;
    goto lbl_35;
lbl_37:
    state->f = state->fbase;
    state->needf = ae_false;
lbl_34:
    trimfunction(&state->f, &state->g, n, state->trimthreshold, _state);
    mcsrch(n, &state->x, &state->f, &state->g, &state->d, &state->stp, state->stpmax, minlbfgs_gtol, &mcinfo, &state->nfev, &state->work, &state->lstate, &state->mcstage, _state);
    goto lbl_31;
lbl_32:
    if( !state->xrep )
    {
        goto lbl_38;
    }
    
    /*
     * report
     */
    minlbfgs_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 16;
    goto lbl_rcomm;
lbl_16:
    state->xupdated = ae_false;
lbl_38:
    state->repnfev = state->repnfev+state->nfev;
    state->repiterationscount = state->repiterationscount+1;
    ae_v_add(&state->sk.ptr.pp_double[state->p][0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_add(&state->yk.ptr.pp_double[state->p][0], 1, &state->g.ptr.p_double[0], 1, ae_v_len(0,n-1));
    
    /*
     * Stopping conditions
     */
    if( state->repiterationscount>=state->maxits&&state->maxits>0 )
    {
        
        /*
         * Too many iterations
         */
        state->repterminationtype = 5;
        result = ae_false;
        return result;
    }
    v = 0;
    for(i=0; i<=n-1; i++)
    {
        v = v+ae_sqr(state->g.ptr.p_double[i]*state->s.ptr.p_double[i], _state);
    }
    if( ae_fp_less_eq(ae_sqrt(v, _state),state->epsg) )
    {
        
        /*
         * Gradient is small enough
         */
        state->repterminationtype = 4;
        result = ae_false;
        return result;
    }
    if( ae_fp_less_eq(state->fold-state->f,state->epsf*ae_maxreal(ae_fabs(state->fold, _state), ae_maxreal(ae_fabs(state->f, _state), 1.0, _state), _state)) )
    {
        
        /*
         * F(k+1)-F(k) is small enough
         */
        state->repterminationtype = 1;
        result = ae_false;
        return result;
    }
    v = 0;
    for(i=0; i<=n-1; i++)
    {
        v = v+ae_sqr(state->sk.ptr.pp_double[state->p][i]/state->s.ptr.p_double[i], _state);
    }
    if( ae_fp_less_eq(ae_sqrt(v, _state),state->epsx) )
    {
        
        /*
         * X(k+1)-X(k) is small enough
         */
        state->repterminationtype = 2;
        result = ae_false;
        return result;
    }
    
    /*
     * If Wolfe conditions are satisfied, we can update
     * limited memory model.
     *
     * However, if conditions are not satisfied (NFEV limit is met,
     * function is too wild, ...), we'll skip L-BFGS update
     */
    if( mcinfo!=1 )
    {
        
        /*
         * Skip update.
         *
         * In such cases we'll initialize search direction by
         * antigradient vector, because it  leads to more
         * transparent code with less number of special cases
         */
        state->fold = state->f;
        ae_v_moveneg(&state->d.ptr.p_double[0], 1, &state->g.ptr.p_double[0], 1, ae_v_len(0,n-1));
    }
    else
    {
        
        /*
         * Calculate Rho[k], GammaK
         */
        v = ae_v_dotproduct(&state->yk.ptr.pp_double[state->p][0], 1, &state->sk.ptr.pp_double[state->p][0], 1, ae_v_len(0,n-1));
        vv = ae_v_dotproduct(&state->yk.ptr.pp_double[state->p][0], 1, &state->yk.ptr.pp_double[state->p][0], 1, ae_v_len(0,n-1));
        if( ae_fp_eq(v,0)||ae_fp_eq(vv,0) )
        {
            
            /*
             * Rounding errors make further iterations impossible.
             */
            state->repterminationtype = -2;
            result = ae_false;
            return result;
        }
        state->rho.ptr.p_double[state->p] = 1/v;
        state->gammak = v/vv;
        
        /*
         *  Calculate d(k+1) = -H(k+1)*g(k+1)
         *
         *  for I:=K downto K-Q do
         *      V = s(i)^T * work(iteration:I)
         *      theta(i) = V
         *      work(iteration:I+1) = work(iteration:I) - V*Rho(i)*y(i)
         *  work(last iteration) = H0*work(last iteration) - preconditioner
         *  for I:=K-Q to K do
         *      V = y(i)^T*work(iteration:I)
         *      work(iteration:I+1) = work(iteration:I) +(-V+theta(i))*Rho(i)*s(i)
         *
         *  NOW WORK CONTAINS d(k+1)
         */
        ae_v_move(&state->work.ptr.p_double[0], 1, &state->g.ptr.p_double[0], 1, ae_v_len(0,n-1));
        for(i=state->k; i>=state->k-state->q; i--)
        {
            ic = i%m;
            v = ae_v_dotproduct(&state->sk.ptr.pp_double[ic][0], 1, &state->work.ptr.p_double[0], 1, ae_v_len(0,n-1));
            state->theta.ptr.p_double[ic] = v;
            vv = v*state->rho.ptr.p_double[ic];
            ae_v_subd(&state->work.ptr.p_double[0], 1, &state->yk.ptr.pp_double[ic][0], 1, ae_v_len(0,n-1), vv);
        }
        if( state->prectype==0 )
        {
            
            /*
             * Simple preconditioner is used
             */
            v = state->gammak;
            ae_v_muld(&state->work.ptr.p_double[0], 1, ae_v_len(0,n-1), v);
        }
        if( state->prectype==1 )
        {
            
            /*
             * Cholesky preconditioner is used
             */
            fblscholeskysolve(&state->denseh, 1, n, ae_true, &state->work, &state->autobuf, _state);
        }
        if( state->prectype==2 )
        {
            
            /*
             * diagonal approximation is used
             */
            for(i=0; i<=n-1; i++)
            {
                state->work.ptr.p_double[i] = state->work.ptr.p_double[i]/state->diagh.ptr.p_double[i];
            }
        }
        if( state->prectype==3 )
        {
            
            /*
             * scale-based preconditioner is used
             */
            for(i=0; i<=n-1; i++)
            {
                state->work.ptr.p_double[i] = state->work.ptr.p_double[i]*state->s.ptr.p_double[i]*state->s.ptr.p_double[i];
            }
        }
        for(i=state->k-state->q; i<=state->k; i++)
        {
            ic = i%m;
            v = ae_v_dotproduct(&state->yk.ptr.pp_double[ic][0], 1, &state->work.ptr.p_double[0], 1, ae_v_len(0,n-1));
            vv = state->rho.ptr.p_double[ic]*(-v+state->theta.ptr.p_double[ic]);
            ae_v_addd(&state->work.ptr.p_double[0], 1, &state->sk.ptr.pp_double[ic][0], 1, ae_v_len(0,n-1), vv);
        }
        ae_v_moveneg(&state->d.ptr.p_double[0], 1, &state->work.ptr.p_double[0], 1, ae_v_len(0,n-1));
        
        /*
         * Next step
         */
        state->fold = state->f;
        state->k = state->k+1;
    }
    goto lbl_29;
lbl_30:
    result = ae_false;
    return result;
    
    /*
     * Saving state
     */
lbl_rcomm:
    result = ae_true;
    state->rstate.ia.ptr.p_int[0] = n;
    state->rstate.ia.ptr.p_int[1] = m;
    state->rstate.ia.ptr.p_int[2] = i;
    state->rstate.ia.ptr.p_int[3] = j;
    state->rstate.ia.ptr.p_int[4] = ic;
    state->rstate.ia.ptr.p_int[5] = mcinfo;
    state->rstate.ra.ptr.p_double[0] = v;
    state->rstate.ra.ptr.p_double[1] = vv;
    return result;
}


/*************************************************************************
L-BFGS algorithm results

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    X       -   array[0..N-1], solution
    Rep     -   optimization report:
                * Rep.TerminationType completetion code:
                    * -7    gradient verification failed.
                            See MinLBFGSSetGradientCheck() for more information.
                    * -2    rounding errors prevent further improvement.
                            X contains best point found.
                    * -1    incorrect parameters were specified
                    *  1    relative function improvement is no more than
                            EpsF.
                    *  2    relative step is no more than EpsX.
                    *  4    gradient norm is no more than EpsG
                    *  5    MaxIts steps was taken
                    *  7    stopping conditions are too stringent,
                            further improvement is impossible
                * Rep.IterationsCount contains iterations count
                * NFEV countains number of function calculations

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgsresults(minlbfgsstate* state,
     /* Real    */ ae_vector* x,
     minlbfgsreport* rep,
     ae_state *_state)
{

    ae_vector_clear(x);
    _minlbfgsreport_clear(rep);

    minlbfgsresultsbuf(state, x, rep, _state);
}


/*************************************************************************
L-BFGS algorithm results

Buffered implementation of MinLBFGSResults which uses pre-allocated buffer
to store X[]. If buffer size is  too  small,  it  resizes  buffer.  It  is
intended to be used in the inner cycles of performance critical algorithms
where array reallocation penalty is too large to be ignored.

  -- ALGLIB --
     Copyright 20.08.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgsresultsbuf(minlbfgsstate* state,
     /* Real    */ ae_vector* x,
     minlbfgsreport* rep,
     ae_state *_state)
{


    if( x->cnt<state->n )
    {
        ae_vector_set_length(x, state->n, _state);
    }
    ae_v_move(&x->ptr.p_double[0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,state->n-1));
    rep->iterationscount = state->repiterationscount;
    rep->nfev = state->repnfev;
    rep->varidx = state->repvaridx;
    rep->terminationtype = state->repterminationtype;
}


/*************************************************************************
This  subroutine restarts LBFGS algorithm from new point. All optimization
parameters are left unchanged.

This  function  allows  to  solve multiple  optimization  problems  (which
must have same number of dimensions) without object reallocation penalty.

INPUT PARAMETERS:
    State   -   structure used to store algorithm state
    X       -   new starting point.

  -- ALGLIB --
     Copyright 30.07.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgsrestartfrom(minlbfgsstate* state,
     /* Real    */ ae_vector* x,
     ae_state *_state)
{


    ae_assert(x->cnt>=state->n, "MinLBFGSRestartFrom: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, state->n, _state), "MinLBFGSRestartFrom: X contains infinite or NaN values!", _state);
    ae_v_move(&state->x.ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,state->n-1));
    ae_vector_set_length(&state->rstate.ia, 5+1, _state);
    ae_vector_set_length(&state->rstate.ra, 1+1, _state);
    state->rstate.stage = -1;
    minlbfgs_clearrequestfields(state, _state);
}


/*************************************************************************
This  subroutine  turns  on  verification  of  the  user-supplied analytic
gradient:
* user calls this subroutine before optimization begins
* MinLBFGSOptimize() is called
* prior to  actual  optimization, for each component  of  parameters being
  optimized X[i] algorithm performs following steps:
  * two trial steps are made to X[i]-TestStep*S[i] and X[i]+TestStep*S[i],
    where X[i] is i-th component of the initial point and S[i] is a  scale
    of i-th parameter
  * if needed, steps are bounded with respect to constraints on X[]
  * F(X) is evaluated at these trial points
  * we perform one more evaluation in the middle point of the interval
  * we  build  cubic  model using function values and derivatives at trial
    points and we compare its prediction with actual value in  the  middle
    point
  * in case difference between prediction and actual value is higher  than
    some predetermined threshold, algorithm stops with completion code -7;
    Rep.VarIdx is set to index of the parameter with incorrect derivative.
* after verification is over, algorithm proceeds to the actual optimization.

NOTE 1: verification  needs  N (parameters count) gradient evaluations. It
        is very costly and you should use  it  only  for  low  dimensional
        problems,  when  you  want  to  be  sure  that  you've   correctly
        calculated  analytic  derivatives.  You  should  not use it in the
        production code (unless you want to check derivatives provided  by
        some third party).

NOTE 2: you  should  carefully  choose  TestStep. Value which is too large
        (so large that function behaviour is significantly non-cubic) will
        lead to false alarms. You may use  different  step  for  different
        parameters by means of setting scale with MinLBFGSSetScale().

NOTE 3: this function may lead to false positives. In case it reports that
        I-th  derivative was calculated incorrectly, you may decrease test
        step  and  try  one  more  time  - maybe your function changes too
        sharply  and  your  step  is  too  large for such rapidly chanding
        function.

INPUT PARAMETERS:
    State       -   structure used to store algorithm state
    TestStep    -   verification step:
                    * TestStep=0 turns verification off
                    * TestStep>0 activates verification

  -- ALGLIB --
     Copyright 24.05.2012 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetgradientcheck(minlbfgsstate* state,
     double teststep,
     ae_state *_state)
{


    ae_assert(ae_isfinite(teststep, _state), "MinLBFGSSetGradientCheck: TestStep contains NaN or Infinite", _state);
    ae_assert(ae_fp_greater_eq(teststep,0), "MinLBFGSSetGradientCheck: invalid argument TestStep(TestStep<0)", _state);
    state->teststep = teststep;
}


/*************************************************************************
Clears request fileds (to be sure that we don't forgot to clear something)
*************************************************************************/
static void minlbfgs_clearrequestfields(minlbfgsstate* state,
     ae_state *) //_state)
{


    state->needf = ae_false;
    state->needfg = ae_false;
    state->xupdated = ae_false;
}


ae_bool _minlbfgsstate_init(minlbfgsstate* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->s, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->rho, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->yk, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->sk, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->theta, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->d, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->work, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->denseh, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->diagh, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->autobuf, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->x, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->g, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !_rcommstate_init(&p->rstate, _state, make_automatic) )
        return ae_false;
    if( !_linminstate_init(&p->lstate, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _minlbfgsstate_init_copy(minlbfgsstate* dst, minlbfgsstate* src, ae_state *_state, ae_bool make_automatic)
{
    dst->n = src->n;
    dst->m = src->m;
    dst->epsg = src->epsg;
    dst->epsf = src->epsf;
    dst->epsx = src->epsx;
    dst->maxits = src->maxits;
    dst->xrep = src->xrep;
    dst->stpmax = src->stpmax;
    if( !ae_vector_init_copy(&dst->s, &src->s, _state, make_automatic) )
        return ae_false;
    dst->diffstep = src->diffstep;
    dst->nfev = src->nfev;
    dst->mcstage = src->mcstage;
    dst->k = src->k;
    dst->q = src->q;
    dst->p = src->p;
    if( !ae_vector_init_copy(&dst->rho, &src->rho, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->yk, &src->yk, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->sk, &src->sk, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->theta, &src->theta, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->d, &src->d, _state, make_automatic) )
        return ae_false;
    dst->stp = src->stp;
    if( !ae_vector_init_copy(&dst->work, &src->work, _state, make_automatic) )
        return ae_false;
    dst->fold = src->fold;
    dst->trimthreshold = src->trimthreshold;
    dst->prectype = src->prectype;
    dst->gammak = src->gammak;
    if( !ae_matrix_init_copy(&dst->denseh, &src->denseh, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->diagh, &src->diagh, _state, make_automatic) )
        return ae_false;
    dst->fbase = src->fbase;
    dst->fm2 = src->fm2;
    dst->fm1 = src->fm1;
    dst->fp1 = src->fp1;
    dst->fp2 = src->fp2;
    if( !ae_vector_init_copy(&dst->autobuf, &src->autobuf, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->x, &src->x, _state, make_automatic) )
        return ae_false;
    dst->f = src->f;
    if( !ae_vector_init_copy(&dst->g, &src->g, _state, make_automatic) )
        return ae_false;
    dst->needf = src->needf;
    dst->needfg = src->needfg;
    dst->xupdated = src->xupdated;
    dst->teststep = src->teststep;
    if( !_rcommstate_init_copy(&dst->rstate, &src->rstate, _state, make_automatic) )
        return ae_false;
    dst->repiterationscount = src->repiterationscount;
    dst->repnfev = src->repnfev;
    dst->repvaridx = src->repvaridx;
    dst->repterminationtype = src->repterminationtype;
    if( !_linminstate_init_copy(&dst->lstate, &src->lstate, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _minlbfgsstate_clear(minlbfgsstate* p)
{
    ae_vector_clear(&p->s);
    ae_vector_clear(&p->rho);
    ae_matrix_clear(&p->yk);
    ae_matrix_clear(&p->sk);
    ae_vector_clear(&p->theta);
    ae_vector_clear(&p->d);
    ae_vector_clear(&p->work);
    ae_matrix_clear(&p->denseh);
    ae_vector_clear(&p->diagh);
    ae_vector_clear(&p->autobuf);
    ae_vector_clear(&p->x);
    ae_vector_clear(&p->g);
    _rcommstate_clear(&p->rstate);
    _linminstate_clear(&p->lstate);
}


ae_bool _minlbfgsreport_init(minlbfgsreport* /*p*/, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    return ae_true;
}


ae_bool _minlbfgsreport_init_copy(minlbfgsreport* dst, minlbfgsreport* src, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    dst->iterationscount = src->iterationscount;
    dst->nfev = src->nfev;
    dst->varidx = src->varidx;
    dst->terminationtype = src->terminationtype;
    return ae_true;
}


void _minlbfgsreport_clear(minlbfgsreport* ) //p)
{
}




/*************************************************************************
This subroutine is used to initialize CQM. By default, empty NxN model  is
generated, with Alpha=Lambda=Theta=0.0 and zero b.

Previously allocated buffer variables are reused as much as possible.

  -- ALGLIB --
     Copyright 12.06.2012 by Bochkanov Sergey
*************************************************************************/
void cqminit(ae_int_t n, convexquadraticmodel* s, ae_state *_state)
{
    ae_int_t i;


    s->n = n;
    s->k = 0;
    s->nfree = n;
    s->ecakind = -1;
    s->alpha = 0.0;
    s->tau = 0.0;
    s->theta = 0.0;
    s->ismaintermchanged = ae_true;
    s->issecondarytermchanged = ae_true;
    s->islineartermchanged = ae_true;
    s->isactivesetchanged = ae_true;
    bvectorsetlengthatleast(&s->activeset, n, _state);
    rvectorsetlengthatleast(&s->xc, n, _state);
    rvectorsetlengthatleast(&s->eb, n, _state);
    rvectorsetlengthatleast(&s->tq1, n, _state);
    rvectorsetlengthatleast(&s->txc, n, _state);
    rvectorsetlengthatleast(&s->tb, n, _state);
    rvectorsetlengthatleast(&s->b, s->n, _state);
    rvectorsetlengthatleast(&s->tk1, s->n, _state);
    for(i=0; i<=n-1; i++)
    {
        s->activeset.ptr.p_bool[i] = ae_false;
        s->xc.ptr.p_double[i] = 0.0;
        s->b.ptr.p_double[i] = 0.0;
    }
}


/*************************************************************************
This subroutine changes main quadratic term of the model.

INPUT PARAMETERS:
    S       -   model
    A       -   NxN matrix, only upper or lower triangle is referenced
    IsUpper -   True, when matrix is stored in upper triangle
    Alpha   -   multiplier; when Alpha=0, A is not referenced at all

  -- ALGLIB --
     Copyright 12.06.2012 by Bochkanov Sergey
*************************************************************************/
void cqmseta(convexquadraticmodel* s,
     /* Real    */ ae_matrix* a,
     ae_bool isupper,
     double alpha,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    double v;


    ae_assert(ae_isfinite(alpha, _state)&&ae_fp_greater_eq(alpha,0), "CQMSetA: Alpha<0 or is not finite number", _state);
    ae_assert(ae_fp_eq(alpha,0)||isfinitertrmatrix(a, s->n, isupper, _state), "CQMSetA: A is not finite NxN matrix", _state);
    s->alpha = alpha;
    if( ae_fp_greater(alpha,0) )
    {
        rmatrixsetlengthatleast(&s->a, s->n, s->n, _state);
        rmatrixsetlengthatleast(&s->ecadense, s->n, s->n, _state);
        rmatrixsetlengthatleast(&s->tq2dense, s->n, s->n, _state);
        for(i=0; i<=s->n-1; i++)
        {
            for(j=i; j<=s->n-1; j++)
            {
                if( isupper )
                {
                    v = a->ptr.pp_double[i][j];
                }
                else
                {
                    v = a->ptr.pp_double[j][i];
                }
                s->a.ptr.pp_double[i][j] = v;
                s->a.ptr.pp_double[j][i] = v;
            }
        }
    }
    s->ismaintermchanged = ae_true;
}


/*************************************************************************
This subroutine rewrites diagonal of the main quadratic term of the  model
(dense  A)  by  vector  Z/Alpha (current value of the Alpha coefficient is
used).

IMPORTANT: in  case  model  has  no  dense  quadratic  term, this function
           allocates N*N dense matrix of zeros, and fills its diagonal  by
           non-zero values.

INPUT PARAMETERS:
    S       -   model
    Z       -   new diagonal, array[N]

  -- ALGLIB --
     Copyright 12.06.2012 by Bochkanov Sergey
*************************************************************************/
void cqmrewritedensediagonal(convexquadraticmodel* s,
     /* Real    */ ae_vector* z,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t i;
    ae_int_t j;


    n = s->n;
    if( ae_fp_eq(s->alpha,0) )
    {
        rmatrixsetlengthatleast(&s->a, s->n, s->n, _state);
        rmatrixsetlengthatleast(&s->ecadense, s->n, s->n, _state);
        rmatrixsetlengthatleast(&s->tq2dense, s->n, s->n, _state);
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                s->a.ptr.pp_double[i][j] = 0.0;
            }
        }
        s->alpha = 1.0;
    }
    for(i=0; i<=s->n-1; i++)
    {
        s->a.ptr.pp_double[i][i] = z->ptr.p_double[i]/s->alpha;
    }
    s->ismaintermchanged = ae_true;
}


/*************************************************************************
This subroutine changes diagonal quadratic term of the model.

INPUT PARAMETERS:
    S       -   model
    D       -   array[N], semidefinite diagonal matrix
    Tau     -   multiplier; when Tau=0, D is not referenced at all

  -- ALGLIB --
     Copyright 12.06.2012 by Bochkanov Sergey
*************************************************************************/
void cqmsetd(convexquadraticmodel* s,
     /* Real    */ ae_vector* d,
     double tau,
     ae_state *_state)
{
    ae_int_t i;


    ae_assert(ae_isfinite(tau, _state)&&ae_fp_greater_eq(tau,0), "CQMSetD: Tau<0 or is not finite number", _state);
    ae_assert(ae_fp_eq(tau,0)||isfinitevector(d, s->n, _state), "CQMSetD: D is not finite Nx1 vector", _state);
    s->tau = tau;
    if( ae_fp_greater(tau,0) )
    {
        rvectorsetlengthatleast(&s->d, s->n, _state);
        rvectorsetlengthatleast(&s->ecadiag, s->n, _state);
        rvectorsetlengthatleast(&s->tq2diag, s->n, _state);
        for(i=0; i<=s->n-1; i++)
        {
            ae_assert(ae_fp_greater_eq(d->ptr.p_double[i],0), "CQMSetD: D[i]<0", _state);
            s->d.ptr.p_double[i] = d->ptr.p_double[i];
        }
    }
    s->ismaintermchanged = ae_true;
}


/*************************************************************************
This subroutine drops main quadratic term A from the model. It is same  as
call  to  CQMSetA()  with  zero  A,   but gives better performance because
algorithm  knows  that  matrix  is  zero  and  can  optimize    subsequent
calculations.

INPUT PARAMETERS:
    S       -   model

  -- ALGLIB --
     Copyright 12.06.2012 by Bochkanov Sergey
*************************************************************************/
void cqmdropa(convexquadraticmodel* s, ae_state *) //_state)
{


    s->alpha = 0.0;
    s->ismaintermchanged = ae_true;
}


/*************************************************************************
This subroutine changes linear term of the model

  -- ALGLIB --
     Copyright 12.06.2012 by Bochkanov Sergey
*************************************************************************/
void cqmsetb(convexquadraticmodel* s,
     /* Real    */ ae_vector* b,
     ae_state *_state)
{
    ae_int_t i;


    ae_assert(isfinitevector(b, s->n, _state), "CQMSetB: B is not finite vector", _state);
    rvectorsetlengthatleast(&s->b, s->n, _state);
    for(i=0; i<=s->n-1; i++)
    {
        s->b.ptr.p_double[i] = b->ptr.p_double[i];
    }
    s->islineartermchanged = ae_true;
}


/*************************************************************************
This subroutine changes linear term of the model

  -- ALGLIB --
     Copyright 12.06.2012 by Bochkanov Sergey
*************************************************************************/
void cqmsetq(convexquadraticmodel* s,
     /* Real    */ ae_matrix* q,
     /* Real    */ ae_vector* r,
     ae_int_t k,
     double theta,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;


    ae_assert(k>=0, "CQMSetQ: K<0", _state);
    ae_assert((k==0||ae_fp_eq(theta,0))||apservisfinitematrix(q, k, s->n, _state), "CQMSetQ: Q is not finite matrix", _state);
    ae_assert((k==0||ae_fp_eq(theta,0))||isfinitevector(r, k, _state), "CQMSetQ: R is not finite vector", _state);
    ae_assert(ae_isfinite(theta, _state)&&ae_fp_greater_eq(theta,0), "CQMSetQ: Theta<0 or is not finite number", _state);
    
    /*
     * degenerate case: K=0 or Theta=0
     */
    if( k==0||ae_fp_eq(theta,0) )
    {
        s->k = 0;
        s->theta = 0;
        s->issecondarytermchanged = ae_true;
        return;
    }
    
    /*
     * General case: both Theta>0 and K>0
     */
    s->k = k;
    s->theta = theta;
    rmatrixsetlengthatleast(&s->q, s->k, s->n, _state);
    rvectorsetlengthatleast(&s->r, s->k, _state);
    rmatrixsetlengthatleast(&s->eq, s->k, s->n, _state);
    rmatrixsetlengthatleast(&s->eccm, s->k, s->k, _state);
    rmatrixsetlengthatleast(&s->tk2, s->k, s->n, _state);
    for(i=0; i<=s->k-1; i++)
    {
        for(j=0; j<=s->n-1; j++)
        {
            s->q.ptr.pp_double[i][j] = q->ptr.pp_double[i][j];
        }
        s->r.ptr.p_double[i] = r->ptr.p_double[i];
    }
    s->issecondarytermchanged = ae_true;
}


/*************************************************************************
This subroutine changes active set

INPUT PARAMETERS
    S       -   model
    X       -   array[N], constraint values
    ActiveSet-  array[N], active set. If ActiveSet[I]=True, then I-th
                variables is constrained to X[I].

  -- ALGLIB --
     Copyright 12.06.2012 by Bochkanov Sergey
*************************************************************************/
void cqmsetactiveset(convexquadraticmodel* s,
     /* Real    */ ae_vector* x,
     /* Boolean */ ae_vector* activeset,
     ae_state *_state)
{
    ae_int_t i;


    ae_assert(x->cnt>=s->n, "CQMSetActiveSet: Length(X)<N", _state);
    ae_assert(activeset->cnt>=s->n, "CQMSetActiveSet: Length(ActiveSet)<N", _state);
    for(i=0; i<=s->n-1; i++)
    {
        s->isactivesetchanged = s->isactivesetchanged||(s->activeset.ptr.p_bool[i]&&!activeset->ptr.p_bool[i]);
        s->isactivesetchanged = s->isactivesetchanged||(activeset->ptr.p_bool[i]&&!s->activeset.ptr.p_bool[i]);
        s->activeset.ptr.p_bool[i] = activeset->ptr.p_bool[i];
        if( activeset->ptr.p_bool[i] )
        {
            ae_assert(ae_isfinite(x->ptr.p_double[i], _state), "CQMSetActiveSet: X[] contains infinite constraints", _state);
            s->isactivesetchanged = s->isactivesetchanged||ae_fp_neq(s->xc.ptr.p_double[i],x->ptr.p_double[i]);
            s->xc.ptr.p_double[i] = x->ptr.p_double[i];
        }
    }
}


/*************************************************************************
This subroutine evaluates model at X. Active constraints are ignored.

  -- ALGLIB --
     Copyright 12.06.2012 by Bochkanov Sergey
*************************************************************************/
double cqmeval(convexquadraticmodel* s,
     /* Real    */ ae_vector* x,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t i;
    ae_int_t j;
    double v;
    double result;


    n = s->n;
    ae_assert(isfinitevector(x, n, _state), "CQMEval: X is not finite vector", _state);
    result = 0.0;
    
    /*
     * main quadratic term
     */
    if( ae_fp_greater(s->alpha,0) )
    {
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                result = result+s->alpha*0.5*x->ptr.p_double[i]*s->a.ptr.pp_double[i][j]*x->ptr.p_double[j];
            }
        }
    }
    if( ae_fp_greater(s->tau,0) )
    {
        for(i=0; i<=n-1; i++)
        {
            result = result+0.5*ae_sqr(x->ptr.p_double[i], _state)*s->tau*s->d.ptr.p_double[i];
        }
    }
    
    /*
     * secondary quadratic term
     */
    if( ae_fp_greater(s->theta,0) )
    {
        for(i=0; i<=s->k-1; i++)
        {
            v = ae_v_dotproduct(&s->q.ptr.pp_double[i][0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,n-1));
            result = result+0.5*s->theta*ae_sqr(v-s->r.ptr.p_double[i], _state);
        }
    }
    
    /*
     * linear term
     */
    for(i=0; i<=s->n-1; i++)
    {
        result = result+x->ptr.p_double[i]*s->b.ptr.p_double[i];
    }
    return result;
}


/*************************************************************************
This subroutine evaluates model at X. Active constraints are ignored.
It returns:
    R   -   model value
    Noise-  estimate of the numerical noise in data

  -- ALGLIB --
     Copyright 12.06.2012 by Bochkanov Sergey
*************************************************************************/
void cqmevalx(convexquadraticmodel* s,
     /* Real    */ ae_vector* x,
     double* r,
     double* noise,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t i;
    ae_int_t j;
    double v;
    double v2;
    double mxq;
    double eps;

    *r = 0;
    *noise = 0;

    n = s->n;
    ae_assert(isfinitevector(x, n, _state), "CQMEval: X is not finite vector", _state);
    *r = 0.0;
    *noise = 0.0;
    eps = 2*ae_machineepsilon;
    mxq = 0.0;
    
    /*
     * Main quadratic term.
     *
     * Noise from the main quadratic term is equal to the
     * maximum summand in the term.
     */
    if( ae_fp_greater(s->alpha,0) )
    {
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                v = s->alpha*0.5*x->ptr.p_double[i]*s->a.ptr.pp_double[i][j]*x->ptr.p_double[j];
                *r = *r+v;
                *noise = ae_maxreal(*noise, eps*ae_fabs(v, _state), _state);
            }
        }
    }
    if( ae_fp_greater(s->tau,0) )
    {
        for(i=0; i<=n-1; i++)
        {
            v = 0.5*ae_sqr(x->ptr.p_double[i], _state)*s->tau*s->d.ptr.p_double[i];
            *r = *r+v;
            *noise = ae_maxreal(*noise, eps*ae_fabs(v, _state), _state);
        }
    }
    
    /*
     * secondary quadratic term
     *
     * Noise from the secondary quadratic term is estimated as follows:
     * * noise in qi*x-r[i] is estimated as
     *   Eps*MXQ = Eps*max(|r[i]|, |q[i,j]*x[j]|)
     * * noise in (qi*x-r[i])^2 is estimated as
     *   NOISE = (|qi*x-r[i]|+Eps*MXQ)^2-(|qi*x-r[i]|)^2
     *         = Eps*MXQ*(2*|qi*x-r[i]|+Eps*MXQ)
     */
    if( ae_fp_greater(s->theta,0) )
    {
        for(i=0; i<=s->k-1; i++)
        {
            v = 0.0;
            mxq = ae_fabs(s->r.ptr.p_double[i], _state);
            for(j=0; j<=n-1; j++)
            {
                v2 = s->q.ptr.pp_double[i][j]*x->ptr.p_double[j];
                v = v+v2;
                mxq = ae_maxreal(mxq, ae_fabs(v2, _state), _state);
            }
            *r = *r+0.5*s->theta*ae_sqr(v-s->r.ptr.p_double[i], _state);
            *noise = ae_maxreal(*noise, eps*mxq*(2*ae_fabs(v-s->r.ptr.p_double[i], _state)+eps*mxq), _state);
        }
    }
    
    /*
     * linear term
     */
    for(i=0; i<=s->n-1; i++)
    {
        *r = *r+x->ptr.p_double[i]*s->b.ptr.p_double[i];
        *noise = ae_maxreal(*noise, eps*ae_fabs(x->ptr.p_double[i]*s->b.ptr.p_double[i], _state), _state);
    }
    
    /*
     * Final update of the noise
     */
    *noise = n*(*noise);
}


/*************************************************************************
This  subroutine  evaluates  gradient of the model; active constraints are
ignored.

INPUT PARAMETERS:
    S       -   convex model
    X       -   point, array[N]
    G       -   possibly preallocated buffer; resized, if too small

  -- ALGLIB --
     Copyright 12.06.2012 by Bochkanov Sergey
*************************************************************************/
void cqmgradunconstrained(convexquadraticmodel* s,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* g,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t i;
    ae_int_t j;
    double v;


    n = s->n;
    ae_assert(isfinitevector(x, n, _state), "CQMEvalGradUnconstrained: X is not finite vector", _state);
    rvectorsetlengthatleast(g, n, _state);
    for(i=0; i<=n-1; i++)
    {
        g->ptr.p_double[i] = 0;
    }
    
    /*
     * main quadratic term
     */
    if( ae_fp_greater(s->alpha,0) )
    {
        for(i=0; i<=n-1; i++)
        {
            v = 0.0;
            for(j=0; j<=n-1; j++)
            {
                v = v+s->alpha*s->a.ptr.pp_double[i][j]*x->ptr.p_double[j];
            }
            g->ptr.p_double[i] = g->ptr.p_double[i]+v;
        }
    }
    if( ae_fp_greater(s->tau,0) )
    {
        for(i=0; i<=n-1; i++)
        {
            g->ptr.p_double[i] = g->ptr.p_double[i]+x->ptr.p_double[i]*s->tau*s->d.ptr.p_double[i];
        }
    }
    
    /*
     * secondary quadratic term
     */
    if( ae_fp_greater(s->theta,0) )
    {
        for(i=0; i<=s->k-1; i++)
        {
            v = ae_v_dotproduct(&s->q.ptr.pp_double[i][0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,n-1));
            v = s->theta*(v-s->r.ptr.p_double[i]);
            ae_v_addd(&g->ptr.p_double[0], 1, &s->q.ptr.pp_double[i][0], 1, ae_v_len(0,n-1), v);
        }
    }
    
    /*
     * linear term
     */
    for(i=0; i<=n-1; i++)
    {
        g->ptr.p_double[i] = g->ptr.p_double[i]+s->b.ptr.p_double[i];
    }
}


/*************************************************************************
This subroutine evaluates x'*(0.5*alpha*A+tau*D)*x

  -- ALGLIB --
     Copyright 12.06.2012 by Bochkanov Sergey
*************************************************************************/
double cqmxtadx2(convexquadraticmodel* s,
     /* Real    */ ae_vector* x,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t i;
    ae_int_t j;
    double result;


    n = s->n;
    ae_assert(isfinitevector(x, n, _state), "CQMEval: X is not finite vector", _state);
    result = 0.0;
    
    /*
     * main quadratic term
     */
    if( ae_fp_greater(s->alpha,0) )
    {
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                result = result+s->alpha*0.5*x->ptr.p_double[i]*s->a.ptr.pp_double[i][j]*x->ptr.p_double[j];
            }
        }
    }
    if( ae_fp_greater(s->tau,0) )
    {
        for(i=0; i<=n-1; i++)
        {
            result = result+0.5*ae_sqr(x->ptr.p_double[i], _state)*s->tau*s->d.ptr.p_double[i];
        }
    }
    return result;
}


/*************************************************************************
This subroutine evaluates (0.5*alpha*A+tau*D)*x

Y is automatically resized if needed

  -- ALGLIB --
     Copyright 12.06.2012 by Bochkanov Sergey
*************************************************************************/
void cqmadx(convexquadraticmodel* s,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t i;
    double v;


    n = s->n;
    ae_assert(isfinitevector(x, n, _state), "CQMEval: X is not finite vector", _state);
    rvectorsetlengthatleast(y, n, _state);
    
    /*
     * main quadratic term
     */
    for(i=0; i<=n-1; i++)
    {
        y->ptr.p_double[i] = 0;
    }
    if( ae_fp_greater(s->alpha,0) )
    {
        for(i=0; i<=n-1; i++)
        {
            v = ae_v_dotproduct(&s->a.ptr.pp_double[i][0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,n-1));
            y->ptr.p_double[i] = y->ptr.p_double[i]+s->alpha*v;
        }
    }
    if( ae_fp_greater(s->tau,0) )
    {
        for(i=0; i<=n-1; i++)
        {
            y->ptr.p_double[i] = y->ptr.p_double[i]+x->ptr.p_double[i]*s->tau*s->d.ptr.p_double[i];
        }
    }
}


/*************************************************************************
This subroutine finds optimum of the model. It returns  False  on  failure
(indefinite/semidefinite matrix).  Optimum  is  found  subject  to  active
constraints.

INPUT PARAMETERS
    S       -   model
    X       -   possibly preallocated buffer; automatically resized, if
                too small enough.

  -- ALGLIB --
     Copyright 12.06.2012 by Bochkanov Sergey
*************************************************************************/
ae_bool cqmconstrainedoptimum(convexquadraticmodel* s,
     /* Real    */ ae_vector* x,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t nfree;
    ae_int_t k;
    ae_int_t i;
    double v;
    ae_int_t cidx0;
    ae_int_t itidx;
    ae_bool result;


    
    /*
     * Rebuild internal structures
     */
    if( !cqmodels_cqmrebuild(s, _state) )
    {
        result = ae_false;
        return result;
    }
    n = s->n;
    k = s->k;
    nfree = s->nfree;
    result = ae_true;
    
    /*
     * Calculate initial point for the iterative refinement:
     * * free components are set to zero
     * * constrained components are set to their constrained values
     */
    rvectorsetlengthatleast(x, n, _state);
    for(i=0; i<=n-1; i++)
    {
        if( s->activeset.ptr.p_bool[i] )
        {
            x->ptr.p_double[i] = s->xc.ptr.p_double[i];
        }
        else
        {
            x->ptr.p_double[i] = 0;
        }
    }
    
    /*
     * Iterative refinement.
     *
     * In an ideal world without numerical errors it would be enough
     * to make just one Newton step from initial point:
     *   x_new = -H^(-1)*grad(x=0)
     * However, roundoff errors can significantly deteriorate quality
     * of the solution. So we have to recalculate gradient and to
     * perform Newton steps several times.
     *
     * Below we perform fixed number of Newton iterations.
     */
    for(itidx=0; itidx<=cqmodels_newtonrefinementits-1; itidx++)
    {
        
        /*
         * Calculate gradient at the current point.
         * Move free components of the gradient in the beginning.
         */
        cqmgradunconstrained(s, x, &s->tmpg, _state);
        cidx0 = 0;
        for(i=0; i<=n-1; i++)
        {
            if( !s->activeset.ptr.p_bool[i] )
            {
                s->tmpg.ptr.p_double[cidx0] = s->tmpg.ptr.p_double[i];
                cidx0 = cidx0+1;
            }
        }
        
        /*
         * Free components of the extrema are calculated in the first NFree elements of TXC.
         *
         * First, we have to calculate original Newton step, without rank-K perturbations
         */
        ae_v_moveneg(&s->txc.ptr.p_double[0], 1, &s->tmpg.ptr.p_double[0], 1, ae_v_len(0,nfree-1));
        cqmodels_cqmsolveea(s, &s->txc, &s->tmp0, _state);
        
        /*
         * Then, we account for rank-K correction.
         * Woodbury matrix identity is used.
         */
        if( s->k>0&&ae_fp_greater(s->theta,0) )
        {
            rvectorsetlengthatleast(&s->tmp0, ae_maxint(nfree, k, _state), _state);
            rvectorsetlengthatleast(&s->tmp1, ae_maxint(nfree, k, _state), _state);
            ae_v_moveneg(&s->tmp1.ptr.p_double[0], 1, &s->tmpg.ptr.p_double[0], 1, ae_v_len(0,nfree-1));
            cqmodels_cqmsolveea(s, &s->tmp1, &s->tmp0, _state);
            for(i=0; i<=k-1; i++)
            {
                v = ae_v_dotproduct(&s->eq.ptr.pp_double[i][0], 1, &s->tmp1.ptr.p_double[0], 1, ae_v_len(0,nfree-1));
                s->tmp0.ptr.p_double[i] = v;
            }
            fblscholeskysolve(&s->eccm, 1.0, k, ae_true, &s->tmp0, &s->tmp1, _state);
            for(i=0; i<=nfree-1; i++)
            {
                s->tmp1.ptr.p_double[i] = 0.0;
            }
            for(i=0; i<=k-1; i++)
            {
                v = s->tmp0.ptr.p_double[i];
                ae_v_addd(&s->tmp1.ptr.p_double[0], 1, &s->eq.ptr.pp_double[i][0], 1, ae_v_len(0,nfree-1), v);
            }
            cqmodels_cqmsolveea(s, &s->tmp1, &s->tmp0, _state);
            ae_v_sub(&s->txc.ptr.p_double[0], 1, &s->tmp1.ptr.p_double[0], 1, ae_v_len(0,nfree-1));
        }
        
        /*
         * Unpack components from TXC into X. We pass through all
         * free components of X and add our step.
         */
        cidx0 = 0;
        for(i=0; i<=n-1; i++)
        {
            if( !s->activeset.ptr.p_bool[i] )
            {
                x->ptr.p_double[i] = x->ptr.p_double[i]+s->txc.ptr.p_double[cidx0];
                cidx0 = cidx0+1;
            }
        }
    }
    return result;
}


/*************************************************************************
This function scales vector  by  multiplying it by inverse of the diagonal
of the Hessian matrix. It should be used to  accelerate  steepest  descent
phase of the QP solver.

Although  it  is  called  "scale-grad",  it  can be called for any vector,
whether it is gradient, anti-gradient, or just some vector.

This function does NOT takes into account current set of  constraints,  it
just performs matrix-vector multiplication  without  taking  into  account
constraints.

INPUT PARAMETERS:
    S       -   model
    X       -   vector to scale

OUTPUT PARAMETERS:
    X       -   scaled vector
    
NOTE:
    when called for non-SPD matrices, it silently skips components of X
    which correspond to zero or negative diagonal elements.
    
NOTE:
    this function uses diagonals of A and D; it ignores Q - rank-K term of
    the quadratic model.

  -- ALGLIB --
     Copyright 12.06.2012 by Bochkanov Sergey
*************************************************************************/
void cqmscalevector(convexquadraticmodel* s,
     /* Real    */ ae_vector* x,
     ae_state *) //_state)
{
    ae_int_t n;
    ae_int_t i;
    double v;


    n = s->n;
    for(i=0; i<=n-1; i++)
    {
        v = 0.0;
        if( ae_fp_greater(s->alpha,0) )
        {
            v = v+s->a.ptr.pp_double[i][i];
        }
        if( ae_fp_greater(s->tau,0) )
        {
            v = v+s->d.ptr.p_double[i];
        }
        if( ae_fp_greater(v,0) )
        {
            x->ptr.p_double[i] = x->ptr.p_double[i]/v;
        }
    }
}


/*************************************************************************
This subroutine calls CQMRebuild() and evaluates model at X subject to
active constraints.

It  is  intended  for  debug  purposes only, because it evaluates model by
means of temporaries, which were calculated  by  CQMRebuild().  The   only
purpose of this function  is  to  check  correctness  of  CQMRebuild()  by
comparing results of this function with ones obtained by CQMEval(),  which
is  used  as  reference  point. The  idea is that significant deviation in
results  of  these  two  functions  is  evidence  of  some  error  in  the
CQMRebuild().

NOTE: suffix T denotes that temporaries marked by T-prefix are used. There
      is one more variant of this function, which uses  "effective"  model
      built by CQMRebuild().

NOTE2: in case CQMRebuild() fails (due to model non-convexity), this
      function returns NAN.

  -- ALGLIB --
     Copyright 12.06.2012 by Bochkanov Sergey
*************************************************************************/
double cqmdebugconstrainedevalt(convexquadraticmodel* s,
     /* Real    */ ae_vector* x,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t nfree;
    ae_int_t i;
    ae_int_t j;
    double v;
    double result;


    n = s->n;
    ae_assert(isfinitevector(x, n, _state), "CQMDebugConstrainedEvalT: X is not finite vector", _state);
    if( !cqmodels_cqmrebuild(s, _state) )
    {
        result = _state->v_nan;
        return result;
    }
    result = 0.0;
    nfree = s->nfree;
    
    /*
     * Reorder variables
     */
    j = 0;
    for(i=0; i<=n-1; i++)
    {
        if( !s->activeset.ptr.p_bool[i] )
        {
            ae_assert(j<nfree, "CQMDebugConstrainedEvalT: internal error", _state);
            s->txc.ptr.p_double[j] = x->ptr.p_double[i];
            j = j+1;
        }
    }
    
    /*
     * TQ2, TQ1, TQ0
     *
     */
    if( ae_fp_greater(s->alpha,0) )
    {
        
        /*
         * Dense TQ2
         */
        for(i=0; i<=nfree-1; i++)
        {
            for(j=0; j<=nfree-1; j++)
            {
                result = result+0.5*s->txc.ptr.p_double[i]*s->tq2dense.ptr.pp_double[i][j]*s->txc.ptr.p_double[j];
            }
        }
    }
    else
    {
        
        /*
         * Diagonal TQ2
         */
        for(i=0; i<=nfree-1; i++)
        {
            result = result+0.5*s->tq2diag.ptr.p_double[i]*ae_sqr(s->txc.ptr.p_double[i], _state);
        }
    }
    for(i=0; i<=nfree-1; i++)
    {
        result = result+s->tq1.ptr.p_double[i]*s->txc.ptr.p_double[i];
    }
    result = result+s->tq0;
    
    /*
     * TK2, TK1, TK0
     */
    if( s->k>0&&ae_fp_greater(s->theta,0) )
    {
        for(i=0; i<=s->k-1; i++)
        {
            v = 0;
            for(j=0; j<=nfree-1; j++)
            {
                v = v+s->tk2.ptr.pp_double[i][j]*s->txc.ptr.p_double[j];
            }
            result = result+0.5*ae_sqr(v, _state);
        }
        for(i=0; i<=nfree-1; i++)
        {
            result = result+s->tk1.ptr.p_double[i]*s->txc.ptr.p_double[i];
        }
        result = result+s->tk0;
    }
    
    /*
     * TB (Bf and Bc parts)
     */
    for(i=0; i<=n-1; i++)
    {
        result = result+s->tb.ptr.p_double[i]*s->txc.ptr.p_double[i];
    }
    return result;
}


/*************************************************************************
This subroutine calls CQMRebuild() and evaluates model at X subject to
active constraints.

It  is  intended  for  debug  purposes only, because it evaluates model by
means of "effective" matrices built by CQMRebuild(). The only  purpose  of
this function is to check correctness of CQMRebuild() by comparing results
of this function with  ones  obtained  by  CQMEval(),  which  is  used  as
reference  point.  The  idea  is  that significant deviation in results of
these two functions is evidence of some error in the CQMRebuild().

NOTE: suffix E denotes that effective matrices. There is one more  variant
      of this function, which uses temporary matrices built by
      CQMRebuild().

NOTE2: in case CQMRebuild() fails (due to model non-convexity), this
      function returns NAN.

  -- ALGLIB --
     Copyright 12.06.2012 by Bochkanov Sergey
*************************************************************************/
double cqmdebugconstrainedevale(convexquadraticmodel* s,
     /* Real    */ ae_vector* x,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t nfree;
    ae_int_t i;
    ae_int_t j;
    double v;
    double result;


    n = s->n;
    ae_assert(isfinitevector(x, n, _state), "CQMDebugConstrainedEvalE: X is not finite vector", _state);
    if( !cqmodels_cqmrebuild(s, _state) )
    {
        result = _state->v_nan;
        return result;
    }
    result = 0.0;
    nfree = s->nfree;
    
    /*
     * Reorder variables
     */
    j = 0;
    for(i=0; i<=n-1; i++)
    {
        if( !s->activeset.ptr.p_bool[i] )
        {
            ae_assert(j<nfree, "CQMDebugConstrainedEvalE: internal error", _state);
            s->txc.ptr.p_double[j] = x->ptr.p_double[i];
            j = j+1;
        }
    }
    
    /*
     * ECA
     */
    ae_assert((s->ecakind==0||s->ecakind==1)||(s->ecakind==-1&&nfree==0), "CQMDebugConstrainedEvalE: unexpected ECAKind", _state);
    if( s->ecakind==0 )
    {
        
        /*
         * Dense ECA
         */
        for(i=0; i<=nfree-1; i++)
        {
            v = 0.0;
            for(j=i; j<=nfree-1; j++)
            {
                v = v+s->ecadense.ptr.pp_double[i][j]*s->txc.ptr.p_double[j];
            }
            result = result+0.5*ae_sqr(v, _state);
        }
    }
    if( s->ecakind==1 )
    {
        
        /*
         * Diagonal ECA
         */
        for(i=0; i<=nfree-1; i++)
        {
            result = result+0.5*ae_sqr(s->ecadiag.ptr.p_double[i]*s->txc.ptr.p_double[i], _state);
        }
    }
    
    /*
     * EQ
     */
    for(i=0; i<=s->k-1; i++)
    {
        v = 0.0;
        for(j=0; j<=nfree-1; j++)
        {
            v = v+s->eq.ptr.pp_double[i][j]*s->txc.ptr.p_double[j];
        }
        result = result+0.5*ae_sqr(v, _state);
    }
    
    /*
     * EB
     */
    for(i=0; i<=nfree-1; i++)
    {
        result = result+s->eb.ptr.p_double[i]*s->txc.ptr.p_double[i];
    }
    
    /*
     * EC
     */
    result = result+s->ec;
    return result;
}


/*************************************************************************
Internal function, rebuilds "effective" model subject to constraints.
Returns False on failure (non-SPD main quadratic term)

  -- ALGLIB --
     Copyright 10.05.2011 by Bochkanov Sergey
*************************************************************************/
static ae_bool cqmodels_cqmrebuild(convexquadraticmodel* s,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t nfree;
    ae_int_t k;
    ae_int_t i;
    ae_int_t j;
    ae_int_t ridx0;
    ae_int_t ridx1;
    ae_int_t cidx0;
    ae_int_t cidx1;
    double v;
    ae_bool result;


    if( ae_fp_eq(s->alpha,0)&&ae_fp_eq(s->tau,0) )
    {
        
        /*
         * Non-SPD model, quick exit
         */
        result = ae_false;
        return result;
    }
    result = ae_true;
    n = s->n;
    k = s->k;
    
    /*
     * Determine number of free variables.
     * Fill TXC - array whose last N-NFree elements store constraints.
     */
    if( s->isactivesetchanged )
    {
        s->nfree = 0;
        for(i=0; i<=n-1; i++)
        {
            if( !s->activeset.ptr.p_bool[i] )
            {
                s->nfree = s->nfree+1;
            }
        }
        j = s->nfree;
        for(i=0; i<=n-1; i++)
        {
            if( s->activeset.ptr.p_bool[i] )
            {
                s->txc.ptr.p_double[j] = s->xc.ptr.p_double[i];
                j = j+1;
            }
        }
    }
    nfree = s->nfree;
    
    /*
     * Re-evaluate TQ2/TQ1/TQ0, if needed
     */
    if( s->isactivesetchanged||s->ismaintermchanged )
    {
        
        /*
         * Handle cases Alpha>0 and Alpha=0 separately:
         * * in the first case we have dense matrix
         * * in the second one we have diagonal matrix, which can be
         *   handled more efficiently
         */
        if( ae_fp_greater(s->alpha,0) )
        {
            
            /*
             * Alpha>0, dense QP
             *
             * Split variables into two groups - free (F) and constrained (C). Reorder
             * variables in such way that free vars come first, constrained are last:
             * x = [xf, xc].
             * 
             * Main quadratic term x'*(alpha*A+tau*D)*x now splits into quadratic part,
             * linear part and constant part:
             *                   ( alpha*Aff+tau*Df  alpha*Afc        ) ( xf )              
             *   0.5*( xf' xc' )*(                                    )*(    ) =
             *                   ( alpha*Acf         alpha*Acc+tau*Dc ) ( xc )
             *
             *   = 0.5*xf'*(alpha*Aff+tau*Df)*xf + (alpha*Afc*xc)'*xf + 0.5*xc'(alpha*Acc+tau*Dc)*xc
             *                    
             * We store these parts into temporary variables:
             * * alpha*Aff+tau*Df, alpha*Afc, alpha*Acc+tau*Dc are stored into upper
             *   triangle of TQ2
             * * alpha*Afc*xc is stored into TQ1
             * * 0.5*xc'(alpha*Acc+tau*Dc)*xc is stored into TQ0
             *
             * Below comes first part of the work - generation of TQ2:
             * * we pass through rows of A and copy I-th row into upper block (Aff/Afc) or
             *   lower one (Acf/Acc) of TQ2, depending on presence of X[i] in the active set.
             *   RIdx0 variable contains current position for insertion into upper block,
             *   RIdx1 contains current position for insertion into lower one.
             * * within each row, we copy J-th element into left half (Aff/Acf) or right
             *   one (Afc/Acc), depending on presence of X[j] in the active set. CIdx0
             *   contains current position for insertion into left block, CIdx1 contains
             *   position for insertion into right one.
             * * during copying, we multiply elements by alpha and add diagonal matrix D.
             */
            ridx0 = 0;
            ridx1 = s->nfree;
            for(i=0; i<=n-1; i++)
            {
                cidx0 = 0;
                cidx1 = s->nfree;
                for(j=0; j<=n-1; j++)
                {
                    if( !s->activeset.ptr.p_bool[i]&&!s->activeset.ptr.p_bool[j] )
                    {
                        
                        /*
                         * Element belongs to Aff
                         */
                        v = s->alpha*s->a.ptr.pp_double[i][j];
                        if( i==j&&ae_fp_greater(s->tau,0) )
                        {
                            v = v+s->tau*s->d.ptr.p_double[i];
                        }
                        s->tq2dense.ptr.pp_double[ridx0][cidx0] = v;
                    }
                    if( !s->activeset.ptr.p_bool[i]&&s->activeset.ptr.p_bool[j] )
                    {
                        
                        /*
                         * Element belongs to Afc
                         */
                        s->tq2dense.ptr.pp_double[ridx0][cidx1] = s->alpha*s->a.ptr.pp_double[i][j];
                    }
                    if( s->activeset.ptr.p_bool[i]&&!s->activeset.ptr.p_bool[j] )
                    {
                        
                        /*
                         * Element belongs to Acf
                         */
                        s->tq2dense.ptr.pp_double[ridx1][cidx0] = s->alpha*s->a.ptr.pp_double[i][j];
                    }
                    if( s->activeset.ptr.p_bool[i]&&s->activeset.ptr.p_bool[j] )
                    {
                        
                        /*
                         * Element belongs to Acc
                         */
                        v = s->alpha*s->a.ptr.pp_double[i][j];
                        if( i==j&&ae_fp_greater(s->tau,0) )
                        {
                            v = v+s->tau*s->d.ptr.p_double[i];
                        }
                        s->tq2dense.ptr.pp_double[ridx1][cidx1] = v;
                    }
                    if( s->activeset.ptr.p_bool[j] )
                    {
                        cidx1 = cidx1+1;
                    }
                    else
                    {
                        cidx0 = cidx0+1;
                    }
                }
                if( s->activeset.ptr.p_bool[i] )
                {
                    ridx1 = ridx1+1;
                }
                else
                {
                    ridx0 = ridx0+1;
                }
            }
            
            /*
             * Now we have TQ2, and we can evaluate TQ1.
             * In the special case when we have Alpha=0, NFree=0 or NFree=N,
             * TQ1 is filled by zeros.
             */
            for(i=0; i<=n-1; i++)
            {
                s->tq1.ptr.p_double[i] = 0.0;
            }
            if( s->nfree>0&&s->nfree<n )
            {
                rmatrixmv(s->nfree, n-s->nfree, &s->tq2dense, 0, s->nfree, 0, &s->txc, s->nfree, &s->tq1, 0, _state);
            }
            
            /*
             * And finally, we evaluate TQ0.
             */
            v = 0.0;
            for(i=s->nfree; i<=n-1; i++)
            {
                for(j=s->nfree; j<=n-1; j++)
                {
                    v = v+0.5*s->txc.ptr.p_double[i]*s->tq2dense.ptr.pp_double[i][j]*s->txc.ptr.p_double[j];
                }
            }
            s->tq0 = v;
        }
        else
        {
            
            /*
             * Alpha=0, diagonal QP
             *
             * Split variables into two groups - free (F) and constrained (C). Reorder
             * variables in such way that free vars come first, constrained are last:
             * x = [xf, xc].
             * 
             * Main quadratic term x'*(tau*D)*x now splits into quadratic and constant
             * parts:
             *                   ( tau*Df        ) ( xf )              
             *   0.5*( xf' xc' )*(               )*(    ) =
             *                   (        tau*Dc ) ( xc )
             *
             *   = 0.5*xf'*(tau*Df)*xf + 0.5*xc'(tau*Dc)*xc
             *                    
             * We store these parts into temporary variables:
             * * tau*Df is stored in TQ2Diag
             * * 0.5*xc'(tau*Dc)*xc is stored into TQ0
             */
            s->tq0 = 0.0;
            ridx0 = 0;
            for(i=0; i<=n-1; i++)
            {
                if( !s->activeset.ptr.p_bool[i] )
                {
                    s->tq2diag.ptr.p_double[ridx0] = s->tau*s->d.ptr.p_double[i];
                    ridx0 = ridx0+1;
                }
                else
                {
                    s->tq0 = s->tq0+0.5*s->tau*s->d.ptr.p_double[i]*ae_sqr(s->xc.ptr.p_double[i], _state);
                }
            }
            for(i=0; i<=n-1; i++)
            {
                s->tq1.ptr.p_double[i] = 0.0;
            }
        }
    }
    
    /*
     * Re-evaluate TK2/TK1/TK0, if needed
     */
    if( s->isactivesetchanged||s->issecondarytermchanged )
    {
        
        /*
         * Split variables into two groups - free (F) and constrained (C). Reorder
         * variables in such way that free vars come first, constrained are last:
         * x = [xf, xc].
         * 
         * Secondary term theta*(Q*x-r)'*(Q*x-r) now splits into quadratic part,
         * linear part and constant part:
         *             (          ( xf )     )'  (          ( xf )     )
         *   0.5*theta*( (Qf Qc)'*(    ) - r ) * ( (Qf Qc)'*(    ) - r ) =
         *             (          ( xc )     )   (          ( xc )     )
         *
         *   = 0.5*theta*xf'*(Qf'*Qf)*xf + theta*((Qc*xc-r)'*Qf)*xf + 
         *     + theta*(-r'*(Qc*xc-r)-0.5*r'*r+0.5*xc'*Qc'*Qc*xc)
         *                    
         * We store these parts into temporary variables:
         * * sqrt(theta)*Qf is stored into TK2
         * * theta*((Qc*xc-r)'*Qf) is stored into TK1
         * * theta*(-r'*(Qc*xc-r)-0.5*r'*r+0.5*xc'*Qc'*Qc*xc) is stored into TK0
         *
         * We use several other temporaries to store intermediate results:
         * * Tmp0 - to store Qc*xc-r
         * * Tmp1 - to store Qc*xc
         *
         * Generation of TK2/TK1/TK0 is performed as follows:
         * * we fill TK2/TK1/TK0 (to handle K=0 or Theta=0)
         * * other steps are performed only for K>0 and Theta>0
         * * we pass through columns of Q and copy I-th column into left block (Qf) or
         *   right one (Qc) of TK2, depending on presence of X[i] in the active set.
         *   CIdx0 variable contains current position for insertion into upper block,
         *   CIdx1 contains current position for insertion into lower one.
         * * we calculate Qc*xc-r and store it into Tmp0
         * * we calculate TK0 and TK1
         * * we multiply leading part of TK2 which stores Qf by sqrt(theta)
         *   it is important to perform this step AFTER calculation of TK0 and TK1,
         *   because we need original (non-modified) Qf to calculate TK0 and TK1.
         */
        for(j=0; j<=n-1; j++)
        {
            for(i=0; i<=k-1; i++)
            {
                s->tk2.ptr.pp_double[i][j] = 0.0;
            }
            s->tk1.ptr.p_double[j] = 0.0;
        }
        s->tk0 = 0.0;
        if( s->k>0&&ae_fp_greater(s->theta,0) )
        {
            
            /*
             * Split Q into Qf and Qc
             * Calculate Qc*xc-r, store in Tmp0
             */
            rvectorsetlengthatleast(&s->tmp0, k, _state);
            rvectorsetlengthatleast(&s->tmp1, k, _state);
            cidx0 = 0;
            cidx1 = nfree;
            for(i=0; i<=k-1; i++)
            {
                s->tmp1.ptr.p_double[i] = 0.0;
            }
            for(j=0; j<=n-1; j++)
            {
                if( s->activeset.ptr.p_bool[j] )
                {
                    for(i=0; i<=k-1; i++)
                    {
                        s->tk2.ptr.pp_double[i][cidx1] = s->q.ptr.pp_double[i][j];
                        s->tmp1.ptr.p_double[i] = s->tmp1.ptr.p_double[i]+s->q.ptr.pp_double[i][j]*s->txc.ptr.p_double[cidx1];
                    }
                    cidx1 = cidx1+1;
                }
                else
                {
                    for(i=0; i<=k-1; i++)
                    {
                        s->tk2.ptr.pp_double[i][cidx0] = s->q.ptr.pp_double[i][j];
                    }
                    cidx0 = cidx0+1;
                }
            }
            for(i=0; i<=k-1; i++)
            {
                s->tmp0.ptr.p_double[i] = s->tmp1.ptr.p_double[i]-s->r.ptr.p_double[i];
            }
            
            /*
             * Calculate TK0
             */
            v = 0.0;
            for(i=0; i<=k-1; i++)
            {
                v = v+s->theta*(0.5*ae_sqr(s->tmp1.ptr.p_double[i], _state)-s->r.ptr.p_double[i]*s->tmp0.ptr.p_double[i]-0.5*ae_sqr(s->r.ptr.p_double[i], _state));
            }
            s->tk0 = v;
            
            /*
             * Calculate TK1
             */
            if( nfree>0 )
            {
                for(i=0; i<=k-1; i++)
                {
                    v = s->theta*s->tmp0.ptr.p_double[i];
                    ae_v_addd(&s->tk1.ptr.p_double[0], 1, &s->tk2.ptr.pp_double[i][0], 1, ae_v_len(0,nfree-1), v);
                }
            }
            
            /*
             * Calculate TK2
             */
            if( nfree>0 )
            {
                v = ae_sqrt(s->theta, _state);
                for(i=0; i<=k-1; i++)
                {
                    ae_v_muld(&s->tk2.ptr.pp_double[i][0], 1, ae_v_len(0,nfree-1), v);
                }
            }
        }
    }
    
    /*
     * Re-evaluate TB
     */
    if( s->isactivesetchanged||s->islineartermchanged )
    {
        ridx0 = 0;
        ridx1 = nfree;
        for(i=0; i<=n-1; i++)
        {
            if( s->activeset.ptr.p_bool[i] )
            {
                s->tb.ptr.p_double[ridx1] = s->b.ptr.p_double[i];
                ridx1 = ridx1+1;
            }
            else
            {
                s->tb.ptr.p_double[ridx0] = s->b.ptr.p_double[i];
                ridx0 = ridx0+1;
            }
        }
    }
    
    /*
     * Compose ECA: either dense ECA or diagonal ECA
     */
    if( (s->isactivesetchanged||s->ismaintermchanged)&&nfree>0 )
    {
        if( ae_fp_greater(s->alpha,0) )
        {
            
            /*
             * Dense ECA
             */
            s->ecakind = 0;
            for(i=0; i<=nfree-1; i++)
            {
                for(j=i; j<=nfree-1; j++)
                {
                    s->ecadense.ptr.pp_double[i][j] = s->tq2dense.ptr.pp_double[i][j];
                }
            }
            if( !spdmatrixcholeskyrec(&s->ecadense, 0, nfree, ae_true, &s->tmp0, _state) )
            {
                result = ae_false;
                return result;
            }
        }
        else
        {
            
            /*
             * Diagonal ECA
             */
            s->ecakind = 1;
            for(i=0; i<=nfree-1; i++)
            {
                if( ae_fp_less(s->tq2diag.ptr.p_double[i],0) )
                {
                    result = ae_false;
                    return result;
                }
                s->ecadiag.ptr.p_double[i] = ae_sqrt(s->tq2diag.ptr.p_double[i], _state);
            }
        }
    }
    
    /*
     * Compose EQ
     */
    if( s->isactivesetchanged||s->issecondarytermchanged )
    {
        for(i=0; i<=k-1; i++)
        {
            for(j=0; j<=nfree-1; j++)
            {
                s->eq.ptr.pp_double[i][j] = s->tk2.ptr.pp_double[i][j];
            }
        }
    }
    
    /*
     * Calculate ECCM
     */
    if( ((((s->isactivesetchanged||s->ismaintermchanged)||s->issecondarytermchanged)&&s->k>0)&&ae_fp_greater(s->theta,0))&&nfree>0 )
    {
        
        /*
         * Calculate ECCM - Cholesky factor of the "effective" capacitance
         * matrix CM = I + EQ*inv(EffectiveA)*EQ'.
         *
         * We calculate CM as follows:
         *   CM = I + EQ*inv(EffectiveA)*EQ'
         *      = I + EQ*ECA^(-1)*ECA^(-T)*EQ'
         *      = I + (EQ*ECA^(-1))*(EQ*ECA^(-1))'
         *
         * Then we perform Cholesky decomposition of CM.
         */
        rmatrixsetlengthatleast(&s->tmp2, k, n, _state);
        rmatrixcopy(k, nfree, &s->eq, 0, 0, &s->tmp2, 0, 0, _state);
        ae_assert(s->ecakind==0||s->ecakind==1, "CQMRebuild: unexpected ECAKind", _state);
        if( s->ecakind==0 )
        {
            rmatrixrighttrsm(k, nfree, &s->ecadense, 0, 0, ae_true, ae_false, 0, &s->tmp2, 0, 0, _state);
        }
        if( s->ecakind==1 )
        {
            for(i=0; i<=k-1; i++)
            {
                for(j=0; j<=nfree-1; j++)
                {
                    s->tmp2.ptr.pp_double[i][j] = s->tmp2.ptr.pp_double[i][j]/s->ecadiag.ptr.p_double[j];
                }
            }
        }
        for(i=0; i<=k-1; i++)
        {
            for(j=0; j<=k-1; j++)
            {
                s->eccm.ptr.pp_double[i][j] = 0.0;
            }
            s->eccm.ptr.pp_double[i][i] = 1.0;
        }
        rmatrixsyrk(k, nfree, 1.0, &s->tmp2, 0, 0, 0, 1.0, &s->eccm, 0, 0, ae_true, _state);
        if( !spdmatrixcholeskyrec(&s->eccm, 0, k, ae_true, &s->tmp0, _state) )
        {
            result = ae_false;
            return result;
        }
    }
    
    /*
     * Compose EB and EC
     *
     * NOTE: because these quantities are cheap to compute, we do not
     * use caching here.
     */
    for(i=0; i<=nfree-1; i++)
    {
        s->eb.ptr.p_double[i] = s->tq1.ptr.p_double[i]+s->tk1.ptr.p_double[i]+s->tb.ptr.p_double[i];
    }
    s->ec = s->tq0+s->tk0;
    for(i=nfree; i<=n-1; i++)
    {
        s->ec = s->ec+s->tb.ptr.p_double[i]*s->txc.ptr.p_double[i];
    }
    
    /*
     * Change cache status - everything is cached 
     */
    s->ismaintermchanged = ae_false;
    s->issecondarytermchanged = ae_false;
    s->islineartermchanged = ae_false;
    s->isactivesetchanged = ae_false;
    return result;
}


/*************************************************************************
Internal function, solves system Effective_A*x = b.
It should be called after successful completion of CQMRebuild().

INPUT PARAMETERS:
    S       -   quadratic model, after call to CQMRebuild()
    X       -   right part B, array[S.NFree]
    Tmp     -   temporary array, automatically reallocated if needed

OUTPUT PARAMETERS:
    X       -   solution, array[S.NFree]
    
NOTE: when called with zero S.NFree, returns silently
NOTE: this function assumes that EA is non-degenerate

  -- ALGLIB --
     Copyright 10.05.2011 by Bochkanov Sergey
*************************************************************************/
static void cqmodels_cqmsolveea(convexquadraticmodel* s,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* tmp,
     ae_state *_state)
{
    ae_int_t i;


    ae_assert((s->ecakind==0||s->ecakind==1)||(s->ecakind==-1&&s->nfree==0), "CQMSolveEA: unexpected ECAKind", _state);
    if( s->ecakind==0 )
    {
        
        /*
         * Dense ECA, use FBLSCholeskySolve() dense solver.
         */
        fblscholeskysolve(&s->ecadense, 1.0, s->nfree, ae_true, x, tmp, _state);
    }
    if( s->ecakind==1 )
    {
        
        /*
         * Diagonal ECA
         */
        for(i=0; i<=s->nfree-1; i++)
        {
            x->ptr.p_double[i] = x->ptr.p_double[i]/ae_sqr(s->ecadiag.ptr.p_double[i], _state);
        }
    }
}


ae_bool _convexquadraticmodel_init(convexquadraticmodel* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_matrix_init(&p->a, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->q, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->b, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->r, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xc, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->d, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->activeset, 0, DT_BOOL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->tq2dense, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->tk2, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tq2diag, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tq1, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tk1, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->txc, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tb, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->ecadense, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->eq, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->eccm, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->ecadiag, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->eb, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tmp0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tmp1, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tmpg, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->tmp2, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _convexquadraticmodel_init_copy(convexquadraticmodel* dst, convexquadraticmodel* src, ae_state *_state, ae_bool make_automatic)
{
    dst->n = src->n;
    dst->k = src->k;
    dst->alpha = src->alpha;
    dst->tau = src->tau;
    dst->theta = src->theta;
    if( !ae_matrix_init_copy(&dst->a, &src->a, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->q, &src->q, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->b, &src->b, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->r, &src->r, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->xc, &src->xc, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->d, &src->d, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->activeset, &src->activeset, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->tq2dense, &src->tq2dense, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->tk2, &src->tk2, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tq2diag, &src->tq2diag, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tq1, &src->tq1, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tk1, &src->tk1, _state, make_automatic) )
        return ae_false;
    dst->tq0 = src->tq0;
    dst->tk0 = src->tk0;
    if( !ae_vector_init_copy(&dst->txc, &src->txc, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tb, &src->tb, _state, make_automatic) )
        return ae_false;
    dst->nfree = src->nfree;
    dst->ecakind = src->ecakind;
    if( !ae_matrix_init_copy(&dst->ecadense, &src->ecadense, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->eq, &src->eq, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->eccm, &src->eccm, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->ecadiag, &src->ecadiag, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->eb, &src->eb, _state, make_automatic) )
        return ae_false;
    dst->ec = src->ec;
    if( !ae_vector_init_copy(&dst->tmp0, &src->tmp0, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tmp1, &src->tmp1, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tmpg, &src->tmpg, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->tmp2, &src->tmp2, _state, make_automatic) )
        return ae_false;
    dst->ismaintermchanged = src->ismaintermchanged;
    dst->issecondarytermchanged = src->issecondarytermchanged;
    dst->islineartermchanged = src->islineartermchanged;
    dst->isactivesetchanged = src->isactivesetchanged;
    return ae_true;
}


void _convexquadraticmodel_clear(convexquadraticmodel* p)
{
    ae_matrix_clear(&p->a);
    ae_matrix_clear(&p->q);
    ae_vector_clear(&p->b);
    ae_vector_clear(&p->r);
    ae_vector_clear(&p->xc);
    ae_vector_clear(&p->d);
    ae_vector_clear(&p->activeset);
    ae_matrix_clear(&p->tq2dense);
    ae_matrix_clear(&p->tk2);
    ae_vector_clear(&p->tq2diag);
    ae_vector_clear(&p->tq1);
    ae_vector_clear(&p->tk1);
    ae_vector_clear(&p->txc);
    ae_vector_clear(&p->tb);
    ae_matrix_clear(&p->ecadense);
    ae_matrix_clear(&p->eq);
    ae_matrix_clear(&p->eccm);
    ae_vector_clear(&p->ecadiag);
    ae_vector_clear(&p->eb);
    ae_vector_clear(&p->tmp0);
    ae_vector_clear(&p->tmp1);
    ae_vector_clear(&p->tmpg);
    ae_matrix_clear(&p->tmp2);
}




/*************************************************************************
                    CONSTRAINED QUADRATIC PROGRAMMING

The subroutine creates QP optimizer. After initial creation,  it  contains
default optimization problem with zero quadratic and linear terms  and  no
constraints. You should set quadratic/linear terms with calls to functions
provided by MinQP subpackage.

INPUT PARAMETERS:
    N       -   problem size
    
OUTPUT PARAMETERS:
    State   -   optimizer with zero quadratic/linear terms
                and no constraints

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpcreate(ae_int_t n, minqpstate* state, ae_state *_state)
{
    ae_int_t i;

    _minqpstate_clear(state);

    ae_assert(n>=1, "MinQPCreate: N<1", _state);
    
    /*
     * initialize QP solver
     */
    state->n = n;
    state->ccount = 0;
    state->nec = 0;
    state->nic = 0;
    state->activelincnt = 0;
    state->repterminationtype = 0;
    state->anorm = 1;
    cqminit(n, &state->a, _state);
    ae_vector_set_length(&state->b, n, _state);
    ae_vector_set_length(&state->bndl, n, _state);
    ae_vector_set_length(&state->bndu, n, _state);
    ae_vector_set_length(&state->workbndl, n, _state);
    ae_vector_set_length(&state->workbndu, n, _state);
    ae_vector_set_length(&state->havebndl, n, _state);
    ae_vector_set_length(&state->havebndu, n, _state);
    ae_vector_set_length(&state->startx, n, _state);
    ae_vector_set_length(&state->xorigin, n, _state);
    ae_vector_set_length(&state->xc, n, _state);
    ae_vector_set_length(&state->xn, n, _state);
    ae_vector_set_length(&state->gc, n, _state);
    ae_vector_set_length(&state->pg, n, _state);
    ae_vector_set_length(&state->activeb, n, _state);
    for(i=0; i<=n-1; i++)
    {
        state->bndl.ptr.p_double[i] = _state->v_neginf;
        state->bndu.ptr.p_double[i] = _state->v_posinf;
        state->havebndl.ptr.p_bool[i] = ae_false;
        state->havebndu.ptr.p_bool[i] = ae_false;
        state->b.ptr.p_double[i] = 0.0;
        state->startx.ptr.p_double[i] = 0.0;
        state->xorigin.ptr.p_double[i] = 0.0;
    }
    state->havex = ae_false;
    minqpsetalgocholesky(state, _state);
    normestimatorcreate(n, n, 5, 5, &state->estimator, _state);
}


/*************************************************************************
This function sets linear term for QP solver.

By default, linear term is zero.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    B       -   linear term, array[N].

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpsetlinearterm(minqpstate* state,
     /* Real    */ ae_vector* b,
     ae_state *_state)
{
    ae_int_t n;


    n = state->n;
    ae_assert(b->cnt>=n, "MinQPSetLinearTerm: Length(B)<N", _state);
    ae_assert(isfinitevector(b, n, _state), "MinQPSetLinearTerm: B contains infinite or NaN elements", _state);
    minqpsetlineartermfast(state, b, _state);
}


/*************************************************************************
This function sets quadratic term for QP solver.

By default quadratic term is zero.

IMPORTANT: this solver minimizes following  function:
    f(x) = 0.5*x'*A*x + b'*x.
Note that quadratic term has 0.5 before it. So if  you  want  to  minimize
    f(x) = x^2 + x
you should rewrite your problem as follows:
    f(x) = 0.5*(2*x^2) + x
and your matrix A will be equal to [[2.0]], not to [[1.0]]

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    A       -   matrix, array[N,N]
    IsUpper -   (optional) storage type:
                * if True, symmetric matrix  A  is  given  by  its  upper
                  triangle, and the lower triangle isn’t used
                * if False, symmetric matrix  A  is  given  by  its lower
                  triangle, and the upper triangle isn’t used
                * if not given, both lower and upper  triangles  must  be
                  filled.

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpsetquadraticterm(minqpstate* state,
     /* Real    */ ae_matrix* a,
     ae_bool isupper,
     ae_state *_state)
{
    ae_int_t n;


    n = state->n;
    ae_assert(a->rows>=n, "MinQPSetQuadraticTerm: Rows(A)<N", _state);
    ae_assert(a->cols>=n, "MinQPSetQuadraticTerm: Cols(A)<N", _state);
    ae_assert(isfinitertrmatrix(a, n, isupper, _state), "MinQPSetQuadraticTerm: A contains infinite or NaN elements", _state);
    minqpsetquadratictermfast(state, a, isupper, 0.0, _state);
}


/*************************************************************************
This function sets starting point for QP solver. It is useful to have
good initial approximation to the solution, because it will increase
speed of convergence and identification of active constraints.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    X       -   starting point, array[N].

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpsetstartingpoint(minqpstate* state,
     /* Real    */ ae_vector* x,
     ae_state *_state)
{
    ae_int_t n;


    n = state->n;
    ae_assert(x->cnt>=n, "MinQPSetStartingPoint: Length(B)<N", _state);
    ae_assert(isfinitevector(x, n, _state), "MinQPSetStartingPoint: X contains infinite or NaN elements", _state);
    minqpsetstartingpointfast(state, x, _state);
}


/*************************************************************************
This  function sets origin for QP solver. By default, following QP program
is solved:

    min(0.5*x'*A*x+b'*x)
    
This function allows to solve different problem:

    min(0.5*(x-x_origin)'*A*(x-x_origin)+b'*(x-x_origin))
    
INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    XOrigin -   origin, array[N].

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpsetorigin(minqpstate* state,
     /* Real    */ ae_vector* xorigin,
     ae_state *_state)
{
    ae_int_t n;


    n = state->n;
    ae_assert(xorigin->cnt>=n, "MinQPSetOrigin: Length(B)<N", _state);
    ae_assert(isfinitevector(xorigin, n, _state), "MinQPSetOrigin: B contains infinite or NaN elements", _state);
    minqpsetoriginfast(state, xorigin, _state);
}


/*************************************************************************
This function tells solver to use Cholesky-based algorithm.

Cholesky-based algorithm can be used when:
* problem is convex
* there is no constraints or only boundary constraints are present

This algorithm has O(N^3) complexity for unconstrained problem and  is  up
to several times slower on bound constrained  problems  (these  additional
iterations are needed to identify active constraints).

INPUT PARAMETERS:
    State   -   structure which stores algorithm state

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpsetalgocholesky(minqpstate* state, ae_state *) //_state)
{


    state->algokind = 1;
}


/*************************************************************************
This function sets boundary constraints for QP solver

Boundary constraints are inactive by default (after initial creation).
After  being  set,  they  are  preserved  until explicitly turned off with
another SetBC() call.

INPUT PARAMETERS:
    State   -   structure stores algorithm state
    BndL    -   lower bounds, array[N].
                If some (all) variables are unbounded, you may specify
                very small number or -INF (latter is recommended because
                it will allow solver to use better algorithm).
    BndU    -   upper bounds, array[N].
                If some (all) variables are unbounded, you may specify
                very large number or +INF (latter is recommended because
                it will allow solver to use better algorithm).
                
NOTE: it is possible to specify BndL[i]=BndU[i]. In this case I-th
variable will be "frozen" at X[i]=BndL[i]=BndU[i].

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpsetbc(minqpstate* state,
     /* Real    */ ae_vector* bndl,
     /* Real    */ ae_vector* bndu,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t n;


    n = state->n;
    ae_assert(bndl->cnt>=n, "MinQPSetBC: Length(BndL)<N", _state);
    ae_assert(bndu->cnt>=n, "MinQPSetBC: Length(BndU)<N", _state);
    for(i=0; i<=n-1; i++)
    {
        ae_assert(ae_isfinite(bndl->ptr.p_double[i], _state)||ae_isneginf(bndl->ptr.p_double[i], _state), "MinQPSetBC: BndL contains NAN or +INF", _state);
        ae_assert(ae_isfinite(bndu->ptr.p_double[i], _state)||ae_isposinf(bndu->ptr.p_double[i], _state), "MinQPSetBC: BndU contains NAN or -INF", _state);
        state->bndl.ptr.p_double[i] = bndl->ptr.p_double[i];
        state->havebndl.ptr.p_bool[i] = ae_isfinite(bndl->ptr.p_double[i], _state);
        state->bndu.ptr.p_double[i] = bndu->ptr.p_double[i];
        state->havebndu.ptr.p_bool[i] = ae_isfinite(bndu->ptr.p_double[i], _state);
    }
}


/*************************************************************************
This function sets linear constraints for QP optimizer.

Linear constraints are inactive by default (after initial creation).

INPUT PARAMETERS:
    State   -   structure previously allocated with MinQPCreate call.
    C       -   linear constraints, array[K,N+1].
                Each row of C represents one constraint, either equality
                or inequality (see below):
                * first N elements correspond to coefficients,
                * last element corresponds to the right part.
                All elements of C (including right part) must be finite.
    CT      -   type of constraints, array[K]:
                * if CT[i]>0, then I-th constraint is C[i,*]*x >= C[i,n+1]
                * if CT[i]=0, then I-th constraint is C[i,*]*x  = C[i,n+1]
                * if CT[i]<0, then I-th constraint is C[i,*]*x <= C[i,n+1]
    K       -   number of equality/inequality constraints, K>=0:
                * if given, only leading K elements of C/CT are used
                * if not given, automatically determined from sizes of C/CT

NOTE 1: linear (non-bound) constraints are satisfied only approximately  -
        there always exists some minor violation (about 10^10...10^13) due
        to rounding errors.

  -- ALGLIB --
     Copyright 19.06.2012 by Bochkanov Sergey
*************************************************************************/
void minqpsetlc(minqpstate* state,
     /* Real    */ ae_matrix* c,
     /* Integer */ ae_vector* ct,
     ae_int_t k,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t i;
    ae_int_t j;
    double v;


    n = state->n;
    
    /*
     * First, check for errors in the inputs
     */
    ae_assert(k>=0, "MinQPSetLC: K<0", _state);
    ae_assert(c->cols>=n+1||k==0, "MinQPSetLC: Cols(C)<N+1", _state);
    ae_assert(c->rows>=k, "MinQPSetLC: Rows(C)<K", _state);
    ae_assert(ct->cnt>=k, "MinQPSetLC: Length(CT)<K", _state);
    ae_assert(apservisfinitematrix(c, k, n+1, _state), "MinQPSetLC: C contains infinite or NaN values!", _state);
    
    /*
     * Handle zero K
     */
    if( k==0 )
    {
        state->nec = 0;
        state->nic = 0;
        return;
    }
    
    /*
     * Equality constraints are stored first, in the upper
     * NEC rows of State.CLEIC matrix. Inequality constraints
     * are stored in the next NIC rows.
     *
     * NOTE: we convert inequality constraints to the form
     * A*x<=b before copying them.
     */
    rmatrixsetlengthatleast(&state->cleic, k, n+1, _state);
    state->nec = 0;
    state->nic = 0;
    for(i=0; i<=k-1; i++)
    {
        if( ct->ptr.p_int[i]==0 )
        {
            ae_v_move(&state->cleic.ptr.pp_double[state->nec][0], 1, &c->ptr.pp_double[i][0], 1, ae_v_len(0,n));
            state->nec = state->nec+1;
        }
    }
    for(i=0; i<=k-1; i++)
    {
        if( ct->ptr.p_int[i]!=0 )
        {
            if( ct->ptr.p_int[i]>0 )
            {
                ae_v_moveneg(&state->cleic.ptr.pp_double[state->nec+state->nic][0], 1, &c->ptr.pp_double[i][0], 1, ae_v_len(0,n));
            }
            else
            {
                ae_v_move(&state->cleic.ptr.pp_double[state->nec+state->nic][0], 1, &c->ptr.pp_double[i][0], 1, ae_v_len(0,n));
            }
            state->nic = state->nic+1;
        }
    }
    
    /*
     * Normalize rows of State.CLEIC: each row must have unit norm.
     * Norm is calculated using first N elements (i.e. right part is
     * not counted when we calculate norm).
     */
    for(i=0; i<=k-1; i++)
    {
        v = 0;
        for(j=0; j<=n-1; j++)
        {
            v = v+ae_sqr(state->cleic.ptr.pp_double[i][j], _state);
        }
        if( ae_fp_eq(v,0) )
        {
            continue;
        }
        v = 1/ae_sqrt(v, _state);
        ae_v_muld(&state->cleic.ptr.pp_double[i][0], 1, ae_v_len(0,n), v);
    }
    
    /*
     * Determine number of constraints,
     * allocate space and copy
     */
    state->ccount = k;
    rmatrixsetlengthatleast(&state->cmatrix, k, n+1, _state);
    ivectorsetlengthatleast(&state->ct, k, _state);
    rvectorsetlengthatleast(&state->cr, k, _state);
    for(i=0; i<=k-1; i++)
    {
        state->ct.ptr.p_int[i] = ct->ptr.p_int[i];
        state->cr.ptr.p_double[i] = c->ptr.pp_double[i][n];
        ae_v_move(&state->cmatrix.ptr.pp_double[i][0], 1, &c->ptr.pp_double[i][0], 1, ae_v_len(0,n));
    }
}


/*************************************************************************
This function solves quadratic programming problem.
You should call it after setting solver options with MinQPSet...() calls.

INPUT PARAMETERS:
    State   -   algorithm state

You should use MinQPResults() function to access results after calls
to this function.

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpoptimize(minqpstate* state, ae_state *_state)
{
    ae_int_t n;
    ae_int_t i;
    ae_int_t j;
    ae_int_t nbc;
    double v0;
    double v1;
    double vn0;
    double v;
    double noisetolerance;
    double fprev;
    double fcur;
    ae_int_t baditscounter;


    noisetolerance = 10;
    n = state->n;
    state->repterminationtype = -5;
    state->repinneriterationscount = 0;
    state->repouteriterationscount = 0;
    state->repncholesky = 0;
    state->repnmv = 0;
    
    /*
     * check correctness of constraints
     */
    for(i=0; i<=n-1; i++)
    {
        if( state->havebndl.ptr.p_bool[i]&&state->havebndu.ptr.p_bool[i] )
        {
            if( ae_fp_greater(state->bndl.ptr.p_double[i],state->bndu.ptr.p_double[i]) )
            {
                state->repterminationtype = -3;
                return;
            }
        }
    }
    
    /*
     * count number of bound and linear constraints
     */
    nbc = 0;
    for(i=0; i<=n-1; i++)
    {
        if( state->havebndl.ptr.p_bool[i] )
        {
            nbc = nbc+1;
        }
        if( state->havebndu.ptr.p_bool[i] )
        {
            nbc = nbc+1;
        }
    }
    
    /*
     * Our formulation of quadratic problem includes origin point,
     * i.e. we have F(x-x_origin) which is minimized subject to
     * constraints on x, instead of having simply F(x).
     *
     * Here we make transition from non-zero origin to zero one.
     * In order to make such transition we have to:
     * 1. subtract x_origin from x_start
     * 2. modify constraints
     * 3. solve problem
     * 4. add x_origin to solution
     *
     * There is alternate solution - to modify quadratic function
     * by expansion of multipliers containing (x-x_origin), but
     * we prefer to modify constraints, because it is a) more precise
     * and b) easier to to.
     *
     * Parts (1)-(2) are done here. After this block is over,
     * we have:
     * * XC, which stores shifted XStart (if we don't have XStart,
     *   value of XC will be ignored later)
     * * WorkBndL, WorkBndU, which store modified boundary constraints.
     */
    for(i=0; i<=n-1; i++)
    {
        if( state->havebndl.ptr.p_bool[i] )
        {
            state->workbndl.ptr.p_double[i] = state->bndl.ptr.p_double[i]-state->xorigin.ptr.p_double[i];
        }
        else
        {
            state->workbndl.ptr.p_double[i] = _state->v_neginf;
        }
        if( state->havebndu.ptr.p_bool[i] )
        {
            state->workbndu.ptr.p_double[i] = state->bndu.ptr.p_double[i]-state->xorigin.ptr.p_double[i];
        }
        else
        {
            state->workbndu.ptr.p_double[i] = _state->v_posinf;
        }
    }
    rmatrixsetlengthatleast(&state->workcleic, state->nec+state->nic, n+1, _state);
    for(i=0; i<=state->nec+state->nic-1; i++)
    {
        v = ae_v_dotproduct(&state->cleic.ptr.pp_double[i][0], 1, &state->xorigin.ptr.p_double[0], 1, ae_v_len(0,n-1));
        ae_v_move(&state->workcleic.ptr.pp_double[i][0], 1, &state->cleic.ptr.pp_double[i][0], 1, ae_v_len(0,n-1));
        state->workcleic.ptr.pp_double[i][n] = state->cleic.ptr.pp_double[i][n]-v;
    }
    
    /*
     * modify starting point XC according to boundary constraints
     */
    if( state->havex )
    {
        
        /*
         * We have starting point in StartX, so we just have to shift and bound it
         */
        for(i=0; i<=n-1; i++)
        {
            state->xc.ptr.p_double[i] = state->startx.ptr.p_double[i]-state->xorigin.ptr.p_double[i];
            if( state->havebndl.ptr.p_bool[i] )
            {
                if( ae_fp_less(state->xc.ptr.p_double[i],state->workbndl.ptr.p_double[i]) )
                {
                    state->xc.ptr.p_double[i] = state->workbndl.ptr.p_double[i];
                }
            }
            if( state->havebndu.ptr.p_bool[i] )
            {
                if( ae_fp_greater(state->xc.ptr.p_double[i],state->workbndu.ptr.p_double[i]) )
                {
                    state->xc.ptr.p_double[i] = state->workbndu.ptr.p_double[i];
                }
            }
        }
    }
    else
    {
        
        /*
         * We don't have starting point, so we deduce it from
         * constraints (if they are present).
         *
         * NOTE: XC contains some meaningless values from previous block
         * which are ignored by code below.
         */
        for(i=0; i<=n-1; i++)
        {
            if( state->havebndl.ptr.p_bool[i]&&state->havebndu.ptr.p_bool[i] )
            {
                state->xc.ptr.p_double[i] = 0.5*(state->workbndl.ptr.p_double[i]+state->workbndu.ptr.p_double[i]);
                if( ae_fp_less(state->xc.ptr.p_double[i],state->workbndl.ptr.p_double[i]) )
                {
                    state->xc.ptr.p_double[i] = state->workbndl.ptr.p_double[i];
                }
                if( ae_fp_greater(state->xc.ptr.p_double[i],state->workbndu.ptr.p_double[i]) )
                {
                    state->xc.ptr.p_double[i] = state->workbndu.ptr.p_double[i];
                }
                continue;
            }
            if( state->havebndl.ptr.p_bool[i] )
            {
                state->xc.ptr.p_double[i] = state->workbndl.ptr.p_double[i];
                continue;
            }
            if( state->havebndu.ptr.p_bool[i] )
            {
                state->xc.ptr.p_double[i] = state->workbndu.ptr.p_double[i];
                continue;
            }
            state->xc.ptr.p_double[i] = 0;
        }
    }
    
    /*
     * Select algo
     */
    if( state->algokind==1 )
    {
        
        /*
         * Cholesky-based algorithm for dense bound constrained problems.
         *
         * This algorithm exists in two variants:
         * * unconstrained one, which can solve problem using only one NxN
         *   double matrix
         * * bound constrained one, which needs two NxN matrices
         *
         * We will try to solve problem using unconstrained algorithm,
         * and will use bound constrained version only when constraints
         * are actually present
         */
        if( nbc==0&&state->ccount==0 )
        {
            
            /*
             * "Simple" unconstrained version
             */
            for(i=0; i<=n-1; i++)
            {
                state->activeb.ptr.p_bool[i] = ae_false;
            }
            state->repncholesky = state->repncholesky+1;
            cqmsetb(&state->a, &state->b, _state);
            cqmsetactiveset(&state->a, &state->xc, &state->activeb, _state);
            if( !cqmconstrainedoptimum(&state->a, &state->xn, _state) )
            {
                state->repterminationtype = -5;
                return;
            }
            ae_v_move(&state->xc.ptr.p_double[0], 1, &state->xn.ptr.p_double[0], 1, ae_v_len(0,n-1));
            ae_v_add(&state->xc.ptr.p_double[0], 1, &state->xorigin.ptr.p_double[0], 1, ae_v_len(0,n-1));
            state->repinneriterationscount = 1;
            state->repouteriterationscount = 1;
            state->repterminationtype = 4;
            return;
        }
        
        /*
         * Determine initial set of active boundary constraints
         */
        if( state->nec+state->nic>0 )
        {
            
            /*
             * General linear constraints are present; general code is used.
             */
            rvectorsetlengthatleast(&state->tmp0, n, _state);
            rvectorsetlengthatleast(&state->tmpfeas, n+state->nic, _state);
            rmatrixsetlengthatleast(&state->tmpm0, state->nec+state->nic, n+state->nic+1, _state);
            for(i=0; i<=state->nec+state->nic-1; i++)
            {
                ae_v_move(&state->tmpm0.ptr.pp_double[i][0], 1, &state->workcleic.ptr.pp_double[i][0], 1, ae_v_len(0,n-1));
                for(j=n; j<=n+state->nic-1; j++)
                {
                    state->tmpm0.ptr.pp_double[i][j] = 0;
                }
                if( i>=state->nec )
                {
                    state->tmpm0.ptr.pp_double[i][n+i-state->nec] = 1.0;
                }
                state->tmpm0.ptr.pp_double[i][n+state->nic] = state->workcleic.ptr.pp_double[i][n];
            }
            ae_v_move(&state->tmpfeas.ptr.p_double[0], 1, &state->xc.ptr.p_double[0], 1, ae_v_len(0,n-1));
            for(i=0; i<=state->nic-1; i++)
            {
                v = ae_v_dotproduct(&state->workcleic.ptr.pp_double[i+state->nec][0], 1, &state->xc.ptr.p_double[0], 1, ae_v_len(0,n-1));
                state->tmpfeas.ptr.p_double[i+n] = ae_maxreal(state->workcleic.ptr.pp_double[i+state->nec][n]-v, 0.0, _state);
            }
            if( !findfeasiblepoint(&state->tmpfeas, &state->workbndl, &state->havebndl, &state->workbndu, &state->havebndu, n, state->nic, &state->tmpm0, state->nec+state->nic, 1.0E-6, &i, &j, _state) )
            {
                state->repterminationtype = -3;
                return;
            }
            ae_v_move(&state->xc.ptr.p_double[0], 1, &state->tmpfeas.ptr.p_double[0], 1, ae_v_len(0,n-1));
            for(i=0; i<=n-1; i++)
            {
                state->activeb.ptr.p_bool[i] = (state->havebndl.ptr.p_bool[i]&&ae_fp_eq(state->xc.ptr.p_double[i],state->workbndl.ptr.p_double[i]))||(state->havebndu.ptr.p_bool[i]&&ae_fp_eq(state->xc.ptr.p_double[i],state->workbndu.ptr.p_double[i]));
            }
            bvectorsetlengthatleast(&state->activelin, state->nec+state->nic, _state);
            for(i=0; i<=state->nec+state->nic-1; i++)
            {
                state->activelin.ptr.p_bool[i] = i<state->nec||ae_fp_eq(state->tmpfeas.ptr.p_double[n+i-state->nec],0);
            }
        }
        else
        {
            
            /*
             * Only bound constraints are present, quick code can be used
             */
            for(i=0; i<=n-1; i++)
            {
                state->activeb.ptr.p_bool[i] = ae_false;
                if( state->havebndl.ptr.p_bool[i]&&ae_fp_less_eq(state->xc.ptr.p_double[i],state->workbndl.ptr.p_double[i]) )
                {
                    state->xc.ptr.p_double[i] = state->workbndl.ptr.p_double[i];
                    state->activeb.ptr.p_bool[i] = ae_true;
                }
                if( state->havebndu.ptr.p_bool[i]&&ae_fp_greater_eq(state->xc.ptr.p_double[i],state->workbndu.ptr.p_double[i]) )
                {
                    state->xc.ptr.p_double[i] = state->workbndu.ptr.p_double[i];
                    state->activeb.ptr.p_bool[i] = ae_true;
                }
            }
        }
        minqp_orthogonalizeactiveconstraints(&state->xc, &state->activeb, n, &state->workcleic, &state->activelin, state->nec+state->nic, &state->activecm, &state->activecr, &state->activelincnt, _state);
        
        /*
         * Main cycle of BLEIC-QP algorithm
         */
        baditscounter = 0;
        state->repterminationtype = 4;
        for(;;)
        {
            
            /*
             * Phase 1.
             * 
             * Active constraints can be added, but can not be removed.
             * This phase allows us to reduce function value without
             * leaving current face.
             *
             * Phase is over when we have no constraints to add. After
             * this phase we have optimal and strictly feasible XC
             * (or hopelessly infeasible - when constraints are inconsistent).
             *
             * Iterations which increase model value instead of decreasing
             * it are called "bad". When number of "bad" iterations exceeds
             * BadQPItsBeforeTermination limit, algorithm is terminated.
             */
            fprev = minqp_minqpmodelvalue(&state->a, &state->b, &state->xc, n, &state->tmp0, _state);
            for(;;)
            {
                
                /*
                 * Calculate optimum subject to presently active constraints
                 */
                state->repncholesky = state->repncholesky+1;
                if( !minqp_minqpconstrainedoptimum(&state->a, state->anorm, &state->b, &state->activeb, &state->xc, n, &state->activecm, &state->activecr, state->activelincnt, &state->xn, &state->tmp0, &state->tmp1, _state) )
                {
                    state->repterminationtype = -5;
                    return;
                }
                
                /*
                 * Add constraints.
                 * If no constraints were added, accept candidate point XN and move to next phase.
                 */
                if( minqp_minqpboundedstepandactivation(state, &state->xc, &state->xn, &state->tmp0, _state)==0 )
                {
                    break;
                }
            }
            fcur = minqp_minqpmodelvalue(&state->a, &state->b, &state->xc, n, &state->tmp0, _state);
            if( ae_fp_greater_eq(fcur,fprev) )
            {
                baditscounter = baditscounter+1;
            }
            if( baditscounter>=minqp_badqpitsbeforetermination )
            {
                state->repterminationtype = 7;
                break;
            }
            
            /*
             * Phase 2 - steepest descent.
             * The only phase where we can terminate algorithm.
             *
             * 1. we calculate projected gradient PG and scale it using
             *    inverse diagonal of Hessian matrix
             * 2. we calculate Cauchy point XN
             * 3. in case difference between f(XN) and f(XC) is dominated by
             *    numerical noise, we terminate algorithm.
             * 4. otherwise, we perform bounded step to XN
             *
             * This is the only phase were we can terminate algorithm
             *
             * NOTE: Cauchy point is calculated as minimizer of
             *       g(t) = f(xc-t*pg) = t^2*(0.5*pg'*A*pg) + t*(-pg*g) + const
             *       where pg is a projected gradient, t is a step length, g is
             *       a non-projected gradient.
             *
             * NOTE 2: here we assume that our problem is convex, hence pg'*A*pg>=0,
             *         and pg'*A*pg=0 if and only if pg=0.
             */
            minqp_minqpdropconstraints(state, &state->xc, &state->tmp0, _state);
            minqp_minqpgradientandconstraineddescent(state, &state->xc, &state->gc, &state->pg, _state);
            cqmevalx(&state->a, &state->xc, &v, &vn0, _state);
            v0 = cqmxtadx2(&state->a, &state->pg, _state);
            if( ae_fp_greater(v0,0) )
            {
                
                /*
                 * Function curvature is positive.
                 * Calculate predicted function decrease, which is equal to (pg'*g)^2/(4*pg'*A*pg).
                 * In case step is dominated by noise, we exit.
                 */
                v1 = ae_v_dotproduct(&state->pg.ptr.p_double[0], 1, &state->gc.ptr.p_double[0], 1, ae_v_len(0,n-1));
                v1 = -v1;
                if( ae_fp_less(ae_sqr(v1, _state)/(4*v0),noisetolerance*vn0) )
                {
                    break;
                }
                v = v1/(2*v0);
                ae_v_move(&state->xn.ptr.p_double[0], 1, &state->xc.ptr.p_double[0], 1, ae_v_len(0,n-1));
                ae_v_addd(&state->xn.ptr.p_double[0], 1, &state->pg.ptr.p_double[0], 1, ae_v_len(0,n-1), v);
                minqp_minqpboundedstepandactivation(state, &state->xc, &state->xn, &state->tmp0, _state);
            }
            else
            {
                
                /*
                 * non-positive curvature detected.
                 * because it is assumed that function is positive definite,
                 * it means that projected gradient is zero (termination is needed).
                 */
                break;
            }
            
            /*
             * Update iterations count
             */
            state->repinneriterationscount = state->repinneriterationscount+1;
        }
        state->repouteriterationscount = 1;
        
        /*
         * Post-process: add XOrigin to XC
         */
        for(i=0; i<=n-1; i++)
        {
            if( state->havebndl.ptr.p_bool[i]&&ae_fp_eq(state->xc.ptr.p_double[i],state->workbndl.ptr.p_double[i]) )
            {
                state->xc.ptr.p_double[i] = state->bndl.ptr.p_double[i];
                continue;
            }
            if( state->havebndu.ptr.p_bool[i]&&ae_fp_eq(state->xc.ptr.p_double[i],state->workbndu.ptr.p_double[i]) )
            {
                state->xc.ptr.p_double[i] = state->bndu.ptr.p_double[i];
                continue;
            }
            state->xc.ptr.p_double[i] = boundval(state->xc.ptr.p_double[i]+state->xorigin.ptr.p_double[i], state->bndl.ptr.p_double[i], state->bndu.ptr.p_double[i], _state);
        }
        return;
    }
}


/*************************************************************************
QP solver results

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    X       -   array[0..N-1], solution
    Rep     -   optimization report. You should check Rep.TerminationType,
                which contains completion code, and you may check  another
                fields which contain another information  about  algorithm
                functioning.

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpresults(minqpstate* state,
     /* Real    */ ae_vector* x,
     minqpreport* rep,
     ae_state *_state)
{

    ae_vector_clear(x);
    _minqpreport_clear(rep);

    minqpresultsbuf(state, x, rep, _state);
}


/*************************************************************************
QP results

Buffered implementation of MinQPResults() which uses pre-allocated  buffer
to store X[]. If buffer size is  too  small,  it  resizes  buffer.  It  is
intended to be used in the inner cycles of performance critical algorithms
where array reallocation penalty is too large to be ignored.

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpresultsbuf(minqpstate* state,
     /* Real    */ ae_vector* x,
     minqpreport* rep,
     ae_state *_state)
{


    if( x->cnt<state->n )
    {
        ae_vector_set_length(x, state->n, _state);
    }
    ae_v_move(&x->ptr.p_double[0], 1, &state->xc.ptr.p_double[0], 1, ae_v_len(0,state->n-1));
    rep->inneriterationscount = state->repinneriterationscount;
    rep->outeriterationscount = state->repouteriterationscount;
    rep->nmv = state->repnmv;
    rep->ncholesky = state->repncholesky;
    rep->terminationtype = state->repterminationtype;
}


/*************************************************************************
Fast version of MinQPSetLinearTerm(), which doesn't check its arguments.
For internal use only.

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpsetlineartermfast(minqpstate* state,
     /* Real    */ ae_vector* b,
     ae_state *) //_state)
{


    ae_v_move(&state->b.ptr.p_double[0], 1, &b->ptr.p_double[0], 1, ae_v_len(0,state->n-1));
}


/*************************************************************************
Fast version of MinQPSetQuadraticTerm(), which doesn't check its arguments.

It accepts additional parameter - shift S, which allows to "shift"  matrix
A by adding s*I to A. S must be positive (although it is not checked).

For internal use only.

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpsetquadratictermfast(minqpstate* state,
     /* Real    */ ae_matrix* a,
     ae_bool isupper,
     double s,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t n;


    n = state->n;
    cqmseta(&state->a, a, isupper, 1.0, _state);
    if( ae_fp_greater(s,0) )
    {
        rvectorsetlengthatleast(&state->tmp0, n, _state);
        for(i=0; i<=n-1; i++)
        {
            state->tmp0.ptr.p_double[i] = a->ptr.pp_double[i][i]+s;
        }
        cqmrewritedensediagonal(&state->a, &state->tmp0, _state);
    }
    
    /*
     * Estimate norm of A
     * (it will be used later in the quadratic penalty function)
     */
    state->anorm = 0;
    for(i=0; i<=n-1; i++)
    {
        if( isupper )
        {
            for(j=i; j<=n-1; j++)
            {
                state->anorm = ae_maxreal(state->anorm, ae_fabs(a->ptr.pp_double[i][j], _state), _state);
            }
        }
        else
        {
            for(j=0; j<=i; j++)
            {
                state->anorm = ae_maxreal(state->anorm, ae_fabs(a->ptr.pp_double[i][j], _state), _state);
            }
        }
    }
    state->anorm = state->anorm*n;
}


/*************************************************************************
Internal function which allows to rewrite diagonal of quadratic term.
For internal use only.

This function can be used only when you have dense A and already made
MinQPSetQuadraticTerm(Fast) call.

  -- ALGLIB --
     Copyright 16.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqprewritediagonal(minqpstate* state,
     /* Real    */ ae_vector* s,
     ae_state *_state)
{


    cqmrewritedensediagonal(&state->a, s, _state);
}


/*************************************************************************
Fast version of MinQPSetStartingPoint(), which doesn't check its arguments.
For internal use only.

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpsetstartingpointfast(minqpstate* state,
     /* Real    */ ae_vector* x,
     ae_state *) //_state)
{
    ae_int_t n;


    n = state->n;
    ae_v_move(&state->startx.ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->havex = ae_true;
}


/*************************************************************************
Fast version of MinQPSetOrigin(), which doesn't check its arguments.
For internal use only.

  -- ALGLIB --
     Copyright 11.01.2011 by Bochkanov Sergey
*************************************************************************/
void minqpsetoriginfast(minqpstate* state,
     /* Real    */ ae_vector* xorigin,
     ae_state *) //_state)
{
    ae_int_t n;


    n = state->n;
    ae_v_move(&state->xorigin.ptr.p_double[0], 1, &xorigin->ptr.p_double[0], 1, ae_v_len(0,n-1));
}


/*************************************************************************
Having current point XC, this function scans active constraints,   detects
ones  which  should  be deactivated by step in the antigradient direction,
deactivates them, and returns constraints count.

For example, if we have:
  XC=1.0,
  antigradient is equal to [-2]
  x>=0, x<=1
then this function will detect that constraint "x<=1" will be  deactivated
by  attempt  to  make  step  in the antigradient direction. Only presently
active constraints are counted and deactivated.

INPUT PARAMETERS:
    State   -   MinQP state.
    XC      -   current point, must be feasible with respect to all
                constrants
    Buf     -   temporary buffer, automatically resized if needed

OUTPUT PARAMETERS:
    State   -   this function changes following fields of State:
                * State.ActiveB     
                * State.ActiveLin   
                * State.ActiveLinCnt
                * State.ActiveCM
                * State.ActiveCR
    
RESULT:

Number of constraints being deactivated.

NOTES:
    this function uses fields of State which are marked as
    "DropConstraints temporaries". These fields are changed during algorithm
    operation.

  -- ALGLIB --
     Copyright 29.02.2012 by Bochkanov Sergey
*************************************************************************/
static ae_int_t minqp_minqpdropconstraints(minqpstate* state,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* buf,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t ncols;
    ae_int_t i;
    ae_int_t j;
    ae_int_t cidx;
    double v;
    ae_bool reactivationflag;
    double norm0;
    double norm1;
    ae_int_t idx;
    ae_int_t nactivebnd;
    ae_int_t orthogonalbasissize;
    ae_int_t result;


    n = state->n;
    
    /*
     * Handle important special case - no active linear constraints,
     * only boundary constraints are present
     */
    if( state->activelincnt==0 )
    {
        
        /*
         * Calculate gradient vector
         */
        rvectorsetlengthatleast(buf, n, _state);
        cqmadx(&state->a, xc, buf, _state);
        ae_v_add(&buf->ptr.p_double[0], 1, &state->b.ptr.p_double[0], 1, ae_v_len(0,n-1));
        result = 0;
        for(i=0; i<=n-1; i++)
        {
            if( state->activeb.ptr.p_bool[i] )
            {
                if( (state->havebndl.ptr.p_bool[i]&&state->havebndu.ptr.p_bool[i])&&ae_fp_eq(state->workbndl.ptr.p_double[i],state->workbndu.ptr.p_double[i]) )
                {
                    continue;
                }
                if( (state->havebndl.ptr.p_bool[i]&&ae_fp_eq(xc->ptr.p_double[i],state->workbndl.ptr.p_double[i]))&&ae_fp_less(buf->ptr.p_double[i],0) )
                {
                    result = result+1;
                    state->activeb.ptr.p_bool[i] = ae_false;
                    continue;
                }
                if( (state->havebndu.ptr.p_bool[i]&&ae_fp_eq(xc->ptr.p_double[i],state->workbndu.ptr.p_double[i]))&&ae_fp_greater(buf->ptr.p_double[i],0) )
                {
                    result = result+1;
                    state->activeb.ptr.p_bool[i] = ae_false;
                    continue;
                }
            }
        }
        return result;
    }
    
    /*
     * Now we have following problem: given ActiveLinCnt>0 active linear
     * constraints and NActiveBnd>=0 active boundary constraints, we want
     * to determine which constraints should be deactivated.
     *
     * We use method of Lagrange multipliers to do so:
     * * active constraints are written as C*x=r, where c[i] is I-th row
     *   of C, and C contains all active constraints - boundary and 
     *   general linear ones.
     * * gradient at x=XC is denoted as g
     * * at x=xc we have following linear system: SUM(lc[i]*C[i])=-g,
     *   where lc[i] is I-th Lagrange multiplier, C[i] is an I-th row of C.
     * * solution of this system is a set of Lagrange coefficients. Inequality
     *   constraints with non-positive Lagrange coefficients are deactivated.
     *   Equality constraints (either boundary constraints of the form a<=x[i]<=a
     *   or general linear ones a'*x=b) are left unchanged. 
     * * after initial deactivation of constraints we perform re-check, because
     *   sometimes algorithm deactivates too many constraints.
     *   We repeatedly calculate descent vector and check its feasibility. In
     *   case step in descent direction activates some constraints, these
     *   constraints are activated and we perform check one more time, until we 
     *   finally get descent direction which does not activate anything.
     *
     * Lagrange system we want to solve has N*(NActiveBnd+ActiveLinCnt) size;
     * it can be pretty large when we have many boundary constraints activated.
     * Straightforward attemp to solve such system gives us O(N^3) worst-case
     * performance.
     *
     * However, we can exploit structure of the matrix to solve this system
     * efficiently. Consider an example problem with N=6, 3 active boundary
     * constraints and 2 active linear constraints. After permutation of rows
     * (ones which correspond to active boundary constraints are moved to the
     * beginning) matrix of the Lagrange system will have following structure:
     *
     *   ( lc0 )T   ( d0                      )      ( -g0 )T
     *   ( lc1 )    (     d1                  )      ( -g1 )
     *   ( lc2 )  * (         d2              )      ( -g2 )
     *   ( lc3 )    ( s00 s01 s02 s03 s04 s05 ) *  = ( -g3 )
     *   ( lc4 )    ( s10 s11 s12 s13 s14 s15 )      ( -g4 )
     *                                               ( -g5 )
     *
     * Furthermore, diagonal elements di are either -1 (in case lower bound
     * is active) or +1 (in case upper bound is active). Such system can be
     * easility orthogonalized by Gram–Schmidt process applied to rows:
     *
     *     lc*C = -g
     *     lc*inv(T)*T*C = -g, where T is transformation matrix
     *     (lc*inv(T))*Q = -g, where Q is orthogonal
     *     lc = -g*Q'*T
     *
     * Special structure of the system matrix C allows us to easily
     * orthogonalize it in O(N*ActiveLinCnt^2) operations. It is important
     * to note that (NActiveBnd+ActiveLinCnt)*(NActiveBnd+ActiveLinCnt) transform
     * matrix T has special structure too:
     *
     *   ( d0                      )
     *   (     d1                  )
     *   (         d2              )
     *   ( t00 t01 t02 t03 t04 t05 )
     *   ( t10 t11 t12 t13 t14 t15 )
     *
     * It allow us to store it using N*ActiveLinCnt storage.
     *
     * Thus, our Lagrangian solver has O(N*ActiveLinCnt^2) time complexity and
     * O(N*ActiveLinCnt) memory complexity, which is significant improvement over
     * naive O(N^3)/O(N^2) approach.
     */
    nactivebnd = 0;
    for(j=0; j<=n-1; j++)
    {
        if( state->activeb.ptr.p_bool[j] )
        {
            nactivebnd = nactivebnd+1;
        }
    }
    ncols = state->activelincnt+nactivebnd;
    
    /*
     * Generate temporary system and right part - ones
     * before permutation of rows is applied:
     * * LagrangeSystemTmp - ActiveLinCnt*N block of the system matrix (right part)
     * * LSRightPartTmp - right part, antigradient
     */
    rvectorsetlengthatleast(&state->lsrightparttmp, n, _state);
    rmatrixsetlengthatleast(&state->lagrangesystemtmp, state->activelincnt, n, _state);
    cqmadx(&state->a, xc, &state->lsrightparttmp, _state);
    ae_v_add(&state->lsrightparttmp.ptr.p_double[0], 1, &state->b.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_muld(&state->lsrightparttmp.ptr.p_double[0], 1, ae_v_len(0,n-1), -1.0);
    cidx = 0;
    for(j=0; j<=state->nec+state->nic-1; j++)
    {
        if( state->activelin.ptr.p_bool[j] )
        {
            ae_v_move(&state->lagrangesystemtmp.ptr.pp_double[cidx][0], 1, &state->workcleic.ptr.pp_double[j][0], 1, ae_v_len(0,n-1));
            cidx = cidx+1;
        }
    }
    cidx = 0;
    for(j=0; j<=n-1; j++)
    {
        if( state->activeb.ptr.p_bool[j] )
        {
            ae_assert((state->havebndl.ptr.p_bool[j]&&ae_fp_eq(xc->ptr.p_double[j],state->workbndl.ptr.p_double[j]))||(state->havebndu.ptr.p_bool[j]&&ae_fp_eq(xc->ptr.p_double[j],state->workbndu.ptr.p_double[j])), "MinQPDropConstraints: internal error (1)", _state);
            if( ae_fp_eq(xc->ptr.p_double[j],state->workbndl.ptr.p_double[j]) )
            {
                
                /*
                 * We normalize our system to have all d[i]=1.0
                 * by multiplying columns corresponding to negative d[i] by -1.
                 */
                ae_v_muld(&state->lagrangesystemtmp.ptr.pp_double[0][j], state->lagrangesystemtmp.stride, ae_v_len(0,state->activelincnt-1), -1);
                state->lsrightparttmp.ptr.p_double[j] = -state->lsrightparttmp.ptr.p_double[j];
            }
            cidx = cidx+1;
        }
    }
    
    /*
     * Initialize right part of the transform matrix (ActiveLinCnt*NCols block)
     * by corresponding part of the identity matrix:
     *
     * LSTransform :=  (          1   )
     *                 (            1 )
     *
     * Initial matrix T is NCols*NCols identity matrix, but top NActiveBnd
     * columns are not changed by orthogonalization process, so we store
     * only its bottom block, which is actually changed during orthogonalization.
     */
    rmatrixsetlengthatleast(&state->lstransform, state->activelincnt, ncols, _state);
    for(i=0; i<=ncols-1; i++)
    {
        for(j=0; j<=state->activelincnt-1; j++)
        {
            state->lstransform.ptr.pp_double[j][i] = 0;
        }
    }
    for(i=0; i<=state->activelincnt-1; i++)
    {
        state->lstransform.ptr.pp_double[i][nactivebnd+i] = 1.0;
    }
    
    /*
     * Reoder columns of the temporary system; columns corresponding to active
     * boundary constraints are moved to the beginning
     */
    rmatrixsetlengthatleast(&state->lagrangesystem, state->activelincnt, n, _state);
    rvectorsetlengthatleast(&state->lsrightpart, n, _state);
    idx = 0;
    for(i=0; i<=n-1; i++)
    {
        if( state->activeb.ptr.p_bool[i] )
        {
            ae_v_move(&state->lagrangesystem.ptr.pp_double[0][idx], state->lagrangesystem.stride, &state->lagrangesystemtmp.ptr.pp_double[0][i], state->lagrangesystemtmp.stride, ae_v_len(0,state->activelincnt-1));
            state->lsrightpart.ptr.p_double[idx] = state->lsrightparttmp.ptr.p_double[i];
            idx = idx+1;
        }
    }
    for(i=0; i<=n-1; i++)
    {
        if( !state->activeb.ptr.p_bool[i] )
        {
            ae_v_move(&state->lagrangesystem.ptr.pp_double[0][idx], state->lagrangesystem.stride, &state->lagrangesystemtmp.ptr.pp_double[0][i], state->lagrangesystemtmp.stride, ae_v_len(0,state->activelincnt-1));
            state->lsrightpart.ptr.p_double[idx] = state->lsrightparttmp.ptr.p_double[i];
            idx = idx+1;
        }
    }
    
    /*
     * Perform orthogonalization of rows of C:
     * * first NActiveBnd rows are already orthogonalized, we just
     *   have to orthogonalize next ActiveLinCnt with respect to each other.
     * * during orthogonalization we monitor size of the current basis
     *   (stored in OrthogonalBasisSize variable). When basis size will become
     *   equal to the dimensionality of the space, we will stop adding new
     *   vectors to basis.
     */
    orthogonalbasissize = nactivebnd;
    for(i=0; i<=state->activelincnt-1; i++)
    {
        
        /*
         * Calculate norm of I-th rows. If it is zero, skip rows.
         */
        v = ae_v_dotproduct(&state->lagrangesystem.ptr.pp_double[i][0], 1, &state->lagrangesystem.ptr.pp_double[i][0], 1, ae_v_len(0,n-1));
        norm0 = ae_sqrt(v, _state);
        if( ae_fp_eq(norm0,0) )
        {
            continue;
        }
        
        /*
         * Orthogonalize with respect to first NActiveBnd rows (correspond
         * to active boundary constraints).
         *
         * Transformation which is applied to LagrangeSystem is also applied
         * to LSTransform, which stores transformation matrix.
         */
        for(j=0; j<=nactivebnd-1; j++)
        {
            state->lstransform.ptr.pp_double[i][j] = -state->lagrangesystem.ptr.pp_double[i][j];
            state->lagrangesystem.ptr.pp_double[i][j] = 0.0;
        }
        
        /*
         * Orthogonalize with respect to active linear constraints
         *
         * NOTE: first NActiveBnd elements of I-th row are already
         *       zero; we perform orthogonalization only when NActiveBnd<N.
         */
        if( nactivebnd<n )
        {
            for(j=0; j<=i-1; j++)
            {
                v = ae_v_dotproduct(&state->lagrangesystem.ptr.pp_double[j][nactivebnd], 1, &state->lagrangesystem.ptr.pp_double[i][nactivebnd], 1, ae_v_len(nactivebnd,n-1));
                ae_v_subd(&state->lagrangesystem.ptr.pp_double[i][nactivebnd], 1, &state->lagrangesystem.ptr.pp_double[j][nactivebnd], 1, ae_v_len(nactivebnd,n-1), v);
                ae_v_subd(&state->lstransform.ptr.pp_double[i][0], 1, &state->lstransform.ptr.pp_double[j][0], 1, ae_v_len(0,j+nactivebnd), v);
            }
        }
        
        /*
         * Determine row norm and decay factor. In case row is almost zero
         * when compared with its norm before orthogonalization, it is considered
         * zero.
         *
         * We perform normalization of the row and apply similar transformation
         * to LSTransform.
         *
         * NOTE: row can be considered zero in one more case - size of the
         *       current basis is already N (dimension of space).
         */
        v = ae_v_dotproduct(&state->lagrangesystem.ptr.pp_double[i][0], 1, &state->lagrangesystem.ptr.pp_double[i][0], 1, ae_v_len(0,n-1));
        norm1 = ae_sqrt(v, _state);
        if( ae_fp_greater_eq(norm1,10*n*ae_machineepsilon*norm0)&&orthogonalbasissize<n )
        {
            v = 1/norm1;
            orthogonalbasissize = orthogonalbasissize+1;
        }
        else
        {
            v = 0.0;
        }
        ae_v_muld(&state->lagrangesystem.ptr.pp_double[i][0], 1, ae_v_len(0,n-1), v);
        ae_v_muld(&state->lstransform.ptr.pp_double[i][0], 1, ae_v_len(0,ncols-1), v);
    }
    
    /*
     * After factorization of the Lagrange system, we have to solve it:
     * * first, we perform multiplication by Q
     * * then, we multiply by T
     */
    rvectorsetlengthatleast(buf, ncols, _state);
    rvectorsetlengthatleast(&state->lagrangecoeffs, ncols, _state);
    for(i=0; i<=nactivebnd-1; i++)
    {
        buf->ptr.p_double[i] = state->lsrightpart.ptr.p_double[i];
    }
    for(i=0; i<=state->activelincnt-1; i++)
    {
        if( nactivebnd<n )
        {
            v = ae_v_dotproduct(&state->lagrangesystem.ptr.pp_double[i][nactivebnd], 1, &state->lsrightpart.ptr.p_double[nactivebnd], 1, ae_v_len(nactivebnd,n-1));
        }
        else
        {
            v = 0;
        }
        buf->ptr.p_double[nactivebnd+i] = v;
    }
    for(i=0; i<=nactivebnd-1; i++)
    {
        v = ae_v_dotproduct(&state->lstransform.ptr.pp_double[0][i], state->lstransform.stride, &buf->ptr.p_double[nactivebnd], 1, ae_v_len(0,state->activelincnt-1));
        state->lagrangecoeffs.ptr.p_double[i] = buf->ptr.p_double[i]+v;
    }
    for(i=0; i<=state->activelincnt-1; i++)
    {
        v = ae_v_dotproduct(&state->lstransform.ptr.pp_double[i][nactivebnd+i], state->lstransform.stride, &buf->ptr.p_double[nactivebnd+i], 1, ae_v_len(i,state->activelincnt-1));
        state->lagrangecoeffs.ptr.p_double[nactivebnd+i] = v;
    }
    
    /*
     * We've solved Lagrange system.
     * Store previous values of ActiveB/ActiveLin
     * Deactivate bound and inequality constraints corresponding to negative Lagrange multipliers.
     * Active equality constraints are left unchanged independently of the multiplier sign.
     */
    bvectorsetlengthatleast(&state->prevactiveb, n, _state);
    bvectorsetlengthatleast(&state->prevactivelin, state->nec+state->nic, _state);
    for(i=0; i<=n-1; i++)
    {
        state->prevactiveb.ptr.p_bool[i] = state->activeb.ptr.p_bool[i];
    }
    for(i=0; i<=state->nec+state->nic-1; i++)
    {
        state->prevactivelin.ptr.p_bool[i] = state->activelin.ptr.p_bool[i];
    }
    cidx = 0;
    for(i=0; i<=n-1; i++)
    {
        if( state->activeb.ptr.p_bool[i] )
        {
            if( (state->havebndl.ptr.p_bool[i]&&state->havebndu.ptr.p_bool[i])&&ae_fp_eq(state->workbndl.ptr.p_double[i],state->workbndu.ptr.p_double[i]) )
            {
                
                /*
                 * Equality constraint a[i]=x[i]=b[i], skip it
                 */
                cidx = cidx+1;
                continue;
            }
            if( ae_fp_less(state->lagrangecoeffs.ptr.p_double[cidx],0) )
            {
                state->activeb.ptr.p_bool[i] = ae_false;
            }
            cidx = cidx+1;
        }
    }
    for(i=0; i<=state->nec+state->nic-1; i++)
    {
        if( state->activelin.ptr.p_bool[i] )
        {
            if( i<state->nec )
            {
                
                /*
                 * Equality constraint, skip it
                 */
                cidx = cidx+1;
                continue;
            }
            if( ae_fp_less(state->lagrangecoeffs.ptr.p_double[cidx],0) )
            {
                state->activelin.ptr.p_bool[i] = ae_false;
            }
            cidx = cidx+1;
        }
    }
    
    /*
     * Post-process constraints, maybe some of them need reactivation.
     * Each iteration of the cycle refreshes ActiveCM and other related structures.
     */
    for(;;)
    {
        minqp_orthogonalizeactiveconstraints(xc, &state->activeb, n, &state->workcleic, &state->activelin, state->nec+state->nic, &state->activecm, &state->activecr, &state->activelincnt, _state);
        minqp_minqpgradientandconstraineddescent(state, xc, buf, &state->lstmp0, _state);
        reactivationflag = ae_false;
        for(i=0; i<=n-1; i++)
        {
            if( !state->activeb.ptr.p_bool[i]&&state->prevactiveb.ptr.p_bool[i] )
            {
                if( ((state->havebndl.ptr.p_bool[i]&&ae_fp_eq(xc->ptr.p_double[i],state->workbndl.ptr.p_double[i]))&&ae_fp_less(state->lstmp0.ptr.p_double[i],0))||((state->havebndu.ptr.p_bool[i]&&ae_fp_eq(xc->ptr.p_double[i],state->workbndu.ptr.p_double[i]))&&ae_fp_greater(state->lstmp0.ptr.p_double[i],0)) )
                {
                    reactivationflag = ae_true;
                    state->activeb.ptr.p_bool[i] = ae_true;
                }
            }
        }
        for(i=state->nec; i<=state->nec+state->nic-1; i++)
        {
            if( !state->activelin.ptr.p_bool[i]&&state->prevactivelin.ptr.p_bool[i] )
            {
                v = ae_v_dotproduct(&state->workcleic.ptr.pp_double[i][0], 1, &state->lstmp0.ptr.p_double[0], 1, ae_v_len(0,n-1));
                if( ae_fp_greater(v,0) )
                {
                    reactivationflag = ae_true;
                    state->activelin.ptr.p_bool[i] = ae_true;
                }
            }
        }
        if( !reactivationflag )
        {
            break;
        }
    }
    
    /*
     * Number of deactivated constraints
     */
    result = 0;
    for(i=0; i<=n-1; i++)
    {
        if( state->prevactiveb.ptr.p_bool[i]&&!state->activeb.ptr.p_bool[i] )
        {
            result = result+1;
        }
    }
    for(i=0; i<=state->nec+state->nic-1; i++)
    {
        if( state->prevactivelin.ptr.p_bool[i]&&!state->activelin.ptr.p_bool[i] )
        {
            result = result+1;
        }
    }
    return result;
}


/*************************************************************************
Having feasible current point XC and possibly infeasible candidate   point
XN,  this  function  performs  longest  step  from  XC to XN which retains
feasibility. In case XN is found to be infeasible, at least one constraint
is activated.

For example, if we have:
  XC=0.5
  XN=1.2
  x>=0, x<=1
then this function will move us to X=1.0 and activate constraint "x<=1".

INPUT PARAMETERS:
    State   -   MinQP state.
    XC      -   current point, must be feasible with respect to
                all constraints
    XN      -   candidate point, can be infeasible with respect to some
                constraints
    Buf     -   temporary buffer, automatically resized if needed

OUTPUT PARAMETERS:
    State   -   this function changes following fields of State:
                * State.ActiveB     -   active boundary constraints
                * State.ActiveC     -   active linear constraints
    XC      -   new position

RESULT:
    number of constraints being activated by step.  Constraints which were
    DEACTIVATED are ignored (do not influence function value).

  -- ALGLIB --
     Copyright 29.02.2012 by Bochkanov Sergey
*************************************************************************/
static ae_int_t minqp_minqpboundedstepandactivation(minqpstate* state,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* xn,
     /* Real    */ ae_vector* buf,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t i;
    ae_int_t varidx;
    ae_int_t linidx;
    double varval;
    double stplen;
    double prevmax;
    double vc;
    double vn;
    ae_int_t result;


    n = state->n;
    rvectorsetlengthatleast(buf, n, _state);
    ae_v_move(&buf->ptr.p_double[0], 1, &xn->ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_sub(&buf->ptr.p_double[0], 1, &xc->ptr.p_double[0], 1, ae_v_len(0,n-1));
    varidx = -1;
    varval = 0.0;
    linidx = -1;
    stplen = ae_maxrealnumber;
    for(i=0; i<=n-1; i++)
    {
        if( !state->activeb.ptr.p_bool[i] )
        {
            ae_assert(!state->havebndl.ptr.p_bool[i]||ae_fp_greater_eq(xc->ptr.p_double[i],state->workbndl.ptr.p_double[i]), "MinQPBoundedStepAndActivation: internal error - infeasible XC", _state);
            ae_assert(!state->havebndu.ptr.p_bool[i]||ae_fp_less_eq(xc->ptr.p_double[i],state->workbndu.ptr.p_double[i]), "MinQPBoundedStepAndActivation: internal error - infeasible XC", _state);
            if( state->havebndl.ptr.p_bool[i]&&ae_fp_less_eq(xn->ptr.p_double[i],state->workbndl.ptr.p_double[i]) )
            {
                prevmax = stplen;
                stplen = safeminposrv(xc->ptr.p_double[i]-state->workbndl.ptr.p_double[i], xc->ptr.p_double[i]-xn->ptr.p_double[i], stplen, _state);
                if( ae_fp_less(stplen,prevmax) )
                {
                    varidx = i;
                    linidx = -1;
                    varval = state->workbndl.ptr.p_double[i];
                }
            }
            if( state->havebndu.ptr.p_bool[i]&&ae_fp_greater_eq(xn->ptr.p_double[i],state->workbndu.ptr.p_double[i]) )
            {
                prevmax = stplen;
                stplen = safeminposrv(state->workbndu.ptr.p_double[i]-xc->ptr.p_double[i], xn->ptr.p_double[i]-xc->ptr.p_double[i], stplen, _state);
                if( ae_fp_less(stplen,prevmax) )
                {
                    varidx = i;
                    linidx = -1;
                    varval = state->workbndu.ptr.p_double[i];
                }
            }
        }
    }
    for(i=state->nec; i<=state->nec+state->nic-1; i++)
    {
        if( !state->activelin.ptr.p_bool[i] )
        {
            vc = ae_v_dotproduct(&state->workcleic.ptr.pp_double[i][0], 1, &xc->ptr.p_double[0], 1, ae_v_len(0,n-1));
            vc = vc-state->workcleic.ptr.pp_double[i][n];
            vn = ae_v_dotproduct(&state->workcleic.ptr.pp_double[i][0], 1, &xn->ptr.p_double[0], 1, ae_v_len(0,n-1));
            vn = vn-state->workcleic.ptr.pp_double[i][n];
            if( ae_fp_greater_eq(vn,0) )
            {
                
                /*
                 * XN is infeasible. Depending on feasibility of XC,
                 * we can perform either non-zero or zero step.
                 */
                if( ae_fp_less(vc,0) )
                {
                    
                    /*
                     * XC is strictly feasible with respect to I-th constraint,
                     * we can perform non-zero step because there is non-zero distance
                     * between XC and bound.
                     */
                    prevmax = stplen;
                    stplen = safeminposrv(-vc, vn-vc, stplen, _state);
                    if( ae_fp_less(stplen,prevmax) )
                    {
                        linidx = i;
                        varidx = -1;
                    }
                }
                else
                {
                    
                    /*
                     * XC is at the boundary (or slightly beyond it), and XN is
                     * strictly infeasible.
                     * The only thing we can do is to perform zero step and activate
                     * I-th constraint.
                     */
                    stplen = 0;
                    linidx = i;
                    varidx = -1;
                }
            }
        }
    }
    stplen = ae_minreal(stplen, 1.0, _state);
    
    /*
     * Perform step
     */
    ae_v_move(&buf->ptr.p_double[0], 1, &xn->ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_sub(&buf->ptr.p_double[0], 1, &xc->ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_muld(&buf->ptr.p_double[0], 1, ae_v_len(0,n-1), stplen);
    ae_v_add(&buf->ptr.p_double[0], 1, &xc->ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_assert(varidx<0||linidx<0, "MinQPBoundedStepAndActivation: internal error (2)", _state);
    if( ae_fp_less(stplen,1.0) )
    {
        if( varidx>=0 )
        {
            state->activeb.ptr.p_bool[varidx] = ae_true;
        }
        if( linidx>=0 )
        {
            state->activelin.ptr.p_bool[linidx] = ae_true;
        }
    }
    result = postprocessboundedstep(buf, xc, &state->workbndl, &state->havebndl, &state->workbndu, &state->havebndu, n, 0, varidx, varval, stplen, stplen, _state);
    if( linidx>=0 )
    {
        result = result+1;
    }
    if( result>0 )
    {
        minqp_orthogonalizeactiveconstraints(buf, &state->activeb, n, &state->workcleic, &state->activelin, state->nec+state->nic, &state->activecm, &state->activecr, &state->activelincnt, _state);
    }
    ae_v_move(&xc->ptr.p_double[0], 1, &buf->ptr.p_double[0], 1, ae_v_len(0,n-1));
    return result;
}


/*************************************************************************
Model value: f = 0.5*x'*A*x + b'*x

INPUT PARAMETERS:
    A       -   convex quadratic model; only main quadratic term is used,
                other parts of the model (D/Q/linear term) are ignored.
                This function does not modify model state.
    B       -   right part
    XC      -   evaluation point
    Tmp     -   temporary buffer, automatically resized if needed

  -- ALGLIB --
     Copyright 20.06.2012 by Bochkanov Sergey
*************************************************************************/
static double minqp_minqpmodelvalue(convexquadraticmodel* a,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* xc,
     ae_int_t n,
     /* Real    */ ae_vector* tmp,
     ae_state *_state)
{
    double v0;
    double v1;
    double result;


    rvectorsetlengthatleast(tmp, n, _state);
    cqmadx(a, xc, tmp, _state);
    v0 = ae_v_dotproduct(&xc->ptr.p_double[0], 1, &tmp->ptr.p_double[0], 1, ae_v_len(0,n-1));
    v1 = ae_v_dotproduct(&xc->ptr.p_double[0], 1, &b->ptr.p_double[0], 1, ae_v_len(0,n-1));
    result = 0.5*v0+v1;
    return result;
}


/*************************************************************************
Constrained descent direction from the feasible point:
* without constraints, equal to the antigradient
* in the presence of active constraints it is equal to the projection of
  antigradient into feasible set

INPUT PARAMETERS:
    State   -   MinQP state.
    XC      -   current point, must be feasible with respect to
                all constraints
    GC      -   possibly preallocated buffer
    D       -   possibly preallocated buffer

OUTPUT PARAMETERS:
    GC      -   gradient of the original (unconstrained model without
                rank-K term) at the current point
    D       -   descent direction chosen with respect to constraints

RESULT:
    length of the direction vector

  -- ALGLIB --
     Copyright 20.06.2012 by Bochkanov Sergey
*************************************************************************/
static double minqp_minqpgradientandconstraineddescent(minqpstate* state,
     /* Real    */ ae_vector* xc,
     /* Real    */ ae_vector* gc,
     /* Real    */ ae_vector* d,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t n;
    double v;
    double norm0;
    double norm1;
    double result;


    n = state->n;
    rvectorsetlengthatleast(gc, n, _state);
    rvectorsetlengthatleast(d, n, _state);
    
    /*
     * Calculate gradient of the original model
     */
    cqmadx(&state->a, xc, gc, _state);
    ae_v_add(&gc->ptr.p_double[0], 1, &state->b.ptr.p_double[0], 1, ae_v_len(0,n-1));
    
    /*
     * Calculate descent direction
     */
    norm0 = 0;
    for(i=0; i<=n-1; i++)
    {
        if( state->activeb.ptr.p_bool[i] )
        {
            d->ptr.p_double[i] = 0;
        }
        else
        {
            d->ptr.p_double[i] = -gc->ptr.p_double[i];
            norm0 = norm0+ae_sqr(d->ptr.p_double[i], _state);
        }
    }
    norm0 = ae_sqrt(norm0, _state);
    for(i=0; i<=state->activelincnt-1; i++)
    {
        v = ae_v_dotproduct(&state->activecm.ptr.pp_double[i][0], 1, &d->ptr.p_double[0], 1, ae_v_len(0,n-1));
        ae_v_subd(&d->ptr.p_double[0], 1, &state->activecm.ptr.pp_double[i][0], 1, ae_v_len(0,n-1), v);
    }
    norm1 = 0;
    for(i=0; i<=n-1; i++)
    {
        norm1 = norm1+ae_sqr(d->ptr.p_double[i], _state);
    }
    norm1 = ae_sqrt(norm1, _state);
    if( ae_fp_less_eq(norm1,1.0E3*ae_machineepsilon*norm0) )
    {
        for(i=0; i<=n-1; i++)
        {
            d->ptr.p_double[i] = 0;
        }
    }
    result = 0.0;
    for(i=0; i<=n-1; i++)
    {
        result = result+ae_sqr(d->ptr.p_double[i], _state);
    }
    result = ae_sqrt(result, _state);
    return result;
}


/*************************************************************************
Optimum of A subject to:
a) active boundary constraints (given by ActiveB[] and corresponding
   elements of XC)
b) active linear constraints (given by C, R, LagrangeC)

INPUT PARAMETERS:
    A       -   main quadratic term of the model;
                although structure may  store  linear  and  rank-K  terms,
                these terms are ignored and rewritten  by  this  function.
    ANorm   -   estimate of ||A|| (2-norm is used)
    B       -   array[N], linear term of the model
    ActiveB -   array[N], active boundary constraints
    XC      -   array[N], constrained values
    C       -   equality constraints, array[K,N]
    R       -   array[K], right part of the equality constraints
    K       -   number of equality constraints
    XN      -   possibly preallocated buffer
    Tmp     -   temporary buffer (automatically resized)
    Tmp1    -   temporary buffer (automatically resized)

OUTPUT PARAMETERS:
    A       -   modified quadratic model (this function changes rank-K
                term and linear term of the model)
    LagrangeC-  current estimate of the Lagrange coefficients
    XN      -   solution

RESULT:
    True on success, False on failure (non-SPD model)

  -- ALGLIB --
     Copyright 20.06.2012 by Bochkanov Sergey
*************************************************************************/
static ae_bool minqp_minqpconstrainedoptimum(convexquadraticmodel* a,
     double anorm,
     /* Real    */ ae_vector* b,
     /* Boolean */ ae_vector* activeb,
     /* Real    */ ae_vector* xc,
     ae_int_t n,
     /* Real    */ ae_matrix* c,
     /* Real    */ ae_vector* r,
     ae_int_t k,
     /* Real    */ ae_vector* xn,
     /* Real    */ ae_vector* tmp,
     /* Real    */ ae_vector* lagrangec,
     ae_state *_state)
{
    ae_int_t itidx;
    ae_int_t i;
    double v;
    double feaserrold;
    double feaserrnew;
    double theta;
    ae_bool result;


    
    /*
     * Prepare model
     */
    theta = 100.0*anorm;
    rvectorsetlengthatleast(tmp, n, _state);
    rvectorsetlengthatleast(lagrangec, k, _state);
    cqmsetactiveset(a, xc, activeb, _state);
    cqmsetq(a, c, r, k, theta, _state);
    
    /*
     * Iterate until optimal values of Lagrange multipliers are found
     */
    for(i=0; i<=k-1; i++)
    {
        lagrangec->ptr.p_double[i] = 0;
    }
    feaserrnew = ae_maxrealnumber;
    result = ae_true;
    for(itidx=1; itidx<=minqp_maxlagrangeits; itidx++)
    {
        
        /*
         * Generate right part B using linear term and current
         * estimate of the Lagrange multipliers.
         */
        ae_v_move(&tmp->ptr.p_double[0], 1, &b->ptr.p_double[0], 1, ae_v_len(0,n-1));
        for(i=0; i<=k-1; i++)
        {
            v = lagrangec->ptr.p_double[i];
            ae_v_subd(&tmp->ptr.p_double[0], 1, &c->ptr.pp_double[i][0], 1, ae_v_len(0,n-1), v);
        }
        cqmsetb(a, tmp, _state);
        
        /*
         * Solve
         */
        result = cqmconstrainedoptimum(a, xn, _state);
        if( !result )
        {
            return result;
        }
        
        /*
         * Compare feasibility errors.
         * Terminate if error decreased too slowly.
         */
        feaserrold = feaserrnew;
        feaserrnew = 0;
        for(i=0; i<=k-1; i++)
        {
            v = ae_v_dotproduct(&c->ptr.pp_double[i][0], 1, &xn->ptr.p_double[0], 1, ae_v_len(0,n-1));
            feaserrnew = feaserrnew+ae_sqr(v-r->ptr.p_double[i], _state);
        }
        feaserrnew = ae_sqrt(feaserrnew, _state);
        if( ae_fp_greater_eq(feaserrnew,0.2*feaserrold) )
        {
            break;
        }
        
        /*
         * Update Lagrange multipliers
         */
        for(i=0; i<=k-1; i++)
        {
            v = ae_v_dotproduct(&c->ptr.pp_double[i][0], 1, &xn->ptr.p_double[0], 1, ae_v_len(0,n-1));
            lagrangec->ptr.p_double[i] = lagrangec->ptr.p_double[i]-theta*(v-r->ptr.p_double[i]);
        }
    }
    return result;
}


/*************************************************************************
This function builds orthogonal constraint matrix using:
* active set given by ActiveB/ActiveLin
* XC, which stores current point (and - values of active boundary constraints)
* original matrix of constraints CLEIC

As result, it fills:
* ActiveCM - by orthogonalized constraints
* ActiveCR - by right parts

INPUT PARAMETERS:
    XC      -   current point, array[N]; used to determine values of active
                boundary constraints.
    ActiveB -   active boundary constraints, array[N]
    N       -   size
    CLEIC   -   matrix of ALL linear constraints, array[CTotal,N+1].
                first N elements of each row store constraint, last
                element contains right part.
    ActiveLin-  active linear constraints
    CTotal  -   number of linear constraints (equality + inequality)
    ActiveCM-   possibly preallocated buffer
    ActiveCR-   possibly preallocated buffer

OUTPUT PARAMETERS:
    ActiveCM-   matrix of active linear constraints, array[...,N].
                Rows of this matrix are:
                * orthogonal to the active boundary constraints (i.e.
                  have zeros at corresponding positions)
                * orthogonal to each other
                Right parts are stored in separate array.
    ActiveCR-   matrix of active right parts
    NActive -   number of active linear constraints

  -- ALGLIB --
     Copyright 20.06.2012 by Bochkanov Sergey
*************************************************************************/
static void minqp_orthogonalizeactiveconstraints(/* Real    */ ae_vector* xc,
     /* Boolean */ ae_vector* activeb,
     ae_int_t n,
     /* Real    */ ae_matrix* cleic,
     /* Boolean */ ae_vector* activelin,
     ae_int_t ctotal,
     /* Real    */ ae_matrix* activecm,
     /* Real    */ ae_vector* activecr,
     ae_int_t* nactive,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t ridx;
    double v;
    double norm0;
    double norm1;


    
    /*
     * Determine number of active constraints
     */
    *nactive = 0;
    for(i=0; i<=ctotal-1; i++)
    {
        if( activelin->ptr.p_bool[i] )
        {
            *nactive = *nactive+1;
        }
    }
    if( *nactive==0 )
    {
        return;
    }
    
    /*
     * Orthogonalize active linear constraints:
     * * with respect to active bound constraints
     * * with respect to each other
     */
    rmatrixsetlengthatleast(activecm, *nactive, n, _state);
    rvectorsetlengthatleast(activecr, *nactive, _state);
    ridx = 0;
    for(i=0; i<=ctotal-1; i++)
    {
        if( activelin->ptr.p_bool[i] )
        {
            ae_v_move(&activecm->ptr.pp_double[ridx][0], 1, &cleic->ptr.pp_double[i][0], 1, ae_v_len(0,n-1));
            activecr->ptr.p_double[ridx] = cleic->ptr.pp_double[i][n];
            ridx = ridx+1;
        }
    }
    for(i=0; i<=*nactive-1; i++)
    {
        
        /*
         * Calculate row norm before transformation.
         * Zero rows are skipped
         */
        norm0 = ae_v_dotproduct(&activecm->ptr.pp_double[i][0], 1, &activecm->ptr.pp_double[i][0], 1, ae_v_len(0,n-1));
        norm0 = ae_sqrt(norm0, _state);
        if( ae_fp_eq(norm0,0) )
        {
            continue;
        }
        
        /*
         * Orthogonalize I-th row with respect to boundary constraints
         */
        for(j=0; j<=n-1; j++)
        {
            if( activeb->ptr.p_bool[j] )
            {
                activecr->ptr.p_double[i] = activecr->ptr.p_double[i]-activecm->ptr.pp_double[i][j]*xc->ptr.p_double[j];
                activecm->ptr.pp_double[i][j] = 0.0;
            }
        }
        
        /*
         * Orthogonalize I-th row with respect to previous linear constraints
         */
        for(k=0; k<=i-1; k++)
        {
            v = ae_v_dotproduct(&activecm->ptr.pp_double[k][0], 1, &activecm->ptr.pp_double[i][0], 1, ae_v_len(0,n-1));
            ae_v_subd(&activecm->ptr.pp_double[i][0], 1, &activecm->ptr.pp_double[k][0], 1, ae_v_len(0,n-1), v);
            activecr->ptr.p_double[i] = activecr->ptr.p_double[i]-v*activecr->ptr.p_double[k];
        }
        
        /*
         * Calculate row norm after transformation.
         * Rows which become almost zero (when compared with their
         * magnitude before orthogonalization) are artificially
         * zeroed.
         */
        norm1 = ae_v_dotproduct(&activecm->ptr.pp_double[i][0], 1, &activecm->ptr.pp_double[i][0], 1, ae_v_len(0,n-1));
        norm1 = ae_sqrt(norm1, _state);
        if( ae_fp_greater(norm1,1.0E4*ae_machineepsilon*norm0) )
        {
            v = 1/norm1;
            ae_v_muld(&activecm->ptr.pp_double[i][0], 1, ae_v_len(0,n-1), v);
            activecr->ptr.p_double[i] = activecr->ptr.p_double[i]*v;
        }
        else
        {
            for(j=0; j<=n-1; j++)
            {
                activecm->ptr.pp_double[i][j] = 0;
            }
            activecr->ptr.p_double[i] = 0;
        }
    }
}


ae_bool _minqpstate_init(minqpstate* p, ae_state *_state, ae_bool make_automatic)
{
    if( !_convexquadraticmodel_init(&p->a, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->b, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->bndl, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->bndu, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->havebndl, 0, DT_BOOL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->havebndu, 0, DT_BOOL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xorigin, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->startx, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->cmatrix, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->cr, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->ct, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->cleic, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xc, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->gc, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xn, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->pg, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->activelin, 0, DT_BOOL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->activeb, 0, DT_BOOL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->activecm, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->activecr, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->workbndl, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->workbndu, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->workcleic, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tmp0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tmp1, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tmpfeas, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->tmpm0, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->lstransform, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->lagrangesystem, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->lagrangesystemtmp, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->prevactiveb, 0, DT_BOOL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->prevactivelin, 0, DT_BOOL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->lsrightpart, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->lsrightparttmp, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->lagrangecoeffs, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->lstmp0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->lstmp1, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !_normestimatorstate_init(&p->estimator, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _minqpstate_init_copy(minqpstate* dst, minqpstate* src, ae_state *_state, ae_bool make_automatic)
{
    dst->n = src->n;
    dst->algokind = src->algokind;
    if( !_convexquadraticmodel_init_copy(&dst->a, &src->a, _state, make_automatic) )
        return ae_false;
    dst->anorm = src->anorm;
    if( !ae_vector_init_copy(&dst->b, &src->b, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->bndl, &src->bndl, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->bndu, &src->bndu, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->havebndl, &src->havebndl, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->havebndu, &src->havebndu, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->xorigin, &src->xorigin, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->startx, &src->startx, _state, make_automatic) )
        return ae_false;
    dst->havex = src->havex;
    if( !ae_matrix_init_copy(&dst->cmatrix, &src->cmatrix, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->cr, &src->cr, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->ct, &src->ct, _state, make_automatic) )
        return ae_false;
    dst->ccount = src->ccount;
    if( !ae_matrix_init_copy(&dst->cleic, &src->cleic, _state, make_automatic) )
        return ae_false;
    dst->nec = src->nec;
    dst->nic = src->nic;
    if( !ae_vector_init_copy(&dst->xc, &src->xc, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->gc, &src->gc, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->xn, &src->xn, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->pg, &src->pg, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->activelin, &src->activelin, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->activeb, &src->activeb, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->activecm, &src->activecm, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->activecr, &src->activecr, _state, make_automatic) )
        return ae_false;
    dst->activelincnt = src->activelincnt;
    if( !ae_vector_init_copy(&dst->workbndl, &src->workbndl, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->workbndu, &src->workbndu, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->workcleic, &src->workcleic, _state, make_automatic) )
        return ae_false;
    dst->repinneriterationscount = src->repinneriterationscount;
    dst->repouteriterationscount = src->repouteriterationscount;
    dst->repncholesky = src->repncholesky;
    dst->repnmv = src->repnmv;
    dst->repterminationtype = src->repterminationtype;
    if( !ae_vector_init_copy(&dst->tmp0, &src->tmp0, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tmp1, &src->tmp1, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tmpfeas, &src->tmpfeas, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->tmpm0, &src->tmpm0, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->lstransform, &src->lstransform, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->lagrangesystem, &src->lagrangesystem, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->lagrangesystemtmp, &src->lagrangesystemtmp, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->prevactiveb, &src->prevactiveb, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->prevactivelin, &src->prevactivelin, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->lsrightpart, &src->lsrightpart, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->lsrightparttmp, &src->lsrightparttmp, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->lagrangecoeffs, &src->lagrangecoeffs, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->lstmp0, &src->lstmp0, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->lstmp1, &src->lstmp1, _state, make_automatic) )
        return ae_false;
    if( !_normestimatorstate_init_copy(&dst->estimator, &src->estimator, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _minqpstate_clear(minqpstate* p)
{
    _convexquadraticmodel_clear(&p->a);
    ae_vector_clear(&p->b);
    ae_vector_clear(&p->bndl);
    ae_vector_clear(&p->bndu);
    ae_vector_clear(&p->havebndl);
    ae_vector_clear(&p->havebndu);
    ae_vector_clear(&p->xorigin);
    ae_vector_clear(&p->startx);
    ae_matrix_clear(&p->cmatrix);
    ae_vector_clear(&p->cr);
    ae_vector_clear(&p->ct);
    ae_matrix_clear(&p->cleic);
    ae_vector_clear(&p->xc);
    ae_vector_clear(&p->gc);
    ae_vector_clear(&p->xn);
    ae_vector_clear(&p->pg);
    ae_vector_clear(&p->activelin);
    ae_vector_clear(&p->activeb);
    ae_matrix_clear(&p->activecm);
    ae_vector_clear(&p->activecr);
    ae_vector_clear(&p->workbndl);
    ae_vector_clear(&p->workbndu);
    ae_matrix_clear(&p->workcleic);
    ae_vector_clear(&p->tmp0);
    ae_vector_clear(&p->tmp1);
    ae_vector_clear(&p->tmpfeas);
    ae_matrix_clear(&p->tmpm0);
    ae_matrix_clear(&p->lstransform);
    ae_matrix_clear(&p->lagrangesystem);
    ae_matrix_clear(&p->lagrangesystemtmp);
    ae_vector_clear(&p->prevactiveb);
    ae_vector_clear(&p->prevactivelin);
    ae_vector_clear(&p->lsrightpart);
    ae_vector_clear(&p->lsrightparttmp);
    ae_vector_clear(&p->lagrangecoeffs);
    ae_vector_clear(&p->lstmp0);
    ae_vector_clear(&p->lstmp1);
    _normestimatorstate_clear(&p->estimator);
}


ae_bool _minqpreport_init(minqpreport* /*p*/, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    return ae_true;
}


ae_bool _minqpreport_init_copy(minqpreport* dst, minqpreport* src, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    dst->inneriterationscount = src->inneriterationscount;
    dst->outeriterationscount = src->outeriterationscount;
    dst->nmv = src->nmv;
    dst->ncholesky = src->ncholesky;
    dst->terminationtype = src->terminationtype;
    return ae_true;
}


void _minqpreport_clear(minqpreport* ) //p)
{
}




/*************************************************************************
                IMPROVED LEVENBERG-MARQUARDT METHOD FOR
                 NON-LINEAR LEAST SQUARES OPTIMIZATION

DESCRIPTION:
This function is used to find minimum of function which is represented  as
sum of squares:
    F(x) = f[0]^2(x[0],...,x[n-1]) + ... + f[m-1]^2(x[0],...,x[n-1])
using value of function vector f[] and Jacobian of f[].


REQUIREMENTS:
This algorithm will request following information during its operation:

* function vector f[] at given point X
* function vector f[] and Jacobian of f[] (simultaneously) at given point

There are several overloaded versions of  MinLMOptimize()  function  which
correspond  to  different LM-like optimization algorithms provided by this
unit. You should choose version which accepts fvec()  and jac() callbacks.
First  one  is used to calculate f[] at given point, second one calculates
f[] and Jacobian df[i]/dx[j].

You can try to initialize MinLMState structure with VJ  function and  then
use incorrect version  of  MinLMOptimize()  (for  example,  version  which
works  with  general  form function and does not provide Jacobian), but it
will  lead  to  exception  being  thrown  after first attempt to calculate
Jacobian.


USAGE:
1. User initializes algorithm state with MinLMCreateVJ() call
2. User tunes solver parameters with MinLMSetCond(),  MinLMSetStpMax() and
   other functions
3. User calls MinLMOptimize() function which  takes algorithm  state   and
   callback functions.
4. User calls MinLMResults() to get solution
5. Optionally, user may call MinLMRestartFrom() to solve  another  problem
   with same N/M but another starting point and/or another function.
   MinLMRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   dimension, N>1
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    M       -   number of functions f[i]
    X       -   initial solution, array[0..N-1]

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTES:
1. you may tune stopping conditions with MinLMSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLMSetStpMax() function to bound algorithm's steps.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatevj(ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_vector* x,
     minlmstate* state,
     ae_state *_state)
{

    _minlmstate_clear(state);

    ae_assert(n>=1, "MinLMCreateVJ: N<1!", _state);
    ae_assert(m>=1, "MinLMCreateVJ: M<1!", _state);
    ae_assert(x->cnt>=n, "MinLMCreateVJ: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "MinLMCreateVJ: X contains infinite or NaN values!", _state);
    
    /*
     * initialize, check parameters
     */
    state->teststep = 0;
    state->n = n;
    state->m = m;
    state->algomode = 1;
    state->hasf = ae_false;
    state->hasfi = ae_true;
    state->hasg = ae_false;
    
    /*
     * second stage of initialization
     */
    minlm_lmprepare(n, m, ae_false, state, _state);
    minlmsetacctype(state, 0, _state);
    minlmsetcond(state, 0, 0, 0, 0, _state);
    minlmsetxrep(state, ae_false, _state);
    minlmsetstpmax(state, 0, _state);
    minlmrestartfrom(state, x, _state);
}


/*************************************************************************
                IMPROVED LEVENBERG-MARQUARDT METHOD FOR
                 NON-LINEAR LEAST SQUARES OPTIMIZATION

DESCRIPTION:
This function is used to find minimum of function which is represented  as
sum of squares:
    F(x) = f[0]^2(x[0],...,x[n-1]) + ... + f[m-1]^2(x[0],...,x[n-1])
using value of function vector f[] only. Finite differences  are  used  to
calculate Jacobian.


REQUIREMENTS:
This algorithm will request following information during its operation:
* function vector f[] at given point X

There are several overloaded versions of  MinLMOptimize()  function  which
correspond  to  different LM-like optimization algorithms provided by this
unit. You should choose version which accepts fvec() callback.

You can try to initialize MinLMState structure with VJ  function and  then
use incorrect version  of  MinLMOptimize()  (for  example,  version  which
works with general form function and does not accept function vector), but
it will  lead  to  exception being thrown after first attempt to calculate
Jacobian.


USAGE:
1. User initializes algorithm state with MinLMCreateV() call
2. User tunes solver parameters with MinLMSetCond(),  MinLMSetStpMax() and
   other functions
3. User calls MinLMOptimize() function which  takes algorithm  state   and
   callback functions.
4. User calls MinLMResults() to get solution
5. Optionally, user may call MinLMRestartFrom() to solve  another  problem
   with same N/M but another starting point and/or another function.
   MinLMRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   dimension, N>1
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    M       -   number of functions f[i]
    X       -   initial solution, array[0..N-1]
    DiffStep-   differentiation step, >0

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

See also MinLMIteration, MinLMResults.

NOTES:
1. you may tune stopping conditions with MinLMSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLMSetStpMax() function to bound algorithm's steps.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatev(ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_vector* x,
     double diffstep,
     minlmstate* state,
     ae_state *_state)
{

    _minlmstate_clear(state);

    ae_assert(ae_isfinite(diffstep, _state), "MinLMCreateV: DiffStep is not finite!", _state);
    ae_assert(ae_fp_greater(diffstep,0), "MinLMCreateV: DiffStep<=0!", _state);
    ae_assert(n>=1, "MinLMCreateV: N<1!", _state);
    ae_assert(m>=1, "MinLMCreateV: M<1!", _state);
    ae_assert(x->cnt>=n, "MinLMCreateV: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "MinLMCreateV: X contains infinite or NaN values!", _state);
    
    /*
     * Initialize
     */
    state->teststep = 0;
    state->n = n;
    state->m = m;
    state->algomode = 0;
    state->hasf = ae_false;
    state->hasfi = ae_true;
    state->hasg = ae_false;
    state->diffstep = diffstep;
    
    /*
     * Second stage of initialization
     */
    minlm_lmprepare(n, m, ae_false, state, _state);
    minlmsetacctype(state, 1, _state);
    minlmsetcond(state, 0, 0, 0, 0, _state);
    minlmsetxrep(state, ae_false, _state);
    minlmsetstpmax(state, 0, _state);
    minlmrestartfrom(state, x, _state);
}


/*************************************************************************
    LEVENBERG-MARQUARDT-LIKE METHOD FOR NON-LINEAR OPTIMIZATION

DESCRIPTION:
This  function  is  used  to  find  minimum  of general form (not "sum-of-
-squares") function
    F = F(x[0], ..., x[n-1])
using  its  gradient  and  Hessian.  Levenberg-Marquardt modification with
L-BFGS pre-optimization and internal pre-conditioned  L-BFGS  optimization
after each Levenberg-Marquardt step is used.


REQUIREMENTS:
This algorithm will request following information during its operation:

* function value F at given point X
* F and gradient G (simultaneously) at given point X
* F, G and Hessian H (simultaneously) at given point X

There are several overloaded versions of  MinLMOptimize()  function  which
correspond  to  different LM-like optimization algorithms provided by this
unit. You should choose version which accepts func(),  grad()  and  hess()
function pointers. First pointer is used to calculate F  at  given  point,
second  one  calculates  F(x)  and  grad F(x),  third one calculates F(x),
grad F(x), hess F(x).

You can try to initialize MinLMState structure with FGH-function and  then
use incorrect version of MinLMOptimize() (for example, version which  does
not provide Hessian matrix), but it will lead to  exception  being  thrown
after first attempt to calculate Hessian.


USAGE:
1. User initializes algorithm state with MinLMCreateFGH() call
2. User tunes solver parameters with MinLMSetCond(),  MinLMSetStpMax() and
   other functions
3. User calls MinLMOptimize() function which  takes algorithm  state   and
   pointers (delegates, etc.) to callback functions.
4. User calls MinLMResults() to get solution
5. Optionally, user may call MinLMRestartFrom() to solve  another  problem
   with same N but another starting point and/or another function.
   MinLMRestartFrom() allows to reuse already initialized structure.


INPUT PARAMETERS:
    N       -   dimension, N>1
                * if given, only leading N elements of X are used
                * if not given, automatically determined from size of X
    X       -   initial solution, array[0..N-1]

OUTPUT PARAMETERS:
    State   -   structure which stores algorithm state

NOTES:
1. you may tune stopping conditions with MinLMSetCond() function
2. if target function contains exp() or other fast growing functions,  and
   optimization algorithm makes too large steps which leads  to  overflow,
   use MinLMSetStpMax() function to bound algorithm's steps.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatefgh(ae_int_t n,
     /* Real    */ ae_vector* x,
     minlmstate* state,
     ae_state *_state)
{

    _minlmstate_clear(state);

    ae_assert(n>=1, "MinLMCreateFGH: N<1!", _state);
    ae_assert(x->cnt>=n, "MinLMCreateFGH: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "MinLMCreateFGH: X contains infinite or NaN values!", _state);
    
    /*
     * initialize
     */
    state->teststep = 0;
    state->n = n;
    state->m = 0;
    state->algomode = 2;
    state->hasf = ae_true;
    state->hasfi = ae_false;
    state->hasg = ae_true;
    
    /*
     * init2
     */
    minlm_lmprepare(n, 0, ae_true, state, _state);
    minlmsetacctype(state, 2, _state);
    minlmsetcond(state, 0, 0, 0, 0, _state);
    minlmsetxrep(state, ae_false, _state);
    minlmsetstpmax(state, 0, _state);
    minlmrestartfrom(state, x, _state);
}


/*************************************************************************
This function sets stopping conditions for Levenberg-Marquardt optimization
algorithm.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    EpsG    -   >=0
                The  subroutine  finishes  its  work   if   the  condition
                |v|<EpsG is satisfied, where:
                * |.| means Euclidian norm
                * v - scaled gradient vector, v[i]=g[i]*s[i]
                * g - gradient
                * s - scaling coefficients set by MinLMSetScale()
    EpsF    -   >=0
                The  subroutine  finishes  its work if on k+1-th iteration
                the  condition  |F(k+1)-F(k)|<=EpsF*max{|F(k)|,|F(k+1)|,1}
                is satisfied.
    EpsX    -   >=0
                The subroutine finishes its work if  on  k+1-th  iteration
                the condition |v|<=EpsX is fulfilled, where:
                * |.| means Euclidian norm
                * v - scaled step vector, v[i]=dx[i]/s[i]
                * dx - ste pvector, dx=X(k+1)-X(k)
                * s - scaling coefficients set by MinLMSetScale()
    MaxIts  -   maximum number of iterations. If MaxIts=0, the  number  of
                iterations   is    unlimited.   Only   Levenberg-Marquardt
                iterations  are  counted  (L-BFGS/CG  iterations  are  NOT
                counted because their cost is very low compared to that of
                LM).

Passing EpsG=0, EpsF=0, EpsX=0 and MaxIts=0 (simultaneously) will lead to
automatic stopping criterion selection (small EpsX).

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlmsetcond(minlmstate* state,
     double epsg,
     double epsf,
     double epsx,
     ae_int_t maxits,
     ae_state *_state)
{


    ae_assert(ae_isfinite(epsg, _state), "MinLMSetCond: EpsG is not finite number!", _state);
    ae_assert(ae_fp_greater_eq(epsg,0), "MinLMSetCond: negative EpsG!", _state);
    ae_assert(ae_isfinite(epsf, _state), "MinLMSetCond: EpsF is not finite number!", _state);
    ae_assert(ae_fp_greater_eq(epsf,0), "MinLMSetCond: negative EpsF!", _state);
    ae_assert(ae_isfinite(epsx, _state), "MinLMSetCond: EpsX is not finite number!", _state);
    ae_assert(ae_fp_greater_eq(epsx,0), "MinLMSetCond: negative EpsX!", _state);
    ae_assert(maxits>=0, "MinLMSetCond: negative MaxIts!", _state);
    if( ((ae_fp_eq(epsg,0)&&ae_fp_eq(epsf,0))&&ae_fp_eq(epsx,0))&&maxits==0 )
    {
        epsx = 1.0E-6;
    }
    state->epsg = epsg;
    state->epsf = epsf;
    state->epsx = epsx;
    state->maxits = maxits;
}


/*************************************************************************
This function turns on/off reporting.

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    NeedXRep-   whether iteration reports are needed or not

If NeedXRep is True, algorithm will call rep() callback function if  it is
provided to MinLMOptimize(). Both Levenberg-Marquardt and internal  L-BFGS
iterations are reported.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlmsetxrep(minlmstate* state, ae_bool needxrep, ae_state *) //_state)
{


    state->xrep = needxrep;
}


/*************************************************************************
This function sets maximum step length

INPUT PARAMETERS:
    State   -   structure which stores algorithm state
    StpMax  -   maximum step length, >=0. Set StpMax to 0.0,  if you don't
                want to limit step length.

Use this subroutine when you optimize target function which contains exp()
or  other  fast  growing  functions,  and optimization algorithm makes too
large  steps  which  leads  to overflow. This function allows us to reject
steps  that  are  too  large  (and  therefore  expose  us  to the possible
overflow) without actually calculating function value at the x+stp*d.

NOTE: non-zero StpMax leads to moderate  performance  degradation  because
intermediate  step  of  preconditioned L-BFGS optimization is incompatible
with limits on step size.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minlmsetstpmax(minlmstate* state, double stpmax, ae_state *_state)
{


    ae_assert(ae_isfinite(stpmax, _state), "MinLMSetStpMax: StpMax is not finite!", _state);
    ae_assert(ae_fp_greater_eq(stpmax,0), "MinLMSetStpMax: StpMax<0!", _state);
    state->stpmax = stpmax;
}


/*************************************************************************
This function sets scaling coefficients for LM optimizer.

ALGLIB optimizers use scaling matrices to test stopping  conditions  (step
size and gradient are scaled before comparison with tolerances).  Scale of
the I-th variable is a translation invariant measure of:
a) "how large" the variable is
b) how large the step should be to make significant changes in the function

Generally, scale is NOT considered to be a form of preconditioner.  But LM
optimizer is unique in that it uses scaling matrix both  in  the  stopping
condition tests and as Marquardt damping factor.

Proper scaling is very important for the algorithm performance. It is less
important for the quality of results, but still has some influence (it  is
easier  to  converge  when  variables  are  properly  scaled, so premature
stopping is possible when very badly scalled variables are  combined  with
relaxed stopping conditions).

INPUT PARAMETERS:
    State   -   structure stores algorithm state
    S       -   array[N], non-zero scaling coefficients
                S[i] may be negative, sign doesn't matter.

  -- ALGLIB --
     Copyright 14.01.2011 by Bochkanov Sergey
*************************************************************************/
void minlmsetscale(minlmstate* state,
     /* Real    */ ae_vector* s,
     ae_state *_state)
{
    ae_int_t i;


    ae_assert(s->cnt>=state->n, "MinLMSetScale: Length(S)<N", _state);
    for(i=0; i<=state->n-1; i++)
    {
        ae_assert(ae_isfinite(s->ptr.p_double[i], _state), "MinLMSetScale: S contains infinite or NAN elements", _state);
        ae_assert(ae_fp_neq(s->ptr.p_double[i],0), "MinLMSetScale: S contains zero elements", _state);
        state->s.ptr.p_double[i] = ae_fabs(s->ptr.p_double[i], _state);
    }
}


/*************************************************************************
This function sets boundary constraints for LM optimizer

Boundary constraints are inactive by default (after initial creation).
They are preserved until explicitly turned off with another SetBC() call.

INPUT PARAMETERS:
    State   -   structure stores algorithm state
    BndL    -   lower bounds, array[N].
                If some (all) variables are unbounded, you may specify
                very small number or -INF (latter is recommended because
                it will allow solver to use better algorithm).
    BndU    -   upper bounds, array[N].
                If some (all) variables are unbounded, you may specify
                very large number or +INF (latter is recommended because
                it will allow solver to use better algorithm).

NOTE 1: it is possible to specify BndL[i]=BndU[i]. In this case I-th
variable will be "frozen" at X[i]=BndL[i]=BndU[i].

NOTE 2: this solver has following useful properties:
* bound constraints are always satisfied exactly
* function is evaluated only INSIDE area specified by bound constraints
  or at its boundary

  -- ALGLIB --
     Copyright 14.01.2011 by Bochkanov Sergey
*************************************************************************/
void minlmsetbc(minlmstate* state,
     /* Real    */ ae_vector* bndl,
     /* Real    */ ae_vector* bndu,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t n;


    n = state->n;
    ae_assert(bndl->cnt>=n, "MinLMSetBC: Length(BndL)<N", _state);
    ae_assert(bndu->cnt>=n, "MinLMSetBC: Length(BndU)<N", _state);
    for(i=0; i<=n-1; i++)
    {
        ae_assert(ae_isfinite(bndl->ptr.p_double[i], _state)||ae_isneginf(bndl->ptr.p_double[i], _state), "MinLMSetBC: BndL contains NAN or +INF", _state);
        ae_assert(ae_isfinite(bndu->ptr.p_double[i], _state)||ae_isposinf(bndu->ptr.p_double[i], _state), "MinLMSetBC: BndU contains NAN or -INF", _state);
        state->bndl.ptr.p_double[i] = bndl->ptr.p_double[i];
        state->havebndl.ptr.p_bool[i] = ae_isfinite(bndl->ptr.p_double[i], _state);
        state->bndu.ptr.p_double[i] = bndu->ptr.p_double[i];
        state->havebndu.ptr.p_bool[i] = ae_isfinite(bndu->ptr.p_double[i], _state);
    }
}


/*************************************************************************
This function is used to change acceleration settings

You can choose between three acceleration strategies:
* AccType=0, no acceleration.
* AccType=1, secant updates are used to update quadratic model after  each
  iteration. After fixed number of iterations (or after  model  breakdown)
  we  recalculate  quadratic  model  using  analytic  Jacobian  or  finite
  differences. Number of secant-based iterations depends  on  optimization
  settings: about 3 iterations - when we have analytic Jacobian, up to 2*N
  iterations - when we use finite differences to calculate Jacobian.

AccType=1 is recommended when Jacobian  calculation  cost  is  prohibitive
high (several Mx1 function vector calculations  followed  by  several  NxN
Cholesky factorizations are faster than calculation of one M*N  Jacobian).
It should also be used when we have no Jacobian, because finite difference
approximation takes too much time to compute.

Table below list  optimization  protocols  (XYZ  protocol  corresponds  to
MinLMCreateXYZ) and acceleration types they support (and use by  default).

ACCELERATION TYPES SUPPORTED BY OPTIMIZATION PROTOCOLS:

protocol    0   1   comment
V           +   +
VJ          +   +
FGH         +

DAFAULT VALUES:

protocol    0   1   comment
V               x   without acceleration it is so slooooooooow
VJ          x
FGH         x

NOTE: this  function should be called before optimization. Attempt to call
it during algorithm iterations may result in unexpected behavior.

NOTE: attempt to call this function with unsupported protocol/acceleration
combination will result in exception being thrown.

  -- ALGLIB --
     Copyright 14.10.2010 by Bochkanov Sergey
*************************************************************************/
void minlmsetacctype(minlmstate* state,
     ae_int_t acctype,
     ae_state *_state)
{


    ae_assert((acctype==0||acctype==1)||acctype==2, "MinLMSetAccType: incorrect AccType!", _state);
    if( acctype==2 )
    {
        acctype = 0;
    }
    if( acctype==0 )
    {
        state->maxmodelage = 0;
        state->makeadditers = ae_false;
        return;
    }
    if( acctype==1 )
    {
        ae_assert(state->hasfi, "MinLMSetAccType: AccType=1 is incompatible with current protocol!", _state);
        if( state->algomode==0 )
        {
            state->maxmodelage = 2*state->n;
        }
        else
        {
            state->maxmodelage = minlm_smallmodelage;
        }
        state->makeadditers = ae_false;
        return;
    }
}


/*************************************************************************
NOTES:

1. Depending on function used to create state  structure,  this  algorithm
   may accept Jacobian and/or Hessian and/or gradient.  According  to  the
   said above, there ase several versions of this function,  which  accept
   different sets of callbacks.

   This flexibility opens way to subtle errors - you may create state with
   MinLMCreateFGH() (optimization using Hessian), but call function  which
   does not accept Hessian. So when algorithm will request Hessian,  there
   will be no callback to call. In this case exception will be thrown.

   Be careful to avoid such errors because there is no way to find them at
   compile time - you can see them at runtime only.

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool minlmiteration(minlmstate* state, ae_state *_state)
{
    ae_int_t n;
    ae_int_t m;
    ae_bool bflag;
    ae_int_t iflag;
    double v;
    double s;
    double t;
    ae_int_t i;
    ae_int_t k;
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
        m = state->rstate.ia.ptr.p_int[1];
        iflag = state->rstate.ia.ptr.p_int[2];
        i = state->rstate.ia.ptr.p_int[3];
        k = state->rstate.ia.ptr.p_int[4];
        bflag = state->rstate.ba.ptr.p_bool[0];
        v = state->rstate.ra.ptr.p_double[0];
        s = state->rstate.ra.ptr.p_double[1];
        t = state->rstate.ra.ptr.p_double[2];
    }
    else
    {
        n = -983;
        m = -989;
        iflag = -834;
        i = 900;
        k = -287;
        bflag = ae_false;
        v = 214;
        s = -338;
        t = -686;
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
    if( state->rstate.stage==4 )
    {
        goto lbl_4;
    }
    if( state->rstate.stage==5 )
    {
        goto lbl_5;
    }
    if( state->rstate.stage==6 )
    {
        goto lbl_6;
    }
    if( state->rstate.stage==7 )
    {
        goto lbl_7;
    }
    if( state->rstate.stage==8 )
    {
        goto lbl_8;
    }
    if( state->rstate.stage==9 )
    {
        goto lbl_9;
    }
    if( state->rstate.stage==10 )
    {
        goto lbl_10;
    }
    if( state->rstate.stage==11 )
    {
        goto lbl_11;
    }
    if( state->rstate.stage==12 )
    {
        goto lbl_12;
    }
    if( state->rstate.stage==13 )
    {
        goto lbl_13;
    }
    if( state->rstate.stage==14 )
    {
        goto lbl_14;
    }
    if( state->rstate.stage==15 )
    {
        goto lbl_15;
    }
    if( state->rstate.stage==16 )
    {
        goto lbl_16;
    }
    if( state->rstate.stage==17 )
    {
        goto lbl_17;
    }
    if( state->rstate.stage==18 )
    {
        goto lbl_18;
    }
    
    /*
     * Routine body
     */
    
    /*
     * prepare
     */
    n = state->n;
    m = state->m;
    state->repiterationscount = 0;
    state->repterminationtype = 0;
    state->repfuncidx = -1;
    state->repvaridx = -1;
    state->repnfunc = 0;
    state->repnjac = 0;
    state->repngrad = 0;
    state->repnhess = 0;
    state->repncholesky = 0;
    
    /*
     * check consistency of constraints,
     * enforce feasibility of the solution
     * set constraints
     */
    if( !enforceboundaryconstraints(&state->xbase, &state->bndl, &state->havebndl, &state->bndu, &state->havebndu, n, 0, _state) )
    {
        state->repterminationtype = -3;
        result = ae_false;
        return result;
    }
    minqpsetbc(&state->qpstate, &state->bndl, &state->bndu, _state);
    
    /*
     *  Check, that transferred derivative value is right
     */
    minlm_clearrequestfields(state, _state);
    if( !(state->algomode==1&&ae_fp_greater(state->teststep,0)) )
    {
        goto lbl_19;
    }
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->needfij = ae_true;
    i = 0;
lbl_21:
    if( i>n-1 )
    {
        goto lbl_23;
    }
    ae_assert((state->havebndl.ptr.p_bool[i]&&ae_fp_less_eq(state->bndl.ptr.p_double[i],state->x.ptr.p_double[i]))||!state->havebndl.ptr.p_bool[i], "MinLM: internal error(State.X is out of bounds)", _state);
    ae_assert((state->havebndu.ptr.p_bool[i]&&ae_fp_less_eq(state->x.ptr.p_double[i],state->bndu.ptr.p_double[i]))||!state->havebndu.ptr.p_bool[i], "MinLMIteration: internal error(State.X is out of bounds)", _state);
    v = state->x.ptr.p_double[i];
    state->x.ptr.p_double[i] = v-state->teststep*state->s.ptr.p_double[i];
    if( state->havebndl.ptr.p_bool[i] )
    {
        state->x.ptr.p_double[i] = ae_maxreal(state->x.ptr.p_double[i], state->bndl.ptr.p_double[i], _state);
    }
    state->xm1 = state->x.ptr.p_double[i];
    state->rstate.stage = 0;
    goto lbl_rcomm;
lbl_0:
    ae_v_move(&state->fm1.ptr.p_double[0], 1, &state->fi.ptr.p_double[0], 1, ae_v_len(0,m-1));
    ae_v_move(&state->gm1.ptr.p_double[0], 1, &state->j.ptr.pp_double[0][i], state->j.stride, ae_v_len(0,m-1));
    state->x.ptr.p_double[i] = v+state->teststep*state->s.ptr.p_double[i];
    if( state->havebndu.ptr.p_bool[i] )
    {
        state->x.ptr.p_double[i] = ae_minreal(state->x.ptr.p_double[i], state->bndu.ptr.p_double[i], _state);
    }
    state->xp1 = state->x.ptr.p_double[i];
    state->rstate.stage = 1;
    goto lbl_rcomm;
lbl_1:
    ae_v_move(&state->fp1.ptr.p_double[0], 1, &state->fi.ptr.p_double[0], 1, ae_v_len(0,m-1));
    ae_v_move(&state->gp1.ptr.p_double[0], 1, &state->j.ptr.pp_double[0][i], state->j.stride, ae_v_len(0,m-1));
    state->x.ptr.p_double[i] = (state->xm1+state->xp1)/2;
    if( state->havebndl.ptr.p_bool[i] )
    {
        state->x.ptr.p_double[i] = ae_maxreal(state->x.ptr.p_double[i], state->bndl.ptr.p_double[i], _state);
    }
    if( state->havebndu.ptr.p_bool[i] )
    {
        state->x.ptr.p_double[i] = ae_minreal(state->x.ptr.p_double[i], state->bndu.ptr.p_double[i], _state);
    }
    state->rstate.stage = 2;
    goto lbl_rcomm;
lbl_2:
    ae_v_move(&state->fc1.ptr.p_double[0], 1, &state->fi.ptr.p_double[0], 1, ae_v_len(0,m-1));
    ae_v_move(&state->gc1.ptr.p_double[0], 1, &state->j.ptr.pp_double[0][i], state->j.stride, ae_v_len(0,m-1));
    state->x.ptr.p_double[i] = v;
    for(k=0; k<=m-1; k++)
    {
        if( !derivativecheck(state->fm1.ptr.p_double[k], state->gm1.ptr.p_double[k], state->fp1.ptr.p_double[k], state->gp1.ptr.p_double[k], state->fc1.ptr.p_double[k], state->gc1.ptr.p_double[k], state->xp1-state->xm1, _state) )
        {
            state->repfuncidx = k;
            state->repvaridx = i;
            state->repterminationtype = -7;
            result = ae_false;
            return result;
        }
    }
    i = i+1;
    goto lbl_21;
lbl_23:
    state->needfij = ae_false;
lbl_19:
    
    /*
     * Initial report of current point
     *
     * Note 1: we rewrite State.X twice because
     * user may accidentally change it after first call.
     *
     * Note 2: we set NeedF or NeedFI depending on what
     * information about function we have.
     */
    if( !state->xrep )
    {
        goto lbl_24;
    }
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minlm_clearrequestfields(state, _state);
    if( !state->hasf )
    {
        goto lbl_26;
    }
    state->needf = ae_true;
    state->rstate.stage = 3;
    goto lbl_rcomm;
lbl_3:
    state->needf = ae_false;
    goto lbl_27;
lbl_26:
    ae_assert(state->hasfi, "MinLM: internal error 2!", _state);
    state->needfi = ae_true;
    state->rstate.stage = 4;
    goto lbl_rcomm;
lbl_4:
    state->needfi = ae_false;
    v = ae_v_dotproduct(&state->fi.ptr.p_double[0], 1, &state->fi.ptr.p_double[0], 1, ae_v_len(0,m-1));
    state->f = v;
lbl_27:
    state->repnfunc = state->repnfunc+1;
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minlm_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 5;
    goto lbl_rcomm;
lbl_5:
    state->xupdated = ae_false;
lbl_24:
    
    /*
     * Prepare control variables
     */
    state->nu = 1;
    state->lambdav = -ae_maxrealnumber;
    state->modelage = state->maxmodelage+1;
    state->deltaxready = ae_false;
    state->deltafready = ae_false;
    
    /*
     * Main cycle.
     *
     * We move through it until either:
     * * one of the stopping conditions is met
     * * we decide that stopping conditions are too stringent
     *   and break from cycle
     *
     */
lbl_28:
    if( ae_false )
    {
        goto lbl_29;
    }
    
    /*
     * First, we have to prepare quadratic model for our function.
     * We use BFlag to ensure that model is prepared;
     * if it is false at the end of this block, something went wrong.
     *
     * We may either calculate brand new model or update old one.
     *
     * Before this block we have:
     * * State.XBase            - current position.
     * * State.DeltaX           - if DeltaXReady is True
     * * State.DeltaF           - if DeltaFReady is True
     *
     * After this block is over, we will have:
     * * State.XBase            - base point (unchanged)
     * * State.FBase            - F(XBase)
     * * State.GBase            - linear term
     * * State.QuadraticModel   - quadratic term
     * * State.LambdaV          - current estimate for lambda
     *
     * We also clear DeltaXReady/DeltaFReady flags
     * after initialization is done.
     */
    bflag = ae_false;
    if( !(state->algomode==0||state->algomode==1) )
    {
        goto lbl_30;
    }
    
    /*
     * Calculate f[] and Jacobian
     */
    if( !(state->modelage>state->maxmodelage||!(state->deltaxready&&state->deltafready)) )
    {
        goto lbl_32;
    }
    
    /*
     * Refresh model (using either finite differences or analytic Jacobian)
     */
    if( state->algomode!=0 )
    {
        goto lbl_34;
    }
    
    /*
     * Optimization using F values only.
     * Use finite differences to estimate Jacobian.
     */
    ae_assert(state->hasfi, "MinLMIteration: internal error when estimating Jacobian (no f[])", _state);
    k = 0;
lbl_36:
    if( k>n-1 )
    {
        goto lbl_38;
    }
    
    /*
     * We guard X[k] from leaving [BndL,BndU].
     * In case BndL=BndU, we assume that derivative in this direction is zero.
     */
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->x.ptr.p_double[k] = state->x.ptr.p_double[k]-state->s.ptr.p_double[k]*state->diffstep;
    if( state->havebndl.ptr.p_bool[k] )
    {
        state->x.ptr.p_double[k] = ae_maxreal(state->x.ptr.p_double[k], state->bndl.ptr.p_double[k], _state);
    }
    if( state->havebndu.ptr.p_bool[k] )
    {
        state->x.ptr.p_double[k] = ae_minreal(state->x.ptr.p_double[k], state->bndu.ptr.p_double[k], _state);
    }
    state->xm1 = state->x.ptr.p_double[k];
    minlm_clearrequestfields(state, _state);
    state->needfi = ae_true;
    state->rstate.stage = 6;
    goto lbl_rcomm;
lbl_6:
    state->repnfunc = state->repnfunc+1;
    ae_v_move(&state->fm1.ptr.p_double[0], 1, &state->fi.ptr.p_double[0], 1, ae_v_len(0,m-1));
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->x.ptr.p_double[k] = state->x.ptr.p_double[k]+state->s.ptr.p_double[k]*state->diffstep;
    if( state->havebndl.ptr.p_bool[k] )
    {
        state->x.ptr.p_double[k] = ae_maxreal(state->x.ptr.p_double[k], state->bndl.ptr.p_double[k], _state);
    }
    if( state->havebndu.ptr.p_bool[k] )
    {
        state->x.ptr.p_double[k] = ae_minreal(state->x.ptr.p_double[k], state->bndu.ptr.p_double[k], _state);
    }
    state->xp1 = state->x.ptr.p_double[k];
    minlm_clearrequestfields(state, _state);
    state->needfi = ae_true;
    state->rstate.stage = 7;
    goto lbl_rcomm;
lbl_7:
    state->repnfunc = state->repnfunc+1;
    ae_v_move(&state->fp1.ptr.p_double[0], 1, &state->fi.ptr.p_double[0], 1, ae_v_len(0,m-1));
    v = state->xp1-state->xm1;
    if( ae_fp_neq(v,0) )
    {
        v = 1/v;
        ae_v_moved(&state->j.ptr.pp_double[0][k], state->j.stride, &state->fp1.ptr.p_double[0], 1, ae_v_len(0,m-1), v);
        ae_v_subd(&state->j.ptr.pp_double[0][k], state->j.stride, &state->fm1.ptr.p_double[0], 1, ae_v_len(0,m-1), v);
    }
    else
    {
        for(i=0; i<=m-1; i++)
        {
            state->j.ptr.pp_double[i][k] = 0;
        }
    }
    k = k+1;
    goto lbl_36;
lbl_38:
    
    /*
     * Calculate F(XBase)
     */
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minlm_clearrequestfields(state, _state);
    state->needfi = ae_true;
    state->rstate.stage = 8;
    goto lbl_rcomm;
lbl_8:
    state->needfi = ae_false;
    state->repnfunc = state->repnfunc+1;
    state->repnjac = state->repnjac+1;
    
    /*
     * New model
     */
    state->modelage = 0;
    goto lbl_35;
lbl_34:
    
    /*
     * Obtain f[] and Jacobian
     */
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minlm_clearrequestfields(state, _state);
    state->needfij = ae_true;
    state->rstate.stage = 9;
    goto lbl_rcomm;
lbl_9:
    state->needfij = ae_false;
    state->repnfunc = state->repnfunc+1;
    state->repnjac = state->repnjac+1;
    
    /*
     * New model
     */
    state->modelage = 0;
lbl_35:
    goto lbl_33;
lbl_32:
    
    /*
     * State.J contains Jacobian or its current approximation;
     * refresh it using secant updates:
     *
     * f(x0+dx) = f(x0) + J*dx,
     * J_new = J_old + u*h'
     * h = x_new-x_old
     * u = (f_new - f_old - J_old*h)/(h'h)
     *
     * We can explicitly generate h and u, but it is
     * preferential to do in-place calculations. Only
     * I-th row of J_old is needed to calculate u[I],
     * so we can update J row by row in one pass.
     *
     * NOTE: we expect that State.XBase contains new point,
     * State.FBase contains old point, State.DeltaX and
     * State.DeltaY contain updates from last step.
     */
    ae_assert(state->deltaxready&&state->deltafready, "MinLMIteration: uninitialized DeltaX/DeltaF", _state);
    t = ae_v_dotproduct(&state->deltax.ptr.p_double[0], 1, &state->deltax.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_assert(ae_fp_neq(t,0), "MinLM: internal error (T=0)", _state);
    for(i=0; i<=m-1; i++)
    {
        v = ae_v_dotproduct(&state->j.ptr.pp_double[i][0], 1, &state->deltax.ptr.p_double[0], 1, ae_v_len(0,n-1));
        v = (state->deltaf.ptr.p_double[i]-v)/t;
        ae_v_addd(&state->j.ptr.pp_double[i][0], 1, &state->deltax.ptr.p_double[0], 1, ae_v_len(0,n-1), v);
    }
    ae_v_move(&state->fi.ptr.p_double[0], 1, &state->fibase.ptr.p_double[0], 1, ae_v_len(0,m-1));
    ae_v_add(&state->fi.ptr.p_double[0], 1, &state->deltaf.ptr.p_double[0], 1, ae_v_len(0,m-1));
    
    /*
     * Increase model age
     */
    state->modelage = state->modelage+1;
lbl_33:
    
    /*
     * Generate quadratic model:
     *     f(xbase+dx) =
     *       = (f0 + J*dx)'(f0 + J*dx)
     *       = f0^2 + dx'J'f0 + f0*J*dx + dx'J'J*dx
     *       = f0^2 + 2*f0*J*dx + dx'J'J*dx
     *
     * Note that we calculate 2*(J'J) instead of J'J because
     * our quadratic model is based on Tailor decomposition,
     * i.e. it has 0.5 before quadratic term.
     */
    rmatrixgemm(n, n, m, 2.0, &state->j, 0, 0, 1, &state->j, 0, 0, 0, 0.0, &state->quadraticmodel, 0, 0, _state);
    rmatrixmv(n, m, &state->j, 0, 0, 1, &state->fi, 0, &state->gbase, 0, _state);
    ae_v_muld(&state->gbase.ptr.p_double[0], 1, ae_v_len(0,n-1), 2);
    v = ae_v_dotproduct(&state->fi.ptr.p_double[0], 1, &state->fi.ptr.p_double[0], 1, ae_v_len(0,m-1));
    state->fbase = v;
    ae_v_move(&state->fibase.ptr.p_double[0], 1, &state->fi.ptr.p_double[0], 1, ae_v_len(0,m-1));
    
    /*
     * set control variables
     */
    bflag = ae_true;
lbl_30:
    if( state->algomode!=2 )
    {
        goto lbl_39;
    }
    ae_assert(!state->hasfi, "MinLMIteration: internal error (HasFI is True in Hessian-based mode)", _state);
    
    /*
     * Obtain F, G, H
     */
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minlm_clearrequestfields(state, _state);
    state->needfgh = ae_true;
    state->rstate.stage = 10;
    goto lbl_rcomm;
lbl_10:
    state->needfgh = ae_false;
    state->repnfunc = state->repnfunc+1;
    state->repngrad = state->repngrad+1;
    state->repnhess = state->repnhess+1;
    rmatrixcopy(n, n, &state->h, 0, 0, &state->quadraticmodel, 0, 0, _state);
    ae_v_move(&state->gbase.ptr.p_double[0], 1, &state->g.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->fbase = state->f;
    
    /*
     * set control variables
     */
    bflag = ae_true;
    state->modelage = 0;
lbl_39:
    ae_assert(bflag, "MinLM: internal integrity check failed!", _state);
    state->deltaxready = ae_false;
    state->deltafready = ae_false;
    
    /*
     * If Lambda is not initialized, initialize it using quadratic model
     */
    if( ae_fp_less(state->lambdav,0) )
    {
        state->lambdav = 0;
        for(i=0; i<=n-1; i++)
        {
            state->lambdav = ae_maxreal(state->lambdav, ae_fabs(state->quadraticmodel.ptr.pp_double[i][i], _state)*ae_sqr(state->s.ptr.p_double[i], _state), _state);
        }
        state->lambdav = 0.001*state->lambdav;
        if( ae_fp_eq(state->lambdav,0) )
        {
            state->lambdav = 1;
        }
    }
    
    /*
     * Test stopping conditions for function gradient
     */
    if( ae_fp_greater(minlm_boundedscaledantigradnorm(state, &state->xbase, &state->gbase, _state),state->epsg) )
    {
        goto lbl_41;
    }
    if( state->modelage!=0 )
    {
        goto lbl_43;
    }
    
    /*
     * Model is fresh, we can rely on it and terminate algorithm
     */
    state->repterminationtype = 4;
    if( !state->xrep )
    {
        goto lbl_45;
    }
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->f = state->fbase;
    minlm_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 11;
    goto lbl_rcomm;
lbl_11:
    state->xupdated = ae_false;
lbl_45:
    result = ae_false;
    return result;
    goto lbl_44;
lbl_43:
    
    /*
     * Model is not fresh, we should refresh it and test
     * conditions once more
     */
    state->modelage = state->maxmodelage+1;
    goto lbl_28;
lbl_44:
lbl_41:
    
    /*
     * Find value of Levenberg-Marquardt damping parameter which:
     * * leads to positive definite damped model
     * * within bounds specified by StpMax
     * * generates step which decreases function value
     *
     * After this block IFlag is set to:
     * * -3, if constraints are infeasible
     * * -2, if model update is needed (either Lambda growth is too large
     *       or step is too short, but we can't rely on model and stop iterations)
     * * -1, if model is fresh, Lambda have grown too large, termination is needed
     * *  0, if everything is OK, continue iterations
     *
     * State.Nu can have any value on enter, but after exit it is set to 1.0
     */
    iflag = -99;
lbl_47:
    if( ae_false )
    {
        goto lbl_48;
    }
    
    /*
     * Do we need model update?
     */
    if( state->modelage>0&&ae_fp_greater_eq(state->nu,minlm_suspiciousnu) )
    {
        iflag = -2;
        goto lbl_48;
    }
    
    /*
     * Setup quadratic solver and solve quadratic programming problem.
     * After problem is solved we'll try to bound step by StpMax
     * (Lambda will be increased if step size is too large).
     *
     * We use BFlag variable to indicate that we have to increase Lambda.
     * If it is False, we will try to increase Lambda and move to new iteration.
     */
    bflag = ae_true;
    minqpsetstartingpointfast(&state->qpstate, &state->xbase, _state);
    minqpsetoriginfast(&state->qpstate, &state->xbase, _state);
    minqpsetlineartermfast(&state->qpstate, &state->gbase, _state);
    minqpsetquadratictermfast(&state->qpstate, &state->quadraticmodel, ae_true, 0.0, _state);
    for(i=0; i<=n-1; i++)
    {
        state->tmp0.ptr.p_double[i] = state->quadraticmodel.ptr.pp_double[i][i]+state->lambdav/ae_sqr(state->s.ptr.p_double[i], _state);
    }
    minqprewritediagonal(&state->qpstate, &state->tmp0, _state);
    minqpoptimize(&state->qpstate, _state);
    minqpresultsbuf(&state->qpstate, &state->xdir, &state->qprep, _state);
    if( state->qprep.terminationtype>0 )
    {
        
        /*
         * successful solution of QP problem
         */
        ae_v_sub(&state->xdir.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
        v = ae_v_dotproduct(&state->xdir.ptr.p_double[0], 1, &state->xdir.ptr.p_double[0], 1, ae_v_len(0,n-1));
        if( ae_isfinite(v, _state) )
        {
            v = ae_sqrt(v, _state);
            if( ae_fp_greater(state->stpmax,0)&&ae_fp_greater(v,state->stpmax) )
            {
                bflag = ae_false;
            }
        }
        else
        {
            bflag = ae_false;
        }
    }
    else
    {
        
        /*
         * Either problem is non-convex (increase LambdaV) or constraints are inconsistent
         */
        ae_assert(state->qprep.terminationtype==-3||state->qprep.terminationtype==-5, "MinLM: unexpected completion code from QP solver", _state);
        if( state->qprep.terminationtype==-3 )
        {
            iflag = -3;
            goto lbl_48;
        }
        bflag = ae_false;
    }
    if( !bflag )
    {
        
        /*
         * Solution failed:
         * try to increase lambda to make matrix positive definite and continue.
         */
        if( !minlm_increaselambda(&state->lambdav, &state->nu, _state) )
        {
            iflag = -1;
            goto lbl_48;
        }
        goto lbl_47;
    }
    
    /*
     * Step in State.XDir and it is bounded by StpMax.
     *
     * We should check stopping conditions on step size here.
     * DeltaX, which is used for secant updates, is initialized here.
     *
     * This code is a bit tricky because sometimes XDir<>0, but
     * it is so small that XDir+XBase==XBase (in finite precision
     * arithmetics). So we set DeltaX to XBase, then
     * add XDir, and then subtract XBase to get exact value of
     * DeltaX.
     *
     * Step length is estimated using DeltaX.
     *
     * NOTE: stopping conditions are tested
     * for fresh models only (ModelAge=0)
     */
    ae_v_move(&state->deltax.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_add(&state->deltax.ptr.p_double[0], 1, &state->xdir.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_sub(&state->deltax.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->deltaxready = ae_true;
    v = 0.0;
    for(i=0; i<=n-1; i++)
    {
        v = v+ae_sqr(state->deltax.ptr.p_double[i]/state->s.ptr.p_double[i], _state);
    }
    v = ae_sqrt(v, _state);
    if( ae_fp_greater(v,state->epsx) )
    {
        goto lbl_49;
    }
    if( state->modelage!=0 )
    {
        goto lbl_51;
    }
    
    /*
     * Step is too short, model is fresh and we can rely on it.
     * Terminating.
     */
    state->repterminationtype = 2;
    if( !state->xrep )
    {
        goto lbl_53;
    }
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->f = state->fbase;
    minlm_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 12;
    goto lbl_rcomm;
lbl_12:
    state->xupdated = ae_false;
lbl_53:
    result = ae_false;
    return result;
    goto lbl_52;
lbl_51:
    
    /*
     * Step is suspiciously short, but model is not fresh
     * and we can't rely on it.
     */
    iflag = -2;
    goto lbl_48;
lbl_52:
lbl_49:
    
    /*
     * Let's evaluate new step:
     * a) if we have Fi vector, we evaluate it using rcomm, and
     *    then we manually calculate State.F as sum of squares of Fi[]
     * b) if we have F value, we just evaluate it through rcomm interface
     *
     * We prefer (a) because we may need Fi vector for additional
     * iterations
     */
    ae_assert(state->hasfi||state->hasf, "MinLM: internal error 2!", _state);
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_add(&state->x.ptr.p_double[0], 1, &state->xdir.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minlm_clearrequestfields(state, _state);
    if( !state->hasfi )
    {
        goto lbl_55;
    }
    state->needfi = ae_true;
    state->rstate.stage = 13;
    goto lbl_rcomm;
lbl_13:
    state->needfi = ae_false;
    v = ae_v_dotproduct(&state->fi.ptr.p_double[0], 1, &state->fi.ptr.p_double[0], 1, ae_v_len(0,m-1));
    state->f = v;
    ae_v_move(&state->deltaf.ptr.p_double[0], 1, &state->fi.ptr.p_double[0], 1, ae_v_len(0,m-1));
    ae_v_sub(&state->deltaf.ptr.p_double[0], 1, &state->fibase.ptr.p_double[0], 1, ae_v_len(0,m-1));
    state->deltafready = ae_true;
    goto lbl_56;
lbl_55:
    state->needf = ae_true;
    state->rstate.stage = 14;
    goto lbl_rcomm;
lbl_14:
    state->needf = ae_false;
lbl_56:
    state->repnfunc = state->repnfunc+1;
    if( ae_fp_greater_eq(state->f,state->fbase) )
    {
        
        /*
         * Increase lambda and continue
         */
        if( !minlm_increaselambda(&state->lambdav, &state->nu, _state) )
        {
            iflag = -1;
            goto lbl_48;
        }
        goto lbl_47;
    }
    
    /*
     * We've found our step!
     */
    iflag = 0;
    goto lbl_48;
    goto lbl_47;
lbl_48:
    state->nu = 1;
    ae_assert(iflag>=-3&&iflag<=0, "MinLM: internal integrity check failed!", _state);
    if( iflag==-3 )
    {
        state->repterminationtype = -3;
        result = ae_false;
        return result;
    }
    if( iflag==-2 )
    {
        state->modelage = state->maxmodelage+1;
        goto lbl_28;
    }
    if( iflag==-1 )
    {
        goto lbl_29;
    }
    
    /*
     * Levenberg-Marquardt step is ready.
     * Compare predicted vs. actual decrease and decide what to do with lambda.
     *
     * NOTE: we expect that State.DeltaX contains direction of step,
     * State.F contains function value at new point.
     */
    ae_assert(state->deltaxready, "MinLM: deltaX is not ready", _state);
    t = 0;
    for(i=0; i<=n-1; i++)
    {
        v = ae_v_dotproduct(&state->quadraticmodel.ptr.pp_double[i][0], 1, &state->deltax.ptr.p_double[0], 1, ae_v_len(0,n-1));
        t = t+state->deltax.ptr.p_double[i]*state->gbase.ptr.p_double[i]+0.5*state->deltax.ptr.p_double[i]*v;
    }
    state->predicteddecrease = -t;
    state->actualdecrease = -(state->f-state->fbase);
    if( ae_fp_less_eq(state->predicteddecrease,0) )
    {
        goto lbl_29;
    }
    v = state->actualdecrease/state->predicteddecrease;
    if( ae_fp_greater_eq(v,0.1) )
    {
        goto lbl_57;
    }
    if( minlm_increaselambda(&state->lambdav, &state->nu, _state) )
    {
        goto lbl_59;
    }
    
    /*
     * Lambda is too large, we have to break iterations.
     */
    state->repterminationtype = 7;
    if( !state->xrep )
    {
        goto lbl_61;
    }
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->f = state->fbase;
    minlm_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 15;
    goto lbl_rcomm;
lbl_15:
    state->xupdated = ae_false;
lbl_61:
    result = ae_false;
    return result;
lbl_59:
lbl_57:
    if( ae_fp_greater(v,0.5) )
    {
        minlm_decreaselambda(&state->lambdav, &state->nu, _state);
    }
    
    /*
     * Accept step, report it and
     * test stopping conditions on iterations count and function decrease.
     *
     * NOTE: we expect that State.DeltaX contains direction of step,
     * State.F contains function value at new point.
     *
     * NOTE2: we should update XBase ONLY. In the beginning of the next
     * iteration we expect that State.FIBase is NOT updated and
     * contains old value of a function vector.
     */
    ae_v_add(&state->xbase.ptr.p_double[0], 1, &state->deltax.ptr.p_double[0], 1, ae_v_len(0,n-1));
    if( !state->xrep )
    {
        goto lbl_63;
    }
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minlm_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 16;
    goto lbl_rcomm;
lbl_16:
    state->xupdated = ae_false;
lbl_63:
    state->repiterationscount = state->repiterationscount+1;
    if( state->repiterationscount>=state->maxits&&state->maxits>0 )
    {
        state->repterminationtype = 5;
    }
    if( state->modelage==0 )
    {
        if( ae_fp_less_eq(ae_fabs(state->f-state->fbase, _state),state->epsf*ae_maxreal(1, ae_maxreal(ae_fabs(state->f, _state), ae_fabs(state->fbase, _state), _state), _state)) )
        {
            state->repterminationtype = 1;
        }
    }
    if( state->repterminationtype<=0 )
    {
        goto lbl_65;
    }
    if( !state->xrep )
    {
        goto lbl_67;
    }
    
    /*
     * Report: XBase contains new point, F contains function value at new point
     */
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    minlm_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 17;
    goto lbl_rcomm;
lbl_17:
    state->xupdated = ae_false;
lbl_67:
    result = ae_false;
    return result;
lbl_65:
    state->modelage = state->modelage+1;
    goto lbl_28;
lbl_29:
    
    /*
     * Lambda is too large, we have to break iterations.
     */
    state->repterminationtype = 7;
    if( !state->xrep )
    {
        goto lbl_69;
    }
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xbase.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->f = state->fbase;
    minlm_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 18;
    goto lbl_rcomm;
lbl_18:
    state->xupdated = ae_false;
lbl_69:
    result = ae_false;
    return result;
    
    /*
     * Saving state
     */
lbl_rcomm:
    result = ae_true;
    state->rstate.ia.ptr.p_int[0] = n;
    state->rstate.ia.ptr.p_int[1] = m;
    state->rstate.ia.ptr.p_int[2] = iflag;
    state->rstate.ia.ptr.p_int[3] = i;
    state->rstate.ia.ptr.p_int[4] = k;
    state->rstate.ba.ptr.p_bool[0] = bflag;
    state->rstate.ra.ptr.p_double[0] = v;
    state->rstate.ra.ptr.p_double[1] = s;
    state->rstate.ra.ptr.p_double[2] = t;
    return result;
}


/*************************************************************************
Levenberg-Marquardt algorithm results

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    X       -   array[0..N-1], solution
    Rep     -   optimization report;
                see comments for this structure for more info.

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmresults(minlmstate* state,
     /* Real    */ ae_vector* x,
     minlmreport* rep,
     ae_state *_state)
{

    ae_vector_clear(x);
    _minlmreport_clear(rep);

    minlmresultsbuf(state, x, rep, _state);
}


/*************************************************************************
Levenberg-Marquardt algorithm results

Buffered implementation of MinLMResults(), which uses pre-allocated buffer
to store X[]. If buffer size is  too  small,  it  resizes  buffer.  It  is
intended to be used in the inner cycles of performance critical algorithms
where array reallocation penalty is too large to be ignored.

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmresultsbuf(minlmstate* state,
     /* Real    */ ae_vector* x,
     minlmreport* rep,
     ae_state *_state)
{


    if( x->cnt<state->n )
    {
        ae_vector_set_length(x, state->n, _state);
    }
    ae_v_move(&x->ptr.p_double[0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,state->n-1));
    rep->iterationscount = state->repiterationscount;
    rep->terminationtype = state->repterminationtype;
    rep->funcidx = state->repfuncidx;
    rep->varidx = state->repvaridx;
    rep->nfunc = state->repnfunc;
    rep->njac = state->repnjac;
    rep->ngrad = state->repngrad;
    rep->nhess = state->repnhess;
    rep->ncholesky = state->repncholesky;
}


/*************************************************************************
This  subroutine  restarts  LM  algorithm from new point. All optimization
parameters are left unchanged.

This  function  allows  to  solve multiple  optimization  problems  (which
must have same number of dimensions) without object reallocation penalty.

INPUT PARAMETERS:
    State   -   structure used for reverse communication previously
                allocated with MinLMCreateXXX call.
    X       -   new starting point.

  -- ALGLIB --
     Copyright 30.07.2010 by Bochkanov Sergey
*************************************************************************/
void minlmrestartfrom(minlmstate* state,
     /* Real    */ ae_vector* x,
     ae_state *_state)
{


    ae_assert(x->cnt>=state->n, "MinLMRestartFrom: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, state->n, _state), "MinLMRestartFrom: X contains infinite or NaN values!", _state);
    ae_v_move(&state->xbase.ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,state->n-1));
    ae_vector_set_length(&state->rstate.ia, 4+1, _state);
    ae_vector_set_length(&state->rstate.ba, 0+1, _state);
    ae_vector_set_length(&state->rstate.ra, 2+1, _state);
    state->rstate.stage = -1;
    minlm_clearrequestfields(state, _state);
}


/*************************************************************************
This is obsolete function.

Since ALGLIB 3.3 it is equivalent to MinLMCreateVJ().

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatevgj(ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_vector* x,
     minlmstate* state,
     ae_state *_state)
{

    _minlmstate_clear(state);

    minlmcreatevj(n, m, x, state, _state);
}


/*************************************************************************
This is obsolete function.

Since ALGLIB 3.3 it is equivalent to MinLMCreateFJ().

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatefgj(ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_vector* x,
     minlmstate* state,
     ae_state *_state)
{

    _minlmstate_clear(state);

    minlmcreatefj(n, m, x, state, _state);
}


/*************************************************************************
This function is considered obsolete since ALGLIB 3.1.0 and is present for
backward  compatibility  only.  We  recommend  to use MinLMCreateVJ, which
provides similar, but more consistent and feature-rich interface.

  -- ALGLIB --
     Copyright 30.03.2009 by Bochkanov Sergey
*************************************************************************/
void minlmcreatefj(ae_int_t n,
     ae_int_t m,
     /* Real    */ ae_vector* x,
     minlmstate* state,
     ae_state *_state)
{

    _minlmstate_clear(state);

    ae_assert(n>=1, "MinLMCreateFJ: N<1!", _state);
    ae_assert(m>=1, "MinLMCreateFJ: M<1!", _state);
    ae_assert(x->cnt>=n, "MinLMCreateFJ: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "MinLMCreateFJ: X contains infinite or NaN values!", _state);
    
    /*
     * initialize
     */
    state->teststep = 0;
    state->n = n;
    state->m = m;
    state->algomode = 1;
    state->hasf = ae_true;
    state->hasfi = ae_false;
    state->hasg = ae_false;
    
    /*
     * init 2
     */
    minlm_lmprepare(n, m, ae_true, state, _state);
    minlmsetacctype(state, 0, _state);
    minlmsetcond(state, 0, 0, 0, 0, _state);
    minlmsetxrep(state, ae_false, _state);
    minlmsetstpmax(state, 0, _state);
    minlmrestartfrom(state, x, _state);
}


/*************************************************************************
This  subroutine  turns  on  verification  of  the  user-supplied analytic
gradient:
* user calls this subroutine before optimization begins
* MinLMOptimize() is called
* prior to actual optimization, for  each  function Fi and each  component
  of parameters  being  optimized X[j] algorithm performs following steps:
  * two trial steps are made to X[j]-TestStep*S[j] and X[j]+TestStep*S[j],
    where X[j] is j-th parameter and S[j] is a scale of j-th parameter
  * if needed, steps are bounded with respect to constraints on X[]
  * Fi(X) is evaluated at these trial points
  * we perform one more evaluation in the middle point of the interval
  * we  build  cubic  model using function values and derivatives at trial
    points and we compare its prediction with actual value in  the  middle
    point
  * in case difference between prediction and actual value is higher  than
    some predetermined threshold, algorithm stops with completion code -7;
    Rep.VarIdx is set to index of the parameter with incorrect derivative,
    Rep.FuncIdx is set to index of the function.
* after verification is over, algorithm proceeds to the actual optimization.

NOTE 1: verification  needs  N (parameters count) Jacobian evaluations. It
        is  very  costly  and  you  should use it only for low dimensional
        problems,  when  you  want  to  be  sure  that  you've   correctly
        calculated  analytic  derivatives.  You should not  use  it in the
        production code  (unless  you  want  to check derivatives provided
        by some third party).

NOTE 2: you  should  carefully  choose  TestStep. Value which is too large
        (so large that function behaviour is significantly non-cubic) will
        lead to false alarms. You may use  different  step  for  different
        parameters by means of setting scale with MinLMSetScale().

NOTE 3: this function may lead to false positives. In case it reports that
        I-th  derivative was calculated incorrectly, you may decrease test
        step  and  try  one  more  time  - maybe your function changes too
        sharply  and  your  step  is  too  large for such rapidly chanding
        function.

INPUT PARAMETERS:
    State       -   structure used to store algorithm state
    TestStep    -   verification step:
                    * TestStep=0 turns verification off
                    * TestStep>0 activates verification

  -- ALGLIB --
     Copyright 15.06.2012 by Bochkanov Sergey
*************************************************************************/
void minlmsetgradientcheck(minlmstate* state,
     double teststep,
     ae_state *_state)
{


    ae_assert(ae_isfinite(teststep, _state), "MinLMSetGradientCheck: TestStep contains NaN or Infinite", _state);
    ae_assert(ae_fp_greater_eq(teststep,0), "MinLMSetGradientCheck: invalid argument TestStep(TestStep<0)", _state);
    state->teststep = teststep;
}


/*************************************************************************
Prepare internal structures (except for RComm).

Note: M must be zero for FGH mode, non-zero for V/VJ/FJ/FGJ mode.
*************************************************************************/
static void minlm_lmprepare(ae_int_t n,
     ae_int_t m,
     ae_bool havegrad,
     minlmstate* state,
     ae_state *_state)
{
    ae_int_t i;


    if( n<=0||m<0 )
    {
        return;
    }
    if( havegrad )
    {
        ae_vector_set_length(&state->g, n, _state);
    }
    if( m!=0 )
    {
        ae_matrix_set_length(&state->j, m, n, _state);
        ae_vector_set_length(&state->fi, m, _state);
        ae_vector_set_length(&state->fibase, m, _state);
        ae_vector_set_length(&state->deltaf, m, _state);
        ae_vector_set_length(&state->fm1, m, _state);
        ae_vector_set_length(&state->fp1, m, _state);
        ae_vector_set_length(&state->fc1, m, _state);
        ae_vector_set_length(&state->gm1, m, _state);
        ae_vector_set_length(&state->gp1, m, _state);
        ae_vector_set_length(&state->gc1, m, _state);
    }
    else
    {
        ae_matrix_set_length(&state->h, n, n, _state);
    }
    ae_vector_set_length(&state->x, n, _state);
    ae_vector_set_length(&state->deltax, n, _state);
    ae_matrix_set_length(&state->quadraticmodel, n, n, _state);
    ae_vector_set_length(&state->xbase, n, _state);
    ae_vector_set_length(&state->gbase, n, _state);
    ae_vector_set_length(&state->xdir, n, _state);
    ae_vector_set_length(&state->tmp0, n, _state);
    
    /*
     * prepare internal L-BFGS
     */
    for(i=0; i<=n-1; i++)
    {
        state->x.ptr.p_double[i] = 0;
    }
    minlbfgscreate(n, ae_minint(minlm_additers, n, _state), &state->x, &state->internalstate, _state);
    minlbfgssetcond(&state->internalstate, 0.0, 0.0, 0.0, ae_minint(minlm_additers, n, _state), _state);
    
    /*
     * Prepare internal QP solver
     */
    minqpcreate(n, &state->qpstate, _state);
    minqpsetalgocholesky(&state->qpstate, _state);
    
    /*
     * Prepare boundary constraints
     */
    ae_vector_set_length(&state->bndl, n, _state);
    ae_vector_set_length(&state->bndu, n, _state);
    ae_vector_set_length(&state->havebndl, n, _state);
    ae_vector_set_length(&state->havebndu, n, _state);
    for(i=0; i<=n-1; i++)
    {
        state->bndl.ptr.p_double[i] = _state->v_neginf;
        state->havebndl.ptr.p_bool[i] = ae_false;
        state->bndu.ptr.p_double[i] = _state->v_posinf;
        state->havebndu.ptr.p_bool[i] = ae_false;
    }
    
    /*
     * Prepare scaling matrix
     */
    ae_vector_set_length(&state->s, n, _state);
    for(i=0; i<=n-1; i++)
    {
        state->s.ptr.p_double[i] = 1.0;
    }
}


/*************************************************************************
Clears request fileds (to be sure that we don't forgot to clear something)
*************************************************************************/
static void minlm_clearrequestfields(minlmstate* state, ae_state *) //_state)
{


    state->needf = ae_false;
    state->needfg = ae_false;
    state->needfgh = ae_false;
    state->needfij = ae_false;
    state->needfi = ae_false;
    state->xupdated = ae_false;
}


/*************************************************************************
Increases lambda, returns False when there is a danger of overflow
*************************************************************************/
static ae_bool minlm_increaselambda(double* lambdav,
     double* nu,
     ae_state *_state)
{
    double lnlambda;
    double lnnu;
    double lnlambdaup;
    double lnmax;
    ae_bool result;


    result = ae_false;
    lnlambda = ae_log(*lambdav, _state);
    lnlambdaup = ae_log(minlm_lambdaup, _state);
    lnnu = ae_log(*nu, _state);
    lnmax = ae_log(ae_maxrealnumber, _state);
    if( ae_fp_greater(lnlambda+lnlambdaup+lnnu,0.25*lnmax) )
    {
        return result;
    }
    if( ae_fp_greater(lnnu+ae_log(2, _state),lnmax) )
    {
        return result;
    }
    *lambdav = *lambdav*minlm_lambdaup*(*nu);
    *nu = *nu*2;
    result = ae_true;
    return result;
}


/*************************************************************************
Decreases lambda, but leaves it unchanged when there is danger of underflow.
*************************************************************************/
static void minlm_decreaselambda(double* lambdav,
     double* nu,
     ae_state *_state)
{


    *nu = 1;
    if( ae_fp_less(ae_log(*lambdav, _state)+ae_log(minlm_lambdadown, _state),ae_log(ae_minrealnumber, _state)) )
    {
        *lambdav = ae_minrealnumber;
    }
    else
    {
        *lambdav = *lambdav*minlm_lambdadown;
    }
}


/*************************************************************************
Returns norm of bounded scaled anti-gradient.

Bounded antigradient is a vector obtained from  anti-gradient  by  zeroing
components which point outwards:
    result = norm(v)
    v[i]=0     if ((-g[i]<0)and(x[i]=bndl[i])) or
                  ((-g[i]>0)and(x[i]=bndu[i]))
    v[i]=-g[i]*s[i] otherwise, where s[i] is a scale for I-th variable

This function may be used to check a stopping criterion.

  -- ALGLIB --
     Copyright 14.01.2011 by Bochkanov Sergey
*************************************************************************/
static double minlm_boundedscaledantigradnorm(minlmstate* state,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* g,
     ae_state *_state)
{
    ae_int_t n;
    ae_int_t i;
    double v;
    double result;


    result = 0;
    n = state->n;
    for(i=0; i<=n-1; i++)
    {
        v = -g->ptr.p_double[i]*state->s.ptr.p_double[i];
        if( state->havebndl.ptr.p_bool[i] )
        {
            if( ae_fp_less_eq(x->ptr.p_double[i],state->bndl.ptr.p_double[i])&&ae_fp_less(-g->ptr.p_double[i],0) )
            {
                v = 0;
            }
        }
        if( state->havebndu.ptr.p_bool[i] )
        {
            if( ae_fp_greater_eq(x->ptr.p_double[i],state->bndu.ptr.p_double[i])&&ae_fp_greater(-g->ptr.p_double[i],0) )
            {
                v = 0;
            }
        }
        result = result+ae_sqr(v, _state);
    }
    result = ae_sqrt(result, _state);
    return result;
}


ae_bool _minlmstate_init(minlmstate* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->x, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->fi, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->j, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->h, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->g, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xbase, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->fibase, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->gbase, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->quadraticmodel, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->bndl, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->bndu, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->havebndl, 0, DT_BOOL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->havebndu, 0, DT_BOOL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->s, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xdir, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->deltax, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->deltaf, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !_rcommstate_init(&p->rstate, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->choleskybuf, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tmp0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->fm1, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->fp1, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->fc1, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->gm1, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->gp1, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->gc1, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !_minlbfgsstate_init(&p->internalstate, _state, make_automatic) )
        return ae_false;
    if( !_minlbfgsreport_init(&p->internalrep, _state, make_automatic) )
        return ae_false;
    if( !_minqpstate_init(&p->qpstate, _state, make_automatic) )
        return ae_false;
    if( !_minqpreport_init(&p->qprep, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _minlmstate_init_copy(minlmstate* dst, minlmstate* src, ae_state *_state, ae_bool make_automatic)
{
    dst->n = src->n;
    dst->m = src->m;
    dst->diffstep = src->diffstep;
    dst->epsg = src->epsg;
    dst->epsf = src->epsf;
    dst->epsx = src->epsx;
    dst->maxits = src->maxits;
    dst->xrep = src->xrep;
    dst->stpmax = src->stpmax;
    dst->maxmodelage = src->maxmodelage;
    dst->makeadditers = src->makeadditers;
    if( !ae_vector_init_copy(&dst->x, &src->x, _state, make_automatic) )
        return ae_false;
    dst->f = src->f;
    if( !ae_vector_init_copy(&dst->fi, &src->fi, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->j, &src->j, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->h, &src->h, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->g, &src->g, _state, make_automatic) )
        return ae_false;
    dst->needf = src->needf;
    dst->needfg = src->needfg;
    dst->needfgh = src->needfgh;
    dst->needfij = src->needfij;
    dst->needfi = src->needfi;
    dst->xupdated = src->xupdated;
    dst->algomode = src->algomode;
    dst->hasf = src->hasf;
    dst->hasfi = src->hasfi;
    dst->hasg = src->hasg;
    if( !ae_vector_init_copy(&dst->xbase, &src->xbase, _state, make_automatic) )
        return ae_false;
    dst->fbase = src->fbase;
    if( !ae_vector_init_copy(&dst->fibase, &src->fibase, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->gbase, &src->gbase, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->quadraticmodel, &src->quadraticmodel, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->bndl, &src->bndl, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->bndu, &src->bndu, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->havebndl, &src->havebndl, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->havebndu, &src->havebndu, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->s, &src->s, _state, make_automatic) )
        return ae_false;
    dst->lambdav = src->lambdav;
    dst->nu = src->nu;
    dst->modelage = src->modelage;
    if( !ae_vector_init_copy(&dst->xdir, &src->xdir, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->deltax, &src->deltax, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->deltaf, &src->deltaf, _state, make_automatic) )
        return ae_false;
    dst->deltaxready = src->deltaxready;
    dst->deltafready = src->deltafready;
    dst->teststep = src->teststep;
    dst->repiterationscount = src->repiterationscount;
    dst->repterminationtype = src->repterminationtype;
    dst->repfuncidx = src->repfuncidx;
    dst->repvaridx = src->repvaridx;
    dst->repnfunc = src->repnfunc;
    dst->repnjac = src->repnjac;
    dst->repngrad = src->repngrad;
    dst->repnhess = src->repnhess;
    dst->repncholesky = src->repncholesky;
    if( !_rcommstate_init_copy(&dst->rstate, &src->rstate, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->choleskybuf, &src->choleskybuf, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tmp0, &src->tmp0, _state, make_automatic) )
        return ae_false;
    dst->actualdecrease = src->actualdecrease;
    dst->predicteddecrease = src->predicteddecrease;
    dst->xm1 = src->xm1;
    dst->xp1 = src->xp1;
    if( !ae_vector_init_copy(&dst->fm1, &src->fm1, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->fp1, &src->fp1, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->fc1, &src->fc1, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->gm1, &src->gm1, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->gp1, &src->gp1, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->gc1, &src->gc1, _state, make_automatic) )
        return ae_false;
    if( !_minlbfgsstate_init_copy(&dst->internalstate, &src->internalstate, _state, make_automatic) )
        return ae_false;
    if( !_minlbfgsreport_init_copy(&dst->internalrep, &src->internalrep, _state, make_automatic) )
        return ae_false;
    if( !_minqpstate_init_copy(&dst->qpstate, &src->qpstate, _state, make_automatic) )
        return ae_false;
    if( !_minqpreport_init_copy(&dst->qprep, &src->qprep, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _minlmstate_clear(minlmstate* p)
{
    ae_vector_clear(&p->x);
    ae_vector_clear(&p->fi);
    ae_matrix_clear(&p->j);
    ae_matrix_clear(&p->h);
    ae_vector_clear(&p->g);
    ae_vector_clear(&p->xbase);
    ae_vector_clear(&p->fibase);
    ae_vector_clear(&p->gbase);
    ae_matrix_clear(&p->quadraticmodel);
    ae_vector_clear(&p->bndl);
    ae_vector_clear(&p->bndu);
    ae_vector_clear(&p->havebndl);
    ae_vector_clear(&p->havebndu);
    ae_vector_clear(&p->s);
    ae_vector_clear(&p->xdir);
    ae_vector_clear(&p->deltax);
    ae_vector_clear(&p->deltaf);
    _rcommstate_clear(&p->rstate);
    ae_vector_clear(&p->choleskybuf);
    ae_vector_clear(&p->tmp0);
    ae_vector_clear(&p->fm1);
    ae_vector_clear(&p->fp1);
    ae_vector_clear(&p->fc1);
    ae_vector_clear(&p->gm1);
    ae_vector_clear(&p->gp1);
    ae_vector_clear(&p->gc1);
    _minlbfgsstate_clear(&p->internalstate);
    _minlbfgsreport_clear(&p->internalrep);
    _minqpstate_clear(&p->qpstate);
    _minqpreport_clear(&p->qprep);
}


ae_bool _minlmreport_init(minlmreport* /*p*/, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    return ae_true;
}


ae_bool _minlmreport_init_copy(minlmreport* dst, minlmreport* src, ae_state */*_state*/, ae_bool) // make_automatic)
{
    dst->iterationscount = src->iterationscount;
    dst->terminationtype = src->terminationtype;
    dst->funcidx = src->funcidx;
    dst->varidx = src->varidx;
    dst->nfunc = src->nfunc;
    dst->njac = src->njac;
    dst->ngrad = src->ngrad;
    dst->nhess = src->nhess;
    dst->ncholesky = src->ncholesky;
    return ae_true;
}


void _minlmreport_clear(minlmreport* ) //p)
{
}




/*************************************************************************
Obsolete function, use MinLBFGSSetPrecDefault() instead.

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetdefaultpreconditioner(minlbfgsstate* state,
     ae_state *_state)
{


    minlbfgssetprecdefault(state, _state);
}


/*************************************************************************
Obsolete function, use MinLBFGSSetCholeskyPreconditioner() instead.

  -- ALGLIB --
     Copyright 13.10.2010 by Bochkanov Sergey
*************************************************************************/
void minlbfgssetcholeskypreconditioner(minlbfgsstate* state,
     /* Real    */ ae_matrix* p,
     ae_bool isupper,
     ae_state *_state)
{


    minlbfgssetpreccholesky(state, p, isupper, _state);
}


/*************************************************************************
This is obsolete function which was used by previous version of the  BLEIC
optimizer. It does nothing in the current version of BLEIC.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetbarrierwidth(minbleicstate* ,//state,
     double ,//mu,
     ae_state *)//_state)
{


}


/*************************************************************************
This is obsolete function which was used by previous version of the  BLEIC
optimizer. It does nothing in the current version of BLEIC.

  -- ALGLIB --
     Copyright 28.11.2010 by Bochkanov Sergey
*************************************************************************/
void minbleicsetbarrierdecay(minbleicstate*,// state,
     double ,//mudecay,
     ae_state *) //_state)
{


}


/*************************************************************************
Obsolete optimization algorithm.
Was replaced by MinBLEIC subpackage.

  -- ALGLIB --
     Copyright 25.03.2010 by Bochkanov Sergey
*************************************************************************/
void minasacreate(ae_int_t n,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* bndl,
     /* Real    */ ae_vector* bndu,
     minasastate* state,
     ae_state *_state)
{
    ae_int_t i;

    _minasastate_clear(state);

    ae_assert(n>=1, "MinASA: N too small!", _state);
    ae_assert(x->cnt>=n, "MinCGCreate: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, n, _state), "MinCGCreate: X contains infinite or NaN values!", _state);
    ae_assert(bndl->cnt>=n, "MinCGCreate: Length(BndL)<N!", _state);
    ae_assert(isfinitevector(bndl, n, _state), "MinCGCreate: BndL contains infinite or NaN values!", _state);
    ae_assert(bndu->cnt>=n, "MinCGCreate: Length(BndU)<N!", _state);
    ae_assert(isfinitevector(bndu, n, _state), "MinCGCreate: BndU contains infinite or NaN values!", _state);
    for(i=0; i<=n-1; i++)
    {
        ae_assert(ae_fp_less_eq(bndl->ptr.p_double[i],bndu->ptr.p_double[i]), "MinASA: inconsistent bounds!", _state);
        ae_assert(ae_fp_less_eq(bndl->ptr.p_double[i],x->ptr.p_double[i]), "MinASA: infeasible X!", _state);
        ae_assert(ae_fp_less_eq(x->ptr.p_double[i],bndu->ptr.p_double[i]), "MinASA: infeasible X!", _state);
    }
    
    /*
     * Initialize
     */
    state->n = n;
    minasasetcond(state, 0, 0, 0, 0, _state);
    minasasetxrep(state, ae_false, _state);
    minasasetstpmax(state, 0, _state);
    minasasetalgorithm(state, -1, _state);
    ae_vector_set_length(&state->bndl, n, _state);
    ae_vector_set_length(&state->bndu, n, _state);
    ae_vector_set_length(&state->ak, n, _state);
    ae_vector_set_length(&state->xk, n, _state);
    ae_vector_set_length(&state->dk, n, _state);
    ae_vector_set_length(&state->an, n, _state);
    ae_vector_set_length(&state->xn, n, _state);
    ae_vector_set_length(&state->dn, n, _state);
    ae_vector_set_length(&state->x, n, _state);
    ae_vector_set_length(&state->d, n, _state);
    ae_vector_set_length(&state->g, n, _state);
    ae_vector_set_length(&state->gc, n, _state);
    ae_vector_set_length(&state->work, n, _state);
    ae_vector_set_length(&state->yk, n, _state);
    minasarestartfrom(state, x, bndl, bndu, _state);
}


/*************************************************************************
Obsolete optimization algorithm.
Was replaced by MinBLEIC subpackage.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minasasetcond(minasastate* state,
     double epsg,
     double epsf,
     double epsx,
     ae_int_t maxits,
     ae_state *_state)
{


    ae_assert(ae_isfinite(epsg, _state), "MinASASetCond: EpsG is not finite number!", _state);
    ae_assert(ae_fp_greater_eq(epsg,0), "MinASASetCond: negative EpsG!", _state);
    ae_assert(ae_isfinite(epsf, _state), "MinASASetCond: EpsF is not finite number!", _state);
    ae_assert(ae_fp_greater_eq(epsf,0), "MinASASetCond: negative EpsF!", _state);
    ae_assert(ae_isfinite(epsx, _state), "MinASASetCond: EpsX is not finite number!", _state);
    ae_assert(ae_fp_greater_eq(epsx,0), "MinASASetCond: negative EpsX!", _state);
    ae_assert(maxits>=0, "MinASASetCond: negative MaxIts!", _state);
    if( ((ae_fp_eq(epsg,0)&&ae_fp_eq(epsf,0))&&ae_fp_eq(epsx,0))&&maxits==0 )
    {
        epsx = 1.0E-6;
    }
    state->epsg = epsg;
    state->epsf = epsf;
    state->epsx = epsx;
    state->maxits = maxits;
}


/*************************************************************************
Obsolete optimization algorithm.
Was replaced by MinBLEIC subpackage.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minasasetxrep(minasastate* state, ae_bool needxrep, ae_state *) //_state)
{


    state->xrep = needxrep;
}


/*************************************************************************
Obsolete optimization algorithm.
Was replaced by MinBLEIC subpackage.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minasasetalgorithm(minasastate* state,
     ae_int_t algotype,
     ae_state *_state)
{


    ae_assert(algotype>=-1&&algotype<=1, "MinASASetAlgorithm: incorrect AlgoType!", _state);
    if( algotype==-1 )
    {
        algotype = 1;
    }
    state->cgtype = algotype;
}


/*************************************************************************
Obsolete optimization algorithm.
Was replaced by MinBLEIC subpackage.

  -- ALGLIB --
     Copyright 02.04.2010 by Bochkanov Sergey
*************************************************************************/
void minasasetstpmax(minasastate* state, double stpmax, ae_state *_state)
{


    ae_assert(ae_isfinite(stpmax, _state), "MinASASetStpMax: StpMax is not finite!", _state);
    ae_assert(ae_fp_greater_eq(stpmax,0), "MinASASetStpMax: StpMax<0!", _state);
    state->stpmax = stpmax;
}


/*************************************************************************

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool minasaiteration(minasastate* state, ae_state *_state)
{
    ae_int_t n;
    ae_int_t i;
    double betak;
    double v;
    double vv;
    ae_int_t mcinfo;
    ae_bool b;
    ae_bool stepfound;
    ae_int_t diffcnt;
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
        i = state->rstate.ia.ptr.p_int[1];
        mcinfo = state->rstate.ia.ptr.p_int[2];
        diffcnt = state->rstate.ia.ptr.p_int[3];
        b = state->rstate.ba.ptr.p_bool[0];
        stepfound = state->rstate.ba.ptr.p_bool[1];
        betak = state->rstate.ra.ptr.p_double[0];
        v = state->rstate.ra.ptr.p_double[1];
        vv = state->rstate.ra.ptr.p_double[2];
    }
    else
    {
        n = -983;
        i = -989;
        mcinfo = -834;
        diffcnt = 900;
        b = ae_true;
        stepfound = ae_false;
        betak = 214;
        v = -338;
        vv = -686;
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
    if( state->rstate.stage==4 )
    {
        goto lbl_4;
    }
    if( state->rstate.stage==5 )
    {
        goto lbl_5;
    }
    if( state->rstate.stage==6 )
    {
        goto lbl_6;
    }
    if( state->rstate.stage==7 )
    {
        goto lbl_7;
    }
    if( state->rstate.stage==8 )
    {
        goto lbl_8;
    }
    if( state->rstate.stage==9 )
    {
        goto lbl_9;
    }
    if( state->rstate.stage==10 )
    {
        goto lbl_10;
    }
    if( state->rstate.stage==11 )
    {
        goto lbl_11;
    }
    if( state->rstate.stage==12 )
    {
        goto lbl_12;
    }
    if( state->rstate.stage==13 )
    {
        goto lbl_13;
    }
    if( state->rstate.stage==14 )
    {
        goto lbl_14;
    }
    
    /*
     * Routine body
     */
    
    /*
     * Prepare
     */
    n = state->n;
    state->repterminationtype = 0;
    state->repiterationscount = 0;
    state->repnfev = 0;
    state->debugrestartscount = 0;
    state->cgtype = 1;
    ae_v_move(&state->xk.ptr.p_double[0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,n-1));
    for(i=0; i<=n-1; i++)
    {
        if( ae_fp_eq(state->xk.ptr.p_double[i],state->bndl.ptr.p_double[i])||ae_fp_eq(state->xk.ptr.p_double[i],state->bndu.ptr.p_double[i]) )
        {
            state->ak.ptr.p_double[i] = 0;
        }
        else
        {
            state->ak.ptr.p_double[i] = 1;
        }
    }
    state->mu = 0.1;
    state->curalgo = 0;
    
    /*
     * Calculate F/G, initialize algorithm
     */
    mincomp_clearrequestfields(state, _state);
    state->needfg = ae_true;
    state->rstate.stage = 0;
    goto lbl_rcomm;
lbl_0:
    state->needfg = ae_false;
    if( !state->xrep )
    {
        goto lbl_15;
    }
    
    /*
     * progress report
     */
    mincomp_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 1;
    goto lbl_rcomm;
lbl_1:
    state->xupdated = ae_false;
lbl_15:
    if( ae_fp_less_eq(mincomp_asaboundedantigradnorm(state, _state),state->epsg) )
    {
        state->repterminationtype = 4;
        result = ae_false;
        return result;
    }
    state->repnfev = state->repnfev+1;
    
    /*
     * Main cycle
     *
     * At the beginning of new iteration:
     * * CurAlgo stores current algorithm selector
     * * State.XK, State.F and State.G store current X/F/G
     * * State.AK stores current set of active constraints
     */
lbl_17:
    if( ae_false )
    {
        goto lbl_18;
    }
    
    /*
     * GPA algorithm
     */
    if( state->curalgo!=0 )
    {
        goto lbl_19;
    }
    state->k = 0;
    state->acount = 0;
lbl_21:
    if( ae_false )
    {
        goto lbl_22;
    }
    
    /*
     * Determine Dk = proj(xk - gk)-xk
     */
    for(i=0; i<=n-1; i++)
    {
        state->d.ptr.p_double[i] = boundval(state->xk.ptr.p_double[i]-state->g.ptr.p_double[i], state->bndl.ptr.p_double[i], state->bndu.ptr.p_double[i], _state)-state->xk.ptr.p_double[i];
    }
    
    /*
     * Armijo line search.
     * * exact search with alpha=1 is tried first,
     *   'exact' means that we evaluate f() EXACTLY at
     *   bound(x-g,bndl,bndu), without intermediate floating
     *   point operations.
     * * alpha<1 are tried if explicit search wasn't successful
     * Result is placed into XN.
     *
     * Two types of search are needed because we can't
     * just use second type with alpha=1 because in finite
     * precision arithmetics (x1-x0)+x0 may differ from x1.
     * So while x1 is correctly bounded (it lie EXACTLY on
     * boundary, if it is active), (x1-x0)+x0 may be
     * not bounded.
     */
    v = ae_v_dotproduct(&state->d.ptr.p_double[0], 1, &state->g.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->dginit = v;
    state->finit = state->f;
    if( !(ae_fp_less_eq(mincomp_asad1norm(state, _state),state->stpmax)||ae_fp_eq(state->stpmax,0)) )
    {
        goto lbl_23;
    }
    
    /*
     * Try alpha=1 step first
     */
    for(i=0; i<=n-1; i++)
    {
        state->x.ptr.p_double[i] = boundval(state->xk.ptr.p_double[i]-state->g.ptr.p_double[i], state->bndl.ptr.p_double[i], state->bndu.ptr.p_double[i], _state);
    }
    mincomp_clearrequestfields(state, _state);
    state->needfg = ae_true;
    state->rstate.stage = 2;
    goto lbl_rcomm;
lbl_2:
    state->needfg = ae_false;
    state->repnfev = state->repnfev+1;
    stepfound = ae_fp_less_eq(state->f,state->finit+mincomp_gpaftol*state->dginit);
    goto lbl_24;
lbl_23:
    stepfound = ae_false;
lbl_24:
    if( !stepfound )
    {
        goto lbl_25;
    }
    
    /*
     * we are at the boundary(ies)
     */
    ae_v_move(&state->xn.ptr.p_double[0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->stp = 1;
    goto lbl_26;
lbl_25:
    
    /*
     * alpha=1 is too large, try smaller values
     */
    state->stp = 1;
    linminnormalized(&state->d, &state->stp, n, _state);
    state->dginit = state->dginit/state->stp;
    state->stp = mincomp_gpadecay*state->stp;
    if( ae_fp_greater(state->stpmax,0) )
    {
        state->stp = ae_minreal(state->stp, state->stpmax, _state);
    }
lbl_27:
    if( ae_false )
    {
        goto lbl_28;
    }
    v = state->stp;
    ae_v_move(&state->x.ptr.p_double[0], 1, &state->xk.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_addd(&state->x.ptr.p_double[0], 1, &state->d.ptr.p_double[0], 1, ae_v_len(0,n-1), v);
    mincomp_clearrequestfields(state, _state);
    state->needfg = ae_true;
    state->rstate.stage = 3;
    goto lbl_rcomm;
lbl_3:
    state->needfg = ae_false;
    state->repnfev = state->repnfev+1;
    if( ae_fp_less_eq(state->stp,mincomp_stpmin) )
    {
        goto lbl_28;
    }
    if( ae_fp_less_eq(state->f,state->finit+state->stp*mincomp_gpaftol*state->dginit) )
    {
        goto lbl_28;
    }
    state->stp = state->stp*mincomp_gpadecay;
    goto lbl_27;
lbl_28:
    ae_v_move(&state->xn.ptr.p_double[0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,n-1));
lbl_26:
    state->repiterationscount = state->repiterationscount+1;
    if( !state->xrep )
    {
        goto lbl_29;
    }
    
    /*
     * progress report
     */
    mincomp_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 4;
    goto lbl_rcomm;
lbl_4:
    state->xupdated = ae_false;
lbl_29:
    
    /*
     * Calculate new set of active constraints.
     * Reset counter if active set was changed.
     * Prepare for the new iteration
     */
    for(i=0; i<=n-1; i++)
    {
        if( ae_fp_eq(state->xn.ptr.p_double[i],state->bndl.ptr.p_double[i])||ae_fp_eq(state->xn.ptr.p_double[i],state->bndu.ptr.p_double[i]) )
        {
            state->an.ptr.p_double[i] = 0;
        }
        else
        {
            state->an.ptr.p_double[i] = 1;
        }
    }
    for(i=0; i<=n-1; i++)
    {
        if( ae_fp_neq(state->ak.ptr.p_double[i],state->an.ptr.p_double[i]) )
        {
            state->acount = -1;
            break;
        }
    }
    state->acount = state->acount+1;
    ae_v_move(&state->xk.ptr.p_double[0], 1, &state->xn.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_move(&state->ak.ptr.p_double[0], 1, &state->an.ptr.p_double[0], 1, ae_v_len(0,n-1));
    
    /*
     * Stopping conditions
     */
    if( !(state->repiterationscount>=state->maxits&&state->maxits>0) )
    {
        goto lbl_31;
    }
    
    /*
     * Too many iterations
     */
    state->repterminationtype = 5;
    if( !state->xrep )
    {
        goto lbl_33;
    }
    mincomp_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 5;
    goto lbl_rcomm;
lbl_5:
    state->xupdated = ae_false;
lbl_33:
    result = ae_false;
    return result;
lbl_31:
    if( ae_fp_greater(mincomp_asaboundedantigradnorm(state, _state),state->epsg) )
    {
        goto lbl_35;
    }
    
    /*
     * Gradient is small enough
     */
    state->repterminationtype = 4;
    if( !state->xrep )
    {
        goto lbl_37;
    }
    mincomp_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 6;
    goto lbl_rcomm;
lbl_6:
    state->xupdated = ae_false;
lbl_37:
    result = ae_false;
    return result;
lbl_35:
    v = ae_v_dotproduct(&state->d.ptr.p_double[0], 1, &state->d.ptr.p_double[0], 1, ae_v_len(0,n-1));
    if( ae_fp_greater(ae_sqrt(v, _state)*state->stp,state->epsx) )
    {
        goto lbl_39;
    }
    
    /*
     * Step size is too small, no further improvement is
     * possible
     */
    state->repterminationtype = 2;
    if( !state->xrep )
    {
        goto lbl_41;
    }
    mincomp_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 7;
    goto lbl_rcomm;
lbl_7:
    state->xupdated = ae_false;
lbl_41:
    result = ae_false;
    return result;
lbl_39:
    if( ae_fp_greater(state->finit-state->f,state->epsf*ae_maxreal(ae_fabs(state->finit, _state), ae_maxreal(ae_fabs(state->f, _state), 1.0, _state), _state)) )
    {
        goto lbl_43;
    }
    
    /*
     * F(k+1)-F(k) is small enough
     */
    state->repterminationtype = 1;
    if( !state->xrep )
    {
        goto lbl_45;
    }
    mincomp_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 8;
    goto lbl_rcomm;
lbl_8:
    state->xupdated = ae_false;
lbl_45:
    result = ae_false;
    return result;
lbl_43:
    
    /*
     * Decide - should we switch algorithm or not
     */
    if( mincomp_asauisempty(state, _state) )
    {
        if( ae_fp_greater_eq(mincomp_asaginorm(state, _state),state->mu*mincomp_asad1norm(state, _state)) )
        {
            state->curalgo = 1;
            goto lbl_22;
        }
        else
        {
            state->mu = state->mu*mincomp_asarho;
        }
    }
    else
    {
        if( state->acount==mincomp_n1 )
        {
            if( ae_fp_greater_eq(mincomp_asaginorm(state, _state),state->mu*mincomp_asad1norm(state, _state)) )
            {
                state->curalgo = 1;
                goto lbl_22;
            }
        }
    }
    
    /*
     * Next iteration
     */
    state->k = state->k+1;
    goto lbl_21;
lbl_22:
lbl_19:
    
    /*
     * CG algorithm
     */
    if( state->curalgo!=1 )
    {
        goto lbl_47;
    }
    
    /*
     * first, check that there are non-active constraints.
     * move to GPA algorithm, if all constraints are active
     */
    b = ae_true;
    for(i=0; i<=n-1; i++)
    {
        if( ae_fp_neq(state->ak.ptr.p_double[i],0) )
        {
            b = ae_false;
            break;
        }
    }
    if( b )
    {
        state->curalgo = 0;
        goto lbl_17;
    }
    
    /*
     * CG iterations
     */
    state->fold = state->f;
    ae_v_move(&state->xk.ptr.p_double[0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,n-1));
    for(i=0; i<=n-1; i++)
    {
        state->dk.ptr.p_double[i] = -state->g.ptr.p_double[i]*state->ak.ptr.p_double[i];
        state->gc.ptr.p_double[i] = state->g.ptr.p_double[i]*state->ak.ptr.p_double[i];
    }
lbl_49:
    if( ae_false )
    {
        goto lbl_50;
    }
    
    /*
     * Store G[k] for later calculation of Y[k]
     */
    for(i=0; i<=n-1; i++)
    {
        state->yk.ptr.p_double[i] = -state->gc.ptr.p_double[i];
    }
    
    /*
     * Make a CG step in direction given by DK[]:
     * * calculate step. Step projection into feasible set
     *   is used. It has several benefits: a) step may be
     *   found with usual line search, b) multiple constraints
     *   may be activated with one step, c) activated constraints
     *   are detected in a natural way - just compare x[i] with
     *   bounds
     * * update active set, set B to True, if there
     *   were changes in the set.
     */
    ae_v_move(&state->d.ptr.p_double[0], 1, &state->dk.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_move(&state->xn.ptr.p_double[0], 1, &state->xk.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->mcstage = 0;
    state->stp = 1;
    linminnormalized(&state->d, &state->stp, n, _state);
    if( ae_fp_neq(state->laststep,0) )
    {
        state->stp = state->laststep;
    }
    mcsrch(n, &state->xn, &state->f, &state->gc, &state->d, &state->stp, state->stpmax, mincomp_gtol, &mcinfo, &state->nfev, &state->work, &state->lstate, &state->mcstage, _state);
lbl_51:
    if( state->mcstage==0 )
    {
        goto lbl_52;
    }
    
    /*
     * preprocess data: bound State.XN so it belongs to the
     * feasible set and store it in the State.X
     */
    for(i=0; i<=n-1; i++)
    {
        state->x.ptr.p_double[i] = boundval(state->xn.ptr.p_double[i], state->bndl.ptr.p_double[i], state->bndu.ptr.p_double[i], _state);
    }
    
    /*
     * RComm
     */
    mincomp_clearrequestfields(state, _state);
    state->needfg = ae_true;
    state->rstate.stage = 9;
    goto lbl_rcomm;
lbl_9:
    state->needfg = ae_false;
    
    /*
     * postprocess data: zero components of G corresponding to
     * the active constraints
     */
    for(i=0; i<=n-1; i++)
    {
        if( ae_fp_eq(state->x.ptr.p_double[i],state->bndl.ptr.p_double[i])||ae_fp_eq(state->x.ptr.p_double[i],state->bndu.ptr.p_double[i]) )
        {
            state->gc.ptr.p_double[i] = 0;
        }
        else
        {
            state->gc.ptr.p_double[i] = state->g.ptr.p_double[i];
        }
    }
    mcsrch(n, &state->xn, &state->f, &state->gc, &state->d, &state->stp, state->stpmax, mincomp_gtol, &mcinfo, &state->nfev, &state->work, &state->lstate, &state->mcstage, _state);
    goto lbl_51;
lbl_52:
    diffcnt = 0;
    for(i=0; i<=n-1; i++)
    {
        
        /*
         * XN contains unprojected result, project it,
         * save copy to X (will be used for progress reporting)
         */
        state->xn.ptr.p_double[i] = boundval(state->xn.ptr.p_double[i], state->bndl.ptr.p_double[i], state->bndu.ptr.p_double[i], _state);
        
        /*
         * update active set
         */
        if( ae_fp_eq(state->xn.ptr.p_double[i],state->bndl.ptr.p_double[i])||ae_fp_eq(state->xn.ptr.p_double[i],state->bndu.ptr.p_double[i]) )
        {
            state->an.ptr.p_double[i] = 0;
        }
        else
        {
            state->an.ptr.p_double[i] = 1;
        }
        if( ae_fp_neq(state->an.ptr.p_double[i],state->ak.ptr.p_double[i]) )
        {
            diffcnt = diffcnt+1;
        }
        state->ak.ptr.p_double[i] = state->an.ptr.p_double[i];
    }
    ae_v_move(&state->xk.ptr.p_double[0], 1, &state->xn.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->repnfev = state->repnfev+state->nfev;
    state->repiterationscount = state->repiterationscount+1;
    if( !state->xrep )
    {
        goto lbl_53;
    }
    
    /*
     * progress report
     */
    mincomp_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 10;
    goto lbl_rcomm;
lbl_10:
    state->xupdated = ae_false;
lbl_53:
    
    /*
     * Update info about step length
     */
    v = ae_v_dotproduct(&state->d.ptr.p_double[0], 1, &state->d.ptr.p_double[0], 1, ae_v_len(0,n-1));
    state->laststep = ae_sqrt(v, _state)*state->stp;
    
    /*
     * Check stopping conditions.
     */
    if( ae_fp_greater(mincomp_asaboundedantigradnorm(state, _state),state->epsg) )
    {
        goto lbl_55;
    }
    
    /*
     * Gradient is small enough
     */
    state->repterminationtype = 4;
    if( !state->xrep )
    {
        goto lbl_57;
    }
    mincomp_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 11;
    goto lbl_rcomm;
lbl_11:
    state->xupdated = ae_false;
lbl_57:
    result = ae_false;
    return result;
lbl_55:
    if( !(state->repiterationscount>=state->maxits&&state->maxits>0) )
    {
        goto lbl_59;
    }
    
    /*
     * Too many iterations
     */
    state->repterminationtype = 5;
    if( !state->xrep )
    {
        goto lbl_61;
    }
    mincomp_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 12;
    goto lbl_rcomm;
lbl_12:
    state->xupdated = ae_false;
lbl_61:
    result = ae_false;
    return result;
lbl_59:
    if( !(ae_fp_greater_eq(mincomp_asaginorm(state, _state),state->mu*mincomp_asad1norm(state, _state))&&diffcnt==0) )
    {
        goto lbl_63;
    }
    
    /*
     * These conditions (EpsF/EpsX) are explicitly or implicitly
     * related to the current step size and influenced
     * by changes in the active constraints.
     *
     * For these reasons they are checked only when we don't
     * want to 'unstick' at the end of the iteration and there
     * were no changes in the active set.
     *
     * NOTE: consition |G|>=Mu*|D1| must be exactly opposite
     * to the condition used to switch back to GPA. At least
     * one inequality must be strict, otherwise infinite cycle
     * may occur when |G|=Mu*|D1| (we DON'T test stopping
     * conditions and we DON'T switch to GPA, so we cycle
     * indefinitely).
     */
    if( ae_fp_greater(state->fold-state->f,state->epsf*ae_maxreal(ae_fabs(state->fold, _state), ae_maxreal(ae_fabs(state->f, _state), 1.0, _state), _state)) )
    {
        goto lbl_65;
    }
    
    /*
     * F(k+1)-F(k) is small enough
     */
    state->repterminationtype = 1;
    if( !state->xrep )
    {
        goto lbl_67;
    }
    mincomp_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 13;
    goto lbl_rcomm;
lbl_13:
    state->xupdated = ae_false;
lbl_67:
    result = ae_false;
    return result;
lbl_65:
    if( ae_fp_greater(state->laststep,state->epsx) )
    {
        goto lbl_69;
    }
    
    /*
     * X(k+1)-X(k) is small enough
     */
    state->repterminationtype = 2;
    if( !state->xrep )
    {
        goto lbl_71;
    }
    mincomp_clearrequestfields(state, _state);
    state->xupdated = ae_true;
    state->rstate.stage = 14;
    goto lbl_rcomm;
lbl_14:
    state->xupdated = ae_false;
lbl_71:
    result = ae_false;
    return result;
lbl_69:
lbl_63:
    
    /*
     * Check conditions for switching
     */
    if( ae_fp_less(mincomp_asaginorm(state, _state),state->mu*mincomp_asad1norm(state, _state)) )
    {
        state->curalgo = 0;
        goto lbl_50;
    }
    if( diffcnt>0 )
    {
        if( mincomp_asauisempty(state, _state)||diffcnt>=mincomp_n2 )
        {
            state->curalgo = 1;
        }
        else
        {
            state->curalgo = 0;
        }
        goto lbl_50;
    }
    
    /*
     * Calculate D(k+1)
     *
     * Line search may result in:
     * * maximum feasible step being taken (already processed)
     * * point satisfying Wolfe conditions
     * * some kind of error (CG is restarted by assigning 0.0 to Beta)
     */
    if( mcinfo==1 )
    {
        
        /*
         * Standard Wolfe conditions are satisfied:
         * * calculate Y[K] and BetaK
         */
        ae_v_add(&state->yk.ptr.p_double[0], 1, &state->gc.ptr.p_double[0], 1, ae_v_len(0,n-1));
        vv = ae_v_dotproduct(&state->yk.ptr.p_double[0], 1, &state->dk.ptr.p_double[0], 1, ae_v_len(0,n-1));
        v = ae_v_dotproduct(&state->gc.ptr.p_double[0], 1, &state->gc.ptr.p_double[0], 1, ae_v_len(0,n-1));
        state->betady = v/vv;
        v = ae_v_dotproduct(&state->gc.ptr.p_double[0], 1, &state->yk.ptr.p_double[0], 1, ae_v_len(0,n-1));
        state->betahs = v/vv;
        if( state->cgtype==0 )
        {
            betak = state->betady;
        }
        if( state->cgtype==1 )
        {
            betak = ae_maxreal(0, ae_minreal(state->betady, state->betahs, _state), _state);
        }
    }
    else
    {
        
        /*
         * Something is wrong (may be function is too wild or too flat).
         *
         * We'll set BetaK=0, which will restart CG algorithm.
         * We can stop later (during normal checks) if stopping conditions are met.
         */
        betak = 0;
        state->debugrestartscount = state->debugrestartscount+1;
    }
    ae_v_moveneg(&state->dn.ptr.p_double[0], 1, &state->gc.ptr.p_double[0], 1, ae_v_len(0,n-1));
    ae_v_addd(&state->dn.ptr.p_double[0], 1, &state->dk.ptr.p_double[0], 1, ae_v_len(0,n-1), betak);
    ae_v_move(&state->dk.ptr.p_double[0], 1, &state->dn.ptr.p_double[0], 1, ae_v_len(0,n-1));
    
    /*
     * update other information
     */
    state->fold = state->f;
    state->k = state->k+1;
    goto lbl_49;
lbl_50:
lbl_47:
    goto lbl_17;
lbl_18:
    result = ae_false;
    return result;
    
    /*
     * Saving state
     */
lbl_rcomm:
    result = ae_true;
    state->rstate.ia.ptr.p_int[0] = n;
    state->rstate.ia.ptr.p_int[1] = i;
    state->rstate.ia.ptr.p_int[2] = mcinfo;
    state->rstate.ia.ptr.p_int[3] = diffcnt;
    state->rstate.ba.ptr.p_bool[0] = b;
    state->rstate.ba.ptr.p_bool[1] = stepfound;
    state->rstate.ra.ptr.p_double[0] = betak;
    state->rstate.ra.ptr.p_double[1] = v;
    state->rstate.ra.ptr.p_double[2] = vv;
    return result;
}


/*************************************************************************
Obsolete optimization algorithm.
Was replaced by MinBLEIC subpackage.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void minasaresults(minasastate* state,
     /* Real    */ ae_vector* x,
     minasareport* rep,
     ae_state *_state)
{

    ae_vector_clear(x);
    _minasareport_clear(rep);

    minasaresultsbuf(state, x, rep, _state);
}


/*************************************************************************
Obsolete optimization algorithm.
Was replaced by MinBLEIC subpackage.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
void minasaresultsbuf(minasastate* state,
     /* Real    */ ae_vector* x,
     minasareport* rep,
     ae_state *_state)
{
    ae_int_t i;


    if( x->cnt<state->n )
    {
        ae_vector_set_length(x, state->n, _state);
    }
    ae_v_move(&x->ptr.p_double[0], 1, &state->x.ptr.p_double[0], 1, ae_v_len(0,state->n-1));
    rep->iterationscount = state->repiterationscount;
    rep->nfev = state->repnfev;
    rep->terminationtype = state->repterminationtype;
    rep->activeconstraints = 0;
    for(i=0; i<=state->n-1; i++)
    {
        if( ae_fp_eq(state->ak.ptr.p_double[i],0) )
        {
            rep->activeconstraints = rep->activeconstraints+1;
        }
    }
}


/*************************************************************************
Obsolete optimization algorithm.
Was replaced by MinBLEIC subpackage.

  -- ALGLIB --
     Copyright 30.07.2010 by Bochkanov Sergey
*************************************************************************/
void minasarestartfrom(minasastate* state,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* bndl,
     /* Real    */ ae_vector* bndu,
     ae_state *_state)
{


    ae_assert(x->cnt>=state->n, "MinASARestartFrom: Length(X)<N!", _state);
    ae_assert(isfinitevector(x, state->n, _state), "MinASARestartFrom: X contains infinite or NaN values!", _state);
    ae_assert(bndl->cnt>=state->n, "MinASARestartFrom: Length(BndL)<N!", _state);
    ae_assert(isfinitevector(bndl, state->n, _state), "MinASARestartFrom: BndL contains infinite or NaN values!", _state);
    ae_assert(bndu->cnt>=state->n, "MinASARestartFrom: Length(BndU)<N!", _state);
    ae_assert(isfinitevector(bndu, state->n, _state), "MinASARestartFrom: BndU contains infinite or NaN values!", _state);
    ae_v_move(&state->x.ptr.p_double[0], 1, &x->ptr.p_double[0], 1, ae_v_len(0,state->n-1));
    ae_v_move(&state->bndl.ptr.p_double[0], 1, &bndl->ptr.p_double[0], 1, ae_v_len(0,state->n-1));
    ae_v_move(&state->bndu.ptr.p_double[0], 1, &bndu->ptr.p_double[0], 1, ae_v_len(0,state->n-1));
    state->laststep = 0;
    ae_vector_set_length(&state->rstate.ia, 3+1, _state);
    ae_vector_set_length(&state->rstate.ba, 1+1, _state);
    ae_vector_set_length(&state->rstate.ra, 2+1, _state);
    state->rstate.stage = -1;
    mincomp_clearrequestfields(state, _state);
}


/*************************************************************************
Returns norm of bounded anti-gradient.

Bounded antigradient is a vector obtained from  anti-gradient  by  zeroing
components which point outwards:
    result = norm(v)
    v[i]=0     if ((-g[i]<0)and(x[i]=bndl[i])) or
                  ((-g[i]>0)and(x[i]=bndu[i]))
    v[i]=-g[i] otherwise

This function may be used to check a stopping criterion.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
static double mincomp_asaboundedantigradnorm(minasastate* state,
     ae_state *_state)
{
    ae_int_t i;
    double v;
    double result;


    result = 0;
    for(i=0; i<=state->n-1; i++)
    {
        v = -state->g.ptr.p_double[i];
        if( ae_fp_eq(state->x.ptr.p_double[i],state->bndl.ptr.p_double[i])&&ae_fp_less(-state->g.ptr.p_double[i],0) )
        {
            v = 0;
        }
        if( ae_fp_eq(state->x.ptr.p_double[i],state->bndu.ptr.p_double[i])&&ae_fp_greater(-state->g.ptr.p_double[i],0) )
        {
            v = 0;
        }
        result = result+ae_sqr(v, _state);
    }
    result = ae_sqrt(result, _state);
    return result;
}


/*************************************************************************
Returns norm of GI(x).

GI(x) is  a  gradient  vector  whose  components  associated  with  active
constraints are zeroed. It  differs  from  bounded  anti-gradient  because
components  of   GI(x)   are   zeroed  independently  of  sign(g[i]),  and
anti-gradient's components are zeroed with respect to both constraint  and
sign.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
static double mincomp_asaginorm(minasastate* state, ae_state *_state)
{
    ae_int_t i;
    double result;


    result = 0;
    for(i=0; i<=state->n-1; i++)
    {
        if( ae_fp_neq(state->x.ptr.p_double[i],state->bndl.ptr.p_double[i])&&ae_fp_neq(state->x.ptr.p_double[i],state->bndu.ptr.p_double[i]) )
        {
            result = result+ae_sqr(state->g.ptr.p_double[i], _state);
        }
    }
    result = ae_sqrt(result, _state);
    return result;
}


/*************************************************************************
Returns norm(D1(State.X))

For a meaning of D1 see 'NEW ACTIVE SET ALGORITHM FOR BOX CONSTRAINED
OPTIMIZATION' by WILLIAM W. HAGER AND HONGCHAO ZHANG.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
static double mincomp_asad1norm(minasastate* state, ae_state *_state)
{
    ae_int_t i;
    double result;


    result = 0;
    for(i=0; i<=state->n-1; i++)
    {
        result = result+ae_sqr(boundval(state->x.ptr.p_double[i]-state->g.ptr.p_double[i], state->bndl.ptr.p_double[i], state->bndu.ptr.p_double[i], _state)-state->x.ptr.p_double[i], _state);
    }
    result = ae_sqrt(result, _state);
    return result;
}


/*************************************************************************
Returns True, if U set is empty.

* State.X is used as point,
* State.G - as gradient,
* D is calculated within function (because State.D may have different
  meaning depending on current optimization algorithm)

For a meaning of U see 'NEW ACTIVE SET ALGORITHM FOR BOX CONSTRAINED
OPTIMIZATION' by WILLIAM W. HAGER AND HONGCHAO ZHANG.

  -- ALGLIB --
     Copyright 20.03.2009 by Bochkanov Sergey
*************************************************************************/
static ae_bool mincomp_asauisempty(minasastate* state, ae_state *_state)
{
    ae_int_t i;
    double d;
    double d2;
    double d32;
    ae_bool result;


    d = mincomp_asad1norm(state, _state);
    d2 = ae_sqrt(d, _state);
    d32 = d*d2;
    result = ae_true;
    for(i=0; i<=state->n-1; i++)
    {
        if( ae_fp_greater_eq(ae_fabs(state->g.ptr.p_double[i], _state),d2)&&ae_fp_greater_eq(ae_minreal(state->x.ptr.p_double[i]-state->bndl.ptr.p_double[i], state->bndu.ptr.p_double[i]-state->x.ptr.p_double[i], _state),d32) )
        {
            result = ae_false;
            return result;
        }
    }
    return result;
}


/*************************************************************************
Clears request fileds (to be sure that we don't forgot to clear something)
*************************************************************************/
static void mincomp_clearrequestfields(minasastate* state,
     ae_state *) //_state)
{


    state->needfg = ae_false;
    state->xupdated = ae_false;
}


ae_bool _minasastate_init(minasastate* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->bndl, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->bndu, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->ak, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xk, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->dk, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->an, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->xn, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->dn, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->d, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->work, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->yk, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->gc, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->x, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->g, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !_rcommstate_init(&p->rstate, _state, make_automatic) )
        return ae_false;
    if( !_linminstate_init(&p->lstate, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _minasastate_init_copy(minasastate* dst, minasastate* src, ae_state *_state, ae_bool make_automatic)
{
    dst->n = src->n;
    dst->epsg = src->epsg;
    dst->epsf = src->epsf;
    dst->epsx = src->epsx;
    dst->maxits = src->maxits;
    dst->xrep = src->xrep;
    dst->stpmax = src->stpmax;
    dst->cgtype = src->cgtype;
    dst->k = src->k;
    dst->nfev = src->nfev;
    dst->mcstage = src->mcstage;
    if( !ae_vector_init_copy(&dst->bndl, &src->bndl, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->bndu, &src->bndu, _state, make_automatic) )
        return ae_false;
    dst->curalgo = src->curalgo;
    dst->acount = src->acount;
    dst->mu = src->mu;
    dst->finit = src->finit;
    dst->dginit = src->dginit;
    if( !ae_vector_init_copy(&dst->ak, &src->ak, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->xk, &src->xk, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->dk, &src->dk, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->an, &src->an, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->xn, &src->xn, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->dn, &src->dn, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->d, &src->d, _state, make_automatic) )
        return ae_false;
    dst->fold = src->fold;
    dst->stp = src->stp;
    if( !ae_vector_init_copy(&dst->work, &src->work, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->yk, &src->yk, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->gc, &src->gc, _state, make_automatic) )
        return ae_false;
    dst->laststep = src->laststep;
    if( !ae_vector_init_copy(&dst->x, &src->x, _state, make_automatic) )
        return ae_false;
    dst->f = src->f;
    if( !ae_vector_init_copy(&dst->g, &src->g, _state, make_automatic) )
        return ae_false;
    dst->needfg = src->needfg;
    dst->xupdated = src->xupdated;
    if( !_rcommstate_init_copy(&dst->rstate, &src->rstate, _state, make_automatic) )
        return ae_false;
    dst->repiterationscount = src->repiterationscount;
    dst->repnfev = src->repnfev;
    dst->repterminationtype = src->repterminationtype;
    dst->debugrestartscount = src->debugrestartscount;
    if( !_linminstate_init_copy(&dst->lstate, &src->lstate, _state, make_automatic) )
        return ae_false;
    dst->betahs = src->betahs;
    dst->betady = src->betady;
    return ae_true;
}


void _minasastate_clear(minasastate* p)
{
    ae_vector_clear(&p->bndl);
    ae_vector_clear(&p->bndu);
    ae_vector_clear(&p->ak);
    ae_vector_clear(&p->xk);
    ae_vector_clear(&p->dk);
    ae_vector_clear(&p->an);
    ae_vector_clear(&p->xn);
    ae_vector_clear(&p->dn);
    ae_vector_clear(&p->d);
    ae_vector_clear(&p->work);
    ae_vector_clear(&p->yk);
    ae_vector_clear(&p->gc);
    ae_vector_clear(&p->x);
    ae_vector_clear(&p->g);
    _rcommstate_clear(&p->rstate);
    _linminstate_clear(&p->lstate);
}


ae_bool _minasareport_init(minasareport* /*p*/, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    return ae_true;
}


ae_bool _minasareport_init_copy(minasareport* dst, minasareport* src, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    dst->iterationscount = src->iterationscount;
    dst->nfev = src->nfev;
    dst->terminationtype = src->terminationtype;
    dst->activeconstraints = src->activeconstraints;
    return ae_true;
}


void _minasareport_clear(minasareport*) //p)
{
}



}

