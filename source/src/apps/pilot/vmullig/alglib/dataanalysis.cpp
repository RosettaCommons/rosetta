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
#include "dataanalysis.h"

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
Optimal binary classification

Algorithms finds optimal (=with minimal cross-entropy) binary partition.
Internal subroutine.

INPUT PARAMETERS:
    A       -   array[0..N-1], variable
    C       -   array[0..N-1], class numbers (0 or 1).
    N       -   array size

OUTPUT PARAMETERS:
    Info    -   completetion code:
                * -3, all values of A[] are same (partition is impossible)
                * -2, one of C[] is incorrect (<0, >1)
                * -1, incorrect pararemets were passed (N<=0).
                *  1, OK
    Threshold-  partiton boundary. Left part contains values which are
                strictly less than Threshold. Right part contains values
                which are greater than or equal to Threshold.
    PAL, PBL-   probabilities P(0|v<Threshold) and P(1|v<Threshold)
    PAR, PBR-   probabilities P(0|v>=Threshold) and P(1|v>=Threshold)
    CVE     -   cross-validation estimate of cross-entropy

  -- ALGLIB --
     Copyright 22.05.2008 by Bochkanov Sergey
*************************************************************************/
void dsoptimalsplit2(const real_1d_array &a, const integer_1d_array &c, const ae_int_t n, ae_int_t &info, double &threshold, double &pal, double &pbl, double &par, double &pbr, double &cve)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::dsoptimalsplit2(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), const_cast<alglib_impl::ae_vector*>(c.c_ptr()), n, &info, &threshold, &pal, &pbl, &par, &pbr, &cve, &_alglib_env_state);
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
Optimal partition, internal subroutine. Fast version.

Accepts:
    A       array[0..N-1]       array of attributes     array[0..N-1]
    C       array[0..N-1]       array of class labels
    TiesBuf array[0..N]         temporaries (ties)
    CntBuf  array[0..2*NC-1]    temporaries (counts)
    Alpha                       centering factor (0<=alpha<=1, recommended value - 0.05)
    BufR    array[0..N-1]       temporaries
    BufI    array[0..N-1]       temporaries

Output:
    Info    error code (">0"=OK, "<0"=bad)
    RMS     training set RMS error
    CVRMS   leave-one-out RMS error

Note:
    content of all arrays is changed by subroutine;
    it doesn't allocate temporaries.

  -- ALGLIB --
     Copyright 11.12.2008 by Bochkanov Sergey
*************************************************************************/
void dsoptimalsplit2fast(real_1d_array &a, integer_1d_array &c, integer_1d_array &tiesbuf, integer_1d_array &cntbuf, real_1d_array &bufr, integer_1d_array &bufi, const ae_int_t n, const ae_int_t nc, const double alpha, ae_int_t &info, double &threshold, double &rms, double &cvrms)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::dsoptimalsplit2fast(const_cast<alglib_impl::ae_vector*>(a.c_ptr()), const_cast<alglib_impl::ae_vector*>(c.c_ptr()), const_cast<alglib_impl::ae_vector*>(tiesbuf.c_ptr()), const_cast<alglib_impl::ae_vector*>(cntbuf.c_ptr()), const_cast<alglib_impl::ae_vector*>(bufr.c_ptr()), const_cast<alglib_impl::ae_vector*>(bufi.c_ptr()), n, nc, alpha, &info, &threshold, &rms, &cvrms, &_alglib_env_state);
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
_decisionforest_owner::_decisionforest_owner()
{
    p_struct = (alglib_impl::decisionforest*)alglib_impl::ae_malloc(sizeof(alglib_impl::decisionforest), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_decisionforest_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_decisionforest_owner::_decisionforest_owner(const _decisionforest_owner &rhs)
{
    p_struct = (alglib_impl::decisionforest*)alglib_impl::ae_malloc(sizeof(alglib_impl::decisionforest), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_decisionforest_init_copy(p_struct, const_cast<alglib_impl::decisionforest*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_decisionforest_owner& _decisionforest_owner::operator=(const _decisionforest_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_decisionforest_clear(p_struct);
    if( !alglib_impl::_decisionforest_init_copy(p_struct, const_cast<alglib_impl::decisionforest*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_decisionforest_owner::~_decisionforest_owner()
{
    alglib_impl::_decisionforest_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::decisionforest* _decisionforest_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::decisionforest* _decisionforest_owner::c_ptr() const
{
    return const_cast<alglib_impl::decisionforest*>(p_struct);
}
decisionforest::decisionforest() : _decisionforest_owner() 
{
}

decisionforest::decisionforest(const decisionforest &rhs):_decisionforest_owner(rhs) 
{
}

decisionforest& decisionforest::operator=(const decisionforest &rhs)
{
    if( this==&rhs )
        return *this;
    _decisionforest_owner::operator=(rhs);
    return *this;
}

decisionforest::~decisionforest()
{
}


/*************************************************************************

*************************************************************************/
_dfreport_owner::_dfreport_owner()
{
    p_struct = (alglib_impl::dfreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::dfreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_dfreport_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_dfreport_owner::_dfreport_owner(const _dfreport_owner &rhs)
{
    p_struct = (alglib_impl::dfreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::dfreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_dfreport_init_copy(p_struct, const_cast<alglib_impl::dfreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_dfreport_owner& _dfreport_owner::operator=(const _dfreport_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_dfreport_clear(p_struct);
    if( !alglib_impl::_dfreport_init_copy(p_struct, const_cast<alglib_impl::dfreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_dfreport_owner::~_dfreport_owner()
{
    alglib_impl::_dfreport_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::dfreport* _dfreport_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::dfreport* _dfreport_owner::c_ptr() const
{
    return const_cast<alglib_impl::dfreport*>(p_struct);
}
dfreport::dfreport() : _dfreport_owner() ,relclserror(p_struct->relclserror),avgce(p_struct->avgce),rmserror(p_struct->rmserror),avgerror(p_struct->avgerror),avgrelerror(p_struct->avgrelerror),oobrelclserror(p_struct->oobrelclserror),oobavgce(p_struct->oobavgce),oobrmserror(p_struct->oobrmserror),oobavgerror(p_struct->oobavgerror),oobavgrelerror(p_struct->oobavgrelerror)
{
}

dfreport::dfreport(const dfreport &rhs):_dfreport_owner(rhs) ,relclserror(p_struct->relclserror),avgce(p_struct->avgce),rmserror(p_struct->rmserror),avgerror(p_struct->avgerror),avgrelerror(p_struct->avgrelerror),oobrelclserror(p_struct->oobrelclserror),oobavgce(p_struct->oobavgce),oobrmserror(p_struct->oobrmserror),oobavgerror(p_struct->oobavgerror),oobavgrelerror(p_struct->oobavgrelerror)
{
}

dfreport& dfreport::operator=(const dfreport &rhs)
{
    if( this==&rhs )
        return *this;
    _dfreport_owner::operator=(rhs);
    return *this;
}

dfreport::~dfreport()
{
}


/*************************************************************************
This function serializes data structure to string.

Important properties of s_out:
* it contains alphanumeric characters, dots, underscores, minus signs
* these symbols are grouped into words, which are separated by spaces
  and Windows-style (CR+LF) newlines
* although  serializer  uses  spaces and CR+LF as separators, you can 
  replace any separator character by arbitrary combination of spaces,
  tabs, Windows or Unix newlines. It allows flexible reformatting  of
  the  string  in  case you want to include it into text or XML file. 
  But you should not insert separators into the middle of the "words"
  nor you should change case of letters.
* s_out can be freely moved between 32-bit and 64-bit systems, little
  and big endian machines, and so on. You can serialize structure  on
  32-bit machine and unserialize it on 64-bit one (or vice versa), or
  serialize  it  on  SPARC  and  unserialize  on  x86.  You  can also 
  serialize  it  in  C++ version of ALGLIB and unserialize in C# one, 
  and vice versa.
*************************************************************************/
void dfserialize(decisionforest &obj, std::string &s_out)
{
    alglib_impl::ae_state state;
    alglib_impl::ae_serializer serializer;
    alglib_impl::ae_int_t ssize;

    alglib_impl::ae_state_init(&state);
    try
    {
        alglib_impl::ae_serializer_init(&serializer);
        alglib_impl::ae_serializer_alloc_start(&serializer);
        alglib_impl::dfalloc(&serializer, obj.c_ptr(), &state);
        ssize = alglib_impl::ae_serializer_get_alloc_size(&serializer);
        s_out.clear();
        s_out.reserve((size_t)(ssize+1));
        alglib_impl::ae_serializer_sstart_str(&serializer, &s_out);
        alglib_impl::dfserialize(&serializer, obj.c_ptr(), &state);
        alglib_impl::ae_serializer_stop(&serializer);
        if( s_out.length()>(size_t)ssize )
            throw ap_error("ALGLIB: serialization integrity error");
        alglib_impl::ae_serializer_clear(&serializer);
        alglib_impl::ae_state_clear(&state);
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}
/*************************************************************************
This function unserializes data structure from string.
*************************************************************************/
void dfunserialize(std::string &s_in, decisionforest &obj)
{
    alglib_impl::ae_state state;
    alglib_impl::ae_serializer serializer;

    alglib_impl::ae_state_init(&state);
    try
    {
        alglib_impl::ae_serializer_init(&serializer);
        alglib_impl::ae_serializer_ustart_str(&serializer, &s_in);
        alglib_impl::dfunserialize(&serializer, obj.c_ptr(), &state);
        alglib_impl::ae_serializer_stop(&serializer);
        alglib_impl::ae_serializer_clear(&serializer);
        alglib_impl::ae_state_clear(&state);
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
This subroutine builds random decision forest.

INPUT PARAMETERS:
    XY          -   training set
    NPoints     -   training set size, NPoints>=1
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   task type:
                    * NClasses=1 - regression task with one
                                   dependent variable
                    * NClasses>1 - classification task with
                                   NClasses classes.
    NTrees      -   number of trees in a forest, NTrees>=1.
                    recommended values: 50-100.
    R           -   percent of a training set used to build
                    individual trees. 0<R<=1.
                    recommended values: 0.1 <= R <= 0.66.

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<1, NVars<1, NClasses<1, NTrees<1, R<=0
                          or R>1).
                    *  1, if task has been solved
    DF          -   model built
    Rep         -   training report, contains error on a training set
                    and out-of-bag estimates of generalization error.

  -- ALGLIB --
     Copyright 19.02.2009 by Bochkanov Sergey
*************************************************************************/
void dfbuildrandomdecisionforest(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, const ae_int_t ntrees, const double r, ae_int_t &info, decisionforest &df, dfreport &rep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::dfbuildrandomdecisionforest(const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, nvars, nclasses, ntrees, r, &info, const_cast<alglib_impl::decisionforest*>(df.c_ptr()), const_cast<alglib_impl::dfreport*>(rep.c_ptr()), &_alglib_env_state);
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
This subroutine builds random decision forest.
This function gives ability to tune number of variables used when choosing
best split.

INPUT PARAMETERS:
    XY          -   training set
    NPoints     -   training set size, NPoints>=1
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   task type:
                    * NClasses=1 - regression task with one
                                   dependent variable
                    * NClasses>1 - classification task with
                                   NClasses classes.
    NTrees      -   number of trees in a forest, NTrees>=1.
                    recommended values: 50-100.
    NRndVars    -   number of variables used when choosing best split
    R           -   percent of a training set used to build
                    individual trees. 0<R<=1.
                    recommended values: 0.1 <= R <= 0.66.

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<1, NVars<1, NClasses<1, NTrees<1, R<=0
                          or R>1).
                    *  1, if task has been solved
    DF          -   model built
    Rep         -   training report, contains error on a training set
                    and out-of-bag estimates of generalization error.

  -- ALGLIB --
     Copyright 19.02.2009 by Bochkanov Sergey
*************************************************************************/
void dfbuildrandomdecisionforestx1(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, const ae_int_t ntrees, const ae_int_t nrndvars, const double r, ae_int_t &info, decisionforest &df, dfreport &rep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::dfbuildrandomdecisionforestx1(const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, nvars, nclasses, ntrees, nrndvars, r, &info, const_cast<alglib_impl::decisionforest*>(df.c_ptr()), const_cast<alglib_impl::dfreport*>(rep.c_ptr()), &_alglib_env_state);
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
Procesing

INPUT PARAMETERS:
    DF      -   decision forest model
    X       -   input vector,  array[0..NVars-1].

OUTPUT PARAMETERS:
    Y       -   result. Regression estimate when solving regression  task,
                vector of posterior probabilities for classification task.

See also DFProcessI.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
void dfprocess(const decisionforest &df, const real_1d_array &x, real_1d_array &y)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::dfprocess(const_cast<alglib_impl::decisionforest*>(df.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::ae_vector*>(y.c_ptr()), &_alglib_env_state);
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
'interactive' variant of DFProcess for languages like Python which support
constructs like "Y = DFProcessI(DF,X)" and interactive mode of interpreter

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void dfprocessi(const decisionforest &df, const real_1d_array &x, real_1d_array &y)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::dfprocessi(const_cast<alglib_impl::decisionforest*>(df.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::ae_vector*>(y.c_ptr()), &_alglib_env_state);
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
Relative classification error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    percent of incorrectly classified cases.
    Zero if model solves regression task.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfrelclserror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::dfrelclserror(const_cast<alglib_impl::decisionforest*>(df.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
Average cross-entropy (in bits per element) on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    CrossEntropy/(NPoints*LN(2)).
    Zero if model solves regression task.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfavgce(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::dfavgce(const_cast<alglib_impl::decisionforest*>(df.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
RMS error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.
    Its meaning for regression task is obvious. As for
    classification task, RMS error means error when estimating posterior
    probabilities.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfrmserror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::dfrmserror(const_cast<alglib_impl::decisionforest*>(df.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
Average error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for
    classification task, it means average error when estimating posterior
    probabilities.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfavgerror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::dfavgerror(const_cast<alglib_impl::decisionforest*>(df.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
Average relative error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for
    classification task, it means average relative error when estimating
    posterior probability of belonging to the correct class.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfavgrelerror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::dfavgrelerror(const_cast<alglib_impl::decisionforest*>(df.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
_linearmodel_owner::_linearmodel_owner()
{
    p_struct = (alglib_impl::linearmodel*)alglib_impl::ae_malloc(sizeof(alglib_impl::linearmodel), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_linearmodel_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_linearmodel_owner::_linearmodel_owner(const _linearmodel_owner &rhs)
{
    p_struct = (alglib_impl::linearmodel*)alglib_impl::ae_malloc(sizeof(alglib_impl::linearmodel), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_linearmodel_init_copy(p_struct, const_cast<alglib_impl::linearmodel*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_linearmodel_owner& _linearmodel_owner::operator=(const _linearmodel_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_linearmodel_clear(p_struct);
    if( !alglib_impl::_linearmodel_init_copy(p_struct, const_cast<alglib_impl::linearmodel*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_linearmodel_owner::~_linearmodel_owner()
{
    alglib_impl::_linearmodel_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::linearmodel* _linearmodel_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::linearmodel* _linearmodel_owner::c_ptr() const
{
    return const_cast<alglib_impl::linearmodel*>(p_struct);
}
linearmodel::linearmodel() : _linearmodel_owner() 
{
}

linearmodel::linearmodel(const linearmodel &rhs):_linearmodel_owner(rhs) 
{
}

linearmodel& linearmodel::operator=(const linearmodel &rhs)
{
    if( this==&rhs )
        return *this;
    _linearmodel_owner::operator=(rhs);
    return *this;
}

linearmodel::~linearmodel()
{
}


/*************************************************************************
LRReport structure contains additional information about linear model:
* C             -   covariation matrix,  array[0..NVars,0..NVars].
                    C[i,j] = Cov(A[i],A[j])
* RMSError      -   root mean square error on a training set
* AvgError      -   average error on a training set
* AvgRelError   -   average relative error on a training set (excluding
                    observations with zero function value).
* CVRMSError    -   leave-one-out cross-validation estimate of
                    generalization error. Calculated using fast algorithm
                    with O(NVars*NPoints) complexity.
* CVAvgError    -   cross-validation estimate of average error
* CVAvgRelError -   cross-validation estimate of average relative error

All other fields of the structure are intended for internal use and should
not be used outside ALGLIB.
*************************************************************************/
_lrreport_owner::_lrreport_owner()
{
    p_struct = (alglib_impl::lrreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::lrreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_lrreport_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_lrreport_owner::_lrreport_owner(const _lrreport_owner &rhs)
{
    p_struct = (alglib_impl::lrreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::lrreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_lrreport_init_copy(p_struct, const_cast<alglib_impl::lrreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_lrreport_owner& _lrreport_owner::operator=(const _lrreport_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_lrreport_clear(p_struct);
    if( !alglib_impl::_lrreport_init_copy(p_struct, const_cast<alglib_impl::lrreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_lrreport_owner::~_lrreport_owner()
{
    alglib_impl::_lrreport_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::lrreport* _lrreport_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::lrreport* _lrreport_owner::c_ptr() const
{
    return const_cast<alglib_impl::lrreport*>(p_struct);
}
lrreport::lrreport() : _lrreport_owner() ,c(&p_struct->c),rmserror(p_struct->rmserror),avgerror(p_struct->avgerror),avgrelerror(p_struct->avgrelerror),cvrmserror(p_struct->cvrmserror),cvavgerror(p_struct->cvavgerror),cvavgrelerror(p_struct->cvavgrelerror),ncvdefects(p_struct->ncvdefects),cvdefects(&p_struct->cvdefects)
{
}

lrreport::lrreport(const lrreport &rhs):_lrreport_owner(rhs) ,c(&p_struct->c),rmserror(p_struct->rmserror),avgerror(p_struct->avgerror),avgrelerror(p_struct->avgrelerror),cvrmserror(p_struct->cvrmserror),cvavgerror(p_struct->cvavgerror),cvavgrelerror(p_struct->cvavgrelerror),ncvdefects(p_struct->ncvdefects),cvdefects(&p_struct->cvdefects)
{
}

lrreport& lrreport::operator=(const lrreport &rhs)
{
    if( this==&rhs )
        return *this;
    _lrreport_owner::operator=(rhs);
    return *this;
}

lrreport::~lrreport()
{
}

/*************************************************************************
Linear regression

Subroutine builds model:

    Y = A(0)*X[0] + ... + A(N-1)*X[N-1] + A(N)

and model found in ALGLIB format, covariation matrix, training set  errors
(rms,  average,  average  relative)   and  leave-one-out  cross-validation
estimate of the generalization error. CV  estimate calculated  using  fast
algorithm with O(NPoints*NVars) complexity.

When  covariation  matrix  is  calculated  standard deviations of function
values are assumed to be equal to RMS error on the training set.

INPUT PARAMETERS:
    XY          -   training set, array [0..NPoints-1,0..NVars]:
                    * NVars columns - independent variables
                    * last column - dependent variable
    NPoints     -   training set size, NPoints>NVars+1
    NVars       -   number of independent variables

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -255, in case of unknown internal error
                    * -4, if internal SVD subroutine haven't converged
                    * -1, if incorrect parameters was passed (NPoints<NVars+2, NVars<1).
                    *  1, if subroutine successfully finished
    LM          -   linear model in the ALGLIB format. Use subroutines of
                    this unit to work with the model.
    AR          -   additional results


  -- ALGLIB --
     Copyright 02.08.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuild(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::lrbuild(const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, nvars, &info, const_cast<alglib_impl::linearmodel*>(lm.c_ptr()), const_cast<alglib_impl::lrreport*>(ar.c_ptr()), &_alglib_env_state);
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
Linear regression

Variant of LRBuild which uses vector of standatd deviations (errors in
function values).

INPUT PARAMETERS:
    XY          -   training set, array [0..NPoints-1,0..NVars]:
                    * NVars columns - independent variables
                    * last column - dependent variable
    S           -   standard deviations (errors in function values)
                    array[0..NPoints-1], S[i]>0.
    NPoints     -   training set size, NPoints>NVars+1
    NVars       -   number of independent variables

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -255, in case of unknown internal error
                    * -4, if internal SVD subroutine haven't converged
                    * -1, if incorrect parameters was passed (NPoints<NVars+2, NVars<1).
                    * -2, if S[I]<=0
                    *  1, if subroutine successfully finished
    LM          -   linear model in the ALGLIB format. Use subroutines of
                    this unit to work with the model.
    AR          -   additional results


  -- ALGLIB --
     Copyright 02.08.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuilds(const real_2d_array &xy, const real_1d_array &s, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::lrbuilds(const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), const_cast<alglib_impl::ae_vector*>(s.c_ptr()), npoints, nvars, &info, const_cast<alglib_impl::linearmodel*>(lm.c_ptr()), const_cast<alglib_impl::lrreport*>(ar.c_ptr()), &_alglib_env_state);
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
Like LRBuildS, but builds model

    Y = A(0)*X[0] + ... + A(N-1)*X[N-1]

i.e. with zero constant term.

  -- ALGLIB --
     Copyright 30.10.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuildzs(const real_2d_array &xy, const real_1d_array &s, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::lrbuildzs(const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), const_cast<alglib_impl::ae_vector*>(s.c_ptr()), npoints, nvars, &info, const_cast<alglib_impl::linearmodel*>(lm.c_ptr()), const_cast<alglib_impl::lrreport*>(ar.c_ptr()), &_alglib_env_state);
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
Like LRBuild but builds model

    Y = A(0)*X[0] + ... + A(N-1)*X[N-1]

i.e. with zero constant term.

  -- ALGLIB --
     Copyright 30.10.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuildz(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::lrbuildz(const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, nvars, &info, const_cast<alglib_impl::linearmodel*>(lm.c_ptr()), const_cast<alglib_impl::lrreport*>(ar.c_ptr()), &_alglib_env_state);
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
Unpacks coefficients of linear model.

INPUT PARAMETERS:
    LM          -   linear model in ALGLIB format

OUTPUT PARAMETERS:
    V           -   coefficients, array[0..NVars]
                    constant term (intercept) is stored in the V[NVars].
    NVars       -   number of independent variables (one less than number
                    of coefficients)

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
void lrunpack(const linearmodel &lm, real_1d_array &v, ae_int_t &nvars)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::lrunpack(const_cast<alglib_impl::linearmodel*>(lm.c_ptr()), const_cast<alglib_impl::ae_vector*>(v.c_ptr()), &nvars, &_alglib_env_state);
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
"Packs" coefficients and creates linear model in ALGLIB format (LRUnpack
reversed).

INPUT PARAMETERS:
    V           -   coefficients, array[0..NVars]
    NVars       -   number of independent variables

OUTPUT PAREMETERS:
    LM          -   linear model.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
void lrpack(const real_1d_array &v, const ae_int_t nvars, linearmodel &lm)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::lrpack(const_cast<alglib_impl::ae_vector*>(v.c_ptr()), nvars, const_cast<alglib_impl::linearmodel*>(lm.c_ptr()), &_alglib_env_state);
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
Procesing

INPUT PARAMETERS:
    LM      -   linear model
    X       -   input vector,  array[0..NVars-1].

Result:
    value of linear model regression estimate

  -- ALGLIB --
     Copyright 03.09.2008 by Bochkanov Sergey
*************************************************************************/
double lrprocess(const linearmodel &lm, const real_1d_array &x)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::lrprocess(const_cast<alglib_impl::linearmodel*>(lm.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
RMS error on the test set

INPUT PARAMETERS:
    LM      -   linear model
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double lrrmserror(const linearmodel &lm, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::lrrmserror(const_cast<alglib_impl::linearmodel*>(lm.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
Average error on the test set

INPUT PARAMETERS:
    LM      -   linear model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average error.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double lravgerror(const linearmodel &lm, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::lravgerror(const_cast<alglib_impl::linearmodel*>(lm.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
RMS error on the test set

INPUT PARAMETERS:
    LM      -   linear model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average relative error.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double lravgrelerror(const linearmodel &lm, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::lravgrelerror(const_cast<alglib_impl::linearmodel*>(lm.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
Filters: simple moving averages (unsymmetric).

This filter replaces array by results of SMA(K) filter. SMA(K) is defined
as filter which averages at most K previous points (previous - not points
AROUND central point) - or less, in case of the first K-1 points.

INPUT PARAMETERS:
    X           -   array[N], array to process. It can be larger than N,
                    in this case only first N points are processed.
    N           -   points count, N>=0
    K           -   K>=1 (K can be larger than N ,  such  cases  will  be
                    correctly handled). Window width. K=1 corresponds  to
                    identity transformation (nothing changes).

OUTPUT PARAMETERS:
    X           -   array, whose first N elements were processed with SMA(K)

NOTE 1: this function uses efficient in-place  algorithm  which  does not
        allocate temporary arrays.

NOTE 2: this algorithm makes only one pass through array and uses running
        sum  to speed-up calculation of the averages. Additional measures
        are taken to ensure that running sum on a long sequence  of  zero
        elements will be correctly reset to zero even in the presence  of
        round-off error.

NOTE 3: this  is  unsymmetric version of the algorithm,  which  does  NOT
        averages points after the current one. Only X[i], X[i-1], ... are
        used when calculating new value of X[i]. We should also note that
        this algorithm uses BOTH previous points and  current  one,  i.e.
        new value of X[i] depends on BOTH previous point and X[i] itself.

  -- ALGLIB --
     Copyright 25.10.2011 by Bochkanov Sergey
*************************************************************************/
void filtersma(real_1d_array &x, const ae_int_t n, const ae_int_t k)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::filtersma(const_cast<alglib_impl::ae_vector*>(x.c_ptr()), n, k, &_alglib_env_state);
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
Filters: simple moving averages (unsymmetric).

This filter replaces array by results of SMA(K) filter. SMA(K) is defined
as filter which averages at most K previous points (previous - not points
AROUND central point) - or less, in case of the first K-1 points.

INPUT PARAMETERS:
    X           -   array[N], array to process. It can be larger than N,
                    in this case only first N points are processed.
    N           -   points count, N>=0
    K           -   K>=1 (K can be larger than N ,  such  cases  will  be
                    correctly handled). Window width. K=1 corresponds  to
                    identity transformation (nothing changes).

OUTPUT PARAMETERS:
    X           -   array, whose first N elements were processed with SMA(K)

NOTE 1: this function uses efficient in-place  algorithm  which  does not
        allocate temporary arrays.

NOTE 2: this algorithm makes only one pass through array and uses running
        sum  to speed-up calculation of the averages. Additional measures
        are taken to ensure that running sum on a long sequence  of  zero
        elements will be correctly reset to zero even in the presence  of
        round-off error.

NOTE 3: this  is  unsymmetric version of the algorithm,  which  does  NOT
        averages points after the current one. Only X[i], X[i-1], ... are
        used when calculating new value of X[i]. We should also note that
        this algorithm uses BOTH previous points and  current  one,  i.e.
        new value of X[i] depends on BOTH previous point and X[i] itself.

  -- ALGLIB --
     Copyright 25.10.2011 by Bochkanov Sergey
*************************************************************************/
void filtersma(real_1d_array &x, const ae_int_t k)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;

    n = x.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::filtersma(const_cast<alglib_impl::ae_vector*>(x.c_ptr()), n, k, &_alglib_env_state);

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
Filters: exponential moving averages.

This filter replaces array by results of EMA(alpha) filter. EMA(alpha) is
defined as filter which replaces X[] by S[]:
    S[0] = X[0]
    S[t] = alpha*X[t] + (1-alpha)*S[t-1]

INPUT PARAMETERS:
    X           -   array[N], array to process. It can be larger than N,
                    in this case only first N points are processed.
    N           -   points count, N>=0
    alpha       -   0<alpha<=1, smoothing parameter.

OUTPUT PARAMETERS:
    X           -   array, whose first N elements were processed
                    with EMA(alpha)

NOTE 1: this function uses efficient in-place  algorithm  which  does not
        allocate temporary arrays.

NOTE 2: this algorithm uses BOTH previous points and  current  one,  i.e.
        new value of X[i] depends on BOTH previous point and X[i] itself.

NOTE 3: technical analytis users quite often work  with  EMA  coefficient
        expressed in DAYS instead of fractions. If you want to  calculate
        EMA(N), where N is a number of days, you can use alpha=2/(N+1).

  -- ALGLIB --
     Copyright 25.10.2011 by Bochkanov Sergey
*************************************************************************/
void filterema(real_1d_array &x, const ae_int_t n, const double alpha)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::filterema(const_cast<alglib_impl::ae_vector*>(x.c_ptr()), n, alpha, &_alglib_env_state);
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
Filters: exponential moving averages.

This filter replaces array by results of EMA(alpha) filter. EMA(alpha) is
defined as filter which replaces X[] by S[]:
    S[0] = X[0]
    S[t] = alpha*X[t] + (1-alpha)*S[t-1]

INPUT PARAMETERS:
    X           -   array[N], array to process. It can be larger than N,
                    in this case only first N points are processed.
    N           -   points count, N>=0
    alpha       -   0<alpha<=1, smoothing parameter.

OUTPUT PARAMETERS:
    X           -   array, whose first N elements were processed
                    with EMA(alpha)

NOTE 1: this function uses efficient in-place  algorithm  which  does not
        allocate temporary arrays.

NOTE 2: this algorithm uses BOTH previous points and  current  one,  i.e.
        new value of X[i] depends on BOTH previous point and X[i] itself.

NOTE 3: technical analytis users quite often work  with  EMA  coefficient
        expressed in DAYS instead of fractions. If you want to  calculate
        EMA(N), where N is a number of days, you can use alpha=2/(N+1).

  -- ALGLIB --
     Copyright 25.10.2011 by Bochkanov Sergey
*************************************************************************/
void filterema(real_1d_array &x, const double alpha)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;

    n = x.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::filterema(const_cast<alglib_impl::ae_vector*>(x.c_ptr()), n, alpha, &_alglib_env_state);

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
Filters: linear regression moving averages.

This filter replaces array by results of LRMA(K) filter.

LRMA(K) is defined as filter which, for each data  point,  builds  linear
regression  model  using  K  prevous  points (point itself is included in
these K points) and calculates value of this linear model at the point in
question.

INPUT PARAMETERS:
    X           -   array[N], array to process. It can be larger than N,
                    in this case only first N points are processed.
    N           -   points count, N>=0
    K           -   K>=1 (K can be larger than N ,  such  cases  will  be
                    correctly handled). Window width. K=1 corresponds  to
                    identity transformation (nothing changes).

OUTPUT PARAMETERS:
    X           -   array, whose first N elements were processed with SMA(K)

NOTE 1: this function uses efficient in-place  algorithm  which  does not
        allocate temporary arrays.

NOTE 2: this algorithm makes only one pass through array and uses running
        sum  to speed-up calculation of the averages. Additional measures
        are taken to ensure that running sum on a long sequence  of  zero
        elements will be correctly reset to zero even in the presence  of
        round-off error.

NOTE 3: this  is  unsymmetric version of the algorithm,  which  does  NOT
        averages points after the current one. Only X[i], X[i-1], ... are
        used when calculating new value of X[i]. We should also note that
        this algorithm uses BOTH previous points and  current  one,  i.e.
        new value of X[i] depends on BOTH previous point and X[i] itself.

  -- ALGLIB --
     Copyright 25.10.2011 by Bochkanov Sergey
*************************************************************************/
void filterlrma(real_1d_array &x, const ae_int_t n, const ae_int_t k)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::filterlrma(const_cast<alglib_impl::ae_vector*>(x.c_ptr()), n, k, &_alglib_env_state);
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
Filters: linear regression moving averages.

This filter replaces array by results of LRMA(K) filter.

LRMA(K) is defined as filter which, for each data  point,  builds  linear
regression  model  using  K  prevous  points (point itself is included in
these K points) and calculates value of this linear model at the point in
question.

INPUT PARAMETERS:
    X           -   array[N], array to process. It can be larger than N,
                    in this case only first N points are processed.
    N           -   points count, N>=0
    K           -   K>=1 (K can be larger than N ,  such  cases  will  be
                    correctly handled). Window width. K=1 corresponds  to
                    identity transformation (nothing changes).

OUTPUT PARAMETERS:
    X           -   array, whose first N elements were processed with SMA(K)

NOTE 1: this function uses efficient in-place  algorithm  which  does not
        allocate temporary arrays.

NOTE 2: this algorithm makes only one pass through array and uses running
        sum  to speed-up calculation of the averages. Additional measures
        are taken to ensure that running sum on a long sequence  of  zero
        elements will be correctly reset to zero even in the presence  of
        round-off error.

NOTE 3: this  is  unsymmetric version of the algorithm,  which  does  NOT
        averages points after the current one. Only X[i], X[i-1], ... are
        used when calculating new value of X[i]. We should also note that
        this algorithm uses BOTH previous points and  current  one,  i.e.
        new value of X[i] depends on BOTH previous point and X[i] itself.

  -- ALGLIB --
     Copyright 25.10.2011 by Bochkanov Sergey
*************************************************************************/
void filterlrma(real_1d_array &x, const ae_int_t k)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t n;

    n = x.length();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::filterlrma(const_cast<alglib_impl::ae_vector*>(x.c_ptr()), n, k, &_alglib_env_state);

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
k-means++ clusterization

INPUT PARAMETERS:
    XY          -   dataset, array [0..NPoints-1,0..NVars-1].
    NPoints     -   dataset size, NPoints>=K
    NVars       -   number of variables, NVars>=1
    K           -   desired number of clusters, K>=1
    Restarts    -   number of restarts, Restarts>=1

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -3, if task is degenerate (number of distinct points is
                          less than K)
                    * -1, if incorrect NPoints/NFeatures/K/Restarts was passed
                    *  1, if subroutine finished successfully
    C           -   array[0..NVars-1,0..K-1].matrix whose columns store
                    cluster's centers
    XYC         -   array[NPoints], which contains cluster indexes

  -- ALGLIB --
     Copyright 21.03.2009 by Bochkanov Sergey
*************************************************************************/
void kmeansgenerate(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t k, const ae_int_t restarts, ae_int_t &info, real_2d_array &c, integer_1d_array &xyc)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::kmeansgenerate(const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, nvars, k, restarts, &info, const_cast<alglib_impl::ae_matrix*>(c.c_ptr()), const_cast<alglib_impl::ae_vector*>(xyc.c_ptr()), &_alglib_env_state);
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
Multiclass Fisher LDA

Subroutine finds coefficients of linear combination which optimally separates
training set on classes.

INPUT PARAMETERS:
    XY          -   training set, array[0..NPoints-1,0..NVars].
                    First NVars columns store values of independent
                    variables, next column stores number of class (from 0
                    to NClasses-1) which dataset element belongs to. Fractional
                    values are rounded to nearest integer.
    NPoints     -   training set size, NPoints>=0
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   number of classes, NClasses>=2


OUTPUT PARAMETERS:
    Info        -   return code:
                    * -4, if internal EVD subroutine hasn't converged
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed (NPoints<0,
                          NVars<1, NClasses<2)
                    *  1, if task has been solved
                    *  2, if there was a multicollinearity in training set,
                          but task has been solved.
    W           -   linear combination coefficients, array[0..NVars-1]

  -- ALGLIB --
     Copyright 31.05.2008 by Bochkanov Sergey
*************************************************************************/
void fisherlda(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, real_1d_array &w)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::fisherlda(const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, nvars, nclasses, &info, const_cast<alglib_impl::ae_vector*>(w.c_ptr()), &_alglib_env_state);
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
N-dimensional multiclass Fisher LDA

Subroutine finds coefficients of linear combinations which optimally separates
training set on classes. It returns N-dimensional basis whose vector are sorted
by quality of training set separation (in descending order).

INPUT PARAMETERS:
    XY          -   training set, array[0..NPoints-1,0..NVars].
                    First NVars columns store values of independent
                    variables, next column stores number of class (from 0
                    to NClasses-1) which dataset element belongs to. Fractional
                    values are rounded to nearest integer.
    NPoints     -   training set size, NPoints>=0
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   number of classes, NClasses>=2


OUTPUT PARAMETERS:
    Info        -   return code:
                    * -4, if internal EVD subroutine hasn't converged
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed (NPoints<0,
                          NVars<1, NClasses<2)
                    *  1, if task has been solved
                    *  2, if there was a multicollinearity in training set,
                          but task has been solved.
    W           -   basis, array[0..NVars-1,0..NVars-1]
                    columns of matrix stores basis vectors, sorted by
                    quality of training set separation (in descending order)

  -- ALGLIB --
     Copyright 31.05.2008 by Bochkanov Sergey
*************************************************************************/
void fisherldan(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, real_2d_array &w)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::fisherldan(const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, nvars, nclasses, &info, const_cast<alglib_impl::ae_matrix*>(w.c_ptr()), &_alglib_env_state);
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
_multilayerperceptron_owner::_multilayerperceptron_owner()
{
    p_struct = (alglib_impl::multilayerperceptron*)alglib_impl::ae_malloc(sizeof(alglib_impl::multilayerperceptron), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_multilayerperceptron_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_multilayerperceptron_owner::_multilayerperceptron_owner(const _multilayerperceptron_owner &rhs)
{
    p_struct = (alglib_impl::multilayerperceptron*)alglib_impl::ae_malloc(sizeof(alglib_impl::multilayerperceptron), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_multilayerperceptron_init_copy(p_struct, const_cast<alglib_impl::multilayerperceptron*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_multilayerperceptron_owner& _multilayerperceptron_owner::operator=(const _multilayerperceptron_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_multilayerperceptron_clear(p_struct);
    if( !alglib_impl::_multilayerperceptron_init_copy(p_struct, const_cast<alglib_impl::multilayerperceptron*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_multilayerperceptron_owner::~_multilayerperceptron_owner()
{
    alglib_impl::_multilayerperceptron_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::multilayerperceptron* _multilayerperceptron_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::multilayerperceptron* _multilayerperceptron_owner::c_ptr() const
{
    return const_cast<alglib_impl::multilayerperceptron*>(p_struct);
}
multilayerperceptron::multilayerperceptron() : _multilayerperceptron_owner() 
{
}

multilayerperceptron::multilayerperceptron(const multilayerperceptron &rhs):_multilayerperceptron_owner(rhs) 
{
}

multilayerperceptron& multilayerperceptron::operator=(const multilayerperceptron &rhs)
{
    if( this==&rhs )
        return *this;
    _multilayerperceptron_owner::operator=(rhs);
    return *this;
}

multilayerperceptron::~multilayerperceptron()
{
}


/*************************************************************************
This function serializes data structure to string.

Important properties of s_out:
* it contains alphanumeric characters, dots, underscores, minus signs
* these symbols are grouped into words, which are separated by spaces
  and Windows-style (CR+LF) newlines
* although  serializer  uses  spaces and CR+LF as separators, you can 
  replace any separator character by arbitrary combination of spaces,
  tabs, Windows or Unix newlines. It allows flexible reformatting  of
  the  string  in  case you want to include it into text or XML file. 
  But you should not insert separators into the middle of the "words"
  nor you should change case of letters.
* s_out can be freely moved between 32-bit and 64-bit systems, little
  and big endian machines, and so on. You can serialize structure  on
  32-bit machine and unserialize it on 64-bit one (or vice versa), or
  serialize  it  on  SPARC  and  unserialize  on  x86.  You  can also 
  serialize  it  in  C++ version of ALGLIB and unserialize in C# one, 
  and vice versa.
*************************************************************************/
void mlpserialize(multilayerperceptron &obj, std::string &s_out)
{
    alglib_impl::ae_state state;
    alglib_impl::ae_serializer serializer;
    alglib_impl::ae_int_t ssize;

    alglib_impl::ae_state_init(&state);
    try
    {
        alglib_impl::ae_serializer_init(&serializer);
        alglib_impl::ae_serializer_alloc_start(&serializer);
        alglib_impl::mlpalloc(&serializer, obj.c_ptr(), &state);
        ssize = alglib_impl::ae_serializer_get_alloc_size(&serializer);
        s_out.clear();
        s_out.reserve((size_t)(ssize+1));
        alglib_impl::ae_serializer_sstart_str(&serializer, &s_out);
        alglib_impl::mlpserialize(&serializer, obj.c_ptr(), &state);
        alglib_impl::ae_serializer_stop(&serializer);
        if( s_out.length()>(size_t)ssize )
            throw ap_error("ALGLIB: serialization integrity error");
        alglib_impl::ae_serializer_clear(&serializer);
        alglib_impl::ae_state_clear(&state);
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}
/*************************************************************************
This function unserializes data structure from string.
*************************************************************************/
void mlpunserialize(std::string &s_in, multilayerperceptron &obj)
{
    alglib_impl::ae_state state;
    alglib_impl::ae_serializer serializer;

    alglib_impl::ae_state_init(&state);
    try
    {
        alglib_impl::ae_serializer_init(&serializer);
        alglib_impl::ae_serializer_ustart_str(&serializer, &s_in);
        alglib_impl::mlpunserialize(&serializer, obj.c_ptr(), &state);
        alglib_impl::ae_serializer_stop(&serializer);
        alglib_impl::ae_serializer_clear(&serializer);
        alglib_impl::ae_state_clear(&state);
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Creates  neural  network  with  NIn  inputs,  NOut outputs, without hidden
layers, with linear output layer. Network weights are  filled  with  small
random values.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreate0(const ae_int_t nin, const ae_int_t nout, multilayerperceptron &network)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpcreate0(nin, nout, const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), &_alglib_env_state);
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
Same  as  MLPCreate0,  but  with  one  hidden  layer  (NHid  neurons) with
non-linear activation function. Output layer is linear.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreate1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, multilayerperceptron &network)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpcreate1(nin, nhid, nout, const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), &_alglib_env_state);
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
Same as MLPCreate0, but with two hidden layers (NHid1 and  NHid2  neurons)
with non-linear activation function. Output layer is linear.
 $ALL

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreate2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, multilayerperceptron &network)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpcreate2(nin, nhid1, nhid2, nout, const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), &_alglib_env_state);
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
Creates  neural  network  with  NIn  inputs,  NOut outputs, without hidden
layers with non-linear output layer. Network weights are filled with small
random values.

Activation function of the output layer takes values:

    (B, +INF), if D>=0

or

    (-INF, B), if D<0.


  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreateb0(const ae_int_t nin, const ae_int_t nout, const double b, const double d, multilayerperceptron &network)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpcreateb0(nin, nout, b, d, const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), &_alglib_env_state);
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
Same as MLPCreateB0 but with non-linear hidden layer.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreateb1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double b, const double d, multilayerperceptron &network)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpcreateb1(nin, nhid, nout, b, d, const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), &_alglib_env_state);
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
Same as MLPCreateB0 but with two non-linear hidden layers.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreateb2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double b, const double d, multilayerperceptron &network)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpcreateb2(nin, nhid1, nhid2, nout, b, d, const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), &_alglib_env_state);
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
Creates  neural  network  with  NIn  inputs,  NOut outputs, without hidden
layers with non-linear output layer. Network weights are filled with small
random values. Activation function of the output layer takes values [A,B].

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreater0(const ae_int_t nin, const ae_int_t nout, const double a, const double b, multilayerperceptron &network)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpcreater0(nin, nout, a, b, const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), &_alglib_env_state);
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
Same as MLPCreateR0, but with non-linear hidden layer.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreater1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double a, const double b, multilayerperceptron &network)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpcreater1(nin, nhid, nout, a, b, const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), &_alglib_env_state);
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
Same as MLPCreateR0, but with two non-linear hidden layers.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreater2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double a, const double b, multilayerperceptron &network)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpcreater2(nin, nhid1, nhid2, nout, a, b, const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), &_alglib_env_state);
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
Creates classifier network with NIn  inputs  and  NOut  possible  classes.
Network contains no hidden layers and linear output  layer  with  SOFTMAX-
normalization  (so  outputs  sums  up  to  1.0  and  converge to posterior
probabilities).

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreatec0(const ae_int_t nin, const ae_int_t nout, multilayerperceptron &network)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpcreatec0(nin, nout, const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), &_alglib_env_state);
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
Same as MLPCreateC0, but with one non-linear hidden layer.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreatec1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, multilayerperceptron &network)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpcreatec1(nin, nhid, nout, const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), &_alglib_env_state);
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
Same as MLPCreateC0, but with two non-linear hidden layers.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreatec2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, multilayerperceptron &network)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpcreatec2(nin, nhid1, nhid2, nout, const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), &_alglib_env_state);
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
Randomization of neural network weights

  -- ALGLIB --
     Copyright 06.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlprandomize(const multilayerperceptron &network)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlprandomize(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), &_alglib_env_state);
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
Randomization of neural network weights and standartisator

  -- ALGLIB --
     Copyright 10.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlprandomizefull(const multilayerperceptron &network)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlprandomizefull(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), &_alglib_env_state);
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
Returns information about initialized network: number of inputs, outputs,
weights.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpproperties(const multilayerperceptron &network, ae_int_t &nin, ae_int_t &nout, ae_int_t &wcount)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpproperties(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), &nin, &nout, &wcount, &_alglib_env_state);
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
Returns number of inputs.

  -- ALGLIB --
     Copyright 19.10.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpgetinputscount(const multilayerperceptron &network)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::ae_int_t result = alglib_impl::mlpgetinputscount(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<ae_int_t*>(&result));
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
Returns number of outputs.

  -- ALGLIB --
     Copyright 19.10.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpgetoutputscount(const multilayerperceptron &network)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::ae_int_t result = alglib_impl::mlpgetoutputscount(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<ae_int_t*>(&result));
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
Returns number of weights.

  -- ALGLIB --
     Copyright 19.10.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpgetweightscount(const multilayerperceptron &network)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::ae_int_t result = alglib_impl::mlpgetweightscount(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<ae_int_t*>(&result));
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
Tells whether network is SOFTMAX-normalized (i.e. classifier) or not.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
bool mlpissoftmax(const multilayerperceptron &network)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        ae_bool result = alglib_impl::mlpissoftmax(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), &_alglib_env_state);
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

/*************************************************************************
This function returns total number of layers (including input, hidden and
output layers).

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpgetlayerscount(const multilayerperceptron &network)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::ae_int_t result = alglib_impl::mlpgetlayerscount(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<ae_int_t*>(&result));
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
This function returns size of K-th layer.

K=0 corresponds to input layer, K=CNT-1 corresponds to output layer.

Size of the output layer is always equal to the number of outputs, although
when we have softmax-normalized network, last neuron doesn't have any
connections - it is just zero.

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpgetlayersize(const multilayerperceptron &network, const ae_int_t k)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::ae_int_t result = alglib_impl::mlpgetlayersize(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), k, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<ae_int_t*>(&result));
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
This function returns offset/scaling coefficients for I-th input of the
network.

INPUT PARAMETERS:
    Network     -   network
    I           -   input index

OUTPUT PARAMETERS:
    Mean        -   mean term
    Sigma       -   sigma term, guaranteed to be nonzero.

I-th input is passed through linear transformation
    IN[i] = (IN[i]-Mean)/Sigma
before feeding to the network

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpgetinputscaling(const multilayerperceptron &network, const ae_int_t i, double &mean, double &sigma)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpgetinputscaling(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), i, &mean, &sigma, &_alglib_env_state);
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
This function returns offset/scaling coefficients for I-th output of the
network.

INPUT PARAMETERS:
    Network     -   network
    I           -   input index

OUTPUT PARAMETERS:
    Mean        -   mean term
    Sigma       -   sigma term, guaranteed to be nonzero.

I-th output is passed through linear transformation
    OUT[i] = OUT[i]*Sigma+Mean
before returning it to user. In case we have SOFTMAX-normalized network,
we return (Mean,Sigma)=(0.0,1.0).

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpgetoutputscaling(const multilayerperceptron &network, const ae_int_t i, double &mean, double &sigma)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpgetoutputscaling(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), i, &mean, &sigma, &_alglib_env_state);
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
This function returns information about Ith neuron of Kth layer

INPUT PARAMETERS:
    Network     -   network
    K           -   layer index
    I           -   neuron index (within layer)

OUTPUT PARAMETERS:
    FKind       -   activation function type (used by MLPActivationFunction())
                    this value is zero for input or linear neurons
    Threshold   -   also called offset, bias
                    zero for input neurons

NOTE: this function throws exception if layer or neuron with  given  index
do not exists.

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpgetneuroninfo(const multilayerperceptron &network, const ae_int_t k, const ae_int_t i, ae_int_t &fkind, double &threshold)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpgetneuroninfo(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), k, i, &fkind, &threshold, &_alglib_env_state);
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
This function returns information about connection from I0-th neuron of
K0-th layer to I1-th neuron of K1-th layer.

INPUT PARAMETERS:
    Network     -   network
    K0          -   layer index
    I0          -   neuron index (within layer)
    K1          -   layer index
    I1          -   neuron index (within layer)

RESULT:
    connection weight (zero for non-existent connections)

This function:
1. throws exception if layer or neuron with given index do not exists.
2. returns zero if neurons exist, but there is no connection between them

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
double mlpgetweight(const multilayerperceptron &network, const ae_int_t k0, const ae_int_t i0, const ae_int_t k1, const ae_int_t i1)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::mlpgetweight(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), k0, i0, k1, i1, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
This function sets offset/scaling coefficients for I-th input of the
network.

INPUT PARAMETERS:
    Network     -   network
    I           -   input index
    Mean        -   mean term
    Sigma       -   sigma term (if zero, will be replaced by 1.0)

NTE: I-th input is passed through linear transformation
    IN[i] = (IN[i]-Mean)/Sigma
before feeding to the network. This function sets Mean and Sigma.

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpsetinputscaling(const multilayerperceptron &network, const ae_int_t i, const double mean, const double sigma)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpsetinputscaling(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), i, mean, sigma, &_alglib_env_state);
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
This function sets offset/scaling coefficients for I-th output of the
network.

INPUT PARAMETERS:
    Network     -   network
    I           -   input index
    Mean        -   mean term
    Sigma       -   sigma term (if zero, will be replaced by 1.0)

OUTPUT PARAMETERS:

NOTE: I-th output is passed through linear transformation
    OUT[i] = OUT[i]*Sigma+Mean
before returning it to user. This function sets Sigma/Mean. In case we
have SOFTMAX-normalized network, you can not set (Sigma,Mean) to anything
other than(0.0,1.0) - this function will throw exception.

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpsetoutputscaling(const multilayerperceptron &network, const ae_int_t i, const double mean, const double sigma)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpsetoutputscaling(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), i, mean, sigma, &_alglib_env_state);
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
This function modifies information about Ith neuron of Kth layer

INPUT PARAMETERS:
    Network     -   network
    K           -   layer index
    I           -   neuron index (within layer)
    FKind       -   activation function type (used by MLPActivationFunction())
                    this value must be zero for input neurons
                    (you can not set activation function for input neurons)
    Threshold   -   also called offset, bias
                    this value must be zero for input neurons
                    (you can not set threshold for input neurons)

NOTES:
1. this function throws exception if layer or neuron with given index do
   not exists.
2. this function also throws exception when you try to set non-linear
   activation function for input neurons (any kind of network) or for output
   neurons of classifier network.
3. this function throws exception when you try to set non-zero threshold for
   input neurons (any kind of network).

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpsetneuroninfo(const multilayerperceptron &network, const ae_int_t k, const ae_int_t i, const ae_int_t fkind, const double threshold)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpsetneuroninfo(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), k, i, fkind, threshold, &_alglib_env_state);
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
This function modifies information about connection from I0-th neuron of
K0-th layer to I1-th neuron of K1-th layer.

INPUT PARAMETERS:
    Network     -   network
    K0          -   layer index
    I0          -   neuron index (within layer)
    K1          -   layer index
    I1          -   neuron index (within layer)
    W           -   connection weight (must be zero for non-existent
                    connections)

This function:
1. throws exception if layer or neuron with given index do not exists.
2. throws exception if you try to set non-zero weight for non-existent
   connection

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpsetweight(const multilayerperceptron &network, const ae_int_t k0, const ae_int_t i0, const ae_int_t k1, const ae_int_t i1, const double w)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpsetweight(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), k0, i0, k1, i1, w, &_alglib_env_state);
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
Neural network activation function

INPUT PARAMETERS:
    NET         -   neuron input
    K           -   function index (zero for linear function)

OUTPUT PARAMETERS:
    F           -   function
    DF          -   its derivative
    D2F         -   its second derivative

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpactivationfunction(const double net, const ae_int_t k, double &f, double &df, double &d2f)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpactivationfunction(net, k, &f, &df, &d2f, &_alglib_env_state);
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
Procesing

INPUT PARAMETERS:
    Network -   neural network
    X       -   input vector,  array[0..NIn-1].

OUTPUT PARAMETERS:
    Y       -   result. Regression estimate when solving regression  task,
                vector of posterior probabilities for classification task.

See also MLPProcessI

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpprocess(const multilayerperceptron &network, const real_1d_array &x, real_1d_array &y)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpprocess(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::ae_vector*>(y.c_ptr()), &_alglib_env_state);
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
'interactive'  variant  of  MLPProcess  for  languages  like  Python which
support constructs like "Y = MLPProcess(NN,X)" and interactive mode of the
interpreter

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 21.09.2010 by Bochkanov Sergey
*************************************************************************/
void mlpprocessi(const multilayerperceptron &network, const real_1d_array &x, real_1d_array &y)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpprocessi(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::ae_vector*>(y.c_ptr()), &_alglib_env_state);
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
Error function for neural network, internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
double mlperror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::mlperror(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), ssize, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
Natural error function for neural network, internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
double mlperrorn(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::mlperrorn(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), ssize, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
Classification error

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpclserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::ae_int_t result = alglib_impl::mlpclserror(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), ssize, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<ae_int_t*>(&result));
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
Relative classification error on the test set

INPUT PARAMETERS:
    Network -   network
    XY      -   test set
    NPoints -   test set size

RESULT:
    percent of incorrectly classified cases. Works both for
    classifier networks and general purpose networks used as
    classifiers.

  -- ALGLIB --
     Copyright 25.12.2008 by Bochkanov Sergey
*************************************************************************/
double mlprelclserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::mlprelclserror(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
Average cross-entropy (in bits per element) on the test set

INPUT PARAMETERS:
    Network -   neural network
    XY      -   test set
    NPoints -   test set size

RESULT:
    CrossEntropy/(NPoints*LN(2)).
    Zero if network solves regression task.

  -- ALGLIB --
     Copyright 08.01.2009 by Bochkanov Sergey
*************************************************************************/
double mlpavgce(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::mlpavgce(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
RMS error on the test set

INPUT PARAMETERS:
    Network -   neural network
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.
    Its meaning for regression task is obvious. As for
    classification task, RMS error means error when estimating posterior
    probabilities.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
double mlprmserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::mlprmserror(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
Average error on the test set

INPUT PARAMETERS:
    Network -   neural network
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for
    classification task, it means average error when estimating posterior
    probabilities.

  -- ALGLIB --
     Copyright 11.03.2008 by Bochkanov Sergey
*************************************************************************/
double mlpavgerror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::mlpavgerror(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
Average relative error on the test set

INPUT PARAMETERS:
    Network -   neural network
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for
    classification task, it means average relative error when estimating
    posterior probability of belonging to the correct class.

  -- ALGLIB --
     Copyright 11.03.2008 by Bochkanov Sergey
*************************************************************************/
double mlpavgrelerror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::mlpavgrelerror(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
Gradient calculation

INPUT PARAMETERS:
    Network -   network initialized with one of the network creation funcs
    X       -   input vector, length of array must be at least NIn
    DesiredY-   desired outputs, length of array must be at least NOut
    Grad    -   possibly preallocated array. If size of array is smaller
                than WCount, it will be reallocated. It is recommended to
                reuse previously allocated array to reduce allocation
                overhead.

OUTPUT PARAMETERS:
    E       -   error function, SUM(sqr(y[i]-desiredy[i])/2,i)
    Grad    -   gradient of E with respect to weights of network, array[WCount]

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgrad(const multilayerperceptron &network, const real_1d_array &x, const real_1d_array &desiredy, double &e, real_1d_array &grad)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpgrad(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::ae_vector*>(desiredy.c_ptr()), &e, const_cast<alglib_impl::ae_vector*>(grad.c_ptr()), &_alglib_env_state);
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
Gradient calculation (natural error function is used)

INPUT PARAMETERS:
    Network -   network initialized with one of the network creation funcs
    X       -   input vector, length of array must be at least NIn
    DesiredY-   desired outputs, length of array must be at least NOut
    Grad    -   possibly preallocated array. If size of array is smaller
                than WCount, it will be reallocated. It is recommended to
                reuse previously allocated array to reduce allocation
                overhead.

OUTPUT PARAMETERS:
    E       -   error function, sum-of-squares for regression networks,
                cross-entropy for classification networks.
    Grad    -   gradient of E with respect to weights of network, array[WCount]

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgradn(const multilayerperceptron &network, const real_1d_array &x, const real_1d_array &desiredy, double &e, real_1d_array &grad)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpgradn(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::ae_vector*>(desiredy.c_ptr()), &e, const_cast<alglib_impl::ae_vector*>(grad.c_ptr()), &_alglib_env_state);
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
Batch gradient calculation for a set of inputs/outputs

INPUT PARAMETERS:
    Network -   network initialized with one of the network creation funcs
    XY      -   set of inputs/outputs; one sample = one row;
                first NIn columns contain inputs,
                next NOut columns - desired outputs.
    SSize   -   number of elements in XY
    Grad    -   possibly preallocated array. If size of array is smaller
                than WCount, it will be reallocated. It is recommended to
                reuse previously allocated array to reduce allocation
                overhead.

OUTPUT PARAMETERS:
    E       -   error function, SUM(sqr(y[i]-desiredy[i])/2,i)
    Grad    -   gradient of E with respect to weights of network, array[WCount]

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgradbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpgradbatch(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), ssize, &e, const_cast<alglib_impl::ae_vector*>(grad.c_ptr()), &_alglib_env_state);
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
Batch gradient calculation for a set of inputs/outputs
(natural error function is used)

INPUT PARAMETERS:
    Network -   network initialized with one of the network creation funcs
    XY      -   set of inputs/outputs; one sample = one row;
                first NIn columns contain inputs,
                next NOut columns - desired outputs.
    SSize   -   number of elements in XY
    Grad    -   possibly preallocated array. If size of array is smaller
                than WCount, it will be reallocated. It is recommended to
                reuse previously allocated array to reduce allocation
                overhead.

OUTPUT PARAMETERS:
    E       -   error function, sum-of-squares for regression networks,
                cross-entropy for classification networks.
    Grad    -   gradient of E with respect to weights of network, array[WCount]

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgradnbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpgradnbatch(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), ssize, &e, const_cast<alglib_impl::ae_vector*>(grad.c_ptr()), &_alglib_env_state);
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
Batch Hessian calculation (natural error function) using R-algorithm.
Internal subroutine.

  -- ALGLIB --
     Copyright 26.01.2008 by Bochkanov Sergey.

     Hessian calculation based on R-algorithm described in
     "Fast Exact Multiplication by the Hessian",
     B. A. Pearlmutter,
     Neural Computation, 1994.
*************************************************************************/
void mlphessiannbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad, real_2d_array &h)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlphessiannbatch(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), ssize, &e, const_cast<alglib_impl::ae_vector*>(grad.c_ptr()), const_cast<alglib_impl::ae_matrix*>(h.c_ptr()), &_alglib_env_state);
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
Batch Hessian calculation using R-algorithm.
Internal subroutine.

  -- ALGLIB --
     Copyright 26.01.2008 by Bochkanov Sergey.

     Hessian calculation based on R-algorithm described in
     "Fast Exact Multiplication by the Hessian",
     B. A. Pearlmutter,
     Neural Computation, 1994.
*************************************************************************/
void mlphessianbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad, real_2d_array &h)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlphessianbatch(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), ssize, &e, const_cast<alglib_impl::ae_vector*>(grad.c_ptr()), const_cast<alglib_impl::ae_matrix*>(h.c_ptr()), &_alglib_env_state);
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
_logitmodel_owner::_logitmodel_owner()
{
    p_struct = (alglib_impl::logitmodel*)alglib_impl::ae_malloc(sizeof(alglib_impl::logitmodel), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_logitmodel_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_logitmodel_owner::_logitmodel_owner(const _logitmodel_owner &rhs)
{
    p_struct = (alglib_impl::logitmodel*)alglib_impl::ae_malloc(sizeof(alglib_impl::logitmodel), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_logitmodel_init_copy(p_struct, const_cast<alglib_impl::logitmodel*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_logitmodel_owner& _logitmodel_owner::operator=(const _logitmodel_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_logitmodel_clear(p_struct);
    if( !alglib_impl::_logitmodel_init_copy(p_struct, const_cast<alglib_impl::logitmodel*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_logitmodel_owner::~_logitmodel_owner()
{
    alglib_impl::_logitmodel_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::logitmodel* _logitmodel_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::logitmodel* _logitmodel_owner::c_ptr() const
{
    return const_cast<alglib_impl::logitmodel*>(p_struct);
}
logitmodel::logitmodel() : _logitmodel_owner() 
{
}

logitmodel::logitmodel(const logitmodel &rhs):_logitmodel_owner(rhs) 
{
}

logitmodel& logitmodel::operator=(const logitmodel &rhs)
{
    if( this==&rhs )
        return *this;
    _logitmodel_owner::operator=(rhs);
    return *this;
}

logitmodel::~logitmodel()
{
}


/*************************************************************************
MNLReport structure contains information about training process:
* NGrad     -   number of gradient calculations
* NHess     -   number of Hessian calculations
*************************************************************************/
_mnlreport_owner::_mnlreport_owner()
{
    p_struct = (alglib_impl::mnlreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::mnlreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_mnlreport_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_mnlreport_owner::_mnlreport_owner(const _mnlreport_owner &rhs)
{
    p_struct = (alglib_impl::mnlreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::mnlreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_mnlreport_init_copy(p_struct, const_cast<alglib_impl::mnlreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_mnlreport_owner& _mnlreport_owner::operator=(const _mnlreport_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_mnlreport_clear(p_struct);
    if( !alglib_impl::_mnlreport_init_copy(p_struct, const_cast<alglib_impl::mnlreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_mnlreport_owner::~_mnlreport_owner()
{
    alglib_impl::_mnlreport_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::mnlreport* _mnlreport_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::mnlreport* _mnlreport_owner::c_ptr() const
{
    return const_cast<alglib_impl::mnlreport*>(p_struct);
}
mnlreport::mnlreport() : _mnlreport_owner() ,ngrad(p_struct->ngrad),nhess(p_struct->nhess)
{
}

mnlreport::mnlreport(const mnlreport &rhs):_mnlreport_owner(rhs) ,ngrad(p_struct->ngrad),nhess(p_struct->nhess)
{
}

mnlreport& mnlreport::operator=(const mnlreport &rhs)
{
    if( this==&rhs )
        return *this;
    _mnlreport_owner::operator=(rhs);
    return *this;
}

mnlreport::~mnlreport()
{
}

/*************************************************************************
This subroutine trains logit model.

INPUT PARAMETERS:
    XY          -   training set, array[0..NPoints-1,0..NVars]
                    First NVars columns store values of independent
                    variables, next column stores number of class (from 0
                    to NClasses-1) which dataset element belongs to. Fractional
                    values are rounded to nearest integer.
    NPoints     -   training set size, NPoints>=1
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   number of classes, NClasses>=2

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<NVars+2, NVars<1, NClasses<2).
                    *  1, if task has been solved
    LM          -   model built
    Rep         -   training report

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
void mnltrainh(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, logitmodel &lm, mnlreport &rep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mnltrainh(const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, nvars, nclasses, &info, const_cast<alglib_impl::logitmodel*>(lm.c_ptr()), const_cast<alglib_impl::mnlreport*>(rep.c_ptr()), &_alglib_env_state);
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
Procesing

INPUT PARAMETERS:
    LM      -   logit model, passed by non-constant reference
                (some fields of structure are used as temporaries
                when calculating model output).
    X       -   input vector,  array[0..NVars-1].
    Y       -   (possibly) preallocated buffer; if size of Y is less than
                NClasses, it will be reallocated.If it is large enough, it
                is NOT reallocated, so we can save some time on reallocation.

OUTPUT PARAMETERS:
    Y       -   result, array[0..NClasses-1]
                Vector of posterior probabilities for classification task.

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
void mnlprocess(const logitmodel &lm, const real_1d_array &x, real_1d_array &y)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mnlprocess(const_cast<alglib_impl::logitmodel*>(lm.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::ae_vector*>(y.c_ptr()), &_alglib_env_state);
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
'interactive'  variant  of  MNLProcess  for  languages  like  Python which
support constructs like "Y = MNLProcess(LM,X)" and interactive mode of the
interpreter

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
void mnlprocessi(const logitmodel &lm, const real_1d_array &x, real_1d_array &y)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mnlprocessi(const_cast<alglib_impl::logitmodel*>(lm.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::ae_vector*>(y.c_ptr()), &_alglib_env_state);
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
Unpacks coefficients of logit model. Logit model have form:

    P(class=i) = S(i) / (S(0) + S(1) + ... +S(M-1))
          S(i) = Exp(A[i,0]*X[0] + ... + A[i,N-1]*X[N-1] + A[i,N]), when i<M-1
        S(M-1) = 1

INPUT PARAMETERS:
    LM          -   logit model in ALGLIB format

OUTPUT PARAMETERS:
    V           -   coefficients, array[0..NClasses-2,0..NVars]
    NVars       -   number of independent variables
    NClasses    -   number of classes

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
void mnlunpack(const logitmodel &lm, real_2d_array &a, ae_int_t &nvars, ae_int_t &nclasses)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mnlunpack(const_cast<alglib_impl::logitmodel*>(lm.c_ptr()), const_cast<alglib_impl::ae_matrix*>(a.c_ptr()), &nvars, &nclasses, &_alglib_env_state);
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
"Packs" coefficients and creates logit model in ALGLIB format (MNLUnpack
reversed).

INPUT PARAMETERS:
    A           -   model (see MNLUnpack)
    NVars       -   number of independent variables
    NClasses    -   number of classes

OUTPUT PARAMETERS:
    LM          -   logit model.

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
void mnlpack(const real_2d_array &a, const ae_int_t nvars, const ae_int_t nclasses, logitmodel &lm)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mnlpack(const_cast<alglib_impl::ae_matrix*>(a.c_ptr()), nvars, nclasses, const_cast<alglib_impl::logitmodel*>(lm.c_ptr()), &_alglib_env_state);
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
Average cross-entropy (in bits per element) on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    CrossEntropy/(NPoints*ln(2)).

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
double mnlavgce(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::mnlavgce(const_cast<alglib_impl::logitmodel*>(lm.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
Relative classification error on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    percent of incorrectly classified cases.

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
double mnlrelclserror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::mnlrelclserror(const_cast<alglib_impl::logitmodel*>(lm.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
RMS error on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error (error when estimating posterior probabilities).

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double mnlrmserror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::mnlrmserror(const_cast<alglib_impl::logitmodel*>(lm.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
Average error on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average error (error when estimating posterior probabilities).

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double mnlavgerror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::mnlavgerror(const_cast<alglib_impl::logitmodel*>(lm.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
Average relative error on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average relative error (error when estimating posterior probabilities).

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double mnlavgrelerror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t ssize)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::mnlavgrelerror(const_cast<alglib_impl::logitmodel*>(lm.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), ssize, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
Classification error on test set = MNLRelClsError*NPoints

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
ae_int_t mnlclserror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::ae_int_t result = alglib_impl::mnlclserror(const_cast<alglib_impl::logitmodel*>(lm.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<ae_int_t*>(&result));
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
This structure is a MCPD (Markov Chains for Population Data) solver.

You should use ALGLIB functions in order to work with this object.

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
_mcpdstate_owner::_mcpdstate_owner()
{
    p_struct = (alglib_impl::mcpdstate*)alglib_impl::ae_malloc(sizeof(alglib_impl::mcpdstate), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_mcpdstate_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_mcpdstate_owner::_mcpdstate_owner(const _mcpdstate_owner &rhs)
{
    p_struct = (alglib_impl::mcpdstate*)alglib_impl::ae_malloc(sizeof(alglib_impl::mcpdstate), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_mcpdstate_init_copy(p_struct, const_cast<alglib_impl::mcpdstate*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_mcpdstate_owner& _mcpdstate_owner::operator=(const _mcpdstate_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_mcpdstate_clear(p_struct);
    if( !alglib_impl::_mcpdstate_init_copy(p_struct, const_cast<alglib_impl::mcpdstate*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_mcpdstate_owner::~_mcpdstate_owner()
{
    alglib_impl::_mcpdstate_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::mcpdstate* _mcpdstate_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::mcpdstate* _mcpdstate_owner::c_ptr() const
{
    return const_cast<alglib_impl::mcpdstate*>(p_struct);
}
mcpdstate::mcpdstate() : _mcpdstate_owner() 
{
}

mcpdstate::mcpdstate(const mcpdstate &rhs):_mcpdstate_owner(rhs) 
{
}

mcpdstate& mcpdstate::operator=(const mcpdstate &rhs)
{
    if( this==&rhs )
        return *this;
    _mcpdstate_owner::operator=(rhs);
    return *this;
}

mcpdstate::~mcpdstate()
{
}


/*************************************************************************
This structure is a MCPD training report:
    InnerIterationsCount    -   number of inner iterations of the
                                underlying optimization algorithm
    OuterIterationsCount    -   number of outer iterations of the
                                underlying optimization algorithm
    NFEV                    -   number of merit function evaluations
    TerminationType         -   termination type
                                (same as for MinBLEIC optimizer, positive
                                values denote success, negative ones -
                                failure)

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
_mcpdreport_owner::_mcpdreport_owner()
{
    p_struct = (alglib_impl::mcpdreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::mcpdreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_mcpdreport_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_mcpdreport_owner::_mcpdreport_owner(const _mcpdreport_owner &rhs)
{
    p_struct = (alglib_impl::mcpdreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::mcpdreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_mcpdreport_init_copy(p_struct, const_cast<alglib_impl::mcpdreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_mcpdreport_owner& _mcpdreport_owner::operator=(const _mcpdreport_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_mcpdreport_clear(p_struct);
    if( !alglib_impl::_mcpdreport_init_copy(p_struct, const_cast<alglib_impl::mcpdreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_mcpdreport_owner::~_mcpdreport_owner()
{
    alglib_impl::_mcpdreport_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::mcpdreport* _mcpdreport_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::mcpdreport* _mcpdreport_owner::c_ptr() const
{
    return const_cast<alglib_impl::mcpdreport*>(p_struct);
}
mcpdreport::mcpdreport() : _mcpdreport_owner() ,inneriterationscount(p_struct->inneriterationscount),outeriterationscount(p_struct->outeriterationscount),nfev(p_struct->nfev),terminationtype(p_struct->terminationtype)
{
}

mcpdreport::mcpdreport(const mcpdreport &rhs):_mcpdreport_owner(rhs) ,inneriterationscount(p_struct->inneriterationscount),outeriterationscount(p_struct->outeriterationscount),nfev(p_struct->nfev),terminationtype(p_struct->terminationtype)
{
}

mcpdreport& mcpdreport::operator=(const mcpdreport &rhs)
{
    if( this==&rhs )
        return *this;
    _mcpdreport_owner::operator=(rhs);
    return *this;
}

mcpdreport::~mcpdreport()
{
}

/*************************************************************************
DESCRIPTION:

This function creates MCPD (Markov Chains for Population Data) solver.

This  solver  can  be  used  to find transition matrix P for N-dimensional
prediction  problem  where transition from X[i] to X[i+1] is  modelled  as
    X[i+1] = P*X[i]
where X[i] and X[i+1] are N-dimensional population vectors (components  of
each X are non-negative), and P is a N*N transition matrix (elements of  P
are non-negative, each column sums to 1.0).

Such models arise when when:
* there is some population of individuals
* individuals can have different states
* individuals can transit from one state to another
* population size is constant, i.e. there is no new individuals and no one
  leaves population
* you want to model transitions of individuals from one state into another

USAGE:

Here we give very brief outline of the MCPD. We strongly recommend you  to
read examples in the ALGLIB Reference Manual and to read ALGLIB User Guide
on data analysis which is available at http://www.alglib.net/dataanalysis/

1. User initializes algorithm state with MCPDCreate() call

2. User  adds  one  or  more  tracks -  sequences of states which describe
   evolution of a system being modelled from different starting conditions

3. User may add optional boundary, equality  and/or  linear constraints on
   the coefficients of P by calling one of the following functions:
   * MCPDSetEC() to set equality constraints
   * MCPDSetBC() to set bound constraints
   * MCPDSetLC() to set linear constraints

4. Optionally,  user  may  set  custom  weights  for prediction errors (by
   default, algorithm assigns non-equal, automatically chosen weights  for
   errors in the prediction of different components of X). It can be  done
   with a call of MCPDSetPredictionWeights() function.

5. User calls MCPDSolve() function which takes algorithm  state and
   pointer (delegate, etc.) to callback function which calculates F/G.

6. User calls MCPDResults() to get solution

INPUT PARAMETERS:
    N       -   problem dimension, N>=1

OUTPUT PARAMETERS:
    State   -   structure stores algorithm state

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdcreate(const ae_int_t n, mcpdstate &s)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mcpdcreate(n, const_cast<alglib_impl::mcpdstate*>(s.c_ptr()), &_alglib_env_state);
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
DESCRIPTION:

This function is a specialized version of MCPDCreate()  function,  and  we
recommend  you  to read comments for this function for general information
about MCPD solver.

This  function  creates  MCPD (Markov Chains for Population  Data)  solver
for "Entry-state" model,  i.e. model  where transition from X[i] to X[i+1]
is modelled as
    X[i+1] = P*X[i]
where
    X[i] and X[i+1] are N-dimensional state vectors
    P is a N*N transition matrix
and  one  selected component of X[] is called "entry" state and is treated
in a special way:
    system state always transits from "entry" state to some another state
    system state can not transit from any state into "entry" state
Such conditions basically mean that row of P which corresponds to  "entry"
state is zero.

Such models arise when:
* there is some population of individuals
* individuals can have different states
* individuals can transit from one state to another
* population size is NOT constant -  at every moment of time there is some
  (unpredictable) amount of "new" individuals, which can transit into  one
  of the states at the next turn, but still no one leaves population
* you want to model transitions of individuals from one state into another
* but you do NOT want to predict amount of "new"  individuals  because  it
  does not depends on individuals already present (hence  system  can  not
  transit INTO entry state - it can only transit FROM it).

This model is discussed  in  more  details  in  the ALGLIB User Guide (see
http://www.alglib.net/dataanalysis/ for more data).

INPUT PARAMETERS:
    N       -   problem dimension, N>=2
    EntryState- index of entry state, in 0..N-1

OUTPUT PARAMETERS:
    State   -   structure stores algorithm state

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdcreateentry(const ae_int_t n, const ae_int_t entrystate, mcpdstate &s)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mcpdcreateentry(n, entrystate, const_cast<alglib_impl::mcpdstate*>(s.c_ptr()), &_alglib_env_state);
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
DESCRIPTION:

This function is a specialized version of MCPDCreate()  function,  and  we
recommend  you  to read comments for this function for general information
about MCPD solver.

This  function  creates  MCPD (Markov Chains for Population  Data)  solver
for "Exit-state" model,  i.e. model  where  transition from X[i] to X[i+1]
is modelled as
    X[i+1] = P*X[i]
where
    X[i] and X[i+1] are N-dimensional state vectors
    P is a N*N transition matrix
and  one  selected component of X[] is called "exit"  state and is treated
in a special way:
    system state can transit from any state into "exit" state
    system state can not transit from "exit" state into any other state
    transition operator discards "exit" state (makes it zero at each turn)
Such  conditions  basically  mean  that  column  of P which corresponds to
"exit" state is zero. Multiplication by such P may decrease sum of  vector
components.

Such models arise when:
* there is some population of individuals
* individuals can have different states
* individuals can transit from one state to another
* population size is NOT constant - individuals can move into "exit" state
  and leave population at the next turn, but there are no new individuals
* amount of individuals which leave population can be predicted
* you want to model transitions of individuals from one state into another
  (including transitions into the "exit" state)

This model is discussed  in  more  details  in  the ALGLIB User Guide (see
http://www.alglib.net/dataanalysis/ for more data).

INPUT PARAMETERS:
    N       -   problem dimension, N>=2
    ExitState-  index of exit state, in 0..N-1

OUTPUT PARAMETERS:
    State   -   structure stores algorithm state

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdcreateexit(const ae_int_t n, const ae_int_t exitstate, mcpdstate &s)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mcpdcreateexit(n, exitstate, const_cast<alglib_impl::mcpdstate*>(s.c_ptr()), &_alglib_env_state);
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
DESCRIPTION:

This function is a specialized version of MCPDCreate()  function,  and  we
recommend  you  to read comments for this function for general information
about MCPD solver.

This  function  creates  MCPD (Markov Chains for Population  Data)  solver
for "Entry-Exit-states" model, i.e. model where  transition  from  X[i] to
X[i+1] is modelled as
    X[i+1] = P*X[i]
where
    X[i] and X[i+1] are N-dimensional state vectors
    P is a N*N transition matrix
one selected component of X[] is called "entry" state and is treated in  a
special way:
    system state always transits from "entry" state to some another state
    system state can not transit from any state into "entry" state
and another one component of X[] is called "exit" state and is treated  in
a special way too:
    system state can transit from any state into "exit" state
    system state can not transit from "exit" state into any other state
    transition operator discards "exit" state (makes it zero at each turn)
Such conditions basically mean that:
    row of P which corresponds to "entry" state is zero
    column of P which corresponds to "exit" state is zero
Multiplication by such P may decrease sum of vector components.

Such models arise when:
* there is some population of individuals
* individuals can have different states
* individuals can transit from one state to another
* population size is NOT constant
* at every moment of time there is some (unpredictable)  amount  of  "new"
  individuals, which can transit into one of the states at the next turn
* some  individuals  can  move  (predictably)  into "exit" state and leave
  population at the next turn
* you want to model transitions of individuals from one state into another,
  including transitions from the "entry" state and into the "exit" state.
* but you do NOT want to predict amount of "new"  individuals  because  it
  does not depends on individuals already present (hence  system  can  not
  transit INTO entry state - it can only transit FROM it).

This model is discussed  in  more  details  in  the ALGLIB User Guide (see
http://www.alglib.net/dataanalysis/ for more data).

INPUT PARAMETERS:
    N       -   problem dimension, N>=2
    EntryState- index of entry state, in 0..N-1
    ExitState-  index of exit state, in 0..N-1

OUTPUT PARAMETERS:
    State   -   structure stores algorithm state

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdcreateentryexit(const ae_int_t n, const ae_int_t entrystate, const ae_int_t exitstate, mcpdstate &s)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mcpdcreateentryexit(n, entrystate, exitstate, const_cast<alglib_impl::mcpdstate*>(s.c_ptr()), &_alglib_env_state);
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
This  function  is  used to add a track - sequence of system states at the
different moments of its evolution.

You  may  add  one  or several tracks to the MCPD solver. In case you have
several tracks, they won't overwrite each other. For example,  if you pass
two tracks, A1-A2-A3 (system at t=A+1, t=A+2 and t=A+3) and B1-B2-B3, then
solver will try to model transitions from t=A+1 to t=A+2, t=A+2 to  t=A+3,
t=B+1 to t=B+2, t=B+2 to t=B+3. But it WONT mix these two tracks - i.e. it
wont try to model transition from t=A+3 to t=B+1.

INPUT PARAMETERS:
    S       -   solver
    XY      -   track, array[K,N]:
                * I-th row is a state at t=I
                * elements of XY must be non-negative (exception will be
                  thrown on negative elements)
    K       -   number of points in a track
                * if given, only leading K rows of XY are used
                * if not given, automatically determined from size of XY

NOTES:

1. Track may contain either proportional or population data:
   * with proportional data all rows of XY must sum to 1.0, i.e. we have
     proportions instead of absolute population values
   * with population data rows of XY contain population counts and generally
     do not sum to 1.0 (although they still must be non-negative)

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdaddtrack(const mcpdstate &s, const real_2d_array &xy, const ae_int_t k)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mcpdaddtrack(const_cast<alglib_impl::mcpdstate*>(s.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), k, &_alglib_env_state);
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
This  function  is  used to add a track - sequence of system states at the
different moments of its evolution.

You  may  add  one  or several tracks to the MCPD solver. In case you have
several tracks, they won't overwrite each other. For example,  if you pass
two tracks, A1-A2-A3 (system at t=A+1, t=A+2 and t=A+3) and B1-B2-B3, then
solver will try to model transitions from t=A+1 to t=A+2, t=A+2 to  t=A+3,
t=B+1 to t=B+2, t=B+2 to t=B+3. But it WONT mix these two tracks - i.e. it
wont try to model transition from t=A+3 to t=B+1.

INPUT PARAMETERS:
    S       -   solver
    XY      -   track, array[K,N]:
                * I-th row is a state at t=I
                * elements of XY must be non-negative (exception will be
                  thrown on negative elements)
    K       -   number of points in a track
                * if given, only leading K rows of XY are used
                * if not given, automatically determined from size of XY

NOTES:

1. Track may contain either proportional or population data:
   * with proportional data all rows of XY must sum to 1.0, i.e. we have
     proportions instead of absolute population values
   * with population data rows of XY contain population counts and generally
     do not sum to 1.0 (although they still must be non-negative)

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdaddtrack(const mcpdstate &s, const real_2d_array &xy)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t k;

    k = xy.rows();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mcpdaddtrack(const_cast<alglib_impl::mcpdstate*>(s.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), k, &_alglib_env_state);

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
This function is used to add equality constraints on the elements  of  the
transition matrix P.

MCPD solver has four types of constraints which can be placed on P:
* user-specified equality constraints (optional)
* user-specified bound constraints (optional)
* user-specified general linear constraints (optional)
* basic constraints (always present):
  * non-negativity: P[i,j]>=0
  * consistency: every column of P sums to 1.0

Final  constraints  which  are  passed  to  the  underlying  optimizer are
calculated  as  intersection  of all present constraints. For example, you
may specify boundary constraint on P[0,0] and equality one:
    0.1<=P[0,0]<=0.9
    P[0,0]=0.5
Such  combination  of  constraints  will  be  silently  reduced  to  their
intersection, which is P[0,0]=0.5.

This  function  can  be  used  to  place equality constraints on arbitrary
subset of elements of P. Set of constraints is specified by EC, which  may
contain either NAN's or finite numbers from [0,1]. NAN denotes absence  of
constraint, finite number denotes equality constraint on specific  element
of P.

You can also  use  MCPDAddEC()  function  which  allows  to  ADD  equality
constraint  for  one  element  of P without changing constraints for other
elements.

These functions (MCPDSetEC and MCPDAddEC) interact as follows:
* there is internal matrix of equality constraints which is stored in  the
  MCPD solver
* MCPDSetEC() replaces this matrix by another one (SET)
* MCPDAddEC() modifies one element of this matrix and  leaves  other  ones
  unchanged (ADD)
* thus  MCPDAddEC()  call  preserves  all  modifications  done by previous
  calls,  while  MCPDSetEC()  completely discards all changes  done to the
  equality constraints.

INPUT PARAMETERS:
    S       -   solver
    EC      -   equality constraints, array[N,N]. Elements of  EC  can  be
                either NAN's or finite  numbers from  [0,1].  NAN  denotes
                absence  of  constraints,  while  finite  value    denotes
                equality constraint on the corresponding element of P.

NOTES:

1. infinite values of EC will lead to exception being thrown. Values  less
than 0.0 or greater than 1.0 will lead to error code being returned  after
call to MCPDSolve().

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsetec(const mcpdstate &s, const real_2d_array &ec)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mcpdsetec(const_cast<alglib_impl::mcpdstate*>(s.c_ptr()), const_cast<alglib_impl::ae_matrix*>(ec.c_ptr()), &_alglib_env_state);
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
This function is used to add equality constraints on the elements  of  the
transition matrix P.

MCPD solver has four types of constraints which can be placed on P:
* user-specified equality constraints (optional)
* user-specified bound constraints (optional)
* user-specified general linear constraints (optional)
* basic constraints (always present):
  * non-negativity: P[i,j]>=0
  * consistency: every column of P sums to 1.0

Final  constraints  which  are  passed  to  the  underlying  optimizer are
calculated  as  intersection  of all present constraints. For example, you
may specify boundary constraint on P[0,0] and equality one:
    0.1<=P[0,0]<=0.9
    P[0,0]=0.5
Such  combination  of  constraints  will  be  silently  reduced  to  their
intersection, which is P[0,0]=0.5.

This function can be used to ADD equality constraint for one element of  P
without changing constraints for other elements.

You  can  also  use  MCPDSetEC()  function  which  allows  you  to specify
arbitrary set of equality constraints in one call.

These functions (MCPDSetEC and MCPDAddEC) interact as follows:
* there is internal matrix of equality constraints which is stored in the
  MCPD solver
* MCPDSetEC() replaces this matrix by another one (SET)
* MCPDAddEC() modifies one element of this matrix and leaves  other  ones
  unchanged (ADD)
* thus  MCPDAddEC()  call  preserves  all  modifications done by previous
  calls,  while  MCPDSetEC()  completely discards all changes done to the
  equality constraints.

INPUT PARAMETERS:
    S       -   solver
    I       -   row index of element being constrained
    J       -   column index of element being constrained
    C       -   value (constraint for P[I,J]).  Can  be  either  NAN  (no
                constraint) or finite value from [0,1].

NOTES:

1. infinite values of C  will lead to exception being thrown. Values  less
than 0.0 or greater than 1.0 will lead to error code being returned  after
call to MCPDSolve().

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdaddec(const mcpdstate &s, const ae_int_t i, const ae_int_t j, const double c)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mcpdaddec(const_cast<alglib_impl::mcpdstate*>(s.c_ptr()), i, j, c, &_alglib_env_state);
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
This function is used to add bound constraints  on  the  elements  of  the
transition matrix P.

MCPD solver has four types of constraints which can be placed on P:
* user-specified equality constraints (optional)
* user-specified bound constraints (optional)
* user-specified general linear constraints (optional)
* basic constraints (always present):
  * non-negativity: P[i,j]>=0
  * consistency: every column of P sums to 1.0

Final  constraints  which  are  passed  to  the  underlying  optimizer are
calculated  as  intersection  of all present constraints. For example, you
may specify boundary constraint on P[0,0] and equality one:
    0.1<=P[0,0]<=0.9
    P[0,0]=0.5
Such  combination  of  constraints  will  be  silently  reduced  to  their
intersection, which is P[0,0]=0.5.

This  function  can  be  used  to  place bound   constraints  on arbitrary
subset  of  elements  of  P.  Set of constraints is specified by BndL/BndU
matrices, which may contain arbitrary combination  of  finite  numbers  or
infinities (like -INF<x<=0.5 or 0.1<=x<+INF).

You can also use MCPDAddBC() function which allows to ADD bound constraint
for one element of P without changing constraints for other elements.

These functions (MCPDSetBC and MCPDAddBC) interact as follows:
* there is internal matrix of bound constraints which is stored in the
  MCPD solver
* MCPDSetBC() replaces this matrix by another one (SET)
* MCPDAddBC() modifies one element of this matrix and  leaves  other  ones
  unchanged (ADD)
* thus  MCPDAddBC()  call  preserves  all  modifications  done by previous
  calls,  while  MCPDSetBC()  completely discards all changes  done to the
  equality constraints.

INPUT PARAMETERS:
    S       -   solver
    BndL    -   lower bounds constraints, array[N,N]. Elements of BndL can
                be finite numbers or -INF.
    BndU    -   upper bounds constraints, array[N,N]. Elements of BndU can
                be finite numbers or +INF.

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsetbc(const mcpdstate &s, const real_2d_array &bndl, const real_2d_array &bndu)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mcpdsetbc(const_cast<alglib_impl::mcpdstate*>(s.c_ptr()), const_cast<alglib_impl::ae_matrix*>(bndl.c_ptr()), const_cast<alglib_impl::ae_matrix*>(bndu.c_ptr()), &_alglib_env_state);
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
This function is used to add bound constraints  on  the  elements  of  the
transition matrix P.

MCPD solver has four types of constraints which can be placed on P:
* user-specified equality constraints (optional)
* user-specified bound constraints (optional)
* user-specified general linear constraints (optional)
* basic constraints (always present):
  * non-negativity: P[i,j]>=0
  * consistency: every column of P sums to 1.0

Final  constraints  which  are  passed  to  the  underlying  optimizer are
calculated  as  intersection  of all present constraints. For example, you
may specify boundary constraint on P[0,0] and equality one:
    0.1<=P[0,0]<=0.9
    P[0,0]=0.5
Such  combination  of  constraints  will  be  silently  reduced  to  their
intersection, which is P[0,0]=0.5.

This  function  can  be  used to ADD bound constraint for one element of P
without changing constraints for other elements.

You  can  also  use  MCPDSetBC()  function  which  allows to  place  bound
constraints  on arbitrary subset of elements of P.   Set of constraints is
specified  by  BndL/BndU matrices, which may contain arbitrary combination
of finite numbers or infinities (like -INF<x<=0.5 or 0.1<=x<+INF).

These functions (MCPDSetBC and MCPDAddBC) interact as follows:
* there is internal matrix of bound constraints which is stored in the
  MCPD solver
* MCPDSetBC() replaces this matrix by another one (SET)
* MCPDAddBC() modifies one element of this matrix and  leaves  other  ones
  unchanged (ADD)
* thus  MCPDAddBC()  call  preserves  all  modifications  done by previous
  calls,  while  MCPDSetBC()  completely discards all changes  done to the
  equality constraints.

INPUT PARAMETERS:
    S       -   solver
    I       -   row index of element being constrained
    J       -   column index of element being constrained
    BndL    -   lower bound
    BndU    -   upper bound

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdaddbc(const mcpdstate &s, const ae_int_t i, const ae_int_t j, const double bndl, const double bndu)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mcpdaddbc(const_cast<alglib_impl::mcpdstate*>(s.c_ptr()), i, j, bndl, bndu, &_alglib_env_state);
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
This function is used to set linear equality/inequality constraints on the
elements of the transition matrix P.

This function can be used to set one or several general linear constraints
on the elements of P. Two types of constraints are supported:
* equality constraints
* inequality constraints (both less-or-equal and greater-or-equal)

Coefficients  of  constraints  are  specified  by  matrix  C (one  of  the
parameters).  One  row  of  C  corresponds  to  one  constraint.   Because
transition  matrix P has N*N elements,  we  need  N*N columns to store all
coefficients  (they  are  stored row by row), and one more column to store
right part - hence C has N*N+1 columns.  Constraint  kind is stored in the
CT array.

Thus, I-th linear constraint is
    P[0,0]*C[I,0] + P[0,1]*C[I,1] + .. + P[0,N-1]*C[I,N-1] +
        + P[1,0]*C[I,N] + P[1,1]*C[I,N+1] + ... +
        + P[N-1,N-1]*C[I,N*N-1]  ?=?  C[I,N*N]
where ?=? can be either "=" (CT[i]=0), "<=" (CT[i]<0) or ">=" (CT[i]>0).

Your constraint may involve only some subset of P (less than N*N elements).
For example it can be something like
    P[0,0] + P[0,1] = 0.5
In this case you still should pass matrix  with N*N+1 columns, but all its
elements (except for C[0,0], C[0,1] and C[0,N*N-1]) will be zero.

INPUT PARAMETERS:
    S       -   solver
    C       -   array[K,N*N+1] - coefficients of constraints
                (see above for complete description)
    CT      -   array[K] - constraint types
                (see above for complete description)
    K       -   number of equality/inequality constraints, K>=0:
                * if given, only leading K elements of C/CT are used
                * if not given, automatically determined from sizes of C/CT

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsetlc(const mcpdstate &s, const real_2d_array &c, const integer_1d_array &ct, const ae_int_t k)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mcpdsetlc(const_cast<alglib_impl::mcpdstate*>(s.c_ptr()), const_cast<alglib_impl::ae_matrix*>(c.c_ptr()), const_cast<alglib_impl::ae_vector*>(ct.c_ptr()), k, &_alglib_env_state);
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
This function is used to set linear equality/inequality constraints on the
elements of the transition matrix P.

This function can be used to set one or several general linear constraints
on the elements of P. Two types of constraints are supported:
* equality constraints
* inequality constraints (both less-or-equal and greater-or-equal)

Coefficients  of  constraints  are  specified  by  matrix  C (one  of  the
parameters).  One  row  of  C  corresponds  to  one  constraint.   Because
transition  matrix P has N*N elements,  we  need  N*N columns to store all
coefficients  (they  are  stored row by row), and one more column to store
right part - hence C has N*N+1 columns.  Constraint  kind is stored in the
CT array.

Thus, I-th linear constraint is
    P[0,0]*C[I,0] + P[0,1]*C[I,1] + .. + P[0,N-1]*C[I,N-1] +
        + P[1,0]*C[I,N] + P[1,1]*C[I,N+1] + ... +
        + P[N-1,N-1]*C[I,N*N-1]  ?=?  C[I,N*N]
where ?=? can be either "=" (CT[i]=0), "<=" (CT[i]<0) or ">=" (CT[i]>0).

Your constraint may involve only some subset of P (less than N*N elements).
For example it can be something like
    P[0,0] + P[0,1] = 0.5
In this case you still should pass matrix  with N*N+1 columns, but all its
elements (except for C[0,0], C[0,1] and C[0,N*N-1]) will be zero.

INPUT PARAMETERS:
    S       -   solver
    C       -   array[K,N*N+1] - coefficients of constraints
                (see above for complete description)
    CT      -   array[K] - constraint types
                (see above for complete description)
    K       -   number of equality/inequality constraints, K>=0:
                * if given, only leading K elements of C/CT are used
                * if not given, automatically determined from sizes of C/CT

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsetlc(const mcpdstate &s, const real_2d_array &c, const integer_1d_array &ct)
{
    alglib_impl::ae_state _alglib_env_state;    
    ae_int_t k;
    if( (c.rows()!=ct.length()))
        throw ap_error("Error while calling 'mcpdsetlc': looks like one of arguments has wrong size");
    k = c.rows();
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mcpdsetlc(const_cast<alglib_impl::mcpdstate*>(s.c_ptr()), const_cast<alglib_impl::ae_matrix*>(c.c_ptr()), const_cast<alglib_impl::ae_vector*>(ct.c_ptr()), k, &_alglib_env_state);

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
This function allows to  tune  amount  of  Tikhonov  regularization  being
applied to your problem.

By default, regularizing term is equal to r*||P-prior_P||^2, where r is  a
small non-zero value,  P is transition matrix, prior_P is identity matrix,
||X||^2 is a sum of squared elements of X.

This  function  allows  you to change coefficient r. You can  also  change
prior values with MCPDSetPrior() function.

INPUT PARAMETERS:
    S       -   solver
    V       -   regularization  coefficient, finite non-negative value. It
                is  not  recommended  to specify zero value unless you are
                pretty sure that you want it.

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsettikhonovregularizer(const mcpdstate &s, const double v)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mcpdsettikhonovregularizer(const_cast<alglib_impl::mcpdstate*>(s.c_ptr()), v, &_alglib_env_state);
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
This  function  allows to set prior values used for regularization of your
problem.

By default, regularizing term is equal to r*||P-prior_P||^2, where r is  a
small non-zero value,  P is transition matrix, prior_P is identity matrix,
||X||^2 is a sum of squared elements of X.

This  function  allows  you to change prior values prior_P. You  can  also
change r with MCPDSetTikhonovRegularizer() function.

INPUT PARAMETERS:
    S       -   solver
    PP      -   array[N,N], matrix of prior values:
                1. elements must be real numbers from [0,1]
                2. columns must sum to 1.0.
                First property is checked (exception is thrown otherwise),
                while second one is not checked/enforced.

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsetprior(const mcpdstate &s, const real_2d_array &pp)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mcpdsetprior(const_cast<alglib_impl::mcpdstate*>(s.c_ptr()), const_cast<alglib_impl::ae_matrix*>(pp.c_ptr()), &_alglib_env_state);
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
This function is used to change prediction weights

MCPD solver scales prediction errors as follows
    Error(P) = ||W*(y-P*x)||^2
where
    x is a system state at time t
    y is a system state at time t+1
    P is a transition matrix
    W is a diagonal scaling matrix

By default, weights are chosen in order  to  minimize  relative prediction
error instead of absolute one. For example, if one component of  state  is
about 0.5 in magnitude and another one is about 0.05, then algorithm  will
make corresponding weights equal to 2.0 and 20.0.

INPUT PARAMETERS:
    S       -   solver
    PW      -   array[N], weights:
                * must be non-negative values (exception will be thrown otherwise)
                * zero values will be replaced by automatically chosen values

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsetpredictionweights(const mcpdstate &s, const real_1d_array &pw)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mcpdsetpredictionweights(const_cast<alglib_impl::mcpdstate*>(s.c_ptr()), const_cast<alglib_impl::ae_vector*>(pw.c_ptr()), &_alglib_env_state);
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
This function is used to start solution of the MCPD problem.

After return from this function, you can use MCPDResults() to get solution
and completion code.

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsolve(const mcpdstate &s)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mcpdsolve(const_cast<alglib_impl::mcpdstate*>(s.c_ptr()), &_alglib_env_state);
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
MCPD results

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    P       -   array[N,N], transition matrix
    Rep     -   optimization report. You should check Rep.TerminationType
                in  order  to  distinguish  successful  termination  from
                unsuccessful one. Speaking short, positive values  denote
                success, negative ones are failures.
                More information about fields of this  structure  can  be
                found in the comments on MCPDReport datatype.


  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdresults(const mcpdstate &s, real_2d_array &p, mcpdreport &rep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mcpdresults(const_cast<alglib_impl::mcpdstate*>(s.c_ptr()), const_cast<alglib_impl::ae_matrix*>(p.c_ptr()), const_cast<alglib_impl::mcpdreport*>(rep.c_ptr()), &_alglib_env_state);
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
Training report:
    * NGrad     - number of gradient calculations
    * NHess     - number of Hessian calculations
    * NCholesky - number of Cholesky decompositions
*************************************************************************/
_mlpreport_owner::_mlpreport_owner()
{
    p_struct = (alglib_impl::mlpreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::mlpreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_mlpreport_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_mlpreport_owner::_mlpreport_owner(const _mlpreport_owner &rhs)
{
    p_struct = (alglib_impl::mlpreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::mlpreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_mlpreport_init_copy(p_struct, const_cast<alglib_impl::mlpreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_mlpreport_owner& _mlpreport_owner::operator=(const _mlpreport_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_mlpreport_clear(p_struct);
    if( !alglib_impl::_mlpreport_init_copy(p_struct, const_cast<alglib_impl::mlpreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_mlpreport_owner::~_mlpreport_owner()
{
    alglib_impl::_mlpreport_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::mlpreport* _mlpreport_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::mlpreport* _mlpreport_owner::c_ptr() const
{
    return const_cast<alglib_impl::mlpreport*>(p_struct);
}
mlpreport::mlpreport() : _mlpreport_owner() ,ngrad(p_struct->ngrad),nhess(p_struct->nhess),ncholesky(p_struct->ncholesky)
{
}

mlpreport::mlpreport(const mlpreport &rhs):_mlpreport_owner(rhs) ,ngrad(p_struct->ngrad),nhess(p_struct->nhess),ncholesky(p_struct->ncholesky)
{
}

mlpreport& mlpreport::operator=(const mlpreport &rhs)
{
    if( this==&rhs )
        return *this;
    _mlpreport_owner::operator=(rhs);
    return *this;
}

mlpreport::~mlpreport()
{
}


/*************************************************************************
Cross-validation estimates of generalization error
*************************************************************************/
_mlpcvreport_owner::_mlpcvreport_owner()
{
    p_struct = (alglib_impl::mlpcvreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::mlpcvreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_mlpcvreport_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_mlpcvreport_owner::_mlpcvreport_owner(const _mlpcvreport_owner &rhs)
{
    p_struct = (alglib_impl::mlpcvreport*)alglib_impl::ae_malloc(sizeof(alglib_impl::mlpcvreport), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_mlpcvreport_init_copy(p_struct, const_cast<alglib_impl::mlpcvreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_mlpcvreport_owner& _mlpcvreport_owner::operator=(const _mlpcvreport_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_mlpcvreport_clear(p_struct);
    if( !alglib_impl::_mlpcvreport_init_copy(p_struct, const_cast<alglib_impl::mlpcvreport*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_mlpcvreport_owner::~_mlpcvreport_owner()
{
    alglib_impl::_mlpcvreport_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::mlpcvreport* _mlpcvreport_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::mlpcvreport* _mlpcvreport_owner::c_ptr() const
{
    return const_cast<alglib_impl::mlpcvreport*>(p_struct);
}
mlpcvreport::mlpcvreport() : _mlpcvreport_owner() ,relclserror(p_struct->relclserror),avgce(p_struct->avgce),rmserror(p_struct->rmserror),avgerror(p_struct->avgerror),avgrelerror(p_struct->avgrelerror)
{
}

mlpcvreport::mlpcvreport(const mlpcvreport &rhs):_mlpcvreport_owner(rhs) ,relclserror(p_struct->relclserror),avgce(p_struct->avgce),rmserror(p_struct->rmserror),avgerror(p_struct->avgerror),avgrelerror(p_struct->avgrelerror)
{
}

mlpcvreport& mlpcvreport::operator=(const mlpcvreport &rhs)
{
    if( this==&rhs )
        return *this;
    _mlpcvreport_owner::operator=(rhs);
    return *this;
}

mlpcvreport::~mlpcvreport()
{
}

/*************************************************************************
Neural network training  using  modified  Levenberg-Marquardt  with  exact
Hessian calculation and regularization. Subroutine trains  neural  network
with restarts from random positions. Algorithm is well  suited  for  small
and medium scale problems (hundreds of weights).

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay constant, >=0.001
                    Decay term 'Decay*||Weights||^2' is added to error
                    function.
                    If you don't know what Decay to choose, use 0.001.
    Restarts    -   number of restarts from random position, >0.
                    If you don't know what Restarts to choose, use 2.

OUTPUT PARAMETERS:
    Network     -   trained neural network.
    Info        -   return code:
                    * -9, if internal matrix inverse subroutine failed
                    * -2, if there is a point with class number
                          outside of [0..NOut-1].
                    * -1, if wrong parameters specified
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void mlptrainlm(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlptrainlm(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, decay, restarts, &info, const_cast<alglib_impl::mlpreport*>(rep.c_ptr()), &_alglib_env_state);
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
Neural  network  training  using  L-BFGS  algorithm  with  regularization.
Subroutine  trains  neural  network  with  restarts from random positions.
Algorithm  is  well  suited  for  problems  of  any dimensionality (memory
requirements and step complexity are linear by weights number).

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay constant, >=0.001
                    Decay term 'Decay*||Weights||^2' is added to error
                    function.
                    If you don't know what Decay to choose, use 0.001.
    Restarts    -   number of restarts from random position, >0.
                    If you don't know what Restarts to choose, use 2.
    WStep       -   stopping criterion. Algorithm stops if  step  size  is
                    less than WStep. Recommended value - 0.01.  Zero  step
                    size means stopping after MaxIts iterations.
    MaxIts      -   stopping   criterion.  Algorithm  stops  after  MaxIts
                    iterations (NOT gradient  calculations).  Zero  MaxIts
                    means stopping when step is sufficiently small.

OUTPUT PARAMETERS:
    Network     -   trained neural network.
    Info        -   return code:
                    * -8, if both WStep=0 and MaxIts=0
                    * -2, if there is a point with class number
                          outside of [0..NOut-1].
                    * -1, if wrong parameters specified
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report

  -- ALGLIB --
     Copyright 09.12.2007 by Bochkanov Sergey
*************************************************************************/
void mlptrainlbfgs(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const double wstep, const ae_int_t maxits, ae_int_t &info, mlpreport &rep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlptrainlbfgs(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, decay, restarts, wstep, maxits, &info, const_cast<alglib_impl::mlpreport*>(rep.c_ptr()), &_alglib_env_state);
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
Neural network training using early stopping (base algorithm - L-BFGS with
regularization).

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry
    TrnXY       -   training set
    TrnSize     -   training set size, TrnSize>0
    ValXY       -   validation set
    ValSize     -   validation set size, ValSize>0
    Decay       -   weight decay constant, >=0.001
                    Decay term 'Decay*||Weights||^2' is added to error
                    function.
                    If you don't know what Decay to choose, use 0.001.
    Restarts    -   number of restarts, either:
                    * strictly positive number - algorithm make specified
                      number of restarts from random position.
                    * -1, in which case algorithm makes exactly one run
                      from the initial state of the network (no randomization).
                    If you don't know what Restarts to choose, choose one
                    one the following:
                    * -1 (deterministic start)
                    * +1 (one random restart)
                    * +5 (moderate amount of random restarts)

OUTPUT PARAMETERS:
    Network     -   trained neural network.
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NOut-1].
                    * -1, if wrong parameters specified
                          (NPoints<0, Restarts<1, ...).
                    *  2, task has been solved, stopping  criterion  met -
                          sufficiently small step size.  Not expected  (we
                          use  EARLY  stopping)  but  possible  and not an
                          error.
                    *  6, task has been solved, stopping  criterion  met -
                          increasing of validation set error.
    Rep         -   training report

NOTE:

Algorithm stops if validation set error increases for  a  long  enough  or
step size is small enought  (there  are  task  where  validation  set  may
decrease for eternity). In any case solution returned corresponds  to  the
minimum of validation set error.

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void mlptraines(const multilayerperceptron &network, const real_2d_array &trnxy, const ae_int_t trnsize, const real_2d_array &valxy, const ae_int_t valsize, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlptraines(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), const_cast<alglib_impl::ae_matrix*>(trnxy.c_ptr()), trnsize, const_cast<alglib_impl::ae_matrix*>(valxy.c_ptr()), valsize, decay, restarts, &info, const_cast<alglib_impl::mlpreport*>(rep.c_ptr()), &_alglib_env_state);
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
Cross-validation estimate of generalization error.

Base algorithm - L-BFGS.

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry.   Network is
                    not changed during cross-validation -  it is used only
                    as a representative of its architecture.
    XY          -   training set.
    SSize       -   training set size
    Decay       -   weight  decay, same as in MLPTrainLBFGS
    Restarts    -   number of restarts, >0.
                    restarts are counted for each partition separately, so
                    total number of restarts will be Restarts*FoldsCount.
    WStep       -   stopping criterion, same as in MLPTrainLBFGS
    MaxIts      -   stopping criterion, same as in MLPTrainLBFGS
    FoldsCount  -   number of folds in k-fold cross-validation,
                    2<=FoldsCount<=SSize.
                    recommended value: 10.

OUTPUT PARAMETERS:
    Info        -   return code, same as in MLPTrainLBFGS
    Rep         -   report, same as in MLPTrainLM/MLPTrainLBFGS
    CVRep       -   generalization error estimates

  -- ALGLIB --
     Copyright 09.12.2007 by Bochkanov Sergey
*************************************************************************/
void mlpkfoldcvlbfgs(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const double wstep, const ae_int_t maxits, const ae_int_t foldscount, ae_int_t &info, mlpreport &rep, mlpcvreport &cvrep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpkfoldcvlbfgs(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, decay, restarts, wstep, maxits, foldscount, &info, const_cast<alglib_impl::mlpreport*>(rep.c_ptr()), const_cast<alglib_impl::mlpcvreport*>(cvrep.c_ptr()), &_alglib_env_state);
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
Cross-validation estimate of generalization error.

Base algorithm - Levenberg-Marquardt.

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry.   Network is
                    not changed during cross-validation -  it is used only
                    as a representative of its architecture.
    XY          -   training set.
    SSize       -   training set size
    Decay       -   weight  decay, same as in MLPTrainLBFGS
    Restarts    -   number of restarts, >0.
                    restarts are counted for each partition separately, so
                    total number of restarts will be Restarts*FoldsCount.
    FoldsCount  -   number of folds in k-fold cross-validation,
                    2<=FoldsCount<=SSize.
                    recommended value: 10.

OUTPUT PARAMETERS:
    Info        -   return code, same as in MLPTrainLBFGS
    Rep         -   report, same as in MLPTrainLM/MLPTrainLBFGS
    CVRep       -   generalization error estimates

  -- ALGLIB --
     Copyright 09.12.2007 by Bochkanov Sergey
*************************************************************************/
void mlpkfoldcvlm(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const ae_int_t foldscount, ae_int_t &info, mlpreport &rep, mlpcvreport &cvrep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpkfoldcvlm(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, decay, restarts, foldscount, &info, const_cast<alglib_impl::mlpreport*>(rep.c_ptr()), const_cast<alglib_impl::mlpcvreport*>(cvrep.c_ptr()), &_alglib_env_state);
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
Neural networks ensemble
*************************************************************************/
_mlpensemble_owner::_mlpensemble_owner()
{
    p_struct = (alglib_impl::mlpensemble*)alglib_impl::ae_malloc(sizeof(alglib_impl::mlpensemble), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_mlpensemble_init(p_struct, NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_mlpensemble_owner::_mlpensemble_owner(const _mlpensemble_owner &rhs)
{
    p_struct = (alglib_impl::mlpensemble*)alglib_impl::ae_malloc(sizeof(alglib_impl::mlpensemble), NULL);
    if( p_struct==NULL )
        throw ap_error("ALGLIB: malloc error");
    if( !alglib_impl::_mlpensemble_init_copy(p_struct, const_cast<alglib_impl::mlpensemble*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
}

_mlpensemble_owner& _mlpensemble_owner::operator=(const _mlpensemble_owner &rhs)
{
    if( this==&rhs )
        return *this;
    alglib_impl::_mlpensemble_clear(p_struct);
    if( !alglib_impl::_mlpensemble_init_copy(p_struct, const_cast<alglib_impl::mlpensemble*>(rhs.p_struct), NULL, ae_false) )
        throw ap_error("ALGLIB: malloc error");
    return *this;
}

_mlpensemble_owner::~_mlpensemble_owner()
{
    alglib_impl::_mlpensemble_clear(p_struct);
    ae_free(p_struct);
}

alglib_impl::mlpensemble* _mlpensemble_owner::c_ptr()
{
    return p_struct;
}

alglib_impl::mlpensemble* _mlpensemble_owner::c_ptr() const
{
    return const_cast<alglib_impl::mlpensemble*>(p_struct);
}
mlpensemble::mlpensemble() : _mlpensemble_owner() 
{
}

mlpensemble::mlpensemble(const mlpensemble &rhs):_mlpensemble_owner(rhs) 
{
}

mlpensemble& mlpensemble::operator=(const mlpensemble &rhs)
{
    if( this==&rhs )
        return *this;
    _mlpensemble_owner::operator=(rhs);
    return *this;
}

mlpensemble::~mlpensemble()
{
}


/*************************************************************************
This function serializes data structure to string.

Important properties of s_out:
* it contains alphanumeric characters, dots, underscores, minus signs
* these symbols are grouped into words, which are separated by spaces
  and Windows-style (CR+LF) newlines
* although  serializer  uses  spaces and CR+LF as separators, you can 
  replace any separator character by arbitrary combination of spaces,
  tabs, Windows or Unix newlines. It allows flexible reformatting  of
  the  string  in  case you want to include it into text or XML file. 
  But you should not insert separators into the middle of the "words"
  nor you should change case of letters.
* s_out can be freely moved between 32-bit and 64-bit systems, little
  and big endian machines, and so on. You can serialize structure  on
  32-bit machine and unserialize it on 64-bit one (or vice versa), or
  serialize  it  on  SPARC  and  unserialize  on  x86.  You  can also 
  serialize  it  in  C++ version of ALGLIB and unserialize in C# one, 
  and vice versa.
*************************************************************************/
void mlpeserialize(mlpensemble &obj, std::string &s_out)
{
    alglib_impl::ae_state state;
    alglib_impl::ae_serializer serializer;
    alglib_impl::ae_int_t ssize;

    alglib_impl::ae_state_init(&state);
    try
    {
        alglib_impl::ae_serializer_init(&serializer);
        alglib_impl::ae_serializer_alloc_start(&serializer);
        alglib_impl::mlpealloc(&serializer, obj.c_ptr(), &state);
        ssize = alglib_impl::ae_serializer_get_alloc_size(&serializer);
        s_out.clear();
        s_out.reserve((size_t)(ssize+1));
        alglib_impl::ae_serializer_sstart_str(&serializer, &s_out);
        alglib_impl::mlpeserialize(&serializer, obj.c_ptr(), &state);
        alglib_impl::ae_serializer_stop(&serializer);
        if( s_out.length()>(size_t)ssize )
            throw ap_error("ALGLIB: serialization integrity error");
        alglib_impl::ae_serializer_clear(&serializer);
        alglib_impl::ae_state_clear(&state);
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}
/*************************************************************************
This function unserializes data structure from string.
*************************************************************************/
void mlpeunserialize(std::string &s_in, mlpensemble &obj)
{
    alglib_impl::ae_state state;
    alglib_impl::ae_serializer serializer;

    alglib_impl::ae_state_init(&state);
    try
    {
        alglib_impl::ae_serializer_init(&serializer);
        alglib_impl::ae_serializer_ustart_str(&serializer, &s_in);
        alglib_impl::mlpeunserialize(&serializer, obj.c_ptr(), &state);
        alglib_impl::ae_serializer_stop(&serializer);
        alglib_impl::ae_serializer_clear(&serializer);
        alglib_impl::ae_state_clear(&state);
    }
    catch(alglib_impl::ae_error_type)
    {
        throw ap_error(state.error_msg);
    }
    catch(...)
    {
        throw;
    }
}

/*************************************************************************
Like MLPCreate0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreate0(const ae_int_t nin, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpecreate0(nin, nout, ensemblesize, const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), &_alglib_env_state);
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
Like MLPCreate1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreate1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpecreate1(nin, nhid, nout, ensemblesize, const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), &_alglib_env_state);
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
Like MLPCreate2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreate2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpecreate2(nin, nhid1, nhid2, nout, ensemblesize, const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), &_alglib_env_state);
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
Like MLPCreateB0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreateb0(const ae_int_t nin, const ae_int_t nout, const double b, const double d, const ae_int_t ensemblesize, mlpensemble &ensemble)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpecreateb0(nin, nout, b, d, ensemblesize, const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), &_alglib_env_state);
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
Like MLPCreateB1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreateb1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double b, const double d, const ae_int_t ensemblesize, mlpensemble &ensemble)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpecreateb1(nin, nhid, nout, b, d, ensemblesize, const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), &_alglib_env_state);
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
Like MLPCreateB2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreateb2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double b, const double d, const ae_int_t ensemblesize, mlpensemble &ensemble)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpecreateb2(nin, nhid1, nhid2, nout, b, d, ensemblesize, const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), &_alglib_env_state);
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
Like MLPCreateR0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreater0(const ae_int_t nin, const ae_int_t nout, const double a, const double b, const ae_int_t ensemblesize, mlpensemble &ensemble)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpecreater0(nin, nout, a, b, ensemblesize, const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), &_alglib_env_state);
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
Like MLPCreateR1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreater1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double a, const double b, const ae_int_t ensemblesize, mlpensemble &ensemble)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpecreater1(nin, nhid, nout, a, b, ensemblesize, const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), &_alglib_env_state);
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
Like MLPCreateR2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreater2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double a, const double b, const ae_int_t ensemblesize, mlpensemble &ensemble)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpecreater2(nin, nhid1, nhid2, nout, a, b, ensemblesize, const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), &_alglib_env_state);
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
Like MLPCreateC0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatec0(const ae_int_t nin, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpecreatec0(nin, nout, ensemblesize, const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), &_alglib_env_state);
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
Like MLPCreateC1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatec1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpecreatec1(nin, nhid, nout, ensemblesize, const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), &_alglib_env_state);
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
Like MLPCreateC2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatec2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpecreatec2(nin, nhid1, nhid2, nout, ensemblesize, const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), &_alglib_env_state);
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
Creates ensemble from network. Only network geometry is copied.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatefromnetwork(const multilayerperceptron &network, const ae_int_t ensemblesize, mlpensemble &ensemble)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpecreatefromnetwork(const_cast<alglib_impl::multilayerperceptron*>(network.c_ptr()), ensemblesize, const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), &_alglib_env_state);
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
Randomization of MLP ensemble

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlperandomize(const mlpensemble &ensemble)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlperandomize(const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), &_alglib_env_state);
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
Return ensemble properties (number of inputs and outputs).

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpeproperties(const mlpensemble &ensemble, ae_int_t &nin, ae_int_t &nout)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpeproperties(const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), &nin, &nout, &_alglib_env_state);
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
Return normalization type (whether ensemble is SOFTMAX-normalized or not).

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
bool mlpeissoftmax(const mlpensemble &ensemble)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        ae_bool result = alglib_impl::mlpeissoftmax(const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), &_alglib_env_state);
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

/*************************************************************************
Procesing

INPUT PARAMETERS:
    Ensemble-   neural networks ensemble
    X       -   input vector,  array[0..NIn-1].
    Y       -   (possibly) preallocated buffer; if size of Y is less than
                NOut, it will be reallocated. If it is large enough, it
                is NOT reallocated, so we can save some time on reallocation.


OUTPUT PARAMETERS:
    Y       -   result. Regression estimate when solving regression  task,
                vector of posterior probabilities for classification task.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpeprocess(const mlpensemble &ensemble, const real_1d_array &x, real_1d_array &y)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpeprocess(const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::ae_vector*>(y.c_ptr()), &_alglib_env_state);
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
'interactive'  variant  of  MLPEProcess  for  languages  like Python which
support constructs like "Y = MLPEProcess(LM,X)" and interactive mode of the
interpreter

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpeprocessi(const mlpensemble &ensemble, const real_1d_array &x, real_1d_array &y)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpeprocessi(const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), const_cast<alglib_impl::ae_vector*>(x.c_ptr()), const_cast<alglib_impl::ae_vector*>(y.c_ptr()), &_alglib_env_state);
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
Relative classification error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    percent of incorrectly classified cases.
    Works both for classifier betwork and for regression networks which
are used as classifiers.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlperelclserror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::mlperelclserror(const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
Average cross-entropy (in bits per element) on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    CrossEntropy/(NPoints*LN(2)).
    Zero if ensemble solves regression task.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlpeavgce(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::mlpeavgce(const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
RMS error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.
    Its meaning for regression task is obvious. As for classification task
RMS error means error when estimating posterior probabilities.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlpermserror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::mlpermserror(const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
Average error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for classification task
it means average error when estimating posterior probabilities.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlpeavgerror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::mlpeavgerror(const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
Average relative error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for classification task
it means average relative error when estimating posterior probabilities.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlpeavgrelerror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        double result = alglib_impl::mlpeavgrelerror(const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, &_alglib_env_state);
        alglib_impl::ae_state_clear(&_alglib_env_state);
        return *(reinterpret_cast<double*>(&result));
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
Training neural networks ensemble using  bootstrap  aggregating (bagging).
Modified Levenberg-Marquardt algorithm is used as base training method.

INPUT PARAMETERS:
    Ensemble    -   model with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay coefficient, >=0.001
    Restarts    -   restarts, >0.

OUTPUT PARAMETERS:
    Ensemble    -   trained model
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report.
    OOBErrors   -   out-of-bag generalization error estimate

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpebagginglm(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep, mlpcvreport &ooberrors)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpebagginglm(const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, decay, restarts, &info, const_cast<alglib_impl::mlpreport*>(rep.c_ptr()), const_cast<alglib_impl::mlpcvreport*>(ooberrors.c_ptr()), &_alglib_env_state);
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
Training neural networks ensemble using  bootstrap  aggregating (bagging).
L-BFGS algorithm is used as base training method.

INPUT PARAMETERS:
    Ensemble    -   model with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay coefficient, >=0.001
    Restarts    -   restarts, >0.
    WStep       -   stopping criterion, same as in MLPTrainLBFGS
    MaxIts      -   stopping criterion, same as in MLPTrainLBFGS

OUTPUT PARAMETERS:
    Ensemble    -   trained model
    Info        -   return code:
                    * -8, if both WStep=0 and MaxIts=0
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report.
    OOBErrors   -   out-of-bag generalization error estimate

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpebagginglbfgs(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const double wstep, const ae_int_t maxits, ae_int_t &info, mlpreport &rep, mlpcvreport &ooberrors)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpebagginglbfgs(const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, decay, restarts, wstep, maxits, &info, const_cast<alglib_impl::mlpreport*>(rep.c_ptr()), const_cast<alglib_impl::mlpcvreport*>(ooberrors.c_ptr()), &_alglib_env_state);
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
Training neural networks ensemble using early stopping.

INPUT PARAMETERS:
    Ensemble    -   model with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay coefficient, >=0.001
    Restarts    -   restarts, >0.

OUTPUT PARAMETERS:
    Ensemble    -   trained model
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<0, Restarts<1).
                    *  6, if task has been solved.
    Rep         -   training report.
    OOBErrors   -   out-of-bag generalization error estimate

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void mlpetraines(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::mlpetraines(const_cast<alglib_impl::mlpensemble*>(ensemble.c_ptr()), const_cast<alglib_impl::ae_matrix*>(xy.c_ptr()), npoints, decay, restarts, &info, const_cast<alglib_impl::mlpreport*>(rep.c_ptr()), &_alglib_env_state);
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
Principal components analysis

Subroutine  builds  orthogonal  basis  where  first  axis  corresponds  to
direction with maximum variance, second axis maximizes variance in subspace
orthogonal to first axis and so on.

It should be noted that, unlike LDA, PCA does not use class labels.

INPUT PARAMETERS:
    X           -   dataset, array[0..NPoints-1,0..NVars-1].
                    matrix contains ONLY INDEPENDENT VARIABLES.
    NPoints     -   dataset size, NPoints>=0
    NVars       -   number of independent variables, NVars>=1

¬€’ŒƒÕ€≈ œ¿–¿Ã≈“–€:
    Info        -   return code:
                    * -4, if SVD subroutine haven't converged
                    * -1, if wrong parameters has been passed (NPoints<0,
                          NVars<1)
                    *  1, if task is solved
    S2          -   array[0..NVars-1]. variance values corresponding
                    to basis vectors.
    V           -   array[0..NVars-1,0..NVars-1]
                    matrix, whose columns store basis vectors.

  -- ALGLIB --
     Copyright 25.08.2008 by Bochkanov Sergey
*************************************************************************/
void pcabuildbasis(const real_2d_array &x, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, real_1d_array &s2, real_2d_array &v)
{
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    try
    {
        alglib_impl::pcabuildbasis(const_cast<alglib_impl::ae_matrix*>(x.c_ptr()), npoints, nvars, &info, const_cast<alglib_impl::ae_vector*>(s2.c_ptr()), const_cast<alglib_impl::ae_matrix*>(v.c_ptr()), &_alglib_env_state);
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
static double bdss_xlny(double x, double y, ae_state *_state);
static double bdss_getcv(/* Integer */ ae_vector* cnt,
     ae_int_t nc,
     ae_state *_state);
static void bdss_tieaddc(/* Integer */ ae_vector* c,
     /* Integer */ ae_vector* ties,
     ae_int_t ntie,
     ae_int_t nc,
     /* Integer */ ae_vector* cnt,
     ae_state *_state);
static void bdss_tiesubc(/* Integer */ ae_vector* c,
     /* Integer */ ae_vector* ties,
     ae_int_t ntie,
     ae_int_t nc,
     /* Integer */ ae_vector* cnt,
     ae_state *_state);


static ae_int_t dforest_innernodewidth = 3;
static ae_int_t dforest_leafnodewidth = 2;
static ae_int_t dforest_dfusestrongsplits = 1;
static ae_int_t dforest_dfuseevs = 2;
static ae_int_t dforest_dffirstversion = 0;
static ae_int_t dforest_dfclserror(decisionforest* df,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
static void dforest_dfprocessinternal(decisionforest* df,
     ae_int_t offs,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
static void dforest_dfbuildtree(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t nfeatures,
     ae_int_t nvarsinpool,
     ae_int_t flags,
     dfinternalbuffers* bufs,
     ae_state *_state);
static void dforest_dfbuildtreerec(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t nfeatures,
     ae_int_t nvarsinpool,
     ae_int_t flags,
     ae_int_t* numprocessed,
     ae_int_t idx1,
     ae_int_t idx2,
     dfinternalbuffers* bufs,
     ae_state *_state);
static void dforest_dfsplitc(/* Real    */ ae_vector* x,
     /* Integer */ ae_vector* c,
     /* Integer */ ae_vector* cntbuf,
     ae_int_t n,
     ae_int_t nc,
     ae_int_t flags,
     ae_int_t* info,
     double* threshold,
     double* e,
     /* Real    */ ae_vector* sortrbuf,
     /* Integer */ ae_vector* sortibuf,
     ae_state *_state);
static void dforest_dfsplitr(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_int_t flags,
     ae_int_t* info,
     double* threshold,
     double* e,
     /* Real    */ ae_vector* sortrbuf,
     /* Real    */ ae_vector* sortrbuf2,
     ae_state *_state);


static ae_int_t linreg_lrvnum = 5;
static void linreg_lrinternal(/* Real    */ ae_matrix* xy,
     /* Real    */ ae_vector* s,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     linearmodel* lm,
     lrreport* ar,
     ae_state *_state);




static ae_bool kmeans_selectcenterpp(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     /* Real    */ ae_matrix* centers,
     /* Boolean */ ae_vector* busycenters,
     ae_int_t ccnt,
     /* Real    */ ae_vector* d2,
     /* Real    */ ae_vector* p,
     /* Real    */ ae_vector* tmp,
     ae_state *_state);




static ae_int_t mlpbase_mlpvnum = 7;
static ae_int_t mlpbase_mlpfirstversion = 0;
static ae_int_t mlpbase_nfieldwidth = 4;
static ae_int_t mlpbase_hlconnfieldwidth = 5;
static ae_int_t mlpbase_hlnfieldwidth = 4;
static ae_int_t mlpbase_chunksize = 32;
static void mlpbase_addinputlayer(ae_int_t ncount,
     /* Integer */ ae_vector* lsizes,
     /* Integer */ ae_vector* ltypes,
     /* Integer */ ae_vector* lconnfirst,
     /* Integer */ ae_vector* lconnlast,
     ae_int_t* lastproc,
     ae_state *_state);
static void mlpbase_addbiasedsummatorlayer(ae_int_t ncount,
     /* Integer */ ae_vector* lsizes,
     /* Integer */ ae_vector* ltypes,
     /* Integer */ ae_vector* lconnfirst,
     /* Integer */ ae_vector* lconnlast,
     ae_int_t* lastproc,
     ae_state *_state);
static void mlpbase_addactivationlayer(ae_int_t functype,
     /* Integer */ ae_vector* lsizes,
     /* Integer */ ae_vector* ltypes,
     /* Integer */ ae_vector* lconnfirst,
     /* Integer */ ae_vector* lconnlast,
     ae_int_t* lastproc,
     ae_state *_state);
static void mlpbase_addzerolayer(/* Integer */ ae_vector* lsizes,
     /* Integer */ ae_vector* ltypes,
     /* Integer */ ae_vector* lconnfirst,
     /* Integer */ ae_vector* lconnlast,
     ae_int_t* lastproc,
     ae_state *_state);
static void mlpbase_hladdinputlayer(multilayerperceptron* network,
     ae_int_t* connidx,
     ae_int_t* neuroidx,
     ae_int_t* structinfoidx,
     ae_int_t nin,
     ae_state *_state);
static void mlpbase_hladdoutputlayer(multilayerperceptron* network,
     ae_int_t* connidx,
     ae_int_t* neuroidx,
     ae_int_t* structinfoidx,
     ae_int_t* weightsidx,
     ae_int_t k,
     ae_int_t nprev,
     ae_int_t nout,
     ae_bool iscls,
     ae_bool islinearout,
     ae_state *_state);
static void mlpbase_hladdhiddenlayer(multilayerperceptron* network,
     ae_int_t* connidx,
     ae_int_t* neuroidx,
     ae_int_t* structinfoidx,
     ae_int_t* weightsidx,
     ae_int_t k,
     ae_int_t nprev,
     ae_int_t ncur,
     ae_state *_state);
static void mlpbase_fillhighlevelinformation(multilayerperceptron* network,
     ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_bool iscls,
     ae_bool islinearout,
     ae_state *_state);
static void mlpbase_mlpcreate(ae_int_t nin,
     ae_int_t nout,
     /* Integer */ ae_vector* lsizes,
     /* Integer */ ae_vector* ltypes,
     /* Integer */ ae_vector* lconnfirst,
     /* Integer */ ae_vector* lconnlast,
     ae_int_t layerscount,
     ae_bool isclsnet,
     multilayerperceptron* network,
     ae_state *_state);
static void mlpbase_mlphessianbatchinternal(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     ae_bool naturalerr,
     double* e,
     /* Real    */ ae_vector* grad,
     /* Real    */ ae_matrix* h,
     ae_state *_state);
static void mlpbase_mlpinternalcalculategradient(multilayerperceptron* network,
     /* Real    */ ae_vector* neurons,
     /* Real    */ ae_vector* weights,
     /* Real    */ ae_vector* derror,
     /* Real    */ ae_vector* grad,
     ae_bool naturalerrorfunc,
     ae_state *_state);
static void mlpbase_mlpchunkedgradient(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t cstart,
     ae_int_t csize,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_bool naturalerrorfunc,
     ae_state *_state);
static double mlpbase_safecrossentropy(double t,
     double z,
     ae_state *_state);


static double logit_xtol = 100*ae_machineepsilon;
static double logit_ftol = 0.0001;
static double logit_gtol = 0.3;
static ae_int_t logit_maxfev = 20;
static double logit_stpmin = 1.0E-2;
static double logit_stpmax = 1.0E5;
static ae_int_t logit_logitvnum = 6;
static void logit_mnliexp(/* Real    */ ae_vector* w,
     /* Real    */ ae_vector* x,
     ae_state *_state);
static void logit_mnlallerrors(logitmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double* relcls,
     double* avgce,
     double* rms,
     double* avg,
     double* avgrel,
     ae_state *_state);
static void logit_mnlmcsrch(ae_int_t n,
     /* Real    */ ae_vector* x,
     double* f,
     /* Real    */ ae_vector* g,
     /* Real    */ ae_vector* s,
     double* stp,
     ae_int_t* info,
     ae_int_t* nfev,
     /* Real    */ ae_vector* wa,
     logitmcstate* state,
     ae_int_t* stage,
     ae_state *_state);
static void logit_mnlmcstep(double* stx,
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


static double mcpd_xtol = 1.0E-8;
static void mcpd_mcpdinit(ae_int_t n,
     ae_int_t entrystate,
     ae_int_t exitstate,
     mcpdstate* s,
     ae_state *_state);


static double mlptrain_mindecay = 0.001;
static void mlptrain_mlpkfoldcvgeneral(multilayerperceptron* n,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     ae_int_t foldscount,
     ae_bool lmalgorithm,
     double wstep,
     ae_int_t maxits,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* cvrep,
     ae_state *_state);
static void mlptrain_mlpkfoldsplit(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nclasses,
     ae_int_t foldscount,
     ae_bool stratifiedsplits,
     /* Integer */ ae_vector* folds,
     ae_state *_state);


//static ae_int_t mlpe_mlpntotaloffset = 3;
//static ae_int_t mlpe_mlpevnum = 9;
static ae_int_t mlpe_mlpefirstversion = 1;
static void mlpe_mlpeallerrors(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double* relcls,
     double* avgce,
     double* rms,
     double* avg,
     double* avgrel,
     ae_state *_state);
static void mlpe_mlpebagginginternal(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     double wstep,
     ae_int_t maxits,
     ae_bool lmalgorithm,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* ooberrors,
     ae_state *_state);







/*************************************************************************
This set of routines (DSErrAllocate, DSErrAccumulate, DSErrFinish)
calculates different error functions (classification error, cross-entropy,
rms, avg, avg.rel errors).

1. DSErrAllocate prepares buffer.
2. DSErrAccumulate accumulates individual errors:
    * Y contains predicted output (posterior probabilities for classification)
    * DesiredY contains desired output (class number for classification)
3. DSErrFinish outputs results:
   * Buf[0] contains relative classification error (zero for regression tasks)
   * Buf[1] contains avg. cross-entropy (zero for regression tasks)
   * Buf[2] contains rms error (regression, classification)
   * Buf[3] contains average error (regression, classification)
   * Buf[4] contains average relative error (regression, classification)
   
NOTES(1):
    "NClasses>0" means that we have classification task.
    "NClasses<0" means regression task with -NClasses real outputs.

NOTES(2):
    rms. avg, avg.rel errors for classification tasks are interpreted as
    errors in posterior probabilities with respect to probabilities given
    by training/test set.

  -- ALGLIB --
     Copyright 11.01.2009 by Bochkanov Sergey
*************************************************************************/
void dserrallocate(ae_int_t nclasses,
     /* Real    */ ae_vector* buf,
     ae_state *_state)
{

    ae_vector_clear(buf);

    ae_vector_set_length(buf, 7+1, _state);
    buf->ptr.p_double[0] = 0;
    buf->ptr.p_double[1] = 0;
    buf->ptr.p_double[2] = 0;
    buf->ptr.p_double[3] = 0;
    buf->ptr.p_double[4] = 0;
    buf->ptr.p_double[5] = nclasses;
    buf->ptr.p_double[6] = 0;
    buf->ptr.p_double[7] = 0;
}


/*************************************************************************
See DSErrAllocate for comments on this routine.

  -- ALGLIB --
     Copyright 11.01.2009 by Bochkanov Sergey
*************************************************************************/
void dserraccumulate(/* Real    */ ae_vector* buf,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* desiredy,
     ae_state *_state)
{
    ae_int_t nclasses;
    ae_int_t nout;
    ae_int_t offs;
    ae_int_t mmax;
    ae_int_t rmax;
    ae_int_t j;
    double v;
    double ev;


    offs = 5;
    nclasses = ae_round(buf->ptr.p_double[offs], _state);
    if( nclasses>0 )
    {
        
        /*
         * Classification
         */
        rmax = ae_round(desiredy->ptr.p_double[0], _state);
        mmax = 0;
        for(j=1; j<=nclasses-1; j++)
        {
            if( ae_fp_greater(y->ptr.p_double[j],y->ptr.p_double[mmax]) )
            {
                mmax = j;
            }
        }
        if( mmax!=rmax )
        {
            buf->ptr.p_double[0] = buf->ptr.p_double[0]+1;
        }
        if( ae_fp_greater(y->ptr.p_double[rmax],0) )
        {
            buf->ptr.p_double[1] = buf->ptr.p_double[1]-ae_log(y->ptr.p_double[rmax], _state);
        }
        else
        {
            buf->ptr.p_double[1] = buf->ptr.p_double[1]+ae_log(ae_maxrealnumber, _state);
        }
        for(j=0; j<=nclasses-1; j++)
        {
            v = y->ptr.p_double[j];
            if( j==rmax )
            {
                ev = 1;
            }
            else
            {
                ev = 0;
            }
            buf->ptr.p_double[2] = buf->ptr.p_double[2]+ae_sqr(v-ev, _state);
            buf->ptr.p_double[3] = buf->ptr.p_double[3]+ae_fabs(v-ev, _state);
            if( ae_fp_neq(ev,0) )
            {
                buf->ptr.p_double[4] = buf->ptr.p_double[4]+ae_fabs((v-ev)/ev, _state);
                buf->ptr.p_double[offs+2] = buf->ptr.p_double[offs+2]+1;
            }
        }
        buf->ptr.p_double[offs+1] = buf->ptr.p_double[offs+1]+1;
    }
    else
    {
        
        /*
         * Regression
         */
        nout = -nclasses;
        rmax = 0;
        for(j=1; j<=nout-1; j++)
        {
            if( ae_fp_greater(desiredy->ptr.p_double[j],desiredy->ptr.p_double[rmax]) )
            {
                rmax = j;
            }
        }
        mmax = 0;
        for(j=1; j<=nout-1; j++)
        {
            if( ae_fp_greater(y->ptr.p_double[j],y->ptr.p_double[mmax]) )
            {
                mmax = j;
            }
        }
        if( mmax!=rmax )
        {
            buf->ptr.p_double[0] = buf->ptr.p_double[0]+1;
        }
        for(j=0; j<=nout-1; j++)
        {
            v = y->ptr.p_double[j];
            ev = desiredy->ptr.p_double[j];
            buf->ptr.p_double[2] = buf->ptr.p_double[2]+ae_sqr(v-ev, _state);
            buf->ptr.p_double[3] = buf->ptr.p_double[3]+ae_fabs(v-ev, _state);
            if( ae_fp_neq(ev,0) )
            {
                buf->ptr.p_double[4] = buf->ptr.p_double[4]+ae_fabs((v-ev)/ev, _state);
                buf->ptr.p_double[offs+2] = buf->ptr.p_double[offs+2]+1;
            }
        }
        buf->ptr.p_double[offs+1] = buf->ptr.p_double[offs+1]+1;
    }
}


/*************************************************************************
See DSErrAllocate for comments on this routine.

  -- ALGLIB --
     Copyright 11.01.2009 by Bochkanov Sergey
*************************************************************************/
void dserrfinish(/* Real    */ ae_vector* buf, ae_state *_state)
{
    ae_int_t nout;
    ae_int_t offs;


    offs = 5;
    nout = ae_iabs(ae_round(buf->ptr.p_double[offs], _state), _state);
    if( ae_fp_neq(buf->ptr.p_double[offs+1],0) )
    {
        buf->ptr.p_double[0] = buf->ptr.p_double[0]/buf->ptr.p_double[offs+1];
        buf->ptr.p_double[1] = buf->ptr.p_double[1]/buf->ptr.p_double[offs+1];
        buf->ptr.p_double[2] = ae_sqrt(buf->ptr.p_double[2]/(nout*buf->ptr.p_double[offs+1]), _state);
        buf->ptr.p_double[3] = buf->ptr.p_double[3]/(nout*buf->ptr.p_double[offs+1]);
    }
    if( ae_fp_neq(buf->ptr.p_double[offs+2],0) )
    {
        buf->ptr.p_double[4] = buf->ptr.p_double[4]/buf->ptr.p_double[offs+2];
    }
}


/*************************************************************************

  -- ALGLIB --
     Copyright 19.05.2008 by Bochkanov Sergey
*************************************************************************/
void dsnormalize(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     /* Real    */ ae_vector* means,
     /* Real    */ ae_vector* sigmas,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_vector tmp;
    double mean;
    double variance;
    double skewness;
    double kurtosis;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    ae_vector_clear(means);
    ae_vector_clear(sigmas);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);

    
    /*
     * Test parameters
     */
    if( npoints<=0||nvars<1 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    *info = 1;
    
    /*
     * Standartization
     */
    ae_vector_set_length(means, nvars-1+1, _state);
    ae_vector_set_length(sigmas, nvars-1+1, _state);
    ae_vector_set_length(&tmp, npoints-1+1, _state);
    for(j=0; j<=nvars-1; j++)
    {
        ae_v_move(&tmp.ptr.p_double[0], 1, &xy->ptr.pp_double[0][j], xy->stride, ae_v_len(0,npoints-1));
        samplemoments(&tmp, npoints, &mean, &variance, &skewness, &kurtosis, _state);
        means->ptr.p_double[j] = mean;
        sigmas->ptr.p_double[j] = ae_sqrt(variance, _state);
        if( ae_fp_eq(sigmas->ptr.p_double[j],0) )
        {
            sigmas->ptr.p_double[j] = 1;
        }
        for(i=0; i<=npoints-1; i++)
        {
            xy->ptr.pp_double[i][j] = (xy->ptr.pp_double[i][j]-means->ptr.p_double[j])/sigmas->ptr.p_double[j];
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************

  -- ALGLIB --
     Copyright 19.05.2008 by Bochkanov Sergey
*************************************************************************/
void dsnormalizec(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     /* Real    */ ae_vector* means,
     /* Real    */ ae_vector* sigmas,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t j;
    ae_vector tmp;
    double mean;
    double variance;
    double skewness;
    double kurtosis;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    ae_vector_clear(means);
    ae_vector_clear(sigmas);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);

    
    /*
     * Test parameters
     */
    if( npoints<=0||nvars<1 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    *info = 1;
    
    /*
     * Standartization
     */
    ae_vector_set_length(means, nvars-1+1, _state);
    ae_vector_set_length(sigmas, nvars-1+1, _state);
    ae_vector_set_length(&tmp, npoints-1+1, _state);
    for(j=0; j<=nvars-1; j++)
    {
        ae_v_move(&tmp.ptr.p_double[0], 1, &xy->ptr.pp_double[0][j], xy->stride, ae_v_len(0,npoints-1));
        samplemoments(&tmp, npoints, &mean, &variance, &skewness, &kurtosis, _state);
        means->ptr.p_double[j] = mean;
        sigmas->ptr.p_double[j] = ae_sqrt(variance, _state);
        if( ae_fp_eq(sigmas->ptr.p_double[j],0) )
        {
            sigmas->ptr.p_double[j] = 1;
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************

  -- ALGLIB --
     Copyright 19.05.2008 by Bochkanov Sergey
*************************************************************************/
double dsgetmeanmindistance(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_vector tmp;
    ae_vector tmp2;
    double v;
    double result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmp2, 0, DT_REAL, _state, ae_true);

    
    /*
     * Test parameters
     */
    if( npoints<=0||nvars<1 )
    {
        result = 0;
        ae_frame_leave(_state);
        return result;
    }
    
    /*
     * Process
     */
    ae_vector_set_length(&tmp, npoints-1+1, _state);
    for(i=0; i<=npoints-1; i++)
    {
        tmp.ptr.p_double[i] = ae_maxrealnumber;
    }
    ae_vector_set_length(&tmp2, nvars-1+1, _state);
    for(i=0; i<=npoints-1; i++)
    {
        for(j=i+1; j<=npoints-1; j++)
        {
            ae_v_move(&tmp2.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
            ae_v_sub(&tmp2.ptr.p_double[0], 1, &xy->ptr.pp_double[j][0], 1, ae_v_len(0,nvars-1));
            v = ae_v_dotproduct(&tmp2.ptr.p_double[0], 1, &tmp2.ptr.p_double[0], 1, ae_v_len(0,nvars-1));
            v = ae_sqrt(v, _state);
            tmp.ptr.p_double[i] = ae_minreal(tmp.ptr.p_double[i], v, _state);
            tmp.ptr.p_double[j] = ae_minreal(tmp.ptr.p_double[j], v, _state);
        }
    }
    result = 0;
    for(i=0; i<=npoints-1; i++)
    {
        result = result+tmp.ptr.p_double[i]/npoints;
    }
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************

  -- ALGLIB --
     Copyright 19.05.2008 by Bochkanov Sergey
*************************************************************************/
void dstie(/* Real    */ ae_vector* a,
     ae_int_t n,
     /* Integer */ ae_vector* ties,
     ae_int_t* tiecount,
     /* Integer */ ae_vector* p1,
     /* Integer */ ae_vector* p2,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t k;
    ae_vector tmp;

    ae_frame_make(_state, &_frame_block);
    ae_vector_clear(ties);
    *tiecount = 0;
    ae_vector_clear(p1);
    ae_vector_clear(p2);
    ae_vector_init(&tmp, 0, DT_INT, _state, ae_true);

    
    /*
     * Special case
     */
    if( n<=0 )
    {
        *tiecount = 0;
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Sort A
     */
    tagsort(a, n, p1, p2, _state);
    
    /*
     * Process ties
     */
    *tiecount = 1;
    for(i=1; i<=n-1; i++)
    {
        if( ae_fp_neq(a->ptr.p_double[i],a->ptr.p_double[i-1]) )
        {
            *tiecount = *tiecount+1;
        }
    }
    ae_vector_set_length(ties, *tiecount+1, _state);
    ties->ptr.p_int[0] = 0;
    k = 1;
    for(i=1; i<=n-1; i++)
    {
        if( ae_fp_neq(a->ptr.p_double[i],a->ptr.p_double[i-1]) )
        {
            ties->ptr.p_int[k] = i;
            k = k+1;
        }
    }
    ties->ptr.p_int[*tiecount] = n;
    ae_frame_leave(_state);
}


/*************************************************************************

  -- ALGLIB --
     Copyright 11.12.2008 by Bochkanov Sergey
*************************************************************************/
void dstiefasti(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     ae_int_t n,
     /* Integer */ ae_vector* ties,
     ae_int_t* tiecount,
     /* Real    */ ae_vector* bufr,
     /* Integer */ ae_vector* bufi,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t k;
    ae_vector tmp;

    ae_frame_make(_state, &_frame_block);
    *tiecount = 0;
    ae_vector_init(&tmp, 0, DT_INT, _state, ae_true);

    
    /*
     * Special case
     */
    if( n<=0 )
    {
        *tiecount = 0;
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Sort A
     */
    tagsortfasti(a, b, bufr, bufi, n, _state);
    
    /*
     * Process ties
     */
    ties->ptr.p_int[0] = 0;
    k = 1;
    for(i=1; i<=n-1; i++)
    {
        if( ae_fp_neq(a->ptr.p_double[i],a->ptr.p_double[i-1]) )
        {
            ties->ptr.p_int[k] = i;
            k = k+1;
        }
    }
    ties->ptr.p_int[k] = n;
    *tiecount = k;
    ae_frame_leave(_state);
}


/*************************************************************************
Optimal binary classification

Algorithms finds optimal (=with minimal cross-entropy) binary partition.
Internal subroutine.

INPUT PARAMETERS:
    A       -   array[0..N-1], variable
    C       -   array[0..N-1], class numbers (0 or 1).
    N       -   array size

OUTPUT PARAMETERS:
    Info    -   completetion code:
                * -3, all values of A[] are same (partition is impossible)
                * -2, one of C[] is incorrect (<0, >1)
                * -1, incorrect pararemets were passed (N<=0).
                *  1, OK
    Threshold-  partiton boundary. Left part contains values which are
                strictly less than Threshold. Right part contains values
                which are greater than or equal to Threshold.
    PAL, PBL-   probabilities P(0|v<Threshold) and P(1|v<Threshold)
    PAR, PBR-   probabilities P(0|v>=Threshold) and P(1|v>=Threshold)
    CVE     -   cross-validation estimate of cross-entropy

  -- ALGLIB --
     Copyright 22.05.2008 by Bochkanov Sergey
*************************************************************************/
void dsoptimalsplit2(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* c,
     ae_int_t n,
     ae_int_t* info,
     double* threshold,
     double* pal,
     double* pbl,
     double* par,
     double* pbr,
     double* cve,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector _a;
    ae_vector _c;
    ae_int_t i;
    ae_int_t t;
    double s;
    ae_vector ties;
    ae_int_t tiecount;
    ae_vector p1;
    ae_vector p2;
    ae_int_t k;
    ae_int_t koptimal;
    double pak;
    double pbk;
    double cvoptimal;
    double cv;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init_copy(&_a, a, _state, ae_true);
    a = &_a;
    ae_vector_init_copy(&_c, c, _state, ae_true);
    c = &_c;
    *info = 0;
    *threshold = 0;
    *pal = 0;
    *pbl = 0;
    *par = 0;
    *pbr = 0;
    *cve = 0;
    ae_vector_init(&ties, 0, DT_INT, _state, ae_true);
    ae_vector_init(&p1, 0, DT_INT, _state, ae_true);
    ae_vector_init(&p2, 0, DT_INT, _state, ae_true);

    
    /*
     * Test for errors in inputs
     */
    if( n<=0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    for(i=0; i<=n-1; i++)
    {
        if( c->ptr.p_int[i]!=0&&c->ptr.p_int[i]!=1 )
        {
            *info = -2;
            ae_frame_leave(_state);
            return;
        }
    }
    *info = 1;
    
    /*
     * Tie
     */
    dstie(a, n, &ties, &tiecount, &p1, &p2, _state);
    for(i=0; i<=n-1; i++)
    {
        if( p2.ptr.p_int[i]!=i )
        {
            t = c->ptr.p_int[i];
            c->ptr.p_int[i] = c->ptr.p_int[p2.ptr.p_int[i]];
            c->ptr.p_int[p2.ptr.p_int[i]] = t;
        }
    }
    
    /*
     * Special case: number of ties is 1.
     *
     * NOTE: we assume that P[i,j] equals to 0 or 1,
     *       intermediate values are not allowed.
     */
    if( tiecount==1 )
    {
        *info = -3;
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * General case, number of ties > 1
     *
     * NOTE: we assume that P[i,j] equals to 0 or 1,
     *       intermediate values are not allowed.
     */
    *pal = 0;
    *pbl = 0;
    *par = 0;
    *pbr = 0;
    for(i=0; i<=n-1; i++)
    {
        if( c->ptr.p_int[i]==0 )
        {
            *par = *par+1;
        }
        if( c->ptr.p_int[i]==1 )
        {
            *pbr = *pbr+1;
        }
    }
    koptimal = -1;
    cvoptimal = ae_maxrealnumber;
    for(k=0; k<=tiecount-2; k++)
    {
        
        /*
         * first, obtain information about K-th tie which is
         * moved from R-part to L-part
         */
        pak = 0;
        pbk = 0;
        for(i=ties.ptr.p_int[k]; i<=ties.ptr.p_int[k+1]-1; i++)
        {
            if( c->ptr.p_int[i]==0 )
            {
                pak = pak+1;
            }
            if( c->ptr.p_int[i]==1 )
            {
                pbk = pbk+1;
            }
        }
        
        /*
         * Calculate cross-validation CE
         */
        cv = 0;
        cv = cv-bdss_xlny(*pal+pak, (*pal+pak)/(*pal+pak+(*pbl)+pbk+1), _state);
        cv = cv-bdss_xlny(*pbl+pbk, (*pbl+pbk)/(*pal+pak+1+(*pbl)+pbk), _state);
        cv = cv-bdss_xlny(*par-pak, (*par-pak)/(*par-pak+(*pbr)-pbk+1), _state);
        cv = cv-bdss_xlny(*pbr-pbk, (*pbr-pbk)/(*par-pak+1+(*pbr)-pbk), _state);
        
        /*
         * Compare with best
         */
        if( ae_fp_less(cv,cvoptimal) )
        {
            cvoptimal = cv;
            koptimal = k;
        }
        
        /*
         * update
         */
        *pal = *pal+pak;
        *pbl = *pbl+pbk;
        *par = *par-pak;
        *pbr = *pbr-pbk;
    }
    *cve = cvoptimal;
    *threshold = 0.5*(a->ptr.p_double[ties.ptr.p_int[koptimal]]+a->ptr.p_double[ties.ptr.p_int[koptimal+1]]);
    *pal = 0;
    *pbl = 0;
    *par = 0;
    *pbr = 0;
    for(i=0; i<=n-1; i++)
    {
        if( ae_fp_less(a->ptr.p_double[i],*threshold) )
        {
            if( c->ptr.p_int[i]==0 )
            {
                *pal = *pal+1;
            }
            else
            {
                *pbl = *pbl+1;
            }
        }
        else
        {
            if( c->ptr.p_int[i]==0 )
            {
                *par = *par+1;
            }
            else
            {
                *pbr = *pbr+1;
            }
        }
    }
    s = *pal+(*pbl);
    *pal = *pal/s;
    *pbl = *pbl/s;
    s = *par+(*pbr);
    *par = *par/s;
    *pbr = *pbr/s;
    ae_frame_leave(_state);
}


/*************************************************************************
Optimal partition, internal subroutine. Fast version.

Accepts:
    A       array[0..N-1]       array of attributes     array[0..N-1]
    C       array[0..N-1]       array of class labels
    TiesBuf array[0..N]         temporaries (ties)
    CntBuf  array[0..2*NC-1]    temporaries (counts)
    Alpha                       centering factor (0<=alpha<=1, recommended value - 0.05)
    BufR    array[0..N-1]       temporaries
    BufI    array[0..N-1]       temporaries

Output:
    Info    error code (">0"=OK, "<0"=bad)
    RMS     training set RMS error
    CVRMS   leave-one-out RMS error
    
Note:
    content of all arrays is changed by subroutine;
    it doesn't allocate temporaries.

  -- ALGLIB --
     Copyright 11.12.2008 by Bochkanov Sergey
*************************************************************************/
void dsoptimalsplit2fast(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* c,
     /* Integer */ ae_vector* tiesbuf,
     /* Integer */ ae_vector* cntbuf,
     /* Real    */ ae_vector* bufr,
     /* Integer */ ae_vector* bufi,
     ae_int_t n,
     ae_int_t nc,
     double alpha,
     ae_int_t* info,
     double* threshold,
     double* rms,
     double* cvrms,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t k;
    ae_int_t cl;
    ae_int_t tiecount;
    double cbest;
    double cc;
    ae_int_t koptimal;
    ae_int_t sl;
    ae_int_t sr;
    double v;
    double w;
    double x;

    *info = 0;
    *threshold = 0;
    *rms = 0;
    *cvrms = 0;

    
    /*
     * Test for errors in inputs
     */
    if( n<=0||nc<2 )
    {
        *info = -1;
        return;
    }
    for(i=0; i<=n-1; i++)
    {
        if( c->ptr.p_int[i]<0||c->ptr.p_int[i]>=nc )
        {
            *info = -2;
            return;
        }
    }
    *info = 1;
    
    /*
     * Tie
     */
    dstiefasti(a, c, n, tiesbuf, &tiecount, bufr, bufi, _state);
    
    /*
     * Special case: number of ties is 1.
     */
    if( tiecount==1 )
    {
        *info = -3;
        return;
    }
    
    /*
     * General case, number of ties > 1
     */
    for(i=0; i<=2*nc-1; i++)
    {
        cntbuf->ptr.p_int[i] = 0;
    }
    for(i=0; i<=n-1; i++)
    {
        cntbuf->ptr.p_int[nc+c->ptr.p_int[i]] = cntbuf->ptr.p_int[nc+c->ptr.p_int[i]]+1;
    }
    koptimal = -1;
    *threshold = a->ptr.p_double[n-1];
    cbest = ae_maxrealnumber;
    sl = 0;
    sr = n;
    for(k=0; k<=tiecount-2; k++)
    {
        
        /*
         * first, move Kth tie from right to left
         */
        for(i=tiesbuf->ptr.p_int[k]; i<=tiesbuf->ptr.p_int[k+1]-1; i++)
        {
            cl = c->ptr.p_int[i];
            cntbuf->ptr.p_int[cl] = cntbuf->ptr.p_int[cl]+1;
            cntbuf->ptr.p_int[nc+cl] = cntbuf->ptr.p_int[nc+cl]-1;
        }
        sl = sl+(tiesbuf->ptr.p_int[k+1]-tiesbuf->ptr.p_int[k]);
        sr = sr-(tiesbuf->ptr.p_int[k+1]-tiesbuf->ptr.p_int[k]);
        
        /*
         * Calculate RMS error
         */
        v = 0;
        for(i=0; i<=nc-1; i++)
        {
            w = cntbuf->ptr.p_int[i];
            v = v+w*ae_sqr(w/sl-1, _state);
            v = v+(sl-w)*ae_sqr(w/sl, _state);
            w = cntbuf->ptr.p_int[nc+i];
            v = v+w*ae_sqr(w/sr-1, _state);
            v = v+(sr-w)*ae_sqr(w/sr, _state);
        }
        v = ae_sqrt(v/(nc*n), _state);
        
        /*
         * Compare with best
         */
        x = (double)(2*sl)/(double)(sl+sr)-1;
        cc = v*(1-alpha+alpha*ae_sqr(x, _state));
        if( ae_fp_less(cc,cbest) )
        {
            
            /*
             * store split
             */
            *rms = v;
            koptimal = k;
            cbest = cc;
            
            /*
             * calculate CVRMS error
             */
            *cvrms = 0;
            for(i=0; i<=nc-1; i++)
            {
                if( sl>1 )
                {
                    w = cntbuf->ptr.p_int[i];
                    *cvrms = *cvrms+w*ae_sqr((w-1)/(sl-1)-1, _state);
                    *cvrms = *cvrms+(sl-w)*ae_sqr(w/(sl-1), _state);
                }
                else
                {
                    w = cntbuf->ptr.p_int[i];
                    *cvrms = *cvrms+w*ae_sqr((double)1/(double)nc-1, _state);
                    *cvrms = *cvrms+(sl-w)*ae_sqr((double)1/(double)nc, _state);
                }
                if( sr>1 )
                {
                    w = cntbuf->ptr.p_int[nc+i];
                    *cvrms = *cvrms+w*ae_sqr((w-1)/(sr-1)-1, _state);
                    *cvrms = *cvrms+(sr-w)*ae_sqr(w/(sr-1), _state);
                }
                else
                {
                    w = cntbuf->ptr.p_int[nc+i];
                    *cvrms = *cvrms+w*ae_sqr((double)1/(double)nc-1, _state);
                    *cvrms = *cvrms+(sr-w)*ae_sqr((double)1/(double)nc, _state);
                }
            }
            *cvrms = ae_sqrt(*cvrms/(nc*n), _state);
        }
    }
    
    /*
     * Calculate threshold.
     * Code is a bit complicated because there can be such
     * numbers that 0.5(A+B) equals to A or B (if A-B=epsilon)
     */
    *threshold = 0.5*(a->ptr.p_double[tiesbuf->ptr.p_int[koptimal]]+a->ptr.p_double[tiesbuf->ptr.p_int[koptimal+1]]);
    if( ae_fp_less_eq(*threshold,a->ptr.p_double[tiesbuf->ptr.p_int[koptimal]]) )
    {
        *threshold = a->ptr.p_double[tiesbuf->ptr.p_int[koptimal+1]];
    }
}


/*************************************************************************
Automatic non-optimal discretization, internal subroutine.

  -- ALGLIB --
     Copyright 22.05.2008 by Bochkanov Sergey
*************************************************************************/
void dssplitk(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* c,
     ae_int_t n,
     ae_int_t nc,
     ae_int_t kmax,
     ae_int_t* info,
     /* Real    */ ae_vector* thresholds,
     ae_int_t* ni,
     double* cve,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector _a;
    ae_vector _c;
    ae_int_t i;
    ae_int_t j;
    ae_int_t j1;
    ae_int_t k;
    ae_vector ties;
    ae_int_t tiecount;
    ae_vector p1;
    ae_vector p2;
    ae_vector cnt;
    double v2;
    ae_int_t bestk;
    double bestcve;
    ae_vector bestsizes;
    double curcve;
    ae_vector cursizes;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init_copy(&_a, a, _state, ae_true);
    a = &_a;
    ae_vector_init_copy(&_c, c, _state, ae_true);
    c = &_c;
    *info = 0;
    ae_vector_clear(thresholds);
    *ni = 0;
    *cve = 0;
    ae_vector_init(&ties, 0, DT_INT, _state, ae_true);
    ae_vector_init(&p1, 0, DT_INT, _state, ae_true);
    ae_vector_init(&p2, 0, DT_INT, _state, ae_true);
    ae_vector_init(&cnt, 0, DT_INT, _state, ae_true);
    ae_vector_init(&bestsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&cursizes, 0, DT_INT, _state, ae_true);

    
    /*
     * Test for errors in inputs
     */
    if( (n<=0||nc<2)||kmax<2 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    for(i=0; i<=n-1; i++)
    {
        if( c->ptr.p_int[i]<0||c->ptr.p_int[i]>=nc )
        {
            *info = -2;
            ae_frame_leave(_state);
            return;
        }
    }
    *info = 1;
    
    /*
     * Tie
     */
    dstie(a, n, &ties, &tiecount, &p1, &p2, _state);
    for(i=0; i<=n-1; i++)
    {
        if( p2.ptr.p_int[i]!=i )
        {
            k = c->ptr.p_int[i];
            c->ptr.p_int[i] = c->ptr.p_int[p2.ptr.p_int[i]];
            c->ptr.p_int[p2.ptr.p_int[i]] = k;
        }
    }
    
    /*
     * Special cases
     */
    if( tiecount==1 )
    {
        *info = -3;
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * General case:
     * 0. allocate arrays
     */
    kmax = ae_minint(kmax, tiecount, _state);
    ae_vector_set_length(&bestsizes, kmax-1+1, _state);
    ae_vector_set_length(&cursizes, kmax-1+1, _state);
    ae_vector_set_length(&cnt, nc-1+1, _state);
    
    /*
     * General case:
     * 1. prepare "weak" solution (two subintervals, divided at median)
     */
    v2 = ae_maxrealnumber;
    j = -1;
    for(i=1; i<=tiecount-1; i++)
    {
        if( ae_fp_less(ae_fabs(ties.ptr.p_int[i]-0.5*(n-1), _state),v2) )
        {
            v2 = ae_fabs(ties.ptr.p_int[i]-0.5*n, _state);
            j = i;
        }
    }
    ae_assert(j>0, "DSSplitK: internal error #1!", _state);
    bestk = 2;
    bestsizes.ptr.p_int[0] = ties.ptr.p_int[j];
    bestsizes.ptr.p_int[1] = n-j;
    bestcve = 0;
    for(i=0; i<=nc-1; i++)
    {
        cnt.ptr.p_int[i] = 0;
    }
    for(i=0; i<=j-1; i++)
    {
        bdss_tieaddc(c, &ties, i, nc, &cnt, _state);
    }
    bestcve = bestcve+bdss_getcv(&cnt, nc, _state);
    for(i=0; i<=nc-1; i++)
    {
        cnt.ptr.p_int[i] = 0;
    }
    for(i=j; i<=tiecount-1; i++)
    {
        bdss_tieaddc(c, &ties, i, nc, &cnt, _state);
    }
    bestcve = bestcve+bdss_getcv(&cnt, nc, _state);
    
    /*
     * General case:
     * 2. Use greedy algorithm to find sub-optimal split in O(KMax*N) time
     */
    for(k=2; k<=kmax; k++)
    {
        
        /*
         * Prepare greedy K-interval split
         */
        for(i=0; i<=k-1; i++)
        {
            cursizes.ptr.p_int[i] = 0;
        }
        i = 0;
        j = 0;
        while(j<=tiecount-1&&i<=k-1)
        {
            
            /*
             * Rule: I-th bin is empty, fill it
             */
            if( cursizes.ptr.p_int[i]==0 )
            {
                cursizes.ptr.p_int[i] = ties.ptr.p_int[j+1]-ties.ptr.p_int[j];
                j = j+1;
                continue;
            }
            
            /*
             * Rule: (K-1-I) bins left, (K-1-I) ties left (1 tie per bin); next bin
             */
            if( tiecount-j==k-1-i )
            {
                i = i+1;
                continue;
            }
            
            /*
             * Rule: last bin, always place in current
             */
            if( i==k-1 )
            {
                cursizes.ptr.p_int[i] = cursizes.ptr.p_int[i]+ties.ptr.p_int[j+1]-ties.ptr.p_int[j];
                j = j+1;
                continue;
            }
            
            /*
             * Place J-th tie in I-th bin, or leave for I+1-th bin.
             */
            if( ae_fp_less(ae_fabs(cursizes.ptr.p_int[i]+ties.ptr.p_int[j+1]-ties.ptr.p_int[j]-(double)n/(double)k, _state),ae_fabs(cursizes.ptr.p_int[i]-(double)n/(double)k, _state)) )
            {
                cursizes.ptr.p_int[i] = cursizes.ptr.p_int[i]+ties.ptr.p_int[j+1]-ties.ptr.p_int[j];
                j = j+1;
            }
            else
            {
                i = i+1;
            }
        }
        ae_assert(cursizes.ptr.p_int[k-1]!=0&&j==tiecount, "DSSplitK: internal error #1", _state);
        
        /*
         * Calculate CVE
         */
        curcve = 0;
        j = 0;
        for(i=0; i<=k-1; i++)
        {
            for(j1=0; j1<=nc-1; j1++)
            {
                cnt.ptr.p_int[j1] = 0;
            }
            for(j1=j; j1<=j+cursizes.ptr.p_int[i]-1; j1++)
            {
                cnt.ptr.p_int[c->ptr.p_int[j1]] = cnt.ptr.p_int[c->ptr.p_int[j1]]+1;
            }
            curcve = curcve+bdss_getcv(&cnt, nc, _state);
            j = j+cursizes.ptr.p_int[i];
        }
        
        /*
         * Choose best variant
         */
        if( ae_fp_less(curcve,bestcve) )
        {
            for(i=0; i<=k-1; i++)
            {
                bestsizes.ptr.p_int[i] = cursizes.ptr.p_int[i];
            }
            bestcve = curcve;
            bestk = k;
        }
    }
    
    /*
     * Transform from sizes to thresholds
     */
    *cve = bestcve;
    *ni = bestk;
    ae_vector_set_length(thresholds, *ni-2+1, _state);
    j = bestsizes.ptr.p_int[0];
    for(i=1; i<=bestk-1; i++)
    {
        thresholds->ptr.p_double[i-1] = 0.5*(a->ptr.p_double[j-1]+a->ptr.p_double[j]);
        j = j+bestsizes.ptr.p_int[i];
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Automatic optimal discretization, internal subroutine.

  -- ALGLIB --
     Copyright 22.05.2008 by Bochkanov Sergey
*************************************************************************/
void dsoptimalsplitk(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* c,
     ae_int_t n,
     ae_int_t nc,
     ae_int_t kmax,
     ae_int_t* info,
     /* Real    */ ae_vector* thresholds,
     ae_int_t* ni,
     double* cve,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector _a;
    ae_vector _c;
    ae_int_t i;
    ae_int_t j;
    ae_int_t s;
    ae_int_t jl;
    ae_int_t jr;
    double v2;
    ae_vector ties;
    ae_int_t tiecount;
    ae_vector p1;
    ae_vector p2;
    double cvtemp;
    ae_vector cnt;
    ae_vector cnt2;
    ae_matrix cv;
    ae_matrix splits;
    ae_int_t k;
    ae_int_t koptimal;
    double cvoptimal;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init_copy(&_a, a, _state, ae_true);
    a = &_a;
    ae_vector_init_copy(&_c, c, _state, ae_true);
    c = &_c;
    *info = 0;
    ae_vector_clear(thresholds);
    *ni = 0;
    *cve = 0;
    ae_vector_init(&ties, 0, DT_INT, _state, ae_true);
    ae_vector_init(&p1, 0, DT_INT, _state, ae_true);
    ae_vector_init(&p2, 0, DT_INT, _state, ae_true);
    ae_vector_init(&cnt, 0, DT_INT, _state, ae_true);
    ae_vector_init(&cnt2, 0, DT_INT, _state, ae_true);
    ae_matrix_init(&cv, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&splits, 0, 0, DT_INT, _state, ae_true);

    
    /*
     * Test for errors in inputs
     */
    if( (n<=0||nc<2)||kmax<2 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    for(i=0; i<=n-1; i++)
    {
        if( c->ptr.p_int[i]<0||c->ptr.p_int[i]>=nc )
        {
            *info = -2;
            ae_frame_leave(_state);
            return;
        }
    }
    *info = 1;
    
    /*
     * Tie
     */
    dstie(a, n, &ties, &tiecount, &p1, &p2, _state);
    for(i=0; i<=n-1; i++)
    {
        if( p2.ptr.p_int[i]!=i )
        {
            k = c->ptr.p_int[i];
            c->ptr.p_int[i] = c->ptr.p_int[p2.ptr.p_int[i]];
            c->ptr.p_int[p2.ptr.p_int[i]] = k;
        }
    }
    
    /*
     * Special cases
     */
    if( tiecount==1 )
    {
        *info = -3;
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * General case
     * Use dynamic programming to find best split in O(KMax*NC*TieCount^2) time
     */
    kmax = ae_minint(kmax, tiecount, _state);
    ae_matrix_set_length(&cv, kmax-1+1, tiecount-1+1, _state);
    ae_matrix_set_length(&splits, kmax-1+1, tiecount-1+1, _state);
    ae_vector_set_length(&cnt, nc-1+1, _state);
    ae_vector_set_length(&cnt2, nc-1+1, _state);
    for(j=0; j<=nc-1; j++)
    {
        cnt.ptr.p_int[j] = 0;
    }
    for(j=0; j<=tiecount-1; j++)
    {
        bdss_tieaddc(c, &ties, j, nc, &cnt, _state);
        splits.ptr.pp_int[0][j] = 0;
        cv.ptr.pp_double[0][j] = bdss_getcv(&cnt, nc, _state);
    }
    for(k=1; k<=kmax-1; k++)
    {
        for(j=0; j<=nc-1; j++)
        {
            cnt.ptr.p_int[j] = 0;
        }
        
        /*
         * Subtask size J in [K..TieCount-1]:
         * optimal K-splitting on ties from 0-th to J-th.
         */
        for(j=k; j<=tiecount-1; j++)
        {
            
            /*
             * Update Cnt - let it contain classes of ties from K-th to J-th
             */
            bdss_tieaddc(c, &ties, j, nc, &cnt, _state);
            
            /*
             * Search for optimal split point S in [K..J]
             */
            for(i=0; i<=nc-1; i++)
            {
                cnt2.ptr.p_int[i] = cnt.ptr.p_int[i];
            }
            cv.ptr.pp_double[k][j] = cv.ptr.pp_double[k-1][j-1]+bdss_getcv(&cnt2, nc, _state);
            splits.ptr.pp_int[k][j] = j;
            for(s=k+1; s<=j; s++)
            {
                
                /*
                 * Update Cnt2 - let it contain classes of ties from S-th to J-th
                 */
                bdss_tiesubc(c, &ties, s-1, nc, &cnt2, _state);
                
                /*
                 * Calculate CVE
                 */
                cvtemp = cv.ptr.pp_double[k-1][s-1]+bdss_getcv(&cnt2, nc, _state);
                if( ae_fp_less(cvtemp,cv.ptr.pp_double[k][j]) )
                {
                    cv.ptr.pp_double[k][j] = cvtemp;
                    splits.ptr.pp_int[k][j] = s;
                }
            }
        }
    }
    
    /*
     * Choose best partition, output result
     */
    koptimal = -1;
    cvoptimal = ae_maxrealnumber;
    for(k=0; k<=kmax-1; k++)
    {
        if( ae_fp_less(cv.ptr.pp_double[k][tiecount-1],cvoptimal) )
        {
            cvoptimal = cv.ptr.pp_double[k][tiecount-1];
            koptimal = k;
        }
    }
    ae_assert(koptimal>=0, "DSOptimalSplitK: internal error #1!", _state);
    if( koptimal==0 )
    {
        
        /*
         * Special case: best partition is one big interval.
         * Even 2-partition is not better.
         * This is possible when dealing with "weak" predictor variables.
         *
         * Make binary split as close to the median as possible.
         */
        v2 = ae_maxrealnumber;
        j = -1;
        for(i=1; i<=tiecount-1; i++)
        {
            if( ae_fp_less(ae_fabs(ties.ptr.p_int[i]-0.5*(n-1), _state),v2) )
            {
                v2 = ae_fabs(ties.ptr.p_int[i]-0.5*(n-1), _state);
                j = i;
            }
        }
        ae_assert(j>0, "DSOptimalSplitK: internal error #2!", _state);
        ae_vector_set_length(thresholds, 0+1, _state);
        thresholds->ptr.p_double[0] = 0.5*(a->ptr.p_double[ties.ptr.p_int[j-1]]+a->ptr.p_double[ties.ptr.p_int[j]]);
        *ni = 2;
        *cve = 0;
        for(i=0; i<=nc-1; i++)
        {
            cnt.ptr.p_int[i] = 0;
        }
        for(i=0; i<=j-1; i++)
        {
            bdss_tieaddc(c, &ties, i, nc, &cnt, _state);
        }
        *cve = *cve+bdss_getcv(&cnt, nc, _state);
        for(i=0; i<=nc-1; i++)
        {
            cnt.ptr.p_int[i] = 0;
        }
        for(i=j; i<=tiecount-1; i++)
        {
            bdss_tieaddc(c, &ties, i, nc, &cnt, _state);
        }
        *cve = *cve+bdss_getcv(&cnt, nc, _state);
    }
    else
    {
        
        /*
         * General case: 2 or more intervals
         */
        ae_vector_set_length(thresholds, koptimal-1+1, _state);
        *ni = koptimal+1;
        *cve = cv.ptr.pp_double[koptimal][tiecount-1];
        jl = splits.ptr.pp_int[koptimal][tiecount-1];
        jr = tiecount-1;
        for(k=koptimal; k>=1; k--)
        {
            thresholds->ptr.p_double[k-1] = 0.5*(a->ptr.p_double[ties.ptr.p_int[jl-1]]+a->ptr.p_double[ties.ptr.p_int[jl]]);
            jr = jl-1;
            jl = splits.ptr.pp_int[k-1][jl-1];
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Internal function
*************************************************************************/
static double bdss_xlny(double x, double y, ae_state *_state)
{
    double result;


    if( ae_fp_eq(x,0) )
    {
        result = 0;
    }
    else
    {
        result = x*ae_log(y, _state);
    }
    return result;
}


/*************************************************************************
Internal function,
returns number of samples of class I in Cnt[I]
*************************************************************************/
static double bdss_getcv(/* Integer */ ae_vector* cnt,
     ae_int_t nc,
     ae_state *_state)
{
    ae_int_t i;
    double s;
    double result;


    s = 0;
    for(i=0; i<=nc-1; i++)
    {
        s = s+cnt->ptr.p_int[i];
    }
    result = 0;
    for(i=0; i<=nc-1; i++)
    {
        result = result-bdss_xlny(cnt->ptr.p_int[i], cnt->ptr.p_int[i]/(s+nc-1), _state);
    }
    return result;
}


/*************************************************************************
Internal function, adds number of samples of class I in tie NTie to Cnt[I]
*************************************************************************/
static void bdss_tieaddc(/* Integer */ ae_vector* c,
     /* Integer */ ae_vector* ties,
     ae_int_t ntie,
     ae_int_t , //nc,
     /* Integer */ ae_vector* cnt,
     ae_state *) //_state)
{
    ae_int_t i;


    for(i=ties->ptr.p_int[ntie]; i<=ties->ptr.p_int[ntie+1]-1; i++)
    {
        cnt->ptr.p_int[c->ptr.p_int[i]] = cnt->ptr.p_int[c->ptr.p_int[i]]+1;
    }
}


/*************************************************************************
Internal function, subtracts number of samples of class I in tie NTie to Cnt[I]
*************************************************************************/
static void bdss_tiesubc(/* Integer */ ae_vector* c,
     /* Integer */ ae_vector* ties,
     ae_int_t ntie,
     ae_int_t , //nc,
     /* Integer */ ae_vector* cnt,
     ae_state *) //_state)
{
    ae_int_t i;


    for(i=ties->ptr.p_int[ntie]; i<=ties->ptr.p_int[ntie+1]-1; i++)
    {
        cnt->ptr.p_int[c->ptr.p_int[i]] = cnt->ptr.p_int[c->ptr.p_int[i]]-1;
    }
}


ae_bool _cvreport_init(cvreport* /*p*/, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    return ae_true;
}


ae_bool _cvreport_init_copy(cvreport* dst, cvreport* src, ae_state */*_state*/, ae_bool) // make_automatic)
{
    dst->relclserror = src->relclserror;
    dst->avgce = src->avgce;
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
    return ae_true;
}


void _cvreport_clear(cvreport* ) //	p)
{
}




/*************************************************************************
This subroutine builds random decision forest.

INPUT PARAMETERS:
    XY          -   training set
    NPoints     -   training set size, NPoints>=1
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   task type:
                    * NClasses=1 - regression task with one
                                   dependent variable
                    * NClasses>1 - classification task with
                                   NClasses classes.
    NTrees      -   number of trees in a forest, NTrees>=1.
                    recommended values: 50-100.
    R           -   percent of a training set used to build
                    individual trees. 0<R<=1.
                    recommended values: 0.1 <= R <= 0.66.

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<1, NVars<1, NClasses<1, NTrees<1, R<=0
                          or R>1).
                    *  1, if task has been solved
    DF          -   model built
    Rep         -   training report, contains error on a training set
                    and out-of-bag estimates of generalization error.

  -- ALGLIB --
     Copyright 19.02.2009 by Bochkanov Sergey
*************************************************************************/
void dfbuildrandomdecisionforest(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t ntrees,
     double r,
     ae_int_t* info,
     decisionforest* df,
     dfreport* rep,
     ae_state *_state)
{
    ae_int_t samplesize;

    *info = 0;
    _decisionforest_clear(df);
    _dfreport_clear(rep);

    if( ae_fp_less_eq(r,0)||ae_fp_greater(r,1) )
    {
        *info = -1;
        return;
    }
    samplesize = ae_maxint(ae_round(r*npoints, _state), 1, _state);
    dfbuildinternal(xy, npoints, nvars, nclasses, ntrees, samplesize, ae_maxint(nvars/2, 1, _state), dforest_dfusestrongsplits+dforest_dfuseevs, info, df, rep, _state);
}


/*************************************************************************
This subroutine builds random decision forest.
This function gives ability to tune number of variables used when choosing
best split.

INPUT PARAMETERS:
    XY          -   training set
    NPoints     -   training set size, NPoints>=1
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   task type:
                    * NClasses=1 - regression task with one
                                   dependent variable
                    * NClasses>1 - classification task with
                                   NClasses classes.
    NTrees      -   number of trees in a forest, NTrees>=1.
                    recommended values: 50-100.
    NRndVars    -   number of variables used when choosing best split
    R           -   percent of a training set used to build
                    individual trees. 0<R<=1.
                    recommended values: 0.1 <= R <= 0.66.

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<1, NVars<1, NClasses<1, NTrees<1, R<=0
                          or R>1).
                    *  1, if task has been solved
    DF          -   model built
    Rep         -   training report, contains error on a training set
                    and out-of-bag estimates of generalization error.

  -- ALGLIB --
     Copyright 19.02.2009 by Bochkanov Sergey
*************************************************************************/
void dfbuildrandomdecisionforestx1(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t ntrees,
     ae_int_t nrndvars,
     double r,
     ae_int_t* info,
     decisionforest* df,
     dfreport* rep,
     ae_state *_state)
{
    ae_int_t samplesize;

    *info = 0;
    _decisionforest_clear(df);
    _dfreport_clear(rep);

    if( ae_fp_less_eq(r,0)||ae_fp_greater(r,1) )
    {
        *info = -1;
        return;
    }
    if( nrndvars<=0||nrndvars>nvars )
    {
        *info = -1;
        return;
    }
    samplesize = ae_maxint(ae_round(r*npoints, _state), 1, _state);
    dfbuildinternal(xy, npoints, nvars, nclasses, ntrees, samplesize, nrndvars, dforest_dfusestrongsplits+dforest_dfuseevs, info, df, rep, _state);
}


void dfbuildinternal(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t ntrees,
     ae_int_t samplesize,
     ae_int_t nfeatures,
     ae_int_t flags,
     ae_int_t* info,
     decisionforest* df,
     dfreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t tmpi;
    ae_int_t lasttreeoffs;
    ae_int_t offs;
    ae_int_t ooboffs;
    ae_int_t treesize;
    ae_int_t nvarsinpool;
    ae_bool useevs;
    dfinternalbuffers bufs;
    ae_vector permbuf;
    ae_vector oobbuf;
    ae_vector oobcntbuf;
    ae_matrix xys;
    ae_vector x;
    ae_vector y;
    ae_int_t oobcnt;
    ae_int_t oobrelcnt;
    double v;
    double vmin;
    double vmax;
    ae_bool bflag;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _decisionforest_clear(df);
    _dfreport_clear(rep);
    _dfinternalbuffers_init(&bufs, _state, ae_true);
    ae_vector_init(&permbuf, 0, DT_INT, _state, ae_true);
    ae_vector_init(&oobbuf, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&oobcntbuf, 0, DT_INT, _state, ae_true);
    ae_matrix_init(&xys, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);

    
    /*
     * Test for inputs
     */
    if( (((((npoints<1||samplesize<1)||samplesize>npoints)||nvars<1)||nclasses<1)||ntrees<1)||nfeatures<1 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    if( nclasses>1 )
    {
        for(i=0; i<=npoints-1; i++)
        {
            if( ae_round(xy->ptr.pp_double[i][nvars], _state)<0||ae_round(xy->ptr.pp_double[i][nvars], _state)>=nclasses )
            {
                *info = -2;
                ae_frame_leave(_state);
                return;
            }
        }
    }
    *info = 1;
    
    /*
     * Flags
     */
    useevs = flags/dforest_dfuseevs%2!=0;
    
    /*
     * Allocate data, prepare header
     */
    treesize = 1+dforest_innernodewidth*(samplesize-1)+dforest_leafnodewidth*samplesize;
    ae_vector_set_length(&permbuf, npoints-1+1, _state);
    ae_vector_set_length(&bufs.treebuf, treesize-1+1, _state);
    ae_vector_set_length(&bufs.idxbuf, npoints-1+1, _state);
    ae_vector_set_length(&bufs.tmpbufr, npoints-1+1, _state);
    ae_vector_set_length(&bufs.tmpbufr2, npoints-1+1, _state);
    ae_vector_set_length(&bufs.tmpbufi, npoints-1+1, _state);
    ae_vector_set_length(&bufs.sortrbuf, npoints, _state);
    ae_vector_set_length(&bufs.sortrbuf2, npoints, _state);
    ae_vector_set_length(&bufs.sortibuf, npoints, _state);
    ae_vector_set_length(&bufs.varpool, nvars-1+1, _state);
    ae_vector_set_length(&bufs.evsbin, nvars-1+1, _state);
    ae_vector_set_length(&bufs.evssplits, nvars-1+1, _state);
    ae_vector_set_length(&bufs.classibuf, 2*nclasses-1+1, _state);
    ae_vector_set_length(&oobbuf, nclasses*npoints-1+1, _state);
    ae_vector_set_length(&oobcntbuf, npoints-1+1, _state);
    ae_vector_set_length(&df->trees, ntrees*treesize-1+1, _state);
    ae_matrix_set_length(&xys, samplesize-1+1, nvars+1, _state);
    ae_vector_set_length(&x, nvars-1+1, _state);
    ae_vector_set_length(&y, nclasses-1+1, _state);
    for(i=0; i<=npoints-1; i++)
    {
        permbuf.ptr.p_int[i] = i;
    }
    for(i=0; i<=npoints*nclasses-1; i++)
    {
        oobbuf.ptr.p_double[i] = 0;
    }
    for(i=0; i<=npoints-1; i++)
    {
        oobcntbuf.ptr.p_int[i] = 0;
    }
    
    /*
     * Prepare variable pool and EVS (extended variable selection/splitting) buffers
     * (whether EVS is turned on or not):
     * 1. detect binary variables and pre-calculate splits for them
     * 2. detect variables with non-distinct values and exclude them from pool
     */
    for(i=0; i<=nvars-1; i++)
    {
        bufs.varpool.ptr.p_int[i] = i;
    }
    nvarsinpool = nvars;
    if( useevs )
    {
        for(j=0; j<=nvars-1; j++)
        {
            vmin = xy->ptr.pp_double[0][j];
            vmax = vmin;
            for(i=0; i<=npoints-1; i++)
            {
                v = xy->ptr.pp_double[i][j];
                vmin = ae_minreal(vmin, v, _state);
                vmax = ae_maxreal(vmax, v, _state);
            }
            if( ae_fp_eq(vmin,vmax) )
            {
                
                /*
                 * exclude variable from pool
                 */
                bufs.varpool.ptr.p_int[j] = bufs.varpool.ptr.p_int[nvarsinpool-1];
                bufs.varpool.ptr.p_int[nvarsinpool-1] = -1;
                nvarsinpool = nvarsinpool-1;
                continue;
            }
            bflag = ae_false;
            for(i=0; i<=npoints-1; i++)
            {
                v = xy->ptr.pp_double[i][j];
                if( ae_fp_neq(v,vmin)&&ae_fp_neq(v,vmax) )
                {
                    bflag = ae_true;
                    break;
                }
            }
            if( bflag )
            {
                
                /*
                 * non-binary variable
                 */
                bufs.evsbin.ptr.p_bool[j] = ae_false;
            }
            else
            {
                
                /*
                 * Prepare
                 */
                bufs.evsbin.ptr.p_bool[j] = ae_true;
                bufs.evssplits.ptr.p_double[j] = 0.5*(vmin+vmax);
                if( ae_fp_less_eq(bufs.evssplits.ptr.p_double[j],vmin) )
                {
                    bufs.evssplits.ptr.p_double[j] = vmax;
                }
            }
        }
    }
    
    /*
     * RANDOM FOREST FORMAT
     * W[0]         -   size of array
     * W[1]         -   version number
     * W[2]         -   NVars
     * W[3]         -   NClasses (1 for regression)
     * W[4]         -   NTrees
     * W[5]         -   trees offset
     *
     *
     * TREE FORMAT
     * W[Offs]      -   size of sub-array
     *     node info:
     * W[K+0]       -   variable number        (-1 for leaf mode)
     * W[K+1]       -   threshold              (class/value for leaf node)
     * W[K+2]       -   ">=" branch index      (absent for leaf node)
     *
     */
    df->nvars = nvars;
    df->nclasses = nclasses;
    df->ntrees = ntrees;
    
    /*
     * Build forest
     */
    offs = 0;
    for(i=0; i<=ntrees-1; i++)
    {
        
        /*
         * Prepare sample
         */
        for(k=0; k<=samplesize-1; k++)
        {
            j = k+ae_randominteger(npoints-k, _state);
            tmpi = permbuf.ptr.p_int[k];
            permbuf.ptr.p_int[k] = permbuf.ptr.p_int[j];
            permbuf.ptr.p_int[j] = tmpi;
            j = permbuf.ptr.p_int[k];
            ae_v_move(&xys.ptr.pp_double[k][0], 1, &xy->ptr.pp_double[j][0], 1, ae_v_len(0,nvars));
        }
        
        /*
         * build tree, copy
         */
        dforest_dfbuildtree(&xys, samplesize, nvars, nclasses, nfeatures, nvarsinpool, flags, &bufs, _state);
        j = ae_round(bufs.treebuf.ptr.p_double[0], _state);
        ae_v_move(&df->trees.ptr.p_double[offs], 1, &bufs.treebuf.ptr.p_double[0], 1, ae_v_len(offs,offs+j-1));
        lasttreeoffs = offs;
        offs = offs+j;
        
        /*
         * OOB estimates
         */
        for(k=samplesize; k<=npoints-1; k++)
        {
            for(j=0; j<=nclasses-1; j++)
            {
                y.ptr.p_double[j] = 0;
            }
            j = permbuf.ptr.p_int[k];
            ae_v_move(&x.ptr.p_double[0], 1, &xy->ptr.pp_double[j][0], 1, ae_v_len(0,nvars-1));
            dforest_dfprocessinternal(df, lasttreeoffs, &x, &y, _state);
            ae_v_add(&oobbuf.ptr.p_double[j*nclasses], 1, &y.ptr.p_double[0], 1, ae_v_len(j*nclasses,(j+1)*nclasses-1));
            oobcntbuf.ptr.p_int[j] = oobcntbuf.ptr.p_int[j]+1;
        }
    }
    df->bufsize = offs;
    
    /*
     * Normalize OOB results
     */
    for(i=0; i<=npoints-1; i++)
    {
        if( oobcntbuf.ptr.p_int[i]!=0 )
        {
            v = (double)1/(double)oobcntbuf.ptr.p_int[i];
            ae_v_muld(&oobbuf.ptr.p_double[i*nclasses], 1, ae_v_len(i*nclasses,i*nclasses+nclasses-1), v);
        }
    }
    
    /*
     * Calculate training set estimates
     */
    rep->relclserror = dfrelclserror(df, xy, npoints, _state);
    rep->avgce = dfavgce(df, xy, npoints, _state);
    rep->rmserror = dfrmserror(df, xy, npoints, _state);
    rep->avgerror = dfavgerror(df, xy, npoints, _state);
    rep->avgrelerror = dfavgrelerror(df, xy, npoints, _state);
    
    /*
     * Calculate OOB estimates.
     */
    rep->oobrelclserror = 0;
    rep->oobavgce = 0;
    rep->oobrmserror = 0;
    rep->oobavgerror = 0;
    rep->oobavgrelerror = 0;
    oobcnt = 0;
    oobrelcnt = 0;
    for(i=0; i<=npoints-1; i++)
    {
        if( oobcntbuf.ptr.p_int[i]!=0 )
        {
            ooboffs = i*nclasses;
            if( nclasses>1 )
            {
                
                /*
                 * classification-specific code
                 */
                k = ae_round(xy->ptr.pp_double[i][nvars], _state);
                tmpi = 0;
                for(j=1; j<=nclasses-1; j++)
                {
                    if( ae_fp_greater(oobbuf.ptr.p_double[ooboffs+j],oobbuf.ptr.p_double[ooboffs+tmpi]) )
                    {
                        tmpi = j;
                    }
                }
                if( tmpi!=k )
                {
                    rep->oobrelclserror = rep->oobrelclserror+1;
                }
                if( ae_fp_neq(oobbuf.ptr.p_double[ooboffs+k],0) )
                {
                    rep->oobavgce = rep->oobavgce-ae_log(oobbuf.ptr.p_double[ooboffs+k], _state);
                }
                else
                {
                    rep->oobavgce = rep->oobavgce-ae_log(ae_minrealnumber, _state);
                }
                for(j=0; j<=nclasses-1; j++)
                {
                    if( j==k )
                    {
                        rep->oobrmserror = rep->oobrmserror+ae_sqr(oobbuf.ptr.p_double[ooboffs+j]-1, _state);
                        rep->oobavgerror = rep->oobavgerror+ae_fabs(oobbuf.ptr.p_double[ooboffs+j]-1, _state);
                        rep->oobavgrelerror = rep->oobavgrelerror+ae_fabs(oobbuf.ptr.p_double[ooboffs+j]-1, _state);
                        oobrelcnt = oobrelcnt+1;
                    }
                    else
                    {
                        rep->oobrmserror = rep->oobrmserror+ae_sqr(oobbuf.ptr.p_double[ooboffs+j], _state);
                        rep->oobavgerror = rep->oobavgerror+ae_fabs(oobbuf.ptr.p_double[ooboffs+j], _state);
                    }
                }
            }
            else
            {
                
                /*
                 * regression-specific code
                 */
                rep->oobrmserror = rep->oobrmserror+ae_sqr(oobbuf.ptr.p_double[ooboffs]-xy->ptr.pp_double[i][nvars], _state);
                rep->oobavgerror = rep->oobavgerror+ae_fabs(oobbuf.ptr.p_double[ooboffs]-xy->ptr.pp_double[i][nvars], _state);
                if( ae_fp_neq(xy->ptr.pp_double[i][nvars],0) )
                {
                    rep->oobavgrelerror = rep->oobavgrelerror+ae_fabs((oobbuf.ptr.p_double[ooboffs]-xy->ptr.pp_double[i][nvars])/xy->ptr.pp_double[i][nvars], _state);
                    oobrelcnt = oobrelcnt+1;
                }
            }
            
            /*
             * update OOB estimates count.
             */
            oobcnt = oobcnt+1;
        }
    }
    if( oobcnt>0 )
    {
        rep->oobrelclserror = rep->oobrelclserror/oobcnt;
        rep->oobavgce = rep->oobavgce/oobcnt;
        rep->oobrmserror = ae_sqrt(rep->oobrmserror/(oobcnt*nclasses), _state);
        rep->oobavgerror = rep->oobavgerror/(oobcnt*nclasses);
        if( oobrelcnt>0 )
        {
            rep->oobavgrelerror = rep->oobavgrelerror/oobrelcnt;
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Procesing

INPUT PARAMETERS:
    DF      -   decision forest model
    X       -   input vector,  array[0..NVars-1].

OUTPUT PARAMETERS:
    Y       -   result. Regression estimate when solving regression  task,
                vector of posterior probabilities for classification task.

See also DFProcessI.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
void dfprocess(decisionforest* df,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t offs;
    ae_int_t i;
    double v;


    
    /*
     * Proceed
     */
    if( y->cnt<df->nclasses )
    {
        ae_vector_set_length(y, df->nclasses, _state);
    }
    offs = 0;
    for(i=0; i<=df->nclasses-1; i++)
    {
        y->ptr.p_double[i] = 0;
    }
    for(i=0; i<=df->ntrees-1; i++)
    {
        
        /*
         * Process basic tree
         */
        dforest_dfprocessinternal(df, offs, x, y, _state);
        
        /*
         * Next tree
         */
        offs = offs+ae_round(df->trees.ptr.p_double[offs], _state);
    }
    v = (double)1/(double)df->ntrees;
    ae_v_muld(&y->ptr.p_double[0], 1, ae_v_len(0,df->nclasses-1), v);
}


/*************************************************************************
'interactive' variant of DFProcess for languages like Python which support
constructs like "Y = DFProcessI(DF,X)" and interactive mode of interpreter

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void dfprocessi(decisionforest* df,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{

    ae_vector_clear(y);

    dfprocess(df, x, y, _state);
}


/*************************************************************************
Relative classification error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    percent of incorrectly classified cases.
    Zero if model solves regression task.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfrelclserror(decisionforest* df,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    double result;


    result = (double)dforest_dfclserror(df, xy, npoints, _state)/(double)npoints;
    return result;
}


/*************************************************************************
Average cross-entropy (in bits per element) on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    CrossEntropy/(NPoints*LN(2)).
    Zero if model solves regression task.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfavgce(decisionforest* df,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector x;
    ae_vector y;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t tmpi;
    double result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);

    ae_vector_set_length(&x, df->nvars-1+1, _state);
    ae_vector_set_length(&y, df->nclasses-1+1, _state);
    result = 0;
    for(i=0; i<=npoints-1; i++)
    {
        ae_v_move(&x.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,df->nvars-1));
        dfprocess(df, &x, &y, _state);
        if( df->nclasses>1 )
        {
            
            /*
             * classification-specific code
             */
            k = ae_round(xy->ptr.pp_double[i][df->nvars], _state);
            tmpi = 0;
            for(j=1; j<=df->nclasses-1; j++)
            {
                if( ae_fp_greater(y.ptr.p_double[j],y.ptr.p_double[tmpi]) )
                {
                    tmpi = j;
                }
            }
            if( ae_fp_neq(y.ptr.p_double[k],0) )
            {
                result = result-ae_log(y.ptr.p_double[k], _state);
            }
            else
            {
                result = result-ae_log(ae_minrealnumber, _state);
            }
        }
    }
    result = result/npoints;
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.
    Its meaning for regression task is obvious. As for
    classification task, RMS error means error when estimating posterior
    probabilities.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfrmserror(decisionforest* df,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector x;
    ae_vector y;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t tmpi;
    double result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);

    ae_vector_set_length(&x, df->nvars-1+1, _state);
    ae_vector_set_length(&y, df->nclasses-1+1, _state);
    result = 0;
    for(i=0; i<=npoints-1; i++)
    {
        ae_v_move(&x.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,df->nvars-1));
        dfprocess(df, &x, &y, _state);
        if( df->nclasses>1 )
        {
            
            /*
             * classification-specific code
             */
            k = ae_round(xy->ptr.pp_double[i][df->nvars], _state);
            tmpi = 0;
            for(j=1; j<=df->nclasses-1; j++)
            {
                if( ae_fp_greater(y.ptr.p_double[j],y.ptr.p_double[tmpi]) )
                {
                    tmpi = j;
                }
            }
            for(j=0; j<=df->nclasses-1; j++)
            {
                if( j==k )
                {
                    result = result+ae_sqr(y.ptr.p_double[j]-1, _state);
                }
                else
                {
                    result = result+ae_sqr(y.ptr.p_double[j], _state);
                }
            }
        }
        else
        {
            
            /*
             * regression-specific code
             */
            result = result+ae_sqr(y.ptr.p_double[0]-xy->ptr.pp_double[i][df->nvars], _state);
        }
    }
    result = ae_sqrt(result/(npoints*df->nclasses), _state);
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for
    classification task, it means average error when estimating posterior
    probabilities.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfavgerror(decisionforest* df,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector x;
    ae_vector y;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    double result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);

    ae_vector_set_length(&x, df->nvars-1+1, _state);
    ae_vector_set_length(&y, df->nclasses-1+1, _state);
    result = 0;
    for(i=0; i<=npoints-1; i++)
    {
        ae_v_move(&x.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,df->nvars-1));
        dfprocess(df, &x, &y, _state);
        if( df->nclasses>1 )
        {
            
            /*
             * classification-specific code
             */
            k = ae_round(xy->ptr.pp_double[i][df->nvars], _state);
            for(j=0; j<=df->nclasses-1; j++)
            {
                if( j==k )
                {
                    result = result+ae_fabs(y.ptr.p_double[j]-1, _state);
                }
                else
                {
                    result = result+ae_fabs(y.ptr.p_double[j], _state);
                }
            }
        }
        else
        {
            
            /*
             * regression-specific code
             */
            result = result+ae_fabs(y.ptr.p_double[0]-xy->ptr.pp_double[i][df->nvars], _state);
        }
    }
    result = result/(npoints*df->nclasses);
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
Average relative error on the test set

INPUT PARAMETERS:
    DF      -   decision forest model
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for
    classification task, it means average relative error when estimating
    posterior probability of belonging to the correct class.

  -- ALGLIB --
     Copyright 16.02.2009 by Bochkanov Sergey
*************************************************************************/
double dfavgrelerror(decisionforest* df,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector x;
    ae_vector y;
    ae_int_t relcnt;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    double result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);

    ae_vector_set_length(&x, df->nvars-1+1, _state);
    ae_vector_set_length(&y, df->nclasses-1+1, _state);
    result = 0;
    relcnt = 0;
    for(i=0; i<=npoints-1; i++)
    {
        ae_v_move(&x.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,df->nvars-1));
        dfprocess(df, &x, &y, _state);
        if( df->nclasses>1 )
        {
            
            /*
             * classification-specific code
             */
            k = ae_round(xy->ptr.pp_double[i][df->nvars], _state);
            for(j=0; j<=df->nclasses-1; j++)
            {
                if( j==k )
                {
                    result = result+ae_fabs(y.ptr.p_double[j]-1, _state);
                    relcnt = relcnt+1;
                }
            }
        }
        else
        {
            
            /*
             * regression-specific code
             */
            if( ae_fp_neq(xy->ptr.pp_double[i][df->nvars],0) )
            {
                result = result+ae_fabs((y.ptr.p_double[0]-xy->ptr.pp_double[i][df->nvars])/xy->ptr.pp_double[i][df->nvars], _state);
                relcnt = relcnt+1;
            }
        }
    }
    if( relcnt>0 )
    {
        result = result/relcnt;
    }
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
Copying of DecisionForest strucure

INPUT PARAMETERS:
    DF1 -   original

OUTPUT PARAMETERS:
    DF2 -   copy

  -- ALGLIB --
     Copyright 13.02.2009 by Bochkanov Sergey
*************************************************************************/
void dfcopy(decisionforest* df1, decisionforest* df2, ae_state *_state)
{

    _decisionforest_clear(df2);

    df2->nvars = df1->nvars;
    df2->nclasses = df1->nclasses;
    df2->ntrees = df1->ntrees;
    df2->bufsize = df1->bufsize;
    ae_vector_set_length(&df2->trees, df1->bufsize-1+1, _state);
    ae_v_move(&df2->trees.ptr.p_double[0], 1, &df1->trees.ptr.p_double[0], 1, ae_v_len(0,df1->bufsize-1));
}


/*************************************************************************
Serializer: allocation

  -- ALGLIB --
     Copyright 14.03.2011 by Bochkanov Sergey
*************************************************************************/
void dfalloc(ae_serializer* s, decisionforest* forest, ae_state *_state)
{


    ae_serializer_alloc_entry(s);
    ae_serializer_alloc_entry(s);
    ae_serializer_alloc_entry(s);
    ae_serializer_alloc_entry(s);
    ae_serializer_alloc_entry(s);
    ae_serializer_alloc_entry(s);
    allocrealarray(s, &forest->trees, forest->bufsize, _state);
}


/*************************************************************************
Serializer: serialization

  -- ALGLIB --
     Copyright 14.03.2011 by Bochkanov Sergey
*************************************************************************/
void dfserialize(ae_serializer* s,
     decisionforest* forest,
     ae_state *_state)
{


    ae_serializer_serialize_int(s, getrdfserializationcode(_state), _state);
    ae_serializer_serialize_int(s, dforest_dffirstversion, _state);
    ae_serializer_serialize_int(s, forest->nvars, _state);
    ae_serializer_serialize_int(s, forest->nclasses, _state);
    ae_serializer_serialize_int(s, forest->ntrees, _state);
    ae_serializer_serialize_int(s, forest->bufsize, _state);
    serializerealarray(s, &forest->trees, forest->bufsize, _state);
}


/*************************************************************************
Serializer: unserialization

  -- ALGLIB --
     Copyright 14.03.2011 by Bochkanov Sergey
*************************************************************************/
void dfunserialize(ae_serializer* s,
     decisionforest* forest,
     ae_state *_state)
{
    ae_int_t i0;
    ae_int_t i1;

    _decisionforest_clear(forest);

    
    /*
     * check correctness of header
     */
    ae_serializer_unserialize_int(s, &i0, _state);
    ae_assert(i0==getrdfserializationcode(_state), "DFUnserialize: stream header corrupted", _state);
    ae_serializer_unserialize_int(s, &i1, _state);
    ae_assert(i1==dforest_dffirstversion, "DFUnserialize: stream header corrupted", _state);
    
    /*
     * Unserialize data
     */
    ae_serializer_unserialize_int(s, &forest->nvars, _state);
    ae_serializer_unserialize_int(s, &forest->nclasses, _state);
    ae_serializer_unserialize_int(s, &forest->ntrees, _state);
    ae_serializer_unserialize_int(s, &forest->bufsize, _state);
    unserializerealarray(s, &forest->trees, _state);
}


/*************************************************************************
Classification error
*************************************************************************/
static ae_int_t dforest_dfclserror(decisionforest* df,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector x;
    ae_vector y;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t tmpi;
    ae_int_t result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);

    if( df->nclasses<=1 )
    {
        result = 0;
        ae_frame_leave(_state);
        return result;
    }
    ae_vector_set_length(&x, df->nvars-1+1, _state);
    ae_vector_set_length(&y, df->nclasses-1+1, _state);
    result = 0;
    for(i=0; i<=npoints-1; i++)
    {
        ae_v_move(&x.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,df->nvars-1));
        dfprocess(df, &x, &y, _state);
        k = ae_round(xy->ptr.pp_double[i][df->nvars], _state);
        tmpi = 0;
        for(j=1; j<=df->nclasses-1; j++)
        {
            if( ae_fp_greater(y.ptr.p_double[j],y.ptr.p_double[tmpi]) )
            {
                tmpi = j;
            }
        }
        if( tmpi!=k )
        {
            result = result+1;
        }
    }
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
Internal subroutine for processing one decision tree starting at Offs
*************************************************************************/
static void dforest_dfprocessinternal(decisionforest* df,
     ae_int_t offs,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t k;
    ae_int_t idx;


    
    /*
     * Set pointer to the root
     */
    k = offs+1;
    
    /*
     * Navigate through the tree
     */
    for(;;)
    {
        if( ae_fp_eq(df->trees.ptr.p_double[k],-1) )
        {
            if( df->nclasses==1 )
            {
                y->ptr.p_double[0] = y->ptr.p_double[0]+df->trees.ptr.p_double[k+1];
            }
            else
            {
                idx = ae_round(df->trees.ptr.p_double[k+1], _state);
                y->ptr.p_double[idx] = y->ptr.p_double[idx]+1;
            }
            break;
        }
        if( ae_fp_less(x->ptr.p_double[ae_round(df->trees.ptr.p_double[k], _state)],df->trees.ptr.p_double[k+1]) )
        {
            k = k+dforest_innernodewidth;
        }
        else
        {
            k = offs+ae_round(df->trees.ptr.p_double[k+2], _state);
        }
    }
}


/*************************************************************************
Builds one decision tree. Just a wrapper for the DFBuildTreeRec.
*************************************************************************/
static void dforest_dfbuildtree(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t nfeatures,
     ae_int_t nvarsinpool,
     ae_int_t flags,
     dfinternalbuffers* bufs,
     ae_state *_state)
{
    ae_int_t numprocessed;
    ae_int_t i;


    ae_assert(npoints>0, "Assertion failed", _state);
    
    /*
     * Prepare IdxBuf. It stores indices of the training set elements.
     * When training set is being split, contents of IdxBuf is
     * correspondingly reordered so we can know which elements belong
     * to which branch of decision tree.
     */
    for(i=0; i<=npoints-1; i++)
    {
        bufs->idxbuf.ptr.p_int[i] = i;
    }
    
    /*
     * Recursive procedure
     */
    numprocessed = 1;
    dforest_dfbuildtreerec(xy, npoints, nvars, nclasses, nfeatures, nvarsinpool, flags, &numprocessed, 0, npoints-1, bufs, _state);
    bufs->treebuf.ptr.p_double[0] = numprocessed;
}


/*************************************************************************
Builds one decision tree (internal recursive subroutine)

Parameters:
    TreeBuf     -   large enough array, at least TreeSize
    IdxBuf      -   at least NPoints elements
    TmpBufR     -   at least NPoints
    TmpBufR2    -   at least NPoints
    TmpBufI     -   at least NPoints
    TmpBufI2    -   at least NPoints+1
*************************************************************************/
static void dforest_dfbuildtreerec(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t nfeatures,
     ae_int_t nvarsinpool,
     ae_int_t flags,
     ae_int_t* numprocessed,
     ae_int_t idx1,
     ae_int_t idx2,
     dfinternalbuffers* bufs,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_bool bflag;
    ae_int_t i1;
    ae_int_t i2;
    ae_int_t info;
    double sl;
    double sr;
    double w;
    ae_int_t idxbest;
    double ebest;
    double tbest;
    ae_int_t varcur;
    double s;
    double v;
    double v1;
    double v2;
    double threshold;
    ae_int_t oldnp;
    double currms;
    ae_bool useevs;


    
    /*
     * these initializers are not really necessary,
     * but without them compiler complains about uninitialized locals
     */
    tbest = 0;
    
    /*
     * Prepare
     */
    ae_assert(npoints>0, "Assertion failed", _state);
    ae_assert(idx2>=idx1, "Assertion failed", _state);
    useevs = flags/dforest_dfuseevs%2!=0;
    
    /*
     * Leaf node
     */
    if( idx2==idx1 )
    {
        bufs->treebuf.ptr.p_double[*numprocessed] = -1;
        bufs->treebuf.ptr.p_double[*numprocessed+1] = xy->ptr.pp_double[bufs->idxbuf.ptr.p_int[idx1]][nvars];
        *numprocessed = *numprocessed+dforest_leafnodewidth;
        return;
    }
    
    /*
     * Non-leaf node.
     * Select random variable, prepare split:
     * 1. prepare default solution - no splitting, class at random
     * 2. investigate possible splits, compare with default/best
     */
    idxbest = -1;
    if( nclasses>1 )
    {
        
        /*
         * default solution for classification
         */
        for(i=0; i<=nclasses-1; i++)
        {
            bufs->classibuf.ptr.p_int[i] = 0;
        }
        s = idx2-idx1+1;
        for(i=idx1; i<=idx2; i++)
        {
            j = ae_round(xy->ptr.pp_double[bufs->idxbuf.ptr.p_int[i]][nvars], _state);
            bufs->classibuf.ptr.p_int[j] = bufs->classibuf.ptr.p_int[j]+1;
        }
        ebest = 0;
        for(i=0; i<=nclasses-1; i++)
        {
            ebest = ebest+bufs->classibuf.ptr.p_int[i]*ae_sqr(1-bufs->classibuf.ptr.p_int[i]/s, _state)+(s-bufs->classibuf.ptr.p_int[i])*ae_sqr(bufs->classibuf.ptr.p_int[i]/s, _state);
        }
        ebest = ae_sqrt(ebest/(nclasses*(idx2-idx1+1)), _state);
    }
    else
    {
        
        /*
         * default solution for regression
         */
        v = 0;
        for(i=idx1; i<=idx2; i++)
        {
            v = v+xy->ptr.pp_double[bufs->idxbuf.ptr.p_int[i]][nvars];
        }
        v = v/(idx2-idx1+1);
        ebest = 0;
        for(i=idx1; i<=idx2; i++)
        {
            ebest = ebest+ae_sqr(xy->ptr.pp_double[bufs->idxbuf.ptr.p_int[i]][nvars]-v, _state);
        }
        ebest = ae_sqrt(ebest/(idx2-idx1+1), _state);
    }
    i = 0;
    while(i<=ae_minint(nfeatures, nvarsinpool, _state)-1)
    {
        
        /*
         * select variables from pool
         */
        j = i+ae_randominteger(nvarsinpool-i, _state);
        k = bufs->varpool.ptr.p_int[i];
        bufs->varpool.ptr.p_int[i] = bufs->varpool.ptr.p_int[j];
        bufs->varpool.ptr.p_int[j] = k;
        varcur = bufs->varpool.ptr.p_int[i];
        
        /*
         * load variable values to working array
         *
         * apply EVS preprocessing: if all variable values are same,
         * variable is excluded from pool.
         *
         * This is necessary for binary pre-splits (see later) to work.
         */
        for(j=idx1; j<=idx2; j++)
        {
            bufs->tmpbufr.ptr.p_double[j-idx1] = xy->ptr.pp_double[bufs->idxbuf.ptr.p_int[j]][varcur];
        }
        if( useevs )
        {
            bflag = ae_false;
            v = bufs->tmpbufr.ptr.p_double[0];
            for(j=0; j<=idx2-idx1; j++)
            {
                if( ae_fp_neq(bufs->tmpbufr.ptr.p_double[j],v) )
                {
                    bflag = ae_true;
                    break;
                }
            }
            if( !bflag )
            {
                
                /*
                 * exclude variable from pool,
                 * go to the next iteration.
                 * I is not increased.
                 */
                k = bufs->varpool.ptr.p_int[i];
                bufs->varpool.ptr.p_int[i] = bufs->varpool.ptr.p_int[nvarsinpool-1];
                bufs->varpool.ptr.p_int[nvarsinpool-1] = k;
                nvarsinpool = nvarsinpool-1;
                continue;
            }
        }
        
        /*
         * load labels to working array
         */
        if( nclasses>1 )
        {
            for(j=idx1; j<=idx2; j++)
            {
                bufs->tmpbufi.ptr.p_int[j-idx1] = ae_round(xy->ptr.pp_double[bufs->idxbuf.ptr.p_int[j]][nvars], _state);
            }
        }
        else
        {
            for(j=idx1; j<=idx2; j++)
            {
                bufs->tmpbufr2.ptr.p_double[j-idx1] = xy->ptr.pp_double[bufs->idxbuf.ptr.p_int[j]][nvars];
            }
        }
        
        /*
         * calculate split
         */
        if( useevs&&bufs->evsbin.ptr.p_bool[varcur] )
        {
            
            /*
             * Pre-calculated splits for binary variables.
             * Threshold is already known, just calculate RMS error
             */
            threshold = bufs->evssplits.ptr.p_double[varcur];
            if( nclasses>1 )
            {
                
                /*
                 * classification-specific code
                 */
                for(j=0; j<=2*nclasses-1; j++)
                {
                    bufs->classibuf.ptr.p_int[j] = 0;
                }
                sl = 0;
                sr = 0;
                for(j=0; j<=idx2-idx1; j++)
                {
                    k = bufs->tmpbufi.ptr.p_int[j];
                    if( ae_fp_less(bufs->tmpbufr.ptr.p_double[j],threshold) )
                    {
                        bufs->classibuf.ptr.p_int[k] = bufs->classibuf.ptr.p_int[k]+1;
                        sl = sl+1;
                    }
                    else
                    {
                        bufs->classibuf.ptr.p_int[k+nclasses] = bufs->classibuf.ptr.p_int[k+nclasses]+1;
                        sr = sr+1;
                    }
                }
                ae_assert(ae_fp_neq(sl,0)&&ae_fp_neq(sr,0), "DFBuildTreeRec: something strange!", _state);
                currms = 0;
                for(j=0; j<=nclasses-1; j++)
                {
                    w = bufs->classibuf.ptr.p_int[j];
                    currms = currms+w*ae_sqr(w/sl-1, _state);
                    currms = currms+(sl-w)*ae_sqr(w/sl, _state);
                    w = bufs->classibuf.ptr.p_int[nclasses+j];
                    currms = currms+w*ae_sqr(w/sr-1, _state);
                    currms = currms+(sr-w)*ae_sqr(w/sr, _state);
                }
                currms = ae_sqrt(currms/(nclasses*(idx2-idx1+1)), _state);
            }
            else
            {
                
                /*
                 * regression-specific code
                 */
                sl = 0;
                sr = 0;
                v1 = 0;
                v2 = 0;
                for(j=0; j<=idx2-idx1; j++)
                {
                    if( ae_fp_less(bufs->tmpbufr.ptr.p_double[j],threshold) )
                    {
                        v1 = v1+bufs->tmpbufr2.ptr.p_double[j];
                        sl = sl+1;
                    }
                    else
                    {
                        v2 = v2+bufs->tmpbufr2.ptr.p_double[j];
                        sr = sr+1;
                    }
                }
                ae_assert(ae_fp_neq(sl,0)&&ae_fp_neq(sr,0), "DFBuildTreeRec: something strange!", _state);
                v1 = v1/sl;
                v2 = v2/sr;
                currms = 0;
                for(j=0; j<=idx2-idx1; j++)
                {
                    if( ae_fp_less(bufs->tmpbufr.ptr.p_double[j],threshold) )
                    {
                        currms = currms+ae_sqr(v1-bufs->tmpbufr2.ptr.p_double[j], _state);
                    }
                    else
                    {
                        currms = currms+ae_sqr(v2-bufs->tmpbufr2.ptr.p_double[j], _state);
                    }
                }
                currms = ae_sqrt(currms/(idx2-idx1+1), _state);
            }
            info = 1;
        }
        else
        {
            
            /*
             * Generic splits
             */
            if( nclasses>1 )
            {
                dforest_dfsplitc(&bufs->tmpbufr, &bufs->tmpbufi, &bufs->classibuf, idx2-idx1+1, nclasses, dforest_dfusestrongsplits, &info, &threshold, &currms, &bufs->sortrbuf, &bufs->sortibuf, _state);
            }
            else
            {
                dforest_dfsplitr(&bufs->tmpbufr, &bufs->tmpbufr2, idx2-idx1+1, dforest_dfusestrongsplits, &info, &threshold, &currms, &bufs->sortrbuf, &bufs->sortrbuf2, _state);
            }
        }
        if( info>0 )
        {
            if( ae_fp_less_eq(currms,ebest) )
            {
                ebest = currms;
                idxbest = varcur;
                tbest = threshold;
            }
        }
        
        /*
         * Next iteration
         */
        i = i+1;
    }
    
    /*
     * to split or not to split
     */
    if( idxbest<0 )
    {
        
        /*
         * All values are same, cannot split.
         */
        bufs->treebuf.ptr.p_double[*numprocessed] = -1;
        if( nclasses>1 )
        {
            
            /*
             * Select random class label (randomness allows us to
             * approximate distribution of the classes)
             */
            bufs->treebuf.ptr.p_double[*numprocessed+1] = ae_round(xy->ptr.pp_double[bufs->idxbuf.ptr.p_int[idx1+ae_randominteger(idx2-idx1+1, _state)]][nvars], _state);
        }
        else
        {
            
            /*
             * Select average (for regression task).
             */
            v = 0;
            for(i=idx1; i<=idx2; i++)
            {
                v = v+xy->ptr.pp_double[bufs->idxbuf.ptr.p_int[i]][nvars]/(idx2-idx1+1);
            }
            bufs->treebuf.ptr.p_double[*numprocessed+1] = v;
        }
        *numprocessed = *numprocessed+dforest_leafnodewidth;
    }
    else
    {
        
        /*
         * we can split
         */
        bufs->treebuf.ptr.p_double[*numprocessed] = idxbest;
        bufs->treebuf.ptr.p_double[*numprocessed+1] = tbest;
        i1 = idx1;
        i2 = idx2;
        while(i1<=i2)
        {
            
            /*
             * Reorder indices so that left partition is in [Idx1..I1-1],
             * and right partition is in [I2+1..Idx2]
             */
            if( ae_fp_less(xy->ptr.pp_double[bufs->idxbuf.ptr.p_int[i1]][idxbest],tbest) )
            {
                i1 = i1+1;
                continue;
            }
            if( ae_fp_greater_eq(xy->ptr.pp_double[bufs->idxbuf.ptr.p_int[i2]][idxbest],tbest) )
            {
                i2 = i2-1;
                continue;
            }
            j = bufs->idxbuf.ptr.p_int[i1];
            bufs->idxbuf.ptr.p_int[i1] = bufs->idxbuf.ptr.p_int[i2];
            bufs->idxbuf.ptr.p_int[i2] = j;
            i1 = i1+1;
            i2 = i2-1;
        }
        oldnp = *numprocessed;
        *numprocessed = *numprocessed+dforest_innernodewidth;
        dforest_dfbuildtreerec(xy, npoints, nvars, nclasses, nfeatures, nvarsinpool, flags, numprocessed, idx1, i1-1, bufs, _state);
        bufs->treebuf.ptr.p_double[oldnp+2] = *numprocessed;
        dforest_dfbuildtreerec(xy, npoints, nvars, nclasses, nfeatures, nvarsinpool, flags, numprocessed, i2+1, idx2, bufs, _state);
    }
}


/*************************************************************************
Makes split on attribute
*************************************************************************/
static void dforest_dfsplitc(/* Real    */ ae_vector* x,
     /* Integer */ ae_vector* c,
     /* Integer */ ae_vector* cntbuf,
     ae_int_t n,
     ae_int_t nc,
     ae_int_t flags,
     ae_int_t* info,
     double* threshold,
     double* e,
     /* Real    */ ae_vector* sortrbuf,
     /* Integer */ ae_vector* sortibuf,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t neq;
    ae_int_t nless;
    ae_int_t ngreater;
    ae_int_t q;
    ae_int_t qmin;
    ae_int_t qmax;
    ae_int_t qcnt;
    double cursplit;
    ae_int_t nleft;
    double v;
    double cure;
    double w;
    double sl;
    double sr;

    *info = 0;
    *threshold = 0;
    *e = 0;

    tagsortfasti(x, c, sortrbuf, sortibuf, n, _state);
    *e = ae_maxrealnumber;
    *threshold = 0.5*(x->ptr.p_double[0]+x->ptr.p_double[n-1]);
    *info = -3;
    if( flags/dforest_dfusestrongsplits%2==0 )
    {
        
        /*
         * weak splits, split at half
         */
        qcnt = 2;
        qmin = 1;
        qmax = 1;
    }
    else
    {
        
        /*
         * strong splits: choose best quartile
         */
        qcnt = 4;
        qmin = 1;
        qmax = 3;
    }
    for(q=qmin; q<=qmax; q++)
    {
        cursplit = x->ptr.p_double[n*q/qcnt];
        neq = 0;
        nless = 0;
        ngreater = 0;
        for(i=0; i<=n-1; i++)
        {
            if( ae_fp_less(x->ptr.p_double[i],cursplit) )
            {
                nless = nless+1;
            }
            if( ae_fp_eq(x->ptr.p_double[i],cursplit) )
            {
                neq = neq+1;
            }
            if( ae_fp_greater(x->ptr.p_double[i],cursplit) )
            {
                ngreater = ngreater+1;
            }
        }
        ae_assert(neq!=0, "DFSplitR: NEq=0, something strange!!!", _state);
        if( nless!=0||ngreater!=0 )
        {
            
            /*
             * set threshold between two partitions, with
             * some tweaking to avoid problems with floating point
             * arithmetics.
             *
             * The problem is that when you calculates C = 0.5*(A+B) there
             * can be no C which lies strictly between A and B (for example,
             * there is no floating point number which is
             * greater than 1 and less than 1+eps). In such situations
             * we choose right side as theshold (remember that
             * points which lie on threshold falls to the right side).
             */
            if( nless<ngreater )
            {
                cursplit = 0.5*(x->ptr.p_double[nless+neq-1]+x->ptr.p_double[nless+neq]);
                nleft = nless+neq;
                if( ae_fp_less_eq(cursplit,x->ptr.p_double[nless+neq-1]) )
                {
                    cursplit = x->ptr.p_double[nless+neq];
                }
            }
            else
            {
                cursplit = 0.5*(x->ptr.p_double[nless-1]+x->ptr.p_double[nless]);
                nleft = nless;
                if( ae_fp_less_eq(cursplit,x->ptr.p_double[nless-1]) )
                {
                    cursplit = x->ptr.p_double[nless];
                }
            }
            *info = 1;
            cure = 0;
            for(i=0; i<=2*nc-1; i++)
            {
                cntbuf->ptr.p_int[i] = 0;
            }
            for(i=0; i<=nleft-1; i++)
            {
                cntbuf->ptr.p_int[c->ptr.p_int[i]] = cntbuf->ptr.p_int[c->ptr.p_int[i]]+1;
            }
            for(i=nleft; i<=n-1; i++)
            {
                cntbuf->ptr.p_int[nc+c->ptr.p_int[i]] = cntbuf->ptr.p_int[nc+c->ptr.p_int[i]]+1;
            }
            sl = nleft;
            sr = n-nleft;
            v = 0;
            for(i=0; i<=nc-1; i++)
            {
                w = cntbuf->ptr.p_int[i];
                v = v+w*ae_sqr(w/sl-1, _state);
                v = v+(sl-w)*ae_sqr(w/sl, _state);
                w = cntbuf->ptr.p_int[nc+i];
                v = v+w*ae_sqr(w/sr-1, _state);
                v = v+(sr-w)*ae_sqr(w/sr, _state);
            }
            cure = ae_sqrt(v/(nc*n), _state);
            if( ae_fp_less(cure,*e) )
            {
                *threshold = cursplit;
                *e = cure;
            }
        }
    }
}


/*************************************************************************
Makes split on attribute
*************************************************************************/
static void dforest_dfsplitr(/* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_int_t n,
     ae_int_t flags,
     ae_int_t* info,
     double* threshold,
     double* e,
     /* Real    */ ae_vector* sortrbuf,
     /* Real    */ ae_vector* sortrbuf2,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t neq;
    ae_int_t nless;
    ae_int_t ngreater;
    ae_int_t q;
    ae_int_t qmin;
    ae_int_t qmax;
    ae_int_t qcnt;
    double cursplit;
    ae_int_t nleft;
    double v;
    double cure;

    *info = 0;
    *threshold = 0;
    *e = 0;

    tagsortfastr(x, y, sortrbuf, sortrbuf2, n, _state);
    *e = ae_maxrealnumber;
    *threshold = 0.5*(x->ptr.p_double[0]+x->ptr.p_double[n-1]);
    *info = -3;
    if( flags/dforest_dfusestrongsplits%2==0 )
    {
        
        /*
         * weak splits, split at half
         */
        qcnt = 2;
        qmin = 1;
        qmax = 1;
    }
    else
    {
        
        /*
         * strong splits: choose best quartile
         */
        qcnt = 4;
        qmin = 1;
        qmax = 3;
    }
    for(q=qmin; q<=qmax; q++)
    {
        cursplit = x->ptr.p_double[n*q/qcnt];
        neq = 0;
        nless = 0;
        ngreater = 0;
        for(i=0; i<=n-1; i++)
        {
            if( ae_fp_less(x->ptr.p_double[i],cursplit) )
            {
                nless = nless+1;
            }
            if( ae_fp_eq(x->ptr.p_double[i],cursplit) )
            {
                neq = neq+1;
            }
            if( ae_fp_greater(x->ptr.p_double[i],cursplit) )
            {
                ngreater = ngreater+1;
            }
        }
        ae_assert(neq!=0, "DFSplitR: NEq=0, something strange!!!", _state);
        if( nless!=0||ngreater!=0 )
        {
            
            /*
             * set threshold between two partitions, with
             * some tweaking to avoid problems with floating point
             * arithmetics.
             *
             * The problem is that when you calculates C = 0.5*(A+B) there
             * can be no C which lies strictly between A and B (for example,
             * there is no floating point number which is
             * greater than 1 and less than 1+eps). In such situations
             * we choose right side as theshold (remember that
             * points which lie on threshold falls to the right side).
             */
            if( nless<ngreater )
            {
                cursplit = 0.5*(x->ptr.p_double[nless+neq-1]+x->ptr.p_double[nless+neq]);
                nleft = nless+neq;
                if( ae_fp_less_eq(cursplit,x->ptr.p_double[nless+neq-1]) )
                {
                    cursplit = x->ptr.p_double[nless+neq];
                }
            }
            else
            {
                cursplit = 0.5*(x->ptr.p_double[nless-1]+x->ptr.p_double[nless]);
                nleft = nless;
                if( ae_fp_less_eq(cursplit,x->ptr.p_double[nless-1]) )
                {
                    cursplit = x->ptr.p_double[nless];
                }
            }
            *info = 1;
            cure = 0;
            v = 0;
            for(i=0; i<=nleft-1; i++)
            {
                v = v+y->ptr.p_double[i];
            }
            v = v/nleft;
            for(i=0; i<=nleft-1; i++)
            {
                cure = cure+ae_sqr(y->ptr.p_double[i]-v, _state);
            }
            v = 0;
            for(i=nleft; i<=n-1; i++)
            {
                v = v+y->ptr.p_double[i];
            }
            v = v/(n-nleft);
            for(i=nleft; i<=n-1; i++)
            {
                cure = cure+ae_sqr(y->ptr.p_double[i]-v, _state);
            }
            cure = ae_sqrt(cure/n, _state);
            if( ae_fp_less(cure,*e) )
            {
                *threshold = cursplit;
                *e = cure;
            }
        }
    }
}


ae_bool _decisionforest_init(decisionforest* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->trees, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _decisionforest_init_copy(decisionforest* dst, decisionforest* src, ae_state *_state, ae_bool make_automatic)
{
    dst->nvars = src->nvars;
    dst->nclasses = src->nclasses;
    dst->ntrees = src->ntrees;
    dst->bufsize = src->bufsize;
    if( !ae_vector_init_copy(&dst->trees, &src->trees, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _decisionforest_clear(decisionforest* p)
{
    ae_vector_clear(&p->trees);
}


ae_bool _dfreport_init(dfreport* /*p*/, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    return ae_true;
}


ae_bool _dfreport_init_copy(dfreport* dst, dfreport* src, ae_state */*_state*/, ae_bool) // make_automatic)
{
    dst->relclserror = src->relclserror;
    dst->avgce = src->avgce;
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
    dst->oobrelclserror = src->oobrelclserror;
    dst->oobavgce = src->oobavgce;
    dst->oobrmserror = src->oobrmserror;
    dst->oobavgerror = src->oobavgerror;
    dst->oobavgrelerror = src->oobavgrelerror;
    return ae_true;
}


void _dfreport_clear(dfreport* ) //p)
{
}


ae_bool _dfinternalbuffers_init(dfinternalbuffers* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->treebuf, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->idxbuf, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tmpbufr, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tmpbufr2, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tmpbufi, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->classibuf, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->sortrbuf, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->sortrbuf2, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->sortibuf, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->varpool, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->evsbin, 0, DT_BOOL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->evssplits, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _dfinternalbuffers_init_copy(dfinternalbuffers* dst, dfinternalbuffers* src, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init_copy(&dst->treebuf, &src->treebuf, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->idxbuf, &src->idxbuf, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tmpbufr, &src->tmpbufr, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tmpbufr2, &src->tmpbufr2, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tmpbufi, &src->tmpbufi, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->classibuf, &src->classibuf, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->sortrbuf, &src->sortrbuf, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->sortrbuf2, &src->sortrbuf2, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->sortibuf, &src->sortibuf, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->varpool, &src->varpool, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->evsbin, &src->evsbin, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->evssplits, &src->evssplits, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _dfinternalbuffers_clear(dfinternalbuffers* p)
{
    ae_vector_clear(&p->treebuf);
    ae_vector_clear(&p->idxbuf);
    ae_vector_clear(&p->tmpbufr);
    ae_vector_clear(&p->tmpbufr2);
    ae_vector_clear(&p->tmpbufi);
    ae_vector_clear(&p->classibuf);
    ae_vector_clear(&p->sortrbuf);
    ae_vector_clear(&p->sortrbuf2);
    ae_vector_clear(&p->sortibuf);
    ae_vector_clear(&p->varpool);
    ae_vector_clear(&p->evsbin);
    ae_vector_clear(&p->evssplits);
}




/*************************************************************************
Linear regression

Subroutine builds model:

    Y = A(0)*X[0] + ... + A(N-1)*X[N-1] + A(N)

and model found in ALGLIB format, covariation matrix, training set  errors
(rms,  average,  average  relative)   and  leave-one-out  cross-validation
estimate of the generalization error. CV  estimate calculated  using  fast
algorithm with O(NPoints*NVars) complexity.

When  covariation  matrix  is  calculated  standard deviations of function
values are assumed to be equal to RMS error on the training set.

INPUT PARAMETERS:
    XY          -   training set, array [0..NPoints-1,0..NVars]:
                    * NVars columns - independent variables
                    * last column - dependent variable
    NPoints     -   training set size, NPoints>NVars+1
    NVars       -   number of independent variables

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -255, in case of unknown internal error
                    * -4, if internal SVD subroutine haven't converged
                    * -1, if incorrect parameters was passed (NPoints<NVars+2, NVars<1).
                    *  1, if subroutine successfully finished
    LM          -   linear model in the ALGLIB format. Use subroutines of
                    this unit to work with the model.
    AR          -   additional results


  -- ALGLIB --
     Copyright 02.08.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuild(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     linearmodel* lm,
     lrreport* ar,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector s;
    ae_int_t i;
    double sigma2;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _linearmodel_clear(lm);
    _lrreport_clear(ar);
    ae_vector_init(&s, 0, DT_REAL, _state, ae_true);

    if( npoints<=nvars+1||nvars<1 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    ae_vector_set_length(&s, npoints-1+1, _state);
    for(i=0; i<=npoints-1; i++)
    {
        s.ptr.p_double[i] = 1;
    }
    lrbuilds(xy, &s, npoints, nvars, info, lm, ar, _state);
    if( *info<0 )
    {
        ae_frame_leave(_state);
        return;
    }
    sigma2 = ae_sqr(ar->rmserror, _state)*npoints/(npoints-nvars-1);
    for(i=0; i<=nvars; i++)
    {
        ae_v_muld(&ar->c.ptr.pp_double[i][0], 1, ae_v_len(0,nvars), sigma2);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Linear regression

Variant of LRBuild which uses vector of standatd deviations (errors in
function values).

INPUT PARAMETERS:
    XY          -   training set, array [0..NPoints-1,0..NVars]:
                    * NVars columns - independent variables
                    * last column - dependent variable
    S           -   standard deviations (errors in function values)
                    array[0..NPoints-1], S[i]>0.
    NPoints     -   training set size, NPoints>NVars+1
    NVars       -   number of independent variables

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -255, in case of unknown internal error
                    * -4, if internal SVD subroutine haven't converged
                    * -1, if incorrect parameters was passed (NPoints<NVars+2, NVars<1).
                    * -2, if S[I]<=0
                    *  1, if subroutine successfully finished
    LM          -   linear model in the ALGLIB format. Use subroutines of
                    this unit to work with the model.
    AR          -   additional results


  -- ALGLIB --
     Copyright 02.08.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuilds(/* Real    */ ae_matrix* xy,
     /* Real    */ ae_vector* s,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     linearmodel* lm,
     lrreport* ar,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix xyi;
    ae_vector x;
    ae_vector means;
    ae_vector sigmas;
    ae_int_t i;
    ae_int_t j;
    double v;
    ae_int_t offs;
    double mean;
    double variance;
    double skewness;
    double kurtosis;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _linearmodel_clear(lm);
    _lrreport_clear(ar);
    ae_matrix_init(&xyi, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&means, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&sigmas, 0, DT_REAL, _state, ae_true);

    
    /*
     * Test parameters
     */
    if( npoints<=nvars+1||nvars<1 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Copy data, add one more column (constant term)
     */
    ae_matrix_set_length(&xyi, npoints-1+1, nvars+1+1, _state);
    for(i=0; i<=npoints-1; i++)
    {
        ae_v_move(&xyi.ptr.pp_double[i][0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
        xyi.ptr.pp_double[i][nvars] = 1;
        xyi.ptr.pp_double[i][nvars+1] = xy->ptr.pp_double[i][nvars];
    }
    
    /*
     * Standartization
     */
    ae_vector_set_length(&x, npoints-1+1, _state);
    ae_vector_set_length(&means, nvars-1+1, _state);
    ae_vector_set_length(&sigmas, nvars-1+1, _state);
    for(j=0; j<=nvars-1; j++)
    {
        ae_v_move(&x.ptr.p_double[0], 1, &xy->ptr.pp_double[0][j], xy->stride, ae_v_len(0,npoints-1));
        samplemoments(&x, npoints, &mean, &variance, &skewness, &kurtosis, _state);
        means.ptr.p_double[j] = mean;
        sigmas.ptr.p_double[j] = ae_sqrt(variance, _state);
        if( ae_fp_eq(sigmas.ptr.p_double[j],0) )
        {
            sigmas.ptr.p_double[j] = 1;
        }
        for(i=0; i<=npoints-1; i++)
        {
            xyi.ptr.pp_double[i][j] = (xyi.ptr.pp_double[i][j]-means.ptr.p_double[j])/sigmas.ptr.p_double[j];
        }
    }
    
    /*
     * Internal processing
     */
    linreg_lrinternal(&xyi, s, npoints, nvars+1, info, lm, ar, _state);
    if( *info<0 )
    {
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Un-standartization
     */
    offs = ae_round(lm->w.ptr.p_double[3], _state);
    for(j=0; j<=nvars-1; j++)
    {
        
        /*
         * Constant term is updated (and its covariance too,
         * since it gets some variance from J-th component)
         */
        lm->w.ptr.p_double[offs+nvars] = lm->w.ptr.p_double[offs+nvars]-lm->w.ptr.p_double[offs+j]*means.ptr.p_double[j]/sigmas.ptr.p_double[j];
        v = means.ptr.p_double[j]/sigmas.ptr.p_double[j];
        ae_v_subd(&ar->c.ptr.pp_double[nvars][0], 1, &ar->c.ptr.pp_double[j][0], 1, ae_v_len(0,nvars), v);
        ae_v_subd(&ar->c.ptr.pp_double[0][nvars], ar->c.stride, &ar->c.ptr.pp_double[0][j], ar->c.stride, ae_v_len(0,nvars), v);
        
        /*
         * J-th term is updated
         */
        lm->w.ptr.p_double[offs+j] = lm->w.ptr.p_double[offs+j]/sigmas.ptr.p_double[j];
        v = 1/sigmas.ptr.p_double[j];
        ae_v_muld(&ar->c.ptr.pp_double[j][0], 1, ae_v_len(0,nvars), v);
        ae_v_muld(&ar->c.ptr.pp_double[0][j], ar->c.stride, ae_v_len(0,nvars), v);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Like LRBuildS, but builds model

    Y = A(0)*X[0] + ... + A(N-1)*X[N-1]

i.e. with zero constant term.

  -- ALGLIB --
     Copyright 30.10.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuildzs(/* Real    */ ae_matrix* xy,
     /* Real    */ ae_vector* s,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     linearmodel* lm,
     lrreport* ar,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix xyi;
    ae_vector x;
    ae_vector c;
    ae_int_t i;
    ae_int_t j;
    double v;
    ae_int_t offs;
    double mean;
    double variance;
    double skewness;
    double kurtosis;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _linearmodel_clear(lm);
    _lrreport_clear(ar);
    ae_matrix_init(&xyi, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&c, 0, DT_REAL, _state, ae_true);

    
    /*
     * Test parameters
     */
    if( npoints<=nvars+1||nvars<1 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Copy data, add one more column (constant term)
     */
    ae_matrix_set_length(&xyi, npoints-1+1, nvars+1+1, _state);
    for(i=0; i<=npoints-1; i++)
    {
        ae_v_move(&xyi.ptr.pp_double[i][0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
        xyi.ptr.pp_double[i][nvars] = 0;
        xyi.ptr.pp_double[i][nvars+1] = xy->ptr.pp_double[i][nvars];
    }
    
    /*
     * Standartization: unusual scaling
     */
    ae_vector_set_length(&x, npoints-1+1, _state);
    ae_vector_set_length(&c, nvars-1+1, _state);
    for(j=0; j<=nvars-1; j++)
    {
        ae_v_move(&x.ptr.p_double[0], 1, &xy->ptr.pp_double[0][j], xy->stride, ae_v_len(0,npoints-1));
        samplemoments(&x, npoints, &mean, &variance, &skewness, &kurtosis, _state);
        if( ae_fp_greater(ae_fabs(mean, _state),ae_sqrt(variance, _state)) )
        {
            
            /*
             * variation is relatively small, it is better to
             * bring mean value to 1
             */
            c.ptr.p_double[j] = mean;
        }
        else
        {
            
            /*
             * variation is large, it is better to bring variance to 1
             */
            if( ae_fp_eq(variance,0) )
            {
                variance = 1;
            }
            c.ptr.p_double[j] = ae_sqrt(variance, _state);
        }
        for(i=0; i<=npoints-1; i++)
        {
            xyi.ptr.pp_double[i][j] = xyi.ptr.pp_double[i][j]/c.ptr.p_double[j];
        }
    }
    
    /*
     * Internal processing
     */
    linreg_lrinternal(&xyi, s, npoints, nvars+1, info, lm, ar, _state);
    if( *info<0 )
    {
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Un-standartization
     */
    offs = ae_round(lm->w.ptr.p_double[3], _state);
    for(j=0; j<=nvars-1; j++)
    {
        
        /*
         * J-th term is updated
         */
        lm->w.ptr.p_double[offs+j] = lm->w.ptr.p_double[offs+j]/c.ptr.p_double[j];
        v = 1/c.ptr.p_double[j];
        ae_v_muld(&ar->c.ptr.pp_double[j][0], 1, ae_v_len(0,nvars), v);
        ae_v_muld(&ar->c.ptr.pp_double[0][j], ar->c.stride, ae_v_len(0,nvars), v);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Like LRBuild but builds model

    Y = A(0)*X[0] + ... + A(N-1)*X[N-1]

i.e. with zero constant term.

  -- ALGLIB --
     Copyright 30.10.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuildz(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     linearmodel* lm,
     lrreport* ar,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector s;
    ae_int_t i;
    double sigma2;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _linearmodel_clear(lm);
    _lrreport_clear(ar);
    ae_vector_init(&s, 0, DT_REAL, _state, ae_true);

    if( npoints<=nvars+1||nvars<1 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    ae_vector_set_length(&s, npoints-1+1, _state);
    for(i=0; i<=npoints-1; i++)
    {
        s.ptr.p_double[i] = 1;
    }
    lrbuildzs(xy, &s, npoints, nvars, info, lm, ar, _state);
    if( *info<0 )
    {
        ae_frame_leave(_state);
        return;
    }
    sigma2 = ae_sqr(ar->rmserror, _state)*npoints/(npoints-nvars-1);
    for(i=0; i<=nvars; i++)
    {
        ae_v_muld(&ar->c.ptr.pp_double[i][0], 1, ae_v_len(0,nvars), sigma2);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Unpacks coefficients of linear model.

INPUT PARAMETERS:
    LM          -   linear model in ALGLIB format

OUTPUT PARAMETERS:
    V           -   coefficients, array[0..NVars]
                    constant term (intercept) is stored in the V[NVars].
    NVars       -   number of independent variables (one less than number
                    of coefficients)

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
void lrunpack(linearmodel* lm,
     /* Real    */ ae_vector* v,
     ae_int_t* nvars,
     ae_state *_state)
{
    ae_int_t offs;

    ae_vector_clear(v);
    *nvars = 0;

    ae_assert(ae_round(lm->w.ptr.p_double[1], _state)==linreg_lrvnum, "LINREG: Incorrect LINREG version!", _state);
    *nvars = ae_round(lm->w.ptr.p_double[2], _state);
    offs = ae_round(lm->w.ptr.p_double[3], _state);
    ae_vector_set_length(v, *nvars+1, _state);
    ae_v_move(&v->ptr.p_double[0], 1, &lm->w.ptr.p_double[offs], 1, ae_v_len(0,*nvars));
}


/*************************************************************************
"Packs" coefficients and creates linear model in ALGLIB format (LRUnpack
reversed).

INPUT PARAMETERS:
    V           -   coefficients, array[0..NVars]
    NVars       -   number of independent variables

OUTPUT PAREMETERS:
    LM          -   linear model.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
void lrpack(/* Real    */ ae_vector* v,
     ae_int_t nvars,
     linearmodel* lm,
     ae_state *_state)
{
    ae_int_t offs;

    _linearmodel_clear(lm);

    ae_vector_set_length(&lm->w, 4+nvars+1, _state);
    offs = 4;
    lm->w.ptr.p_double[0] = 4+nvars+1;
    lm->w.ptr.p_double[1] = linreg_lrvnum;
    lm->w.ptr.p_double[2] = nvars;
    lm->w.ptr.p_double[3] = offs;
    ae_v_move(&lm->w.ptr.p_double[offs], 1, &v->ptr.p_double[0], 1, ae_v_len(offs,offs+nvars));
}


/*************************************************************************
Procesing

INPUT PARAMETERS:
    LM      -   linear model
    X       -   input vector,  array[0..NVars-1].

Result:
    value of linear model regression estimate

  -- ALGLIB --
     Copyright 03.09.2008 by Bochkanov Sergey
*************************************************************************/
double lrprocess(linearmodel* lm,
     /* Real    */ ae_vector* x,
     ae_state *_state)
{
    double v;
    ae_int_t offs;
    ae_int_t nvars;
    double result;


    ae_assert(ae_round(lm->w.ptr.p_double[1], _state)==linreg_lrvnum, "LINREG: Incorrect LINREG version!", _state);
    nvars = ae_round(lm->w.ptr.p_double[2], _state);
    offs = ae_round(lm->w.ptr.p_double[3], _state);
    v = ae_v_dotproduct(&x->ptr.p_double[0], 1, &lm->w.ptr.p_double[offs], 1, ae_v_len(0,nvars-1));
    result = v+lm->w.ptr.p_double[offs+nvars];
    return result;
}


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    LM      -   linear model
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double lrrmserror(linearmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    ae_int_t i;
    double v;
    ae_int_t offs;
    ae_int_t nvars;
    double result;


    ae_assert(ae_round(lm->w.ptr.p_double[1], _state)==linreg_lrvnum, "LINREG: Incorrect LINREG version!", _state);
    nvars = ae_round(lm->w.ptr.p_double[2], _state);
    offs = ae_round(lm->w.ptr.p_double[3], _state);
    result = 0;
    for(i=0; i<=npoints-1; i++)
    {
        v = ae_v_dotproduct(&xy->ptr.pp_double[i][0], 1, &lm->w.ptr.p_double[offs], 1, ae_v_len(0,nvars-1));
        v = v+lm->w.ptr.p_double[offs+nvars];
        result = result+ae_sqr(v-xy->ptr.pp_double[i][nvars], _state);
    }
    result = ae_sqrt(result/npoints, _state);
    return result;
}


/*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    LM      -   linear model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average error.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double lravgerror(linearmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    ae_int_t i;
    double v;
    ae_int_t offs;
    ae_int_t nvars;
    double result;


    ae_assert(ae_round(lm->w.ptr.p_double[1], _state)==linreg_lrvnum, "LINREG: Incorrect LINREG version!", _state);
    nvars = ae_round(lm->w.ptr.p_double[2], _state);
    offs = ae_round(lm->w.ptr.p_double[3], _state);
    result = 0;
    for(i=0; i<=npoints-1; i++)
    {
        v = ae_v_dotproduct(&xy->ptr.pp_double[i][0], 1, &lm->w.ptr.p_double[offs], 1, ae_v_len(0,nvars-1));
        v = v+lm->w.ptr.p_double[offs+nvars];
        result = result+ae_fabs(v-xy->ptr.pp_double[i][nvars], _state);
    }
    result = result/npoints;
    return result;
}


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    LM      -   linear model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average relative error.

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double lravgrelerror(linearmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t k;
    double v;
    ae_int_t offs;
    ae_int_t nvars;
    double result;


    ae_assert(ae_round(lm->w.ptr.p_double[1], _state)==linreg_lrvnum, "LINREG: Incorrect LINREG version!", _state);
    nvars = ae_round(lm->w.ptr.p_double[2], _state);
    offs = ae_round(lm->w.ptr.p_double[3], _state);
    result = 0;
    k = 0;
    for(i=0; i<=npoints-1; i++)
    {
        if( ae_fp_neq(xy->ptr.pp_double[i][nvars],0) )
        {
            v = ae_v_dotproduct(&xy->ptr.pp_double[i][0], 1, &lm->w.ptr.p_double[offs], 1, ae_v_len(0,nvars-1));
            v = v+lm->w.ptr.p_double[offs+nvars];
            result = result+ae_fabs((v-xy->ptr.pp_double[i][nvars])/xy->ptr.pp_double[i][nvars], _state);
            k = k+1;
        }
    }
    if( k!=0 )
    {
        result = result/k;
    }
    return result;
}


/*************************************************************************
Copying of LinearModel strucure

INPUT PARAMETERS:
    LM1 -   original

OUTPUT PARAMETERS:
    LM2 -   copy

  -- ALGLIB --
     Copyright 15.03.2009 by Bochkanov Sergey
*************************************************************************/
void lrcopy(linearmodel* lm1, linearmodel* lm2, ae_state *_state)
{
    ae_int_t k;

    _linearmodel_clear(lm2);

    k = ae_round(lm1->w.ptr.p_double[0], _state);
    ae_vector_set_length(&lm2->w, k-1+1, _state);
    ae_v_move(&lm2->w.ptr.p_double[0], 1, &lm1->w.ptr.p_double[0], 1, ae_v_len(0,k-1));
}


void lrlines(/* Real    */ ae_matrix* xy,
     /* Real    */ ae_vector* s,
     ae_int_t n,
     ae_int_t* info,
     double* a,
     double* b,
     double* vara,
     double* varb,
     double* covab,
     double* corrab,
     double* p,
     ae_state *_state)
{
    ae_int_t i;
    double ss;
    double sx;
    double sxx;
    double sy;
    double stt;
    double e1;
    double e2;
    double t;
    double chi2;

    *info = 0;
    *a = 0;
    *b = 0;
    *vara = 0;
    *varb = 0;
    *covab = 0;
    *corrab = 0;
    *p = 0;

    if( n<2 )
    {
        *info = -1;
        return;
    }
    for(i=0; i<=n-1; i++)
    {
        if( ae_fp_less_eq(s->ptr.p_double[i],0) )
        {
            *info = -2;
            return;
        }
    }
    *info = 1;
    
    /*
     * Calculate S, SX, SY, SXX
     */
    ss = 0;
    sx = 0;
    sy = 0;
    sxx = 0;
    for(i=0; i<=n-1; i++)
    {
        t = ae_sqr(s->ptr.p_double[i], _state);
        ss = ss+1/t;
        sx = sx+xy->ptr.pp_double[i][0]/t;
        sy = sy+xy->ptr.pp_double[i][1]/t;
        sxx = sxx+ae_sqr(xy->ptr.pp_double[i][0], _state)/t;
    }
    
    /*
     * Test for condition number
     */
    t = ae_sqrt(4*ae_sqr(sx, _state)+ae_sqr(ss-sxx, _state), _state);
    e1 = 0.5*(ss+sxx+t);
    e2 = 0.5*(ss+sxx-t);
    if( ae_fp_less_eq(ae_minreal(e1, e2, _state),1000*ae_machineepsilon*ae_maxreal(e1, e2, _state)) )
    {
        *info = -3;
        return;
    }
    
    /*
     * Calculate A, B
     */
    *a = 0;
    *b = 0;
    stt = 0;
    for(i=0; i<=n-1; i++)
    {
        t = (xy->ptr.pp_double[i][0]-sx/ss)/s->ptr.p_double[i];
        *b = *b+t*xy->ptr.pp_double[i][1]/s->ptr.p_double[i];
        stt = stt+ae_sqr(t, _state);
    }
    *b = *b/stt;
    *a = (sy-sx*(*b))/ss;
    
    /*
     * Calculate goodness-of-fit
     */
    if( n>2 )
    {
        chi2 = 0;
        for(i=0; i<=n-1; i++)
        {
            chi2 = chi2+ae_sqr((xy->ptr.pp_double[i][1]-(*a)-*b*xy->ptr.pp_double[i][0])/s->ptr.p_double[i], _state);
        }
        *p = incompletegammac((double)(n-2)/(double)2, chi2/2, _state);
    }
    else
    {
        *p = 1;
    }
    
    /*
     * Calculate other parameters
     */
    *vara = (1+ae_sqr(sx, _state)/(ss*stt))/ss;
    *varb = 1/stt;
    *covab = -sx/(ss*stt);
    *corrab = *covab/ae_sqrt(*vara*(*varb), _state);
}


void lrline(/* Real    */ ae_matrix* xy,
     ae_int_t n,
     ae_int_t* info,
     double* a,
     double* b,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector s;
    ae_int_t i;
    double vara;
    double varb;
    double covab;
    double corrab;
    double p;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    *a = 0;
    *b = 0;
    ae_vector_init(&s, 0, DT_REAL, _state, ae_true);

    if( n<2 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    ae_vector_set_length(&s, n-1+1, _state);
    for(i=0; i<=n-1; i++)
    {
        s.ptr.p_double[i] = 1;
    }
    lrlines(xy, &s, n, info, a, b, &vara, &varb, &covab, &corrab, &p, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Internal linear regression subroutine
*************************************************************************/
static void linreg_lrinternal(/* Real    */ ae_matrix* xy,
     /* Real    */ ae_vector* s,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     linearmodel* lm,
     lrreport* ar,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix a;
    ae_matrix u;
    ae_matrix vt;
    ae_matrix vm;
    ae_matrix xym;
    ae_vector b;
    ae_vector sv;
    ae_vector t;
    ae_vector svi;
    ae_vector work;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t ncv;
    ae_int_t na;
    ae_int_t nacv;
    double r;
    double p;
    double epstol;
    lrreport ar2;
    ae_int_t offs;
    linearmodel tlm;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _linearmodel_clear(lm);
    _lrreport_clear(ar);
    ae_matrix_init(&a, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&u, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&vt, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&vm, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&xym, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&b, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&sv, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&t, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&svi, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&work, 0, DT_REAL, _state, ae_true);
    _lrreport_init(&ar2, _state, ae_true);
    _linearmodel_init(&tlm, _state, ae_true);

    epstol = 1000;
    
    /*
     * Check for errors in data
     */
    if( npoints<nvars||nvars<1 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    for(i=0; i<=npoints-1; i++)
    {
        if( ae_fp_less_eq(s->ptr.p_double[i],0) )
        {
            *info = -2;
            ae_frame_leave(_state);
            return;
        }
    }
    *info = 1;
    
    /*
     * Create design matrix
     */
    ae_matrix_set_length(&a, npoints-1+1, nvars-1+1, _state);
    ae_vector_set_length(&b, npoints-1+1, _state);
    for(i=0; i<=npoints-1; i++)
    {
        r = 1/s->ptr.p_double[i];
        ae_v_moved(&a.ptr.pp_double[i][0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1), r);
        b.ptr.p_double[i] = xy->ptr.pp_double[i][nvars]/s->ptr.p_double[i];
    }
    
    /*
     * Allocate W:
     * W[0]     array size
     * W[1]     version number, 0
     * W[2]     NVars (minus 1, to be compatible with external representation)
     * W[3]     coefficients offset
     */
    ae_vector_set_length(&lm->w, 4+nvars-1+1, _state);
    offs = 4;
    lm->w.ptr.p_double[0] = 4+nvars;
    lm->w.ptr.p_double[1] = linreg_lrvnum;
    lm->w.ptr.p_double[2] = nvars-1;
    lm->w.ptr.p_double[3] = offs;
    
    /*
     * Solve problem using SVD:
     *
     * 0. check for degeneracy (different types)
     * 1. A = U*diag(sv)*V'
     * 2. T = b'*U
     * 3. w = SUM((T[i]/sv[i])*V[..,i])
     * 4. cov(wi,wj) = SUM(Vji*Vjk/sv[i]^2,K=1..M)
     *
     * see $15.4 of "Numerical Recipes in C" for more information
     */
    ae_vector_set_length(&t, nvars-1+1, _state);
    ae_vector_set_length(&svi, nvars-1+1, _state);
    ae_matrix_set_length(&ar->c, nvars-1+1, nvars-1+1, _state);
    ae_matrix_set_length(&vm, nvars-1+1, nvars-1+1, _state);
    if( !rmatrixsvd(&a, npoints, nvars, 1, 1, 2, &sv, &u, &vt, _state) )
    {
        *info = -4;
        ae_frame_leave(_state);
        return;
    }
    if( ae_fp_less_eq(sv.ptr.p_double[0],0) )
    {
        
        /*
         * Degenerate case: zero design matrix.
         */
        for(i=offs; i<=offs+nvars-1; i++)
        {
            lm->w.ptr.p_double[i] = 0;
        }
        ar->rmserror = lrrmserror(lm, xy, npoints, _state);
        ar->avgerror = lravgerror(lm, xy, npoints, _state);
        ar->avgrelerror = lravgrelerror(lm, xy, npoints, _state);
        ar->cvrmserror = ar->rmserror;
        ar->cvavgerror = ar->avgerror;
        ar->cvavgrelerror = ar->avgrelerror;
        ar->ncvdefects = 0;
        ae_vector_set_length(&ar->cvdefects, nvars-1+1, _state);
        ae_matrix_set_length(&ar->c, nvars-1+1, nvars-1+1, _state);
        for(i=0; i<=nvars-1; i++)
        {
            for(j=0; j<=nvars-1; j++)
            {
                ar->c.ptr.pp_double[i][j] = 0;
            }
        }
        ae_frame_leave(_state);
        return;
    }
    if( ae_fp_less_eq(sv.ptr.p_double[nvars-1],epstol*ae_machineepsilon*sv.ptr.p_double[0]) )
    {
        
        /*
         * Degenerate case, non-zero design matrix.
         *
         * We can leave it and solve task in SVD least squares fashion.
         * Solution and covariance matrix will be obtained correctly,
         * but CV error estimates - will not. It is better to reduce
         * it to non-degenerate task and to obtain correct CV estimates.
         */
        for(k=nvars; k>=1; k--)
        {
            if( ae_fp_greater(sv.ptr.p_double[k-1],epstol*ae_machineepsilon*sv.ptr.p_double[0]) )
            {
                
                /*
                 * Reduce
                 */
                ae_matrix_set_length(&xym, npoints-1+1, k+1, _state);
                for(i=0; i<=npoints-1; i++)
                {
                    for(j=0; j<=k-1; j++)
                    {
                        r = ae_v_dotproduct(&xy->ptr.pp_double[i][0], 1, &vt.ptr.pp_double[j][0], 1, ae_v_len(0,nvars-1));
                        xym.ptr.pp_double[i][j] = r;
                    }
                    xym.ptr.pp_double[i][k] = xy->ptr.pp_double[i][nvars];
                }
                
                /*
                 * Solve
                 */
                linreg_lrinternal(&xym, s, npoints, k, info, &tlm, &ar2, _state);
                if( *info!=1 )
                {
                    ae_frame_leave(_state);
                    return;
                }
                
                /*
                 * Convert back to un-reduced format
                 */
                for(j=0; j<=nvars-1; j++)
                {
                    lm->w.ptr.p_double[offs+j] = 0;
                }
                for(j=0; j<=k-1; j++)
                {
                    r = tlm.w.ptr.p_double[offs+j];
                    ae_v_addd(&lm->w.ptr.p_double[offs], 1, &vt.ptr.pp_double[j][0], 1, ae_v_len(offs,offs+nvars-1), r);
                }
                ar->rmserror = ar2.rmserror;
                ar->avgerror = ar2.avgerror;
                ar->avgrelerror = ar2.avgrelerror;
                ar->cvrmserror = ar2.cvrmserror;
                ar->cvavgerror = ar2.cvavgerror;
                ar->cvavgrelerror = ar2.cvavgrelerror;
                ar->ncvdefects = ar2.ncvdefects;
                ae_vector_set_length(&ar->cvdefects, nvars-1+1, _state);
                for(j=0; j<=ar->ncvdefects-1; j++)
                {
                    ar->cvdefects.ptr.p_int[j] = ar2.cvdefects.ptr.p_int[j];
                }
                ae_matrix_set_length(&ar->c, nvars-1+1, nvars-1+1, _state);
                ae_vector_set_length(&work, nvars+1, _state);
                matrixmatrixmultiply(&ar2.c, 0, k-1, 0, k-1, ae_false, &vt, 0, k-1, 0, nvars-1, ae_false, 1.0, &vm, 0, k-1, 0, nvars-1, 0.0, &work, _state);
                matrixmatrixmultiply(&vt, 0, k-1, 0, nvars-1, ae_true, &vm, 0, k-1, 0, nvars-1, ae_false, 1.0, &ar->c, 0, nvars-1, 0, nvars-1, 0.0, &work, _state);
                ae_frame_leave(_state);
                return;
            }
        }
        *info = -255;
        ae_frame_leave(_state);
        return;
    }
    for(i=0; i<=nvars-1; i++)
    {
        if( ae_fp_greater(sv.ptr.p_double[i],epstol*ae_machineepsilon*sv.ptr.p_double[0]) )
        {
            svi.ptr.p_double[i] = 1/sv.ptr.p_double[i];
        }
        else
        {
            svi.ptr.p_double[i] = 0;
        }
    }
    for(i=0; i<=nvars-1; i++)
    {
        t.ptr.p_double[i] = 0;
    }
    for(i=0; i<=npoints-1; i++)
    {
        r = b.ptr.p_double[i];
        ae_v_addd(&t.ptr.p_double[0], 1, &u.ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1), r);
    }
    for(i=0; i<=nvars-1; i++)
    {
        lm->w.ptr.p_double[offs+i] = 0;
    }
    for(i=0; i<=nvars-1; i++)
    {
        r = t.ptr.p_double[i]*svi.ptr.p_double[i];
        ae_v_addd(&lm->w.ptr.p_double[offs], 1, &vt.ptr.pp_double[i][0], 1, ae_v_len(offs,offs+nvars-1), r);
    }
    for(j=0; j<=nvars-1; j++)
    {
        r = svi.ptr.p_double[j];
        ae_v_moved(&vm.ptr.pp_double[0][j], vm.stride, &vt.ptr.pp_double[j][0], 1, ae_v_len(0,nvars-1), r);
    }
    for(i=0; i<=nvars-1; i++)
    {
        for(j=i; j<=nvars-1; j++)
        {
            r = ae_v_dotproduct(&vm.ptr.pp_double[i][0], 1, &vm.ptr.pp_double[j][0], 1, ae_v_len(0,nvars-1));
            ar->c.ptr.pp_double[i][j] = r;
            ar->c.ptr.pp_double[j][i] = r;
        }
    }
    
    /*
     * Leave-1-out cross-validation error.
     *
     * NOTATIONS:
     * A            design matrix
     * A*x = b      original linear least squares task
     * U*S*V'       SVD of A
     * ai           i-th row of the A
     * bi           i-th element of the b
     * xf           solution of the original LLS task
     *
     * Cross-validation error of i-th element from a sample is
     * calculated using following formula:
     *
     *     ERRi = ai*xf - (ai*xf-bi*(ui*ui'))/(1-ui*ui')     (1)
     *
     * This formula can be derived from normal equations of the
     * original task
     *
     *     (A'*A)x = A'*b                                    (2)
     *
     * by applying modification (zeroing out i-th row of A) to (2):
     *
     *     (A-ai)'*(A-ai) = (A-ai)'*b
     *
     * and using Sherman-Morrison formula for updating matrix inverse
     *
     * NOTE 1: b is not zeroed out since it is much simpler and
     * does not influence final result.
     *
     * NOTE 2: some design matrices A have such ui that 1-ui*ui'=0.
     * Formula (1) can't be applied for such cases and they are skipped
     * from CV calculation (which distorts resulting CV estimate).
     * But from the properties of U we can conclude that there can
     * be no more than NVars such vectors. Usually
     * NVars << NPoints, so in a normal case it only slightly
     * influences result.
     */
    ncv = 0;
    na = 0;
    nacv = 0;
    ar->rmserror = 0;
    ar->avgerror = 0;
    ar->avgrelerror = 0;
    ar->cvrmserror = 0;
    ar->cvavgerror = 0;
    ar->cvavgrelerror = 0;
    ar->ncvdefects = 0;
    ae_vector_set_length(&ar->cvdefects, nvars-1+1, _state);
    for(i=0; i<=npoints-1; i++)
    {
        
        /*
         * Error on a training set
         */
        r = ae_v_dotproduct(&xy->ptr.pp_double[i][0], 1, &lm->w.ptr.p_double[offs], 1, ae_v_len(0,nvars-1));
        ar->rmserror = ar->rmserror+ae_sqr(r-xy->ptr.pp_double[i][nvars], _state);
        ar->avgerror = ar->avgerror+ae_fabs(r-xy->ptr.pp_double[i][nvars], _state);
        if( ae_fp_neq(xy->ptr.pp_double[i][nvars],0) )
        {
            ar->avgrelerror = ar->avgrelerror+ae_fabs((r-xy->ptr.pp_double[i][nvars])/xy->ptr.pp_double[i][nvars], _state);
            na = na+1;
        }
        
        /*
         * Error using fast leave-one-out cross-validation
         */
        p = ae_v_dotproduct(&u.ptr.pp_double[i][0], 1, &u.ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
        if( ae_fp_greater(p,1-epstol*ae_machineepsilon) )
        {
            ar->cvdefects.ptr.p_int[ar->ncvdefects] = i;
            ar->ncvdefects = ar->ncvdefects+1;
            continue;
        }
        r = s->ptr.p_double[i]*(r/s->ptr.p_double[i]-b.ptr.p_double[i]*p)/(1-p);
        ar->cvrmserror = ar->cvrmserror+ae_sqr(r-xy->ptr.pp_double[i][nvars], _state);
        ar->cvavgerror = ar->cvavgerror+ae_fabs(r-xy->ptr.pp_double[i][nvars], _state);
        if( ae_fp_neq(xy->ptr.pp_double[i][nvars],0) )
        {
            ar->cvavgrelerror = ar->cvavgrelerror+ae_fabs((r-xy->ptr.pp_double[i][nvars])/xy->ptr.pp_double[i][nvars], _state);
            nacv = nacv+1;
        }
        ncv = ncv+1;
    }
    if( ncv==0 )
    {
        
        /*
         * Something strange: ALL ui are degenerate.
         * Unexpected...
         */
        *info = -255;
        ae_frame_leave(_state);
        return;
    }
    ar->rmserror = ae_sqrt(ar->rmserror/npoints, _state);
    ar->avgerror = ar->avgerror/npoints;
    if( na!=0 )
    {
        ar->avgrelerror = ar->avgrelerror/na;
    }
    ar->cvrmserror = ae_sqrt(ar->cvrmserror/ncv, _state);
    ar->cvavgerror = ar->cvavgerror/ncv;
    if( nacv!=0 )
    {
        ar->cvavgrelerror = ar->cvavgrelerror/nacv;
    }
    ae_frame_leave(_state);
}


ae_bool _linearmodel_init(linearmodel* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->w, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _linearmodel_init_copy(linearmodel* dst, linearmodel* src, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init_copy(&dst->w, &src->w, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _linearmodel_clear(linearmodel* p)
{
    ae_vector_clear(&p->w);
}


ae_bool _lrreport_init(lrreport* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_matrix_init(&p->c, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->cvdefects, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _lrreport_init_copy(lrreport* dst, lrreport* src, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_matrix_init_copy(&dst->c, &src->c, _state, make_automatic) )
        return ae_false;
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
    dst->cvrmserror = src->cvrmserror;
    dst->cvavgerror = src->cvavgerror;
    dst->cvavgrelerror = src->cvavgrelerror;
    dst->ncvdefects = src->ncvdefects;
    if( !ae_vector_init_copy(&dst->cvdefects, &src->cvdefects, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _lrreport_clear(lrreport* p)
{
    ae_matrix_clear(&p->c);
    ae_vector_clear(&p->cvdefects);
}




/*************************************************************************
Filters: simple moving averages (unsymmetric).

This filter replaces array by results of SMA(K) filter. SMA(K) is defined
as filter which averages at most K previous points (previous - not points
AROUND central point) - or less, in case of the first K-1 points.

INPUT PARAMETERS:
    X           -   array[N], array to process. It can be larger than N,
                    in this case only first N points are processed.
    N           -   points count, N>=0
    K           -   K>=1 (K can be larger than N ,  such  cases  will  be
                    correctly handled). Window width. K=1 corresponds  to
                    identity transformation (nothing changes).

OUTPUT PARAMETERS:
    X           -   array, whose first N elements were processed with SMA(K)

NOTE 1: this function uses efficient in-place  algorithm  which  does not
        allocate temporary arrays.

NOTE 2: this algorithm makes only one pass through array and uses running
        sum  to speed-up calculation of the averages. Additional measures
        are taken to ensure that running sum on a long sequence  of  zero
        elements will be correctly reset to zero even in the presence  of
        round-off error.

NOTE 3: this  is  unsymmetric version of the algorithm,  which  does  NOT
        averages points after the current one. Only X[i], X[i-1], ... are
        used when calculating new value of X[i]. We should also note that
        this algorithm uses BOTH previous points and  current  one,  i.e.
        new value of X[i] depends on BOTH previous point and X[i] itself.

  -- ALGLIB --
     Copyright 25.10.2011 by Bochkanov Sergey
*************************************************************************/
void filtersma(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_int_t k,
     ae_state *_state)
{
    ae_int_t i;
    double runningsum;
    double termsinsum;
    ae_int_t zeroprefix;
    double v;


    ae_assert(n>=0, "FilterSMA: N<0", _state);
    ae_assert(x->cnt>=n, "FilterSMA: Length(X)<N", _state);
    ae_assert(isfinitevector(x, n, _state), "FilterSMA: X contains INF or NAN", _state);
    ae_assert(k>=1, "FilterSMA: K<1", _state);
    
    /*
     * Quick exit, if necessary
     */
    if( n<=1||k==1 )
    {
        return;
    }
    
    /*
     * Prepare variables (see below for explanation)
     */
    runningsum = 0.0;
    termsinsum = 0;
    for(i=ae_maxint(n-k, 0, _state); i<=n-1; i++)
    {
        runningsum = runningsum+x->ptr.p_double[i];
        termsinsum = termsinsum+1;
    }
    i = ae_maxint(n-k, 0, _state);
    zeroprefix = 0;
    while(i<=n-1&&ae_fp_eq(x->ptr.p_double[i],0))
    {
        zeroprefix = zeroprefix+1;
        i = i+1;
    }
    
    /*
     * General case: we assume that N>1 and K>1
     *
     * Make one pass through all elements. At the beginning of
     * the iteration we have:
     * * I              element being processed
     * * RunningSum     current value of the running sum
     *                  (including I-th element)
     * * TermsInSum     number of terms in sum, 0<=TermsInSum<=K
     * * ZeroPrefix     length of the sequence of zero elements
     *                  which starts at X[I-K+1] and continues towards X[I].
     *                  Equal to zero in case X[I-K+1] is non-zero.
     *                  This value is used to make RunningSum exactly zero
     *                  when it follows from the problem properties.
     */
    for(i=n-1; i>=0; i--)
    {
        
        /*
         * Store new value of X[i], save old value in V
         */
        v = x->ptr.p_double[i];
        x->ptr.p_double[i] = runningsum/termsinsum;
        
        /*
         * Update RunningSum and TermsInSum
         */
        if( i-k>=0 )
        {
            runningsum = runningsum-v+x->ptr.p_double[i-k];
        }
        else
        {
            runningsum = runningsum-v;
            termsinsum = termsinsum-1;
        }
        
        /*
         * Update ZeroPrefix.
         * In case we have ZeroPrefix=TermsInSum,
         * RunningSum is reset to zero.
         */
        if( i-k>=0 )
        {
            if( ae_fp_neq(x->ptr.p_double[i-k],0) )
            {
                zeroprefix = 0;
            }
            else
            {
                zeroprefix = ae_minint(zeroprefix+1, k, _state);
            }
        }
        else
        {
            zeroprefix = ae_minint(zeroprefix, i+1, _state);
        }
        if( ae_fp_eq(zeroprefix,termsinsum) )
        {
            runningsum = 0;
        }
    }
}


/*************************************************************************
Filters: exponential moving averages.

This filter replaces array by results of EMA(alpha) filter. EMA(alpha) is
defined as filter which replaces X[] by S[]:
    S[0] = X[0]
    S[t] = alpha*X[t] + (1-alpha)*S[t-1]

INPUT PARAMETERS:
    X           -   array[N], array to process. It can be larger than N,
                    in this case only first N points are processed.
    N           -   points count, N>=0
    alpha       -   0<alpha<=1, smoothing parameter.

OUTPUT PARAMETERS:
    X           -   array, whose first N elements were processed
                    with EMA(alpha)

NOTE 1: this function uses efficient in-place  algorithm  which  does not
        allocate temporary arrays.

NOTE 2: this algorithm uses BOTH previous points and  current  one,  i.e.
        new value of X[i] depends on BOTH previous point and X[i] itself.

NOTE 3: technical analytis users quite often work  with  EMA  coefficient
        expressed in DAYS instead of fractions. If you want to  calculate
        EMA(N), where N is a number of days, you can use alpha=2/(N+1).

  -- ALGLIB --
     Copyright 25.10.2011 by Bochkanov Sergey
*************************************************************************/
void filterema(/* Real    */ ae_vector* x,
     ae_int_t n,
     double alpha,
     ae_state *_state)
{
    ae_int_t i;


    ae_assert(n>=0, "FilterEMA: N<0", _state);
    ae_assert(x->cnt>=n, "FilterEMA: Length(X)<N", _state);
    ae_assert(isfinitevector(x, n, _state), "FilterEMA: X contains INF or NAN", _state);
    ae_assert(ae_fp_greater(alpha,0), "FilterEMA: Alpha<=0", _state);
    ae_assert(ae_fp_less_eq(alpha,1), "FilterEMA: Alpha>1", _state);
    
    /*
     * Quick exit, if necessary
     */
    if( n<=1||ae_fp_eq(alpha,1) )
    {
        return;
    }
    
    /*
     * Process
     */
    for(i=1; i<=n-1; i++)
    {
        x->ptr.p_double[i] = alpha*x->ptr.p_double[i]+(1-alpha)*x->ptr.p_double[i-1];
    }
}


/*************************************************************************
Filters: linear regression moving averages.

This filter replaces array by results of LRMA(K) filter.

LRMA(K) is defined as filter which, for each data  point,  builds  linear
regression  model  using  K  prevous  points (point itself is included in
these K points) and calculates value of this linear model at the point in
question.

INPUT PARAMETERS:
    X           -   array[N], array to process. It can be larger than N,
                    in this case only first N points are processed.
    N           -   points count, N>=0
    K           -   K>=1 (K can be larger than N ,  such  cases  will  be
                    correctly handled). Window width. K=1 corresponds  to
                    identity transformation (nothing changes).

OUTPUT PARAMETERS:
    X           -   array, whose first N elements were processed with SMA(K)

NOTE 1: this function uses efficient in-place  algorithm  which  does not
        allocate temporary arrays.

NOTE 2: this algorithm makes only one pass through array and uses running
        sum  to speed-up calculation of the averages. Additional measures
        are taken to ensure that running sum on a long sequence  of  zero
        elements will be correctly reset to zero even in the presence  of
        round-off error.

NOTE 3: this  is  unsymmetric version of the algorithm,  which  does  NOT
        averages points after the current one. Only X[i], X[i-1], ... are
        used when calculating new value of X[i]. We should also note that
        this algorithm uses BOTH previous points and  current  one,  i.e.
        new value of X[i] depends on BOTH previous point and X[i] itself.

  -- ALGLIB --
     Copyright 25.10.2011 by Bochkanov Sergey
*************************************************************************/
void filterlrma(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_int_t k,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t m;
    ae_matrix xy;
    ae_vector s;
    ae_int_t info;
    double a;
    double b;
    double vara;
    double varb;
    double covab;
    double corrab;
    double p;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init(&xy, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&s, 0, DT_REAL, _state, ae_true);

    ae_assert(n>=0, "FilterLRMA: N<0", _state);
    ae_assert(x->cnt>=n, "FilterLRMA: Length(X)<N", _state);
    ae_assert(isfinitevector(x, n, _state), "FilterLRMA: X contains INF or NAN", _state);
    ae_assert(k>=1, "FilterLRMA: K<1", _state);
    
    /*
     * Quick exit, if necessary:
     * * either N is equal to 1 (nothing to average)
     * * or K is 1 (only point itself is used) or 2 (model is too simple,
     *   we will always get identity transformation)
     */
    if( n<=1||k<=2 )
    {
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * General case: K>2, N>1.
     * We do not process points with I<2 because first two points (I=0 and I=1) will be
     * left unmodified by LRMA filter in any case.
     */
    ae_matrix_set_length(&xy, k, 2, _state);
    ae_vector_set_length(&s, k, _state);
    for(i=0; i<=k-1; i++)
    {
        xy.ptr.pp_double[i][0] = i;
        s.ptr.p_double[i] = 1.0;
    }
    for(i=n-1; i>=2; i--)
    {
        m = ae_minint(i+1, k, _state);
        ae_v_move(&xy.ptr.pp_double[0][1], xy.stride, &x->ptr.p_double[i-m+1], 1, ae_v_len(0,m-1));
        lrlines(&xy, &s, m, &info, &a, &b, &vara, &varb, &covab, &corrab, &p, _state);
        ae_assert(info==1, "FilterLRMA: internal error", _state);
        x->ptr.p_double[i] = a+b*(m-1);
    }
    ae_frame_leave(_state);
}




/*************************************************************************
k-means++ clusterization

INPUT PARAMETERS:
    XY          -   dataset, array [0..NPoints-1,0..NVars-1].
    NPoints     -   dataset size, NPoints>=K
    NVars       -   number of variables, NVars>=1
    K           -   desired number of clusters, K>=1
    Restarts    -   number of restarts, Restarts>=1

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -3, if task is degenerate (number of distinct points is
                          less than K)
                    * -1, if incorrect NPoints/NFeatures/K/Restarts was passed
                    *  1, if subroutine finished successfully
    C           -   array[0..NVars-1,0..K-1].matrix whose columns store
                    cluster's centers
    XYC         -   array[NPoints], which contains cluster indexes

  -- ALGLIB --
     Copyright 21.03.2009 by Bochkanov Sergey
*************************************************************************/
void kmeansgenerate(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t k,
     ae_int_t restarts,
     ae_int_t* info,
     /* Real    */ ae_matrix* c,
     /* Integer */ ae_vector* xyc,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_matrix ct;
    ae_matrix ctbest;
    ae_vector xycbest;
    double e;
    double ebest;
    ae_vector x;
    ae_vector tmp;
    ae_vector d2;
    ae_vector p;
    ae_vector csizes;
    ae_vector cbusy;
    double v;
    ae_int_t cclosest;
    double dclosest;
    ae_vector work;
    ae_bool waschanges;
    ae_bool zerosizeclusters;
    ae_int_t pass;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    ae_matrix_clear(c);
    ae_vector_clear(xyc);
    ae_matrix_init(&ct, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&ctbest, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&xycbest, 0, DT_INT, _state, ae_true);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tmp, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&d2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&p, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&csizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&cbusy, 0, DT_BOOL, _state, ae_true);
    ae_vector_init(&work, 0, DT_REAL, _state, ae_true);

    
    /*
     * Test parameters
     */
    if( ((npoints<k||nvars<1)||k<1)||restarts<1 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * TODO: special case K=1
     * TODO: special case K=NPoints
     */
    *info = 1;
    
    /*
     * Multiple passes of k-means++ algorithm
     */
    ae_matrix_set_length(&ct, k, nvars, _state);
    ae_matrix_set_length(&ctbest, k, nvars, _state);
    ae_vector_set_length(xyc, npoints, _state);
    ae_vector_set_length(&xycbest, npoints, _state);
    ae_vector_set_length(&d2, npoints, _state);
    ae_vector_set_length(&p, npoints, _state);
    ae_vector_set_length(&tmp, nvars, _state);
    ae_vector_set_length(&csizes, k, _state);
    ae_vector_set_length(&cbusy, k, _state);
    ebest = ae_maxrealnumber;
    for(pass=1; pass<=restarts; pass++)
    {
        
        /*
         * Select initial centers  using k-means++ algorithm
         * 1. Choose first center at random
         * 2. Choose next centers using their distance from centers already chosen
         *
         * Note that for performance reasons centers are stored in ROWS of CT, not
         * in columns. We'll transpose CT in the end and store it in the C.
         */
        i = ae_randominteger(npoints, _state);
        ae_v_move(&ct.ptr.pp_double[0][0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
        cbusy.ptr.p_bool[0] = ae_true;
        for(i=1; i<=k-1; i++)
        {
            cbusy.ptr.p_bool[i] = ae_false;
        }
        if( !kmeans_selectcenterpp(xy, npoints, nvars, &ct, &cbusy, k, &d2, &p, &tmp, _state) )
        {
            *info = -3;
            ae_frame_leave(_state);
            return;
        }
        
        /*
         * Update centers:
         * 2. update center positions
         */
        for(i=0; i<=npoints-1; i++)
        {
            xyc->ptr.p_int[i] = -1;
        }
        for(;;)
        {
            
            /*
             * fill XYC with center numbers
             */
            waschanges = ae_false;
            for(i=0; i<=npoints-1; i++)
            {
                cclosest = -1;
                dclosest = ae_maxrealnumber;
                for(j=0; j<=k-1; j++)
                {
                    ae_v_move(&tmp.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
                    ae_v_sub(&tmp.ptr.p_double[0], 1, &ct.ptr.pp_double[j][0], 1, ae_v_len(0,nvars-1));
                    v = ae_v_dotproduct(&tmp.ptr.p_double[0], 1, &tmp.ptr.p_double[0], 1, ae_v_len(0,nvars-1));
                    if( ae_fp_less(v,dclosest) )
                    {
                        cclosest = j;
                        dclosest = v;
                    }
                }
                if( xyc->ptr.p_int[i]!=cclosest )
                {
                    waschanges = ae_true;
                }
                xyc->ptr.p_int[i] = cclosest;
            }
            
            /*
             * Update centers
             */
            for(j=0; j<=k-1; j++)
            {
                csizes.ptr.p_int[j] = 0;
            }
            for(i=0; i<=k-1; i++)
            {
                for(j=0; j<=nvars-1; j++)
                {
                    ct.ptr.pp_double[i][j] = 0;
                }
            }
            for(i=0; i<=npoints-1; i++)
            {
                csizes.ptr.p_int[xyc->ptr.p_int[i]] = csizes.ptr.p_int[xyc->ptr.p_int[i]]+1;
                ae_v_add(&ct.ptr.pp_double[xyc->ptr.p_int[i]][0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
            }
            zerosizeclusters = ae_false;
            for(i=0; i<=k-1; i++)
            {
                cbusy.ptr.p_bool[i] = csizes.ptr.p_int[i]!=0;
                zerosizeclusters = zerosizeclusters||csizes.ptr.p_int[i]==0;
            }
            if( zerosizeclusters )
            {
                
                /*
                 * Some clusters have zero size - rare, but possible.
                 * We'll choose new centers for such clusters using k-means++ rule
                 * and restart algorithm
                 */
                if( !kmeans_selectcenterpp(xy, npoints, nvars, &ct, &cbusy, k, &d2, &p, &tmp, _state) )
                {
                    *info = -3;
                    ae_frame_leave(_state);
                    return;
                }
                continue;
            }
            for(j=0; j<=k-1; j++)
            {
                v = (double)1/(double)csizes.ptr.p_int[j];
                ae_v_muld(&ct.ptr.pp_double[j][0], 1, ae_v_len(0,nvars-1), v);
            }
            
            /*
             * if nothing has changed during iteration
             */
            if( !waschanges )
            {
                break;
            }
        }
        
        /*
         * 3. Calculate E, compare with best centers found so far
         */
        e = 0;
        for(i=0; i<=npoints-1; i++)
        {
            ae_v_move(&tmp.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
            ae_v_sub(&tmp.ptr.p_double[0], 1, &ct.ptr.pp_double[xyc->ptr.p_int[i]][0], 1, ae_v_len(0,nvars-1));
            v = ae_v_dotproduct(&tmp.ptr.p_double[0], 1, &tmp.ptr.p_double[0], 1, ae_v_len(0,nvars-1));
            e = e+v;
        }
        if( ae_fp_less(e,ebest) )
        {
            
            /*
             * store partition.
             */
            ebest = e;
            copymatrix(&ct, 0, k-1, 0, nvars-1, &ctbest, 0, k-1, 0, nvars-1, _state);
            for(i=0; i<=npoints-1; i++)
            {
                xycbest.ptr.p_int[i] = xyc->ptr.p_int[i];
            }
        }
    }
    
    /*
     * Copy and transpose
     */
    ae_matrix_set_length(c, nvars-1+1, k-1+1, _state);
    copyandtranspose(&ctbest, 0, k-1, 0, nvars-1, c, 0, nvars-1, 0, k-1, _state);
    for(i=0; i<=npoints-1; i++)
    {
        xyc->ptr.p_int[i] = xycbest.ptr.p_int[i];
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Select center for a new cluster using k-means++ rule
*************************************************************************/
static ae_bool kmeans_selectcenterpp(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     /* Real    */ ae_matrix* centers,
     /* Boolean */ ae_vector* busycenters,
     ae_int_t ccnt,
     /* Real    */ ae_vector* d2,
     /* Real    */ ae_vector* p,
     /* Real    */ ae_vector* tmp,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector _busycenters;
    ae_int_t i;
    ae_int_t j;
    ae_int_t cc;
    double v;
    double s;
    ae_bool result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init_copy(&_busycenters, busycenters, _state, ae_true);
    busycenters = &_busycenters;

    result = ae_true;
    for(cc=0; cc<=ccnt-1; cc++)
    {
        if( !busycenters->ptr.p_bool[cc] )
        {
            
            /*
             * fill D2
             */
            for(i=0; i<=npoints-1; i++)
            {
                d2->ptr.p_double[i] = ae_maxrealnumber;
                for(j=0; j<=ccnt-1; j++)
                {
                    if( busycenters->ptr.p_bool[j] )
                    {
                        ae_v_move(&tmp->ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
                        ae_v_sub(&tmp->ptr.p_double[0], 1, &centers->ptr.pp_double[j][0], 1, ae_v_len(0,nvars-1));
                        v = ae_v_dotproduct(&tmp->ptr.p_double[0], 1, &tmp->ptr.p_double[0], 1, ae_v_len(0,nvars-1));
                        if( ae_fp_less(v,d2->ptr.p_double[i]) )
                        {
                            d2->ptr.p_double[i] = v;
                        }
                    }
                }
            }
            
            /*
             * calculate P (non-cumulative)
             */
            s = 0;
            for(i=0; i<=npoints-1; i++)
            {
                s = s+d2->ptr.p_double[i];
            }
            if( ae_fp_eq(s,0) )
            {
                result = ae_false;
                ae_frame_leave(_state);
                return result;
            }
            s = 1/s;
            ae_v_moved(&p->ptr.p_double[0], 1, &d2->ptr.p_double[0], 1, ae_v_len(0,npoints-1), s);
            
            /*
             * choose one of points with probability P
             * random number within (0,1) is generated and
             * inverse empirical CDF is used to randomly choose a point.
             */
            s = 0;
            v = ae_randomreal(_state);
            for(i=0; i<=npoints-1; i++)
            {
                s = s+p->ptr.p_double[i];
                if( ae_fp_less_eq(v,s)||i==npoints-1 )
                {
                    ae_v_move(&centers->ptr.pp_double[cc][0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
                    busycenters->ptr.p_bool[cc] = ae_true;
                    break;
                }
            }
        }
    }
    ae_frame_leave(_state);
    return result;
}




/*************************************************************************
Multiclass Fisher LDA

Subroutine finds coefficients of linear combination which optimally separates
training set on classes.

INPUT PARAMETERS:
    XY          -   training set, array[0..NPoints-1,0..NVars].
                    First NVars columns store values of independent
                    variables, next column stores number of class (from 0
                    to NClasses-1) which dataset element belongs to. Fractional
                    values are rounded to nearest integer.
    NPoints     -   training set size, NPoints>=0
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   number of classes, NClasses>=2


OUTPUT PARAMETERS:
    Info        -   return code:
                    * -4, if internal EVD subroutine hasn't converged
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed (NPoints<0,
                          NVars<1, NClasses<2)
                    *  1, if task has been solved
                    *  2, if there was a multicollinearity in training set,
                          but task has been solved.
    W           -   linear combination coefficients, array[0..NVars-1]

  -- ALGLIB --
     Copyright 31.05.2008 by Bochkanov Sergey
*************************************************************************/
void fisherlda(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t* info,
     /* Real    */ ae_vector* w,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix w2;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    ae_vector_clear(w);
    ae_matrix_init(&w2, 0, 0, DT_REAL, _state, ae_true);

    fisherldan(xy, npoints, nvars, nclasses, info, &w2, _state);
    if( *info>0 )
    {
        ae_vector_set_length(w, nvars-1+1, _state);
        ae_v_move(&w->ptr.p_double[0], 1, &w2.ptr.pp_double[0][0], w2.stride, ae_v_len(0,nvars-1));
    }
    ae_frame_leave(_state);
}


/*************************************************************************
N-dimensional multiclass Fisher LDA

Subroutine finds coefficients of linear combinations which optimally separates
training set on classes. It returns N-dimensional basis whose vector are sorted
by quality of training set separation (in descending order).

INPUT PARAMETERS:
    XY          -   training set, array[0..NPoints-1,0..NVars].
                    First NVars columns store values of independent
                    variables, next column stores number of class (from 0
                    to NClasses-1) which dataset element belongs to. Fractional
                    values are rounded to nearest integer.
    NPoints     -   training set size, NPoints>=0
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   number of classes, NClasses>=2


OUTPUT PARAMETERS:
    Info        -   return code:
                    * -4, if internal EVD subroutine hasn't converged
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed (NPoints<0,
                          NVars<1, NClasses<2)
                    *  1, if task has been solved
                    *  2, if there was a multicollinearity in training set,
                          but task has been solved.
    W           -   basis, array[0..NVars-1,0..NVars-1]
                    columns of matrix stores basis vectors, sorted by
                    quality of training set separation (in descending order)

  -- ALGLIB --
     Copyright 31.05.2008 by Bochkanov Sergey
*************************************************************************/
void fisherldan(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t* info,
     /* Real    */ ae_matrix* w,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t m;
    double v;
    ae_vector c;
    ae_vector mu;
    ae_matrix muc;
    ae_vector nc;
    ae_matrix sw;
    ae_matrix st;
    ae_matrix z;
    ae_matrix z2;
    ae_matrix tm;
    ae_matrix sbroot;
    ae_matrix a;
    ae_matrix xyproj;
    ae_matrix wproj;
    ae_vector tf;
    ae_vector d;
    ae_vector d2;
    ae_vector work;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    ae_matrix_clear(w);
    ae_vector_init(&c, 0, DT_INT, _state, ae_true);
    ae_vector_init(&mu, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&muc, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&nc, 0, DT_INT, _state, ae_true);
    ae_matrix_init(&sw, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&st, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&z, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&z2, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&tm, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&sbroot, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&a, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&xyproj, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&wproj, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&tf, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&d, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&d2, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&work, 0, DT_REAL, _state, ae_true);

    
    /*
     * Test data
     */
    if( (npoints<0||nvars<1)||nclasses<2 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    for(i=0; i<=npoints-1; i++)
    {
        if( ae_round(xy->ptr.pp_double[i][nvars], _state)<0||ae_round(xy->ptr.pp_double[i][nvars], _state)>=nclasses )
        {
            *info = -2;
            ae_frame_leave(_state);
            return;
        }
    }
    *info = 1;
    
    /*
     * Special case: NPoints<=1
     * Degenerate task.
     */
    if( npoints<=1 )
    {
        *info = 2;
        ae_matrix_set_length(w, nvars-1+1, nvars-1+1, _state);
        for(i=0; i<=nvars-1; i++)
        {
            for(j=0; j<=nvars-1; j++)
            {
                if( i==j )
                {
                    w->ptr.pp_double[i][j] = 1;
                }
                else
                {
                    w->ptr.pp_double[i][j] = 0;
                }
            }
        }
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Prepare temporaries
     */
    ae_vector_set_length(&tf, nvars-1+1, _state);
    ae_vector_set_length(&work, ae_maxint(nvars, npoints, _state)+1, _state);
    
    /*
     * Convert class labels from reals to integers (just for convenience)
     */
    ae_vector_set_length(&c, npoints-1+1, _state);
    for(i=0; i<=npoints-1; i++)
    {
        c.ptr.p_int[i] = ae_round(xy->ptr.pp_double[i][nvars], _state);
    }
    
    /*
     * Calculate class sizes and means
     */
    ae_vector_set_length(&mu, nvars-1+1, _state);
    ae_matrix_set_length(&muc, nclasses-1+1, nvars-1+1, _state);
    ae_vector_set_length(&nc, nclasses-1+1, _state);
    for(j=0; j<=nvars-1; j++)
    {
        mu.ptr.p_double[j] = 0;
    }
    for(i=0; i<=nclasses-1; i++)
    {
        nc.ptr.p_int[i] = 0;
        for(j=0; j<=nvars-1; j++)
        {
            muc.ptr.pp_double[i][j] = 0;
        }
    }
    for(i=0; i<=npoints-1; i++)
    {
        ae_v_add(&mu.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
        ae_v_add(&muc.ptr.pp_double[c.ptr.p_int[i]][0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
        nc.ptr.p_int[c.ptr.p_int[i]] = nc.ptr.p_int[c.ptr.p_int[i]]+1;
    }
    for(i=0; i<=nclasses-1; i++)
    {
        v = (double)1/(double)nc.ptr.p_int[i];
        ae_v_muld(&muc.ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1), v);
    }
    v = (double)1/(double)npoints;
    ae_v_muld(&mu.ptr.p_double[0], 1, ae_v_len(0,nvars-1), v);
    
    /*
     * Create ST matrix
     */
    ae_matrix_set_length(&st, nvars-1+1, nvars-1+1, _state);
    for(i=0; i<=nvars-1; i++)
    {
        for(j=0; j<=nvars-1; j++)
        {
            st.ptr.pp_double[i][j] = 0;
        }
    }
    for(k=0; k<=npoints-1; k++)
    {
        ae_v_move(&tf.ptr.p_double[0], 1, &xy->ptr.pp_double[k][0], 1, ae_v_len(0,nvars-1));
        ae_v_sub(&tf.ptr.p_double[0], 1, &mu.ptr.p_double[0], 1, ae_v_len(0,nvars-1));
        for(i=0; i<=nvars-1; i++)
        {
            v = tf.ptr.p_double[i];
            ae_v_addd(&st.ptr.pp_double[i][0], 1, &tf.ptr.p_double[0], 1, ae_v_len(0,nvars-1), v);
        }
    }
    
    /*
     * Create SW matrix
     */
    ae_matrix_set_length(&sw, nvars-1+1, nvars-1+1, _state);
    for(i=0; i<=nvars-1; i++)
    {
        for(j=0; j<=nvars-1; j++)
        {
            sw.ptr.pp_double[i][j] = 0;
        }
    }
    for(k=0; k<=npoints-1; k++)
    {
        ae_v_move(&tf.ptr.p_double[0], 1, &xy->ptr.pp_double[k][0], 1, ae_v_len(0,nvars-1));
        ae_v_sub(&tf.ptr.p_double[0], 1, &muc.ptr.pp_double[c.ptr.p_int[k]][0], 1, ae_v_len(0,nvars-1));
        for(i=0; i<=nvars-1; i++)
        {
            v = tf.ptr.p_double[i];
            ae_v_addd(&sw.ptr.pp_double[i][0], 1, &tf.ptr.p_double[0], 1, ae_v_len(0,nvars-1), v);
        }
    }
    
    /*
     * Maximize ratio J=(w'*ST*w)/(w'*SW*w).
     *
     * First, make transition from w to v such that w'*ST*w becomes v'*v:
     *    v  = root(ST)*w = R*w
     *    R  = root(D)*Z'
     *    w  = (root(ST)^-1)*v = RI*v
     *    RI = Z*inv(root(D))
     *    J  = (v'*v)/(v'*(RI'*SW*RI)*v)
     *    ST = Z*D*Z'
     *
     *    so we have
     *
     *    J = (v'*v) / (v'*(inv(root(D))*Z'*SW*Z*inv(root(D)))*v)  =
     *      = (v'*v) / (v'*A*v)
     */
    if( !smatrixevd(&st, nvars, 1, ae_true, &d, &z, _state) )
    {
        *info = -4;
        ae_frame_leave(_state);
        return;
    }
    ae_matrix_set_length(w, nvars-1+1, nvars-1+1, _state);
    if( ae_fp_less_eq(d.ptr.p_double[nvars-1],0)||ae_fp_less_eq(d.ptr.p_double[0],1000*ae_machineepsilon*d.ptr.p_double[nvars-1]) )
    {
        
        /*
         * Special case: D[NVars-1]<=0
         * Degenerate task (all variables takes the same value).
         */
        if( ae_fp_less_eq(d.ptr.p_double[nvars-1],0) )
        {
            *info = 2;
            for(i=0; i<=nvars-1; i++)
            {
                for(j=0; j<=nvars-1; j++)
                {
                    if( i==j )
                    {
                        w->ptr.pp_double[i][j] = 1;
                    }
                    else
                    {
                        w->ptr.pp_double[i][j] = 0;
                    }
                }
            }
            ae_frame_leave(_state);
            return;
        }
        
        /*
         * Special case: degenerate ST matrix, multicollinearity found.
         * Since we know ST eigenvalues/vectors we can translate task to
         * non-degenerate form.
         *
         * Let WG is orthogonal basis of the non zero variance subspace
         * of the ST and let WZ is orthogonal basis of the zero variance
         * subspace.
         *
         * Projection on WG allows us to use LDA on reduced M-dimensional
         * subspace, N-M vectors of WZ allows us to update reduced LDA
         * factors to full N-dimensional subspace.
         */
        m = 0;
        for(k=0; k<=nvars-1; k++)
        {
            if( ae_fp_less_eq(d.ptr.p_double[k],1000*ae_machineepsilon*d.ptr.p_double[nvars-1]) )
            {
                m = k+1;
            }
        }
        ae_assert(m!=0, "FisherLDAN: internal error #1", _state);
        ae_matrix_set_length(&xyproj, npoints-1+1, nvars-m+1, _state);
        matrixmatrixmultiply(xy, 0, npoints-1, 0, nvars-1, ae_false, &z, 0, nvars-1, m, nvars-1, ae_false, 1.0, &xyproj, 0, npoints-1, 0, nvars-m-1, 0.0, &work, _state);
        for(i=0; i<=npoints-1; i++)
        {
            xyproj.ptr.pp_double[i][nvars-m] = xy->ptr.pp_double[i][nvars];
        }
        fisherldan(&xyproj, npoints, nvars-m, nclasses, info, &wproj, _state);
        if( *info<0 )
        {
            ae_frame_leave(_state);
            return;
        }
        matrixmatrixmultiply(&z, 0, nvars-1, m, nvars-1, ae_false, &wproj, 0, nvars-m-1, 0, nvars-m-1, ae_false, 1.0, w, 0, nvars-1, 0, nvars-m-1, 0.0, &work, _state);
        for(k=nvars-m; k<=nvars-1; k++)
        {
            ae_v_move(&w->ptr.pp_double[0][k], w->stride, &z.ptr.pp_double[0][k-(nvars-m)], z.stride, ae_v_len(0,nvars-1));
        }
        *info = 2;
    }
    else
    {
        
        /*
         * General case: no multicollinearity
         */
        ae_matrix_set_length(&tm, nvars-1+1, nvars-1+1, _state);
        ae_matrix_set_length(&a, nvars-1+1, nvars-1+1, _state);
        matrixmatrixmultiply(&sw, 0, nvars-1, 0, nvars-1, ae_false, &z, 0, nvars-1, 0, nvars-1, ae_false, 1.0, &tm, 0, nvars-1, 0, nvars-1, 0.0, &work, _state);
        matrixmatrixmultiply(&z, 0, nvars-1, 0, nvars-1, ae_true, &tm, 0, nvars-1, 0, nvars-1, ae_false, 1.0, &a, 0, nvars-1, 0, nvars-1, 0.0, &work, _state);
        for(i=0; i<=nvars-1; i++)
        {
            for(j=0; j<=nvars-1; j++)
            {
                a.ptr.pp_double[i][j] = a.ptr.pp_double[i][j]/ae_sqrt(d.ptr.p_double[i]*d.ptr.p_double[j], _state);
            }
        }
        if( !smatrixevd(&a, nvars, 1, ae_true, &d2, &z2, _state) )
        {
            *info = -4;
            ae_frame_leave(_state);
            return;
        }
        for(k=0; k<=nvars-1; k++)
        {
            for(i=0; i<=nvars-1; i++)
            {
                tf.ptr.p_double[i] = z2.ptr.pp_double[i][k]/ae_sqrt(d.ptr.p_double[i], _state);
            }
            for(i=0; i<=nvars-1; i++)
            {
                v = ae_v_dotproduct(&z.ptr.pp_double[i][0], 1, &tf.ptr.p_double[0], 1, ae_v_len(0,nvars-1));
                w->ptr.pp_double[i][k] = v;
            }
        }
    }
    
    /*
     * Post-processing:
     * * normalization
     * * converting to non-negative form, if possible
     */
    for(k=0; k<=nvars-1; k++)
    {
        v = ae_v_dotproduct(&w->ptr.pp_double[0][k], w->stride, &w->ptr.pp_double[0][k], w->stride, ae_v_len(0,nvars-1));
        v = 1/ae_sqrt(v, _state);
        ae_v_muld(&w->ptr.pp_double[0][k], w->stride, ae_v_len(0,nvars-1), v);
        v = 0;
        for(i=0; i<=nvars-1; i++)
        {
            v = v+w->ptr.pp_double[i][k];
        }
        if( ae_fp_less(v,0) )
        {
            ae_v_muld(&w->ptr.pp_double[0][k], w->stride, ae_v_len(0,nvars-1), -1);
        }
    }
    ae_frame_leave(_state);
}




/*************************************************************************
Creates  neural  network  with  NIn  inputs,  NOut outputs, without hidden
layers, with linear output layer. Network weights are  filled  with  small
random values.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreate0(ae_int_t nin,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    layerscount = 1+3;
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(-5, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_false, network, _state);
    mlpbase_fillhighlevelinformation(network, nin, 0, 0, nout, ae_false, ae_true, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Same  as  MLPCreate0,  but  with  one  hidden  layer  (NHid  neurons) with
non-linear activation function. Output layer is linear.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreate1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    layerscount = 1+3+3;
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(-5, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_false, network, _state);
    mlpbase_fillhighlevelinformation(network, nin, nhid, 0, nout, ae_false, ae_true, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Same as MLPCreate0, but with two hidden layers (NHid1 and  NHid2  neurons)
with non-linear activation function. Output layer is linear.
 $ALL

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreate2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    layerscount = 1+3+3+3;
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid2, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(-5, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_false, network, _state);
    mlpbase_fillhighlevelinformation(network, nin, nhid1, nhid2, nout, ae_false, ae_true, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Creates  neural  network  with  NIn  inputs,  NOut outputs, without hidden
layers with non-linear output layer. Network weights are filled with small
random values.

Activation function of the output layer takes values:

    (B, +INF), if D>=0

or

    (-INF, B), if D<0.


  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreateb0(ae_int_t nin,
     ae_int_t nout,
     double b,
     double d,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    layerscount = 1+3;
    if( ae_fp_greater_eq(d,0) )
    {
        d = 1;
    }
    else
    {
        d = -1;
    }
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(3, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_false, network, _state);
    mlpbase_fillhighlevelinformation(network, nin, 0, 0, nout, ae_false, ae_false, _state);
    
    /*
     * Turn on ouputs shift/scaling.
     */
    for(i=nin; i<=nin+nout-1; i++)
    {
        network->columnmeans.ptr.p_double[i] = b;
        network->columnsigmas.ptr.p_double[i] = d;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Same as MLPCreateB0 but with non-linear hidden layer.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreateb1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     double b,
     double d,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    layerscount = 1+3+3;
    if( ae_fp_greater_eq(d,0) )
    {
        d = 1;
    }
    else
    {
        d = -1;
    }
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(3, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_false, network, _state);
    mlpbase_fillhighlevelinformation(network, nin, nhid, 0, nout, ae_false, ae_false, _state);
    
    /*
     * Turn on ouputs shift/scaling.
     */
    for(i=nin; i<=nin+nout-1; i++)
    {
        network->columnmeans.ptr.p_double[i] = b;
        network->columnsigmas.ptr.p_double[i] = d;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Same as MLPCreateB0 but with two non-linear hidden layers.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreateb2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     double b,
     double d,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    layerscount = 1+3+3+3;
    if( ae_fp_greater_eq(d,0) )
    {
        d = 1;
    }
    else
    {
        d = -1;
    }
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid2, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(3, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_false, network, _state);
    mlpbase_fillhighlevelinformation(network, nin, nhid1, nhid2, nout, ae_false, ae_false, _state);
    
    /*
     * Turn on ouputs shift/scaling.
     */
    for(i=nin; i<=nin+nout-1; i++)
    {
        network->columnmeans.ptr.p_double[i] = b;
        network->columnsigmas.ptr.p_double[i] = d;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Creates  neural  network  with  NIn  inputs,  NOut outputs, without hidden
layers with non-linear output layer. Network weights are filled with small
random values. Activation function of the output layer takes values [A,B].

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreater0(ae_int_t nin,
     ae_int_t nout,
     double a,
     double b,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    layerscount = 1+3;
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_false, network, _state);
    mlpbase_fillhighlevelinformation(network, nin, 0, 0, nout, ae_false, ae_false, _state);
    
    /*
     * Turn on outputs shift/scaling.
     */
    for(i=nin; i<=nin+nout-1; i++)
    {
        network->columnmeans.ptr.p_double[i] = 0.5*(a+b);
        network->columnsigmas.ptr.p_double[i] = 0.5*(a-b);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Same as MLPCreateR0, but with non-linear hidden layer.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreater1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     double a,
     double b,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    layerscount = 1+3+3;
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_false, network, _state);
    mlpbase_fillhighlevelinformation(network, nin, nhid, 0, nout, ae_false, ae_false, _state);
    
    /*
     * Turn on outputs shift/scaling.
     */
    for(i=nin; i<=nin+nout-1; i++)
    {
        network->columnmeans.ptr.p_double[i] = 0.5*(a+b);
        network->columnsigmas.ptr.p_double[i] = 0.5*(a-b);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Same as MLPCreateR0, but with two non-linear hidden layers.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreater2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     double a,
     double b,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;
    ae_int_t i;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    layerscount = 1+3+3+3;
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid2, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_false, network, _state);
    mlpbase_fillhighlevelinformation(network, nin, nhid1, nhid2, nout, ae_false, ae_false, _state);
    
    /*
     * Turn on outputs shift/scaling.
     */
    for(i=nin; i<=nin+nout-1; i++)
    {
        network->columnmeans.ptr.p_double[i] = 0.5*(a+b);
        network->columnsigmas.ptr.p_double[i] = 0.5*(a-b);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Creates classifier network with NIn  inputs  and  NOut  possible  classes.
Network contains no hidden layers and linear output  layer  with  SOFTMAX-
normalization  (so  outputs  sums  up  to  1.0  and  converge to posterior
probabilities).

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreatec0(ae_int_t nin,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    ae_assert(nout>=2, "MLPCreateC0: NOut<2!", _state);
    layerscount = 1+2+1;
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout-1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addzerolayer(&lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_true, network, _state);
    mlpbase_fillhighlevelinformation(network, nin, 0, 0, nout, ae_true, ae_true, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Same as MLPCreateC0, but with one non-linear hidden layer.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreatec1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    ae_assert(nout>=2, "MLPCreateC1: NOut<2!", _state);
    layerscount = 1+3+2+1;
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout-1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addzerolayer(&lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_true, network, _state);
    mlpbase_fillhighlevelinformation(network, nin, nhid, 0, nout, ae_true, ae_true, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Same as MLPCreateC0, but with two non-linear hidden layers.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreatec2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_vector lsizes;
    ae_vector ltypes;
    ae_vector lconnfirst;
    ae_vector lconnlast;
    ae_int_t layerscount;
    ae_int_t lastproc;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&lsizes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&ltypes, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lconnlast, 0, DT_INT, _state, ae_true);

    ae_assert(nout>=2, "MLPCreateC2: NOut<2!", _state);
    layerscount = 1+3+3+2+1;
    
    /*
     * Allocate arrays
     */
    ae_vector_set_length(&lsizes, layerscount-1+1, _state);
    ae_vector_set_length(&ltypes, layerscount-1+1, _state);
    ae_vector_set_length(&lconnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lconnlast, layerscount-1+1, _state);
    
    /*
     * Layers
     */
    mlpbase_addinputlayer(nin, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nhid2, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addactivationlayer(1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addbiasedsummatorlayer(nout-1, &lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    mlpbase_addzerolayer(&lsizes, &ltypes, &lconnfirst, &lconnlast, &lastproc, _state);
    
    /*
     * Create
     */
    mlpbase_mlpcreate(nin, nout, &lsizes, &ltypes, &lconnfirst, &lconnlast, layerscount, ae_true, network, _state);
    mlpbase_fillhighlevelinformation(network, nin, nhid1, nhid2, nout, ae_true, ae_true, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Copying of neural network

INPUT PARAMETERS:
    Network1 -   original

OUTPUT PARAMETERS:
    Network2 -   copy

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcopy(multilayerperceptron* network1,
     multilayerperceptron* network2,
     ae_state *_state)
{

    _multilayerperceptron_clear(network2);

    network2->hlnetworktype = network1->hlnetworktype;
    network2->hlnormtype = network1->hlnormtype;
    copyintegerarray(&network1->hllayersizes, &network2->hllayersizes, _state);
    copyintegerarray(&network1->hlconnections, &network2->hlconnections, _state);
    copyintegerarray(&network1->hlneurons, &network2->hlneurons, _state);
    copyintegerarray(&network1->structinfo, &network2->structinfo, _state);
    copyrealarray(&network1->weights, &network2->weights, _state);
    copyrealarray(&network1->columnmeans, &network2->columnmeans, _state);
    copyrealarray(&network1->columnsigmas, &network2->columnsigmas, _state);
    copyrealarray(&network1->neurons, &network2->neurons, _state);
    copyrealarray(&network1->dfdnet, &network2->dfdnet, _state);
    copyrealarray(&network1->derror, &network2->derror, _state);
    copyrealarray(&network1->x, &network2->x, _state);
    copyrealarray(&network1->y, &network2->y, _state);
    copyrealmatrix(&network1->chunks, &network2->chunks, _state);
    copyrealarray(&network1->nwbuf, &network2->nwbuf, _state);
    copyintegerarray(&network1->integerbuf, &network2->integerbuf, _state);
}


/*************************************************************************
Serialization of MultiLayerPerceptron strucure

INPUT PARAMETERS:
    Network -   original

OUTPUT PARAMETERS:
    RA      -   array of real numbers which stores network,
                array[0..RLen-1]
    RLen    -   RA lenght

  -- ALGLIB --
     Copyright 29.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpserializeold(multilayerperceptron* network,
     /* Real    */ ae_vector* ra,
     ae_int_t* rlen,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t ssize;
    ae_int_t ntotal;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    ae_int_t sigmalen;
    ae_int_t offs;

    ae_vector_clear(ra);
    *rlen = 0;

    
    /*
     * Unload info
     */
    ssize = network->structinfo.ptr.p_int[0];
    nin = network->structinfo.ptr.p_int[1];
    nout = network->structinfo.ptr.p_int[2];
    ntotal = network->structinfo.ptr.p_int[3];
    wcount = network->structinfo.ptr.p_int[4];
    if( mlpissoftmax(network, _state) )
    {
        sigmalen = nin;
    }
    else
    {
        sigmalen = nin+nout;
    }
    
    /*
     *  RA format:
     *      LEN         DESRC.
     *      1           RLen
     *      1           version (MLPVNum)
     *      1           StructInfo size
     *      SSize       StructInfo
     *      WCount      Weights
     *      SigmaLen    ColumnMeans
     *      SigmaLen    ColumnSigmas
     */
    *rlen = 3+ssize+wcount+2*sigmalen;
    ae_vector_set_length(ra, *rlen-1+1, _state);
    ra->ptr.p_double[0] = *rlen;
    ra->ptr.p_double[1] = mlpbase_mlpvnum;
    ra->ptr.p_double[2] = ssize;
    offs = 3;
    for(i=0; i<=ssize-1; i++)
    {
        ra->ptr.p_double[offs+i] = network->structinfo.ptr.p_int[i];
    }
    offs = offs+ssize;
    ae_v_move(&ra->ptr.p_double[offs], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(offs,offs+wcount-1));
    offs = offs+wcount;
    ae_v_move(&ra->ptr.p_double[offs], 1, &network->columnmeans.ptr.p_double[0], 1, ae_v_len(offs,offs+sigmalen-1));
    offs = offs+sigmalen;
    ae_v_move(&ra->ptr.p_double[offs], 1, &network->columnsigmas.ptr.p_double[0], 1, ae_v_len(offs,offs+sigmalen-1));
    offs = offs+sigmalen;
}


/*************************************************************************
Unserialization of MultiLayerPerceptron strucure

INPUT PARAMETERS:
    RA      -   real array which stores network

OUTPUT PARAMETERS:
    Network -   restored network

  -- ALGLIB --
     Copyright 29.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpunserializeold(/* Real    */ ae_vector* ra,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t ssize;
    ae_int_t ntotal;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    ae_int_t sigmalen;
    ae_int_t offs;

    _multilayerperceptron_clear(network);

    ae_assert(ae_round(ra->ptr.p_double[1], _state)==mlpbase_mlpvnum, "MLPUnserialize: incorrect array!", _state);
    
    /*
     * Unload StructInfo from IA
     */
    offs = 3;
    ssize = ae_round(ra->ptr.p_double[2], _state);
    ae_vector_set_length(&network->structinfo, ssize-1+1, _state);
    for(i=0; i<=ssize-1; i++)
    {
        network->structinfo.ptr.p_int[i] = ae_round(ra->ptr.p_double[offs+i], _state);
    }
    offs = offs+ssize;
    
    /*
     * Unload info from StructInfo
     */
    ssize = network->structinfo.ptr.p_int[0];
    nin = network->structinfo.ptr.p_int[1];
    nout = network->structinfo.ptr.p_int[2];
    ntotal = network->structinfo.ptr.p_int[3];
    wcount = network->structinfo.ptr.p_int[4];
    if( network->structinfo.ptr.p_int[6]==0 )
    {
        sigmalen = nin+nout;
    }
    else
    {
        sigmalen = nin;
    }
    
    /*
     * Allocate space for other fields
     */
    ae_vector_set_length(&network->weights, wcount-1+1, _state);
    ae_vector_set_length(&network->columnmeans, sigmalen-1+1, _state);
    ae_vector_set_length(&network->columnsigmas, sigmalen-1+1, _state);
    ae_vector_set_length(&network->neurons, ntotal-1+1, _state);
    ae_matrix_set_length(&network->chunks, 3*ntotal+1, mlpbase_chunksize-1+1, _state);
    ae_vector_set_length(&network->nwbuf, ae_maxint(wcount, 2*nout, _state)-1+1, _state);
    ae_vector_set_length(&network->dfdnet, ntotal-1+1, _state);
    ae_vector_set_length(&network->x, nin-1+1, _state);
    ae_vector_set_length(&network->y, nout-1+1, _state);
    ae_vector_set_length(&network->derror, ntotal-1+1, _state);
    
    /*
     * Copy parameters from RA
     */
    ae_v_move(&network->weights.ptr.p_double[0], 1, &ra->ptr.p_double[offs], 1, ae_v_len(0,wcount-1));
    offs = offs+wcount;
    ae_v_move(&network->columnmeans.ptr.p_double[0], 1, &ra->ptr.p_double[offs], 1, ae_v_len(0,sigmalen-1));
    offs = offs+sigmalen;
    ae_v_move(&network->columnsigmas.ptr.p_double[0], 1, &ra->ptr.p_double[offs], 1, ae_v_len(0,sigmalen-1));
    offs = offs+sigmalen;
}


/*************************************************************************
Randomization of neural network weights

  -- ALGLIB --
     Copyright 06.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlprandomize(multilayerperceptron* network, ae_state *_state)
{
    ae_int_t i;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;


    mlpproperties(network, &nin, &nout, &wcount, _state);
    for(i=0; i<=wcount-1; i++)
    {
        network->weights.ptr.p_double[i] = ae_randomreal(_state)-0.5;
    }
}


/*************************************************************************
Randomization of neural network weights and standartisator

  -- ALGLIB --
     Copyright 10.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlprandomizefull(multilayerperceptron* network, ae_state *_state)
{
    ae_int_t i;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    ae_int_t ntotal;
    ae_int_t istart;
    ae_int_t offs;
    ae_int_t ntype;


    mlpproperties(network, &nin, &nout, &wcount, _state);
    ntotal = network->structinfo.ptr.p_int[3];
    istart = network->structinfo.ptr.p_int[5];
    
    /*
     * Process network
     */
    for(i=0; i<=wcount-1; i++)
    {
        network->weights.ptr.p_double[i] = ae_randomreal(_state)-0.5;
    }
    for(i=0; i<=nin-1; i++)
    {
        network->columnmeans.ptr.p_double[i] = 2*ae_randomreal(_state)-1;
        network->columnsigmas.ptr.p_double[i] = 1.5*ae_randomreal(_state)+0.5;
    }
    if( !mlpissoftmax(network, _state) )
    {
        for(i=0; i<=nout-1; i++)
        {
            offs = istart+(ntotal-nout+i)*mlpbase_nfieldwidth;
            ntype = network->structinfo.ptr.p_int[offs+0];
            if( ntype==0 )
            {
                
                /*
                 * Shifts are changed only for linear outputs neurons
                 */
                network->columnmeans.ptr.p_double[nin+i] = 2*ae_randomreal(_state)-1;
            }
            if( ntype==0||ntype==3 )
            {
                
                /*
                 * Scales are changed only for linear or bounded outputs neurons.
                 * Note that scale randomization preserves sign.
                 */
                network->columnsigmas.ptr.p_double[nin+i] = ae_sign(network->columnsigmas.ptr.p_double[nin+i], _state)*(1.5*ae_randomreal(_state)+0.5);
            }
        }
    }
}


/*************************************************************************
Internal subroutine.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpinitpreprocessor(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_int_t jmax;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    ae_int_t ntotal;
    ae_int_t istart;
    ae_int_t offs;
    ae_int_t ntype;
    ae_vector means;
    ae_vector sigmas;
    double s;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&means, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&sigmas, 0, DT_REAL, _state, ae_true);

    mlpproperties(network, &nin, &nout, &wcount, _state);
    ntotal = network->structinfo.ptr.p_int[3];
    istart = network->structinfo.ptr.p_int[5];
    
    /*
     * Means/Sigmas
     */
    if( mlpissoftmax(network, _state) )
    {
        jmax = nin-1;
    }
    else
    {
        jmax = nin+nout-1;
    }
    ae_vector_set_length(&means, jmax+1, _state);
    ae_vector_set_length(&sigmas, jmax+1, _state);
    for(j=0; j<=jmax; j++)
    {
        means.ptr.p_double[j] = 0;
        for(i=0; i<=ssize-1; i++)
        {
            means.ptr.p_double[j] = means.ptr.p_double[j]+xy->ptr.pp_double[i][j];
        }
        means.ptr.p_double[j] = means.ptr.p_double[j]/ssize;
        sigmas.ptr.p_double[j] = 0;
        for(i=0; i<=ssize-1; i++)
        {
            sigmas.ptr.p_double[j] = sigmas.ptr.p_double[j]+ae_sqr(xy->ptr.pp_double[i][j]-means.ptr.p_double[j], _state);
        }
        sigmas.ptr.p_double[j] = ae_sqrt(sigmas.ptr.p_double[j]/ssize, _state);
    }
    
    /*
     * Inputs
     */
    for(i=0; i<=nin-1; i++)
    {
        network->columnmeans.ptr.p_double[i] = means.ptr.p_double[i];
        network->columnsigmas.ptr.p_double[i] = sigmas.ptr.p_double[i];
        if( ae_fp_eq(network->columnsigmas.ptr.p_double[i],0) )
        {
            network->columnsigmas.ptr.p_double[i] = 1;
        }
    }
    
    /*
     * Outputs
     */
    if( !mlpissoftmax(network, _state) )
    {
        for(i=0; i<=nout-1; i++)
        {
            offs = istart+(ntotal-nout+i)*mlpbase_nfieldwidth;
            ntype = network->structinfo.ptr.p_int[offs+0];
            
            /*
             * Linear outputs
             */
            if( ntype==0 )
            {
                network->columnmeans.ptr.p_double[nin+i] = means.ptr.p_double[nin+i];
                network->columnsigmas.ptr.p_double[nin+i] = sigmas.ptr.p_double[nin+i];
                if( ae_fp_eq(network->columnsigmas.ptr.p_double[nin+i],0) )
                {
                    network->columnsigmas.ptr.p_double[nin+i] = 1;
                }
            }
            
            /*
             * Bounded outputs (half-interval)
             */
            if( ntype==3 )
            {
                s = means.ptr.p_double[nin+i]-network->columnmeans.ptr.p_double[nin+i];
                if( ae_fp_eq(s,0) )
                {
                    s = ae_sign(network->columnsigmas.ptr.p_double[nin+i], _state);
                }
                if( ae_fp_eq(s,0) )
                {
                    s = 1.0;
                }
                network->columnsigmas.ptr.p_double[nin+i] = ae_sign(network->columnsigmas.ptr.p_double[nin+i], _state)*ae_fabs(s, _state);
                if( ae_fp_eq(network->columnsigmas.ptr.p_double[nin+i],0) )
                {
                    network->columnsigmas.ptr.p_double[nin+i] = 1;
                }
            }
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Returns information about initialized network: number of inputs, outputs,
weights.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpproperties(multilayerperceptron* network,
     ae_int_t* nin,
     ae_int_t* nout,
     ae_int_t* wcount,
     ae_state *) //_state)
{

    *nin = 0;
    *nout = 0;
    *wcount = 0;

    *nin = network->structinfo.ptr.p_int[1];
    *nout = network->structinfo.ptr.p_int[2];
    *wcount = network->structinfo.ptr.p_int[4];
}


/*************************************************************************
Returns number of inputs.

  -- ALGLIB --
     Copyright 19.10.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpgetinputscount(multilayerperceptron* network,
     ae_state *) //_state)
{
    ae_int_t result;


    result = network->structinfo.ptr.p_int[1];
    return result;
}


/*************************************************************************
Returns number of outputs.

  -- ALGLIB --
     Copyright 19.10.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpgetoutputscount(multilayerperceptron* network,
     ae_state *) //_state)
{
    ae_int_t result;


    result = network->structinfo.ptr.p_int[2];
    return result;
}


/*************************************************************************
Returns number of weights.

  -- ALGLIB --
     Copyright 19.10.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpgetweightscount(multilayerperceptron* network,
     ae_state *) //_state)
{
    ae_int_t result;


    result = network->structinfo.ptr.p_int[4];
    return result;
}


/*************************************************************************
Tells whether network is SOFTMAX-normalized (i.e. classifier) or not.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
ae_bool mlpissoftmax(multilayerperceptron* network, ae_state *) //_state)
{
    ae_bool result;


    result = network->structinfo.ptr.p_int[6]==1;
    return result;
}


/*************************************************************************
This function returns total number of layers (including input, hidden and
output layers).

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpgetlayerscount(multilayerperceptron* network,
     ae_state *) //_state)
{
    ae_int_t result;


    result = network->hllayersizes.cnt;
    return result;
}


/*************************************************************************
This function returns size of K-th layer.

K=0 corresponds to input layer, K=CNT-1 corresponds to output layer.

Size of the output layer is always equal to the number of outputs, although
when we have softmax-normalized network, last neuron doesn't have any
connections - it is just zero.

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpgetlayersize(multilayerperceptron* network,
     ae_int_t k,
     ae_state *_state)
{
    ae_int_t result;


    ae_assert(k>=0&&k<network->hllayersizes.cnt, "MLPGetLayerSize: incorrect layer index", _state);
    result = network->hllayersizes.ptr.p_int[k];
    return result;
}


/*************************************************************************
This function returns offset/scaling coefficients for I-th input of the
network.

INPUT PARAMETERS:
    Network     -   network
    I           -   input index

OUTPUT PARAMETERS:
    Mean        -   mean term
    Sigma       -   sigma term, guaranteed to be nonzero.

I-th input is passed through linear transformation
    IN[i] = (IN[i]-Mean)/Sigma
before feeding to the network

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpgetinputscaling(multilayerperceptron* network,
     ae_int_t i,
     double* mean,
     double* sigma,
     ae_state *_state)
{

    *mean = 0;
    *sigma = 0;

    ae_assert(i>=0&&i<network->hllayersizes.ptr.p_int[0], "MLPGetInputScaling: incorrect (nonexistent) I", _state);
    *mean = network->columnmeans.ptr.p_double[i];
    *sigma = network->columnsigmas.ptr.p_double[i];
    if( ae_fp_eq(*sigma,0) )
    {
        *sigma = 1;
    }
}


/*************************************************************************
This function returns offset/scaling coefficients for I-th output of the
network.

INPUT PARAMETERS:
    Network     -   network
    I           -   input index

OUTPUT PARAMETERS:
    Mean        -   mean term
    Sigma       -   sigma term, guaranteed to be nonzero.

I-th output is passed through linear transformation
    OUT[i] = OUT[i]*Sigma+Mean
before returning it to user. In case we have SOFTMAX-normalized network,
we return (Mean,Sigma)=(0.0,1.0).

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpgetoutputscaling(multilayerperceptron* network,
     ae_int_t i,
     double* mean,
     double* sigma,
     ae_state* _state)
{

    *mean = 0;
    *sigma = 0;

    ae_assert(i>=0&&i<network->hllayersizes.ptr.p_int[network->hllayersizes.cnt-1], "MLPGetOutputScaling: incorrect (nonexistent) I", _state);
    if( network->structinfo.ptr.p_int[6]==1 )
    {
        *mean = 0;
        *sigma = 1;
    }
    else
    {
        *mean = network->columnmeans.ptr.p_double[network->hllayersizes.ptr.p_int[0]+i];
        *sigma = network->columnsigmas.ptr.p_double[network->hllayersizes.ptr.p_int[0]+i];
    }
}


/*************************************************************************
This function returns information about Ith neuron of Kth layer

INPUT PARAMETERS:
    Network     -   network
    K           -   layer index
    I           -   neuron index (within layer)

OUTPUT PARAMETERS:
    FKind       -   activation function type (used by MLPActivationFunction())
                    this value is zero for input or linear neurons
    Threshold   -   also called offset, bias
                    zero for input neurons
                    
NOTE: this function throws exception if layer or neuron with  given  index
do not exists.

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpgetneuroninfo(multilayerperceptron* network,
     ae_int_t k,
     ae_int_t i,
     ae_int_t* fkind,
     double* threshold,
     ae_state *_state)
{
    ae_int_t ncnt;
    ae_int_t istart;
    ae_int_t highlevelidx;
    ae_int_t activationoffset;

    *fkind = 0;
    *threshold = 0;

    ncnt = network->hlneurons.cnt/mlpbase_hlnfieldwidth;
    istart = network->structinfo.ptr.p_int[5];
    
    /*
     * search
     */
    network->integerbuf.ptr.p_int[0] = k;
    network->integerbuf.ptr.p_int[1] = i;
    highlevelidx = recsearch(&network->hlneurons, mlpbase_hlnfieldwidth, 2, 0, ncnt, &network->integerbuf, _state);
    ae_assert(highlevelidx>=0, "MLPGetNeuronInfo: incorrect (nonexistent) layer or neuron index", _state);
    
    /*
     * 1. find offset of the activation function record in the
     */
    if( network->hlneurons.ptr.p_int[highlevelidx*mlpbase_hlnfieldwidth+2]>=0 )
    {
        activationoffset = istart+network->hlneurons.ptr.p_int[highlevelidx*mlpbase_hlnfieldwidth+2]*mlpbase_nfieldwidth;
        *fkind = network->structinfo.ptr.p_int[activationoffset+0];
    }
    else
    {
        *fkind = 0;
    }
    if( network->hlneurons.ptr.p_int[highlevelidx*mlpbase_hlnfieldwidth+3]>=0 )
    {
        *threshold = network->weights.ptr.p_double[network->hlneurons.ptr.p_int[highlevelidx*mlpbase_hlnfieldwidth+3]];
    }
    else
    {
        *threshold = 0;
    }
}


/*************************************************************************
This function returns information about connection from I0-th neuron of
K0-th layer to I1-th neuron of K1-th layer.

INPUT PARAMETERS:
    Network     -   network
    K0          -   layer index
    I0          -   neuron index (within layer)
    K1          -   layer index
    I1          -   neuron index (within layer)

RESULT:
    connection weight (zero for non-existent connections)

This function:
1. throws exception if layer or neuron with given index do not exists.
2. returns zero if neurons exist, but there is no connection between them

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
double mlpgetweight(multilayerperceptron* network,
     ae_int_t k0,
     ae_int_t i0,
     ae_int_t k1,
     ae_int_t i1,
     ae_state *_state)
{
    ae_int_t ccnt;
    ae_int_t highlevelidx;
    double result;


    ccnt = network->hlconnections.cnt/mlpbase_hlconnfieldwidth;
    
    /*
     * check params
     */
    ae_assert(k0>=0&&k0<network->hllayersizes.cnt, "MLPGetWeight: incorrect (nonexistent) K0", _state);
    ae_assert(i0>=0&&i0<network->hllayersizes.ptr.p_int[k0], "MLPGetWeight: incorrect (nonexistent) I0", _state);
    ae_assert(k1>=0&&k1<network->hllayersizes.cnt, "MLPGetWeight: incorrect (nonexistent) K1", _state);
    ae_assert(i1>=0&&i1<network->hllayersizes.ptr.p_int[k1], "MLPGetWeight: incorrect (nonexistent) I1", _state);
    
    /*
     * search
     */
    network->integerbuf.ptr.p_int[0] = k0;
    network->integerbuf.ptr.p_int[1] = i0;
    network->integerbuf.ptr.p_int[2] = k1;
    network->integerbuf.ptr.p_int[3] = i1;
    highlevelidx = recsearch(&network->hlconnections, mlpbase_hlconnfieldwidth, 4, 0, ccnt, &network->integerbuf, _state);
    if( highlevelidx>=0 )
    {
        result = network->weights.ptr.p_double[network->hlconnections.ptr.p_int[highlevelidx*mlpbase_hlconnfieldwidth+4]];
    }
    else
    {
        result = 0;
    }
    return result;
}


/*************************************************************************
This function sets offset/scaling coefficients for I-th input of the
network.

INPUT PARAMETERS:
    Network     -   network
    I           -   input index
    Mean        -   mean term
    Sigma       -   sigma term (if zero, will be replaced by 1.0)

NTE: I-th input is passed through linear transformation
    IN[i] = (IN[i]-Mean)/Sigma
before feeding to the network. This function sets Mean and Sigma.

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpsetinputscaling(multilayerperceptron* network,
     ae_int_t i,
     double mean,
     double sigma,
     ae_state *_state)
{


    ae_assert(i>=0&&i<network->hllayersizes.ptr.p_int[0], "MLPSetInputScaling: incorrect (nonexistent) I", _state);
    ae_assert(ae_isfinite(mean, _state), "MLPSetInputScaling: infinite or NAN Mean", _state);
    ae_assert(ae_isfinite(sigma, _state), "MLPSetInputScaling: infinite or NAN Sigma", _state);
    if( ae_fp_eq(sigma,0) )
    {
        sigma = 1;
    }
    network->columnmeans.ptr.p_double[i] = mean;
    network->columnsigmas.ptr.p_double[i] = sigma;
}


/*************************************************************************
This function sets offset/scaling coefficients for I-th output of the
network.

INPUT PARAMETERS:
    Network     -   network
    I           -   input index
    Mean        -   mean term
    Sigma       -   sigma term (if zero, will be replaced by 1.0)

OUTPUT PARAMETERS:

NOTE: I-th output is passed through linear transformation
    OUT[i] = OUT[i]*Sigma+Mean
before returning it to user. This function sets Sigma/Mean. In case we
have SOFTMAX-normalized network, you can not set (Sigma,Mean) to anything
other than(0.0,1.0) - this function will throw exception.

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpsetoutputscaling(multilayerperceptron* network,
     ae_int_t i,
     double mean,
     double sigma,
     ae_state *_state)
{


    ae_assert(i>=0&&i<network->hllayersizes.ptr.p_int[network->hllayersizes.cnt-1], "MLPSetOutputScaling: incorrect (nonexistent) I", _state);
    ae_assert(ae_isfinite(mean, _state), "MLPSetOutputScaling: infinite or NAN Mean", _state);
    ae_assert(ae_isfinite(sigma, _state), "MLPSetOutputScaling: infinite or NAN Sigma", _state);
    if( network->structinfo.ptr.p_int[6]==1 )
    {
        ae_assert(ae_fp_eq(mean,0), "MLPSetOutputScaling: you can not set non-zero Mean term for classifier network", _state);
        ae_assert(ae_fp_eq(sigma,1), "MLPSetOutputScaling: you can not set non-unit Sigma term for classifier network", _state);
    }
    else
    {
        if( ae_fp_eq(sigma,0) )
        {
            sigma = 1;
        }
        network->columnmeans.ptr.p_double[network->hllayersizes.ptr.p_int[0]+i] = mean;
        network->columnsigmas.ptr.p_double[network->hllayersizes.ptr.p_int[0]+i] = sigma;
    }
}


/*************************************************************************
This function modifies information about Ith neuron of Kth layer

INPUT PARAMETERS:
    Network     -   network
    K           -   layer index
    I           -   neuron index (within layer)
    FKind       -   activation function type (used by MLPActivationFunction())
                    this value must be zero for input neurons
                    (you can not set activation function for input neurons)
    Threshold   -   also called offset, bias
                    this value must be zero for input neurons
                    (you can not set threshold for input neurons)

NOTES:
1. this function throws exception if layer or neuron with given index do
   not exists.
2. this function also throws exception when you try to set non-linear
   activation function for input neurons (any kind of network) or for output
   neurons of classifier network.
3. this function throws exception when you try to set non-zero threshold for
   input neurons (any kind of network).

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpsetneuroninfo(multilayerperceptron* network,
     ae_int_t k,
     ae_int_t i,
     ae_int_t fkind,
     double threshold,
     ae_state *_state)
{
    ae_int_t ncnt;
    ae_int_t istart;
    ae_int_t highlevelidx;
    ae_int_t activationoffset;


    ae_assert(ae_isfinite(threshold, _state), "MLPSetNeuronInfo: infinite or NAN Threshold", _state);
    
    /*
     * convenience vars
     */
    ncnt = network->hlneurons.cnt/mlpbase_hlnfieldwidth;
    istart = network->structinfo.ptr.p_int[5];
    
    /*
     * search
     */
    network->integerbuf.ptr.p_int[0] = k;
    network->integerbuf.ptr.p_int[1] = i;
    highlevelidx = recsearch(&network->hlneurons, mlpbase_hlnfieldwidth, 2, 0, ncnt, &network->integerbuf, _state);
    ae_assert(highlevelidx>=0, "MLPSetNeuronInfo: incorrect (nonexistent) layer or neuron index", _state);
    
    /*
     * activation function
     */
    if( network->hlneurons.ptr.p_int[highlevelidx*mlpbase_hlnfieldwidth+2]>=0 )
    {
        activationoffset = istart+network->hlneurons.ptr.p_int[highlevelidx*mlpbase_hlnfieldwidth+2]*mlpbase_nfieldwidth;
        network->structinfo.ptr.p_int[activationoffset+0] = fkind;
    }
    else
    {
        ae_assert(fkind==0, "MLPSetNeuronInfo: you try to set activation function for neuron which can not have one", _state);
    }
    
    /*
     * Threshold
     */
    if( network->hlneurons.ptr.p_int[highlevelidx*mlpbase_hlnfieldwidth+3]>=0 )
    {
        network->weights.ptr.p_double[network->hlneurons.ptr.p_int[highlevelidx*mlpbase_hlnfieldwidth+3]] = threshold;
    }
    else
    {
        ae_assert(ae_fp_eq(threshold,0), "MLPSetNeuronInfo: you try to set non-zero threshold for neuron which can not have one", _state);
    }
}


/*************************************************************************
This function modifies information about connection from I0-th neuron of
K0-th layer to I1-th neuron of K1-th layer.

INPUT PARAMETERS:
    Network     -   network
    K0          -   layer index
    I0          -   neuron index (within layer)
    K1          -   layer index
    I1          -   neuron index (within layer)
    W           -   connection weight (must be zero for non-existent
                    connections)

This function:
1. throws exception if layer or neuron with given index do not exists.
2. throws exception if you try to set non-zero weight for non-existent
   connection

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpsetweight(multilayerperceptron* network,
     ae_int_t k0,
     ae_int_t i0,
     ae_int_t k1,
     ae_int_t i1,
     double w,
     ae_state *_state)
{
    ae_int_t ccnt;
    ae_int_t highlevelidx;


    ccnt = network->hlconnections.cnt/mlpbase_hlconnfieldwidth;
    
    /*
     * check params
     */
    ae_assert(k0>=0&&k0<network->hllayersizes.cnt, "MLPSetWeight: incorrect (nonexistent) K0", _state);
    ae_assert(i0>=0&&i0<network->hllayersizes.ptr.p_int[k0], "MLPSetWeight: incorrect (nonexistent) I0", _state);
    ae_assert(k1>=0&&k1<network->hllayersizes.cnt, "MLPSetWeight: incorrect (nonexistent) K1", _state);
    ae_assert(i1>=0&&i1<network->hllayersizes.ptr.p_int[k1], "MLPSetWeight: incorrect (nonexistent) I1", _state);
    ae_assert(ae_isfinite(w, _state), "MLPSetWeight: infinite or NAN weight", _state);
    
    /*
     * search
     */
    network->integerbuf.ptr.p_int[0] = k0;
    network->integerbuf.ptr.p_int[1] = i0;
    network->integerbuf.ptr.p_int[2] = k1;
    network->integerbuf.ptr.p_int[3] = i1;
    highlevelidx = recsearch(&network->hlconnections, mlpbase_hlconnfieldwidth, 4, 0, ccnt, &network->integerbuf, _state);
    if( highlevelidx>=0 )
    {
        network->weights.ptr.p_double[network->hlconnections.ptr.p_int[highlevelidx*mlpbase_hlconnfieldwidth+4]] = w;
    }
    else
    {
        ae_assert(ae_fp_eq(w,0), "MLPSetWeight: you try to set non-zero weight for non-existent connection", _state);
    }
}


/*************************************************************************
Neural network activation function

INPUT PARAMETERS:
    NET         -   neuron input
    K           -   function index (zero for linear function)

OUTPUT PARAMETERS:
    F           -   function
    DF          -   its derivative
    D2F         -   its second derivative

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpactivationfunction(double net,
     ae_int_t k,
     double* f,
     double* df,
     double* d2f,
     ae_state *_state)
{
    double net2;
    double arg;
    double root;
    double r;

    *f = 0;
    *df = 0;
    *d2f = 0;

    if( k==0||k==-5 )
    {
        *f = net;
        *df = 1;
        *d2f = 0;
        return;
    }
    if( k==1 )
    {
        
        /*
         * TanH activation function
         */
        if( ae_fp_less(ae_fabs(net, _state),100) )
        {
            *f = ae_tanh(net, _state);
        }
        else
        {
            *f = ae_sign(net, _state);
        }
        *df = 1-ae_sqr(*f, _state);
        *d2f = -2*(*f)*(*df);
        return;
    }
    if( k==3 )
    {
        
        /*
         * EX activation function
         */
        if( ae_fp_greater_eq(net,0) )
        {
            net2 = net*net;
            arg = net2+1;
            root = ae_sqrt(arg, _state);
            *f = net+root;
            r = net/root;
            *df = 1+r;
            *d2f = (root-net*r)/arg;
        }
        else
        {
            *f = ae_exp(net, _state);
            *df = *f;
            *d2f = *f;
        }
        return;
    }
    if( k==2 )
    {
        *f = ae_exp(-ae_sqr(net, _state), _state);
        *df = -2*net*(*f);
        *d2f = -2*(*f+*df*net);
        return;
    }
    *f = 0;
    *df = 0;
    *d2f = 0;
}


/*************************************************************************
Procesing

INPUT PARAMETERS:
    Network -   neural network
    X       -   input vector,  array[0..NIn-1].

OUTPUT PARAMETERS:
    Y       -   result. Regression estimate when solving regression  task,
                vector of posterior probabilities for classification task.

See also MLPProcessI

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpprocess(multilayerperceptron* network,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{


    if( y->cnt<network->structinfo.ptr.p_int[2] )
    {
        ae_vector_set_length(y, network->structinfo.ptr.p_int[2], _state);
    }
    mlpinternalprocessvector(&network->structinfo, &network->weights, &network->columnmeans, &network->columnsigmas, &network->neurons, &network->dfdnet, x, y, _state);
}


/*************************************************************************
'interactive'  variant  of  MLPProcess  for  languages  like  Python which
support constructs like "Y = MLPProcess(NN,X)" and interactive mode of the
interpreter

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 21.09.2010 by Bochkanov Sergey
*************************************************************************/
void mlpprocessi(multilayerperceptron* network,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{

    ae_vector_clear(y);

    mlpprocess(network, x, y, _state);
}


/*************************************************************************
Error function for neural network, internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
double mlperror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t k;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    double e;
    double result;


    mlpproperties(network, &nin, &nout, &wcount, _state);
    result = 0;
    for(i=0; i<=ssize-1; i++)
    {
        ae_v_move(&network->x.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nin-1));
        mlpprocess(network, &network->x, &network->y, _state);
        if( mlpissoftmax(network, _state) )
        {
            
            /*
             * class labels outputs
             */
            k = ae_round(xy->ptr.pp_double[i][nin], _state);
            if( k>=0&&k<nout )
            {
                network->y.ptr.p_double[k] = network->y.ptr.p_double[k]-1;
            }
        }
        else
        {
            
            /*
             * real outputs
             */
            ae_v_sub(&network->y.ptr.p_double[0], 1, &xy->ptr.pp_double[i][nin], 1, ae_v_len(0,nout-1));
        }
        e = ae_v_dotproduct(&network->y.ptr.p_double[0], 1, &network->y.ptr.p_double[0], 1, ae_v_len(0,nout-1));
        result = result+e/2;
    }
    return result;
}


/*************************************************************************
Natural error function for neural network, internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
double mlperrorn(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t k;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    double e;
    double result;


    mlpproperties(network, &nin, &nout, &wcount, _state);
    result = 0;
    for(i=0; i<=ssize-1; i++)
    {
        
        /*
         * Process vector
         */
        ae_v_move(&network->x.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nin-1));
        mlpprocess(network, &network->x, &network->y, _state);
        
        /*
         * Update error function
         */
        if( network->structinfo.ptr.p_int[6]==0 )
        {
            
            /*
             * Least squares error function
             */
            ae_v_sub(&network->y.ptr.p_double[0], 1, &xy->ptr.pp_double[i][nin], 1, ae_v_len(0,nout-1));
            e = ae_v_dotproduct(&network->y.ptr.p_double[0], 1, &network->y.ptr.p_double[0], 1, ae_v_len(0,nout-1));
            result = result+e/2;
        }
        else
        {
            
            /*
             * Cross-entropy error function
             */
            k = ae_round(xy->ptr.pp_double[i][nin], _state);
            if( k>=0&&k<nout )
            {
                result = result+mlpbase_safecrossentropy(1, network->y.ptr.p_double[k], _state);
            }
        }
    }
    return result;
}


/*************************************************************************
Classification error

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpclserror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    ae_vector workx;
    ae_vector worky;
    ae_int_t nn;
    ae_int_t ns;
    ae_int_t nmax;
    ae_int_t result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&workx, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&worky, 0, DT_REAL, _state, ae_true);

    mlpproperties(network, &nin, &nout, &wcount, _state);
    ae_vector_set_length(&workx, nin-1+1, _state);
    ae_vector_set_length(&worky, nout-1+1, _state);
    result = 0;
    for(i=0; i<=ssize-1; i++)
    {
        
        /*
         * Process
         */
        ae_v_move(&workx.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nin-1));
        mlpprocess(network, &workx, &worky, _state);
        
        /*
         * Network version of the answer
         */
        nmax = 0;
        for(j=0; j<=nout-1; j++)
        {
            if( ae_fp_greater(worky.ptr.p_double[j],worky.ptr.p_double[nmax]) )
            {
                nmax = j;
            }
        }
        nn = nmax;
        
        /*
         * Right answer
         */
        if( mlpissoftmax(network, _state) )
        {
            ns = ae_round(xy->ptr.pp_double[i][nin], _state);
        }
        else
        {
            nmax = 0;
            for(j=0; j<=nout-1; j++)
            {
                if( ae_fp_greater(xy->ptr.pp_double[i][nin+j],xy->ptr.pp_double[i][nin+nmax]) )
                {
                    nmax = j;
                }
            }
            ns = nmax;
        }
        
        /*
         * compare
         */
        if( nn!=ns )
        {
            result = result+1;
        }
    }
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
Relative classification error on the test set

INPUT PARAMETERS:
    Network -   network
    XY      -   test set
    NPoints -   test set size

RESULT:
    percent of incorrectly classified cases. Works both for
    classifier networks and general purpose networks used as
    classifiers.

  -- ALGLIB --
     Copyright 25.12.2008 by Bochkanov Sergey
*************************************************************************/
double mlprelclserror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    double result;


    result = (double)mlpclserror(network, xy, npoints, _state)/(double)npoints;
    return result;
}


/*************************************************************************
Average cross-entropy (in bits per element) on the test set

INPUT PARAMETERS:
    Network -   neural network
    XY      -   test set
    NPoints -   test set size

RESULT:
    CrossEntropy/(NPoints*LN(2)).
    Zero if network solves regression task.

  -- ALGLIB --
     Copyright 08.01.2009 by Bochkanov Sergey
*************************************************************************/
double mlpavgce(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    double result;


    if( mlpissoftmax(network, _state) )
    {
        mlpproperties(network, &nin, &nout, &wcount, _state);
        result = mlperrorn(network, xy, npoints, _state)/(npoints*ae_log(2, _state));
    }
    else
    {
        result = 0;
    }
    return result;
}


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    Network -   neural network
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.
    Its meaning for regression task is obvious. As for
    classification task, RMS error means error when estimating posterior
    probabilities.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
double mlprmserror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    double result;


    mlpproperties(network, &nin, &nout, &wcount, _state);
    result = ae_sqrt(2*mlperror(network, xy, npoints, _state)/(npoints*nout), _state);
    return result;
}


/*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    Network -   neural network
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for
    classification task, it means average error when estimating posterior
    probabilities.

  -- ALGLIB --
     Copyright 11.03.2008 by Bochkanov Sergey
*************************************************************************/
double mlpavgerror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    double result;


    mlpproperties(network, &nin, &nout, &wcount, _state);
    result = 0;
    for(i=0; i<=npoints-1; i++)
    {
        ae_v_move(&network->x.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nin-1));
        mlpprocess(network, &network->x, &network->y, _state);
        if( mlpissoftmax(network, _state) )
        {
            
            /*
             * class labels
             */
            k = ae_round(xy->ptr.pp_double[i][nin], _state);
            for(j=0; j<=nout-1; j++)
            {
                if( j==k )
                {
                    result = result+ae_fabs(1-network->y.ptr.p_double[j], _state);
                }
                else
                {
                    result = result+ae_fabs(network->y.ptr.p_double[j], _state);
                }
            }
        }
        else
        {
            
            /*
             * real outputs
             */
            for(j=0; j<=nout-1; j++)
            {
                result = result+ae_fabs(xy->ptr.pp_double[i][nin+j]-network->y.ptr.p_double[j], _state);
            }
        }
    }
    result = result/(npoints*nout);
    return result;
}


/*************************************************************************
Average relative error on the test set

INPUT PARAMETERS:
    Network -   neural network
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for
    classification task, it means average relative error when estimating
    posterior probability of belonging to the correct class.

  -- ALGLIB --
     Copyright 11.03.2008 by Bochkanov Sergey
*************************************************************************/
double mlpavgrelerror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t lk;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    double result;


    mlpproperties(network, &nin, &nout, &wcount, _state);
    result = 0;
    k = 0;
    for(i=0; i<=npoints-1; i++)
    {
        ae_v_move(&network->x.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nin-1));
        mlpprocess(network, &network->x, &network->y, _state);
        if( mlpissoftmax(network, _state) )
        {
            
            /*
             * class labels
             */
            lk = ae_round(xy->ptr.pp_double[i][nin], _state);
            for(j=0; j<=nout-1; j++)
            {
                if( j==lk )
                {
                    result = result+ae_fabs(1-network->y.ptr.p_double[j], _state);
                    k = k+1;
                }
            }
        }
        else
        {
            
            /*
             * real outputs
             */
            for(j=0; j<=nout-1; j++)
            {
                if( ae_fp_neq(xy->ptr.pp_double[i][nin+j],0) )
                {
                    result = result+ae_fabs(xy->ptr.pp_double[i][nin+j]-network->y.ptr.p_double[j], _state)/ae_fabs(xy->ptr.pp_double[i][nin+j], _state);
                    k = k+1;
                }
            }
        }
    }
    if( k!=0 )
    {
        result = result/k;
    }
    return result;
}


/*************************************************************************
Gradient calculation

INPUT PARAMETERS:
    Network -   network initialized with one of the network creation funcs
    X       -   input vector, length of array must be at least NIn
    DesiredY-   desired outputs, length of array must be at least NOut
    Grad    -   possibly preallocated array. If size of array is smaller
                than WCount, it will be reallocated. It is recommended to
                reuse previously allocated array to reduce allocation
                overhead.

OUTPUT PARAMETERS:
    E       -   error function, SUM(sqr(y[i]-desiredy[i])/2,i)
    Grad    -   gradient of E with respect to weights of network, array[WCount]
    
  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgrad(multilayerperceptron* network,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* desiredy,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t nout;
    ae_int_t ntotal;

    *e = 0;

    
    /*
     * Alloc
     */
    if( grad->cnt<network->structinfo.ptr.p_int[4] )
    {
        ae_vector_set_length(grad, network->structinfo.ptr.p_int[4], _state);
    }
    
    /*
     * Prepare dError/dOut, internal structures
     */
    mlpprocess(network, x, &network->y, _state);
    nout = network->structinfo.ptr.p_int[2];
    ntotal = network->structinfo.ptr.p_int[3];
    *e = 0;
    for(i=0; i<=ntotal-1; i++)
    {
        network->derror.ptr.p_double[i] = 0;
    }
    for(i=0; i<=nout-1; i++)
    {
        network->derror.ptr.p_double[ntotal-nout+i] = network->y.ptr.p_double[i]-desiredy->ptr.p_double[i];
        *e = *e+ae_sqr(network->y.ptr.p_double[i]-desiredy->ptr.p_double[i], _state)/2;
    }
    
    /*
     * gradient
     */
    mlpbase_mlpinternalcalculategradient(network, &network->neurons, &network->weights, &network->derror, grad, ae_false, _state);
}


/*************************************************************************
Gradient calculation (natural error function is used)

INPUT PARAMETERS:
    Network -   network initialized with one of the network creation funcs
    X       -   input vector, length of array must be at least NIn
    DesiredY-   desired outputs, length of array must be at least NOut
    Grad    -   possibly preallocated array. If size of array is smaller
                than WCount, it will be reallocated. It is recommended to
                reuse previously allocated array to reduce allocation
                overhead.

OUTPUT PARAMETERS:
    E       -   error function, sum-of-squares for regression networks,
                cross-entropy for classification networks.
    Grad    -   gradient of E with respect to weights of network, array[WCount]

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgradn(multilayerperceptron* network,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* desiredy,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_state *_state)
{
    double s;
    ae_int_t i;
    ae_int_t nout;
    ae_int_t ntotal;

    *e = 0;

    
    /*
     * Alloc
     */
    if( grad->cnt<network->structinfo.ptr.p_int[4] )
    {
        ae_vector_set_length(grad, network->structinfo.ptr.p_int[4], _state);
    }
    
    /*
     * Prepare dError/dOut, internal structures
     */
    mlpprocess(network, x, &network->y, _state);
    nout = network->structinfo.ptr.p_int[2];
    ntotal = network->structinfo.ptr.p_int[3];
    for(i=0; i<=ntotal-1; i++)
    {
        network->derror.ptr.p_double[i] = 0;
    }
    *e = 0;
    if( network->structinfo.ptr.p_int[6]==0 )
    {
        
        /*
         * Regression network, least squares
         */
        for(i=0; i<=nout-1; i++)
        {
            network->derror.ptr.p_double[ntotal-nout+i] = network->y.ptr.p_double[i]-desiredy->ptr.p_double[i];
            *e = *e+ae_sqr(network->y.ptr.p_double[i]-desiredy->ptr.p_double[i], _state)/2;
        }
    }
    else
    {
        
        /*
         * Classification network, cross-entropy
         */
        s = 0;
        for(i=0; i<=nout-1; i++)
        {
            s = s+desiredy->ptr.p_double[i];
        }
        for(i=0; i<=nout-1; i++)
        {
            network->derror.ptr.p_double[ntotal-nout+i] = s*network->y.ptr.p_double[i]-desiredy->ptr.p_double[i];
            *e = *e+mlpbase_safecrossentropy(desiredy->ptr.p_double[i], network->y.ptr.p_double[i], _state);
        }
    }
    
    /*
     * gradient
     */
    mlpbase_mlpinternalcalculategradient(network, &network->neurons, &network->weights, &network->derror, grad, ae_true, _state);
}


/*************************************************************************
Batch gradient calculation for a set of inputs/outputs

INPUT PARAMETERS:
    Network -   network initialized with one of the network creation funcs
    XY      -   set of inputs/outputs; one sample = one row;
                first NIn columns contain inputs,
                next NOut columns - desired outputs.
    SSize   -   number of elements in XY
    Grad    -   possibly preallocated array. If size of array is smaller
                than WCount, it will be reallocated. It is recommended to
                reuse previously allocated array to reduce allocation
                overhead.

OUTPUT PARAMETERS:
    E       -   error function, SUM(sqr(y[i]-desiredy[i])/2,i)
    Grad    -   gradient of E with respect to weights of network, array[WCount]

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgradbatch(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;

    *e = 0;

    mlpproperties(network, &nin, &nout, &wcount, _state);
    for(i=0; i<=wcount-1; i++)
    {
        grad->ptr.p_double[i] = 0;
    }
    *e = 0;
    i = 0;
    while(i<=ssize-1)
    {
        mlpbase_mlpchunkedgradient(network, xy, i, ae_minint(ssize, i+mlpbase_chunksize, _state)-i, e, grad, ae_false, _state);
        i = i+mlpbase_chunksize;
    }
}


/*************************************************************************
Batch gradient calculation for a set of inputs/outputs
(natural error function is used)

INPUT PARAMETERS:
    Network -   network initialized with one of the network creation funcs
    XY      -   set of inputs/outputs; one sample = one row;
                first NIn columns contain inputs,
                next NOut columns - desired outputs.
    SSize   -   number of elements in XY
    Grad    -   possibly preallocated array. If size of array is smaller
                than WCount, it will be reallocated. It is recommended to
                reuse previously allocated array to reduce allocation
                overhead.

OUTPUT PARAMETERS:
    E       -   error function, sum-of-squares for regression networks,
                cross-entropy for classification networks.
    Grad    -   gradient of E with respect to weights of network, array[WCount]

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpgradnbatch(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;

    *e = 0;

    mlpproperties(network, &nin, &nout, &wcount, _state);
    for(i=0; i<=wcount-1; i++)
    {
        grad->ptr.p_double[i] = 0;
    }
    *e = 0;
    i = 0;
    while(i<=ssize-1)
    {
        mlpbase_mlpchunkedgradient(network, xy, i, ae_minint(ssize, i+mlpbase_chunksize, _state)-i, e, grad, ae_true, _state);
        i = i+mlpbase_chunksize;
    }
}


/*************************************************************************
Batch Hessian calculation (natural error function) using R-algorithm.
Internal subroutine.

  -- ALGLIB --
     Copyright 26.01.2008 by Bochkanov Sergey.
     
     Hessian calculation based on R-algorithm described in
     "Fast Exact Multiplication by the Hessian",
     B. A. Pearlmutter,
     Neural Computation, 1994.
*************************************************************************/
void mlphessiannbatch(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     double* e,
     /* Real    */ ae_vector* grad,
     /* Real    */ ae_matrix* h,
     ae_state *_state)
{

    *e = 0;

    mlpbase_mlphessianbatchinternal(network, xy, ssize, ae_true, e, grad, h, _state);
}


/*************************************************************************
Batch Hessian calculation using R-algorithm.
Internal subroutine.

  -- ALGLIB --
     Copyright 26.01.2008 by Bochkanov Sergey.

     Hessian calculation based on R-algorithm described in
     "Fast Exact Multiplication by the Hessian",
     B. A. Pearlmutter,
     Neural Computation, 1994.
*************************************************************************/
void mlphessianbatch(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     double* e,
     /* Real    */ ae_vector* grad,
     /* Real    */ ae_matrix* h,
     ae_state *_state)
{

    *e = 0;

    mlpbase_mlphessianbatchinternal(network, xy, ssize, ae_false, e, grad, h, _state);
}


/*************************************************************************
Internal subroutine, shouldn't be called by user.
*************************************************************************/
void mlpinternalprocessvector(/* Integer */ ae_vector* structinfo,
     /* Real    */ ae_vector* weights,
     /* Real    */ ae_vector* columnmeans,
     /* Real    */ ae_vector* columnsigmas,
     /* Real    */ ae_vector* neurons,
     /* Real    */ ae_vector* dfdnet,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t n1;
    ae_int_t n2;
    ae_int_t w1;
    ae_int_t w2;
    ae_int_t ntotal;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t istart;
    ae_int_t offs;
    double net;
    double f;
    double df;
    double d2f;
    double mx;
    ae_bool perr;


    
    /*
     * Read network geometry
     */
    nin = structinfo->ptr.p_int[1];
    nout = structinfo->ptr.p_int[2];
    ntotal = structinfo->ptr.p_int[3];
    istart = structinfo->ptr.p_int[5];
    
    /*
     * Inputs standartisation and putting in the network
     */
    for(i=0; i<=nin-1; i++)
    {
        if( ae_fp_neq(columnsigmas->ptr.p_double[i],0) )
        {
            neurons->ptr.p_double[i] = (x->ptr.p_double[i]-columnmeans->ptr.p_double[i])/columnsigmas->ptr.p_double[i];
        }
        else
        {
            neurons->ptr.p_double[i] = x->ptr.p_double[i]-columnmeans->ptr.p_double[i];
        }
    }
    
    /*
     * Process network
     */
    for(i=0; i<=ntotal-1; i++)
    {
        offs = istart+i*mlpbase_nfieldwidth;
        if( structinfo->ptr.p_int[offs+0]>0||structinfo->ptr.p_int[offs+0]==-5 )
        {
            
            /*
             * Activation function
             */
            mlpactivationfunction(neurons->ptr.p_double[structinfo->ptr.p_int[offs+2]], structinfo->ptr.p_int[offs+0], &f, &df, &d2f, _state);
            neurons->ptr.p_double[i] = f;
            dfdnet->ptr.p_double[i] = df;
            continue;
        }
        if( structinfo->ptr.p_int[offs+0]==0 )
        {
            
            /*
             * Adaptive summator
             */
            n1 = structinfo->ptr.p_int[offs+2];
            n2 = n1+structinfo->ptr.p_int[offs+1]-1;
            w1 = structinfo->ptr.p_int[offs+3];
            w2 = w1+structinfo->ptr.p_int[offs+1]-1;
            net = ae_v_dotproduct(&weights->ptr.p_double[w1], 1, &neurons->ptr.p_double[n1], 1, ae_v_len(w1,w2));
            neurons->ptr.p_double[i] = net;
            dfdnet->ptr.p_double[i] = 1.0;
            continue;
        }
        if( structinfo->ptr.p_int[offs+0]<0 )
        {
            perr = ae_true;
            if( structinfo->ptr.p_int[offs+0]==-2 )
            {
                
                /*
                 * input neuron, left unchanged
                 */
                perr = ae_false;
            }
            if( structinfo->ptr.p_int[offs+0]==-3 )
            {
                
                /*
                 * "-1" neuron
                 */
                neurons->ptr.p_double[i] = -1;
                perr = ae_false;
            }
            if( structinfo->ptr.p_int[offs+0]==-4 )
            {
                
                /*
                 * "0" neuron
                 */
                neurons->ptr.p_double[i] = 0;
                perr = ae_false;
            }
            ae_assert(!perr, "MLPInternalProcessVector: internal error - unknown neuron type!", _state);
            continue;
        }
    }
    
    /*
     * Extract result
     */
    ae_v_move(&y->ptr.p_double[0], 1, &neurons->ptr.p_double[ntotal-nout], 1, ae_v_len(0,nout-1));
    
    /*
     * Softmax post-processing or standardisation if needed
     */
    ae_assert(structinfo->ptr.p_int[6]==0||structinfo->ptr.p_int[6]==1, "MLPInternalProcessVector: unknown normalization type!", _state);
    if( structinfo->ptr.p_int[6]==1 )
    {
        
        /*
         * Softmax
         */
        mx = y->ptr.p_double[0];
        for(i=1; i<=nout-1; i++)
        {
            mx = ae_maxreal(mx, y->ptr.p_double[i], _state);
        }
        net = 0;
        for(i=0; i<=nout-1; i++)
        {
            y->ptr.p_double[i] = ae_exp(y->ptr.p_double[i]-mx, _state);
            net = net+y->ptr.p_double[i];
        }
        for(i=0; i<=nout-1; i++)
        {
            y->ptr.p_double[i] = y->ptr.p_double[i]/net;
        }
    }
    else
    {
        
        /*
         * Standardisation
         */
        for(i=0; i<=nout-1; i++)
        {
            y->ptr.p_double[i] = y->ptr.p_double[i]*columnsigmas->ptr.p_double[nin+i]+columnmeans->ptr.p_double[nin+i];
        }
    }
}


/*************************************************************************
Serializer: allocation

  -- ALGLIB --
     Copyright 14.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpalloc(ae_serializer* s,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t fkind;
    double threshold;
    double v0;
    double v1;
    ae_int_t nin;
    ae_int_t nout;


    nin = network->hllayersizes.ptr.p_int[0];
    nout = network->hllayersizes.ptr.p_int[network->hllayersizes.cnt-1];
    ae_serializer_alloc_entry(s);
    ae_serializer_alloc_entry(s);
    ae_serializer_alloc_entry(s);
    allocintegerarray(s, &network->hllayersizes, -1, _state);
    for(i=1; i<=network->hllayersizes.cnt-1; i++)
    {
        for(j=0; j<=network->hllayersizes.ptr.p_int[i]-1; j++)
        {
            mlpgetneuroninfo(network, i, j, &fkind, &threshold, _state);
            ae_serializer_alloc_entry(s);
            ae_serializer_alloc_entry(s);
            for(k=0; k<=network->hllayersizes.ptr.p_int[i-1]-1; k++)
            {
                ae_serializer_alloc_entry(s);
            }
        }
    }
    for(j=0; j<=nin-1; j++)
    {
        mlpgetinputscaling(network, j, &v0, &v1, _state);
        ae_serializer_alloc_entry(s);
        ae_serializer_alloc_entry(s);
    }
    for(j=0; j<=nout-1; j++)
    {
        mlpgetoutputscaling(network, j, &v0, &v1, _state);
        ae_serializer_alloc_entry(s);
        ae_serializer_alloc_entry(s);
    }
}


/*************************************************************************
Serializer: serialization

  -- ALGLIB --
     Copyright 14.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpserialize(ae_serializer* s,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t fkind;
    double threshold;
    double v0;
    double v1;
    ae_int_t nin;
    ae_int_t nout;


    nin = network->hllayersizes.ptr.p_int[0];
    nout = network->hllayersizes.ptr.p_int[network->hllayersizes.cnt-1];
    ae_serializer_serialize_int(s, getmlpserializationcode(_state), _state);
    ae_serializer_serialize_int(s, mlpbase_mlpfirstversion, _state);
    ae_serializer_serialize_bool(s, mlpissoftmax(network, _state), _state);
    serializeintegerarray(s, &network->hllayersizes, -1, _state);
    for(i=1; i<=network->hllayersizes.cnt-1; i++)
    {
        for(j=0; j<=network->hllayersizes.ptr.p_int[i]-1; j++)
        {
            mlpgetneuroninfo(network, i, j, &fkind, &threshold, _state);
            ae_serializer_serialize_int(s, fkind, _state);
            ae_serializer_serialize_double(s, threshold, _state);
            for(k=0; k<=network->hllayersizes.ptr.p_int[i-1]-1; k++)
            {
                ae_serializer_serialize_double(s, mlpgetweight(network, i-1, k, i, j, _state), _state);
            }
        }
    }
    for(j=0; j<=nin-1; j++)
    {
        mlpgetinputscaling(network, j, &v0, &v1, _state);
        ae_serializer_serialize_double(s, v0, _state);
        ae_serializer_serialize_double(s, v1, _state);
    }
    for(j=0; j<=nout-1; j++)
    {
        mlpgetoutputscaling(network, j, &v0, &v1, _state);
        ae_serializer_serialize_double(s, v0, _state);
        ae_serializer_serialize_double(s, v1, _state);
    }
}


/*************************************************************************
Serializer: unserialization

  -- ALGLIB --
     Copyright 14.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpunserialize(ae_serializer* s,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i0;
    ae_int_t i1;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t fkind;
    double threshold;
    double v0;
    double v1;
    ae_int_t nin;
    ae_int_t nout;
    ae_bool issoftmax;
    ae_vector layersizes;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&layersizes, 0, DT_INT, _state, ae_true);

    
    /*
     * check correctness of header
     */
    ae_serializer_unserialize_int(s, &i0, _state);
    ae_assert(i0==getmlpserializationcode(_state), "MLPUnserialize: stream header corrupted", _state);
    ae_serializer_unserialize_int(s, &i1, _state);
    ae_assert(i1==mlpbase_mlpfirstversion, "MLPUnserialize: stream header corrupted", _state);
    
    /*
     * Create network
     */
    ae_serializer_unserialize_bool(s, &issoftmax, _state);
    unserializeintegerarray(s, &layersizes, _state);
    ae_assert((layersizes.cnt==2||layersizes.cnt==3)||layersizes.cnt==4, "MLPUnserialize: too many hidden layers!", _state);
    nin = layersizes.ptr.p_int[0];
    nout = layersizes.ptr.p_int[layersizes.cnt-1];
    if( layersizes.cnt==2 )
    {
        if( issoftmax )
        {
            mlpcreatec0(layersizes.ptr.p_int[0], layersizes.ptr.p_int[1], network, _state);
        }
        else
        {
            mlpcreate0(layersizes.ptr.p_int[0], layersizes.ptr.p_int[1], network, _state);
        }
    }
    if( layersizes.cnt==3 )
    {
        if( issoftmax )
        {
            mlpcreatec1(layersizes.ptr.p_int[0], layersizes.ptr.p_int[1], layersizes.ptr.p_int[2], network, _state);
        }
        else
        {
            mlpcreate1(layersizes.ptr.p_int[0], layersizes.ptr.p_int[1], layersizes.ptr.p_int[2], network, _state);
        }
    }
    if( layersizes.cnt==4 )
    {
        if( issoftmax )
        {
            mlpcreatec2(layersizes.ptr.p_int[0], layersizes.ptr.p_int[1], layersizes.ptr.p_int[2], layersizes.ptr.p_int[3], network, _state);
        }
        else
        {
            mlpcreate2(layersizes.ptr.p_int[0], layersizes.ptr.p_int[1], layersizes.ptr.p_int[2], layersizes.ptr.p_int[3], network, _state);
        }
    }
    
    /*
     * Load neurons and weights
     */
    for(i=1; i<=layersizes.cnt-1; i++)
    {
        for(j=0; j<=layersizes.ptr.p_int[i]-1; j++)
        {
            ae_serializer_unserialize_int(s, &fkind, _state);
            ae_serializer_unserialize_double(s, &threshold, _state);
            mlpsetneuroninfo(network, i, j, fkind, threshold, _state);
            for(k=0; k<=layersizes.ptr.p_int[i-1]-1; k++)
            {
                ae_serializer_unserialize_double(s, &v0, _state);
                mlpsetweight(network, i-1, k, i, j, v0, _state);
            }
        }
    }
    
    /*
     * Load standartizator
     */
    for(j=0; j<=nin-1; j++)
    {
        ae_serializer_unserialize_double(s, &v0, _state);
        ae_serializer_unserialize_double(s, &v1, _state);
        mlpsetinputscaling(network, j, v0, v1, _state);
    }
    for(j=0; j<=nout-1; j++)
    {
        ae_serializer_unserialize_double(s, &v0, _state);
        ae_serializer_unserialize_double(s, &v1, _state);
        mlpsetoutputscaling(network, j, v0, v1, _state);
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Internal subroutine: adding new input layer to network
*************************************************************************/
static void mlpbase_addinputlayer(ae_int_t ncount,
     /* Integer */ ae_vector* lsizes,
     /* Integer */ ae_vector* ltypes,
     /* Integer */ ae_vector* lconnfirst,
     /* Integer */ ae_vector* lconnlast,
     ae_int_t* lastproc,
     ae_state *) //_state)
{


    lsizes->ptr.p_int[0] = ncount;
    ltypes->ptr.p_int[0] = -2;
    lconnfirst->ptr.p_int[0] = 0;
    lconnlast->ptr.p_int[0] = 0;
    *lastproc = 0;
}


/*************************************************************************
Internal subroutine: adding new summator layer to network
*************************************************************************/
static void mlpbase_addbiasedsummatorlayer(ae_int_t ncount,
     /* Integer */ ae_vector* lsizes,
     /* Integer */ ae_vector* ltypes,
     /* Integer */ ae_vector* lconnfirst,
     /* Integer */ ae_vector* lconnlast,
     ae_int_t* lastproc,
     ae_state *) //_state)
{


    lsizes->ptr.p_int[*lastproc+1] = 1;
    ltypes->ptr.p_int[*lastproc+1] = -3;
    lconnfirst->ptr.p_int[*lastproc+1] = 0;
    lconnlast->ptr.p_int[*lastproc+1] = 0;
    lsizes->ptr.p_int[*lastproc+2] = ncount;
    ltypes->ptr.p_int[*lastproc+2] = 0;
    lconnfirst->ptr.p_int[*lastproc+2] = *lastproc;
    lconnlast->ptr.p_int[*lastproc+2] = *lastproc+1;
    *lastproc = *lastproc+2;
}


/*************************************************************************
Internal subroutine: adding new summator layer to network
*************************************************************************/
static void mlpbase_addactivationlayer(ae_int_t functype,
     /* Integer */ ae_vector* lsizes,
     /* Integer */ ae_vector* ltypes,
     /* Integer */ ae_vector* lconnfirst,
     /* Integer */ ae_vector* lconnlast,
     ae_int_t* lastproc,
     ae_state *_state)
{


    ae_assert(functype>0||functype==-5, "AddActivationLayer: incorrect function type", _state);
    lsizes->ptr.p_int[*lastproc+1] = lsizes->ptr.p_int[*lastproc];
    ltypes->ptr.p_int[*lastproc+1] = functype;
    lconnfirst->ptr.p_int[*lastproc+1] = *lastproc;
    lconnlast->ptr.p_int[*lastproc+1] = *lastproc;
    *lastproc = *lastproc+1;
}


/*************************************************************************
Internal subroutine: adding new zero layer to network
*************************************************************************/
static void mlpbase_addzerolayer(/* Integer */ ae_vector* lsizes,
     /* Integer */ ae_vector* ltypes,
     /* Integer */ ae_vector* lconnfirst,
     /* Integer */ ae_vector* lconnlast,
     ae_int_t* lastproc,
     ae_state *) //_state)
{


    lsizes->ptr.p_int[*lastproc+1] = 1;
    ltypes->ptr.p_int[*lastproc+1] = -4;
    lconnfirst->ptr.p_int[*lastproc+1] = 0;
    lconnlast->ptr.p_int[*lastproc+1] = 0;
    *lastproc = *lastproc+1;
}


/*************************************************************************
This routine adds input layer to the high-level description of the network.

It modifies Network.HLConnections and Network.HLNeurons  and  assumes that
these  arrays  have  enough  place  to  store  data.  It accepts following
parameters:
    Network     -   network
    ConnIdx     -   index of the first free entry in the HLConnections
    NeuroIdx    -   index of the first free entry in the HLNeurons
    StructInfoIdx-  index of the first entry in the low level description
                    of the current layer (in the StructInfo array)
    NIn         -   number of inputs
                    
It modified Network and indices.
*************************************************************************/
static void mlpbase_hladdinputlayer(multilayerperceptron* network,
     ae_int_t* ,//connidx,
     ae_int_t* neuroidx,
     ae_int_t* structinfoidx,
     ae_int_t nin,
     ae_state *) //_state)
{
    ae_int_t i;
    ae_int_t offs;


    offs = mlpbase_hlnfieldwidth*(*neuroidx);
    for(i=0; i<=nin-1; i++)
    {
        network->hlneurons.ptr.p_int[offs+0] = 0;
        network->hlneurons.ptr.p_int[offs+1] = i;
        network->hlneurons.ptr.p_int[offs+2] = -1;
        network->hlneurons.ptr.p_int[offs+3] = -1;
        offs = offs+mlpbase_hlnfieldwidth;
    }
    *neuroidx = *neuroidx+nin;
    *structinfoidx = *structinfoidx+nin;
}


/*************************************************************************
This routine adds output layer to the high-level description of
the network.

It modifies Network.HLConnections and Network.HLNeurons  and  assumes that
these  arrays  have  enough  place  to  store  data.  It accepts following
parameters:
    Network     -   network
    ConnIdx     -   index of the first free entry in the HLConnections
    NeuroIdx    -   index of the first free entry in the HLNeurons
    StructInfoIdx-  index of the first entry in the low level description
                    of the current layer (in the StructInfo array)
    WeightsIdx  -   index of the first entry in the Weights array which
                    corresponds to the current layer
    K           -   current layer index
    NPrev       -   number of neurons in the previous layer
    NOut        -   number of outputs
    IsCls       -   is it classifier network?
    IsLinear    -   is it network with linear output?

It modified Network and ConnIdx/NeuroIdx/StructInfoIdx/WeightsIdx.
*************************************************************************/
static void mlpbase_hladdoutputlayer(multilayerperceptron* network,
     ae_int_t* connidx,
     ae_int_t* neuroidx,
     ae_int_t* structinfoidx,
     ae_int_t* weightsidx,
     ae_int_t k,
     ae_int_t nprev,
     ae_int_t nout,
     ae_bool iscls,
     ae_bool islinearout,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t neurooffs;
    ae_int_t connoffs;


    ae_assert((iscls&&islinearout)||!iscls, "HLAddOutputLayer: internal error", _state);
    neurooffs = mlpbase_hlnfieldwidth*(*neuroidx);
    connoffs = mlpbase_hlconnfieldwidth*(*connidx);
    if( !iscls )
    {
        
        /*
         * Regression network
         */
        for(i=0; i<=nout-1; i++)
        {
            network->hlneurons.ptr.p_int[neurooffs+0] = k;
            network->hlneurons.ptr.p_int[neurooffs+1] = i;
            network->hlneurons.ptr.p_int[neurooffs+2] = *structinfoidx+1+nout+i;
            network->hlneurons.ptr.p_int[neurooffs+3] = *weightsidx+nprev+(nprev+1)*i;
            neurooffs = neurooffs+mlpbase_hlnfieldwidth;
        }
        for(i=0; i<=nprev-1; i++)
        {
            for(j=0; j<=nout-1; j++)
            {
                network->hlconnections.ptr.p_int[connoffs+0] = k-1;
                network->hlconnections.ptr.p_int[connoffs+1] = i;
                network->hlconnections.ptr.p_int[connoffs+2] = k;
                network->hlconnections.ptr.p_int[connoffs+3] = j;
                network->hlconnections.ptr.p_int[connoffs+4] = *weightsidx+i+j*(nprev+1);
                connoffs = connoffs+mlpbase_hlconnfieldwidth;
            }
        }
        *connidx = *connidx+nprev*nout;
        *neuroidx = *neuroidx+nout;
        *structinfoidx = *structinfoidx+2*nout+1;
        *weightsidx = *weightsidx+nout*(nprev+1);
    }
    else
    {
        
        /*
         * Classification network
         */
        for(i=0; i<=nout-2; i++)
        {
            network->hlneurons.ptr.p_int[neurooffs+0] = k;
            network->hlneurons.ptr.p_int[neurooffs+1] = i;
            network->hlneurons.ptr.p_int[neurooffs+2] = -1;
            network->hlneurons.ptr.p_int[neurooffs+3] = *weightsidx+nprev+(nprev+1)*i;
            neurooffs = neurooffs+mlpbase_hlnfieldwidth;
        }
        network->hlneurons.ptr.p_int[neurooffs+0] = k;
        network->hlneurons.ptr.p_int[neurooffs+1] = i;
        network->hlneurons.ptr.p_int[neurooffs+2] = -1;
        network->hlneurons.ptr.p_int[neurooffs+3] = -1;
        for(i=0; i<=nprev-1; i++)
        {
            for(j=0; j<=nout-2; j++)
            {
                network->hlconnections.ptr.p_int[connoffs+0] = k-1;
                network->hlconnections.ptr.p_int[connoffs+1] = i;
                network->hlconnections.ptr.p_int[connoffs+2] = k;
                network->hlconnections.ptr.p_int[connoffs+3] = j;
                network->hlconnections.ptr.p_int[connoffs+4] = *weightsidx+i+j*(nprev+1);
                connoffs = connoffs+mlpbase_hlconnfieldwidth;
            }
        }
        *connidx = *connidx+nprev*(nout-1);
        *neuroidx = *neuroidx+nout;
        *structinfoidx = *structinfoidx+nout+2;
        *weightsidx = *weightsidx+(nout-1)*(nprev+1);
    }
}


/*************************************************************************
This routine adds hidden layer to the high-level description of
the network.

It modifies Network.HLConnections and Network.HLNeurons  and  assumes that
these  arrays  have  enough  place  to  store  data.  It accepts following
parameters:
    Network     -   network
    ConnIdx     -   index of the first free entry in the HLConnections
    NeuroIdx    -   index of the first free entry in the HLNeurons
    StructInfoIdx-  index of the first entry in the low level description
                    of the current layer (in the StructInfo array)
    WeightsIdx  -   index of the first entry in the Weights array which
                    corresponds to the current layer
    K           -   current layer index
    NPrev       -   number of neurons in the previous layer
    NCur        -   number of neurons in the current layer

It modified Network and ConnIdx/NeuroIdx/StructInfoIdx/WeightsIdx.
*************************************************************************/
static void mlpbase_hladdhiddenlayer(multilayerperceptron* network,
     ae_int_t* connidx,
     ae_int_t* neuroidx,
     ae_int_t* structinfoidx,
     ae_int_t* weightsidx,
     ae_int_t k,
     ae_int_t nprev,
     ae_int_t ncur,
     ae_state *) //_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t neurooffs;
    ae_int_t connoffs;


    neurooffs = mlpbase_hlnfieldwidth*(*neuroidx);
    connoffs = mlpbase_hlconnfieldwidth*(*connidx);
    for(i=0; i<=ncur-1; i++)
    {
        network->hlneurons.ptr.p_int[neurooffs+0] = k;
        network->hlneurons.ptr.p_int[neurooffs+1] = i;
        network->hlneurons.ptr.p_int[neurooffs+2] = *structinfoidx+1+ncur+i;
        network->hlneurons.ptr.p_int[neurooffs+3] = *weightsidx+nprev+(nprev+1)*i;
        neurooffs = neurooffs+mlpbase_hlnfieldwidth;
    }
    for(i=0; i<=nprev-1; i++)
    {
        for(j=0; j<=ncur-1; j++)
        {
            network->hlconnections.ptr.p_int[connoffs+0] = k-1;
            network->hlconnections.ptr.p_int[connoffs+1] = i;
            network->hlconnections.ptr.p_int[connoffs+2] = k;
            network->hlconnections.ptr.p_int[connoffs+3] = j;
            network->hlconnections.ptr.p_int[connoffs+4] = *weightsidx+i+j*(nprev+1);
            connoffs = connoffs+mlpbase_hlconnfieldwidth;
        }
    }
    *connidx = *connidx+nprev*ncur;
    *neuroidx = *neuroidx+ncur;
    *structinfoidx = *structinfoidx+2*ncur+1;
    *weightsidx = *weightsidx+ncur*(nprev+1);
}


/*************************************************************************
This function fills high level information about network created using
internal MLPCreate() function.

This function does NOT examine StructInfo for low level information, it
just expects that network has following structure:

    input neuron            \
    ...                      | input layer
    input neuron            /
    
    "-1" neuron             \
    biased summator          |
    ...                      |
    biased summator          | hidden layer(s), if there are exists any
    activation function      |
    ...                      |
    activation function     /
    
    "-1" neuron            \
    biased summator         | output layer:
    ...                     |
    biased summator         | * we have NOut summators/activators for regression networks
    activation function     | * we have only NOut-1 summators and no activators for classifiers
    ...                     | * we have "0" neuron only when we have classifier
    activation function     |
    "0" neuron              /


  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
static void mlpbase_fillhighlevelinformation(multilayerperceptron* network,
     ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_bool iscls,
     ae_bool islinearout,
     ae_state *_state)
{
    ae_int_t idxweights;
    ae_int_t idxstruct;
    ae_int_t idxneuro;
    ae_int_t idxconn;


    ae_assert((iscls&&islinearout)||!iscls, "FillHighLevelInformation: internal error", _state);
    
    /*
     * Preparations common to all types of networks
     */
    idxweights = 0;
    idxneuro = 0;
    idxstruct = 0;
    idxconn = 0;
    network->hlnetworktype = 0;
    
    /*
     * network without hidden layers
     */
    if( nhid1==0 )
    {
        ae_vector_set_length(&network->hllayersizes, 2, _state);
        network->hllayersizes.ptr.p_int[0] = nin;
        network->hllayersizes.ptr.p_int[1] = nout;
        if( !iscls )
        {
            ae_vector_set_length(&network->hlconnections, mlpbase_hlconnfieldwidth*nin*nout, _state);
            ae_vector_set_length(&network->hlneurons, mlpbase_hlnfieldwidth*(nin+nout), _state);
            network->hlnormtype = 0;
        }
        else
        {
            ae_vector_set_length(&network->hlconnections, mlpbase_hlconnfieldwidth*nin*(nout-1), _state);
            ae_vector_set_length(&network->hlneurons, mlpbase_hlnfieldwidth*(nin+nout), _state);
            network->hlnormtype = 1;
        }
        mlpbase_hladdinputlayer(network, &idxconn, &idxneuro, &idxstruct, nin, _state);
        mlpbase_hladdoutputlayer(network, &idxconn, &idxneuro, &idxstruct, &idxweights, 1, nin, nout, iscls, islinearout, _state);
        return;
    }
    
    /*
     * network with one hidden layers
     */
    if( nhid2==0 )
    {
        ae_vector_set_length(&network->hllayersizes, 3, _state);
        network->hllayersizes.ptr.p_int[0] = nin;
        network->hllayersizes.ptr.p_int[1] = nhid1;
        network->hllayersizes.ptr.p_int[2] = nout;
        if( !iscls )
        {
            ae_vector_set_length(&network->hlconnections, mlpbase_hlconnfieldwidth*(nin*nhid1+nhid1*nout), _state);
            ae_vector_set_length(&network->hlneurons, mlpbase_hlnfieldwidth*(nin+nhid1+nout), _state);
            network->hlnormtype = 0;
        }
        else
        {
            ae_vector_set_length(&network->hlconnections, mlpbase_hlconnfieldwidth*(nin*nhid1+nhid1*(nout-1)), _state);
            ae_vector_set_length(&network->hlneurons, mlpbase_hlnfieldwidth*(nin+nhid1+nout), _state);
            network->hlnormtype = 1;
        }
        mlpbase_hladdinputlayer(network, &idxconn, &idxneuro, &idxstruct, nin, _state);
        mlpbase_hladdhiddenlayer(network, &idxconn, &idxneuro, &idxstruct, &idxweights, 1, nin, nhid1, _state);
        mlpbase_hladdoutputlayer(network, &idxconn, &idxneuro, &idxstruct, &idxweights, 2, nhid1, nout, iscls, islinearout, _state);
        return;
    }
    
    /*
     * Two hidden layers
     */
    ae_vector_set_length(&network->hllayersizes, 4, _state);
    network->hllayersizes.ptr.p_int[0] = nin;
    network->hllayersizes.ptr.p_int[1] = nhid1;
    network->hllayersizes.ptr.p_int[2] = nhid2;
    network->hllayersizes.ptr.p_int[3] = nout;
    if( !iscls )
    {
        ae_vector_set_length(&network->hlconnections, mlpbase_hlconnfieldwidth*(nin*nhid1+nhid1*nhid2+nhid2*nout), _state);
        ae_vector_set_length(&network->hlneurons, mlpbase_hlnfieldwidth*(nin+nhid1+nhid2+nout), _state);
        network->hlnormtype = 0;
    }
    else
    {
        ae_vector_set_length(&network->hlconnections, mlpbase_hlconnfieldwidth*(nin*nhid1+nhid1*nhid2+nhid2*(nout-1)), _state);
        ae_vector_set_length(&network->hlneurons, mlpbase_hlnfieldwidth*(nin+nhid1+nhid2+nout), _state);
        network->hlnormtype = 1;
    }
    mlpbase_hladdinputlayer(network, &idxconn, &idxneuro, &idxstruct, nin, _state);
    mlpbase_hladdhiddenlayer(network, &idxconn, &idxneuro, &idxstruct, &idxweights, 1, nin, nhid1, _state);
    mlpbase_hladdhiddenlayer(network, &idxconn, &idxneuro, &idxstruct, &idxweights, 2, nhid1, nhid2, _state);
    mlpbase_hladdoutputlayer(network, &idxconn, &idxneuro, &idxstruct, &idxweights, 3, nhid2, nout, iscls, islinearout, _state);
}


/*************************************************************************
Internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
static void mlpbase_mlpcreate(ae_int_t nin,
     ae_int_t nout,
     /* Integer */ ae_vector* lsizes,
     /* Integer */ ae_vector* ltypes,
     /* Integer */ ae_vector* lconnfirst,
     /* Integer */ ae_vector* lconnlast,
     ae_int_t layerscount,
     ae_bool isclsnet,
     multilayerperceptron* network,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_int_t ssize;
    ae_int_t ntotal;
    ae_int_t wcount;
    ae_int_t offs;
    ae_int_t nprocessed;
    ae_int_t wallocated;
    ae_vector localtemp;
    ae_vector lnfirst;
    ae_vector lnsyn;

    ae_frame_make(_state, &_frame_block);
    _multilayerperceptron_clear(network);
    ae_vector_init(&localtemp, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lnfirst, 0, DT_INT, _state, ae_true);
    ae_vector_init(&lnsyn, 0, DT_INT, _state, ae_true);

    
    /*
     * Check
     */
    ae_assert(layerscount>0, "MLPCreate: wrong parameters!", _state);
    ae_assert(ltypes->ptr.p_int[0]==-2, "MLPCreate: wrong LTypes[0] (must be -2)!", _state);
    for(i=0; i<=layerscount-1; i++)
    {
        ae_assert(lsizes->ptr.p_int[i]>0, "MLPCreate: wrong LSizes!", _state);
        ae_assert(lconnfirst->ptr.p_int[i]>=0&&(lconnfirst->ptr.p_int[i]<i||i==0), "MLPCreate: wrong LConnFirst!", _state);
        ae_assert(lconnlast->ptr.p_int[i]>=lconnfirst->ptr.p_int[i]&&(lconnlast->ptr.p_int[i]<i||i==0), "MLPCreate: wrong LConnLast!", _state);
    }
    
    /*
     * Build network geometry
     */
    ae_vector_set_length(&lnfirst, layerscount-1+1, _state);
    ae_vector_set_length(&lnsyn, layerscount-1+1, _state);
    ntotal = 0;
    wcount = 0;
    for(i=0; i<=layerscount-1; i++)
    {
        
        /*
         * Analyze connections.
         * This code must throw an assertion in case of unknown LTypes[I]
         */
        lnsyn.ptr.p_int[i] = -1;
        if( ltypes->ptr.p_int[i]>=0||ltypes->ptr.p_int[i]==-5 )
        {
            lnsyn.ptr.p_int[i] = 0;
            for(j=lconnfirst->ptr.p_int[i]; j<=lconnlast->ptr.p_int[i]; j++)
            {
                lnsyn.ptr.p_int[i] = lnsyn.ptr.p_int[i]+lsizes->ptr.p_int[j];
            }
        }
        else
        {
            if( (ltypes->ptr.p_int[i]==-2||ltypes->ptr.p_int[i]==-3)||ltypes->ptr.p_int[i]==-4 )
            {
                lnsyn.ptr.p_int[i] = 0;
            }
        }
        ae_assert(lnsyn.ptr.p_int[i]>=0, "MLPCreate: internal error #0!", _state);
        
        /*
         * Other info
         */
        lnfirst.ptr.p_int[i] = ntotal;
        ntotal = ntotal+lsizes->ptr.p_int[i];
        if( ltypes->ptr.p_int[i]==0 )
        {
            wcount = wcount+lnsyn.ptr.p_int[i]*lsizes->ptr.p_int[i];
        }
    }
    ssize = 7+ntotal*mlpbase_nfieldwidth;
    
    /*
     * Allocate
     */
    ae_vector_set_length(&network->structinfo, ssize-1+1, _state);
    ae_vector_set_length(&network->weights, wcount-1+1, _state);
    if( isclsnet )
    {
        ae_vector_set_length(&network->columnmeans, nin-1+1, _state);
        ae_vector_set_length(&network->columnsigmas, nin-1+1, _state);
    }
    else
    {
        ae_vector_set_length(&network->columnmeans, nin+nout-1+1, _state);
        ae_vector_set_length(&network->columnsigmas, nin+nout-1+1, _state);
    }
    ae_vector_set_length(&network->neurons, ntotal-1+1, _state);
    ae_matrix_set_length(&network->chunks, 3*ntotal+1, mlpbase_chunksize-1+1, _state);
    ae_vector_set_length(&network->nwbuf, ae_maxint(wcount, 2*nout, _state)-1+1, _state);
    ae_vector_set_length(&network->integerbuf, 3+1, _state);
    ae_vector_set_length(&network->dfdnet, ntotal-1+1, _state);
    ae_vector_set_length(&network->x, nin-1+1, _state);
    ae_vector_set_length(&network->y, nout-1+1, _state);
    ae_vector_set_length(&network->derror, ntotal-1+1, _state);
    
    /*
     * Fill structure: global info
     */
    network->structinfo.ptr.p_int[0] = ssize;
    network->structinfo.ptr.p_int[1] = nin;
    network->structinfo.ptr.p_int[2] = nout;
    network->structinfo.ptr.p_int[3] = ntotal;
    network->structinfo.ptr.p_int[4] = wcount;
    network->structinfo.ptr.p_int[5] = 7;
    if( isclsnet )
    {
        network->structinfo.ptr.p_int[6] = 1;
    }
    else
    {
        network->structinfo.ptr.p_int[6] = 0;
    }
    
    /*
     * Fill structure: neuron connections
     */
    nprocessed = 0;
    wallocated = 0;
    for(i=0; i<=layerscount-1; i++)
    {
        for(j=0; j<=lsizes->ptr.p_int[i]-1; j++)
        {
            offs = network->structinfo.ptr.p_int[5]+nprocessed*mlpbase_nfieldwidth;
            network->structinfo.ptr.p_int[offs+0] = ltypes->ptr.p_int[i];
            if( ltypes->ptr.p_int[i]==0 )
            {
                
                /*
                 * Adaptive summator:
                 * * connections with weights to previous neurons
                 */
                network->structinfo.ptr.p_int[offs+1] = lnsyn.ptr.p_int[i];
                network->structinfo.ptr.p_int[offs+2] = lnfirst.ptr.p_int[lconnfirst->ptr.p_int[i]];
                network->structinfo.ptr.p_int[offs+3] = wallocated;
                wallocated = wallocated+lnsyn.ptr.p_int[i];
                nprocessed = nprocessed+1;
            }
            if( ltypes->ptr.p_int[i]>0||ltypes->ptr.p_int[i]==-5 )
            {
                
                /*
                 * Activation layer:
                 * * each neuron connected to one (only one) of previous neurons.
                 * * no weights
                 */
                network->structinfo.ptr.p_int[offs+1] = 1;
                network->structinfo.ptr.p_int[offs+2] = lnfirst.ptr.p_int[lconnfirst->ptr.p_int[i]]+j;
                network->structinfo.ptr.p_int[offs+3] = -1;
                nprocessed = nprocessed+1;
            }
            if( (ltypes->ptr.p_int[i]==-2||ltypes->ptr.p_int[i]==-3)||ltypes->ptr.p_int[i]==-4 )
            {
                nprocessed = nprocessed+1;
            }
        }
    }
    ae_assert(wallocated==wcount, "MLPCreate: internal error #1!", _state);
    ae_assert(nprocessed==ntotal, "MLPCreate: internal error #2!", _state);
    
    /*
     * Fill weights by small random values
     * Initialize means and sigmas
     */
    for(i=0; i<=wcount-1; i++)
    {
        network->weights.ptr.p_double[i] = ae_randomreal(_state)-0.5;
    }
    for(i=0; i<=nin-1; i++)
    {
        network->columnmeans.ptr.p_double[i] = 0;
        network->columnsigmas.ptr.p_double[i] = 1;
    }
    if( !isclsnet )
    {
        for(i=0; i<=nout-1; i++)
        {
            network->columnmeans.ptr.p_double[nin+i] = 0;
            network->columnsigmas.ptr.p_double[nin+i] = 1;
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Internal subroutine for Hessian calculation.

WARNING!!! Unspeakable math far beyong human capabilities :)
*************************************************************************/
static void mlpbase_mlphessianbatchinternal(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     ae_bool naturalerr,
     double* e,
     /* Real    */ ae_vector* grad,
     /* Real    */ ae_matrix* h,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    ae_int_t ntotal;
    ae_int_t istart;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t kl;
    ae_int_t offs;
    ae_int_t n1;
    ae_int_t n2;
    ae_int_t w1;
    ae_int_t w2;
    double s;
    double t;
    double v;
    double et;
    ae_bool bflag;
    double f;
    double df;
    double d2f;
    double deidyj;
    double mx;
    double q;
    double z;
    double s2;
    double expi;
    double expj;
    ae_vector x;
    ae_vector desiredy;
    ae_vector gt;
    ae_vector zeros;
    ae_matrix rx;
    ae_matrix ry;
    ae_matrix rdx;
    ae_matrix rdy;

    ae_frame_make(_state, &_frame_block);
    *e = 0;
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&desiredy, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&gt, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&zeros, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&rx, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&ry, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&rdx, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&rdy, 0, 0, DT_REAL, _state, ae_true);

    mlpproperties(network, &nin, &nout, &wcount, _state);
    ntotal = network->structinfo.ptr.p_int[3];
    istart = network->structinfo.ptr.p_int[5];
    
    /*
     * Prepare
     */
    ae_vector_set_length(&x, nin-1+1, _state);
    ae_vector_set_length(&desiredy, nout-1+1, _state);
    ae_vector_set_length(&zeros, wcount-1+1, _state);
    ae_vector_set_length(&gt, wcount-1+1, _state);
    ae_matrix_set_length(&rx, ntotal+nout-1+1, wcount-1+1, _state);
    ae_matrix_set_length(&ry, ntotal+nout-1+1, wcount-1+1, _state);
    ae_matrix_set_length(&rdx, ntotal+nout-1+1, wcount-1+1, _state);
    ae_matrix_set_length(&rdy, ntotal+nout-1+1, wcount-1+1, _state);
    *e = 0;
    for(i=0; i<=wcount-1; i++)
    {
        zeros.ptr.p_double[i] = 0;
    }
    ae_v_move(&grad->ptr.p_double[0], 1, &zeros.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
    for(i=0; i<=wcount-1; i++)
    {
        ae_v_move(&h->ptr.pp_double[i][0], 1, &zeros.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
    }
    
    /*
     * Process
     */
    for(k=0; k<=ssize-1; k++)
    {
        
        /*
         * Process vector with MLPGradN.
         * Now Neurons, DFDNET and DError contains results of the last run.
         */
        ae_v_move(&x.ptr.p_double[0], 1, &xy->ptr.pp_double[k][0], 1, ae_v_len(0,nin-1));
        if( mlpissoftmax(network, _state) )
        {
            
            /*
             * class labels outputs
             */
            kl = ae_round(xy->ptr.pp_double[k][nin], _state);
            for(i=0; i<=nout-1; i++)
            {
                if( i==kl )
                {
                    desiredy.ptr.p_double[i] = 1;
                }
                else
                {
                    desiredy.ptr.p_double[i] = 0;
                }
            }
        }
        else
        {
            
            /*
             * real outputs
             */
            ae_v_move(&desiredy.ptr.p_double[0], 1, &xy->ptr.pp_double[k][nin], 1, ae_v_len(0,nout-1));
        }
        if( naturalerr )
        {
            mlpgradn(network, &x, &desiredy, &et, &gt, _state);
        }
        else
        {
            mlpgrad(network, &x, &desiredy, &et, &gt, _state);
        }
        
        /*
         * grad, error
         */
        *e = *e+et;
        ae_v_add(&grad->ptr.p_double[0], 1, &gt.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        
        /*
         * Hessian.
         * Forward pass of the R-algorithm
         */
        for(i=0; i<=ntotal-1; i++)
        {
            offs = istart+i*mlpbase_nfieldwidth;
            ae_v_move(&rx.ptr.pp_double[i][0], 1, &zeros.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            ae_v_move(&ry.ptr.pp_double[i][0], 1, &zeros.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            if( network->structinfo.ptr.p_int[offs+0]>0||network->structinfo.ptr.p_int[offs+0]==-5 )
            {
                
                /*
                 * Activation function
                 */
                n1 = network->structinfo.ptr.p_int[offs+2];
                ae_v_move(&rx.ptr.pp_double[i][0], 1, &ry.ptr.pp_double[n1][0], 1, ae_v_len(0,wcount-1));
                v = network->dfdnet.ptr.p_double[i];
                ae_v_moved(&ry.ptr.pp_double[i][0], 1, &rx.ptr.pp_double[i][0], 1, ae_v_len(0,wcount-1), v);
                continue;
            }
            if( network->structinfo.ptr.p_int[offs+0]==0 )
            {
                
                /*
                 * Adaptive summator
                 */
                n1 = network->structinfo.ptr.p_int[offs+2];
                n2 = n1+network->structinfo.ptr.p_int[offs+1]-1;
                w1 = network->structinfo.ptr.p_int[offs+3];
                w2 = w1+network->structinfo.ptr.p_int[offs+1]-1;
                for(j=n1; j<=n2; j++)
                {
                    v = network->weights.ptr.p_double[w1+j-n1];
                    ae_v_addd(&rx.ptr.pp_double[i][0], 1, &ry.ptr.pp_double[j][0], 1, ae_v_len(0,wcount-1), v);
                    rx.ptr.pp_double[i][w1+j-n1] = rx.ptr.pp_double[i][w1+j-n1]+network->neurons.ptr.p_double[j];
                }
                ae_v_move(&ry.ptr.pp_double[i][0], 1, &rx.ptr.pp_double[i][0], 1, ae_v_len(0,wcount-1));
                continue;
            }
            if( network->structinfo.ptr.p_int[offs+0]<0 )
            {
                bflag = ae_true;
                if( network->structinfo.ptr.p_int[offs+0]==-2 )
                {
                    
                    /*
                     * input neuron, left unchanged
                     */
                    bflag = ae_false;
                }
                if( network->structinfo.ptr.p_int[offs+0]==-3 )
                {
                    
                    /*
                     * "-1" neuron, left unchanged
                     */
                    bflag = ae_false;
                }
                if( network->structinfo.ptr.p_int[offs+0]==-4 )
                {
                    
                    /*
                     * "0" neuron, left unchanged
                     */
                    bflag = ae_false;
                }
                ae_assert(!bflag, "MLPHessianNBatch: internal error - unknown neuron type!", _state);
                continue;
            }
        }
        
        /*
         * Hessian. Backward pass of the R-algorithm.
         *
         * Stage 1. Initialize RDY
         */
        for(i=0; i<=ntotal+nout-1; i++)
        {
            ae_v_move(&rdy.ptr.pp_double[i][0], 1, &zeros.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        }
        if( network->structinfo.ptr.p_int[6]==0 )
        {
            
            /*
             * Standardisation.
             *
             * In context of the Hessian calculation standardisation
             * is considered as additional layer with weightless
             * activation function:
             *
             * F(NET) := Sigma*NET
             *
             * So we add one more layer to forward pass, and
             * make forward/backward pass through this layer.
             */
            for(i=0; i<=nout-1; i++)
            {
                n1 = ntotal-nout+i;
                n2 = ntotal+i;
                
                /*
                 * Forward pass from N1 to N2
                 */
                ae_v_move(&rx.ptr.pp_double[n2][0], 1, &ry.ptr.pp_double[n1][0], 1, ae_v_len(0,wcount-1));
                v = network->columnsigmas.ptr.p_double[nin+i];
                ae_v_moved(&ry.ptr.pp_double[n2][0], 1, &rx.ptr.pp_double[n2][0], 1, ae_v_len(0,wcount-1), v);
                
                /*
                 * Initialization of RDY
                 */
                ae_v_move(&rdy.ptr.pp_double[n2][0], 1, &ry.ptr.pp_double[n2][0], 1, ae_v_len(0,wcount-1));
                
                /*
                 * Backward pass from N2 to N1:
                 * 1. Calculate R(dE/dX).
                 * 2. No R(dE/dWij) is needed since weight of activation neuron
                 *    is fixed to 1. So we can update R(dE/dY) for
                 *    the connected neuron (note that Vij=0, Wij=1)
                 */
                df = network->columnsigmas.ptr.p_double[nin+i];
                ae_v_moved(&rdx.ptr.pp_double[n2][0], 1, &rdy.ptr.pp_double[n2][0], 1, ae_v_len(0,wcount-1), df);
                ae_v_add(&rdy.ptr.pp_double[n1][0], 1, &rdx.ptr.pp_double[n2][0], 1, ae_v_len(0,wcount-1));
            }
        }
        else
        {
            
            /*
             * Softmax.
             *
             * Initialize RDY using generalized expression for ei'(yi)
             * (see expression (9) from p. 5 of "Fast Exact Multiplication by the Hessian").
             *
             * When we are working with softmax network, generalized
             * expression for ei'(yi) is used because softmax
             * normalization leads to ei, which depends on all y's
             */
            if( naturalerr )
            {
                
                /*
                 * softmax + cross-entropy.
                 * We have:
                 *
                 * S = sum(exp(yk)),
                 * ei = sum(trn)*exp(yi)/S-trn_i
                 *
                 * j=i:   d(ei)/d(yj) = T*exp(yi)*(S-exp(yi))/S^2
                 * j<>i:  d(ei)/d(yj) = -T*exp(yi)*exp(yj)/S^2
                 */
                t = 0;
                for(i=0; i<=nout-1; i++)
                {
                    t = t+desiredy.ptr.p_double[i];
                }
                mx = network->neurons.ptr.p_double[ntotal-nout];
                for(i=0; i<=nout-1; i++)
                {
                    mx = ae_maxreal(mx, network->neurons.ptr.p_double[ntotal-nout+i], _state);
                }
                s = 0;
                for(i=0; i<=nout-1; i++)
                {
                    network->nwbuf.ptr.p_double[i] = ae_exp(network->neurons.ptr.p_double[ntotal-nout+i]-mx, _state);
                    s = s+network->nwbuf.ptr.p_double[i];
                }
                for(i=0; i<=nout-1; i++)
                {
                    for(j=0; j<=nout-1; j++)
                    {
                        if( j==i )
                        {
                            deidyj = t*network->nwbuf.ptr.p_double[i]*(s-network->nwbuf.ptr.p_double[i])/ae_sqr(s, _state);
                            ae_v_addd(&rdy.ptr.pp_double[ntotal-nout+i][0], 1, &ry.ptr.pp_double[ntotal-nout+i][0], 1, ae_v_len(0,wcount-1), deidyj);
                        }
                        else
                        {
                            deidyj = -t*network->nwbuf.ptr.p_double[i]*network->nwbuf.ptr.p_double[j]/ae_sqr(s, _state);
                            ae_v_addd(&rdy.ptr.pp_double[ntotal-nout+i][0], 1, &ry.ptr.pp_double[ntotal-nout+j][0], 1, ae_v_len(0,wcount-1), deidyj);
                        }
                    }
                }
            }
            else
            {
                
                /*
                 * For a softmax + squared error we have expression
                 * far beyond human imagination so we dont even try
                 * to comment on it. Just enjoy the code...
                 *
                 * P.S. That's why "natural error" is called "natural" -
                 * compact beatiful expressions, fast code....
                 */
                mx = network->neurons.ptr.p_double[ntotal-nout];
                for(i=0; i<=nout-1; i++)
                {
                    mx = ae_maxreal(mx, network->neurons.ptr.p_double[ntotal-nout+i], _state);
                }
                s = 0;
                s2 = 0;
                for(i=0; i<=nout-1; i++)
                {
                    network->nwbuf.ptr.p_double[i] = ae_exp(network->neurons.ptr.p_double[ntotal-nout+i]-mx, _state);
                    s = s+network->nwbuf.ptr.p_double[i];
                    s2 = s2+ae_sqr(network->nwbuf.ptr.p_double[i], _state);
                }
                q = 0;
                for(i=0; i<=nout-1; i++)
                {
                    q = q+(network->y.ptr.p_double[i]-desiredy.ptr.p_double[i])*network->nwbuf.ptr.p_double[i];
                }
                for(i=0; i<=nout-1; i++)
                {
                    z = -q+(network->y.ptr.p_double[i]-desiredy.ptr.p_double[i])*s;
                    expi = network->nwbuf.ptr.p_double[i];
                    for(j=0; j<=nout-1; j++)
                    {
                        expj = network->nwbuf.ptr.p_double[j];
                        if( j==i )
                        {
                            deidyj = expi/ae_sqr(s, _state)*((z+expi)*(s-2*expi)/s+expi*s2/ae_sqr(s, _state));
                        }
                        else
                        {
                            deidyj = expi*expj/ae_sqr(s, _state)*(s2/ae_sqr(s, _state)-2*z/s-(expi+expj)/s+(network->y.ptr.p_double[i]-desiredy.ptr.p_double[i])-(network->y.ptr.p_double[j]-desiredy.ptr.p_double[j]));
                        }
                        ae_v_addd(&rdy.ptr.pp_double[ntotal-nout+i][0], 1, &ry.ptr.pp_double[ntotal-nout+j][0], 1, ae_v_len(0,wcount-1), deidyj);
                    }
                }
            }
        }
        
        /*
         * Hessian. Backward pass of the R-algorithm
         *
         * Stage 2. Process.
         */
        for(i=ntotal-1; i>=0; i--)
        {
            
            /*
             * Possible variants:
             * 1. Activation function
             * 2. Adaptive summator
             * 3. Special neuron
             */
            offs = istart+i*mlpbase_nfieldwidth;
            if( network->structinfo.ptr.p_int[offs+0]>0||network->structinfo.ptr.p_int[offs+0]==-5 )
            {
                n1 = network->structinfo.ptr.p_int[offs+2];
                
                /*
                 * First, calculate R(dE/dX).
                 */
                mlpactivationfunction(network->neurons.ptr.p_double[n1], network->structinfo.ptr.p_int[offs+0], &f, &df, &d2f, _state);
                v = d2f*network->derror.ptr.p_double[i];
                ae_v_moved(&rdx.ptr.pp_double[i][0], 1, &rdy.ptr.pp_double[i][0], 1, ae_v_len(0,wcount-1), df);
                ae_v_addd(&rdx.ptr.pp_double[i][0], 1, &rx.ptr.pp_double[i][0], 1, ae_v_len(0,wcount-1), v);
                
                /*
                 * No R(dE/dWij) is needed since weight of activation neuron
                 * is fixed to 1.
                 *
                 * So we can update R(dE/dY) for the connected neuron.
                 * (note that Vij=0, Wij=1)
                 */
                ae_v_add(&rdy.ptr.pp_double[n1][0], 1, &rdx.ptr.pp_double[i][0], 1, ae_v_len(0,wcount-1));
                continue;
            }
            if( network->structinfo.ptr.p_int[offs+0]==0 )
            {
                
                /*
                 * Adaptive summator
                 */
                n1 = network->structinfo.ptr.p_int[offs+2];
                n2 = n1+network->structinfo.ptr.p_int[offs+1]-1;
                w1 = network->structinfo.ptr.p_int[offs+3];
                w2 = w1+network->structinfo.ptr.p_int[offs+1]-1;
                
                /*
                 * First, calculate R(dE/dX).
                 */
                ae_v_move(&rdx.ptr.pp_double[i][0], 1, &rdy.ptr.pp_double[i][0], 1, ae_v_len(0,wcount-1));
                
                /*
                 * Then, calculate R(dE/dWij)
                 */
                for(j=w1; j<=w2; j++)
                {
                    v = network->neurons.ptr.p_double[n1+j-w1];
                    ae_v_addd(&h->ptr.pp_double[j][0], 1, &rdx.ptr.pp_double[i][0], 1, ae_v_len(0,wcount-1), v);
                    v = network->derror.ptr.p_double[i];
                    ae_v_addd(&h->ptr.pp_double[j][0], 1, &ry.ptr.pp_double[n1+j-w1][0], 1, ae_v_len(0,wcount-1), v);
                }
                
                /*
                 * And finally, update R(dE/dY) for connected neurons.
                 */
                for(j=w1; j<=w2; j++)
                {
                    v = network->weights.ptr.p_double[j];
                    ae_v_addd(&rdy.ptr.pp_double[n1+j-w1][0], 1, &rdx.ptr.pp_double[i][0], 1, ae_v_len(0,wcount-1), v);
                    rdy.ptr.pp_double[n1+j-w1][j] = rdy.ptr.pp_double[n1+j-w1][j]+network->derror.ptr.p_double[i];
                }
                continue;
            }
            if( network->structinfo.ptr.p_int[offs+0]<0 )
            {
                bflag = ae_false;
                if( (network->structinfo.ptr.p_int[offs+0]==-2||network->structinfo.ptr.p_int[offs+0]==-3)||network->structinfo.ptr.p_int[offs+0]==-4 )
                {
                    
                    /*
                     * Special neuron type, no back-propagation required
                     */
                    bflag = ae_true;
                }
                ae_assert(bflag, "MLPHessianNBatch: unknown neuron type!", _state);
                continue;
            }
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Internal subroutine

Network must be processed by MLPProcess on X
*************************************************************************/
static void mlpbase_mlpinternalcalculategradient(multilayerperceptron* network,
     /* Real    */ ae_vector* neurons,
     /* Real    */ ae_vector* weights,
     /* Real    */ ae_vector* derror,
     /* Real    */ ae_vector* grad,
     ae_bool naturalerrorfunc,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t n1;
    ae_int_t n2;
    ae_int_t w1;
    ae_int_t w2;
    ae_int_t ntotal;
    ae_int_t istart;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t offs;
    double dedf;
    double dfdnet;
    double v;
    double fown;
    double deown;
    double net;
    double mx;
    ae_bool bflag;


    
    /*
     * Read network geometry
     */
    nin = network->structinfo.ptr.p_int[1];
    nout = network->structinfo.ptr.p_int[2];
    ntotal = network->structinfo.ptr.p_int[3];
    istart = network->structinfo.ptr.p_int[5];
    
    /*
     * Pre-processing of dError/dOut:
     * from dError/dOut(normalized) to dError/dOut(non-normalized)
     */
    ae_assert(network->structinfo.ptr.p_int[6]==0||network->structinfo.ptr.p_int[6]==1, "MLPInternalCalculateGradient: unknown normalization type!", _state);
    if( network->structinfo.ptr.p_int[6]==1 )
    {
        
        /*
         * Softmax
         */
        if( !naturalerrorfunc )
        {
            mx = network->neurons.ptr.p_double[ntotal-nout];
            for(i=0; i<=nout-1; i++)
            {
                mx = ae_maxreal(mx, network->neurons.ptr.p_double[ntotal-nout+i], _state);
            }
            net = 0;
            for(i=0; i<=nout-1; i++)
            {
                network->nwbuf.ptr.p_double[i] = ae_exp(network->neurons.ptr.p_double[ntotal-nout+i]-mx, _state);
                net = net+network->nwbuf.ptr.p_double[i];
            }
            v = ae_v_dotproduct(&network->derror.ptr.p_double[ntotal-nout], 1, &network->nwbuf.ptr.p_double[0], 1, ae_v_len(ntotal-nout,ntotal-1));
            for(i=0; i<=nout-1; i++)
            {
                fown = network->nwbuf.ptr.p_double[i];
                deown = network->derror.ptr.p_double[ntotal-nout+i];
                network->nwbuf.ptr.p_double[nout+i] = (-v+deown*fown+deown*(net-fown))*fown/ae_sqr(net, _state);
            }
            for(i=0; i<=nout-1; i++)
            {
                network->derror.ptr.p_double[ntotal-nout+i] = network->nwbuf.ptr.p_double[nout+i];
            }
        }
    }
    else
    {
        
        /*
         * Un-standardisation
         */
        for(i=0; i<=nout-1; i++)
        {
            network->derror.ptr.p_double[ntotal-nout+i] = network->derror.ptr.p_double[ntotal-nout+i]*network->columnsigmas.ptr.p_double[nin+i];
        }
    }
    
    /*
     * Backpropagation
     */
    for(i=ntotal-1; i>=0; i--)
    {
        
        /*
         * Extract info
         */
        offs = istart+i*mlpbase_nfieldwidth;
        if( network->structinfo.ptr.p_int[offs+0]>0||network->structinfo.ptr.p_int[offs+0]==-5 )
        {
            
            /*
             * Activation function
             */
            dedf = network->derror.ptr.p_double[i];
            dfdnet = network->dfdnet.ptr.p_double[i];
            derror->ptr.p_double[network->structinfo.ptr.p_int[offs+2]] = derror->ptr.p_double[network->structinfo.ptr.p_int[offs+2]]+dedf*dfdnet;
            continue;
        }
        if( network->structinfo.ptr.p_int[offs+0]==0 )
        {
            
            /*
             * Adaptive summator
             */
            n1 = network->structinfo.ptr.p_int[offs+2];
            n2 = n1+network->structinfo.ptr.p_int[offs+1]-1;
            w1 = network->structinfo.ptr.p_int[offs+3];
            w2 = w1+network->structinfo.ptr.p_int[offs+1]-1;
            dedf = network->derror.ptr.p_double[i];
            dfdnet = 1.0;
            v = dedf*dfdnet;
            ae_v_moved(&grad->ptr.p_double[w1], 1, &neurons->ptr.p_double[n1], 1, ae_v_len(w1,w2), v);
            ae_v_addd(&derror->ptr.p_double[n1], 1, &weights->ptr.p_double[w1], 1, ae_v_len(n1,n2), v);
            continue;
        }
        if( network->structinfo.ptr.p_int[offs+0]<0 )
        {
            bflag = ae_false;
            if( (network->structinfo.ptr.p_int[offs+0]==-2||network->structinfo.ptr.p_int[offs+0]==-3)||network->structinfo.ptr.p_int[offs+0]==-4 )
            {
                
                /*
                 * Special neuron type, no back-propagation required
                 */
                bflag = ae_true;
            }
            ae_assert(bflag, "MLPInternalCalculateGradient: unknown neuron type!", _state);
            continue;
        }
    }
}


/*************************************************************************
Internal subroutine, chunked gradient
*************************************************************************/
static void mlpbase_mlpchunkedgradient(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t cstart,
     ae_int_t csize,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_bool naturalerrorfunc,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t kl;
    ae_int_t n1;
    ae_int_t n2;
    ae_int_t w1;
    ae_int_t w2;
    ae_int_t c1;
    ae_int_t c2;
    ae_int_t ntotal;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t offs;
    double f;
    double df;
    double d2f;
    double v;
    double s;
    double fown;
    double deown;
    double net;
    double lnnet;
    double mx;
    ae_bool bflag;
    ae_int_t istart;
    ae_int_t ineurons;
    ae_int_t idfdnet;
    ae_int_t iderror;
    ae_int_t izeros;


    
    /*
     * Read network geometry, prepare data
     */
    nin = network->structinfo.ptr.p_int[1];
    nout = network->structinfo.ptr.p_int[2];
    ntotal = network->structinfo.ptr.p_int[3];
    istart = network->structinfo.ptr.p_int[5];
    c1 = cstart;
    c2 = cstart+csize-1;
    ineurons = 0;
    idfdnet = ntotal;
    iderror = 2*ntotal;
    izeros = 3*ntotal;
    for(j=0; j<=csize-1; j++)
    {
        network->chunks.ptr.pp_double[izeros][j] = 0;
    }
    
    /*
     * Forward pass:
     * 1. Load inputs from XY to Chunks[0:NIn-1,0:CSize-1]
     * 2. Forward pass
     */
    for(i=0; i<=nin-1; i++)
    {
        for(j=0; j<=csize-1; j++)
        {
            if( ae_fp_neq(network->columnsigmas.ptr.p_double[i],0) )
            {
                network->chunks.ptr.pp_double[i][j] = (xy->ptr.pp_double[c1+j][i]-network->columnmeans.ptr.p_double[i])/network->columnsigmas.ptr.p_double[i];
            }
            else
            {
                network->chunks.ptr.pp_double[i][j] = xy->ptr.pp_double[c1+j][i]-network->columnmeans.ptr.p_double[i];
            }
        }
    }
    for(i=0; i<=ntotal-1; i++)
    {
        offs = istart+i*mlpbase_nfieldwidth;
        if( network->structinfo.ptr.p_int[offs+0]>0||network->structinfo.ptr.p_int[offs+0]==-5 )
        {
            
            /*
             * Activation function:
             * * calculate F vector, F(i) = F(NET(i))
             */
            n1 = network->structinfo.ptr.p_int[offs+2];
            ae_v_move(&network->chunks.ptr.pp_double[i][0], 1, &network->chunks.ptr.pp_double[n1][0], 1, ae_v_len(0,csize-1));
            for(j=0; j<=csize-1; j++)
            {
                mlpactivationfunction(network->chunks.ptr.pp_double[i][j], network->structinfo.ptr.p_int[offs+0], &f, &df, &d2f, _state);
                network->chunks.ptr.pp_double[i][j] = f;
                network->chunks.ptr.pp_double[idfdnet+i][j] = df;
            }
            continue;
        }
        if( network->structinfo.ptr.p_int[offs+0]==0 )
        {
            
            /*
             * Adaptive summator:
             * * calculate NET vector, NET(i) = SUM(W(j,i)*Neurons(j),j=N1..N2)
             */
            n1 = network->structinfo.ptr.p_int[offs+2];
            n2 = n1+network->structinfo.ptr.p_int[offs+1]-1;
            w1 = network->structinfo.ptr.p_int[offs+3];
            w2 = w1+network->structinfo.ptr.p_int[offs+1]-1;
            ae_v_move(&network->chunks.ptr.pp_double[i][0], 1, &network->chunks.ptr.pp_double[izeros][0], 1, ae_v_len(0,csize-1));
            for(j=n1; j<=n2; j++)
            {
                v = network->weights.ptr.p_double[w1+j-n1];
                ae_v_addd(&network->chunks.ptr.pp_double[i][0], 1, &network->chunks.ptr.pp_double[j][0], 1, ae_v_len(0,csize-1), v);
            }
            continue;
        }
        if( network->structinfo.ptr.p_int[offs+0]<0 )
        {
            bflag = ae_false;
            if( network->structinfo.ptr.p_int[offs+0]==-2 )
            {
                
                /*
                 * input neuron, left unchanged
                 */
                bflag = ae_true;
            }
            if( network->structinfo.ptr.p_int[offs+0]==-3 )
            {
                
                /*
                 * "-1" neuron
                 */
                for(k=0; k<=csize-1; k++)
                {
                    network->chunks.ptr.pp_double[i][k] = -1;
                }
                bflag = ae_true;
            }
            if( network->structinfo.ptr.p_int[offs+0]==-4 )
            {
                
                /*
                 * "0" neuron
                 */
                for(k=0; k<=csize-1; k++)
                {
                    network->chunks.ptr.pp_double[i][k] = 0;
                }
                bflag = ae_true;
            }
            ae_assert(bflag, "MLPChunkedGradient: internal error - unknown neuron type!", _state);
            continue;
        }
    }
    
    /*
     * Post-processing, error, dError/dOut
     */
    for(i=0; i<=ntotal-1; i++)
    {
        ae_v_move(&network->chunks.ptr.pp_double[iderror+i][0], 1, &network->chunks.ptr.pp_double[izeros][0], 1, ae_v_len(0,csize-1));
    }
    ae_assert(network->structinfo.ptr.p_int[6]==0||network->structinfo.ptr.p_int[6]==1, "MLPChunkedGradient: unknown normalization type!", _state);
    if( network->structinfo.ptr.p_int[6]==1 )
    {
        
        /*
         * Softmax output, classification network.
         *
         * For each K = 0..CSize-1 do:
         * 1. place exp(outputs[k]) to NWBuf[0:NOut-1]
         * 2. place sum(exp(..)) to NET
         * 3. calculate dError/dOut and place it to the second block of Chunks
         */
        for(k=0; k<=csize-1; k++)
        {
            
            /*
             * Normalize
             */
            mx = network->chunks.ptr.pp_double[ntotal-nout][k];
            for(i=1; i<=nout-1; i++)
            {
                mx = ae_maxreal(mx, network->chunks.ptr.pp_double[ntotal-nout+i][k], _state);
            }
            net = 0;
            for(i=0; i<=nout-1; i++)
            {
                network->nwbuf.ptr.p_double[i] = ae_exp(network->chunks.ptr.pp_double[ntotal-nout+i][k]-mx, _state);
                net = net+network->nwbuf.ptr.p_double[i];
            }
            
            /*
             * Calculate error function and dError/dOut
             */
            if( naturalerrorfunc )
            {
                
                /*
                 * Natural error func.
                 *
                 */
                s = 1;
                lnnet = ae_log(net, _state);
                kl = ae_round(xy->ptr.pp_double[cstart+k][nin], _state);
                for(i=0; i<=nout-1; i++)
                {
                    if( i==kl )
                    {
                        v = 1;
                    }
                    else
                    {
                        v = 0;
                    }
                    network->chunks.ptr.pp_double[iderror+ntotal-nout+i][k] = s*network->nwbuf.ptr.p_double[i]/net-v;
                    *e = *e+mlpbase_safecrossentropy(v, network->nwbuf.ptr.p_double[i]/net, _state);
                }
            }
            else
            {
                
                /*
                 * Least squares error func
                 * Error, dError/dOut(normalized)
                 */
                kl = ae_round(xy->ptr.pp_double[cstart+k][nin], _state);
                for(i=0; i<=nout-1; i++)
                {
                    if( i==kl )
                    {
                        v = network->nwbuf.ptr.p_double[i]/net-1;
                    }
                    else
                    {
                        v = network->nwbuf.ptr.p_double[i]/net;
                    }
                    network->nwbuf.ptr.p_double[nout+i] = v;
                    *e = *e+ae_sqr(v, _state)/2;
                }
                
                /*
                 * From dError/dOut(normalized) to dError/dOut(non-normalized)
                 */
                v = ae_v_dotproduct(&network->nwbuf.ptr.p_double[nout], 1, &network->nwbuf.ptr.p_double[0], 1, ae_v_len(nout,2*nout-1));
                for(i=0; i<=nout-1; i++)
                {
                    fown = network->nwbuf.ptr.p_double[i];
                    deown = network->nwbuf.ptr.p_double[nout+i];
                    network->chunks.ptr.pp_double[iderror+ntotal-nout+i][k] = (-v+deown*fown+deown*(net-fown))*fown/ae_sqr(net, _state);
                }
            }
        }
    }
    else
    {
        
        /*
         * Normal output, regression network
         *
         * For each K = 0..CSize-1 do:
         * 1. calculate dError/dOut and place it to the second block of Chunks
         */
        for(i=0; i<=nout-1; i++)
        {
            for(j=0; j<=csize-1; j++)
            {
                v = network->chunks.ptr.pp_double[ntotal-nout+i][j]*network->columnsigmas.ptr.p_double[nin+i]+network->columnmeans.ptr.p_double[nin+i]-xy->ptr.pp_double[cstart+j][nin+i];
                network->chunks.ptr.pp_double[iderror+ntotal-nout+i][j] = v*network->columnsigmas.ptr.p_double[nin+i];
                *e = *e+ae_sqr(v, _state)/2;
            }
        }
    }
    
    /*
     * Backpropagation
     */
    for(i=ntotal-1; i>=0; i--)
    {
        
        /*
         * Extract info
         */
        offs = istart+i*mlpbase_nfieldwidth;
        if( network->structinfo.ptr.p_int[offs+0]>0||network->structinfo.ptr.p_int[offs+0]==-5 )
        {
            
            /*
             * Activation function
             */
            n1 = network->structinfo.ptr.p_int[offs+2];
            for(k=0; k<=csize-1; k++)
            {
                network->chunks.ptr.pp_double[iderror+i][k] = network->chunks.ptr.pp_double[iderror+i][k]*network->chunks.ptr.pp_double[idfdnet+i][k];
            }
            ae_v_add(&network->chunks.ptr.pp_double[iderror+n1][0], 1, &network->chunks.ptr.pp_double[iderror+i][0], 1, ae_v_len(0,csize-1));
            continue;
        }
        if( network->structinfo.ptr.p_int[offs+0]==0 )
        {
            
            /*
             * "Normal" activation function
             */
            n1 = network->structinfo.ptr.p_int[offs+2];
            n2 = n1+network->structinfo.ptr.p_int[offs+1]-1;
            w1 = network->structinfo.ptr.p_int[offs+3];
            w2 = w1+network->structinfo.ptr.p_int[offs+1]-1;
            for(j=w1; j<=w2; j++)
            {
                v = ae_v_dotproduct(&network->chunks.ptr.pp_double[n1+j-w1][0], 1, &network->chunks.ptr.pp_double[iderror+i][0], 1, ae_v_len(0,csize-1));
                grad->ptr.p_double[j] = grad->ptr.p_double[j]+v;
            }
            for(j=n1; j<=n2; j++)
            {
                v = network->weights.ptr.p_double[w1+j-n1];
                ae_v_addd(&network->chunks.ptr.pp_double[iderror+j][0], 1, &network->chunks.ptr.pp_double[iderror+i][0], 1, ae_v_len(0,csize-1), v);
            }
            continue;
        }
        if( network->structinfo.ptr.p_int[offs+0]<0 )
        {
            bflag = ae_false;
            if( (network->structinfo.ptr.p_int[offs+0]==-2||network->structinfo.ptr.p_int[offs+0]==-3)||network->structinfo.ptr.p_int[offs+0]==-4 )
            {
                
                /*
                 * Special neuron type, no back-propagation required
                 */
                bflag = ae_true;
            }
            ae_assert(bflag, "MLPInternalCalculateGradient: unknown neuron type!", _state);
            continue;
        }
    }
}


/*************************************************************************
Returns T*Ln(T/Z), guarded against overflow/underflow.
Internal subroutine.
*************************************************************************/
static double mlpbase_safecrossentropy(double t,
     double z,
     ae_state *_state)
{
    double r;
    double result;


    if( ae_fp_eq(t,0) )
    {
        result = 0;
    }
    else
    {
        if( ae_fp_greater(ae_fabs(z, _state),1) )
        {
            
            /*
             * Shouldn't be the case with softmax,
             * but we just want to be sure.
             */
            if( ae_fp_eq(t/z,0) )
            {
                r = ae_minrealnumber;
            }
            else
            {
                r = t/z;
            }
        }
        else
        {
            
            /*
             * Normal case
             */
            if( ae_fp_eq(z,0)||ae_fp_greater_eq(ae_fabs(t, _state),ae_maxrealnumber*ae_fabs(z, _state)) )
            {
                r = ae_maxrealnumber;
            }
            else
            {
                r = t/z;
            }
        }
        result = t*ae_log(r, _state);
    }
    return result;
}


ae_bool _multilayerperceptron_init(multilayerperceptron* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->hllayersizes, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->hlconnections, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->hlneurons, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->structinfo, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->weights, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->columnmeans, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->columnsigmas, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->neurons, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->dfdnet, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->derror, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->x, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->y, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->chunks, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->nwbuf, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->integerbuf, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _multilayerperceptron_init_copy(multilayerperceptron* dst, multilayerperceptron* src, ae_state *_state, ae_bool make_automatic)
{
    dst->hlnetworktype = src->hlnetworktype;
    dst->hlnormtype = src->hlnormtype;
    if( !ae_vector_init_copy(&dst->hllayersizes, &src->hllayersizes, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->hlconnections, &src->hlconnections, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->hlneurons, &src->hlneurons, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->structinfo, &src->structinfo, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->weights, &src->weights, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->columnmeans, &src->columnmeans, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->columnsigmas, &src->columnsigmas, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->neurons, &src->neurons, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->dfdnet, &src->dfdnet, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->derror, &src->derror, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->x, &src->x, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->y, &src->y, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->chunks, &src->chunks, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->nwbuf, &src->nwbuf, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->integerbuf, &src->integerbuf, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _multilayerperceptron_clear(multilayerperceptron* p)
{
    ae_vector_clear(&p->hllayersizes);
    ae_vector_clear(&p->hlconnections);
    ae_vector_clear(&p->hlneurons);
    ae_vector_clear(&p->structinfo);
    ae_vector_clear(&p->weights);
    ae_vector_clear(&p->columnmeans);
    ae_vector_clear(&p->columnsigmas);
    ae_vector_clear(&p->neurons);
    ae_vector_clear(&p->dfdnet);
    ae_vector_clear(&p->derror);
    ae_vector_clear(&p->x);
    ae_vector_clear(&p->y);
    ae_matrix_clear(&p->chunks);
    ae_vector_clear(&p->nwbuf);
    ae_vector_clear(&p->integerbuf);
}




/*************************************************************************
This subroutine trains logit model.

INPUT PARAMETERS:
    XY          -   training set, array[0..NPoints-1,0..NVars]
                    First NVars columns store values of independent
                    variables, next column stores number of class (from 0
                    to NClasses-1) which dataset element belongs to. Fractional
                    values are rounded to nearest integer.
    NPoints     -   training set size, NPoints>=1
    NVars       -   number of independent variables, NVars>=1
    NClasses    -   number of classes, NClasses>=2

OUTPUT PARAMETERS:
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<NVars+2, NVars<1, NClasses<2).
                    *  1, if task has been solved
    LM          -   model built
    Rep         -   training report

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
void mnltrainh(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t* info,
     logitmodel* lm,
     mnlreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t ssize;
    ae_bool allsame;
    ae_int_t offs;
    double threshold;
    double wminstep;
    double decay;
    ae_int_t wdim;
    ae_int_t expoffs;
    double v;
    double s;
    multilayerperceptron network;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    double e;
    ae_vector g;
    ae_matrix h;
    ae_bool spd;
    ae_vector x;
    ae_vector y;
    ae_vector wbase;
    double wstep;
    ae_vector wdir;
    ae_vector work;
    ae_int_t mcstage;
    logitmcstate mcstate;
    ae_int_t mcinfo;
    ae_int_t mcnfev;
    ae_int_t solverinfo;
    densesolverreport solverrep;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _logitmodel_clear(lm);
    _mnlreport_clear(rep);
    _multilayerperceptron_init(&network, _state, ae_true);
    ae_vector_init(&g, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&h, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&wbase, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&wdir, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&work, 0, DT_REAL, _state, ae_true);
    _logitmcstate_init(&mcstate, _state, ae_true);
    _densesolverreport_init(&solverrep, _state, ae_true);

    threshold = 1000*ae_machineepsilon;
    wminstep = 0.001;
    decay = 0.001;
    
    /*
     * Test for inputs
     */
    if( (npoints<nvars+2||nvars<1)||nclasses<2 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    for(i=0; i<=npoints-1; i++)
    {
        if( ae_round(xy->ptr.pp_double[i][nvars], _state)<0||ae_round(xy->ptr.pp_double[i][nvars], _state)>=nclasses )
        {
            *info = -2;
            ae_frame_leave(_state);
            return;
        }
    }
    *info = 1;
    
    /*
     * Initialize data
     */
    rep->ngrad = 0;
    rep->nhess = 0;
    
    /*
     * Allocate array
     */
    wdim = (nvars+1)*(nclasses-1);
    offs = 5;
    expoffs = offs+wdim;
    ssize = 5+(nvars+1)*(nclasses-1)+nclasses;
    ae_vector_set_length(&lm->w, ssize-1+1, _state);
    lm->w.ptr.p_double[0] = ssize;
    lm->w.ptr.p_double[1] = logit_logitvnum;
    lm->w.ptr.p_double[2] = nvars;
    lm->w.ptr.p_double[3] = nclasses;
    lm->w.ptr.p_double[4] = offs;
    
    /*
     * Degenerate case: all outputs are equal
     */
    allsame = ae_true;
    for(i=1; i<=npoints-1; i++)
    {
        if( ae_round(xy->ptr.pp_double[i][nvars], _state)!=ae_round(xy->ptr.pp_double[i-1][nvars], _state) )
        {
            allsame = ae_false;
        }
    }
    if( allsame )
    {
        for(i=0; i<=(nvars+1)*(nclasses-1)-1; i++)
        {
            lm->w.ptr.p_double[offs+i] = 0;
        }
        v = -2*ae_log(ae_minrealnumber, _state);
        k = ae_round(xy->ptr.pp_double[0][nvars], _state);
        if( k==nclasses-1 )
        {
            for(i=0; i<=nclasses-2; i++)
            {
                lm->w.ptr.p_double[offs+i*(nvars+1)+nvars] = -v;
            }
        }
        else
        {
            for(i=0; i<=nclasses-2; i++)
            {
                if( i==k )
                {
                    lm->w.ptr.p_double[offs+i*(nvars+1)+nvars] = v;
                }
                else
                {
                    lm->w.ptr.p_double[offs+i*(nvars+1)+nvars] = 0;
                }
            }
        }
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * General case.
     * Prepare task and network. Allocate space.
     */
    mlpcreatec0(nvars, nclasses, &network, _state);
    mlpinitpreprocessor(&network, xy, npoints, _state);
    mlpproperties(&network, &nin, &nout, &wcount, _state);
    for(i=0; i<=wcount-1; i++)
    {
        network.weights.ptr.p_double[i] = (2*ae_randomreal(_state)-1)/nvars;
    }
    ae_vector_set_length(&g, wcount-1+1, _state);
    ae_matrix_set_length(&h, wcount-1+1, wcount-1+1, _state);
    ae_vector_set_length(&wbase, wcount-1+1, _state);
    ae_vector_set_length(&wdir, wcount-1+1, _state);
    ae_vector_set_length(&work, wcount-1+1, _state);
    
    /*
     * First stage: optimize in gradient direction.
     */
    for(k=0; k<=wcount/3+10; k++)
    {
        
        /*
         * Calculate gradient in starting point
         */
        mlpgradnbatch(&network, xy, npoints, &e, &g, _state);
        v = ae_v_dotproduct(&network.weights.ptr.p_double[0], 1, &network.weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        e = e+0.5*decay*v;
        ae_v_addd(&g.ptr.p_double[0], 1, &network.weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1), decay);
        rep->ngrad = rep->ngrad+1;
        
        /*
         * Setup optimization scheme
         */
        ae_v_moveneg(&wdir.ptr.p_double[0], 1, &g.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        v = ae_v_dotproduct(&wdir.ptr.p_double[0], 1, &wdir.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        wstep = ae_sqrt(v, _state);
        v = 1/ae_sqrt(v, _state);
        ae_v_muld(&wdir.ptr.p_double[0], 1, ae_v_len(0,wcount-1), v);
        mcstage = 0;
        logit_mnlmcsrch(wcount, &network.weights, &e, &g, &wdir, &wstep, &mcinfo, &mcnfev, &work, &mcstate, &mcstage, _state);
        while(mcstage!=0)
        {
            mlpgradnbatch(&network, xy, npoints, &e, &g, _state);
            v = ae_v_dotproduct(&network.weights.ptr.p_double[0], 1, &network.weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            e = e+0.5*decay*v;
            ae_v_addd(&g.ptr.p_double[0], 1, &network.weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1), decay);
            rep->ngrad = rep->ngrad+1;
            logit_mnlmcsrch(wcount, &network.weights, &e, &g, &wdir, &wstep, &mcinfo, &mcnfev, &work, &mcstate, &mcstage, _state);
        }
    }
    
    /*
     * Second stage: use Hessian when we are close to the minimum
     */
    for(;;)
    {
        
        /*
         * Calculate and update E/G/H
         */
        mlphessiannbatch(&network, xy, npoints, &e, &g, &h, _state);
        v = ae_v_dotproduct(&network.weights.ptr.p_double[0], 1, &network.weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        e = e+0.5*decay*v;
        ae_v_addd(&g.ptr.p_double[0], 1, &network.weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1), decay);
        for(k=0; k<=wcount-1; k++)
        {
            h.ptr.pp_double[k][k] = h.ptr.pp_double[k][k]+decay;
        }
        rep->nhess = rep->nhess+1;
        
        /*
         * Select step direction
         * NOTE: it is important to use lower-triangle Cholesky
         * factorization since it is much faster than higher-triangle version.
         */
        spd = spdmatrixcholesky(&h, wcount, ae_false, _state);
        spdmatrixcholeskysolve(&h, wcount, ae_false, &g, &solverinfo, &solverrep, &wdir, _state);
        spd = solverinfo>0;
        if( spd )
        {
            
            /*
             * H is positive definite.
             * Step in Newton direction.
             */
            ae_v_muld(&wdir.ptr.p_double[0], 1, ae_v_len(0,wcount-1), -1);
            spd = ae_true;
        }
        else
        {
            
            /*
             * H is indefinite.
             * Step in gradient direction.
             */
            ae_v_moveneg(&wdir.ptr.p_double[0], 1, &g.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            spd = ae_false;
        }
        
        /*
         * Optimize in WDir direction
         */
        v = ae_v_dotproduct(&wdir.ptr.p_double[0], 1, &wdir.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        wstep = ae_sqrt(v, _state);
        v = 1/ae_sqrt(v, _state);
        ae_v_muld(&wdir.ptr.p_double[0], 1, ae_v_len(0,wcount-1), v);
        mcstage = 0;
        logit_mnlmcsrch(wcount, &network.weights, &e, &g, &wdir, &wstep, &mcinfo, &mcnfev, &work, &mcstate, &mcstage, _state);
        while(mcstage!=0)
        {
            mlpgradnbatch(&network, xy, npoints, &e, &g, _state);
            v = ae_v_dotproduct(&network.weights.ptr.p_double[0], 1, &network.weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            e = e+0.5*decay*v;
            ae_v_addd(&g.ptr.p_double[0], 1, &network.weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1), decay);
            rep->ngrad = rep->ngrad+1;
            logit_mnlmcsrch(wcount, &network.weights, &e, &g, &wdir, &wstep, &mcinfo, &mcnfev, &work, &mcstate, &mcstage, _state);
        }
        if( spd&&((mcinfo==2||mcinfo==4)||mcinfo==6) )
        {
            break;
        }
    }
    
    /*
     * Convert from NN format to MNL format
     */
    ae_v_move(&lm->w.ptr.p_double[offs], 1, &network.weights.ptr.p_double[0], 1, ae_v_len(offs,offs+wcount-1));
    for(k=0; k<=nvars-1; k++)
    {
        for(i=0; i<=nclasses-2; i++)
        {
            s = network.columnsigmas.ptr.p_double[k];
            if( ae_fp_eq(s,0) )
            {
                s = 1;
            }
            j = offs+(nvars+1)*i;
            v = lm->w.ptr.p_double[j+k];
            lm->w.ptr.p_double[j+k] = v/s;
            lm->w.ptr.p_double[j+nvars] = lm->w.ptr.p_double[j+nvars]+v*network.columnmeans.ptr.p_double[k]/s;
        }
    }
    for(k=0; k<=nclasses-2; k++)
    {
        lm->w.ptr.p_double[offs+(nvars+1)*k+nvars] = -lm->w.ptr.p_double[offs+(nvars+1)*k+nvars];
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Procesing

INPUT PARAMETERS:
    LM      -   logit model, passed by non-constant reference
                (some fields of structure are used as temporaries
                when calculating model output).
    X       -   input vector,  array[0..NVars-1].
    Y       -   (possibly) preallocated buffer; if size of Y is less than
                NClasses, it will be reallocated.If it is large enough, it
                is NOT reallocated, so we can save some time on reallocation.

OUTPUT PARAMETERS:
    Y       -   result, array[0..NClasses-1]
                Vector of posterior probabilities for classification task.

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
void mnlprocess(logitmodel* lm,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t nvars;
    ae_int_t nclasses;
    ae_int_t offs;
    ae_int_t i;
    ae_int_t i1;
    double s;


    ae_assert(ae_fp_eq(lm->w.ptr.p_double[1],logit_logitvnum), "MNLProcess: unexpected model version", _state);
    nvars = ae_round(lm->w.ptr.p_double[2], _state);
    nclasses = ae_round(lm->w.ptr.p_double[3], _state);
    offs = ae_round(lm->w.ptr.p_double[4], _state);
    logit_mnliexp(&lm->w, x, _state);
    s = 0;
    i1 = offs+(nvars+1)*(nclasses-1);
    for(i=i1; i<=i1+nclasses-1; i++)
    {
        s = s+lm->w.ptr.p_double[i];
    }
    if( y->cnt<nclasses )
    {
        ae_vector_set_length(y, nclasses, _state);
    }
    for(i=0; i<=nclasses-1; i++)
    {
        y->ptr.p_double[i] = lm->w.ptr.p_double[i1+i]/s;
    }
}


/*************************************************************************
'interactive'  variant  of  MNLProcess  for  languages  like  Python which
support constructs like "Y = MNLProcess(LM,X)" and interactive mode of the
interpreter

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
void mnlprocessi(logitmodel* lm,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{

    ae_vector_clear(y);

    mnlprocess(lm, x, y, _state);
}


/*************************************************************************
Unpacks coefficients of logit model. Logit model have form:

    P(class=i) = S(i) / (S(0) + S(1) + ... +S(M-1))
          S(i) = Exp(A[i,0]*X[0] + ... + A[i,N-1]*X[N-1] + A[i,N]), when i<M-1
        S(M-1) = 1

INPUT PARAMETERS:
    LM          -   logit model in ALGLIB format

OUTPUT PARAMETERS:
    V           -   coefficients, array[0..NClasses-2,0..NVars]
    NVars       -   number of independent variables
    NClasses    -   number of classes

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
void mnlunpack(logitmodel* lm,
     /* Real    */ ae_matrix* a,
     ae_int_t* nvars,
     ae_int_t* nclasses,
     ae_state *_state)
{
    ae_int_t offs;
    ae_int_t i;

    ae_matrix_clear(a);
    *nvars = 0;
    *nclasses = 0;

    ae_assert(ae_fp_eq(lm->w.ptr.p_double[1],logit_logitvnum), "MNLUnpack: unexpected model version", _state);
    *nvars = ae_round(lm->w.ptr.p_double[2], _state);
    *nclasses = ae_round(lm->w.ptr.p_double[3], _state);
    offs = ae_round(lm->w.ptr.p_double[4], _state);
    ae_matrix_set_length(a, *nclasses-2+1, *nvars+1, _state);
    for(i=0; i<=*nclasses-2; i++)
    {
        ae_v_move(&a->ptr.pp_double[i][0], 1, &lm->w.ptr.p_double[offs+i*(*nvars+1)], 1, ae_v_len(0,*nvars));
    }
}


/*************************************************************************
"Packs" coefficients and creates logit model in ALGLIB format (MNLUnpack
reversed).

INPUT PARAMETERS:
    A           -   model (see MNLUnpack)
    NVars       -   number of independent variables
    NClasses    -   number of classes

OUTPUT PARAMETERS:
    LM          -   logit model.

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
void mnlpack(/* Real    */ ae_matrix* a,
     ae_int_t nvars,
     ae_int_t nclasses,
     logitmodel* lm,
     ae_state *_state)
{
    ae_int_t offs;
    ae_int_t i;
    ae_int_t wdim;
    ae_int_t ssize;

    _logitmodel_clear(lm);

    wdim = (nvars+1)*(nclasses-1);
    offs = 5;
    ssize = 5+(nvars+1)*(nclasses-1)+nclasses;
    ae_vector_set_length(&lm->w, ssize-1+1, _state);
    lm->w.ptr.p_double[0] = ssize;
    lm->w.ptr.p_double[1] = logit_logitvnum;
    lm->w.ptr.p_double[2] = nvars;
    lm->w.ptr.p_double[3] = nclasses;
    lm->w.ptr.p_double[4] = offs;
    for(i=0; i<=nclasses-2; i++)
    {
        ae_v_move(&lm->w.ptr.p_double[offs+i*(nvars+1)], 1, &a->ptr.pp_double[i][0], 1, ae_v_len(offs+i*(nvars+1),offs+i*(nvars+1)+nvars));
    }
}


/*************************************************************************
Copying of LogitModel strucure

INPUT PARAMETERS:
    LM1 -   original

OUTPUT PARAMETERS:
    LM2 -   copy

  -- ALGLIB --
     Copyright 15.03.2009 by Bochkanov Sergey
*************************************************************************/
void mnlcopy(logitmodel* lm1, logitmodel* lm2, ae_state *_state)
{
    ae_int_t k;

    _logitmodel_clear(lm2);

    k = ae_round(lm1->w.ptr.p_double[0], _state);
    ae_vector_set_length(&lm2->w, k-1+1, _state);
    ae_v_move(&lm2->w.ptr.p_double[0], 1, &lm1->w.ptr.p_double[0], 1, ae_v_len(0,k-1));
}


/*************************************************************************
Average cross-entropy (in bits per element) on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    CrossEntropy/(NPoints*ln(2)).

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
double mnlavgce(logitmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t nvars;
    ae_int_t nclasses;
    ae_int_t i;
    ae_vector workx;
    ae_vector worky;
    double result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&workx, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&worky, 0, DT_REAL, _state, ae_true);

    ae_assert(ae_fp_eq(lm->w.ptr.p_double[1],logit_logitvnum), "MNLClsError: unexpected model version", _state);
    nvars = ae_round(lm->w.ptr.p_double[2], _state);
    nclasses = ae_round(lm->w.ptr.p_double[3], _state);
    ae_vector_set_length(&workx, nvars-1+1, _state);
    ae_vector_set_length(&worky, nclasses-1+1, _state);
    result = 0;
    for(i=0; i<=npoints-1; i++)
    {
        ae_assert(ae_round(xy->ptr.pp_double[i][nvars], _state)>=0&&ae_round(xy->ptr.pp_double[i][nvars], _state)<nclasses, "MNLAvgCE: incorrect class number!", _state);
        
        /*
         * Process
         */
        ae_v_move(&workx.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
        mnlprocess(lm, &workx, &worky, _state);
        if( ae_fp_greater(worky.ptr.p_double[ae_round(xy->ptr.pp_double[i][nvars], _state)],0) )
        {
            result = result-ae_log(worky.ptr.p_double[ae_round(xy->ptr.pp_double[i][nvars], _state)], _state);
        }
        else
        {
            result = result-ae_log(ae_minrealnumber, _state);
        }
    }
    result = result/(npoints*ae_log(2, _state));
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
Relative classification error on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    percent of incorrectly classified cases.

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
double mnlrelclserror(logitmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    double result;


    result = (double)mnlclserror(lm, xy, npoints, _state)/(double)npoints;
    return result;
}


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error (error when estimating posterior probabilities).

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double mnlrmserror(logitmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    double relcls;
    double avgce;
    double rms;
    double avg;
    double avgrel;
    double result;


    ae_assert(ae_round(lm->w.ptr.p_double[1], _state)==logit_logitvnum, "MNLRMSError: Incorrect MNL version!", _state);
    logit_mnlallerrors(lm, xy, npoints, &relcls, &avgce, &rms, &avg, &avgrel, _state);
    result = rms;
    return result;
}


/*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average error (error when estimating posterior probabilities).

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double mnlavgerror(logitmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    double relcls;
    double avgce;
    double rms;
    double avg;
    double avgrel;
    double result;


    ae_assert(ae_round(lm->w.ptr.p_double[1], _state)==logit_logitvnum, "MNLRMSError: Incorrect MNL version!", _state);
    logit_mnlallerrors(lm, xy, npoints, &relcls, &avgce, &rms, &avg, &avgrel, _state);
    result = avg;
    return result;
}


/*************************************************************************
Average relative error on the test set

INPUT PARAMETERS:
    LM      -   logit model
    XY      -   test set
    NPoints -   test set size

RESULT:
    average relative error (error when estimating posterior probabilities).

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
double mnlavgrelerror(logitmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     ae_state *_state)
{
    double relcls;
    double avgce;
    double rms;
    double avg;
    double avgrel;
    double result;


    ae_assert(ae_round(lm->w.ptr.p_double[1], _state)==logit_logitvnum, "MNLRMSError: Incorrect MNL version!", _state);
    logit_mnlallerrors(lm, xy, ssize, &relcls, &avgce, &rms, &avg, &avgrel, _state);
    result = avgrel;
    return result;
}


/*************************************************************************
Classification error on test set = MNLRelClsError*NPoints

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
ae_int_t mnlclserror(logitmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t nvars;
    ae_int_t nclasses;
    ae_int_t i;
    ae_int_t j;
    ae_vector workx;
    ae_vector worky;
    ae_int_t nmax;
    ae_int_t result;

    ae_frame_make(_state, &_frame_block);
    ae_vector_init(&workx, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&worky, 0, DT_REAL, _state, ae_true);

    ae_assert(ae_fp_eq(lm->w.ptr.p_double[1],logit_logitvnum), "MNLClsError: unexpected model version", _state);
    nvars = ae_round(lm->w.ptr.p_double[2], _state);
    nclasses = ae_round(lm->w.ptr.p_double[3], _state);
    ae_vector_set_length(&workx, nvars-1+1, _state);
    ae_vector_set_length(&worky, nclasses-1+1, _state);
    result = 0;
    for(i=0; i<=npoints-1; i++)
    {
        
        /*
         * Process
         */
        ae_v_move(&workx.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
        mnlprocess(lm, &workx, &worky, _state);
        
        /*
         * Logit version of the answer
         */
        nmax = 0;
        for(j=0; j<=nclasses-1; j++)
        {
            if( ae_fp_greater(worky.ptr.p_double[j],worky.ptr.p_double[nmax]) )
            {
                nmax = j;
            }
        }
        
        /*
         * compare
         */
        if( nmax!=ae_round(xy->ptr.pp_double[i][nvars], _state) )
        {
            result = result+1;
        }
    }
    ae_frame_leave(_state);
    return result;
}


/*************************************************************************
Internal subroutine. Places exponents of the anti-overflow shifted
internal linear outputs into the service part of the W array.
*************************************************************************/
static void logit_mnliexp(/* Real    */ ae_vector* w,
     /* Real    */ ae_vector* x,
     ae_state *_state)
{
    ae_int_t nvars;
    ae_int_t nclasses;
    ae_int_t offs;
    ae_int_t i;
    ae_int_t i1;
    double v;
    double mx;


    ae_assert(ae_fp_eq(w->ptr.p_double[1],logit_logitvnum), "LOGIT: unexpected model version", _state);
    nvars = ae_round(w->ptr.p_double[2], _state);
    nclasses = ae_round(w->ptr.p_double[3], _state);
    offs = ae_round(w->ptr.p_double[4], _state);
    i1 = offs+(nvars+1)*(nclasses-1);
    for(i=0; i<=nclasses-2; i++)
    {
        v = ae_v_dotproduct(&w->ptr.p_double[offs+i*(nvars+1)], 1, &x->ptr.p_double[0], 1, ae_v_len(offs+i*(nvars+1),offs+i*(nvars+1)+nvars-1));
        w->ptr.p_double[i1+i] = v+w->ptr.p_double[offs+i*(nvars+1)+nvars];
    }
    w->ptr.p_double[i1+nclasses-1] = 0;
    mx = 0;
    for(i=i1; i<=i1+nclasses-1; i++)
    {
        mx = ae_maxreal(mx, w->ptr.p_double[i], _state);
    }
    for(i=i1; i<=i1+nclasses-1; i++)
    {
        w->ptr.p_double[i] = ae_exp(w->ptr.p_double[i]-mx, _state);
    }
}


/*************************************************************************
Calculation of all types of errors

  -- ALGLIB --
     Copyright 30.08.2008 by Bochkanov Sergey
*************************************************************************/
static void logit_mnlallerrors(logitmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double* relcls,
     double* avgce,
     double* rms,
     double* avg,
     double* avgrel,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t nvars;
    ae_int_t nclasses;
    ae_int_t i;
    ae_vector buf;
    ae_vector workx;
    ae_vector y;
    ae_vector dy;

    ae_frame_make(_state, &_frame_block);
    *relcls = 0;
    *avgce = 0;
    *rms = 0;
    *avg = 0;
    *avgrel = 0;
    ae_vector_init(&buf, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&workx, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&dy, 0, DT_REAL, _state, ae_true);

    ae_assert(ae_round(lm->w.ptr.p_double[1], _state)==logit_logitvnum, "MNL unit: Incorrect MNL version!", _state);
    nvars = ae_round(lm->w.ptr.p_double[2], _state);
    nclasses = ae_round(lm->w.ptr.p_double[3], _state);
    ae_vector_set_length(&workx, nvars-1+1, _state);
    ae_vector_set_length(&y, nclasses-1+1, _state);
    ae_vector_set_length(&dy, 0+1, _state);
    dserrallocate(nclasses, &buf, _state);
    for(i=0; i<=npoints-1; i++)
    {
        ae_v_move(&workx.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
        mnlprocess(lm, &workx, &y, _state);
        dy.ptr.p_double[0] = xy->ptr.pp_double[i][nvars];
        dserraccumulate(&buf, &y, &dy, _state);
    }
    dserrfinish(&buf, _state);
    *relcls = buf.ptr.p_double[0];
    *avgce = buf.ptr.p_double[1];
    *rms = buf.ptr.p_double[2];
    *avg = buf.ptr.p_double[3];
    *avgrel = buf.ptr.p_double[4];
    ae_frame_leave(_state);
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

PARAMETERS DESCRIPRION

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
static void logit_mnlmcsrch(ae_int_t n,
     /* Real    */ ae_vector* x,
     double* f,
     /* Real    */ ae_vector* g,
     /* Real    */ ae_vector* s,
     double* stp,
     ae_int_t* info,
     ae_int_t* nfev,
     /* Real    */ ae_vector* wa,
     logitmcstate* state,
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
            if( ((((((n<=0||ae_fp_less_eq(*stp,0))||ae_fp_less(logit_ftol,0))||ae_fp_less(logit_gtol,zero))||ae_fp_less(logit_xtol,zero))||ae_fp_less(logit_stpmin,zero))||ae_fp_less(logit_stpmax,logit_stpmin))||logit_maxfev<=0 )
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
            state->dgtest = logit_ftol*state->dginit;
            state->width = logit_stpmax-logit_stpmin;
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
            if( ae_fp_greater(*stp,logit_stpmax) )
            {
                *stp = logit_stpmax;
            }
            if( ae_fp_less(*stp,logit_stpmin) )
            {
                *stp = logit_stpmin;
            }
            
            /*
             *        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
             *        STP BE THE LOWEST POINT OBTAINED SO FAR.
             */
            if( (((state->brackt&&(ae_fp_less_eq(*stp,state->stmin)||ae_fp_greater_eq(*stp,state->stmax)))||*nfev>=logit_maxfev-1)||state->infoc==0)||(state->brackt&&ae_fp_less_eq(state->stmax-state->stmin,logit_xtol*state->stmax)) )
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
            if( (ae_fp_eq(*stp,logit_stpmax)&&ae_fp_less_eq(*f,state->ftest1))&&ae_fp_less_eq(state->dg,state->dgtest) )
            {
                *info = 5;
            }
            if( ae_fp_eq(*stp,logit_stpmin)&&(ae_fp_greater(*f,state->ftest1)||ae_fp_greater_eq(state->dg,state->dgtest)) )
            {
                *info = 4;
            }
            if( *nfev>=logit_maxfev )
            {
                *info = 3;
            }
            if( state->brackt&&ae_fp_less_eq(state->stmax-state->stmin,logit_xtol*state->stmax) )
            {
                *info = 2;
            }
            if( ae_fp_less_eq(*f,state->ftest1)&&ae_fp_less_eq(ae_fabs(state->dg, _state),-logit_gtol*state->dginit) )
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
            if( (state->stage1&&ae_fp_less_eq(*f,state->ftest1))&&ae_fp_greater_eq(state->dg,ae_minreal(logit_ftol, logit_gtol, _state)*state->dginit) )
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
                logit_mnlmcstep(&state->stx, &state->fxm, &state->dgxm, &state->sty, &state->fym, &state->dgym, stp, state->fm, state->dgm, &state->brackt, state->stmin, state->stmax, &state->infoc, _state);
                
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
                logit_mnlmcstep(&state->stx, &state->fx, &state->dgx, &state->sty, &state->fy, &state->dgy, stp, *f, state->dg, &state->brackt, state->stmin, state->stmax, &state->infoc, _state);
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


static void logit_mnlmcstep(double* stx,
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


ae_bool _logitmodel_init(logitmodel* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->w, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _logitmodel_init_copy(logitmodel* dst, logitmodel* src, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init_copy(&dst->w, &src->w, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _logitmodel_clear(logitmodel* p)
{
    ae_vector_clear(&p->w);
}


ae_bool _logitmcstate_init(logitmcstate* /*p*/, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    return ae_true;
}


ae_bool _logitmcstate_init_copy(logitmcstate* dst, logitmcstate* src, ae_state */*_state*/, ae_bool ) //make_automatic)
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


void _logitmcstate_clear(logitmcstate* ) //p)
{
}


ae_bool _mnlreport_init(mnlreport* /*p*/, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    return ae_true;
}


ae_bool _mnlreport_init_copy(mnlreport* dst, mnlreport* src, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    dst->ngrad = src->ngrad;
    dst->nhess = src->nhess;
    return ae_true;
}


void _mnlreport_clear(mnlreport* ) //p)
{
}




/*************************************************************************
DESCRIPTION:

This function creates MCPD (Markov Chains for Population Data) solver.

This  solver  can  be  used  to find transition matrix P for N-dimensional
prediction  problem  where transition from X[i] to X[i+1] is  modelled  as
    X[i+1] = P*X[i]
where X[i] and X[i+1] are N-dimensional population vectors (components  of
each X are non-negative), and P is a N*N transition matrix (elements of  P
are non-negative, each column sums to 1.0).

Such models arise when when:
* there is some population of individuals
* individuals can have different states
* individuals can transit from one state to another
* population size is constant, i.e. there is no new individuals and no one
  leaves population
* you want to model transitions of individuals from one state into another

USAGE:

Here we give very brief outline of the MCPD. We strongly recommend you  to
read examples in the ALGLIB Reference Manual and to read ALGLIB User Guide
on data analysis which is available at http://www.alglib.net/dataanalysis/

1. User initializes algorithm state with MCPDCreate() call

2. User  adds  one  or  more  tracks -  sequences of states which describe
   evolution of a system being modelled from different starting conditions

3. User may add optional boundary, equality  and/or  linear constraints on
   the coefficients of P by calling one of the following functions:
   * MCPDSetEC() to set equality constraints
   * MCPDSetBC() to set bound constraints
   * MCPDSetLC() to set linear constraints

4. Optionally,  user  may  set  custom  weights  for prediction errors (by
   default, algorithm assigns non-equal, automatically chosen weights  for
   errors in the prediction of different components of X). It can be  done
   with a call of MCPDSetPredictionWeights() function.

5. User calls MCPDSolve() function which takes algorithm  state and
   pointer (delegate, etc.) to callback function which calculates F/G.

6. User calls MCPDResults() to get solution

INPUT PARAMETERS:
    N       -   problem dimension, N>=1

OUTPUT PARAMETERS:
    State   -   structure stores algorithm state

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdcreate(ae_int_t n, mcpdstate* s, ae_state *_state)
{

    _mcpdstate_clear(s);

    ae_assert(n>=1, "MCPDCreate: N<1", _state);
    mcpd_mcpdinit(n, -1, -1, s, _state);
}


/*************************************************************************
DESCRIPTION:

This function is a specialized version of MCPDCreate()  function,  and  we
recommend  you  to read comments for this function for general information
about MCPD solver.

This  function  creates  MCPD (Markov Chains for Population  Data)  solver
for "Entry-state" model,  i.e. model  where transition from X[i] to X[i+1]
is modelled as
    X[i+1] = P*X[i]
where
    X[i] and X[i+1] are N-dimensional state vectors
    P is a N*N transition matrix
and  one  selected component of X[] is called "entry" state and is treated
in a special way:
    system state always transits from "entry" state to some another state
    system state can not transit from any state into "entry" state
Such conditions basically mean that row of P which corresponds to  "entry"
state is zero.

Such models arise when:
* there is some population of individuals
* individuals can have different states
* individuals can transit from one state to another
* population size is NOT constant -  at every moment of time there is some
  (unpredictable) amount of "new" individuals, which can transit into  one
  of the states at the next turn, but still no one leaves population
* you want to model transitions of individuals from one state into another
* but you do NOT want to predict amount of "new"  individuals  because  it
  does not depends on individuals already present (hence  system  can  not
  transit INTO entry state - it can only transit FROM it).

This model is discussed  in  more  details  in  the ALGLIB User Guide (see
http://www.alglib.net/dataanalysis/ for more data).

INPUT PARAMETERS:
    N       -   problem dimension, N>=2
    EntryState- index of entry state, in 0..N-1

OUTPUT PARAMETERS:
    State   -   structure stores algorithm state

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdcreateentry(ae_int_t n,
     ae_int_t entrystate,
     mcpdstate* s,
     ae_state *_state)
{

    _mcpdstate_clear(s);

    ae_assert(n>=2, "MCPDCreateEntry: N<2", _state);
    ae_assert(entrystate>=0, "MCPDCreateEntry: EntryState<0", _state);
    ae_assert(entrystate<n, "MCPDCreateEntry: EntryState>=N", _state);
    mcpd_mcpdinit(n, entrystate, -1, s, _state);
}


/*************************************************************************
DESCRIPTION:

This function is a specialized version of MCPDCreate()  function,  and  we
recommend  you  to read comments for this function for general information
about MCPD solver.

This  function  creates  MCPD (Markov Chains for Population  Data)  solver
for "Exit-state" model,  i.e. model  where  transition from X[i] to X[i+1]
is modelled as
    X[i+1] = P*X[i]
where
    X[i] and X[i+1] are N-dimensional state vectors
    P is a N*N transition matrix
and  one  selected component of X[] is called "exit"  state and is treated
in a special way:
    system state can transit from any state into "exit" state
    system state can not transit from "exit" state into any other state
    transition operator discards "exit" state (makes it zero at each turn)
Such  conditions  basically  mean  that  column  of P which corresponds to
"exit" state is zero. Multiplication by such P may decrease sum of  vector
components.

Such models arise when:
* there is some population of individuals
* individuals can have different states
* individuals can transit from one state to another
* population size is NOT constant - individuals can move into "exit" state
  and leave population at the next turn, but there are no new individuals
* amount of individuals which leave population can be predicted
* you want to model transitions of individuals from one state into another
  (including transitions into the "exit" state)

This model is discussed  in  more  details  in  the ALGLIB User Guide (see
http://www.alglib.net/dataanalysis/ for more data).

INPUT PARAMETERS:
    N       -   problem dimension, N>=2
    ExitState-  index of exit state, in 0..N-1

OUTPUT PARAMETERS:
    State   -   structure stores algorithm state

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdcreateexit(ae_int_t n,
     ae_int_t exitstate,
     mcpdstate* s,
     ae_state *_state)
{

    _mcpdstate_clear(s);

    ae_assert(n>=2, "MCPDCreateExit: N<2", _state);
    ae_assert(exitstate>=0, "MCPDCreateExit: ExitState<0", _state);
    ae_assert(exitstate<n, "MCPDCreateExit: ExitState>=N", _state);
    mcpd_mcpdinit(n, -1, exitstate, s, _state);
}


/*************************************************************************
DESCRIPTION:

This function is a specialized version of MCPDCreate()  function,  and  we
recommend  you  to read comments for this function for general information
about MCPD solver.

This  function  creates  MCPD (Markov Chains for Population  Data)  solver
for "Entry-Exit-states" model, i.e. model where  transition  from  X[i] to
X[i+1] is modelled as
    X[i+1] = P*X[i]
where
    X[i] and X[i+1] are N-dimensional state vectors
    P is a N*N transition matrix
one selected component of X[] is called "entry" state and is treated in  a
special way:
    system state always transits from "entry" state to some another state
    system state can not transit from any state into "entry" state
and another one component of X[] is called "exit" state and is treated  in
a special way too:
    system state can transit from any state into "exit" state
    system state can not transit from "exit" state into any other state
    transition operator discards "exit" state (makes it zero at each turn)
Such conditions basically mean that:
    row of P which corresponds to "entry" state is zero
    column of P which corresponds to "exit" state is zero
Multiplication by such P may decrease sum of vector components.

Such models arise when:
* there is some population of individuals
* individuals can have different states
* individuals can transit from one state to another
* population size is NOT constant
* at every moment of time there is some (unpredictable)  amount  of  "new"
  individuals, which can transit into one of the states at the next turn
* some  individuals  can  move  (predictably)  into "exit" state and leave
  population at the next turn
* you want to model transitions of individuals from one state into another,
  including transitions from the "entry" state and into the "exit" state.
* but you do NOT want to predict amount of "new"  individuals  because  it
  does not depends on individuals already present (hence  system  can  not
  transit INTO entry state - it can only transit FROM it).

This model is discussed  in  more  details  in  the ALGLIB User Guide (see
http://www.alglib.net/dataanalysis/ for more data).

INPUT PARAMETERS:
    N       -   problem dimension, N>=2
    EntryState- index of entry state, in 0..N-1
    ExitState-  index of exit state, in 0..N-1

OUTPUT PARAMETERS:
    State   -   structure stores algorithm state

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdcreateentryexit(ae_int_t n,
     ae_int_t entrystate,
     ae_int_t exitstate,
     mcpdstate* s,
     ae_state *_state)
{

    _mcpdstate_clear(s);

    ae_assert(n>=2, "MCPDCreateEntryExit: N<2", _state);
    ae_assert(entrystate>=0, "MCPDCreateEntryExit: EntryState<0", _state);
    ae_assert(entrystate<n, "MCPDCreateEntryExit: EntryState>=N", _state);
    ae_assert(exitstate>=0, "MCPDCreateEntryExit: ExitState<0", _state);
    ae_assert(exitstate<n, "MCPDCreateEntryExit: ExitState>=N", _state);
    ae_assert(entrystate!=exitstate, "MCPDCreateEntryExit: EntryState=ExitState", _state);
    mcpd_mcpdinit(n, entrystate, exitstate, s, _state);
}


/*************************************************************************
This  function  is  used to add a track - sequence of system states at the
different moments of its evolution.

You  may  add  one  or several tracks to the MCPD solver. In case you have
several tracks, they won't overwrite each other. For example,  if you pass
two tracks, A1-A2-A3 (system at t=A+1, t=A+2 and t=A+3) and B1-B2-B3, then
solver will try to model transitions from t=A+1 to t=A+2, t=A+2 to  t=A+3,
t=B+1 to t=B+2, t=B+2 to t=B+3. But it WONT mix these two tracks - i.e. it
wont try to model transition from t=A+3 to t=B+1.

INPUT PARAMETERS:
    S       -   solver
    XY      -   track, array[K,N]:
                * I-th row is a state at t=I
                * elements of XY must be non-negative (exception will be
                  thrown on negative elements)
    K       -   number of points in a track
                * if given, only leading K rows of XY are used
                * if not given, automatically determined from size of XY

NOTES:

1. Track may contain either proportional or population data:
   * with proportional data all rows of XY must sum to 1.0, i.e. we have
     proportions instead of absolute population values
   * with population data rows of XY contain population counts and generally
     do not sum to 1.0 (although they still must be non-negative)

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdaddtrack(mcpdstate* s,
     /* Real    */ ae_matrix* xy,
     ae_int_t k,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t n;
    double s0;
    double s1;


    n = s->n;
    ae_assert(k>=0, "MCPDAddTrack: K<0", _state);
    ae_assert(xy->cols>=n, "MCPDAddTrack: Cols(XY)<N", _state);
    ae_assert(xy->rows>=k, "MCPDAddTrack: Rows(XY)<K", _state);
    ae_assert(apservisfinitematrix(xy, k, n, _state), "MCPDAddTrack: XY contains infinite or NaN elements", _state);
    for(i=0; i<=k-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            ae_assert(ae_fp_greater_eq(xy->ptr.pp_double[i][j],0), "MCPDAddTrack: XY contains negative elements", _state);
        }
    }
    if( k<2 )
    {
        return;
    }
    if( s->data.rows<s->npairs+k-1 )
    {
        rmatrixresize(&s->data, ae_maxint(2*s->data.rows, s->npairs+k-1, _state), 2*n, _state);
    }
    for(i=0; i<=k-2; i++)
    {
        s0 = 0;
        s1 = 0;
        for(j=0; j<=n-1; j++)
        {
            if( s->states.ptr.p_int[j]>=0 )
            {
                s0 = s0+xy->ptr.pp_double[i][j];
            }
            if( s->states.ptr.p_int[j]<=0 )
            {
                s1 = s1+xy->ptr.pp_double[i+1][j];
            }
        }
        if( ae_fp_greater(s0,0)&&ae_fp_greater(s1,0) )
        {
            for(j=0; j<=n-1; j++)
            {
                if( s->states.ptr.p_int[j]>=0 )
                {
                    s->data.ptr.pp_double[s->npairs][j] = xy->ptr.pp_double[i][j]/s0;
                }
                else
                {
                    s->data.ptr.pp_double[s->npairs][j] = 0.0;
                }
                if( s->states.ptr.p_int[j]<=0 )
                {
                    s->data.ptr.pp_double[s->npairs][n+j] = xy->ptr.pp_double[i+1][j]/s1;
                }
                else
                {
                    s->data.ptr.pp_double[s->npairs][n+j] = 0.0;
                }
            }
            s->npairs = s->npairs+1;
        }
    }
}


/*************************************************************************
This function is used to add equality constraints on the elements  of  the
transition matrix P.

MCPD solver has four types of constraints which can be placed on P:
* user-specified equality constraints (optional)
* user-specified bound constraints (optional)
* user-specified general linear constraints (optional)
* basic constraints (always present):
  * non-negativity: P[i,j]>=0
  * consistency: every column of P sums to 1.0

Final  constraints  which  are  passed  to  the  underlying  optimizer are
calculated  as  intersection  of all present constraints. For example, you
may specify boundary constraint on P[0,0] and equality one:
    0.1<=P[0,0]<=0.9
    P[0,0]=0.5
Such  combination  of  constraints  will  be  silently  reduced  to  their
intersection, which is P[0,0]=0.5.

This  function  can  be  used  to  place equality constraints on arbitrary
subset of elements of P. Set of constraints is specified by EC, which  may
contain either NAN's or finite numbers from [0,1]. NAN denotes absence  of
constraint, finite number denotes equality constraint on specific  element
of P.

You can also  use  MCPDAddEC()  function  which  allows  to  ADD  equality
constraint  for  one  element  of P without changing constraints for other
elements.

These functions (MCPDSetEC and MCPDAddEC) interact as follows:
* there is internal matrix of equality constraints which is stored in  the
  MCPD solver
* MCPDSetEC() replaces this matrix by another one (SET)
* MCPDAddEC() modifies one element of this matrix and  leaves  other  ones
  unchanged (ADD)
* thus  MCPDAddEC()  call  preserves  all  modifications  done by previous
  calls,  while  MCPDSetEC()  completely discards all changes  done to the
  equality constraints.

INPUT PARAMETERS:
    S       -   solver
    EC      -   equality constraints, array[N,N]. Elements of  EC  can  be
                either NAN's or finite  numbers from  [0,1].  NAN  denotes
                absence  of  constraints,  while  finite  value    denotes
                equality constraint on the corresponding element of P.

NOTES:

1. infinite values of EC will lead to exception being thrown. Values  less
than 0.0 or greater than 1.0 will lead to error code being returned  after
call to MCPDSolve().

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsetec(mcpdstate* s,
     /* Real    */ ae_matrix* ec,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t n;


    n = s->n;
    ae_assert(ec->cols>=n, "MCPDSetEC: Cols(EC)<N", _state);
    ae_assert(ec->rows>=n, "MCPDSetEC: Rows(EC)<N", _state);
    for(i=0; i<=n-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            ae_assert(ae_isfinite(ec->ptr.pp_double[i][j], _state)||ae_isnan(ec->ptr.pp_double[i][j], _state), "MCPDSetEC: EC containts infinite elements", _state);
            s->ec.ptr.pp_double[i][j] = ec->ptr.pp_double[i][j];
        }
    }
}


/*************************************************************************
This function is used to add equality constraints on the elements  of  the
transition matrix P.

MCPD solver has four types of constraints which can be placed on P:
* user-specified equality constraints (optional)
* user-specified bound constraints (optional)
* user-specified general linear constraints (optional)
* basic constraints (always present):
  * non-negativity: P[i,j]>=0
  * consistency: every column of P sums to 1.0

Final  constraints  which  are  passed  to  the  underlying  optimizer are
calculated  as  intersection  of all present constraints. For example, you
may specify boundary constraint on P[0,0] and equality one:
    0.1<=P[0,0]<=0.9
    P[0,0]=0.5
Such  combination  of  constraints  will  be  silently  reduced  to  their
intersection, which is P[0,0]=0.5.

This function can be used to ADD equality constraint for one element of  P
without changing constraints for other elements.

You  can  also  use  MCPDSetEC()  function  which  allows  you  to specify
arbitrary set of equality constraints in one call.

These functions (MCPDSetEC and MCPDAddEC) interact as follows:
* there is internal matrix of equality constraints which is stored in the
  MCPD solver
* MCPDSetEC() replaces this matrix by another one (SET)
* MCPDAddEC() modifies one element of this matrix and leaves  other  ones
  unchanged (ADD)
* thus  MCPDAddEC()  call  preserves  all  modifications done by previous
  calls,  while  MCPDSetEC()  completely discards all changes done to the
  equality constraints.

INPUT PARAMETERS:
    S       -   solver
    I       -   row index of element being constrained
    J       -   column index of element being constrained
    C       -   value (constraint for P[I,J]).  Can  be  either  NAN  (no
                constraint) or finite value from [0,1].
                
NOTES:

1. infinite values of C  will lead to exception being thrown. Values  less
than 0.0 or greater than 1.0 will lead to error code being returned  after
call to MCPDSolve().

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdaddec(mcpdstate* s,
     ae_int_t i,
     ae_int_t j,
     double c,
     ae_state *_state)
{


    ae_assert(i>=0, "MCPDAddEC: I<0", _state);
    ae_assert(i<s->n, "MCPDAddEC: I>=N", _state);
    ae_assert(j>=0, "MCPDAddEC: J<0", _state);
    ae_assert(j<s->n, "MCPDAddEC: J>=N", _state);
    ae_assert(ae_isnan(c, _state)||ae_isfinite(c, _state), "MCPDAddEC: C is not finite number or NAN", _state);
    s->ec.ptr.pp_double[i][j] = c;
}


/*************************************************************************
This function is used to add bound constraints  on  the  elements  of  the
transition matrix P.

MCPD solver has four types of constraints which can be placed on P:
* user-specified equality constraints (optional)
* user-specified bound constraints (optional)
* user-specified general linear constraints (optional)
* basic constraints (always present):
  * non-negativity: P[i,j]>=0
  * consistency: every column of P sums to 1.0

Final  constraints  which  are  passed  to  the  underlying  optimizer are
calculated  as  intersection  of all present constraints. For example, you
may specify boundary constraint on P[0,0] and equality one:
    0.1<=P[0,0]<=0.9
    P[0,0]=0.5
Such  combination  of  constraints  will  be  silently  reduced  to  their
intersection, which is P[0,0]=0.5.

This  function  can  be  used  to  place bound   constraints  on arbitrary
subset  of  elements  of  P.  Set of constraints is specified by BndL/BndU
matrices, which may contain arbitrary combination  of  finite  numbers  or
infinities (like -INF<x<=0.5 or 0.1<=x<+INF).

You can also use MCPDAddBC() function which allows to ADD bound constraint
for one element of P without changing constraints for other elements.

These functions (MCPDSetBC and MCPDAddBC) interact as follows:
* there is internal matrix of bound constraints which is stored in the
  MCPD solver
* MCPDSetBC() replaces this matrix by another one (SET)
* MCPDAddBC() modifies one element of this matrix and  leaves  other  ones
  unchanged (ADD)
* thus  MCPDAddBC()  call  preserves  all  modifications  done by previous
  calls,  while  MCPDSetBC()  completely discards all changes  done to the
  equality constraints.

INPUT PARAMETERS:
    S       -   solver
    BndL    -   lower bounds constraints, array[N,N]. Elements of BndL can
                be finite numbers or -INF.
    BndU    -   upper bounds constraints, array[N,N]. Elements of BndU can
                be finite numbers or +INF.

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsetbc(mcpdstate* s,
     /* Real    */ ae_matrix* bndl,
     /* Real    */ ae_matrix* bndu,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t n;


    n = s->n;
    ae_assert(bndl->cols>=n, "MCPDSetBC: Cols(BndL)<N", _state);
    ae_assert(bndl->rows>=n, "MCPDSetBC: Rows(BndL)<N", _state);
    ae_assert(bndu->cols>=n, "MCPDSetBC: Cols(BndU)<N", _state);
    ae_assert(bndu->rows>=n, "MCPDSetBC: Rows(BndU)<N", _state);
    for(i=0; i<=n-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            ae_assert(ae_isfinite(bndl->ptr.pp_double[i][j], _state)||ae_isneginf(bndl->ptr.pp_double[i][j], _state), "MCPDSetBC: BndL containts NAN or +INF", _state);
            ae_assert(ae_isfinite(bndu->ptr.pp_double[i][j], _state)||ae_isposinf(bndu->ptr.pp_double[i][j], _state), "MCPDSetBC: BndU containts NAN or -INF", _state);
            s->bndl.ptr.pp_double[i][j] = bndl->ptr.pp_double[i][j];
            s->bndu.ptr.pp_double[i][j] = bndu->ptr.pp_double[i][j];
        }
    }
}


/*************************************************************************
This function is used to add bound constraints  on  the  elements  of  the
transition matrix P.

MCPD solver has four types of constraints which can be placed on P:
* user-specified equality constraints (optional)
* user-specified bound constraints (optional)
* user-specified general linear constraints (optional)
* basic constraints (always present):
  * non-negativity: P[i,j]>=0
  * consistency: every column of P sums to 1.0

Final  constraints  which  are  passed  to  the  underlying  optimizer are
calculated  as  intersection  of all present constraints. For example, you
may specify boundary constraint on P[0,0] and equality one:
    0.1<=P[0,0]<=0.9
    P[0,0]=0.5
Such  combination  of  constraints  will  be  silently  reduced  to  their
intersection, which is P[0,0]=0.5.

This  function  can  be  used to ADD bound constraint for one element of P
without changing constraints for other elements.

You  can  also  use  MCPDSetBC()  function  which  allows to  place  bound
constraints  on arbitrary subset of elements of P.   Set of constraints is
specified  by  BndL/BndU matrices, which may contain arbitrary combination
of finite numbers or infinities (like -INF<x<=0.5 or 0.1<=x<+INF).

These functions (MCPDSetBC and MCPDAddBC) interact as follows:
* there is internal matrix of bound constraints which is stored in the
  MCPD solver
* MCPDSetBC() replaces this matrix by another one (SET)
* MCPDAddBC() modifies one element of this matrix and  leaves  other  ones
  unchanged (ADD)
* thus  MCPDAddBC()  call  preserves  all  modifications  done by previous
  calls,  while  MCPDSetBC()  completely discards all changes  done to the
  equality constraints.

INPUT PARAMETERS:
    S       -   solver
    I       -   row index of element being constrained
    J       -   column index of element being constrained
    BndL    -   lower bound
    BndU    -   upper bound

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdaddbc(mcpdstate* s,
     ae_int_t i,
     ae_int_t j,
     double bndl,
     double bndu,
     ae_state *_state)
{


    ae_assert(i>=0, "MCPDAddBC: I<0", _state);
    ae_assert(i<s->n, "MCPDAddBC: I>=N", _state);
    ae_assert(j>=0, "MCPDAddBC: J<0", _state);
    ae_assert(j<s->n, "MCPDAddBC: J>=N", _state);
    ae_assert(ae_isfinite(bndl, _state)||ae_isneginf(bndl, _state), "MCPDAddBC: BndL is NAN or +INF", _state);
    ae_assert(ae_isfinite(bndu, _state)||ae_isposinf(bndu, _state), "MCPDAddBC: BndU is NAN or -INF", _state);
    s->bndl.ptr.pp_double[i][j] = bndl;
    s->bndu.ptr.pp_double[i][j] = bndu;
}


/*************************************************************************
This function is used to set linear equality/inequality constraints on the
elements of the transition matrix P.

This function can be used to set one or several general linear constraints
on the elements of P. Two types of constraints are supported:
* equality constraints
* inequality constraints (both less-or-equal and greater-or-equal)

Coefficients  of  constraints  are  specified  by  matrix  C (one  of  the
parameters).  One  row  of  C  corresponds  to  one  constraint.   Because
transition  matrix P has N*N elements,  we  need  N*N columns to store all
coefficients  (they  are  stored row by row), and one more column to store
right part - hence C has N*N+1 columns.  Constraint  kind is stored in the
CT array.

Thus, I-th linear constraint is
    P[0,0]*C[I,0] + P[0,1]*C[I,1] + .. + P[0,N-1]*C[I,N-1] +
        + P[1,0]*C[I,N] + P[1,1]*C[I,N+1] + ... +
        + P[N-1,N-1]*C[I,N*N-1]  ?=?  C[I,N*N]
where ?=? can be either "=" (CT[i]=0), "<=" (CT[i]<0) or ">=" (CT[i]>0).

Your constraint may involve only some subset of P (less than N*N elements).
For example it can be something like
    P[0,0] + P[0,1] = 0.5
In this case you still should pass matrix  with N*N+1 columns, but all its
elements (except for C[0,0], C[0,1] and C[0,N*N-1]) will be zero.

INPUT PARAMETERS:
    S       -   solver
    C       -   array[K,N*N+1] - coefficients of constraints
                (see above for complete description)
    CT      -   array[K] - constraint types
                (see above for complete description)
    K       -   number of equality/inequality constraints, K>=0:
                * if given, only leading K elements of C/CT are used
                * if not given, automatically determined from sizes of C/CT

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsetlc(mcpdstate* s,
     /* Real    */ ae_matrix* c,
     /* Integer */ ae_vector* ct,
     ae_int_t k,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t n;


    n = s->n;
    ae_assert(c->cols>=n*n+1, "MCPDSetLC: Cols(C)<N*N+1", _state);
    ae_assert(c->rows>=k, "MCPDSetLC: Rows(C)<K", _state);
    ae_assert(ct->cnt>=k, "MCPDSetLC: Len(CT)<K", _state);
    ae_assert(apservisfinitematrix(c, k, n*n+1, _state), "MCPDSetLC: C contains infinite or NaN values!", _state);
    rmatrixsetlengthatleast(&s->c, k, n*n+1, _state);
    ivectorsetlengthatleast(&s->ct, k, _state);
    for(i=0; i<=k-1; i++)
    {
        for(j=0; j<=n*n; j++)
        {
            s->c.ptr.pp_double[i][j] = c->ptr.pp_double[i][j];
        }
        s->ct.ptr.p_int[i] = ct->ptr.p_int[i];
    }
    s->ccnt = k;
}


/*************************************************************************
This function allows to  tune  amount  of  Tikhonov  regularization  being
applied to your problem.

By default, regularizing term is equal to r*||P-prior_P||^2, where r is  a
small non-zero value,  P is transition matrix, prior_P is identity matrix,
||X||^2 is a sum of squared elements of X.

This  function  allows  you to change coefficient r. You can  also  change
prior values with MCPDSetPrior() function.

INPUT PARAMETERS:
    S       -   solver
    V       -   regularization  coefficient, finite non-negative value. It
                is  not  recommended  to specify zero value unless you are
                pretty sure that you want it.

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsettikhonovregularizer(mcpdstate* s, double v, ae_state *_state)
{


    ae_assert(ae_isfinite(v, _state), "MCPDSetTikhonovRegularizer: V is infinite or NAN", _state);
    ae_assert(ae_fp_greater_eq(v,0.0), "MCPDSetTikhonovRegularizer: V is less than zero", _state);
    s->regterm = v;
}


/*************************************************************************
This  function  allows to set prior values used for regularization of your
problem.

By default, regularizing term is equal to r*||P-prior_P||^2, where r is  a
small non-zero value,  P is transition matrix, prior_P is identity matrix,
||X||^2 is a sum of squared elements of X.

This  function  allows  you to change prior values prior_P. You  can  also
change r with MCPDSetTikhonovRegularizer() function.

INPUT PARAMETERS:
    S       -   solver
    PP      -   array[N,N], matrix of prior values:
                1. elements must be real numbers from [0,1]
                2. columns must sum to 1.0.
                First property is checked (exception is thrown otherwise),
                while second one is not checked/enforced.

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsetprior(mcpdstate* s,
     /* Real    */ ae_matrix* pp,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix _pp;
    ae_int_t i;
    ae_int_t j;
    ae_int_t n;

    ae_frame_make(_state, &_frame_block);
    ae_matrix_init_copy(&_pp, pp, _state, ae_true);
    pp = &_pp;

    n = s->n;
    ae_assert(pp->cols>=n, "MCPDSetPrior: Cols(PP)<N", _state);
    ae_assert(pp->rows>=n, "MCPDSetPrior: Rows(PP)<K", _state);
    for(i=0; i<=n-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            ae_assert(ae_isfinite(pp->ptr.pp_double[i][j], _state), "MCPDSetPrior: PP containts infinite elements", _state);
            ae_assert(ae_fp_greater_eq(pp->ptr.pp_double[i][j],0.0)&&ae_fp_less_eq(pp->ptr.pp_double[i][j],1.0), "MCPDSetPrior: PP[i,j] is less than 0.0 or greater than 1.0", _state);
            s->priorp.ptr.pp_double[i][j] = pp->ptr.pp_double[i][j];
        }
    }
    ae_frame_leave(_state);
}


/*************************************************************************
This function is used to change prediction weights

MCPD solver scales prediction errors as follows
    Error(P) = ||W*(y-P*x)||^2
where
    x is a system state at time t
    y is a system state at time t+1
    P is a transition matrix
    W is a diagonal scaling matrix

By default, weights are chosen in order  to  minimize  relative prediction
error instead of absolute one. For example, if one component of  state  is
about 0.5 in magnitude and another one is about 0.05, then algorithm  will
make corresponding weights equal to 2.0 and 20.0.

INPUT PARAMETERS:
    S       -   solver
    PW      -   array[N], weights:
                * must be non-negative values (exception will be thrown otherwise)
                * zero values will be replaced by automatically chosen values

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsetpredictionweights(mcpdstate* s,
     /* Real    */ ae_vector* pw,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t n;


    n = s->n;
    ae_assert(pw->cnt>=n, "MCPDSetPredictionWeights: Length(PW)<N", _state);
    for(i=0; i<=n-1; i++)
    {
        ae_assert(ae_isfinite(pw->ptr.p_double[i], _state), "MCPDSetPredictionWeights: PW containts infinite or NAN elements", _state);
        ae_assert(ae_fp_greater_eq(pw->ptr.p_double[i],0), "MCPDSetPredictionWeights: PW containts negative elements", _state);
        s->pw.ptr.p_double[i] = pw->ptr.p_double[i];
    }
}


/*************************************************************************
This function is used to start solution of the MCPD problem.

After return from this function, you can use MCPDResults() to get solution
and completion code.

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsolve(mcpdstate* s, ae_state *_state)
{
    ae_int_t n;
    ae_int_t npairs;
    ae_int_t ccnt;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    ae_int_t k2;
    double v;
    double vv;


    n = s->n;
    npairs = s->npairs;
    
    /*
     * init fields of S
     */
    s->repterminationtype = 0;
    s->repinneriterationscount = 0;
    s->repouteriterationscount = 0;
    s->repnfev = 0;
    for(k=0; k<=n-1; k++)
    {
        for(k2=0; k2<=n-1; k2++)
        {
            s->p.ptr.pp_double[k][k2] = _state->v_nan;
        }
    }
    
    /*
     * Generate "effective" weights for prediction and calculate preconditioner
     */
    for(i=0; i<=n-1; i++)
    {
        if( ae_fp_eq(s->pw.ptr.p_double[i],0) )
        {
            v = 0;
            k = 0;
            for(j=0; j<=npairs-1; j++)
            {
                if( ae_fp_neq(s->data.ptr.pp_double[j][n+i],0) )
                {
                    v = v+s->data.ptr.pp_double[j][n+i];
                    k = k+1;
                }
            }
            if( k!=0 )
            {
                s->effectivew.ptr.p_double[i] = k/v;
            }
            else
            {
                s->effectivew.ptr.p_double[i] = 1.0;
            }
        }
        else
        {
            s->effectivew.ptr.p_double[i] = s->pw.ptr.p_double[i];
        }
    }
    for(i=0; i<=n-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            s->h.ptr.p_double[i*n+j] = 2*s->regterm;
        }
    }
    for(k=0; k<=npairs-1; k++)
    {
        for(i=0; i<=n-1; i++)
        {
            for(j=0; j<=n-1; j++)
            {
                s->h.ptr.p_double[i*n+j] = s->h.ptr.p_double[i*n+j]+2*ae_sqr(s->effectivew.ptr.p_double[i], _state)*ae_sqr(s->data.ptr.pp_double[k][j], _state);
            }
        }
    }
    for(i=0; i<=n-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            if( ae_fp_eq(s->h.ptr.p_double[i*n+j],0) )
            {
                s->h.ptr.p_double[i*n+j] = 1;
            }
        }
    }
    
    /*
     * Generate "effective" BndL/BndU
     */
    for(i=0; i<=n-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            
            /*
             * Set default boundary constraints.
             * Lower bound is always zero, upper bound is calculated
             * with respect to entry/exit states.
             */
            s->effectivebndl.ptr.p_double[i*n+j] = 0.0;
            if( s->states.ptr.p_int[i]>0||s->states.ptr.p_int[j]<0 )
            {
                s->effectivebndu.ptr.p_double[i*n+j] = 0.0;
            }
            else
            {
                s->effectivebndu.ptr.p_double[i*n+j] = 1.0;
            }
            
            /*
             * Calculate intersection of the default and user-specified bound constraints.
             * This code checks consistency of such combination.
             */
            if( ae_isfinite(s->bndl.ptr.pp_double[i][j], _state)&&ae_fp_greater(s->bndl.ptr.pp_double[i][j],s->effectivebndl.ptr.p_double[i*n+j]) )
            {
                s->effectivebndl.ptr.p_double[i*n+j] = s->bndl.ptr.pp_double[i][j];
            }
            if( ae_isfinite(s->bndu.ptr.pp_double[i][j], _state)&&ae_fp_less(s->bndu.ptr.pp_double[i][j],s->effectivebndu.ptr.p_double[i*n+j]) )
            {
                s->effectivebndu.ptr.p_double[i*n+j] = s->bndu.ptr.pp_double[i][j];
            }
            if( ae_fp_greater(s->effectivebndl.ptr.p_double[i*n+j],s->effectivebndu.ptr.p_double[i*n+j]) )
            {
                s->repterminationtype = -3;
                return;
            }
            
            /*
             * Calculate intersection of the effective bound constraints
             * and user-specified equality constraints.
             * This code checks consistency of such combination.
             */
            if( ae_isfinite(s->ec.ptr.pp_double[i][j], _state) )
            {
                if( ae_fp_less(s->ec.ptr.pp_double[i][j],s->effectivebndl.ptr.p_double[i*n+j])||ae_fp_greater(s->ec.ptr.pp_double[i][j],s->effectivebndu.ptr.p_double[i*n+j]) )
                {
                    s->repterminationtype = -3;
                    return;
                }
                s->effectivebndl.ptr.p_double[i*n+j] = s->ec.ptr.pp_double[i][j];
                s->effectivebndu.ptr.p_double[i*n+j] = s->ec.ptr.pp_double[i][j];
            }
        }
    }
    
    /*
     * Generate linear constraints:
     * * "default" sums-to-one constraints (not generated for "exit" states)
     */
    rmatrixsetlengthatleast(&s->effectivec, s->ccnt+n, n*n+1, _state);
    ivectorsetlengthatleast(&s->effectivect, s->ccnt+n, _state);
    ccnt = s->ccnt;
    for(i=0; i<=s->ccnt-1; i++)
    {
        for(j=0; j<=n*n; j++)
        {
            s->effectivec.ptr.pp_double[i][j] = s->c.ptr.pp_double[i][j];
        }
        s->effectivect.ptr.p_int[i] = s->ct.ptr.p_int[i];
    }
    for(i=0; i<=n-1; i++)
    {
        if( s->states.ptr.p_int[i]>=0 )
        {
            for(k=0; k<=n*n-1; k++)
            {
                s->effectivec.ptr.pp_double[ccnt][k] = 0;
            }
            for(k=0; k<=n-1; k++)
            {
                s->effectivec.ptr.pp_double[ccnt][k*n+i] = 1;
            }
            s->effectivec.ptr.pp_double[ccnt][n*n] = 1.0;
            s->effectivect.ptr.p_int[ccnt] = 0;
            ccnt = ccnt+1;
        }
    }
    
    /*
     * create optimizer
     */
    for(i=0; i<=n-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            s->tmpp.ptr.p_double[i*n+j] = (double)1/(double)n;
        }
    }
    minbleicrestartfrom(&s->bs, &s->tmpp, _state);
    minbleicsetbc(&s->bs, &s->effectivebndl, &s->effectivebndu, _state);
    minbleicsetlc(&s->bs, &s->effectivec, &s->effectivect, ccnt, _state);
    minbleicsetinnercond(&s->bs, 0, 0, mcpd_xtol, _state);
    minbleicsetoutercond(&s->bs, mcpd_xtol, 1.0E-5, _state);
    minbleicsetprecdiag(&s->bs, &s->h, _state);
    
    /*
     * solve problem
     */
    while(minbleiciteration(&s->bs, _state))
    {
        ae_assert(s->bs.needfg, "MCPDSolve: internal error", _state);
        if( s->bs.needfg )
        {
            
            /*
             * Calculate regularization term
             */
            s->bs.f = 0.0;
            vv = s->regterm;
            for(i=0; i<=n-1; i++)
            {
                for(j=0; j<=n-1; j++)
                {
                    s->bs.f = s->bs.f+vv*ae_sqr(s->bs.x.ptr.p_double[i*n+j]-s->priorp.ptr.pp_double[i][j], _state);
                    s->bs.g.ptr.p_double[i*n+j] = 2*vv*(s->bs.x.ptr.p_double[i*n+j]-s->priorp.ptr.pp_double[i][j]);
                }
            }
            
            /*
             * calculate prediction error/gradient for K-th pair
             */
            for(k=0; k<=npairs-1; k++)
            {
                for(i=0; i<=n-1; i++)
                {
                    v = ae_v_dotproduct(&s->bs.x.ptr.p_double[i*n], 1, &s->data.ptr.pp_double[k][0], 1, ae_v_len(i*n,i*n+n-1));
                    vv = s->effectivew.ptr.p_double[i];
                    s->bs.f = s->bs.f+ae_sqr(vv*(v-s->data.ptr.pp_double[k][n+i]), _state);
                    for(j=0; j<=n-1; j++)
                    {
                        s->bs.g.ptr.p_double[i*n+j] = s->bs.g.ptr.p_double[i*n+j]+2*vv*vv*(v-s->data.ptr.pp_double[k][n+i])*s->data.ptr.pp_double[k][j];
                    }
                }
            }
            
            /*
             * continue
             */
            continue;
        }
    }
    minbleicresultsbuf(&s->bs, &s->tmpp, &s->br, _state);
    for(i=0; i<=n-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            s->p.ptr.pp_double[i][j] = s->tmpp.ptr.p_double[i*n+j];
        }
    }
    s->repterminationtype = s->br.terminationtype;
    s->repinneriterationscount = s->br.inneriterationscount;
    s->repouteriterationscount = s->br.outeriterationscount;
    s->repnfev = s->br.nfev;
}


/*************************************************************************
MCPD results

INPUT PARAMETERS:
    State   -   algorithm state

OUTPUT PARAMETERS:
    P       -   array[N,N], transition matrix
    Rep     -   optimization report. You should check Rep.TerminationType
                in  order  to  distinguish  successful  termination  from
                unsuccessful one. Speaking short, positive values  denote
                success, negative ones are failures.
                More information about fields of this  structure  can  be
                found in the comments on MCPDReport datatype.


  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdresults(mcpdstate* s,
     /* Real    */ ae_matrix* p,
     mcpdreport* rep,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;

    ae_matrix_clear(p);
    _mcpdreport_clear(rep);

    ae_matrix_set_length(p, s->n, s->n, _state);
    for(i=0; i<=s->n-1; i++)
    {
        for(j=0; j<=s->n-1; j++)
        {
            p->ptr.pp_double[i][j] = s->p.ptr.pp_double[i][j];
        }
    }
    rep->terminationtype = s->repterminationtype;
    rep->inneriterationscount = s->repinneriterationscount;
    rep->outeriterationscount = s->repouteriterationscount;
    rep->nfev = s->repnfev;
}


/*************************************************************************
Internal initialization function

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
static void mcpd_mcpdinit(ae_int_t n,
     ae_int_t entrystate,
     ae_int_t exitstate,
     mcpdstate* s,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;


    ae_assert(n>=1, "MCPDCreate: N<1", _state);
    s->n = n;
    ae_vector_set_length(&s->states, n, _state);
    for(i=0; i<=n-1; i++)
    {
        s->states.ptr.p_int[i] = 0;
    }
    if( entrystate>=0 )
    {
        s->states.ptr.p_int[entrystate] = 1;
    }
    if( exitstate>=0 )
    {
        s->states.ptr.p_int[exitstate] = -1;
    }
    s->npairs = 0;
    s->regterm = 1.0E-8;
    s->ccnt = 0;
    ae_matrix_set_length(&s->p, n, n, _state);
    ae_matrix_set_length(&s->ec, n, n, _state);
    ae_matrix_set_length(&s->bndl, n, n, _state);
    ae_matrix_set_length(&s->bndu, n, n, _state);
    ae_vector_set_length(&s->pw, n, _state);
    ae_matrix_set_length(&s->priorp, n, n, _state);
    ae_vector_set_length(&s->tmpp, n*n, _state);
    ae_vector_set_length(&s->effectivew, n, _state);
    ae_vector_set_length(&s->effectivebndl, n*n, _state);
    ae_vector_set_length(&s->effectivebndu, n*n, _state);
    ae_vector_set_length(&s->h, n*n, _state);
    for(i=0; i<=n-1; i++)
    {
        for(j=0; j<=n-1; j++)
        {
            s->p.ptr.pp_double[i][j] = 0.0;
            s->priorp.ptr.pp_double[i][j] = 0.0;
            s->bndl.ptr.pp_double[i][j] = _state->v_neginf;
            s->bndu.ptr.pp_double[i][j] = _state->v_posinf;
            s->ec.ptr.pp_double[i][j] = _state->v_nan;
        }
        s->pw.ptr.p_double[i] = 0.0;
        s->priorp.ptr.pp_double[i][i] = 1.0;
    }
    ae_matrix_set_length(&s->data, 1, 2*n, _state);
    for(i=0; i<=2*n-1; i++)
    {
        s->data.ptr.pp_double[0][i] = 0.0;
    }
    for(i=0; i<=n*n-1; i++)
    {
        s->tmpp.ptr.p_double[i] = 0.0;
    }
    minbleiccreate(n*n, &s->tmpp, &s->bs, _state);
}


ae_bool _mcpdstate_init(mcpdstate* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->states, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->data, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->ec, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->bndl, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->bndu, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->c, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->ct, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->pw, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->priorp, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !_minbleicstate_init(&p->bs, _state, make_automatic) )
        return ae_false;
    if( !_minbleicreport_init(&p->br, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->tmpp, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->effectivew, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->effectivebndl, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->effectivebndu, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->effectivec, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->effectivect, 0, DT_INT, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->h, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init(&p->p, 0, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _mcpdstate_init_copy(mcpdstate* dst, mcpdstate* src, ae_state *_state, ae_bool make_automatic)
{
    dst->n = src->n;
    if( !ae_vector_init_copy(&dst->states, &src->states, _state, make_automatic) )
        return ae_false;
    dst->npairs = src->npairs;
    if( !ae_matrix_init_copy(&dst->data, &src->data, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->ec, &src->ec, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->bndl, &src->bndl, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->bndu, &src->bndu, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->c, &src->c, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->ct, &src->ct, _state, make_automatic) )
        return ae_false;
    dst->ccnt = src->ccnt;
    if( !ae_vector_init_copy(&dst->pw, &src->pw, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->priorp, &src->priorp, _state, make_automatic) )
        return ae_false;
    dst->regterm = src->regterm;
    if( !_minbleicstate_init_copy(&dst->bs, &src->bs, _state, make_automatic) )
        return ae_false;
    dst->repinneriterationscount = src->repinneriterationscount;
    dst->repouteriterationscount = src->repouteriterationscount;
    dst->repnfev = src->repnfev;
    dst->repterminationtype = src->repterminationtype;
    if( !_minbleicreport_init_copy(&dst->br, &src->br, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->tmpp, &src->tmpp, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->effectivew, &src->effectivew, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->effectivebndl, &src->effectivebndl, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->effectivebndu, &src->effectivebndu, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->effectivec, &src->effectivec, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->effectivect, &src->effectivect, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->h, &src->h, _state, make_automatic) )
        return ae_false;
    if( !ae_matrix_init_copy(&dst->p, &src->p, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _mcpdstate_clear(mcpdstate* p)
{
    ae_vector_clear(&p->states);
    ae_matrix_clear(&p->data);
    ae_matrix_clear(&p->ec);
    ae_matrix_clear(&p->bndl);
    ae_matrix_clear(&p->bndu);
    ae_matrix_clear(&p->c);
    ae_vector_clear(&p->ct);
    ae_vector_clear(&p->pw);
    ae_matrix_clear(&p->priorp);
    _minbleicstate_clear(&p->bs);
    _minbleicreport_clear(&p->br);
    ae_vector_clear(&p->tmpp);
    ae_vector_clear(&p->effectivew);
    ae_vector_clear(&p->effectivebndl);
    ae_vector_clear(&p->effectivebndu);
    ae_matrix_clear(&p->effectivec);
    ae_vector_clear(&p->effectivect);
    ae_vector_clear(&p->h);
    ae_matrix_clear(&p->p);
}


ae_bool _mcpdreport_init(mcpdreport* /*p*/, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    return ae_true;
}


ae_bool _mcpdreport_init_copy(mcpdreport* dst, mcpdreport* src, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    dst->inneriterationscount = src->inneriterationscount;
    dst->outeriterationscount = src->outeriterationscount;
    dst->nfev = src->nfev;
    dst->terminationtype = src->terminationtype;
    return ae_true;
}


void _mcpdreport_clear(mcpdreport*) //p)
{
}




/*************************************************************************
Neural network training  using  modified  Levenberg-Marquardt  with  exact
Hessian calculation and regularization. Subroutine trains  neural  network
with restarts from random positions. Algorithm is well  suited  for  small
and medium scale problems (hundreds of weights).

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay constant, >=0.001
                    Decay term 'Decay*||Weights||^2' is added to error
                    function.
                    If you don't know what Decay to choose, use 0.001.
    Restarts    -   number of restarts from random position, >0.
                    If you don't know what Restarts to choose, use 2.

OUTPUT PARAMETERS:
    Network     -   trained neural network.
    Info        -   return code:
                    * -9, if internal matrix inverse subroutine failed
                    * -2, if there is a point with class number
                          outside of [0..NOut-1].
                    * -1, if wrong parameters specified
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void mlptrainlm(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     ae_int_t* info,
     mlpreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    double lmftol;
    double lmsteptol;
    ae_int_t i;
    ae_int_t k;
    double v;
    double e;
    double enew;
    double xnorm2;
    double stepnorm;
    ae_vector g;
    ae_vector d;
    ae_matrix h;
    ae_matrix hmod;
    ae_matrix z;
    ae_bool spd;
    double nu;
    double lambdav;
    double lambdaup;
    double lambdadown;
    minlbfgsreport internalrep;
    minlbfgsstate state;
    ae_vector x;
    ae_vector y;
    ae_vector wbase;
    ae_vector wdir;
    ae_vector wt;
    ae_vector wx;
    ae_int_t pass;
    ae_vector wbest;
    double ebest;
    ae_int_t invinfo;
    matinvreport invrep;
    ae_int_t solverinfo;
    densesolverreport solverrep;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _mlpreport_clear(rep);
    ae_vector_init(&g, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&d, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&h, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&hmod, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&z, 0, 0, DT_REAL, _state, ae_true);
    _minlbfgsreport_init(&internalrep, _state, ae_true);
    _minlbfgsstate_init(&state, _state, ae_true);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&wbase, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&wdir, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&wt, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&wx, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&wbest, 0, DT_REAL, _state, ae_true);
    _matinvreport_init(&invrep, _state, ae_true);
    _densesolverreport_init(&solverrep, _state, ae_true);

    mlpproperties(network, &nin, &nout, &wcount, _state);
    lambdaup = 10;
    lambdadown = 0.3;
    lmftol = 0.001;
    lmsteptol = 0.001;
    
    /*
     * Test for inputs
     */
    if( npoints<=0||restarts<1 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    if( mlpissoftmax(network, _state) )
    {
        for(i=0; i<=npoints-1; i++)
        {
            if( ae_round(xy->ptr.pp_double[i][nin], _state)<0||ae_round(xy->ptr.pp_double[i][nin], _state)>=nout )
            {
                *info = -2;
                ae_frame_leave(_state);
                return;
            }
        }
    }
    decay = ae_maxreal(decay, mlptrain_mindecay, _state);
    *info = 2;
    
    /*
     * Initialize data
     */
    rep->ngrad = 0;
    rep->nhess = 0;
    rep->ncholesky = 0;
    
    /*
     * General case.
     * Prepare task and network. Allocate space.
     */
    mlpinitpreprocessor(network, xy, npoints, _state);
    ae_vector_set_length(&g, wcount-1+1, _state);
    ae_matrix_set_length(&h, wcount-1+1, wcount-1+1, _state);
    ae_matrix_set_length(&hmod, wcount-1+1, wcount-1+1, _state);
    ae_vector_set_length(&wbase, wcount-1+1, _state);
    ae_vector_set_length(&wdir, wcount-1+1, _state);
    ae_vector_set_length(&wbest, wcount-1+1, _state);
    ae_vector_set_length(&wt, wcount-1+1, _state);
    ae_vector_set_length(&wx, wcount-1+1, _state);
    ebest = ae_maxrealnumber;
    
    /*
     * Multiple passes
     */
    for(pass=1; pass<=restarts; pass++)
    {
        
        /*
         * Initialize weights
         */
        mlprandomize(network, _state);
        
        /*
         * First stage of the hybrid algorithm: LBFGS
         */
        ae_v_move(&wbase.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        minlbfgscreate(wcount, ae_minint(wcount, 5, _state), &wbase, &state, _state);
        minlbfgssetcond(&state, 0, 0, 0, ae_maxint(25, wcount, _state), _state);
        while(minlbfgsiteration(&state, _state))
        {
            
            /*
             * gradient
             */
            ae_v_move(&network->weights.ptr.p_double[0], 1, &state.x.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            mlpgradbatch(network, xy, npoints, &state.f, &state.g, _state);
            
            /*
             * weight decay
             */
            v = ae_v_dotproduct(&network->weights.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            state.f = state.f+0.5*decay*v;
            ae_v_addd(&state.g.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1), decay);
            
            /*
             * next iteration
             */
            rep->ngrad = rep->ngrad+1;
        }
        minlbfgsresults(&state, &wbase, &internalrep, _state);
        ae_v_move(&network->weights.ptr.p_double[0], 1, &wbase.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        
        /*
         * Second stage of the hybrid algorithm: LM
         *
         * Initialize H with identity matrix,
         * G with gradient,
         * E with regularized error.
         */
        mlphessianbatch(network, xy, npoints, &e, &g, &h, _state);
        v = ae_v_dotproduct(&network->weights.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        e = e+0.5*decay*v;
        ae_v_addd(&g.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1), decay);
        for(k=0; k<=wcount-1; k++)
        {
            h.ptr.pp_double[k][k] = h.ptr.pp_double[k][k]+decay;
        }
        rep->nhess = rep->nhess+1;
        lambdav = 0.001;
        nu = 2;
        for(;;)
        {
            
            /*
             * 1. HMod = H+lambda*I
             * 2. Try to solve (H+Lambda*I)*dx = -g.
             *    Increase lambda if left part is not positive definite.
             */
            for(i=0; i<=wcount-1; i++)
            {
                ae_v_move(&hmod.ptr.pp_double[i][0], 1, &h.ptr.pp_double[i][0], 1, ae_v_len(0,wcount-1));
                hmod.ptr.pp_double[i][i] = hmod.ptr.pp_double[i][i]+lambdav;
            }
            spd = spdmatrixcholesky(&hmod, wcount, ae_true, _state);
            rep->ncholesky = rep->ncholesky+1;
            if( !spd )
            {
                lambdav = lambdav*lambdaup*nu;
                nu = nu*2;
                continue;
            }
            spdmatrixcholeskysolve(&hmod, wcount, ae_true, &g, &solverinfo, &solverrep, &wdir, _state);
            if( solverinfo<0 )
            {
                lambdav = lambdav*lambdaup*nu;
                nu = nu*2;
                continue;
            }
            ae_v_muld(&wdir.ptr.p_double[0], 1, ae_v_len(0,wcount-1), -1);
            
            /*
             * Lambda found.
             * 1. Save old w in WBase
             * 1. Test some stopping criterions
             * 2. If error(w+wdir)>error(w), increase lambda
             */
            ae_v_add(&network->weights.ptr.p_double[0], 1, &wdir.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            xnorm2 = ae_v_dotproduct(&network->weights.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            stepnorm = ae_v_dotproduct(&wdir.ptr.p_double[0], 1, &wdir.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            stepnorm = ae_sqrt(stepnorm, _state);
            enew = mlperror(network, xy, npoints, _state)+0.5*decay*xnorm2;
            if( ae_fp_less(stepnorm,lmsteptol*(1+ae_sqrt(xnorm2, _state))) )
            {
                break;
            }
            if( ae_fp_greater(enew,e) )
            {
                lambdav = lambdav*lambdaup*nu;
                nu = nu*2;
                continue;
            }
            
            /*
             * Optimize using inv(cholesky(H)) as preconditioner
             */
            rmatrixtrinverse(&hmod, wcount, ae_true, ae_false, &invinfo, &invrep, _state);
            if( invinfo<=0 )
            {
                
                /*
                 * if matrix can't be inverted then exit with errors
                 * TODO: make WCount steps in direction suggested by HMod
                 */
                *info = -9;
                ae_frame_leave(_state);
                return;
            }
            ae_v_move(&wbase.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            for(i=0; i<=wcount-1; i++)
            {
                wt.ptr.p_double[i] = 0;
            }
            minlbfgscreatex(wcount, wcount, &wt, 1, 0.0, &state, _state);
            minlbfgssetcond(&state, 0, 0, 0, 5, _state);
            while(minlbfgsiteration(&state, _state))
            {
                
                /*
                 * gradient
                 */
                for(i=0; i<=wcount-1; i++)
                {
                    v = ae_v_dotproduct(&state.x.ptr.p_double[i], 1, &hmod.ptr.pp_double[i][i], 1, ae_v_len(i,wcount-1));
                    network->weights.ptr.p_double[i] = wbase.ptr.p_double[i]+v;
                }
                mlpgradbatch(network, xy, npoints, &state.f, &g, _state);
                for(i=0; i<=wcount-1; i++)
                {
                    state.g.ptr.p_double[i] = 0;
                }
                for(i=0; i<=wcount-1; i++)
                {
                    v = g.ptr.p_double[i];
                    ae_v_addd(&state.g.ptr.p_double[i], 1, &hmod.ptr.pp_double[i][i], 1, ae_v_len(i,wcount-1), v);
                }
                
                /*
                 * weight decay
                 * grad(x'*x) = A'*(x0+A*t)
                 */
                v = ae_v_dotproduct(&network->weights.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
                state.f = state.f+0.5*decay*v;
                for(i=0; i<=wcount-1; i++)
                {
                    v = decay*network->weights.ptr.p_double[i];
                    ae_v_addd(&state.g.ptr.p_double[i], 1, &hmod.ptr.pp_double[i][i], 1, ae_v_len(i,wcount-1), v);
                }
                
                /*
                 * next iteration
                 */
                rep->ngrad = rep->ngrad+1;
            }
            minlbfgsresults(&state, &wt, &internalrep, _state);
            
            /*
             * Accept new position.
             * Calculate Hessian
             */
            for(i=0; i<=wcount-1; i++)
            {
                v = ae_v_dotproduct(&wt.ptr.p_double[i], 1, &hmod.ptr.pp_double[i][i], 1, ae_v_len(i,wcount-1));
                network->weights.ptr.p_double[i] = wbase.ptr.p_double[i]+v;
            }
            mlphessianbatch(network, xy, npoints, &e, &g, &h, _state);
            v = ae_v_dotproduct(&network->weights.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            e = e+0.5*decay*v;
            ae_v_addd(&g.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1), decay);
            for(k=0; k<=wcount-1; k++)
            {
                h.ptr.pp_double[k][k] = h.ptr.pp_double[k][k]+decay;
            }
            rep->nhess = rep->nhess+1;
            
            /*
             * Update lambda
             */
            lambdav = lambdav*lambdadown;
            nu = 2;
        }
        
        /*
         * update WBest
         */
        v = ae_v_dotproduct(&network->weights.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        e = 0.5*decay*v+mlperror(network, xy, npoints, _state);
        if( ae_fp_less(e,ebest) )
        {
            ebest = e;
            ae_v_move(&wbest.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        }
    }
    
    /*
     * copy WBest to output
     */
    ae_v_move(&network->weights.ptr.p_double[0], 1, &wbest.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
    ae_frame_leave(_state);
}


/*************************************************************************
Neural  network  training  using  L-BFGS  algorithm  with  regularization.
Subroutine  trains  neural  network  with  restarts from random positions.
Algorithm  is  well  suited  for  problems  of  any dimensionality (memory
requirements and step complexity are linear by weights number).

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay constant, >=0.001
                    Decay term 'Decay*||Weights||^2' is added to error
                    function.
                    If you don't know what Decay to choose, use 0.001.
    Restarts    -   number of restarts from random position, >0.
                    If you don't know what Restarts to choose, use 2.
    WStep       -   stopping criterion. Algorithm stops if  step  size  is
                    less than WStep. Recommended value - 0.01.  Zero  step
                    size means stopping after MaxIts iterations.
    MaxIts      -   stopping   criterion.  Algorithm  stops  after  MaxIts
                    iterations (NOT gradient  calculations).  Zero  MaxIts
                    means stopping when step is sufficiently small.

OUTPUT PARAMETERS:
    Network     -   trained neural network.
    Info        -   return code:
                    * -8, if both WStep=0 and MaxIts=0
                    * -2, if there is a point with class number
                          outside of [0..NOut-1].
                    * -1, if wrong parameters specified
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report

  -- ALGLIB --
     Copyright 09.12.2007 by Bochkanov Sergey
*************************************************************************/
void mlptrainlbfgs(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     double wstep,
     ae_int_t maxits,
     ae_int_t* info,
     mlpreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t pass;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    ae_vector w;
    ae_vector wbest;
    double e;
    double v;
    double ebest;
    minlbfgsreport internalrep;
    minlbfgsstate state;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _mlpreport_clear(rep);
    ae_vector_init(&w, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&wbest, 0, DT_REAL, _state, ae_true);
    _minlbfgsreport_init(&internalrep, _state, ae_true);
    _minlbfgsstate_init(&state, _state, ae_true);

    
    /*
     * Test inputs, parse flags, read network geometry
     */
    if( ae_fp_eq(wstep,0)&&maxits==0 )
    {
        *info = -8;
        ae_frame_leave(_state);
        return;
    }
    if( ((npoints<=0||restarts<1)||ae_fp_less(wstep,0))||maxits<0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    mlpproperties(network, &nin, &nout, &wcount, _state);
    if( mlpissoftmax(network, _state) )
    {
        for(i=0; i<=npoints-1; i++)
        {
            if( ae_round(xy->ptr.pp_double[i][nin], _state)<0||ae_round(xy->ptr.pp_double[i][nin], _state)>=nout )
            {
                *info = -2;
                ae_frame_leave(_state);
                return;
            }
        }
    }
    decay = ae_maxreal(decay, mlptrain_mindecay, _state);
    *info = 2;
    
    /*
     * Prepare
     */
    mlpinitpreprocessor(network, xy, npoints, _state);
    ae_vector_set_length(&w, wcount-1+1, _state);
    ae_vector_set_length(&wbest, wcount-1+1, _state);
    ebest = ae_maxrealnumber;
    
    /*
     * Multiple starts
     */
    rep->ncholesky = 0;
    rep->nhess = 0;
    rep->ngrad = 0;
    for(pass=1; pass<=restarts; pass++)
    {
        
        /*
         * Process
         */
        mlprandomize(network, _state);
        ae_v_move(&w.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        minlbfgscreate(wcount, ae_minint(wcount, 10, _state), &w, &state, _state);
        minlbfgssetcond(&state, 0.0, 0.0, wstep, maxits, _state);
        while(minlbfgsiteration(&state, _state))
        {
            ae_v_move(&network->weights.ptr.p_double[0], 1, &state.x.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            mlpgradnbatch(network, xy, npoints, &state.f, &state.g, _state);
            v = ae_v_dotproduct(&network->weights.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            state.f = state.f+0.5*decay*v;
            ae_v_addd(&state.g.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1), decay);
            rep->ngrad = rep->ngrad+1;
        }
        minlbfgsresults(&state, &w, &internalrep, _state);
        ae_v_move(&network->weights.ptr.p_double[0], 1, &w.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        
        /*
         * Compare with best
         */
        v = ae_v_dotproduct(&network->weights.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        e = mlperrorn(network, xy, npoints, _state)+0.5*decay*v;
        if( ae_fp_less(e,ebest) )
        {
            ae_v_move(&wbest.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            ebest = e;
        }
    }
    
    /*
     * The best network
     */
    ae_v_move(&network->weights.ptr.p_double[0], 1, &wbest.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
    ae_frame_leave(_state);
}


/*************************************************************************
Neural network training using early stopping (base algorithm - L-BFGS with
regularization).

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry
    TrnXY       -   training set
    TrnSize     -   training set size, TrnSize>0
    ValXY       -   validation set
    ValSize     -   validation set size, ValSize>0
    Decay       -   weight decay constant, >=0.001
                    Decay term 'Decay*||Weights||^2' is added to error
                    function.
                    If you don't know what Decay to choose, use 0.001.
    Restarts    -   number of restarts, either:
                    * strictly positive number - algorithm make specified
                      number of restarts from random position.
                    * -1, in which case algorithm makes exactly one run
                      from the initial state of the network (no randomization).
                    If you don't know what Restarts to choose, choose one
                    one the following:
                    * -1 (deterministic start)
                    * +1 (one random restart)
                    * +5 (moderate amount of random restarts)

OUTPUT PARAMETERS:
    Network     -   trained neural network.
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NOut-1].
                    * -1, if wrong parameters specified
                          (NPoints<0, Restarts<1, ...).
                    *  2, task has been solved, stopping  criterion  met -
                          sufficiently small step size.  Not expected  (we
                          use  EARLY  stopping)  but  possible  and not an
                          error.
                    *  6, task has been solved, stopping  criterion  met -
                          increasing of validation set error.
    Rep         -   training report

NOTE:

Algorithm stops if validation set error increases for  a  long  enough  or
step size is small enought  (there  are  task  where  validation  set  may
decrease for eternity). In any case solution returned corresponds  to  the
minimum of validation set error.

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void mlptraines(multilayerperceptron* network,
     /* Real    */ ae_matrix* trnxy,
     ae_int_t trnsize,
     /* Real    */ ae_matrix* valxy,
     ae_int_t valsize,
     double decay,
     ae_int_t restarts,
     ae_int_t* info,
     mlpreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t pass;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;
    ae_vector w;
    ae_vector wbest;
    double e;
    double v;
    double ebest;
    ae_vector wfinal;
    double efinal;
    ae_int_t itcnt;
    ae_int_t itbest;
    minlbfgsreport internalrep;
    minlbfgsstate state;
    double wstep;
    ae_bool needrandomization;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _mlpreport_clear(rep);
    ae_vector_init(&w, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&wbest, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&wfinal, 0, DT_REAL, _state, ae_true);
    _minlbfgsreport_init(&internalrep, _state, ae_true);
    _minlbfgsstate_init(&state, _state, ae_true);

    wstep = 0.001;
    
    /*
     * Test inputs, parse flags, read network geometry
     */
    if( ((trnsize<=0||valsize<=0)||(restarts<1&&restarts!=-1))||ae_fp_less(decay,0) )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    if( restarts==-1 )
    {
        needrandomization = ae_false;
        restarts = 1;
    }
    else
    {
        needrandomization = ae_true;
    }
    mlpproperties(network, &nin, &nout, &wcount, _state);
    if( mlpissoftmax(network, _state) )
    {
        for(i=0; i<=trnsize-1; i++)
        {
            if( ae_round(trnxy->ptr.pp_double[i][nin], _state)<0||ae_round(trnxy->ptr.pp_double[i][nin], _state)>=nout )
            {
                *info = -2;
                ae_frame_leave(_state);
                return;
            }
        }
        for(i=0; i<=valsize-1; i++)
        {
            if( ae_round(valxy->ptr.pp_double[i][nin], _state)<0||ae_round(valxy->ptr.pp_double[i][nin], _state)>=nout )
            {
                *info = -2;
                ae_frame_leave(_state);
                return;
            }
        }
    }
    *info = 2;
    
    /*
     * Prepare
     */
    mlpinitpreprocessor(network, trnxy, trnsize, _state);
    ae_vector_set_length(&w, wcount-1+1, _state);
    ae_vector_set_length(&wbest, wcount-1+1, _state);
    ae_vector_set_length(&wfinal, wcount-1+1, _state);
    efinal = ae_maxrealnumber;
    for(i=0; i<=wcount-1; i++)
    {
        wfinal.ptr.p_double[i] = 0;
    }
    
    /*
     * Multiple starts
     */
    rep->ncholesky = 0;
    rep->nhess = 0;
    rep->ngrad = 0;
    for(pass=1; pass<=restarts; pass++)
    {
        
        /*
         * Process
         */
        if( needrandomization )
        {
            mlprandomize(network, _state);
        }
        ebest = mlperror(network, valxy, valsize, _state);
        ae_v_move(&wbest.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        itbest = 0;
        itcnt = 0;
        ae_v_move(&w.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
        minlbfgscreate(wcount, ae_minint(wcount, 10, _state), &w, &state, _state);
        minlbfgssetcond(&state, 0.0, 0.0, wstep, 0, _state);
        minlbfgssetxrep(&state, ae_true, _state);
        while(minlbfgsiteration(&state, _state))
        {
            
            /*
             * Calculate gradient
             */
            if( state.needfg )
            {
                ae_v_move(&network->weights.ptr.p_double[0], 1, &state.x.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
                mlpgradnbatch(network, trnxy, trnsize, &state.f, &state.g, _state);
                v = ae_v_dotproduct(&network->weights.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
                state.f = state.f+0.5*decay*v;
                ae_v_addd(&state.g.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1), decay);
                rep->ngrad = rep->ngrad+1;
            }
            
            /*
             * Validation set
             */
            if( state.xupdated )
            {
                ae_v_move(&network->weights.ptr.p_double[0], 1, &state.x.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
                e = mlperror(network, valxy, valsize, _state);
                if( ae_fp_less(e,ebest) )
                {
                    ebest = e;
                    ae_v_move(&wbest.ptr.p_double[0], 1, &network->weights.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
                    itbest = itcnt;
                }
                if( itcnt>30&&ae_fp_greater(itcnt,1.5*itbest) )
                {
                    *info = 6;
                    break;
                }
                itcnt = itcnt+1;
            }
        }
        minlbfgsresults(&state, &w, &internalrep, _state);
        
        /*
         * Compare with final answer
         */
        if( ae_fp_less(ebest,efinal) )
        {
            ae_v_move(&wfinal.ptr.p_double[0], 1, &wbest.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
            efinal = ebest;
        }
    }
    
    /*
     * The best network
     */
    ae_v_move(&network->weights.ptr.p_double[0], 1, &wfinal.ptr.p_double[0], 1, ae_v_len(0,wcount-1));
    ae_frame_leave(_state);
}


/*************************************************************************
Cross-validation estimate of generalization error.

Base algorithm - L-BFGS.

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry.   Network is
                    not changed during cross-validation -  it is used only
                    as a representative of its architecture.
    XY          -   training set.
    SSize       -   training set size
    Decay       -   weight  decay, same as in MLPTrainLBFGS
    Restarts    -   number of restarts, >0.
                    restarts are counted for each partition separately, so
                    total number of restarts will be Restarts*FoldsCount.
    WStep       -   stopping criterion, same as in MLPTrainLBFGS
    MaxIts      -   stopping criterion, same as in MLPTrainLBFGS
    FoldsCount  -   number of folds in k-fold cross-validation,
                    2<=FoldsCount<=SSize.
                    recommended value: 10.

OUTPUT PARAMETERS:
    Info        -   return code, same as in MLPTrainLBFGS
    Rep         -   report, same as in MLPTrainLM/MLPTrainLBFGS
    CVRep       -   generalization error estimates

  -- ALGLIB --
     Copyright 09.12.2007 by Bochkanov Sergey
*************************************************************************/
void mlpkfoldcvlbfgs(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     double wstep,
     ae_int_t maxits,
     ae_int_t foldscount,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* cvrep,
     ae_state *_state)
{

    *info = 0;
    _mlpreport_clear(rep);
    _mlpcvreport_clear(cvrep);

    mlptrain_mlpkfoldcvgeneral(network, xy, npoints, decay, restarts, foldscount, ae_false, wstep, maxits, info, rep, cvrep, _state);
}


/*************************************************************************
Cross-validation estimate of generalization error.

Base algorithm - Levenberg-Marquardt.

INPUT PARAMETERS:
    Network     -   neural network with initialized geometry.   Network is
                    not changed during cross-validation -  it is used only
                    as a representative of its architecture.
    XY          -   training set.
    SSize       -   training set size
    Decay       -   weight  decay, same as in MLPTrainLBFGS
    Restarts    -   number of restarts, >0.
                    restarts are counted for each partition separately, so
                    total number of restarts will be Restarts*FoldsCount.
    FoldsCount  -   number of folds in k-fold cross-validation,
                    2<=FoldsCount<=SSize.
                    recommended value: 10.

OUTPUT PARAMETERS:
    Info        -   return code, same as in MLPTrainLBFGS
    Rep         -   report, same as in MLPTrainLM/MLPTrainLBFGS
    CVRep       -   generalization error estimates

  -- ALGLIB --
     Copyright 09.12.2007 by Bochkanov Sergey
*************************************************************************/
void mlpkfoldcvlm(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     ae_int_t foldscount,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* cvrep,
     ae_state *_state)
{

    *info = 0;
    _mlpreport_clear(rep);
    _mlpcvreport_clear(cvrep);

    mlptrain_mlpkfoldcvgeneral(network, xy, npoints, decay, restarts, foldscount, ae_true, 0.0, 0, info, rep, cvrep, _state);
}


/*************************************************************************
Internal cross-validation subroutine
*************************************************************************/
static void mlptrain_mlpkfoldcvgeneral(multilayerperceptron* n,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     ae_int_t foldscount,
     ae_bool lmalgorithm,
     double wstep,
     ae_int_t maxits,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* cvrep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t fold;
    ae_int_t j;
    ae_int_t k;
    multilayerperceptron network;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t rowlen;
    ae_int_t wcount;
    ae_int_t nclasses;
    ae_int_t tssize;
    ae_int_t cvssize;
    ae_matrix cvset;
    ae_matrix testset;
    ae_vector folds;
    ae_int_t relcnt;
    mlpreport internalrep;
    ae_vector x;
    ae_vector y;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _mlpreport_clear(rep);
    _mlpcvreport_clear(cvrep);
    _multilayerperceptron_init(&network, _state, ae_true);
    ae_matrix_init(&cvset, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&testset, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&folds, 0, DT_INT, _state, ae_true);
    _mlpreport_init(&internalrep, _state, ae_true);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);

    
    /*
     * Read network geometry, test parameters
     */
    mlpproperties(n, &nin, &nout, &wcount, _state);
    if( mlpissoftmax(n, _state) )
    {
        nclasses = nout;
        rowlen = nin+1;
    }
    else
    {
        nclasses = -nout;
        rowlen = nin+nout;
    }
    if( (npoints<=0||foldscount<2)||foldscount>npoints )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    mlpcopy(n, &network, _state);
    
    /*
     * K-fold out cross-validation.
     * First, estimate generalization error
     */
    ae_matrix_set_length(&testset, npoints-1+1, rowlen-1+1, _state);
    ae_matrix_set_length(&cvset, npoints-1+1, rowlen-1+1, _state);
    ae_vector_set_length(&x, nin-1+1, _state);
    ae_vector_set_length(&y, nout-1+1, _state);
    mlptrain_mlpkfoldsplit(xy, npoints, nclasses, foldscount, ae_false, &folds, _state);
    cvrep->relclserror = 0;
    cvrep->avgce = 0;
    cvrep->rmserror = 0;
    cvrep->avgerror = 0;
    cvrep->avgrelerror = 0;
    rep->ngrad = 0;
    rep->nhess = 0;
    rep->ncholesky = 0;
    relcnt = 0;
    for(fold=0; fold<=foldscount-1; fold++)
    {
        
        /*
         * Separate set
         */
        tssize = 0;
        cvssize = 0;
        for(i=0; i<=npoints-1; i++)
        {
            if( folds.ptr.p_int[i]==fold )
            {
                ae_v_move(&testset.ptr.pp_double[tssize][0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,rowlen-1));
                tssize = tssize+1;
            }
            else
            {
                ae_v_move(&cvset.ptr.pp_double[cvssize][0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,rowlen-1));
                cvssize = cvssize+1;
            }
        }
        
        /*
         * Train on CV training set
         */
        if( lmalgorithm )
        {
            mlptrainlm(&network, &cvset, cvssize, decay, restarts, info, &internalrep, _state);
        }
        else
        {
            mlptrainlbfgs(&network, &cvset, cvssize, decay, restarts, wstep, maxits, info, &internalrep, _state);
        }
        if( *info<0 )
        {
            cvrep->relclserror = 0;
            cvrep->avgce = 0;
            cvrep->rmserror = 0;
            cvrep->avgerror = 0;
            cvrep->avgrelerror = 0;
            ae_frame_leave(_state);
            return;
        }
        rep->ngrad = rep->ngrad+internalrep.ngrad;
        rep->nhess = rep->nhess+internalrep.nhess;
        rep->ncholesky = rep->ncholesky+internalrep.ncholesky;
        
        /*
         * Estimate error using CV test set
         */
        if( mlpissoftmax(&network, _state) )
        {
            
            /*
             * classification-only code
             */
            cvrep->relclserror = cvrep->relclserror+mlpclserror(&network, &testset, tssize, _state);
            cvrep->avgce = cvrep->avgce+mlperrorn(&network, &testset, tssize, _state);
        }
        for(i=0; i<=tssize-1; i++)
        {
            ae_v_move(&x.ptr.p_double[0], 1, &testset.ptr.pp_double[i][0], 1, ae_v_len(0,nin-1));
            mlpprocess(&network, &x, &y, _state);
            if( mlpissoftmax(&network, _state) )
            {
                
                /*
                 * Classification-specific code
                 */
                k = ae_round(testset.ptr.pp_double[i][nin], _state);
                for(j=0; j<=nout-1; j++)
                {
                    if( j==k )
                    {
                        cvrep->rmserror = cvrep->rmserror+ae_sqr(y.ptr.p_double[j]-1, _state);
                        cvrep->avgerror = cvrep->avgerror+ae_fabs(y.ptr.p_double[j]-1, _state);
                        cvrep->avgrelerror = cvrep->avgrelerror+ae_fabs(y.ptr.p_double[j]-1, _state);
                        relcnt = relcnt+1;
                    }
                    else
                    {
                        cvrep->rmserror = cvrep->rmserror+ae_sqr(y.ptr.p_double[j], _state);
                        cvrep->avgerror = cvrep->avgerror+ae_fabs(y.ptr.p_double[j], _state);
                    }
                }
            }
            else
            {
                
                /*
                 * Regression-specific code
                 */
                for(j=0; j<=nout-1; j++)
                {
                    cvrep->rmserror = cvrep->rmserror+ae_sqr(y.ptr.p_double[j]-testset.ptr.pp_double[i][nin+j], _state);
                    cvrep->avgerror = cvrep->avgerror+ae_fabs(y.ptr.p_double[j]-testset.ptr.pp_double[i][nin+j], _state);
                    if( ae_fp_neq(testset.ptr.pp_double[i][nin+j],0) )
                    {
                        cvrep->avgrelerror = cvrep->avgrelerror+ae_fabs((y.ptr.p_double[j]-testset.ptr.pp_double[i][nin+j])/testset.ptr.pp_double[i][nin+j], _state);
                        relcnt = relcnt+1;
                    }
                }
            }
        }
    }
    if( mlpissoftmax(&network, _state) )
    {
        cvrep->relclserror = cvrep->relclserror/npoints;
        cvrep->avgce = cvrep->avgce/(ae_log(2, _state)*npoints);
    }
    cvrep->rmserror = ae_sqrt(cvrep->rmserror/(npoints*nout), _state);
    cvrep->avgerror = cvrep->avgerror/(npoints*nout);
    cvrep->avgrelerror = cvrep->avgrelerror/relcnt;
    *info = 1;
    ae_frame_leave(_state);
}


/*************************************************************************
Subroutine prepares K-fold split of the training set.

NOTES:
    "NClasses>0" means that we have classification task.
    "NClasses<0" means regression task with -NClasses real outputs.
*************************************************************************/
static void mlptrain_mlpkfoldsplit(/* Real    */ ae_matrix* ,//xy,
     ae_int_t npoints,
     ae_int_t nclasses,
     ae_int_t foldscount,
     ae_bool stratifiedsplits,
     /* Integer */ ae_vector* folds,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;

    ae_vector_clear(folds);

    
    /*
     * test parameters
     */
    ae_assert(npoints>0, "MLPKFoldSplit: wrong NPoints!", _state);
    ae_assert(nclasses>1||nclasses<0, "MLPKFoldSplit: wrong NClasses!", _state);
    ae_assert(foldscount>=2&&foldscount<=npoints, "MLPKFoldSplit: wrong FoldsCount!", _state);
    ae_assert(!stratifiedsplits, "MLPKFoldSplit: stratified splits are not supported!", _state);
    
    /*
     * Folds
     */
    ae_vector_set_length(folds, npoints-1+1, _state);
    for(i=0; i<=npoints-1; i++)
    {
        folds->ptr.p_int[i] = i*foldscount/npoints;
    }
    for(i=0; i<=npoints-2; i++)
    {
        j = i+ae_randominteger(npoints-i, _state);
        if( j!=i )
        {
            k = folds->ptr.p_int[i];
            folds->ptr.p_int[i] = folds->ptr.p_int[j];
            folds->ptr.p_int[j] = k;
        }
    }
}


ae_bool _mlpreport_init(mlpreport* /*p*/, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    return ae_true;
}


ae_bool _mlpreport_init_copy(mlpreport* dst, mlpreport* src, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    dst->ngrad = src->ngrad;
    dst->nhess = src->nhess;
    dst->ncholesky = src->ncholesky;
    return ae_true;
}


void _mlpreport_clear(mlpreport*) //p)
{
}


ae_bool _mlpcvreport_init(mlpcvreport* /*p*/, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    return ae_true;
}


ae_bool _mlpcvreport_init_copy(mlpcvreport* dst, mlpcvreport* src, ae_state */*_state*/, ae_bool ) //make_automatic)
{
    dst->relclserror = src->relclserror;
    dst->avgce = src->avgce;
    dst->rmserror = src->rmserror;
    dst->avgerror = src->avgerror;
    dst->avgrelerror = src->avgrelerror;
    return ae_true;
}


void _mlpcvreport_clear(mlpcvreport* ) //p)
{
}




/*************************************************************************
Like MLPCreate0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreate0(ae_int_t nin,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreate0(nin, nout, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreate1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreate1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreate1(nin, nhid, nout, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreate2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreate2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreate2(nin, nhid1, nhid2, nout, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreateB0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreateb0(ae_int_t nin,
     ae_int_t nout,
     double b,
     double d,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreateb0(nin, nout, b, d, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreateB1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreateb1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     double b,
     double d,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreateb1(nin, nhid, nout, b, d, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreateB2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreateb2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     double b,
     double d,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreateb2(nin, nhid1, nhid2, nout, b, d, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreateR0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreater0(ae_int_t nin,
     ae_int_t nout,
     double a,
     double b,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreater0(nin, nout, a, b, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreateR1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreater1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     double a,
     double b,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreater1(nin, nhid, nout, a, b, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreateR2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreater2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     double a,
     double b,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreater2(nin, nhid1, nhid2, nout, a, b, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreateC0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatec0(ae_int_t nin,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreatec0(nin, nout, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreateC1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatec1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreatec1(nin, nhid, nout, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Like MLPCreateC2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatec2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_frame _frame_block;
    multilayerperceptron net;

    ae_frame_make(_state, &_frame_block);
    _mlpensemble_clear(ensemble);
    _multilayerperceptron_init(&net, _state, ae_true);

    mlpcreatec2(nin, nhid1, nhid2, nout, &net, _state);
    mlpecreatefromnetwork(&net, ensemblesize, ensemble, _state);
    ae_frame_leave(_state);
}


/*************************************************************************
Creates ensemble from network. Only network geometry is copied.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatefromnetwork(multilayerperceptron* network,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t ccount;
    ae_int_t wcount;

    _mlpensemble_clear(ensemble);

    ae_assert(ensemblesize>0, "MLPECreate: incorrect ensemble size!", _state);
    
    /*
     * Copy network
     */
    mlpcopy(network, &ensemble->network, _state);
    
    /*
     * network properties
     */
    if( mlpissoftmax(network, _state) )
    {
        ccount = mlpgetinputscount(&ensemble->network, _state);
    }
    else
    {
        ccount = mlpgetinputscount(&ensemble->network, _state)+mlpgetoutputscount(&ensemble->network, _state);
    }
    wcount = mlpgetweightscount(&ensemble->network, _state);
    ensemble->ensemblesize = ensemblesize;
    
    /*
     * weights, means, sigmas
     */
    ae_vector_set_length(&ensemble->weights, ensemblesize*wcount, _state);
    ae_vector_set_length(&ensemble->columnmeans, ensemblesize*ccount, _state);
    ae_vector_set_length(&ensemble->columnsigmas, ensemblesize*ccount, _state);
    for(i=0; i<=ensemblesize*wcount-1; i++)
    {
        ensemble->weights.ptr.p_double[i] = ae_randomreal(_state)-0.5;
    }
    for(i=0; i<=ensemblesize-1; i++)
    {
        ae_v_move(&ensemble->columnmeans.ptr.p_double[i*ccount], 1, &network->columnmeans.ptr.p_double[0], 1, ae_v_len(i*ccount,(i+1)*ccount-1));
        ae_v_move(&ensemble->columnsigmas.ptr.p_double[i*ccount], 1, &network->columnsigmas.ptr.p_double[0], 1, ae_v_len(i*ccount,(i+1)*ccount-1));
    }
    
    /*
     * temporaries, internal buffers
     */
    ae_vector_set_length(&ensemble->y, mlpgetoutputscount(&ensemble->network, _state), _state);
}


/*************************************************************************
Copying of MLPEnsemble strucure

INPUT PARAMETERS:
    Ensemble1 -   original

OUTPUT PARAMETERS:
    Ensemble2 -   copy

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecopy(mlpensemble* ensemble1,
     mlpensemble* ensemble2,
     ae_state *_state)
{
    ae_int_t ccount;
    ae_int_t wcount;

    _mlpensemble_clear(ensemble2);

    
    /*
     * Unload info
     */
    if( mlpissoftmax(&ensemble1->network, _state) )
    {
        ccount = mlpgetinputscount(&ensemble1->network, _state);
    }
    else
    {
        ccount = mlpgetinputscount(&ensemble1->network, _state)+mlpgetoutputscount(&ensemble1->network, _state);
    }
    wcount = mlpgetweightscount(&ensemble1->network, _state);
    
    /*
     * Allocate space
     */
    ae_vector_set_length(&ensemble2->weights, ensemble1->ensemblesize*wcount, _state);
    ae_vector_set_length(&ensemble2->columnmeans, ensemble1->ensemblesize*ccount, _state);
    ae_vector_set_length(&ensemble2->columnsigmas, ensemble1->ensemblesize*ccount, _state);
    ae_vector_set_length(&ensemble2->y, mlpgetoutputscount(&ensemble1->network, _state), _state);
    
    /*
     * Copy
     */
    ensemble2->ensemblesize = ensemble1->ensemblesize;
    ae_v_move(&ensemble2->weights.ptr.p_double[0], 1, &ensemble1->weights.ptr.p_double[0], 1, ae_v_len(0,ensemble1->ensemblesize*wcount-1));
    ae_v_move(&ensemble2->columnmeans.ptr.p_double[0], 1, &ensemble1->columnmeans.ptr.p_double[0], 1, ae_v_len(0,ensemble1->ensemblesize*ccount-1));
    ae_v_move(&ensemble2->columnsigmas.ptr.p_double[0], 1, &ensemble1->columnsigmas.ptr.p_double[0], 1, ae_v_len(0,ensemble1->ensemblesize*ccount-1));
    mlpcopy(&ensemble1->network, &ensemble2->network, _state);
}


/*************************************************************************
Randomization of MLP ensemble

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlperandomize(mlpensemble* ensemble, ae_state *_state)
{
    ae_int_t i;
    ae_int_t wcount;


    wcount = mlpgetweightscount(&ensemble->network, _state);
    for(i=0; i<=ensemble->ensemblesize*wcount-1; i++)
    {
        ensemble->weights.ptr.p_double[i] = ae_randomreal(_state)-0.5;
    }
}


/*************************************************************************
Return ensemble properties (number of inputs and outputs).

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpeproperties(mlpensemble* ensemble,
     ae_int_t* nin,
     ae_int_t* nout,
     ae_state *_state)
{

    *nin = 0;
    *nout = 0;

    *nin = mlpgetinputscount(&ensemble->network, _state);
    *nout = mlpgetoutputscount(&ensemble->network, _state);
}


/*************************************************************************
Return normalization type (whether ensemble is SOFTMAX-normalized or not).

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
ae_bool mlpeissoftmax(mlpensemble* ensemble, ae_state *_state)
{
    ae_bool result;


    result = mlpissoftmax(&ensemble->network, _state);
    return result;
}


/*************************************************************************
Procesing

INPUT PARAMETERS:
    Ensemble-   neural networks ensemble
    X       -   input vector,  array[0..NIn-1].
    Y       -   (possibly) preallocated buffer; if size of Y is less than
                NOut, it will be reallocated. If it is large enough, it
                is NOT reallocated, so we can save some time on reallocation.


OUTPUT PARAMETERS:
    Y       -   result. Regression estimate when solving regression  task,
                vector of posterior probabilities for classification task.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpeprocess(mlpensemble* ensemble,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{
    ae_int_t i;
    ae_int_t es;
    ae_int_t wc;
    ae_int_t cc;
    double v;
    ae_int_t nout;


    if( y->cnt<mlpgetoutputscount(&ensemble->network, _state) )
    {
        ae_vector_set_length(y, mlpgetoutputscount(&ensemble->network, _state), _state);
    }
    es = ensemble->ensemblesize;
    wc = mlpgetweightscount(&ensemble->network, _state);
    if( mlpissoftmax(&ensemble->network, _state) )
    {
        cc = mlpgetinputscount(&ensemble->network, _state);
    }
    else
    {
        cc = mlpgetinputscount(&ensemble->network, _state)+mlpgetoutputscount(&ensemble->network, _state);
    }
    v = (double)1/(double)es;
    nout = mlpgetoutputscount(&ensemble->network, _state);
    for(i=0; i<=nout-1; i++)
    {
        y->ptr.p_double[i] = 0;
    }
    for(i=0; i<=es-1; i++)
    {
        ae_v_move(&ensemble->network.weights.ptr.p_double[0], 1, &ensemble->weights.ptr.p_double[i*wc], 1, ae_v_len(0,wc-1));
        ae_v_move(&ensemble->network.columnmeans.ptr.p_double[0], 1, &ensemble->columnmeans.ptr.p_double[i*cc], 1, ae_v_len(0,cc-1));
        ae_v_move(&ensemble->network.columnsigmas.ptr.p_double[0], 1, &ensemble->columnsigmas.ptr.p_double[i*cc], 1, ae_v_len(0,cc-1));
        mlpprocess(&ensemble->network, x, &ensemble->y, _state);
        ae_v_addd(&y->ptr.p_double[0], 1, &ensemble->y.ptr.p_double[0], 1, ae_v_len(0,nout-1), v);
    }
}


/*************************************************************************
'interactive'  variant  of  MLPEProcess  for  languages  like Python which
support constructs like "Y = MLPEProcess(LM,X)" and interactive mode of the
interpreter

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpeprocessi(mlpensemble* ensemble,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state)
{

    ae_vector_clear(y);

    mlpeprocess(ensemble, x, y, _state);
}


/*************************************************************************
Relative classification error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    percent of incorrectly classified cases.
    Works both for classifier betwork and for regression networks which
are used as classifiers.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlperelclserror(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    double relcls;
    double avgce;
    double rms;
    double avg;
    double avgrel;
    double result;


    mlpe_mlpeallerrors(ensemble, xy, npoints, &relcls, &avgce, &rms, &avg, &avgrel, _state);
    result = relcls;
    return result;
}


/*************************************************************************
Average cross-entropy (in bits per element) on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    CrossEntropy/(NPoints*LN(2)).
    Zero if ensemble solves regression task.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlpeavgce(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    double relcls;
    double avgce;
    double rms;
    double avg;
    double avgrel;
    double result;


    mlpe_mlpeallerrors(ensemble, xy, npoints, &relcls, &avgce, &rms, &avg, &avgrel, _state);
    result = avgce;
    return result;
}


/*************************************************************************
RMS error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    root mean square error.
    Its meaning for regression task is obvious. As for classification task
RMS error means error when estimating posterior probabilities.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlpermserror(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    double relcls;
    double avgce;
    double rms;
    double avg;
    double avgrel;
    double result;


    mlpe_mlpeallerrors(ensemble, xy, npoints, &relcls, &avgce, &rms, &avg, &avgrel, _state);
    result = rms;
    return result;
}


/*************************************************************************
Average error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for classification task
it means average error when estimating posterior probabilities.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlpeavgerror(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    double relcls;
    double avgce;
    double rms;
    double avg;
    double avgrel;
    double result;


    mlpe_mlpeallerrors(ensemble, xy, npoints, &relcls, &avgce, &rms, &avg, &avgrel, _state);
    result = avg;
    return result;
}


/*************************************************************************
Average relative error on the test set

INPUT PARAMETERS:
    Ensemble-   ensemble
    XY      -   test set
    NPoints -   test set size

RESULT:
    Its meaning for regression task is obvious. As for classification task
it means average relative error when estimating posterior probabilities.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
double mlpeavgrelerror(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state)
{
    double relcls;
    double avgce;
    double rms;
    double avg;
    double avgrel;
    double result;


    mlpe_mlpeallerrors(ensemble, xy, npoints, &relcls, &avgce, &rms, &avg, &avgrel, _state);
    result = avgrel;
    return result;
}


/*************************************************************************
Training neural networks ensemble using  bootstrap  aggregating (bagging).
Modified Levenberg-Marquardt algorithm is used as base training method.

INPUT PARAMETERS:
    Ensemble    -   model with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay coefficient, >=0.001
    Restarts    -   restarts, >0.

OUTPUT PARAMETERS:
    Ensemble    -   trained model
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report.
    OOBErrors   -   out-of-bag generalization error estimate

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpebagginglm(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* ooberrors,
     ae_state *_state)
{

    *info = 0;
    _mlpreport_clear(rep);
    _mlpcvreport_clear(ooberrors);

    mlpe_mlpebagginginternal(ensemble, xy, npoints, decay, restarts, 0.0, 0, ae_true, info, rep, ooberrors, _state);
}


/*************************************************************************
Training neural networks ensemble using  bootstrap  aggregating (bagging).
L-BFGS algorithm is used as base training method.

INPUT PARAMETERS:
    Ensemble    -   model with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay coefficient, >=0.001
    Restarts    -   restarts, >0.
    WStep       -   stopping criterion, same as in MLPTrainLBFGS
    MaxIts      -   stopping criterion, same as in MLPTrainLBFGS

OUTPUT PARAMETERS:
    Ensemble    -   trained model
    Info        -   return code:
                    * -8, if both WStep=0 and MaxIts=0
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<0, Restarts<1).
                    *  2, if task has been solved.
    Rep         -   training report.
    OOBErrors   -   out-of-bag generalization error estimate

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpebagginglbfgs(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     double wstep,
     ae_int_t maxits,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* ooberrors,
     ae_state *_state)
{

    *info = 0;
    _mlpreport_clear(rep);
    _mlpcvreport_clear(ooberrors);

    mlpe_mlpebagginginternal(ensemble, xy, npoints, decay, restarts, wstep, maxits, ae_false, info, rep, ooberrors, _state);
}


/*************************************************************************
Training neural networks ensemble using early stopping.

INPUT PARAMETERS:
    Ensemble    -   model with initialized geometry
    XY          -   training set
    NPoints     -   training set size
    Decay       -   weight decay coefficient, >=0.001
    Restarts    -   restarts, >0.

OUTPUT PARAMETERS:
    Ensemble    -   trained model
    Info        -   return code:
                    * -2, if there is a point with class number
                          outside of [0..NClasses-1].
                    * -1, if incorrect parameters was passed
                          (NPoints<0, Restarts<1).
                    *  6, if task has been solved.
    Rep         -   training report.
    OOBErrors   -   out-of-bag generalization error estimate

  -- ALGLIB --
     Copyright 10.03.2009 by Bochkanov Sergey
*************************************************************************/
void mlpetraines(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     ae_int_t* info,
     mlpreport* rep,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_int_t k;
    ae_int_t ccount;
    ae_int_t pcount;
    ae_matrix trnxy;
    ae_matrix valxy;
    ae_int_t trnsize;
    ae_int_t valsize;
    ae_int_t tmpinfo;
    mlpreport tmprep;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _mlpreport_clear(rep);
    ae_matrix_init(&trnxy, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&valxy, 0, 0, DT_REAL, _state, ae_true);
    _mlpreport_init(&tmprep, _state, ae_true);

    nin = mlpgetinputscount(&ensemble->network, _state);
    nout = mlpgetoutputscount(&ensemble->network, _state);
    wcount = mlpgetweightscount(&ensemble->network, _state);
    if( (npoints<2||restarts<1)||ae_fp_less(decay,0) )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    if( mlpissoftmax(&ensemble->network, _state) )
    {
        for(i=0; i<=npoints-1; i++)
        {
            if( ae_round(xy->ptr.pp_double[i][nin], _state)<0||ae_round(xy->ptr.pp_double[i][nin], _state)>=nout )
            {
                *info = -2;
                ae_frame_leave(_state);
                return;
            }
        }
    }
    *info = 6;
    
    /*
     * allocate
     */
    if( mlpissoftmax(&ensemble->network, _state) )
    {
        ccount = nin+1;
        pcount = nin;
    }
    else
    {
        ccount = nin+nout;
        pcount = nin+nout;
    }
    ae_matrix_set_length(&trnxy, npoints, ccount, _state);
    ae_matrix_set_length(&valxy, npoints, ccount, _state);
    rep->ngrad = 0;
    rep->nhess = 0;
    rep->ncholesky = 0;
    
    /*
     * train networks
     */
    for(k=0; k<=ensemble->ensemblesize-1; k++)
    {
        
        /*
         * Split set
         */
        do
        {
            trnsize = 0;
            valsize = 0;
            for(i=0; i<=npoints-1; i++)
            {
                if( ae_fp_less(ae_randomreal(_state),0.66) )
                {
                    
                    /*
                     * Assign sample to training set
                     */
                    ae_v_move(&trnxy.ptr.pp_double[trnsize][0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,ccount-1));
                    trnsize = trnsize+1;
                }
                else
                {
                    
                    /*
                     * Assign sample to validation set
                     */
                    ae_v_move(&valxy.ptr.pp_double[valsize][0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,ccount-1));
                    valsize = valsize+1;
                }
            }
        }
        while(!(trnsize!=0&&valsize!=0));
        
        /*
         * Train
         */
        mlptraines(&ensemble->network, &trnxy, trnsize, &valxy, valsize, decay, restarts, &tmpinfo, &tmprep, _state);
        if( tmpinfo<0 )
        {
            *info = tmpinfo;
            ae_frame_leave(_state);
            return;
        }
        
        /*
         * save results
         */
        ae_v_move(&ensemble->weights.ptr.p_double[k*wcount], 1, &ensemble->network.weights.ptr.p_double[0], 1, ae_v_len(k*wcount,(k+1)*wcount-1));
        ae_v_move(&ensemble->columnmeans.ptr.p_double[k*pcount], 1, &ensemble->network.columnmeans.ptr.p_double[0], 1, ae_v_len(k*pcount,(k+1)*pcount-1));
        ae_v_move(&ensemble->columnsigmas.ptr.p_double[k*pcount], 1, &ensemble->network.columnsigmas.ptr.p_double[0], 1, ae_v_len(k*pcount,(k+1)*pcount-1));
        rep->ngrad = rep->ngrad+tmprep.ngrad;
        rep->nhess = rep->nhess+tmprep.nhess;
        rep->ncholesky = rep->ncholesky+tmprep.ncholesky;
    }
    ae_frame_leave(_state);
}


/*************************************************************************
Serializer: allocation

  -- ALGLIB --
     Copyright 19.10.2011 by Bochkanov Sergey
*************************************************************************/
void mlpealloc(ae_serializer* s, mlpensemble* ensemble, ae_state *_state)
{


    ae_serializer_alloc_entry(s);
    ae_serializer_alloc_entry(s);
    ae_serializer_alloc_entry(s);
    allocrealarray(s, &ensemble->weights, -1, _state);
    allocrealarray(s, &ensemble->columnmeans, -1, _state);
    allocrealarray(s, &ensemble->columnsigmas, -1, _state);
    mlpalloc(s, &ensemble->network, _state);
}


/*************************************************************************
Serializer: serialization

  -- ALGLIB --
     Copyright 14.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpeserialize(ae_serializer* s,
     mlpensemble* ensemble,
     ae_state *_state)
{


    ae_serializer_serialize_int(s, getmlpeserializationcode(_state), _state);
    ae_serializer_serialize_int(s, mlpe_mlpefirstversion, _state);
    ae_serializer_serialize_int(s, ensemble->ensemblesize, _state);
    serializerealarray(s, &ensemble->weights, -1, _state);
    serializerealarray(s, &ensemble->columnmeans, -1, _state);
    serializerealarray(s, &ensemble->columnsigmas, -1, _state);
    mlpserialize(s, &ensemble->network, _state);
}


/*************************************************************************
Serializer: unserialization

  -- ALGLIB --
     Copyright 14.03.2011 by Bochkanov Sergey
*************************************************************************/
void mlpeunserialize(ae_serializer* s,
     mlpensemble* ensemble,
     ae_state *_state)
{
    ae_int_t i0;
    ae_int_t i1;

    _mlpensemble_clear(ensemble);

    
    /*
     * check correctness of header
     */
    ae_serializer_unserialize_int(s, &i0, _state);
    ae_assert(i0==getmlpeserializationcode(_state), "MLPEUnserialize: stream header corrupted", _state);
    ae_serializer_unserialize_int(s, &i1, _state);
    ae_assert(i1==mlpe_mlpefirstversion, "MLPEUnserialize: stream header corrupted", _state);
    
    /*
     * Create network
     */
    ae_serializer_unserialize_int(s, &ensemble->ensemblesize, _state);
    unserializerealarray(s, &ensemble->weights, _state);
    unserializerealarray(s, &ensemble->columnmeans, _state);
    unserializerealarray(s, &ensemble->columnsigmas, _state);
    mlpunserialize(s, &ensemble->network, _state);
    
    /*
     * Allocate termoraries
     */
    ae_vector_set_length(&ensemble->y, mlpgetoutputscount(&ensemble->network, _state), _state);
}


/*************************************************************************
Calculation of all types of errors

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
static void mlpe_mlpeallerrors(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double* relcls,
     double* avgce,
     double* rms,
     double* avg,
     double* avgrel,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_int_t i;
    ae_vector buf;
    ae_vector workx;
    ae_vector y;
    ae_vector dy;
    ae_int_t nin;
    ae_int_t nout;

    ae_frame_make(_state, &_frame_block);
    *relcls = 0;
    *avgce = 0;
    *rms = 0;
    *avg = 0;
    *avgrel = 0;
    ae_vector_init(&buf, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&workx, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&dy, 0, DT_REAL, _state, ae_true);

    nin = mlpgetinputscount(&ensemble->network, _state);
    nout = mlpgetoutputscount(&ensemble->network, _state);
    ae_vector_set_length(&workx, nin, _state);
    ae_vector_set_length(&y, nout, _state);
    if( mlpissoftmax(&ensemble->network, _state) )
    {
        ae_vector_set_length(&dy, 1, _state);
        dserrallocate(nout, &buf, _state);
    }
    else
    {
        ae_vector_set_length(&dy, nout, _state);
        dserrallocate(-nout, &buf, _state);
    }
    for(i=0; i<=npoints-1; i++)
    {
        ae_v_move(&workx.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nin-1));
        mlpeprocess(ensemble, &workx, &y, _state);
        if( mlpissoftmax(&ensemble->network, _state) )
        {
            dy.ptr.p_double[0] = xy->ptr.pp_double[i][nin];
        }
        else
        {
            ae_v_move(&dy.ptr.p_double[0], 1, &xy->ptr.pp_double[i][nin], 1, ae_v_len(0,nout-1));
        }
        dserraccumulate(&buf, &y, &dy, _state);
    }
    dserrfinish(&buf, _state);
    *relcls = buf.ptr.p_double[0];
    *avgce = buf.ptr.p_double[1];
    *rms = buf.ptr.p_double[2];
    *avg = buf.ptr.p_double[3];
    *avgrel = buf.ptr.p_double[4];
    ae_frame_leave(_state);
}


/*************************************************************************
Internal bagging subroutine.

  -- ALGLIB --
     Copyright 19.02.2009 by Bochkanov Sergey
*************************************************************************/
static void mlpe_mlpebagginginternal(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     double wstep,
     ae_int_t maxits,
     ae_bool lmalgorithm,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* ooberrors,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix xys;
    ae_vector s;
    ae_matrix oobbuf;
    ae_vector oobcntbuf;
    ae_vector x;
    ae_vector y;
    ae_vector dy;
    ae_vector dsbuf;
    ae_int_t ccnt;
    ae_int_t pcnt;
    ae_int_t i;
    ae_int_t j;
    ae_int_t k;
    double v;
    mlpreport tmprep;
    ae_int_t nin;
    ae_int_t nout;
    ae_int_t wcount;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    _mlpreport_clear(rep);
    _mlpcvreport_clear(ooberrors);
    ae_matrix_init(&xys, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&s, 0, DT_BOOL, _state, ae_true);
    ae_matrix_init(&oobbuf, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&oobcntbuf, 0, DT_INT, _state, ae_true);
    ae_vector_init(&x, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&y, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&dy, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&dsbuf, 0, DT_REAL, _state, ae_true);
    _mlpreport_init(&tmprep, _state, ae_true);

    nin = mlpgetinputscount(&ensemble->network, _state);
    nout = mlpgetoutputscount(&ensemble->network, _state);
    wcount = mlpgetweightscount(&ensemble->network, _state);
    
    /*
     * Test for inputs
     */
    if( (!lmalgorithm&&ae_fp_eq(wstep,0))&&maxits==0 )
    {
        *info = -8;
        ae_frame_leave(_state);
        return;
    }
    if( ((npoints<=0||restarts<1)||ae_fp_less(wstep,0))||maxits<0 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    if( mlpissoftmax(&ensemble->network, _state) )
    {
        for(i=0; i<=npoints-1; i++)
        {
            if( ae_round(xy->ptr.pp_double[i][nin], _state)<0||ae_round(xy->ptr.pp_double[i][nin], _state)>=nout )
            {
                *info = -2;
                ae_frame_leave(_state);
                return;
            }
        }
    }
    
    /*
     * allocate temporaries
     */
    *info = 2;
    rep->ngrad = 0;
    rep->nhess = 0;
    rep->ncholesky = 0;
    ooberrors->relclserror = 0;
    ooberrors->avgce = 0;
    ooberrors->rmserror = 0;
    ooberrors->avgerror = 0;
    ooberrors->avgrelerror = 0;
    if( mlpissoftmax(&ensemble->network, _state) )
    {
        ccnt = nin+1;
        pcnt = nin;
    }
    else
    {
        ccnt = nin+nout;
        pcnt = nin+nout;
    }
    ae_matrix_set_length(&xys, npoints, ccnt, _state);
    ae_vector_set_length(&s, npoints, _state);
    ae_matrix_set_length(&oobbuf, npoints, nout, _state);
    ae_vector_set_length(&oobcntbuf, npoints, _state);
    ae_vector_set_length(&x, nin, _state);
    ae_vector_set_length(&y, nout, _state);
    if( mlpissoftmax(&ensemble->network, _state) )
    {
        ae_vector_set_length(&dy, 1, _state);
    }
    else
    {
        ae_vector_set_length(&dy, nout, _state);
    }
    for(i=0; i<=npoints-1; i++)
    {
        for(j=0; j<=nout-1; j++)
        {
            oobbuf.ptr.pp_double[i][j] = 0;
        }
    }
    for(i=0; i<=npoints-1; i++)
    {
        oobcntbuf.ptr.p_int[i] = 0;
    }
    
    /*
     * main bagging cycle
     */
    for(k=0; k<=ensemble->ensemblesize-1; k++)
    {
        
        /*
         * prepare dataset
         */
        for(i=0; i<=npoints-1; i++)
        {
            s.ptr.p_bool[i] = ae_false;
        }
        for(i=0; i<=npoints-1; i++)
        {
            j = ae_randominteger(npoints, _state);
            s.ptr.p_bool[j] = ae_true;
            ae_v_move(&xys.ptr.pp_double[i][0], 1, &xy->ptr.pp_double[j][0], 1, ae_v_len(0,ccnt-1));
        }
        
        /*
         * train
         */
        if( lmalgorithm )
        {
            mlptrainlm(&ensemble->network, &xys, npoints, decay, restarts, info, &tmprep, _state);
        }
        else
        {
            mlptrainlbfgs(&ensemble->network, &xys, npoints, decay, restarts, wstep, maxits, info, &tmprep, _state);
        }
        if( *info<0 )
        {
            ae_frame_leave(_state);
            return;
        }
        
        /*
         * save results
         */
        rep->ngrad = rep->ngrad+tmprep.ngrad;
        rep->nhess = rep->nhess+tmprep.nhess;
        rep->ncholesky = rep->ncholesky+tmprep.ncholesky;
        ae_v_move(&ensemble->weights.ptr.p_double[k*wcount], 1, &ensemble->network.weights.ptr.p_double[0], 1, ae_v_len(k*wcount,(k+1)*wcount-1));
        ae_v_move(&ensemble->columnmeans.ptr.p_double[k*pcnt], 1, &ensemble->network.columnmeans.ptr.p_double[0], 1, ae_v_len(k*pcnt,(k+1)*pcnt-1));
        ae_v_move(&ensemble->columnsigmas.ptr.p_double[k*pcnt], 1, &ensemble->network.columnsigmas.ptr.p_double[0], 1, ae_v_len(k*pcnt,(k+1)*pcnt-1));
        
        /*
         * OOB estimates
         */
        for(i=0; i<=npoints-1; i++)
        {
            if( !s.ptr.p_bool[i] )
            {
                ae_v_move(&x.ptr.p_double[0], 1, &xy->ptr.pp_double[i][0], 1, ae_v_len(0,nin-1));
                mlpprocess(&ensemble->network, &x, &y, _state);
                ae_v_add(&oobbuf.ptr.pp_double[i][0], 1, &y.ptr.p_double[0], 1, ae_v_len(0,nout-1));
                oobcntbuf.ptr.p_int[i] = oobcntbuf.ptr.p_int[i]+1;
            }
        }
    }
    
    /*
     * OOB estimates
     */
    if( mlpissoftmax(&ensemble->network, _state) )
    {
        dserrallocate(nout, &dsbuf, _state);
    }
    else
    {
        dserrallocate(-nout, &dsbuf, _state);
    }
    for(i=0; i<=npoints-1; i++)
    {
        if( oobcntbuf.ptr.p_int[i]!=0 )
        {
            v = (double)1/(double)oobcntbuf.ptr.p_int[i];
            ae_v_moved(&y.ptr.p_double[0], 1, &oobbuf.ptr.pp_double[i][0], 1, ae_v_len(0,nout-1), v);
            if( mlpissoftmax(&ensemble->network, _state) )
            {
                dy.ptr.p_double[0] = xy->ptr.pp_double[i][nin];
            }
            else
            {
                ae_v_moved(&dy.ptr.p_double[0], 1, &xy->ptr.pp_double[i][nin], 1, ae_v_len(0,nout-1), v);
            }
            dserraccumulate(&dsbuf, &y, &dy, _state);
        }
    }
    dserrfinish(&dsbuf, _state);
    ooberrors->relclserror = dsbuf.ptr.p_double[0];
    ooberrors->avgce = dsbuf.ptr.p_double[1];
    ooberrors->rmserror = dsbuf.ptr.p_double[2];
    ooberrors->avgerror = dsbuf.ptr.p_double[3];
    ooberrors->avgrelerror = dsbuf.ptr.p_double[4];
    ae_frame_leave(_state);
}


ae_bool _mlpensemble_init(mlpensemble* p, ae_state *_state, ae_bool make_automatic)
{
    if( !ae_vector_init(&p->weights, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->columnmeans, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->columnsigmas, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    if( !_multilayerperceptron_init(&p->network, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init(&p->y, 0, DT_REAL, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


ae_bool _mlpensemble_init_copy(mlpensemble* dst, mlpensemble* src, ae_state *_state, ae_bool make_automatic)
{
    dst->ensemblesize = src->ensemblesize;
    if( !ae_vector_init_copy(&dst->weights, &src->weights, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->columnmeans, &src->columnmeans, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->columnsigmas, &src->columnsigmas, _state, make_automatic) )
        return ae_false;
    if( !_multilayerperceptron_init_copy(&dst->network, &src->network, _state, make_automatic) )
        return ae_false;
    if( !ae_vector_init_copy(&dst->y, &src->y, _state, make_automatic) )
        return ae_false;
    return ae_true;
}


void _mlpensemble_clear(mlpensemble* p)
{
    ae_vector_clear(&p->weights);
    ae_vector_clear(&p->columnmeans);
    ae_vector_clear(&p->columnsigmas);
    _multilayerperceptron_clear(&p->network);
    ae_vector_clear(&p->y);
}




/*************************************************************************
Principal components analysis

Subroutine  builds  orthogonal  basis  where  first  axis  corresponds  to
direction with maximum variance, second axis maximizes variance in subspace
orthogonal to first axis and so on.

It should be noted that, unlike LDA, PCA does not use class labels.

INPUT PARAMETERS:
    X           -   dataset, array[0..NPoints-1,0..NVars-1].
                    matrix contains ONLY INDEPENDENT VARIABLES.
    NPoints     -   dataset size, NPoints>=0
    NVars       -   number of independent variables, NVars>=1

¬€’ŒƒÕ€≈ œ¿–¿Ã≈“–€:
    Info        -   return code:
                    * -4, if SVD subroutine haven't converged
                    * -1, if wrong parameters has been passed (NPoints<0,
                          NVars<1)
                    *  1, if task is solved
    S2          -   array[0..NVars-1]. variance values corresponding
                    to basis vectors.
    V           -   array[0..NVars-1,0..NVars-1]
                    matrix, whose columns store basis vectors.

  -- ALGLIB --
     Copyright 25.08.2008 by Bochkanov Sergey
*************************************************************************/
void pcabuildbasis(/* Real    */ ae_matrix* x,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     /* Real    */ ae_vector* s2,
     /* Real    */ ae_matrix* v,
     ae_state *_state)
{
    ae_frame _frame_block;
    ae_matrix a;
    ae_matrix u;
    ae_matrix vt;
    ae_vector m;
    ae_vector t;
    ae_int_t i;
    ae_int_t j;
    double mean;
    double variance;
    double skewness;
    double kurtosis;

    ae_frame_make(_state, &_frame_block);
    *info = 0;
    ae_vector_clear(s2);
    ae_matrix_clear(v);
    ae_matrix_init(&a, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&u, 0, 0, DT_REAL, _state, ae_true);
    ae_matrix_init(&vt, 0, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&m, 0, DT_REAL, _state, ae_true);
    ae_vector_init(&t, 0, DT_REAL, _state, ae_true);

    
    /*
     * Check input data
     */
    if( npoints<0||nvars<1 )
    {
        *info = -1;
        ae_frame_leave(_state);
        return;
    }
    *info = 1;
    
    /*
     * Special case: NPoints=0
     */
    if( npoints==0 )
    {
        ae_vector_set_length(s2, nvars-1+1, _state);
        ae_matrix_set_length(v, nvars-1+1, nvars-1+1, _state);
        for(i=0; i<=nvars-1; i++)
        {
            s2->ptr.p_double[i] = 0;
        }
        for(i=0; i<=nvars-1; i++)
        {
            for(j=0; j<=nvars-1; j++)
            {
                if( i==j )
                {
                    v->ptr.pp_double[i][j] = 1;
                }
                else
                {
                    v->ptr.pp_double[i][j] = 0;
                }
            }
        }
        ae_frame_leave(_state);
        return;
    }
    
    /*
     * Calculate means
     */
    ae_vector_set_length(&m, nvars-1+1, _state);
    ae_vector_set_length(&t, npoints-1+1, _state);
    for(j=0; j<=nvars-1; j++)
    {
        ae_v_move(&t.ptr.p_double[0], 1, &x->ptr.pp_double[0][j], x->stride, ae_v_len(0,npoints-1));
        samplemoments(&t, npoints, &mean, &variance, &skewness, &kurtosis, _state);
        m.ptr.p_double[j] = mean;
    }
    
    /*
     * Center, apply SVD, prepare output
     */
    ae_matrix_set_length(&a, ae_maxint(npoints, nvars, _state)-1+1, nvars-1+1, _state);
    for(i=0; i<=npoints-1; i++)
    {
        ae_v_move(&a.ptr.pp_double[i][0], 1, &x->ptr.pp_double[i][0], 1, ae_v_len(0,nvars-1));
        ae_v_sub(&a.ptr.pp_double[i][0], 1, &m.ptr.p_double[0], 1, ae_v_len(0,nvars-1));
    }
    for(i=npoints; i<=nvars-1; i++)
    {
        for(j=0; j<=nvars-1; j++)
        {
            a.ptr.pp_double[i][j] = 0;
        }
    }
    if( !rmatrixsvd(&a, ae_maxint(npoints, nvars, _state), nvars, 0, 1, 2, s2, &u, &vt, _state) )
    {
        *info = -4;
        ae_frame_leave(_state);
        return;
    }
    if( npoints!=1 )
    {
        for(i=0; i<=nvars-1; i++)
        {
            s2->ptr.p_double[i] = ae_sqr(s2->ptr.p_double[i], _state)/(npoints-1);
        }
    }
    ae_matrix_set_length(v, nvars-1+1, nvars-1+1, _state);
    copyandtranspose(&vt, 0, nvars-1, 0, nvars-1, v, 0, nvars-1, 0, nvars-1, _state);
    ae_frame_leave(_state);
}



}

