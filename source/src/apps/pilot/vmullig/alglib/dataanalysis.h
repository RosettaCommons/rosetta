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
#ifndef _dataanalysis_pkg_h
#define _dataanalysis_pkg_h
#include "ap.h"
#include "alglibinternal.h"
#include "linalg.h"
#include "statistics.h"
#include "specialfunctions.h"
#include "alglibmisc.h"
#include "solvers.h"
#include "optimization.h"

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS COMPUTATIONAL CORE DECLARATIONS (DATATYPES)
//
/////////////////////////////////////////////////////////////////////////
namespace alglib_impl
{
typedef struct
{
    double relclserror;
    double avgce;
    double rmserror;
    double avgerror;
    double avgrelerror;
} cvreport;
typedef struct
{
    ae_int_t nvars;
    ae_int_t nclasses;
    ae_int_t ntrees;
    ae_int_t bufsize;
    ae_vector trees;
} decisionforest;
typedef struct
{
    double relclserror;
    double avgce;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double oobrelclserror;
    double oobavgce;
    double oobrmserror;
    double oobavgerror;
    double oobavgrelerror;
} dfreport;
typedef struct
{
    ae_vector treebuf;
    ae_vector idxbuf;
    ae_vector tmpbufr;
    ae_vector tmpbufr2;
    ae_vector tmpbufi;
    ae_vector classibuf;
    ae_vector sortrbuf;
    ae_vector sortrbuf2;
    ae_vector sortibuf;
    ae_vector varpool;
    ae_vector evsbin;
    ae_vector evssplits;
} dfinternalbuffers;
typedef struct
{
    ae_vector w;
} linearmodel;
typedef struct
{
    ae_matrix c;
    double rmserror;
    double avgerror;
    double avgrelerror;
    double cvrmserror;
    double cvavgerror;
    double cvavgrelerror;
    ae_int_t ncvdefects;
    ae_vector cvdefects;
} lrreport;
typedef struct
{
    ae_int_t hlnetworktype;
    ae_int_t hlnormtype;
    ae_vector hllayersizes;
    ae_vector hlconnections;
    ae_vector hlneurons;
    ae_vector structinfo;
    ae_vector weights;
    ae_vector columnmeans;
    ae_vector columnsigmas;
    ae_vector neurons;
    ae_vector dfdnet;
    ae_vector derror;
    ae_vector x;
    ae_vector y;
    ae_matrix chunks;
    ae_vector nwbuf;
    ae_vector integerbuf;
} multilayerperceptron;
typedef struct
{
    ae_vector w;
} logitmodel;
typedef struct
{
    ae_bool brackt;
    ae_bool stage1;
    ae_int_t infoc;
    double dg;
    double dgm;
    double dginit;
    double dgtest;
    double dgx;
    double dgxm;
    double dgy;
    double dgym;
    double finit;
    double ftest1;
    double fm;
    double fx;
    double fxm;
    double fy;
    double fym;
    double stx;
    double sty;
    double stmin;
    double stmax;
    double width;
    double width1;
    double xtrapf;
} logitmcstate;
typedef struct
{
    ae_int_t ngrad;
    ae_int_t nhess;
} mnlreport;
typedef struct
{
    ae_int_t n;
    ae_vector states;
    ae_int_t npairs;
    ae_matrix data;
    ae_matrix ec;
    ae_matrix bndl;
    ae_matrix bndu;
    ae_matrix c;
    ae_vector ct;
    ae_int_t ccnt;
    ae_vector pw;
    ae_matrix priorp;
    double regterm;
    minbleicstate bs;
    ae_int_t repinneriterationscount;
    ae_int_t repouteriterationscount;
    ae_int_t repnfev;
    ae_int_t repterminationtype;
    minbleicreport br;
    ae_vector tmpp;
    ae_vector effectivew;
    ae_vector effectivebndl;
    ae_vector effectivebndu;
    ae_matrix effectivec;
    ae_vector effectivect;
    ae_vector h;
    ae_matrix p;
} mcpdstate;
typedef struct
{
    ae_int_t inneriterationscount;
    ae_int_t outeriterationscount;
    ae_int_t nfev;
    ae_int_t terminationtype;
} mcpdreport;
typedef struct
{
    ae_int_t ngrad;
    ae_int_t nhess;
    ae_int_t ncholesky;
} mlpreport;
typedef struct
{
    double relclserror;
    double avgce;
    double rmserror;
    double avgerror;
    double avgrelerror;
} mlpcvreport;
typedef struct
{
    ae_int_t ensemblesize;
    ae_vector weights;
    ae_vector columnmeans;
    ae_vector columnsigmas;
    multilayerperceptron network;
    ae_vector y;
} mlpensemble;

}

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS C++ INTERFACE
//
/////////////////////////////////////////////////////////////////////////
namespace alglib
{



/*************************************************************************

*************************************************************************/
class _decisionforest_owner
{
public:
    _decisionforest_owner();
    _decisionforest_owner(const _decisionforest_owner &rhs);
    _decisionforest_owner& operator=(const _decisionforest_owner &rhs);
    virtual ~_decisionforest_owner();
    alglib_impl::decisionforest* c_ptr();
    alglib_impl::decisionforest* c_ptr() const;
protected:
    alglib_impl::decisionforest *p_struct;
};
class decisionforest : public _decisionforest_owner
{
public:
    decisionforest();
    decisionforest(const decisionforest &rhs);
    decisionforest& operator=(const decisionforest &rhs);
    virtual ~decisionforest();

};


/*************************************************************************

*************************************************************************/
class _dfreport_owner
{
public:
    _dfreport_owner();
    _dfreport_owner(const _dfreport_owner &rhs);
    _dfreport_owner& operator=(const _dfreport_owner &rhs);
    virtual ~_dfreport_owner();
    alglib_impl::dfreport* c_ptr();
    alglib_impl::dfreport* c_ptr() const;
protected:
    alglib_impl::dfreport *p_struct;
};
class dfreport : public _dfreport_owner
{
public:
    dfreport();
    dfreport(const dfreport &rhs);
    dfreport& operator=(const dfreport &rhs);
    virtual ~dfreport();
    double &relclserror;
    double &avgce;
    double &rmserror;
    double &avgerror;
    double &avgrelerror;
    double &oobrelclserror;
    double &oobavgce;
    double &oobrmserror;
    double &oobavgerror;
    double &oobavgrelerror;

};

/*************************************************************************

*************************************************************************/
class _linearmodel_owner
{
public:
    _linearmodel_owner();
    _linearmodel_owner(const _linearmodel_owner &rhs);
    _linearmodel_owner& operator=(const _linearmodel_owner &rhs);
    virtual ~_linearmodel_owner();
    alglib_impl::linearmodel* c_ptr();
    alglib_impl::linearmodel* c_ptr() const;
protected:
    alglib_impl::linearmodel *p_struct;
};
class linearmodel : public _linearmodel_owner
{
public:
    linearmodel();
    linearmodel(const linearmodel &rhs);
    linearmodel& operator=(const linearmodel &rhs);
    virtual ~linearmodel();

};


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
class _lrreport_owner
{
public:
    _lrreport_owner();
    _lrreport_owner(const _lrreport_owner &rhs);
    _lrreport_owner& operator=(const _lrreport_owner &rhs);
    virtual ~_lrreport_owner();
    alglib_impl::lrreport* c_ptr();
    alglib_impl::lrreport* c_ptr() const;
protected:
    alglib_impl::lrreport *p_struct;
};
class lrreport : public _lrreport_owner
{
public:
    lrreport();
    lrreport(const lrreport &rhs);
    lrreport& operator=(const lrreport &rhs);
    virtual ~lrreport();
    real_2d_array c;
    double &rmserror;
    double &avgerror;
    double &avgrelerror;
    double &cvrmserror;
    double &cvavgerror;
    double &cvavgrelerror;
    ae_int_t &ncvdefects;
    integer_1d_array cvdefects;

};







/*************************************************************************

*************************************************************************/
class _multilayerperceptron_owner
{
public:
    _multilayerperceptron_owner();
    _multilayerperceptron_owner(const _multilayerperceptron_owner &rhs);
    _multilayerperceptron_owner& operator=(const _multilayerperceptron_owner &rhs);
    virtual ~_multilayerperceptron_owner();
    alglib_impl::multilayerperceptron* c_ptr();
    alglib_impl::multilayerperceptron* c_ptr() const;
protected:
    alglib_impl::multilayerperceptron *p_struct;
};
class multilayerperceptron : public _multilayerperceptron_owner
{
public:
    multilayerperceptron();
    multilayerperceptron(const multilayerperceptron &rhs);
    multilayerperceptron& operator=(const multilayerperceptron &rhs);
    virtual ~multilayerperceptron();

};

/*************************************************************************

*************************************************************************/
class _logitmodel_owner
{
public:
    _logitmodel_owner();
    _logitmodel_owner(const _logitmodel_owner &rhs);
    _logitmodel_owner& operator=(const _logitmodel_owner &rhs);
    virtual ~_logitmodel_owner();
    alglib_impl::logitmodel* c_ptr();
    alglib_impl::logitmodel* c_ptr() const;
protected:
    alglib_impl::logitmodel *p_struct;
};
class logitmodel : public _logitmodel_owner
{
public:
    logitmodel();
    logitmodel(const logitmodel &rhs);
    logitmodel& operator=(const logitmodel &rhs);
    virtual ~logitmodel();

};


/*************************************************************************
MNLReport structure contains information about training process:
* NGrad     -   number of gradient calculations
* NHess     -   number of Hessian calculations
*************************************************************************/
class _mnlreport_owner
{
public:
    _mnlreport_owner();
    _mnlreport_owner(const _mnlreport_owner &rhs);
    _mnlreport_owner& operator=(const _mnlreport_owner &rhs);
    virtual ~_mnlreport_owner();
    alglib_impl::mnlreport* c_ptr();
    alglib_impl::mnlreport* c_ptr() const;
protected:
    alglib_impl::mnlreport *p_struct;
};
class mnlreport : public _mnlreport_owner
{
public:
    mnlreport();
    mnlreport(const mnlreport &rhs);
    mnlreport& operator=(const mnlreport &rhs);
    virtual ~mnlreport();
    ae_int_t &ngrad;
    ae_int_t &nhess;

};

/*************************************************************************
This structure is a MCPD (Markov Chains for Population Data) solver.

You should use ALGLIB functions in order to work with this object.

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
class _mcpdstate_owner
{
public:
    _mcpdstate_owner();
    _mcpdstate_owner(const _mcpdstate_owner &rhs);
    _mcpdstate_owner& operator=(const _mcpdstate_owner &rhs);
    virtual ~_mcpdstate_owner();
    alglib_impl::mcpdstate* c_ptr();
    alglib_impl::mcpdstate* c_ptr() const;
protected:
    alglib_impl::mcpdstate *p_struct;
};
class mcpdstate : public _mcpdstate_owner
{
public:
    mcpdstate();
    mcpdstate(const mcpdstate &rhs);
    mcpdstate& operator=(const mcpdstate &rhs);
    virtual ~mcpdstate();

};


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
class _mcpdreport_owner
{
public:
    _mcpdreport_owner();
    _mcpdreport_owner(const _mcpdreport_owner &rhs);
    _mcpdreport_owner& operator=(const _mcpdreport_owner &rhs);
    virtual ~_mcpdreport_owner();
    alglib_impl::mcpdreport* c_ptr();
    alglib_impl::mcpdreport* c_ptr() const;
protected:
    alglib_impl::mcpdreport *p_struct;
};
class mcpdreport : public _mcpdreport_owner
{
public:
    mcpdreport();
    mcpdreport(const mcpdreport &rhs);
    mcpdreport& operator=(const mcpdreport &rhs);
    virtual ~mcpdreport();
    ae_int_t &inneriterationscount;
    ae_int_t &outeriterationscount;
    ae_int_t &nfev;
    ae_int_t &terminationtype;

};

/*************************************************************************
Training report:
    * NGrad     - number of gradient calculations
    * NHess     - number of Hessian calculations
    * NCholesky - number of Cholesky decompositions
*************************************************************************/
class _mlpreport_owner
{
public:
    _mlpreport_owner();
    _mlpreport_owner(const _mlpreport_owner &rhs);
    _mlpreport_owner& operator=(const _mlpreport_owner &rhs);
    virtual ~_mlpreport_owner();
    alglib_impl::mlpreport* c_ptr();
    alglib_impl::mlpreport* c_ptr() const;
protected:
    alglib_impl::mlpreport *p_struct;
};
class mlpreport : public _mlpreport_owner
{
public:
    mlpreport();
    mlpreport(const mlpreport &rhs);
    mlpreport& operator=(const mlpreport &rhs);
    virtual ~mlpreport();
    ae_int_t &ngrad;
    ae_int_t &nhess;
    ae_int_t &ncholesky;

};


/*************************************************************************
Cross-validation estimates of generalization error
*************************************************************************/
class _mlpcvreport_owner
{
public:
    _mlpcvreport_owner();
    _mlpcvreport_owner(const _mlpcvreport_owner &rhs);
    _mlpcvreport_owner& operator=(const _mlpcvreport_owner &rhs);
    virtual ~_mlpcvreport_owner();
    alglib_impl::mlpcvreport* c_ptr();
    alglib_impl::mlpcvreport* c_ptr() const;
protected:
    alglib_impl::mlpcvreport *p_struct;
};
class mlpcvreport : public _mlpcvreport_owner
{
public:
    mlpcvreport();
    mlpcvreport(const mlpcvreport &rhs);
    mlpcvreport& operator=(const mlpcvreport &rhs);
    virtual ~mlpcvreport();
    double &relclserror;
    double &avgce;
    double &rmserror;
    double &avgerror;
    double &avgrelerror;

};

/*************************************************************************
Neural networks ensemble
*************************************************************************/
class _mlpensemble_owner
{
public:
    _mlpensemble_owner();
    _mlpensemble_owner(const _mlpensemble_owner &rhs);
    _mlpensemble_owner& operator=(const _mlpensemble_owner &rhs);
    virtual ~_mlpensemble_owner();
    alglib_impl::mlpensemble* c_ptr();
    alglib_impl::mlpensemble* c_ptr() const;
protected:
    alglib_impl::mlpensemble *p_struct;
};
class mlpensemble : public _mlpensemble_owner
{
public:
    mlpensemble();
    mlpensemble(const mlpensemble &rhs);
    mlpensemble& operator=(const mlpensemble &rhs);
    virtual ~mlpensemble();

};

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
void dsoptimalsplit2(const real_1d_array &a, const integer_1d_array &c, const ae_int_t n, ae_int_t &info, double &threshold, double &pal, double &pbl, double &par, double &pbr, double &cve);


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
void dsoptimalsplit2fast(real_1d_array &a, integer_1d_array &c, integer_1d_array &tiesbuf, integer_1d_array &cntbuf, real_1d_array &bufr, integer_1d_array &bufi, const ae_int_t n, const ae_int_t nc, const double alpha, ae_int_t &info, double &threshold, double &rms, double &cvrms);

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
void dfserialize(decisionforest &obj, std::string &s_out);


/*************************************************************************
This function unserializes data structure from string.
*************************************************************************/
void dfunserialize(std::string &s_in, decisionforest &obj);


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
void dfbuildrandomdecisionforest(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, const ae_int_t ntrees, const double r, ae_int_t &info, decisionforest &df, dfreport &rep);


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
void dfbuildrandomdecisionforestx1(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, const ae_int_t ntrees, const ae_int_t nrndvars, const double r, ae_int_t &info, decisionforest &df, dfreport &rep);


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
void dfprocess(const decisionforest &df, const real_1d_array &x, real_1d_array &y);


/*************************************************************************
'interactive' variant of DFProcess for languages like Python which support
constructs like "Y = DFProcessI(DF,X)" and interactive mode of interpreter

This function allocates new array on each call,  so  it  is  significantly
slower than its 'non-interactive' counterpart, but it is  more  convenient
when you call it from command line.

  -- ALGLIB --
     Copyright 28.02.2010 by Bochkanov Sergey
*************************************************************************/
void dfprocessi(const decisionforest &df, const real_1d_array &x, real_1d_array &y);


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
double dfrelclserror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);


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
double dfavgce(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);


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
double dfrmserror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);


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
double dfavgerror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);


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
double dfavgrelerror(const decisionforest &df, const real_2d_array &xy, const ae_int_t npoints);

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
void lrbuild(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar);


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
void lrbuilds(const real_2d_array &xy, const real_1d_array &s, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar);


/*************************************************************************
Like LRBuildS, but builds model

    Y = A(0)*X[0] + ... + A(N-1)*X[N-1]

i.e. with zero constant term.

  -- ALGLIB --
     Copyright 30.10.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuildzs(const real_2d_array &xy, const real_1d_array &s, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar);


/*************************************************************************
Like LRBuild but builds model

    Y = A(0)*X[0] + ... + A(N-1)*X[N-1]

i.e. with zero constant term.

  -- ALGLIB --
     Copyright 30.10.2008 by Bochkanov Sergey
*************************************************************************/
void lrbuildz(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, linearmodel &lm, lrreport &ar);


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
void lrunpack(const linearmodel &lm, real_1d_array &v, ae_int_t &nvars);


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
void lrpack(const real_1d_array &v, const ae_int_t nvars, linearmodel &lm);


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
double lrprocess(const linearmodel &lm, const real_1d_array &x);


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
double lrrmserror(const linearmodel &lm, const real_2d_array &xy, const ae_int_t npoints);


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
double lravgerror(const linearmodel &lm, const real_2d_array &xy, const ae_int_t npoints);


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
double lravgrelerror(const linearmodel &lm, const real_2d_array &xy, const ae_int_t npoints);

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
void filtersma(real_1d_array &x, const ae_int_t n, const ae_int_t k);
void filtersma(real_1d_array &x, const ae_int_t k);


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
void filterema(real_1d_array &x, const ae_int_t n, const double alpha);
void filterema(real_1d_array &x, const double alpha);


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
void filterlrma(real_1d_array &x, const ae_int_t n, const ae_int_t k);
void filterlrma(real_1d_array &x, const ae_int_t k);

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
void kmeansgenerate(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t k, const ae_int_t restarts, ae_int_t &info, real_2d_array &c, integer_1d_array &xyc);

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
void fisherlda(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, real_1d_array &w);


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
void fisherldan(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, real_2d_array &w);

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
void mlpserialize(multilayerperceptron &obj, std::string &s_out);


/*************************************************************************
This function unserializes data structure from string.
*************************************************************************/
void mlpunserialize(std::string &s_in, multilayerperceptron &obj);


/*************************************************************************
Creates  neural  network  with  NIn  inputs,  NOut outputs, without hidden
layers, with linear output layer. Network weights are  filled  with  small
random values.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreate0(const ae_int_t nin, const ae_int_t nout, multilayerperceptron &network);


/*************************************************************************
Same  as  MLPCreate0,  but  with  one  hidden  layer  (NHid  neurons) with
non-linear activation function. Output layer is linear.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreate1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, multilayerperceptron &network);


/*************************************************************************
Same as MLPCreate0, but with two hidden layers (NHid1 and  NHid2  neurons)
with non-linear activation function. Output layer is linear.
 $ALL

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreate2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, multilayerperceptron &network);


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
void mlpcreateb0(const ae_int_t nin, const ae_int_t nout, const double b, const double d, multilayerperceptron &network);


/*************************************************************************
Same as MLPCreateB0 but with non-linear hidden layer.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreateb1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double b, const double d, multilayerperceptron &network);


/*************************************************************************
Same as MLPCreateB0 but with two non-linear hidden layers.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreateb2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double b, const double d, multilayerperceptron &network);


/*************************************************************************
Creates  neural  network  with  NIn  inputs,  NOut outputs, without hidden
layers with non-linear output layer. Network weights are filled with small
random values. Activation function of the output layer takes values [A,B].

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreater0(const ae_int_t nin, const ae_int_t nout, const double a, const double b, multilayerperceptron &network);


/*************************************************************************
Same as MLPCreateR0, but with non-linear hidden layer.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreater1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double a, const double b, multilayerperceptron &network);


/*************************************************************************
Same as MLPCreateR0, but with two non-linear hidden layers.

  -- ALGLIB --
     Copyright 30.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlpcreater2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double a, const double b, multilayerperceptron &network);


/*************************************************************************
Creates classifier network with NIn  inputs  and  NOut  possible  classes.
Network contains no hidden layers and linear output  layer  with  SOFTMAX-
normalization  (so  outputs  sums  up  to  1.0  and  converge to posterior
probabilities).

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreatec0(const ae_int_t nin, const ae_int_t nout, multilayerperceptron &network);


/*************************************************************************
Same as MLPCreateC0, but with one non-linear hidden layer.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreatec1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, multilayerperceptron &network);


/*************************************************************************
Same as MLPCreateC0, but with two non-linear hidden layers.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpcreatec2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, multilayerperceptron &network);


/*************************************************************************
Randomization of neural network weights

  -- ALGLIB --
     Copyright 06.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlprandomize(const multilayerperceptron &network);


/*************************************************************************
Randomization of neural network weights and standartisator

  -- ALGLIB --
     Copyright 10.03.2008 by Bochkanov Sergey
*************************************************************************/
void mlprandomizefull(const multilayerperceptron &network);


/*************************************************************************
Returns information about initialized network: number of inputs, outputs,
weights.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
void mlpproperties(const multilayerperceptron &network, ae_int_t &nin, ae_int_t &nout, ae_int_t &wcount);


/*************************************************************************
Returns number of inputs.

  -- ALGLIB --
     Copyright 19.10.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpgetinputscount(const multilayerperceptron &network);


/*************************************************************************
Returns number of outputs.

  -- ALGLIB --
     Copyright 19.10.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpgetoutputscount(const multilayerperceptron &network);


/*************************************************************************
Returns number of weights.

  -- ALGLIB --
     Copyright 19.10.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpgetweightscount(const multilayerperceptron &network);


/*************************************************************************
Tells whether network is SOFTMAX-normalized (i.e. classifier) or not.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
bool mlpissoftmax(const multilayerperceptron &network);


/*************************************************************************
This function returns total number of layers (including input, hidden and
output layers).

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpgetlayerscount(const multilayerperceptron &network);


/*************************************************************************
This function returns size of K-th layer.

K=0 corresponds to input layer, K=CNT-1 corresponds to output layer.

Size of the output layer is always equal to the number of outputs, although
when we have softmax-normalized network, last neuron doesn't have any
connections - it is just zero.

  -- ALGLIB --
     Copyright 25.03.2011 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpgetlayersize(const multilayerperceptron &network, const ae_int_t k);


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
void mlpgetinputscaling(const multilayerperceptron &network, const ae_int_t i, double &mean, double &sigma);


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
void mlpgetoutputscaling(const multilayerperceptron &network, const ae_int_t i, double &mean, double &sigma);


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
void mlpgetneuroninfo(const multilayerperceptron &network, const ae_int_t k, const ae_int_t i, ae_int_t &fkind, double &threshold);


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
double mlpgetweight(const multilayerperceptron &network, const ae_int_t k0, const ae_int_t i0, const ae_int_t k1, const ae_int_t i1);


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
void mlpsetinputscaling(const multilayerperceptron &network, const ae_int_t i, const double mean, const double sigma);


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
void mlpsetoutputscaling(const multilayerperceptron &network, const ae_int_t i, const double mean, const double sigma);


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
void mlpsetneuroninfo(const multilayerperceptron &network, const ae_int_t k, const ae_int_t i, const ae_int_t fkind, const double threshold);


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
void mlpsetweight(const multilayerperceptron &network, const ae_int_t k0, const ae_int_t i0, const ae_int_t k1, const ae_int_t i1, const double w);


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
void mlpactivationfunction(const double net, const ae_int_t k, double &f, double &df, double &d2f);


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
void mlpprocess(const multilayerperceptron &network, const real_1d_array &x, real_1d_array &y);


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
void mlpprocessi(const multilayerperceptron &network, const real_1d_array &x, real_1d_array &y);


/*************************************************************************
Error function for neural network, internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
double mlperror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize);


/*************************************************************************
Natural error function for neural network, internal subroutine.

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
double mlperrorn(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize);


/*************************************************************************
Classification error

  -- ALGLIB --
     Copyright 04.11.2007 by Bochkanov Sergey
*************************************************************************/
ae_int_t mlpclserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize);


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
double mlprelclserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);


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
double mlpavgce(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);


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
double mlprmserror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);


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
double mlpavgerror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);


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
double mlpavgrelerror(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints);


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
void mlpgrad(const multilayerperceptron &network, const real_1d_array &x, const real_1d_array &desiredy, double &e, real_1d_array &grad);


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
void mlpgradn(const multilayerperceptron &network, const real_1d_array &x, const real_1d_array &desiredy, double &e, real_1d_array &grad);


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
void mlpgradbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad);


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
void mlpgradnbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad);


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
void mlphessiannbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad, real_2d_array &h);


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
void mlphessianbatch(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t ssize, double &e, real_1d_array &grad, real_2d_array &h);

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
void mnltrainh(const real_2d_array &xy, const ae_int_t npoints, const ae_int_t nvars, const ae_int_t nclasses, ae_int_t &info, logitmodel &lm, mnlreport &rep);


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
void mnlprocess(const logitmodel &lm, const real_1d_array &x, real_1d_array &y);


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
void mnlprocessi(const logitmodel &lm, const real_1d_array &x, real_1d_array &y);


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
void mnlunpack(const logitmodel &lm, real_2d_array &a, ae_int_t &nvars, ae_int_t &nclasses);


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
void mnlpack(const real_2d_array &a, const ae_int_t nvars, const ae_int_t nclasses, logitmodel &lm);


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
double mnlavgce(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);


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
double mnlrelclserror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);


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
double mnlrmserror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);


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
double mnlavgerror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);


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
double mnlavgrelerror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t ssize);


/*************************************************************************
Classification error on test set = MNLRelClsError*NPoints

  -- ALGLIB --
     Copyright 10.09.2008 by Bochkanov Sergey
*************************************************************************/
ae_int_t mnlclserror(const logitmodel &lm, const real_2d_array &xy, const ae_int_t npoints);

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
void mcpdcreate(const ae_int_t n, mcpdstate &s);


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
void mcpdcreateentry(const ae_int_t n, const ae_int_t entrystate, mcpdstate &s);


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
void mcpdcreateexit(const ae_int_t n, const ae_int_t exitstate, mcpdstate &s);


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
void mcpdcreateentryexit(const ae_int_t n, const ae_int_t entrystate, const ae_int_t exitstate, mcpdstate &s);


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
void mcpdaddtrack(const mcpdstate &s, const real_2d_array &xy, const ae_int_t k);
void mcpdaddtrack(const mcpdstate &s, const real_2d_array &xy);


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
void mcpdsetec(const mcpdstate &s, const real_2d_array &ec);


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
void mcpdaddec(const mcpdstate &s, const ae_int_t i, const ae_int_t j, const double c);


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
void mcpdsetbc(const mcpdstate &s, const real_2d_array &bndl, const real_2d_array &bndu);


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
void mcpdaddbc(const mcpdstate &s, const ae_int_t i, const ae_int_t j, const double bndl, const double bndu);


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
void mcpdsetlc(const mcpdstate &s, const real_2d_array &c, const integer_1d_array &ct, const ae_int_t k);
void mcpdsetlc(const mcpdstate &s, const real_2d_array &c, const integer_1d_array &ct);


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
void mcpdsettikhonovregularizer(const mcpdstate &s, const double v);


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
void mcpdsetprior(const mcpdstate &s, const real_2d_array &pp);


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
void mcpdsetpredictionweights(const mcpdstate &s, const real_1d_array &pw);


/*************************************************************************
This function is used to start solution of the MCPD problem.

After return from this function, you can use MCPDResults() to get solution
and completion code.

  -- ALGLIB --
     Copyright 23.05.2010 by Bochkanov Sergey
*************************************************************************/
void mcpdsolve(const mcpdstate &s);


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
void mcpdresults(const mcpdstate &s, real_2d_array &p, mcpdreport &rep);

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
void mlptrainlm(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep);


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
void mlptrainlbfgs(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const double wstep, const ae_int_t maxits, ae_int_t &info, mlpreport &rep);


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
void mlptraines(const multilayerperceptron &network, const real_2d_array &trnxy, const ae_int_t trnsize, const real_2d_array &valxy, const ae_int_t valsize, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep);


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
void mlpkfoldcvlbfgs(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const double wstep, const ae_int_t maxits, const ae_int_t foldscount, ae_int_t &info, mlpreport &rep, mlpcvreport &cvrep);


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
void mlpkfoldcvlm(const multilayerperceptron &network, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const ae_int_t foldscount, ae_int_t &info, mlpreport &rep, mlpcvreport &cvrep);

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
void mlpeserialize(mlpensemble &obj, std::string &s_out);


/*************************************************************************
This function unserializes data structure from string.
*************************************************************************/
void mlpeunserialize(std::string &s_in, mlpensemble &obj);


/*************************************************************************
Like MLPCreate0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreate0(const ae_int_t nin, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreate1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreate1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreate2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreate2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreateB0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreateb0(const ae_int_t nin, const ae_int_t nout, const double b, const double d, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreateB1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreateb1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double b, const double d, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreateB2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreateb2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double b, const double d, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreateR0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreater0(const ae_int_t nin, const ae_int_t nout, const double a, const double b, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreateR1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreater1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const double a, const double b, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreateR2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreater2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const double a, const double b, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreateC0, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatec0(const ae_int_t nin, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreateC1, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatec1(const ae_int_t nin, const ae_int_t nhid, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Like MLPCreateC2, but for ensembles.

  -- ALGLIB --
     Copyright 18.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatec2(const ae_int_t nin, const ae_int_t nhid1, const ae_int_t nhid2, const ae_int_t nout, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Creates ensemble from network. Only network geometry is copied.

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpecreatefromnetwork(const multilayerperceptron &network, const ae_int_t ensemblesize, mlpensemble &ensemble);


/*************************************************************************
Randomization of MLP ensemble

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlperandomize(const mlpensemble &ensemble);


/*************************************************************************
Return ensemble properties (number of inputs and outputs).

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
void mlpeproperties(const mlpensemble &ensemble, ae_int_t &nin, ae_int_t &nout);


/*************************************************************************
Return normalization type (whether ensemble is SOFTMAX-normalized or not).

  -- ALGLIB --
     Copyright 17.02.2009 by Bochkanov Sergey
*************************************************************************/
bool mlpeissoftmax(const mlpensemble &ensemble);


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
void mlpeprocess(const mlpensemble &ensemble, const real_1d_array &x, real_1d_array &y);


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
void mlpeprocessi(const mlpensemble &ensemble, const real_1d_array &x, real_1d_array &y);


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
double mlperelclserror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);


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
double mlpeavgce(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);


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
double mlpermserror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);


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
double mlpeavgerror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);


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
double mlpeavgrelerror(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints);


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
void mlpebagginglm(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep, mlpcvreport &ooberrors);


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
void mlpebagginglbfgs(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, const double wstep, const ae_int_t maxits, ae_int_t &info, mlpreport &rep, mlpcvreport &ooberrors);


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
void mlpetraines(const mlpensemble &ensemble, const real_2d_array &xy, const ae_int_t npoints, const double decay, const ae_int_t restarts, ae_int_t &info, mlpreport &rep);

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
void pcabuildbasis(const real_2d_array &x, const ae_int_t npoints, const ae_int_t nvars, ae_int_t &info, real_1d_array &s2, real_2d_array &v);
}

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS COMPUTATIONAL CORE DECLARATIONS (FUNCTIONS)
//
/////////////////////////////////////////////////////////////////////////
namespace alglib_impl
{
void dserrallocate(ae_int_t nclasses,
     /* Real    */ ae_vector* buf,
     ae_state *_state);
void dserraccumulate(/* Real    */ ae_vector* buf,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* desiredy,
     ae_state *_state);
void dserrfinish(/* Real    */ ae_vector* buf, ae_state *_state);
void dsnormalize(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     /* Real    */ ae_vector* means,
     /* Real    */ ae_vector* sigmas,
     ae_state *_state);
void dsnormalizec(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     /* Real    */ ae_vector* means,
     /* Real    */ ae_vector* sigmas,
     ae_state *_state);
double dsgetmeanmindistance(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_state *_state);
void dstie(/* Real    */ ae_vector* a,
     ae_int_t n,
     /* Integer */ ae_vector* ties,
     ae_int_t* tiecount,
     /* Integer */ ae_vector* p1,
     /* Integer */ ae_vector* p2,
     ae_state *_state);
void dstiefasti(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     ae_int_t n,
     /* Integer */ ae_vector* ties,
     ae_int_t* tiecount,
     /* Real    */ ae_vector* bufr,
     /* Integer */ ae_vector* bufi,
     ae_state *_state);
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
     ae_state *_state);
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
     ae_state *_state);
void dssplitk(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* c,
     ae_int_t n,
     ae_int_t nc,
     ae_int_t kmax,
     ae_int_t* info,
     /* Real    */ ae_vector* thresholds,
     ae_int_t* ni,
     double* cve,
     ae_state *_state);
void dsoptimalsplitk(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* c,
     ae_int_t n,
     ae_int_t nc,
     ae_int_t kmax,
     ae_int_t* info,
     /* Real    */ ae_vector* thresholds,
     ae_int_t* ni,
     double* cve,
     ae_state *_state);
ae_bool _cvreport_init(cvreport* p, ae_state *_state, ae_bool make_automatic);
ae_bool _cvreport_init_copy(cvreport* dst, cvreport* src, ae_state *_state, ae_bool make_automatic);
void _cvreport_clear(cvreport* p);
void dfbuildrandomdecisionforest(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t ntrees,
     double r,
     ae_int_t* info,
     decisionforest* df,
     dfreport* rep,
     ae_state *_state);
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
     ae_state *_state);
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
     ae_state *_state);
void dfprocess(decisionforest* df,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
void dfprocessi(decisionforest* df,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
double dfrelclserror(decisionforest* df,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double dfavgce(decisionforest* df,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double dfrmserror(decisionforest* df,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double dfavgerror(decisionforest* df,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double dfavgrelerror(decisionforest* df,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
void dfcopy(decisionforest* df1, decisionforest* df2, ae_state *_state);
void dfalloc(ae_serializer* s, decisionforest* forest, ae_state *_state);
void dfserialize(ae_serializer* s,
     decisionforest* forest,
     ae_state *_state);
void dfunserialize(ae_serializer* s,
     decisionforest* forest,
     ae_state *_state);
ae_bool _decisionforest_init(decisionforest* p, ae_state *_state, ae_bool make_automatic);
ae_bool _decisionforest_init_copy(decisionforest* dst, decisionforest* src, ae_state *_state, ae_bool make_automatic);
void _decisionforest_clear(decisionforest* p);
ae_bool _dfreport_init(dfreport* p, ae_state *_state, ae_bool make_automatic);
ae_bool _dfreport_init_copy(dfreport* dst, dfreport* src, ae_state *_state, ae_bool make_automatic);
void _dfreport_clear(dfreport* p);
ae_bool _dfinternalbuffers_init(dfinternalbuffers* p, ae_state *_state, ae_bool make_automatic);
ae_bool _dfinternalbuffers_init_copy(dfinternalbuffers* dst, dfinternalbuffers* src, ae_state *_state, ae_bool make_automatic);
void _dfinternalbuffers_clear(dfinternalbuffers* p);
void lrbuild(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     linearmodel* lm,
     lrreport* ar,
     ae_state *_state);
void lrbuilds(/* Real    */ ae_matrix* xy,
     /* Real    */ ae_vector* s,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     linearmodel* lm,
     lrreport* ar,
     ae_state *_state);
void lrbuildzs(/* Real    */ ae_matrix* xy,
     /* Real    */ ae_vector* s,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     linearmodel* lm,
     lrreport* ar,
     ae_state *_state);
void lrbuildz(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     linearmodel* lm,
     lrreport* ar,
     ae_state *_state);
void lrunpack(linearmodel* lm,
     /* Real    */ ae_vector* v,
     ae_int_t* nvars,
     ae_state *_state);
void lrpack(/* Real    */ ae_vector* v,
     ae_int_t nvars,
     linearmodel* lm,
     ae_state *_state);
double lrprocess(linearmodel* lm,
     /* Real    */ ae_vector* x,
     ae_state *_state);
double lrrmserror(linearmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double lravgerror(linearmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double lravgrelerror(linearmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
void lrcopy(linearmodel* lm1, linearmodel* lm2, ae_state *_state);
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
     ae_state *_state);
void lrline(/* Real    */ ae_matrix* xy,
     ae_int_t n,
     ae_int_t* info,
     double* a,
     double* b,
     ae_state *_state);
ae_bool _linearmodel_init(linearmodel* p, ae_state *_state, ae_bool make_automatic);
ae_bool _linearmodel_init_copy(linearmodel* dst, linearmodel* src, ae_state *_state, ae_bool make_automatic);
void _linearmodel_clear(linearmodel* p);
ae_bool _lrreport_init(lrreport* p, ae_state *_state, ae_bool make_automatic);
ae_bool _lrreport_init_copy(lrreport* dst, lrreport* src, ae_state *_state, ae_bool make_automatic);
void _lrreport_clear(lrreport* p);
void filtersma(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_int_t k,
     ae_state *_state);
void filterema(/* Real    */ ae_vector* x,
     ae_int_t n,
     double alpha,
     ae_state *_state);
void filterlrma(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_int_t k,
     ae_state *_state);
void kmeansgenerate(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t k,
     ae_int_t restarts,
     ae_int_t* info,
     /* Real    */ ae_matrix* c,
     /* Integer */ ae_vector* xyc,
     ae_state *_state);
void fisherlda(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t* info,
     /* Real    */ ae_vector* w,
     ae_state *_state);
void fisherldan(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t* info,
     /* Real    */ ae_matrix* w,
     ae_state *_state);
void mlpcreate0(ae_int_t nin,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreate1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreate2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreateb0(ae_int_t nin,
     ae_int_t nout,
     double b,
     double d,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreateb1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     double b,
     double d,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreateb2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     double b,
     double d,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreater0(ae_int_t nin,
     ae_int_t nout,
     double a,
     double b,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreater1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     double a,
     double b,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreater2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     double a,
     double b,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreatec0(ae_int_t nin,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreatec1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcreatec2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     multilayerperceptron* network,
     ae_state *_state);
void mlpcopy(multilayerperceptron* network1,
     multilayerperceptron* network2,
     ae_state *_state);
void mlpserializeold(multilayerperceptron* network,
     /* Real    */ ae_vector* ra,
     ae_int_t* rlen,
     ae_state *_state);
void mlpunserializeold(/* Real    */ ae_vector* ra,
     multilayerperceptron* network,
     ae_state *_state);
void mlprandomize(multilayerperceptron* network, ae_state *_state);
void mlprandomizefull(multilayerperceptron* network, ae_state *_state);
void mlpinitpreprocessor(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     ae_state *_state);
void mlpproperties(multilayerperceptron* network,
     ae_int_t* nin,
     ae_int_t* nout,
     ae_int_t* wcount,
     ae_state *_state);
ae_int_t mlpgetinputscount(multilayerperceptron* network,
     ae_state *_state);
ae_int_t mlpgetoutputscount(multilayerperceptron* network,
     ae_state *_state);
ae_int_t mlpgetweightscount(multilayerperceptron* network,
     ae_state *_state);
ae_bool mlpissoftmax(multilayerperceptron* network, ae_state *_state);
ae_int_t mlpgetlayerscount(multilayerperceptron* network,
     ae_state *_state);
ae_int_t mlpgetlayersize(multilayerperceptron* network,
     ae_int_t k,
     ae_state *_state);
void mlpgetinputscaling(multilayerperceptron* network,
     ae_int_t i,
     double* mean,
     double* sigma,
     ae_state *_state);
void mlpgetoutputscaling(multilayerperceptron* network,
     ae_int_t i,
     double* mean,
     double* sigma,
     ae_state *_state);
void mlpgetneuroninfo(multilayerperceptron* network,
     ae_int_t k,
     ae_int_t i,
     ae_int_t* fkind,
     double* threshold,
     ae_state *_state);
double mlpgetweight(multilayerperceptron* network,
     ae_int_t k0,
     ae_int_t i0,
     ae_int_t k1,
     ae_int_t i1,
     ae_state *_state);
void mlpsetinputscaling(multilayerperceptron* network,
     ae_int_t i,
     double mean,
     double sigma,
     ae_state *_state);
void mlpsetoutputscaling(multilayerperceptron* network,
     ae_int_t i,
     double mean,
     double sigma,
     ae_state *_state);
void mlpsetneuroninfo(multilayerperceptron* network,
     ae_int_t k,
     ae_int_t i,
     ae_int_t fkind,
     double threshold,
     ae_state *_state);
void mlpsetweight(multilayerperceptron* network,
     ae_int_t k0,
     ae_int_t i0,
     ae_int_t k1,
     ae_int_t i1,
     double w,
     ae_state *_state);
void mlpactivationfunction(double net,
     ae_int_t k,
     double* f,
     double* df,
     double* d2f,
     ae_state *_state);
void mlpprocess(multilayerperceptron* network,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
void mlpprocessi(multilayerperceptron* network,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
double mlperror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     ae_state *_state);
double mlperrorn(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     ae_state *_state);
ae_int_t mlpclserror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     ae_state *_state);
double mlprelclserror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double mlpavgce(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double mlprmserror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double mlpavgerror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double mlpavgrelerror(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
void mlpgrad(multilayerperceptron* network,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* desiredy,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_state *_state);
void mlpgradn(multilayerperceptron* network,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* desiredy,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_state *_state);
void mlpgradbatch(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_state *_state);
void mlpgradnbatch(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     double* e,
     /* Real    */ ae_vector* grad,
     ae_state *_state);
void mlphessiannbatch(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     double* e,
     /* Real    */ ae_vector* grad,
     /* Real    */ ae_matrix* h,
     ae_state *_state);
void mlphessianbatch(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     double* e,
     /* Real    */ ae_vector* grad,
     /* Real    */ ae_matrix* h,
     ae_state *_state);
void mlpinternalprocessvector(/* Integer */ ae_vector* structinfo,
     /* Real    */ ae_vector* weights,
     /* Real    */ ae_vector* columnmeans,
     /* Real    */ ae_vector* columnsigmas,
     /* Real    */ ae_vector* neurons,
     /* Real    */ ae_vector* dfdnet,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
void mlpalloc(ae_serializer* s,
     multilayerperceptron* network,
     ae_state *_state);
void mlpserialize(ae_serializer* s,
     multilayerperceptron* network,
     ae_state *_state);
void mlpunserialize(ae_serializer* s,
     multilayerperceptron* network,
     ae_state *_state);
ae_bool _multilayerperceptron_init(multilayerperceptron* p, ae_state *_state, ae_bool make_automatic);
ae_bool _multilayerperceptron_init_copy(multilayerperceptron* dst, multilayerperceptron* src, ae_state *_state, ae_bool make_automatic);
void _multilayerperceptron_clear(multilayerperceptron* p);
void mnltrainh(/* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t nclasses,
     ae_int_t* info,
     logitmodel* lm,
     mnlreport* rep,
     ae_state *_state);
void mnlprocess(logitmodel* lm,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
void mnlprocessi(logitmodel* lm,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
void mnlunpack(logitmodel* lm,
     /* Real    */ ae_matrix* a,
     ae_int_t* nvars,
     ae_int_t* nclasses,
     ae_state *_state);
void mnlpack(/* Real    */ ae_matrix* a,
     ae_int_t nvars,
     ae_int_t nclasses,
     logitmodel* lm,
     ae_state *_state);
void mnlcopy(logitmodel* lm1, logitmodel* lm2, ae_state *_state);
double mnlavgce(logitmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double mnlrelclserror(logitmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double mnlrmserror(logitmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double mnlavgerror(logitmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double mnlavgrelerror(logitmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t ssize,
     ae_state *_state);
ae_int_t mnlclserror(logitmodel* lm,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
ae_bool _logitmodel_init(logitmodel* p, ae_state *_state, ae_bool make_automatic);
ae_bool _logitmodel_init_copy(logitmodel* dst, logitmodel* src, ae_state *_state, ae_bool make_automatic);
void _logitmodel_clear(logitmodel* p);
ae_bool _logitmcstate_init(logitmcstate* p, ae_state *_state, ae_bool make_automatic);
ae_bool _logitmcstate_init_copy(logitmcstate* dst, logitmcstate* src, ae_state *_state, ae_bool make_automatic);
void _logitmcstate_clear(logitmcstate* p);
ae_bool _mnlreport_init(mnlreport* p, ae_state *_state, ae_bool make_automatic);
ae_bool _mnlreport_init_copy(mnlreport* dst, mnlreport* src, ae_state *_state, ae_bool make_automatic);
void _mnlreport_clear(mnlreport* p);
void mcpdcreate(ae_int_t n, mcpdstate* s, ae_state *_state);
void mcpdcreateentry(ae_int_t n,
     ae_int_t entrystate,
     mcpdstate* s,
     ae_state *_state);
void mcpdcreateexit(ae_int_t n,
     ae_int_t exitstate,
     mcpdstate* s,
     ae_state *_state);
void mcpdcreateentryexit(ae_int_t n,
     ae_int_t entrystate,
     ae_int_t exitstate,
     mcpdstate* s,
     ae_state *_state);
void mcpdaddtrack(mcpdstate* s,
     /* Real    */ ae_matrix* xy,
     ae_int_t k,
     ae_state *_state);
void mcpdsetec(mcpdstate* s,
     /* Real    */ ae_matrix* ec,
     ae_state *_state);
void mcpdaddec(mcpdstate* s,
     ae_int_t i,
     ae_int_t j,
     double c,
     ae_state *_state);
void mcpdsetbc(mcpdstate* s,
     /* Real    */ ae_matrix* bndl,
     /* Real    */ ae_matrix* bndu,
     ae_state *_state);
void mcpdaddbc(mcpdstate* s,
     ae_int_t i,
     ae_int_t j,
     double bndl,
     double bndu,
     ae_state *_state);
void mcpdsetlc(mcpdstate* s,
     /* Real    */ ae_matrix* c,
     /* Integer */ ae_vector* ct,
     ae_int_t k,
     ae_state *_state);
void mcpdsettikhonovregularizer(mcpdstate* s, double v, ae_state *_state);
void mcpdsetprior(mcpdstate* s,
     /* Real    */ ae_matrix* pp,
     ae_state *_state);
void mcpdsetpredictionweights(mcpdstate* s,
     /* Real    */ ae_vector* pw,
     ae_state *_state);
void mcpdsolve(mcpdstate* s, ae_state *_state);
void mcpdresults(mcpdstate* s,
     /* Real    */ ae_matrix* p,
     mcpdreport* rep,
     ae_state *_state);
ae_bool _mcpdstate_init(mcpdstate* p, ae_state *_state, ae_bool make_automatic);
ae_bool _mcpdstate_init_copy(mcpdstate* dst, mcpdstate* src, ae_state *_state, ae_bool make_automatic);
void _mcpdstate_clear(mcpdstate* p);
ae_bool _mcpdreport_init(mcpdreport* p, ae_state *_state, ae_bool make_automatic);
ae_bool _mcpdreport_init_copy(mcpdreport* dst, mcpdreport* src, ae_state *_state, ae_bool make_automatic);
void _mcpdreport_clear(mcpdreport* p);
void mlptrainlm(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     ae_int_t* info,
     mlpreport* rep,
     ae_state *_state);
void mlptrainlbfgs(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     double wstep,
     ae_int_t maxits,
     ae_int_t* info,
     mlpreport* rep,
     ae_state *_state);
void mlptraines(multilayerperceptron* network,
     /* Real    */ ae_matrix* trnxy,
     ae_int_t trnsize,
     /* Real    */ ae_matrix* valxy,
     ae_int_t valsize,
     double decay,
     ae_int_t restarts,
     ae_int_t* info,
     mlpreport* rep,
     ae_state *_state);
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
     ae_state *_state);
void mlpkfoldcvlm(multilayerperceptron* network,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     ae_int_t foldscount,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* cvrep,
     ae_state *_state);
ae_bool _mlpreport_init(mlpreport* p, ae_state *_state, ae_bool make_automatic);
ae_bool _mlpreport_init_copy(mlpreport* dst, mlpreport* src, ae_state *_state, ae_bool make_automatic);
void _mlpreport_clear(mlpreport* p);
ae_bool _mlpcvreport_init(mlpcvreport* p, ae_state *_state, ae_bool make_automatic);
ae_bool _mlpcvreport_init_copy(mlpcvreport* dst, mlpcvreport* src, ae_state *_state, ae_bool make_automatic);
void _mlpcvreport_clear(mlpcvreport* p);
void mlpecreate0(ae_int_t nin,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreate1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreate2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreateb0(ae_int_t nin,
     ae_int_t nout,
     double b,
     double d,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreateb1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     double b,
     double d,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreateb2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     double b,
     double d,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreater0(ae_int_t nin,
     ae_int_t nout,
     double a,
     double b,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreater1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     double a,
     double b,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreater2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     double a,
     double b,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreatec0(ae_int_t nin,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreatec1(ae_int_t nin,
     ae_int_t nhid,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreatec2(ae_int_t nin,
     ae_int_t nhid1,
     ae_int_t nhid2,
     ae_int_t nout,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecreatefromnetwork(multilayerperceptron* network,
     ae_int_t ensemblesize,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpecopy(mlpensemble* ensemble1,
     mlpensemble* ensemble2,
     ae_state *_state);
void mlperandomize(mlpensemble* ensemble, ae_state *_state);
void mlpeproperties(mlpensemble* ensemble,
     ae_int_t* nin,
     ae_int_t* nout,
     ae_state *_state);
ae_bool mlpeissoftmax(mlpensemble* ensemble, ae_state *_state);
void mlpeprocess(mlpensemble* ensemble,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
void mlpeprocessi(mlpensemble* ensemble,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
double mlperelclserror(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double mlpeavgce(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double mlpermserror(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double mlpeavgerror(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
double mlpeavgrelerror(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     ae_state *_state);
void mlpebagginglm(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     ae_int_t* info,
     mlpreport* rep,
     mlpcvreport* ooberrors,
     ae_state *_state);
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
     ae_state *_state);
void mlpetraines(mlpensemble* ensemble,
     /* Real    */ ae_matrix* xy,
     ae_int_t npoints,
     double decay,
     ae_int_t restarts,
     ae_int_t* info,
     mlpreport* rep,
     ae_state *_state);
void mlpealloc(ae_serializer* s, mlpensemble* ensemble, ae_state *_state);
void mlpeserialize(ae_serializer* s,
     mlpensemble* ensemble,
     ae_state *_state);
void mlpeunserialize(ae_serializer* s,
     mlpensemble* ensemble,
     ae_state *_state);
ae_bool _mlpensemble_init(mlpensemble* p, ae_state *_state, ae_bool make_automatic);
ae_bool _mlpensemble_init_copy(mlpensemble* dst, mlpensemble* src, ae_state *_state, ae_bool make_automatic);
void _mlpensemble_clear(mlpensemble* p);
void pcabuildbasis(/* Real    */ ae_matrix* x,
     ae_int_t npoints,
     ae_int_t nvars,
     ae_int_t* info,
     /* Real    */ ae_vector* s2,
     /* Real    */ ae_matrix* v,
     ae_state *_state);

}
#endif

