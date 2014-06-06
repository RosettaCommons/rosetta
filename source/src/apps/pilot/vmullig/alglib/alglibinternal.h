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
#ifndef _alglibinternal_pkg_h
#define _alglibinternal_pkg_h
#include "ap.h"


/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS COMPUTATIONAL CORE DECLARATIONS (DATATYPES)
//
/////////////////////////////////////////////////////////////////////////
namespace alglib_impl
{
typedef struct
{
    ae_vector ia0;
    ae_vector ia1;
    ae_vector ia2;
    ae_vector ia3;
    ae_vector ra0;
    ae_vector ra1;
    ae_vector ra2;
    ae_vector ra3;
} apbuffers;
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
} linminstate;
typedef struct
{
    ae_bool needf;
    ae_vector x;
    double f;
    ae_int_t n;
    ae_vector xbase;
    ae_vector s;
    double stplen;
    double fcur;
    double stpmax;
    ae_int_t fmax;
    ae_int_t nfev;
    ae_int_t info;
    rcommstate rstate;
} armijostate;
typedef struct
{
    ae_vector plan;
    ae_vector precomputed;
    ae_vector tmpbuf;
    ae_vector stackbuf;
} ftplan;

}

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS C++ INTERFACE
//
/////////////////////////////////////////////////////////////////////////
namespace alglib
{


}

/////////////////////////////////////////////////////////////////////////
//
// THIS SECTION CONTAINS COMPUTATIONAL CORE DECLARATIONS (FUNCTIONS)
//
/////////////////////////////////////////////////////////////////////////
namespace alglib_impl
{
ae_int_t getrdfserializationcode(ae_state *_state);
ae_int_t getkdtreeserializationcode(ae_state *_state);
ae_int_t getmlpserializationcode(ae_state *_state);
ae_int_t getmlpeserializationcode(ae_state *_state);
ae_int_t getrbfserializationcode(ae_state *_state);
ae_bool approxequalrel(double a, double b, double tol, ae_state *_state);
void taskgenint1d(double a,
     double b,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
void taskgenint1dequidist(double a,
     double b,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
void taskgenint1dcheb1(double a,
     double b,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
void taskgenint1dcheb2(double a,
     double b,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     ae_state *_state);
ae_bool aredistinct(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_state *_state);
ae_bool aresameboolean(ae_bool v1, ae_bool v2, ae_state *_state);
void bvectorsetlengthatleast(/* Boolean */ ae_vector* x,
     ae_int_t n,
     ae_state *_state);
void ivectorsetlengthatleast(/* Integer */ ae_vector* x,
     ae_int_t n,
     ae_state *_state);
void rvectorsetlengthatleast(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_state *_state);
void rmatrixsetlengthatleast(/* Real    */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state);
void rmatrixresize(/* Real    */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state);
ae_bool isfinitevector(/* Real    */ ae_vector* x,
     ae_int_t n,
     ae_state *_state);
ae_bool isfinitecvector(/* Complex */ ae_vector* z,
     ae_int_t n,
     ae_state *_state);
ae_bool apservisfinitematrix(/* Real    */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state);
ae_bool apservisfinitecmatrix(/* Complex */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state);
ae_bool isfinitertrmatrix(/* Real    */ ae_matrix* x,
     ae_int_t n,
     ae_bool isupper,
     ae_state *_state);
ae_bool apservisfinitectrmatrix(/* Complex */ ae_matrix* x,
     ae_int_t n,
     ae_bool isupper,
     ae_state *_state);
ae_bool apservisfiniteornanmatrix(/* Real    */ ae_matrix* x,
     ae_int_t m,
     ae_int_t n,
     ae_state *_state);
double safepythag2(double x, double y, ae_state *_state);
double safepythag3(double x, double y, double z, ae_state *_state);
ae_int_t saferdiv(double x, double y, double* r, ae_state *_state);
double safeminposrv(double x, double y, double v, ae_state *_state);
void apperiodicmap(double* x,
     double a,
     double b,
     double* k,
     ae_state *_state);
double boundval(double x, double b1, double b2, ae_state *_state);
void alloccomplex(ae_serializer* s, ae_complex v, ae_state *_state);
void serializecomplex(ae_serializer* s, ae_complex v, ae_state *_state);
ae_complex unserializecomplex(ae_serializer* s, ae_state *_state);
void allocrealarray(ae_serializer* s,
     /* Real    */ ae_vector* v,
     ae_int_t n,
     ae_state *_state);
void serializerealarray(ae_serializer* s,
     /* Real    */ ae_vector* v,
     ae_int_t n,
     ae_state *_state);
void unserializerealarray(ae_serializer* s,
     /* Real    */ ae_vector* v,
     ae_state *_state);
void allocintegerarray(ae_serializer* s,
     /* Integer */ ae_vector* v,
     ae_int_t n,
     ae_state *_state);
void serializeintegerarray(ae_serializer* s,
     /* Integer */ ae_vector* v,
     ae_int_t n,
     ae_state *_state);
void unserializeintegerarray(ae_serializer* s,
     /* Integer */ ae_vector* v,
     ae_state *_state);
void allocrealmatrix(ae_serializer* s,
     /* Real    */ ae_matrix* v,
     ae_int_t n0,
     ae_int_t n1,
     ae_state *_state);
void serializerealmatrix(ae_serializer* s,
     /* Real    */ ae_matrix* v,
     ae_int_t n0,
     ae_int_t n1,
     ae_state *_state);
void unserializerealmatrix(ae_serializer* s,
     /* Real    */ ae_matrix* v,
     ae_state *_state);
void copyintegerarray(/* Integer */ ae_vector* src,
     /* Integer */ ae_vector* dst,
     ae_state *_state);
void copyrealarray(/* Real    */ ae_vector* src,
     /* Real    */ ae_vector* dst,
     ae_state *_state);
void copyrealmatrix(/* Real    */ ae_matrix* src,
     /* Real    */ ae_matrix* dst,
     ae_state *_state);
ae_int_t recsearch(/* Integer */ ae_vector* a,
     ae_int_t nrec,
     ae_int_t nheader,
     ae_int_t i0,
     ae_int_t i1,
     /* Integer */ ae_vector* b,
     ae_state *_state);
ae_bool _apbuffers_init(apbuffers* p, ae_state *_state, ae_bool make_automatic);
ae_bool _apbuffers_init_copy(apbuffers* dst, apbuffers* src, ae_state *_state, ae_bool make_automatic);
void _apbuffers_clear(apbuffers* p);
void tagsort(/* Real    */ ae_vector* a,
     ae_int_t n,
     /* Integer */ ae_vector* p1,
     /* Integer */ ae_vector* p2,
     ae_state *_state);
void tagsortbuf(/* Real    */ ae_vector* a,
     ae_int_t n,
     /* Integer */ ae_vector* p1,
     /* Integer */ ae_vector* p2,
     apbuffers* buf,
     ae_state *_state);
void tagsortfasti(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     /* Real    */ ae_vector* bufa,
     /* Integer */ ae_vector* bufb,
     ae_int_t n,
     ae_state *_state);
void tagsortfastr(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* b,
     /* Real    */ ae_vector* bufa,
     /* Real    */ ae_vector* bufb,
     ae_int_t n,
     ae_state *_state);
void tagsortfast(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* bufa,
     ae_int_t n,
     ae_state *_state);
void tagsortmiddleir(/* Integer */ ae_vector* a,
     /* Real    */ ae_vector* b,
     ae_int_t offset,
     ae_int_t n,
     ae_state *_state);
void tagheappushi(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     ae_int_t* n,
     double va,
     ae_int_t vb,
     ae_state *_state);
void tagheapreplacetopi(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     ae_int_t n,
     double va,
     ae_int_t vb,
     ae_state *_state);
void tagheappopi(/* Real    */ ae_vector* a,
     /* Integer */ ae_vector* b,
     ae_int_t* n,
     ae_state *_state);
ae_int_t lowerbound(/* Real    */ ae_vector* a,
     ae_int_t n,
     double t,
     ae_state *_state);
ae_int_t upperbound(/* Real    */ ae_vector* a,
     ae_int_t n,
     double t,
     ae_state *_state);
void rankx(/* Real    */ ae_vector* x,
     ae_int_t n,
     apbuffers* buf,
     ae_state *_state);
ae_bool cmatrixrank1f(ae_int_t m,
     ae_int_t n,
     /* Complex */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     /* Complex */ ae_vector* u,
     ae_int_t iu,
     /* Complex */ ae_vector* v,
     ae_int_t iv,
     ae_state *_state);
ae_bool rmatrixrank1f(ae_int_t m,
     ae_int_t n,
     /* Real    */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     /* Real    */ ae_vector* u,
     ae_int_t iu,
     /* Real    */ ae_vector* v,
     ae_int_t iv,
     ae_state *_state);
ae_bool cmatrixmvf(ae_int_t m,
     ae_int_t n,
     /* Complex */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     ae_int_t opa,
     /* Complex */ ae_vector* x,
     ae_int_t ix,
     /* Complex */ ae_vector* y,
     ae_int_t iy,
     ae_state *_state);
ae_bool rmatrixmvf(ae_int_t m,
     ae_int_t n,
     /* Real    */ ae_matrix* a,
     ae_int_t ia,
     ae_int_t ja,
     ae_int_t opa,
     /* Real    */ ae_vector* x,
     ae_int_t ix,
     /* Real    */ ae_vector* y,
     ae_int_t iy,
     ae_state *_state);
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
     ae_state *_state);
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
     ae_state *_state);
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
     ae_state *_state);
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
     ae_state *_state);
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
     ae_state *_state);
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
     ae_state *_state);
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
     ae_state *_state);
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
     ae_state *_state);
void hermitianmatrixvectormultiply(/* Complex */ ae_matrix* a,
     ae_bool isupper,
     ae_int_t i1,
     ae_int_t i2,
     /* Complex */ ae_vector* x,
     ae_complex alpha,
     /* Complex */ ae_vector* y,
     ae_state *_state);
void hermitianrank2update(/* Complex */ ae_matrix* a,
     ae_bool isupper,
     ae_int_t i1,
     ae_int_t i2,
     /* Complex */ ae_vector* x,
     /* Complex */ ae_vector* y,
     /* Complex */ ae_vector* t,
     ae_complex alpha,
     ae_state *_state);
void generatereflection(/* Real    */ ae_vector* x,
     ae_int_t n,
     double* tau,
     ae_state *_state);
void applyreflectionfromtheleft(/* Real    */ ae_matrix* c,
     double tau,
     /* Real    */ ae_vector* v,
     ae_int_t m1,
     ae_int_t m2,
     ae_int_t n1,
     ae_int_t n2,
     /* Real    */ ae_vector* work,
     ae_state *_state);
void applyreflectionfromtheright(/* Real    */ ae_matrix* c,
     double tau,
     /* Real    */ ae_vector* v,
     ae_int_t m1,
     ae_int_t m2,
     ae_int_t n1,
     ae_int_t n2,
     /* Real    */ ae_vector* work,
     ae_state *_state);
void complexgeneratereflection(/* Complex */ ae_vector* x,
     ae_int_t n,
     ae_complex* tau,
     ae_state *_state);
void complexapplyreflectionfromtheleft(/* Complex */ ae_matrix* c,
     ae_complex tau,
     /* Complex */ ae_vector* v,
     ae_int_t m1,
     ae_int_t m2,
     ae_int_t n1,
     ae_int_t n2,
     /* Complex */ ae_vector* work,
     ae_state *_state);
void complexapplyreflectionfromtheright(/* Complex */ ae_matrix* c,
     ae_complex tau,
     /* Complex */ ae_vector* v,
     ae_int_t m1,
     ae_int_t m2,
     ae_int_t n1,
     ae_int_t n2,
     /* Complex */ ae_vector* work,
     ae_state *_state);
void symmetricmatrixvectormultiply(/* Real    */ ae_matrix* a,
     ae_bool isupper,
     ae_int_t i1,
     ae_int_t i2,
     /* Real    */ ae_vector* x,
     double alpha,
     /* Real    */ ae_vector* y,
     ae_state *_state);
void symmetricrank2update(/* Real    */ ae_matrix* a,
     ae_bool isupper,
     ae_int_t i1,
     ae_int_t i2,
     /* Real    */ ae_vector* x,
     /* Real    */ ae_vector* y,
     /* Real    */ ae_vector* t,
     double alpha,
     ae_state *_state);
double vectornorm2(/* Real    */ ae_vector* x,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state);
ae_int_t vectoridxabsmax(/* Real    */ ae_vector* x,
     ae_int_t i1,
     ae_int_t i2,
     ae_state *_state);
ae_int_t columnidxabsmax(/* Real    */ ae_matrix* x,
     ae_int_t i1,
     ae_int_t i2,
     ae_int_t j,
     ae_state *_state);
ae_int_t rowidxabsmax(/* Real    */ ae_matrix* x,
     ae_int_t j1,
     ae_int_t j2,
     ae_int_t i,
     ae_state *_state);
double upperhessenberg1norm(/* Real    */ ae_matrix* a,
     ae_int_t i1,
     ae_int_t i2,
     ae_int_t j1,
     ae_int_t j2,
     /* Real    */ ae_vector* work,
     ae_state *_state);
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
     ae_state *_state);
void inplacetranspose(/* Real    */ ae_matrix* a,
     ae_int_t i1,
     ae_int_t i2,
     ae_int_t j1,
     ae_int_t j2,
     /* Real    */ ae_vector* work,
     ae_state *_state);
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
     ae_state *_state);
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
     ae_state *_state);
double pythag2(double x, double y, ae_state *_state);
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
     ae_state *_state);
void applyrotationsfromtheleft(ae_bool isforward,
     ae_int_t m1,
     ae_int_t m2,
     ae_int_t n1,
     ae_int_t n2,
     /* Real    */ ae_vector* c,
     /* Real    */ ae_vector* s,
     /* Real    */ ae_matrix* a,
     /* Real    */ ae_vector* work,
     ae_state *_state);
void applyrotationsfromtheright(ae_bool isforward,
     ae_int_t m1,
     ae_int_t m2,
     ae_int_t n1,
     ae_int_t n2,
     /* Real    */ ae_vector* c,
     /* Real    */ ae_vector* s,
     /* Real    */ ae_matrix* a,
     /* Real    */ ae_vector* work,
     ae_state *_state);
void generaterotation(double f,
     double g,
     double* cs,
     double* sn,
     double* r,
     ae_state *_state);
ae_bool upperhessenbergschurdecomposition(/* Real    */ ae_matrix* h,
     ae_int_t n,
     /* Real    */ ae_matrix* s,
     ae_state *_state);
void internalschurdecomposition(/* Real    */ ae_matrix* h,
     ae_int_t n,
     ae_int_t tneeded,
     ae_int_t zneeded,
     /* Real    */ ae_vector* wr,
     /* Real    */ ae_vector* wi,
     /* Real    */ ae_matrix* z,
     ae_int_t* info,
     ae_state *_state);
void rmatrixtrsafesolve(/* Real    */ ae_matrix* a,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     double* s,
     ae_bool isupper,
     ae_bool istrans,
     ae_bool isunit,
     ae_state *_state);
void safesolvetriangular(/* Real    */ ae_matrix* a,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     double* s,
     ae_bool isupper,
     ae_bool istrans,
     ae_bool isunit,
     ae_bool normin,
     /* Real    */ ae_vector* cnorm,
     ae_state *_state);
ae_bool rmatrixscaledtrsafesolve(/* Real    */ ae_matrix* a,
     double sa,
     ae_int_t n,
     /* Real    */ ae_vector* x,
     ae_bool isupper,
     ae_int_t trans,
     ae_bool isunit,
     double maxgrowth,
     ae_state *_state);
ae_bool cmatrixscaledtrsafesolve(/* Complex */ ae_matrix* a,
     double sa,
     ae_int_t n,
     /* Complex */ ae_vector* x,
     ae_bool isupper,
     ae_int_t trans,
     ae_bool isunit,
     double maxgrowth,
     ae_state *_state);
void xdot(/* Real    */ ae_vector* a,
     /* Real    */ ae_vector* b,
     ae_int_t n,
     /* Real    */ ae_vector* temp,
     double* r,
     double* rerr,
     ae_state *_state);
void xcdot(/* Complex */ ae_vector* a,
     /* Complex */ ae_vector* b,
     ae_int_t n,
     /* Real    */ ae_vector* temp,
     ae_complex* r,
     double* rerr,
     ae_state *_state);
void linminnormalized(/* Real    */ ae_vector* d,
     double* stp,
     ae_int_t n,
     ae_state *_state);
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
     ae_state *_state);
void armijocreate(ae_int_t n,
     /* Real    */ ae_vector* x,
     double f,
     /* Real    */ ae_vector* s,
     double stp,
     double stpmax,
     ae_int_t fmax,
     armijostate* state,
     ae_state *_state);
ae_bool armijoiteration(armijostate* state, ae_state *_state);
void armijoresults(armijostate* state,
     ae_int_t* info,
     double* stp,
     double* f,
     ae_state *_state);
ae_bool _linminstate_init(linminstate* p, ae_state *_state, ae_bool make_automatic);
ae_bool _linminstate_init_copy(linminstate* dst, linminstate* src, ae_state *_state, ae_bool make_automatic);
void _linminstate_clear(linminstate* p);
ae_bool _armijostate_init(armijostate* p, ae_state *_state, ae_bool make_automatic);
ae_bool _armijostate_init_copy(armijostate* dst, armijostate* src, ae_state *_state, ae_bool make_automatic);
void _armijostate_clear(armijostate* p);
void ftbasegeneratecomplexfftplan(ae_int_t n,
     ftplan* plan,
     ae_state *_state);
void ftbasegeneraterealfftplan(ae_int_t n, ftplan* plan, ae_state *_state);
void ftbasegeneraterealfhtplan(ae_int_t n, ftplan* plan, ae_state *_state);
void ftbaseexecuteplan(/* Real    */ ae_vector* a,
     ae_int_t aoffset,
     ae_int_t n,
     ftplan* plan,
     ae_state *_state);
void ftbaseexecuteplanrec(/* Real    */ ae_vector* a,
     ae_int_t aoffset,
     ftplan* plan,
     ae_int_t entryoffset,
     ae_int_t stackptr,
     ae_state *_state);
void ftbasefactorize(ae_int_t n,
     ae_int_t tasktype,
     ae_int_t* n1,
     ae_int_t* n2,
     ae_state *_state);
ae_bool ftbaseissmooth(ae_int_t n, ae_state *_state);
ae_int_t ftbasefindsmooth(ae_int_t n, ae_state *_state);
ae_int_t ftbasefindsmootheven(ae_int_t n, ae_state *_state);
double ftbasegetflopestimate(ae_int_t n, ae_state *_state);
ae_bool _ftplan_init(ftplan* p, ae_state *_state, ae_bool make_automatic);
ae_bool _ftplan_init_copy(ftplan* dst, ftplan* src, ae_state *_state, ae_bool make_automatic);
void _ftplan_clear(ftplan* p);
double nulog1p(double x, ae_state *_state);
double nuexpm1(double x, ae_state *_state);
double nucosm1(double x, ae_state *_state);

}
#endif

