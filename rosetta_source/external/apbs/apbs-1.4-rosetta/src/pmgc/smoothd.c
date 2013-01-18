/**
 *  @ingroup PMGC
 *  @author  Tucker Beck [fortran ->c translation], Michael Holst [original]
 *  @brief
 *  @version $Id:
 *
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 * Nathan A. Baker (nathan.baker@pnl.gov)
 * Pacific Northwest National Laboratory
 *
 * Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2010-2012 Battelle Memorial Institute. Developed at the Pacific Northwest National Laboratory, operated by Battelle Memorial Institute, Pacific Northwest Division for the U.S. Department Energy.  Portions Copyright (c) 2002-2010, Washington University in St. Louis.  Portions Copyright (c) 2002-2010, Nathan A. Baker.  Portions Copyright (c) 1999-2002, The Regents of the University of California. Portions Copyright (c) 1995, Michael Holst.
 * All rights reserved.
 *
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * -  Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * - Neither the name of Washington University in St. Louis nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @endverbatim
 */

#include "smoothd.h"

VEXTERNC void Vsmooth(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double *ac, double *cc, double *fc,
        double *x, double *w1, double *w2, double *r,
        int *itmax, int *iters,
        double *errtol, double *omega,
        int *iresid, int *iadjoint, int *meth) {

    // Do in one step
    if (*meth == 0) {
        VABORT_MSG0( "wjac not yet translated" );
        //wjac(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,itmax,iters,errtol,omega,iresid,iadjoint);
    } else if (*meth == 1) {
        Vgsrb(nx, ny, nz,
                ipc, rpc,
                ac, cc, fc,
                x, w1, w2, r,
                itmax, iters,
                errtol, omega,
                iresid, iadjoint);
    } else if (*meth == 2) {
        VABORT_MSG0( "sor not yet translated" );
        //sor(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,itmax,iters,errtol,omega,iresid,iadjoint);
    } else if (*meth == 3) {
        VABORT_MSG0( "rich not yet translated" );
        //rich(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,itmax,iters,errtol,omega,iresid,iadjoint);
    } else if (*meth == 4) {
        Vcghs(nx, ny, nz,
                ipc, rpc,
                ac, cc, fc,
                x, w1, w2, r,
                itmax, iters,
                errtol, omega,
                iresid, iadjoint);
    } else {
        VABORT_MSG1("Bad smoothing routine specified = %d", *meth);
    }
}


VEXTERNC void Vnsmooth(int *nx, int *ny, int *nz,
        int *ipc, double *rpc,
        double *ac, double *cc, double *fc,
        double *x, double *w1, double *w2, double *r,
        int *itmax, int *iters,
        double *errtol, double *omega,
        int *iresid, int *iadjoint, int *meth) {

    WARN_UNTESTED;

    // Do in one step
    if (*meth == 0) {
        VABORT_MSG0( "nwjac not yet translated" );
        //nwjac(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,itmax,iters,errtol,omega,iresid,iadjoint)
    } else if (*meth == 1) {
        VABORT_MSG0( "ngsrb not yet translated" );
        //ngsrb(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,itmax,iters,errtol,omega,iresid,iadjoint)
    } else if (*meth == 2) {
        VABORT_MSG0( "nsor not yet translated" );
        //nsor(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,itmax,iters,errtol,omega,iresid,iadjoint)
    } else if (*meth == 3) {
        VABORT_MSG0( "nrich not yet translated" );
        //nrich(nx,ny,nz,ipc,rpc,ac,cc,fc,x,w1,w2,r,itmax,iters,errtol,omega,iresid,iadjoint)
    } else {
        VABORT_MSG1("Bad smoothing routine specified: %d", *meth );
    }
}
