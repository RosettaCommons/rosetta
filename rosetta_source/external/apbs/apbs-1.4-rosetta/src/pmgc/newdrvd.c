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

#include "newdrvd.h"

VEXTERNC void Vnewdriv(
        int *iparm, double *rparm,
        int *iwork, double *rwork,
        double *u,
        double *xf, double *yf, double *zf,
        double *gxcf, double *gycf, double *gzcf,
        double *a1cf, double *a2cf, double *a3cf,
        double *ccf, double *fcf, double *tcf) {

    int nxc;    /// @todo: Doc
    int nyc;    /// @todo: Doc
    int nzc;    /// @todo: Doc
    int nf;     /// @todo: Doc
    int nc;     /// @todo: Doc
    int narr;   /// @todo: Doc
    int narrc;  /// @todo: Doc
    int n_rpc;  /// @todo: Doc
    int n_iz;   /// @todo: Doc
    int n_ipc; /// @todo: Doc
    int iretot; /// @todo: Doc
    int iintot; /// @todo: Doc

    int nrwk;   /// @todo: Doc
    int niwk;   /// @todo: Doc
    int nx;     /// @todo: Doc
    int ny;     /// @todo: Doc
    int nz;     /// @todo: Doc
    int nlev;   /// @todo: Doc
    int ierror; /// @todo: Doc
    int maxlev; /// @todo: Doc
    int mxlv;   /// @todo: Doc
    int mgcoar; /// @todo: Doc
    int mgdisc; /// @todo: Doc
    int mgsolv; /// @todo: Doc
    int k_iz;   /// @todo: Doc
    int k_w1;   /// @todo: Doc
    int k_w2;   /// @todo: Doc
    int k_ipc;  /// @todo: Doc
    int k_rpc;  /// @todo: Doc
    int k_ac;   /// @todo: Doc
    int k_cc;   /// @todo: Doc
    int k_fc;   /// @todo: Doc
    int k_pc;   /// @todo: Doc

    // Decode some parameters
    nrwk   = VAT(iparm, 1);
    niwk   = VAT(iparm, 2);
    nx     = VAT(iparm, 3);
    ny     = VAT(iparm, 4);
    nz     = VAT(iparm, 5);
    nlev   = VAT(iparm, 6);

    // Some checks on input ***
    VASSERT_MSG0(nlev > 0, "The nlev parameter must be positive");
    VASSERT_MSG0(nx > 0, "The nx parameter must be positive");
    VASSERT_MSG0(ny > 0, "The ny parameter must be positive");
    VASSERT_MSG0(nz > 0, "The nz parameter must be positive");

    mxlv = Vmaxlev(nx,ny,nz);

    VASSERT_MSG1(nlev <= mxlv, "Max lev for your grid size is: %d", mxlv);

    // Basic grid sizes, etc.
    mgcoar = VAT(iparm, 18);
    mgdisc = VAT(iparm, 19);
    mgsolv = VAT(iparm, 21);

    Vmgsz(&mgcoar, &mgdisc, &mgsolv,
            &nx, &ny, &nz,
            &nlev,
            &nxc, &nyc, &nzc,
            &nf, &nc,
            &narr, &narrc,
            &n_rpc, &n_iz, &n_ipc,
            &iretot, &iintot);

    // Allocate space for two additional work vectors ***
    iretot = iretot + 2 * nf;

    // Some more checks on input
    VASSERT_MSG1( nrwk >= iretot, "Real work space must be: %d", iretot );
    VASSERT_MSG1( niwk >= iintot, "Integer work space must be: %d", iintot );

    // Split up the integer work array
    k_iz   = 1;
    k_ipc  = k_iz   + n_iz;

    // Split up the real work array
    k_rpc  = 1;
    k_cc   = k_rpc  + n_rpc;
    k_fc   = k_cc   + narr;
    k_w1   = k_fc   + narr;
    k_w2   = k_w1   + nf;
    k_pc   = k_w2   + nf;
    k_ac   = k_pc   + 27 * narrc;
    // k_ac_after = 4*nf + 14*narrc;

    // Call the Newton Driver
    Vnewdriv2(iparm, rparm,
            &nx, &ny, &nz,
            u, RAT(iwork, k_iz),
            RAT(rwork, k_w1),  RAT(rwork, k_w2),
            RAT(iwork, k_ipc), RAT(rwork, k_rpc),
            RAT(rwork, k_pc),  RAT(rwork, k_ac), RAT(rwork, k_cc), RAT(rwork, k_fc),
            xf, yf, zf,
            gxcf, gycf, gzcf,
            a1cf, a2cf, a3cf,
            ccf, fcf, tcf);
}



VPUBLIC void Vnewdriv2(int *iparm, double *rparm,
        int *nx, int *ny, int *nz,
        double *u, int *iz,
        double *w1, double *w2,
        int *ipc, double *rpc,
        double *pc, double *ac, double *cc, double *fc,
        double *xf, double *yf, double *zf,
        double *gxcf, double *gycf, double *gzcf,
        double *a1cf, double *a2cf, double *a3cf,
        double *ccf, double *fcf, double *tcf) {

    int mgkey;      /// @todo:  Doc
    int nlev;       /// @todo:  Doc
    int itmax;      /// @todo:  Doc
    int iok;        /// @todo:  Doc
    int iinfo;      /// @todo:  Doc
    int istop;      /// @todo:  Doc
    int ipkey;      /// @todo:  Doc
    int nu1;        /// @todo:  Doc
    int nu2;        /// @todo:  Doc
    int ilev;       /// @todo:  Doc
    int ido;        /// @todo:  Doc
    int iters;      /// @todo:  Doc
    int ierror;     /// @todo:  Doc
    int nlev_real;  /// @todo:  Doc
    int ibound;     /// @todo:  Doc
    int mgprol;     /// @todo:  Doc
    int mgcoar;     /// @todo:  Doc
    int mgsolv;     /// @todo:  Doc
    int mgdisc;     /// @todo:  Doc
    int mgsmoo;     /// @todo:  Doc
    int mode;       /// @todo:  Doc
    double epsiln;  /// @todo:  Doc
    double epsmac;  /// @todo:  Doc
    double errtol;  /// @todo:  Doc
    double omegal;  /// @todo:  Doc
    double omegan;  /// @todo:  Doc
    double bf;      /// @todo:  Doc
    double oh;      /// @todo:  Doc
    double tsetupf; /// @todo:  Doc
    double tsetupc; /// @todo:  Doc
    double tsolve;  /// @todo:  Doc



    // Utility variables
    int numlev;

    int iok_t;
    int iters_t;
    double rsnrm_t;
    double rsden_t;
    double orsnrm_t;

    int i;



    // Decode the iparm array
    nlev   = VAT(iparm, 6);
    nu1    = VAT(iparm, 7);
    nu2    = VAT(iparm, 8);
    mgkey  = VAT(iparm, 9);
    itmax  = VAT(iparm, 10);
    istop  = VAT(iparm, 11);
    iinfo  = VAT(iparm, 12);
    ipkey  = VAT(iparm, 14);
    mgprol = VAT(iparm, 17);
    mgcoar = VAT(iparm, 18);
    mgdisc = VAT(iparm, 19);
    mgsmoo = VAT(iparm, 20);
    mgsolv = VAT(iparm, 21);

    errtol = VAT(rparm,  1);
    omegal = VAT(rparm,  9);
    omegan = VAT(rparm, 10);

    Vprtstp(0, -99, 0.0, 0.0, 0.0);

    // Build the multigrid data structure in iz
    Vbuildstr(nx, ny, nz, &nlev, iz);

    // Start the timer
    Vnm_tstart(30, "Vnewdrv2: fine problem setup");

    // Build op and rhs on fine grid ***
    ido = 0;
    Vbuildops(nx, ny, nz,
            &nlev, &ipkey, &iinfo, &ido, iz,
            &mgprol, &mgcoar, &mgsolv, &mgdisc,
            ipc, rpc,
            pc, ac, cc, fc,
            xf, yf, zf,
            gxcf, gycf, gzcf,
            a1cf, a2cf, a3cf,
            ccf, fcf, tcf);

    // Stop the timer
    Vnm_tstop(30, "Vnewdrv2: fine problem setup");

    // Start the timer
    Vnm_tstart(30, "Vnewdrv2: coarse problem setup");

    // Build op and rhs on all coarse grids
    ido = 1;
    Vbuildops(nx, ny, nz,
            &nlev, &ipkey, &iinfo, &ido, iz,
            &mgprol, &mgcoar, &mgsolv, &mgdisc,
            ipc, rpc,
            pc, ac, cc, fc,
            xf, yf, zf,
            gxcf, gycf, gzcf,
            a1cf, a2cf, a3cf,
            ccf, fcf, tcf);

    // Stop the timer
    Vnm_tstop(30, "Vnewdrv2: coarse problem setup");



    /* ******************************************************************* *
     * *** this overwrites the rhs array provided by pde specification *** *
     * ****** compute an algebraically produced rhs for the given tcf  *** *
     * ******************************************************************* */
    mode = 1;

    if (istop == 4 || istop == 5) {
        WARN_UNTESTED;
        Vbuildalg(nx, ny, nz,
                 &mode, &nlev, iz,
                 ipc, rpc, ac, cc, ccf, tcf, fc, fcf);
    }



    // Determine machine epsilon
    epsiln = Vnm_epsmac();

    // Impose zero dirichlet boundary conditions (now in source fcn)
    VfboundPMG00(nx, ny, nz, u);

    // Start the timer
    Vnm_tstart(30, "Vnewdrv2: solve");

    // Call specified multigrid method
    nlev_real = nlev;
    iok  = 1;
    ilev = 1;
    if (mgkey == 0) {
        Vnewton(nx, ny, nz,
                u, iz,
                ccf, fcf, w1, w2,
                &istop, &itmax, &iters, &ierror,
                &nlev, &ilev, &nlev_real, &mgsolv,
                &iok, &iinfo,
                &epsiln, &errtol, &omegan,
                &nu1, &nu2, &mgsmoo,
                a1cf, a2cf, a3cf,
                ipc, rpc,
                pc, ac, cc, fc, tcf);
    } else if (mgkey == 1) {
        Vfnewton(nx, ny, nz,
                u, iz, ccf, fcf, w1, w2,
                &istop, &itmax, &iters, &ierror,
                &nlev, &ilev, &nlev_real, &mgsolv,
                &iok, &iinfo,
                &epsiln, &errtol, &omegan,
                &nu1, &nu2, &mgsmoo,
                a1cf, a2cf, a3cf,
                ipc, rpc,
                pc, ac, cc, fc, tcf);
    } else {
        VABORT_MSG1("Bad mgkey given: %d", mgkey);
    }

    // Stop the timer
    Vnm_tstop(30, "Vnewdrv2: solve");

    // Restore boundary conditions
    ibound = 1;
    VfboundPMG(&ibound, nx, ny, nz, u, gxcf, gycf, gzcf);
}
