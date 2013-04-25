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

#include "mgdrvd.h"

VPUBLIC void Vmgdriv(int* iparm, double* rparm,
        int* iwork, double* rwork, double* u,
        double* xf, double* yf, double* zf,
        double* gxcf, double* gycf, double* gzcf,
        double* a1cf, double* a2cf, double* a3cf,
        double* ccf, double* fcf, double* tcf) {

    // The following variables will be returned from mgsz
    int nxc    = 0;
    int nyc    = 0;
    int nzc    = 0;
    int nf     = 0;
    int nc     = 0;
    int narr   = 0;
    int narrc  = 0;
    int n_rpc  = 0;
    int n_iz   = 0;
    int n_ipc  = 0;
    int iretot = 0;
    int iintot = 0;

    // Miscellaneous variables
    int nrwk   = 0;
    int niwk   = 0;
    int nx     = 0;
    int ny     = 0;
    int nz     = 0;
    int nlev   = 0;
    int ierror = 0;
    int mxlv   = 0;
    int mgcoar = 0;
    int mgdisc = 0;
    int mgsolv = 0;
    int k_iz   = 0;
    int k_ipc  = 0;
    int k_rpc  = 0;
    int k_ac   = 0;
    int k_cc   = 0;
    int k_fc   = 0;
    int k_pc   = 0;

    // Utility pointers to help in passing values
    int *iz     = VNULL;
    int *ipc    = VNULL;
    double *rpc = VNULL;
    double *pc  = VNULL;
    double *ac  = VNULL;
    double *cc  = VNULL;
    double *fc  = VNULL;

    // Decode some parameters
    nrwk   = VAT(iparm, 1);
    niwk   = VAT(iparm, 2);
    nx     = VAT(iparm, 3);
    ny     = VAT(iparm, 4);
    nz     = VAT(iparm, 5);
    nlev   = VAT(iparm, 6);

    // Perform some checks on input
    VASSERT_MSG1(nlev > 0, "nlev must be positive: %d", nlev);
    VASSERT_MSG1(  nx > 0, "nx must be positive: %d", nx);
    VASSERT_MSG1(  ny > 0, "nv must be positive: %d", ny);
    VASSERT_MSG1(  nz > 0, "nz must be positive: %d", nz);

    mxlv = Vmaxlev(nx, ny, nz);
    VASSERT_MSG2(
        nlev <= mxlv,
        "number of levels exceeds maximum: %d > %d",
        nlev, mxlv
        );


    // Extract basic grid sizes, etc.
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

    // Perform some more checks on input
    VASSERT_MSG2(
        iretot >= nrwk,
        "real workspace exceeds maximum size: %d > %d",
        nrwk, iretot
        );
    VASSERT_MSG2(
        iintot >= niwk,
        "integer workspace exceeds maximum size: %d > %d",
        niwk, iintot
        );

    // Split up the integer work array
    k_iz  = 1;
    k_ipc = k_iz + n_iz;

    // Split up the real work array ***
    k_rpc = 1;
    k_cc  = k_rpc + n_rpc;
    k_fc  = k_cc  + narr;
    k_pc  = k_fc  + narr;
    k_ac  = k_pc  + 27 * narrc;
    // k_ac_after =  4 * nf +  4 * narrc;
    // k_ac_after =  4 * nf + 14 * narrc;
    // k_ac_after = 14 * nf + 14 * narrc;

    iz  = RAT(iwork, k_iz);
    ipc = RAT(iwork, k_ipc);

    rpc = RAT(rwork, k_rpc);
    pc  = RAT(rwork, k_pc);
    ac  = RAT(rwork, k_ac);
    cc  = RAT(rwork, k_cc);
    fc  = RAT(rwork, k_fc);

    // Call the multigrid driver
    Vmgdriv2(iparm, rparm,
                &nx, &ny, &nz,
                u,
                iz, ipc, rpc,
                pc, ac, cc, fc,
                xf, yf, zf,
                gxcf, gycf, gzcf,
                a1cf, a2cf, a3cf,
                ccf, fcf, tcf);
}

VPUBLIC void Vmgdriv2(int *iparm, double *rparm,
        int *nx, int *ny, int *nz,
        double *u,
        int *iz, int *ipc, double *rpc,
        double *pc, double *ac, double *cc, double *fc,
        double *xf, double *yf, double *zf,
        double *gxcf, double *gycf, double *gzcf,
        double *a1cf, double *a2cf, double *a3cf,
        double *ccf, double *fcf, double *tcf) {

    // @todo Document this function

    // Miscellaneous Variables
    int mgkey     = 0;
    int itmax     = 0;
    int iok       = 0;
    int iinfo     = 0;
    int istop     = 0;
    int ipkey     = 0;
    int nu1       = 0;
    int nu2       = 0;
    int ilev      = 0;
    int ido       = 0;
    int iters     = 0;
    int ierror    = 0;
    int nlev_real = 0;
    int ibound    = 0;
    int mgprol    = 0;
    int mgcoar    = 0;
    int mgsolv    = 0;
    int mgdisc    = 0;
    int mgsmoo    = 0;
    int iperf     = 0;
    int mode      = 0;

    double epsiln  = 0.0;
    double epsmac  = 0.0;
    double errtol  = 0.0;
    double omegal  = 0.0;
    double omegan  = 0.0;
    double bf      = 0.0;
    double oh      = 0.0;
    double tsetupf = 0.0;
    double tsetupc = 0.0;
    double tsolve  = 0.0;



    // More miscellaneous variables
    int itmax_p = 0;
    int iters_p = 0;
    int iok_p   = 0;
    int iinfo_p = 0;

    double errtol_p    = 0.0;
    double rho_p       = 0.0;
    double rho_min     = 0.0;
    double rho_max     = 0.0;
    double rho_min_mod = 0.0;
    double rho_max_mod = 0.0;

    int nxf   = 0;
    int nyf   = 0;
    int nzf   = 0;
    int nxc   = 0;
    int nyc   = 0;
    int nzc   = 0;
    int level = 0;
    int nlevd = 0;



    // Utility variables
    int numlev = 0;

    // Get the value of nlev here because it is needed for the iz matrix
    int nlev   = VAT(iparm,  6);
    MAT2(iz, 50, nlev);

    // Decode integer parameters from the iparm array
    nu1    = VAT(iparm,  7);
    nu2    = VAT(iparm,  8);
    mgkey  = VAT(iparm,  9);
    itmax  = VAT(iparm, 10);
    istop  = VAT(iparm, 11);
    iinfo  = VAT(iparm, 12);
    ipkey  = VAT(iparm, 14);
    mode   = VAT(iparm, 16);
    mgprol = VAT(iparm, 17);
    mgcoar = VAT(iparm, 18);
    mgdisc = VAT(iparm, 19);
    mgsmoo = VAT(iparm, 20);
    mgsolv = VAT(iparm, 21);
    iperf  = VAT(iparm, 22);

    // Decode real parameters from the rparm array
    errtol = VAT(rparm,  1);
    omegal = VAT(rparm,  9);
    omegan = VAT(rparm, 10);

    /// @todo replace timer setup
    Vprtstp(0, -99, 0.0, 0.0, 0.0);

    // Build the multigrid data structure in iz
    Vbuildstr(nx, ny, nz, &nlev, iz);

    // Start the timer
    Vnm_tstart(30, "Vmgdrv2: fine problem setup");

    // Build operator and rhs on fine grid
    ido = 0;
    Vbuildops(nx, ny, nz,
            &nlev, &ipkey, &iinfo, &ido, iz,
            &mgprol, &mgcoar, &mgsolv, &mgdisc,
            ipc, rpc, pc, ac, cc, fc,
            xf, yf, zf,
            gxcf, gycf, gzcf,
            a1cf, a2cf, a3cf,
            ccf, fcf, tcf);

    // Stop the timer
    Vnm_tstop(30, "Vmgdrv2: fine problem setup");

    // Start the timer
    Vnm_tstart(30, "Vmgdrv2: coarse problem setup");

    // Build operator and rhs on all coarse grids
    ido = 1;
    Vbuildops(nx, ny, nz,
            &nlev, &ipkey, &iinfo, &ido, iz,
            &mgprol, &mgcoar, &mgsolv, &mgdisc,
            ipc, rpc, pc, ac, cc, fc,
            xf, yf, zf,
            gxcf, gycf, gzcf,
            a1cf, a2cf, a3cf,
            ccf, fcf, tcf);

    // Stop the timer
    Vnm_tstop(30, "Vmgdrv2: coarse problem setup");

    // Determine Machine Epsilon
    epsiln = Vnm_epsmac();

    /******************************************************************
        *** analysis ***
        *** note: we destroy the rhs function "fc" here in "mpower" ***
        ******************************************************************/

    // errtol and itmax
    itmax_p   = 1000;
    iok_p     = 0;
    nlev_real = nlev;
    nlevd     = nlev_real;

    // Finest level initialization
    nxf  = *nx;
    nyf  = *ny;
    nzf  = *nz;

    // Go down grids: compute max/min eigenvalues of all operators
    for (level=1; level <= nlev_real; level++) {
        nlevd = nlev_real - level + 1;

        // Move down the grids
        if (level != 1) {

            // Find new grid size
            numlev = 1;
            Vmkcors(&numlev, &nxf, &nyf, &nzf, &nxc, &nyc, &nzc);

            // New grid size ***
            nxf = nxc;
            nyf = nyc;
            nzf = nzc;
        }

        if (iinfo > 1) {
            VMESSAGE3("Analysis ==> (%3d, %3d, %3d)", nxf, nyf, nzf);
        }


        // Largest eigenvalue of the system matrix A
        if (iperf == 1 || iperf == 3) {

            if (iinfo > 1) {
                VMESSAGE0("Power calculating rho(A)");
            }

            iters_p   = 0;
            iinfo_p   = iinfo;
            errtol_p  = 1.0e-4;

            Vpower(&nxf, &nyf, &nzf,
                    iz, &level,
                    ipc, rpc, ac, cc,
                    a1cf, a2cf, a3cf, ccf,
                    &rho_max, &rho_max_mod, &errtol_p,
                    &itmax_p, &iters_p, &iinfo_p);

            if (iinfo > 1) {
                VMESSAGE1("Power iters   = %d", iters_p);
                VMESSAGE1("Power eigmax  = %f", rho_max);
                VMESSAGE1("Power (MODEL) = %f", rho_max_mod);
            }

            // Smallest eigenvalue of the system matrix A
            if (iinfo > 1) {
                VMESSAGE0("Ipower calculating lambda_min(A)...");
            }

            iters_p   = 0;
            iinfo_p   = iinfo;
            errtol_p  = 1.0e-4;

            Vazeros(&nxf, &nyf, &nzf, u);

            Vipower(&nxf, &nyf, &nzf, u, iz,
                    a1cf, a2cf, a3cf, ccf, fcf,
                    &rho_min, &rho_min_mod, &errtol_p, &itmax_p, &iters_p,
                    &nlevd, &level, &nlev_real, &mgsolv,
                    &iok_p, &iinfo_p, &epsiln, &errtol, &omegal,
                    &nu1, &nu2, &mgsmoo,
                    ipc, rpc, pc, ac, cc, tcf);

            if (iinfo > 1) {
                VMESSAGE1("Ipower iters   = %d", iters_p);
                VMESSAGE1("Ipower eigmin  = %f", rho_min);
                VMESSAGE1("Ipower (MODEL) = %f", rho_min_mod);

                // Condition number estimate
                VMESSAGE1("Condition number  = %f", rho_max / rho_min);
                VMESSAGE1("Condition (MODEL) = %f", rho_max_mod / rho_min_mod);
            }
        }

        // Spectral radius of the multigrid operator M
        // NOTE: due to lack of vectors, we destroy "fc" in mpower...
        if (iperf == 2 || iperf == 3) {

            if (iinfo > 1) {
                VMESSAGE0("Mpower calculating rho(M)");
            }

            iters_p   = 0;
            iinfo_p   = iinfo;
            errtol_p  = epsiln;

            Vazeros(&nxf, &nyf, &nzf, RAT(u, VAT2(iz, 1, level)));

            WARN_UNTESTED;
            Vmpower(&nxf, &nyf, &nzf, u, iz,
                            a1cf, a2cf, a3cf, ccf, fcf,
                            &rho_p, &errtol_p, &itmax_p, &iters_p,
                            &nlevd, &level, &nlev_real, &mgsolv,
                            &iok_p, &iinfo_p, &epsiln,
                            &errtol, &omegal, &nu1, &nu2, &mgsmoo,
                            ipc, rpc, pc, ac, cc, fc, tcf);

            if (iinfo > 1) {
                VMESSAGE1("Mpower iters  = %d", iters_p);
                VMESSAGE1("Mpower rho(M) = %f", rho_p);
            }
        }

        // Reinitialize the solution function

        Vazeros(&nxf, &nyf, &nzf, RAT(u, VAT2(iz, 1, level)));

        // Next grid
    }

    // Reinitialize the solution function
        Vazeros(nx, ny, nz, u);

    /*******************************************************************
     *** this overwrites the rhs array provided by pde specification ***
     ***** compute an algebraically produced rhs for the given tcf *****/

    if (istop == 4 || istop == 5 || iperf != 0 ) {

            if (iinfo > 1) {
                VMESSAGE0("Generating algebraic RHS from your soln...");
            }


            WARN_UNTESTED;
            Vbuildalg(nx, ny, nz, &mode, &nlev, iz,
                    ipc, rpc, ac, cc, ccf, tcf, fc, fcf);
    }

    /*******************************************************************/

    // Impose zero dirichlet boundary conditions (now in source fcn)
        VfboundPMG00(nx, ny, nz, u);

    // Start the timer
    Vnm_tstart(30, "Vmgdrv2: solve");

    // Call specified multigrid method
    if (mode == 0 || mode == 2) {
            nlev_real = nlev;
            iok  = 1;
            ilev = 1;

            if (mgkey == 0) {

                Vmvcs(nx, ny, nz,
                        u, iz, a1cf, a2cf, a3cf, ccf,
                        &istop, &itmax, &iters, &ierror, &nlev,
                        &ilev, &nlev_real, &mgsolv,
                        &iok, &iinfo, &epsiln, &errtol, &omegal,
                        &nu1, &nu2, &mgsmoo,
                        ipc, rpc, pc, ac, cc, fc, tcf);

            } else if (mgkey == 1) {

                Vmvcs(nx, ny, nz,
                        u, iz, a1cf, a2cf, a3cf, ccf,
                        &istop, &itmax, &iters, &ierror, &nlev,
                        &ilev, &nlev_real, &mgsolv,
                        &iok, &iinfo, &epsiln, &errtol, &omegal,
                        &nu1, &nu2, &mgsmoo,
                        ipc, rpc, pc, ac, cc, fc, tcf);

            } else {
                VABORT_MSG1("Bad mgkey given: %d", mgkey);
            }
    }

    if (mode == 1 || mode == 2) {

            nlev_real = nlev;
            iok  = 1;
            ilev = 1;

            if (mgkey == 0) {

                Vmvfas(nx, ny, nz,
                        u, iz, a1cf, a2cf, a3cf, ccf, fcf,
                        &istop, &itmax, &iters, &ierror, &nlev,
                        &ilev, &nlev_real, &mgsolv,
                        &iok, &iinfo, &epsiln, &errtol, &omegan,
                        &nu1, &nu2, &mgsmoo,
                        ipc, rpc, pc, ac, cc, fc, tcf);

            } else if (mgkey == 1) {

                Vfmvfas(nx, ny, nz,
                        u, iz,
                        a1cf, a2cf, a3cf, ccf, fcf,
                        &istop, &itmax, &iters, &ierror, &nlev,
                        &ilev, &nlev_real, &mgsolv,
                        &iok, &iinfo, &epsiln, &errtol, &omegan,
                        &nu1, &nu2, &mgsmoo,
                        ipc, rpc, pc, ac, cc, fc, tcf);

        } else {
            VABORT_MSG1("Bad mgkey given: %d", mgkey);
        }
    }

    // Stop the timer
    Vnm_tstop(30, "Vmgdrv2: solve");

    // Restore boundary conditions
    ibound = 1;

    VfboundPMG(&ibound, nx, ny, nz, u, gxcf, gycf, gzcf);
}



VPUBLIC void Vmgsz(int *mgcoar, int *mgdisc, int *mgsolv,
        int *nx, int *ny, int *nz,
        int *nlev,
        int *nxc, int *nyc, int *nzc,
        int *nf, int *nc,
        int *narr, int *narrc,
        int *n_rpc, int *n_iz, int *n_ipc,
        int *iretot, int *iintot) {

    // Constants: num of different types of arrays in mg code
    int num_nf = 0;
    int num_narr = 2;
    int num_narrc = 27;

    // Misc variables
    int nc_band, num_band, n_band;
    int nxf, nyf, nzf;
    int level;
    int num_nf_oper, num_narrc_oper;

    // Utility variables
    int numlev;

    // Go down grids: compute max/min eigenvalues of all operators
    *nf   = *nx * *ny * *nz;

    *narr = *nf;

    nxf  = *nx;
    nyf  = *ny;
    nzf  = *nz;

    *nxc  = *nx;
    *nyc  = *ny;
    *nzc  = *nz;

    for (level=2; level<=*nlev; level++) {

        //find new grid size ***

        numlev = 1;
        Vmkcors(&numlev, &nxf, &nyf, &nzf, nxc, nyc, nzc);

        // New grid size
        nxf = *nxc;
        nyf = *nyc;
        nzf = *nzc;

        // Add the unknowns on this level to the total
        *narr += nxf * nyf * nzf;
    }
    *nc = *nxc * *nyc * *nzc;
    *narrc = *narr - *nf;

    // Box or fem on fine grid?
    if (*mgdisc == 0) {
        num_nf_oper = 4;
    } else if (*mgdisc == 1) {
        num_nf_oper = 14;
    } else {
        Vnm_print(2, "Vmgsz: invalid mgdisc parameter: %d\n", *mgdisc);
    }

    // Galerkin or standard coarsening?
    if ((*mgcoar == 0 || *mgcoar == 1) && *mgdisc == 0) {
        num_narrc_oper = 4;
    } else if (*mgcoar == 2) {
        num_narrc_oper = 14;
    } else {
        Vnm_print(2, "Vmgsz: invalid mgcoar parameter: %d\n", *mgcoar);
    }

    // Symmetric banded linpack storage on coarse grid
    if (*mgsolv == 0) {
        n_band = 0;
    } else if (*mgsolv == 1) {
        if ((*mgcoar == 0 || *mgcoar == 1) && *mgdisc == 0) {
            num_band = 1 + (*nxc - 2) * (*nyc - 2);
        } else {
            num_band = 1 + (*nxc - 2) * (*nyc - 2) + (*nxc - 2) + 1;
        }
        nc_band = (*nxc - 2) * (*nyc - 2) * (*nzc - 2);
        n_band  = nc_band * num_band;
    } else {
        Vnm_print(2, "Vmgsz: invalid mgsolv parameter: %d\n", *mgsolv);
    }

    // Info work array required storage
    *n_rpc = 100 * (*nlev + 1);

    // Resulting total required real storage for method
    *iretot = num_narr * *narr
            + (num_nf    + num_nf_oper) * *nf
            + (num_narrc + num_narrc_oper) * *narrc
            + n_band
            + *n_rpc;

    // The integer storage parameters ***
    *n_iz  = 50  * (*nlev + 1);
    *n_ipc = 100 * (*nlev + 1);

    // Resulting total required integer storage for method
    *iintot = *n_iz + *n_ipc;
}
