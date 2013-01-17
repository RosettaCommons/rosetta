/**
*  @file    routines.c
 *  @author  Nathan Baker
 *  @brief   Supporting routines for APBS front end
 *  @version $Id$
 *  @attention
 *  @verbatim
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 *  Nathan A. Baker (nathan.baker@pnnl.gov)
 *  Pacific Northwest National Laboratory
 *
 *  Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2010-2012 Battelle Memorial Institute. Developed at the
 * Pacific Northwest National Laboratory, operated by Battelle Memorial
 * Institute, Pacific Northwest Division for the U.S. Department of Energy.
 *
 * Portions Copyright (c) 2002-2010, Washington University in St. Louis.
 * Portions Copyright (c) 2002-2010, Nathan A. Baker.
 * Portions Copyright (c) 1999-2002, The Regents of the University of
 * California.
 * Portions Copyright (c) 1995, Michael Holst.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * Neither the name of the developer nor the names of its contributors may be
 * used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 * @endverbatim
 */

#include "routines.h"

VEMBED(rcsid="$Id$")

VPUBLIC void startVio() { Vio_start(); }

VPUBLIC Vparam* loadParameter(NOsh *nosh) {

    Vparam *param = VNULL;

    if (nosh->gotparm) {
        param = Vparam_ctor();
        switch (nosh->parmfmt) {
            case NPF_FLAT:
                Vnm_tprint( 1, "Reading parameter data from %s.\n",
                            nosh->parmpath);
                if (Vparam_readFlatFile(param, "FILE", "ASC", VNULL,
                                        nosh->parmpath) != 1) {
                    Vnm_tprint(2, "Error reading parameter file (%s)!\n", nosh->parmpath);
                    return VNULL;
                }
                    break;
            case NPF_XML:
                Vnm_tprint( 1, "Reading parameter data from %s.\n",
                            nosh->parmpath);
                if (Vparam_readXMLFile(param, "FILE", "ASC", VNULL,
                                       nosh->parmpath) != 1) {
                    Vnm_tprint(2, "Error reading parameter file (%s)!\n", nosh->parmpath);
                    return VNULL;
                }
                    break;
            default:
                Vnm_tprint(2, "Error! Undefined parameter file type (%d)!\n", nosh->parmfmt);
                return VNULL;
        } /* switch parmfmt */
    }

    return param;
}


VPUBLIC int loadMolecules(NOsh *nosh, Vparam *param, Valist *alist[NOSH_MAXMOL]) {

    int i;
    int use_params = 0;
    Vrc_Codes rc;

    Vio *sock = VNULL;

    Vnm_tprint( 1, "Got paths for %d molecules\n", nosh->nmol);
    if (nosh->nmol <= 0) {
        Vnm_tprint(2, "You didn't specify any molecules (correctly)!\n");
        Vnm_tprint(2, "Bailing out!\n");
        return 0;
    }

    if (nosh->gotparm) {
        if (param == VNULL) {
            Vnm_tprint(2, "Error!  You don't have a valid parameter object!\n");
            return 0;
        }
        use_params = 1;
    }

    for (i=0; i<nosh->nmol; i++) {
        if(alist[i] == VNULL){
            alist[i] = Valist_ctor();
        }else{
            alist[i] = VNULL;
            alist[i] = Valist_ctor();
        }

        switch (nosh->molfmt[i]) {
            case NMF_PQR:
                    /* Print out a warning to the user letting them know that we are overriding PQR
                    values for charge, radius and epsilon */
                    if (use_params) {
                        Vnm_print(2, "\nWARNING!!  Radius/charge information from PQR file %s\n", nosh->molpath[i]);
                        Vnm_print(2, "will be replaced with data from parameter file (%s)!\n", nosh->parmpath);
                    }
                    Vnm_tprint( 1, "Reading PQR-format atom data from %s.\n",
                            nosh->molpath[i]);
                    sock = Vio_ctor("FILE", "ASC", VNULL, nosh->molpath[i], "r");
                    if (sock == VNULL) {
                        Vnm_print(2, "Problem opening virtual socket %s!\n",
                              nosh->molpath[i]);
                        return 0;
                    }
                    if (Vio_accept(sock, 0) < 0) {
                        Vnm_print(2, "Problem accepting virtual socket %s!\n",
                                  nosh->molpath[i]);
                        return 0;
                    }
                    if(use_params){
                        rc = Valist_readPQR(alist[i], param, sock);
                    }else{
                        rc = Valist_readPQR(alist[i], VNULL, sock);
                    }
                    if(rc == 0) return 0;

                    Vio_acceptFree(sock);
                    Vio_dtor(&sock);
                    break;
            case NMF_PDB:
                    /* Load parameters */
                    if (!nosh->gotparm) {
                        Vnm_tprint(2, "NOsh:  Error!  Can't read PDB without specifying PARM file!\n");
                        return 0;
                    }
                    Vnm_tprint( 1, "Reading PDB-format atom data from %s.\n",
                            nosh->molpath[i]);
                    sock = Vio_ctor("FILE", "ASC", VNULL, nosh->molpath[i], "r");
                    if (sock == VNULL) {
                        Vnm_print(2, "Problem opening virtual socket %s!\n",
                              nosh->molpath[i]);
                        return 0;
                    }
                    if (Vio_accept(sock, 0) < 0) {
                        Vnm_print(2, "Problem accepting virtual socket %s!\n", nosh->molpath[i]);
                        return 0;
                    }
                    rc = Valist_readPDB(alist[i], param, sock);
                    /* If we are looking for an atom/residue that does not exist
                     * then abort and return 0 */
                    if(rc == 0)
                        return 0;

                    Vio_acceptFree(sock);
                    Vio_dtor(&sock);
                    break;
            case NMF_XML:
                Vnm_tprint( 1, "Reading XML-format atom data from %s.\n",
                            nosh->molpath[i]);
                sock = Vio_ctor("FILE", "ASC", VNULL, nosh->molpath[i], "r");
                if (sock == VNULL) {
                    Vnm_print(2, "Problem opening virtual socket %s!\n",
                              nosh->molpath[i]);
                    return 0;
                }
                    if (Vio_accept(sock, 0) < 0) {
                        Vnm_print(2, "Problem accepting virtual socket %s!\n",
                                  nosh->molpath[i]);
                        return 0;
                    }
                    if(use_params){
                        rc = Valist_readXML(alist[i], param, sock);
                    }else{
                        rc = Valist_readXML(alist[i], VNULL, sock);
                    }
                    if(rc == 0)
                        return 0;

                Vio_acceptFree(sock);
                Vio_dtor(&sock);
                break;
            default:
                Vnm_tprint(2, "NOsh:  Error!  Undefined molecule file type \
(%d)!\n", nosh->molfmt[i]);
                return 0;
        } /* switch molfmt */

        if (rc != 1) {
            Vnm_tprint( 2, "Error while reading molecule from %s\n",
                        nosh->molpath[i]);
            return 0;
        }

        Vnm_tprint( 1, "  %d atoms\n", Valist_getNumberAtoms(alist[i]));
        Vnm_tprint( 1, "  Centered at (%4.3e, %4.3e, %4.3e)\n",
                    alist[i]->center[0], alist[i]->center[1],
                    alist[i]->center[2]);
        Vnm_tprint( 1, "  Net charge %3.2e e\n", alist[i]->charge);

    }

    return 1;

}

VPUBLIC void killMolecules(NOsh *nosh, Valist *alist[NOSH_MAXMOL]) {

    int i;

#ifndef VAPBSQUIET
    Vnm_tprint( 1, "Destroying %d molecules\n", nosh->nmol);
#endif

    for (i=0; i<nosh->nmol; i++)
        Valist_dtor(&(alist[i]));

}

/**
 * Loads dielectric map path data into NOsh object
 * @return 1 on success, 0 on error
 */
VPUBLIC int loadDielMaps(NOsh *nosh,
                         Vgrid *dielXMap[NOSH_MAXMOL],
                         Vgrid *dielYMap[NOSH_MAXMOL],
                         Vgrid *dielZMap[NOSH_MAXMOL]
                        ) {

    int i,
        ii,
        nx,
        ny,
        nz;
    double sum,
           hx,
           hy,
           hzed,
           xmin,
           ymin,
           zmin;

    // Check to be sure we have dieletric map paths; if not, return.
    if (nosh->ndiel > 0)
        Vnm_tprint( 1, "Got paths for %d dielectric map sets\n",
                    nosh->ndiel);
    else
        return 1;

    // For each dielectric map path, read the data and calculate needed values.
    for (i=0; i<nosh->ndiel; i++) {
        Vnm_tprint( 1, "Reading x-shifted dielectric map data from \
%s:\n", nosh->dielXpath[i]);
        dielXMap[i] = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);

        // Determine the format and read data if the format is valid.
        switch (nosh->dielfmt[i]) {
            // OpenDX (Data Explorer) format
            case VDF_DX:
                if (Vgrid_readDX(dielXMap[i], "FILE", "ASC", VNULL,
                                 nosh->dielXpath[i]) != 1) {
                    Vnm_tprint( 2, "Fatal error while reading from %s\n",
                                nosh->dielXpath[i]);
                    return 0;
                }

                // Set grid sizes
                nx = dielXMap[i]->nx;
                ny = dielXMap[i]->ny;
                nz = dielXMap[i]->nz;

                // Set spacings
                hx = dielXMap[i]->hx;
                hy = dielXMap[i]->hy;
                hzed = dielXMap[i]->hzed;

                // Set minimum lower corner
                xmin = dielXMap[i]->xmin;
                ymin = dielXMap[i]->ymin;
                zmin = dielXMap[i]->zmin;
                Vnm_tprint(1, "  %d x %d x %d grid\n", nx, ny, nz);
                Vnm_tprint(1, "  (%g, %g, %g) A spacings\n", hx, hy, hzed);
                Vnm_tprint(1, "  (%g, %g, %g) A lower corner\n",
                           xmin, ymin, zmin);
                sum = 0;
                for (ii=0; ii<(nx*ny*nz); ii++)
                    sum += (dielXMap[i]->data[ii]);
                    sum = sum*hx*hy*hzed;
                Vnm_tprint(1, "  Volume integral = %3.2e A^3\n", sum);
                break;
            // Binary file (GZip)
            case VDF_GZ:
                if (Vgrid_readGZ(dielXMap[i], nosh->dielXpath[i]) != 1) {
                    Vnm_tprint( 2, "Fatal error while reading from %s\n",
                               nosh->dielXpath[i]);
                    return 0;
                }

                // Set grid sizes
                nx = dielXMap[i]->nx;
                ny = dielXMap[i]->ny;
                nz = dielXMap[i]->nz;

                // Set spacings
                hx = dielXMap[i]->hx;
                hy = dielXMap[i]->hy;
                hzed = dielXMap[i]->hzed;

                // Set minimum lower corner
                xmin = dielXMap[i]->xmin;
                ymin = dielXMap[i]->ymin;
                zmin = dielXMap[i]->zmin;
                Vnm_tprint(1, "  %d x %d x %d grid\n", nx, ny, nz);
                Vnm_tprint(1, "  (%g, %g, %g) A spacings\n", hx, hy, hzed);
                Vnm_tprint(1, "  (%g, %g, %g) A lower corner\n",
                           xmin, ymin, zmin);
                sum = 0;
                for (ii=0; ii<(nx*ny*nz); ii++)
                    sum += (dielXMap[i]->data[ii]);
                sum = sum*hx*hy*hzed;
                Vnm_tprint(1, "  Volume integral = %3.2e A^3\n", sum);
                break;
            // UHBD format
            case VDF_UHBD:
                Vnm_tprint( 2, "UHBD input not supported yet!\n");
                return 0;
            // AVS UCD format
            case VDF_AVS:
                Vnm_tprint( 2, "AVS input not supported yet!\n");
                return 0;
            // FEtk MC Simplex Format (MCSF)
            case VDF_MCSF:
                Vnm_tprint( 2, "MCSF input not supported yet!\n");
                return 0;
            default:
                Vnm_tprint( 2, "Invalid data format (%d)!\n",
                            nosh->dielfmt[i]);
                return 0;
        }
        Vnm_tprint( 1, "Reading y-shifted dielectric map data from \
%s:\n", nosh->dielYpath[i]);
        dielYMap[i] = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);

        // Determine the format and read data if the format is valid.
        switch (nosh->dielfmt[i]) {
            // OpenDX (Data Explorer) format
            case VDF_DX:
                if (Vgrid_readDX(dielYMap[i], "FILE", "ASC", VNULL,
                                 nosh->dielYpath[i]) != 1) {
                    Vnm_tprint( 2, "Fatal error while reading from %s\n",
                                nosh->dielYpath[i]);
                    return 0;
                }

                // Read grid
                nx = dielYMap[i]->nx;
                ny = dielYMap[i]->ny;
                nz = dielYMap[i]->nz;

                // Read spacings
                hx = dielYMap[i]->hx;
                hy = dielYMap[i]->hy;
                hzed = dielYMap[i]->hzed;

                // Read minimum lower corner
                xmin = dielYMap[i]->xmin;
                ymin = dielYMap[i]->ymin;
                zmin = dielYMap[i]->zmin;
                Vnm_tprint(1, "  %d x %d x %d grid\n", nx, ny, nz);
                Vnm_tprint(1, "  (%g, %g, %g) A spacings\n", hx, hy, hzed);
                Vnm_tprint(1, "  (%g, %g, %g) A lower corner\n",
                           xmin, ymin, zmin);
                sum = 0;
                for (ii=0; ii<(nx*ny*nz); ii++)
                    sum += (dielYMap[i]->data[ii]);
                    sum = sum*hx*hy*hzed;
                Vnm_tprint(1, "  Volume integral = %3.2e A^3\n", sum);
                break;
            // Binary file (GZip) format
            case VDF_GZ:
                if (Vgrid_readGZ(dielYMap[i], nosh->dielYpath[i]) != 1) {
                    Vnm_tprint( 2, "Fatal error while reading from %s\n",
                               nosh->dielYpath[i]);
                    return 0;
                }

                // Read grid
                nx = dielYMap[i]->nx;
                ny = dielYMap[i]->ny;
                nz = dielYMap[i]->nz;

                // Read spacings
                hx = dielYMap[i]->hx;
                hy = dielYMap[i]->hy;
                hzed = dielYMap[i]->hzed;

                // Read minimum lower corner
                xmin = dielYMap[i]->xmin;
                ymin = dielYMap[i]->ymin;
                zmin = dielYMap[i]->zmin;
                Vnm_tprint(1, "  %d x %d x %d grid\n", nx, ny, nz);
                Vnm_tprint(1, "  (%g, %g, %g) A spacings\n", hx, hy, hzed);
                Vnm_tprint(1, "  (%g, %g, %g) A lower corner\n",
                           xmin, ymin, zmin);
                sum = 0;
                for (ii=0; ii<(nx*ny*nz); ii++)
                    sum += (dielYMap[i]->data[ii]);
                sum = sum*hx*hy*hzed;
                Vnm_tprint(1, "  Volume integral = %3.2e A^3\n", sum);
                break;
            // UHBD format
            case VDF_UHBD:
                Vnm_tprint( 2, "UHBD input not supported yet!\n");
                return 0;
            // AVS UCD format
            case VDF_AVS:
                Vnm_tprint( 2, "AVS input not supported yet!\n");
                return 0;
            // FEtk MC Simplex Format (MCSF)
            case VDF_MCSF:
                Vnm_tprint( 2, "MCSF input not supported yet!\n");
                return 0;
            default:
                Vnm_tprint( 2, "Invalid data format (%d)!\n",
                            nosh->dielfmt[i]);
                return 0;
        }

        Vnm_tprint( 1, "Reading z-shifted dielectric map data from \
%s:\n", nosh->dielZpath[i]);
        dielZMap[i] = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);

        // Determine the format and read data if the format is valid.
        switch (nosh->dielfmt[i]) {
            // OpenDX (Data Explorer) format
            case VDF_DX:
                if (Vgrid_readDX(dielZMap[i], "FILE", "ASC", VNULL,
                                 nosh->dielZpath[i]) != 1) {
                    Vnm_tprint( 2, "Fatal error while reading from %s\n",
                                nosh->dielZpath[i]);
                    return 0;
                }

                // Read grid
                nx = dielZMap[i]->nx;
                ny = dielZMap[i]->ny;
                nz = dielZMap[i]->nz;

                // Read spacings
                hx = dielZMap[i]->hx;
                hy = dielZMap[i]->hy;
                hzed = dielZMap[i]->hzed;

                // Read minimum lower corner
                xmin = dielZMap[i]->xmin;
                ymin = dielZMap[i]->ymin;
                zmin = dielZMap[i]->zmin;
                Vnm_tprint(1, "  %d x %d x %d grid\n",
                           nx, ny, nz);
                Vnm_tprint(1, "  (%g, %g, %g) A spacings\n",
                           hx, hy, hzed);
                Vnm_tprint(1, "  (%g, %g, %g) A lower corner\n",
                           xmin, ymin, zmin);
                sum = 0;
                for (ii=0; ii<(nx*ny*nz); ii++) sum += (dielZMap[i]->data[ii]);
                    sum = sum*hx*hy*hzed;
                Vnm_tprint(1, "  Volume integral = %3.2e A^3\n", sum);
                break;
            // Binary file (GZip) format
            case VDF_GZ:
                if (Vgrid_readGZ(dielZMap[i], nosh->dielZpath[i]) != 1) {
                    Vnm_tprint( 2, "Fatal error while reading from %s\n",
                               nosh->dielZpath[i]);
                    return 0;
                }

                // Read grid
                nx = dielZMap[i]->nx;
                ny = dielZMap[i]->ny;
                nz = dielZMap[i]->nz;

                // Read spacings
                hx = dielZMap[i]->hx;
                hy = dielZMap[i]->hy;
                hzed = dielZMap[i]->hzed;

                // Read minimum lower corner
                xmin = dielZMap[i]->xmin;
                ymin = dielZMap[i]->ymin;
                zmin = dielZMap[i]->zmin;
                Vnm_tprint(1, "  %d x %d x %d grid\n",
                           nx, ny, nz);
                Vnm_tprint(1, "  (%g, %g, %g) A spacings\n",
                           hx, hy, hzed);
                Vnm_tprint(1, "  (%g, %g, %g) A lower corner\n",
                           xmin, ymin, zmin);
                sum = 0;
                for (ii=0; ii<(nx*ny*nz); ii++) sum += (dielZMap[i]->data[ii]);
                sum = sum*hx*hy*hzed;
                Vnm_tprint(1, "  Volume integral = %3.2e A^3\n", sum);
                break;
            // UHBD format
            case VDF_UHBD:
                Vnm_tprint( 2, "UHBD input not supported yet!\n");
                return 0;
            // AVS UCD format
            case VDF_AVS:
                Vnm_tprint( 2, "AVS input not supported yet!\n");
                return 0;
            // FEtk MC Simplex Format (MCSF)
            case VDF_MCSF:
                Vnm_tprint( 2, "MCSF input not supported yet!\n");
                return 0;
            default:
                Vnm_tprint( 2, "Invalid data format (%d)!\n",
                            nosh->dielfmt[i]);
                return 0;
        }
    }

    return 1;
}

VPUBLIC void killDielMaps(NOsh *nosh,
                          Vgrid *dielXMap[NOSH_MAXMOL],
                          Vgrid *dielYMap[NOSH_MAXMOL],
                          Vgrid *dielZMap[NOSH_MAXMOL]) {

    int i;

    if (nosh->ndiel > 0) {
#ifndef VAPBSQUIET
        Vnm_tprint( 1, "Destroying %d dielectric map sets\n",
                    nosh->ndiel);
#endif
        for (i=0; i<nosh->ndiel; i++) {
            Vgrid_dtor(&(dielXMap[i]));
            Vgrid_dtor(&(dielYMap[i]));
            Vgrid_dtor(&(dielZMap[i]));
        }
    }
    else return;

}

/**
 * @return 0 on failure, 1 on success
 */
VPUBLIC int loadKappaMaps(NOsh *nosh,
                          Vgrid *map[NOSH_MAXMOL]
                         ) {

    int i,
        ii,
        len;
    double sum;

    if (nosh->nkappa > 0)
        Vnm_tprint( 1, "Got paths for %d kappa maps\n", nosh->nkappa);
    else return 1;

    for (i=0; i<nosh->nkappa; i++) {
        Vnm_tprint( 1, "Reading kappa map data from %s:\n",
                    nosh->kappapath[i]);
        map[i] = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);

        // Determine the format and read data if the format is valid.
        switch (nosh->kappafmt[i]) {
            // OpenDX (Data Explorer) format
            case VDF_DX:
                if (Vgrid_readDX(map[i], "FILE", "ASC", VNULL,
                                 nosh->kappapath[i]) != 1) {
                    Vnm_tprint( 2, "Fatal error while reading from %s\n",
                                nosh->kappapath[i]);
                    return 0;
                }
                Vnm_tprint(1, "  %d x %d x %d grid\n",
                           map[i]->nx, map[i]->ny, map[i]->nz);
                Vnm_tprint(1, "  (%g, %g, %g) A spacings\n",
                           map[i]->hx, map[i]->hy, map[i]->hzed);
                Vnm_tprint(1, "  (%g, %g, %g) A lower corner\n",
                           map[i]->xmin, map[i]->ymin, map[i]->zmin);
                sum = 0;
                for (ii = 0, len = map[i]->nx * map[i]->ny * map[i]->nz;
                     ii < len;
                     ii++
                    ) {
                    sum += (map[i]->data[ii]);
                }
                    sum = sum*map[i]->hx*map[i]->hy*map[i]->hzed;
                Vnm_tprint(1, "  Volume integral = %3.2e A^3\n", sum);
                break;
            // UHBD format
            case VDF_UHBD:
                Vnm_tprint( 2, "UHBD input not supported yet!\n");
                return 0;
            // FEtk MC Simplex Format (MCSF)
            case VDF_MCSF:
                Vnm_tprint( 2, "MCSF input not supported yet!\n");
                return 0;
            // AVS UCD format
            case VDF_AVS:
                Vnm_tprint( 2, "AVS input not supported yet!\n");
                return 0;
            // Binary file (GZip) format
            case VDF_GZ:
                if (Vgrid_readGZ(map[i], nosh->kappapath[i]) != 1) {
                    Vnm_tprint( 2, "Fatal error while reading from %s\n",
                               nosh->kappapath[i]);
                    return 0;
                }
                Vnm_tprint(1, "  %d x %d x %d grid\n",
                           map[i]->nx, map[i]->ny, map[i]->nz);
                Vnm_tprint(1, "  (%g, %g, %g) A spacings\n",
                           map[i]->hx, map[i]->hy, map[i]->hzed);
                Vnm_tprint(1, "  (%g, %g, %g) A lower corner\n",
                           map[i]->xmin, map[i]->ymin, map[i]->zmin);
                sum = 0;
                for (ii=0, len=map[i]->nx*map[i]->ny*map[i]->nz; ii<len; ii++) {
                    sum += (map[i]->data[ii]);
                }
                sum = sum*map[i]->hx*map[i]->hy*map[i]->hzed;
                Vnm_tprint(1, "  Volume integral = %3.2e A^3\n", sum);
                break;
            default:
                Vnm_tprint( 2, "Invalid data format (%d)!\n",
                            nosh->kappafmt[i]);
                return 0;
        }
    }

    return 1;

}

VPUBLIC void killKappaMaps(NOsh *nosh, Vgrid *map[NOSH_MAXMOL]) {

    int i;

    if (nosh->nkappa > 0) {
#ifndef VAPBSQUIET
        Vnm_tprint( 1, "Destroying %d kappa maps\n", nosh->nkappa);
#endif
        for (i=0; i<nosh->nkappa; i++) Vgrid_dtor(&(map[i]));
    }
    else return;

}

/**
 * @return 0 on failure, 1 on success
 */
VPUBLIC int loadPotMaps(NOsh *nosh,
                        Vgrid *map[NOSH_MAXMOL]
                       ) {

    int i,
        ii,
        len;
    double sum;

    if (nosh->npot > 0)
        Vnm_tprint( 1, "Got paths for %d potential maps\n", nosh->npot);
    else return 1;

    for (i=0; i<nosh->npot; i++) {
        Vnm_tprint( 1, "Reading potential map data from %s:\n",
                   nosh->potpath[i]);
        map[i] = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);
        switch (nosh->potfmt[i]) {
            // OpenDX (Data Explorer) format
            case VDF_DX:
            // Binary file (GZip) format
            case VDF_GZ:
                if (nosh->potfmt[i] == VDF_DX) {
                    if (Vgrid_readDX(map[i], "FILE", "ASC", VNULL,
                                     nosh->potpath[i]) != 1) {
                        Vnm_tprint( 2, "Fatal error while reading from %s\n",
                                   nosh->potpath[i]);
                        return 0;
                    }
                }else {
                    if (Vgrid_readGZ(map[i], nosh->potpath[i]) != 1) {
                        Vnm_tprint( 2, "Fatal error while reading from %s\n",
                                   nosh->potpath[i]);
                        return 0;
                    }
                }
                Vnm_tprint(1, "  %d x %d x %d grid\n",
                           map[i]->nx, map[i]->ny, map[i]->nz);
                Vnm_tprint(1, "  (%g, %g, %g) A spacings\n",
                           map[i]->hx, map[i]->hy, map[i]->hzed);
                Vnm_tprint(1, "  (%g, %g, %g) A lower corner\n",
                           map[i]->xmin, map[i]->ymin, map[i]->zmin);
                sum = 0;
                for (ii=0,len=map[i]->nx*map[i]->ny*map[i]->nz; ii<len; ii++) {
                    sum += (map[i]->data[ii]);
                }
                sum = sum*map[i]->hx*map[i]->hy*map[i]->hzed;
                Vnm_tprint(1, "  Volume integral = %3.2e A^3\n", sum);
                break;
            // UHBD format
            case VDF_UHBD:
                Vnm_tprint( 2, "UHBD input not supported yet!\n");
                return 0;
            // FEtk MC Simplex Format (MCSF)
            case VDF_MCSF:
                Vnm_tprint( 2, "MCSF input not supported yet!\n");
                return 0;
            // AVS UCD format
            case VDF_AVS:
                Vnm_tprint( 2, "AVS input not supported yet!\n");
                return 0;
            default:
                Vnm_tprint( 2, "Invalid data format (%d)!\n",
                           nosh->potfmt[i]);
                return 0;
        }
    }

    return 1;

}

VPUBLIC void killPotMaps(NOsh *nosh,
                         Vgrid *map[NOSH_MAXMOL]
                        ) {

    int i;

    if (nosh->npot > 0) {
#ifndef VAPBSQUIET
        Vnm_tprint( 1, "Destroying %d potential maps\n", nosh->npot);
#endif
        for (i=0; i<nosh->npot; i++) Vgrid_dtor(&(map[i]));
    }
    else return;

}

/**
 * @return 0 on failure, 1 on success
 */
VPUBLIC int loadChargeMaps(NOsh *nosh,
                           Vgrid *map[NOSH_MAXMOL]
                          ) {

    int i,
        ii,
        len;
    double sum;

    if (nosh->ncharge > 0)
        Vnm_tprint( 1, "Got paths for %d charge maps\n", nosh->ncharge);
    else return 1;

    for (i=0; i<nosh->ncharge; i++) {
        Vnm_tprint( 1, "Reading charge map data from %s:\n",
                    nosh->chargepath[i]);
        map[i] = Vgrid_ctor(0, 0, 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, VNULL);

        // Determine data format and read data
        switch (nosh->chargefmt[i]) {
            case VDF_DX:
                if (Vgrid_readDX(map[i], "FILE", "ASC", VNULL,
                                 nosh->chargepath[i]) != 1) {
                    Vnm_tprint( 2, "Fatal error while reading from %s\n",
                                nosh->chargepath[i]);
                    return 0;
                }
                Vnm_tprint(1, "  %d x %d x %d grid\n",
                           map[i]->nx, map[i]->ny, map[i]->nz);
                Vnm_tprint(1, "  (%g, %g, %g) A spacings\n",
                           map[i]->hx, map[i]->hy, map[i]->hzed);
                Vnm_tprint(1, "  (%g, %g, %g) A lower corner\n",
                           map[i]->xmin, map[i]->ymin, map[i]->zmin);
                sum = 0;
                for (ii=0,len=map[i]->nx*map[i]->ny*map[i]->nz; ii<len; ii++) {
                    sum += (map[i]->data[ii]);
                }
                    sum = sum*map[i]->hx*map[i]->hy*map[i]->hzed;
                Vnm_tprint(1, "  Charge map integral = %3.2e e\n", sum);
                break;
            case VDF_UHBD:
                Vnm_tprint( 2, "UHBD input not supported yet!\n");
                return 0;
            case VDF_AVS:
                Vnm_tprint( 2, "AVS input not supported yet!\n");
                return 0;
            case VDF_MCSF:
                Vnm_tprint(2, "MCSF input not supported yet!\n");
                return 0;
            case VDF_GZ:
                if (Vgrid_readGZ(map[i], nosh->chargepath[i]) != 1) {
                    Vnm_tprint( 2, "Fatal error while reading from %s\n",
                               nosh->chargepath[i]);
                    return 0;
                }
                Vnm_tprint(1, "  %d x %d x %d grid\n",
                           map[i]->nx, map[i]->ny, map[i]->nz);
                Vnm_tprint(1, "  (%g, %g, %g) A spacings\n",
                           map[i]->hx, map[i]->hy, map[i]->hzed);
                Vnm_tprint(1, "  (%g, %g, %g) A lower corner\n",
                           map[i]->xmin, map[i]->ymin, map[i]->zmin);
                sum = 0;
                for (ii=0,len=map[i]->nx*map[i]->ny*map[i]->nz; ii<len; ii++) {
                    sum += (map[i]->data[ii]);
                }
                sum = sum*map[i]->hx*map[i]->hy*map[i]->hzed;
                Vnm_tprint(1, "  Charge map integral = %3.2e e\n", sum);
                break;
            default:
                Vnm_tprint( 2, "Invalid data format (%d)!\n",
                            nosh->kappafmt[i]);
                return 0;
        }
    }

    return 1;

}

VPUBLIC void killChargeMaps(NOsh *nosh,
                            Vgrid *map[NOSH_MAXMOL]
                           ) {

    int i;

    if (nosh->ncharge > 0) {
#ifndef VAPBSQUIET
        Vnm_tprint( 1, "Destroying %d charge maps\n", nosh->ncharge);
#endif

        for (i=0; i<nosh->ncharge; i++) Vgrid_dtor(&(map[i]));
    }

    else return;

}

VPUBLIC void printPBEPARM(PBEparm *pbeparm) {

    int i;
    double ionstr = 0.0;

    for (i=0; i<pbeparm->nion; i++)
        ionstr += 0.5*(VSQR(pbeparm->ionq[i])*pbeparm->ionc[i]);

    Vnm_tprint( 1, "  Molecule ID: %d\n", pbeparm->molid);
    switch (pbeparm->pbetype) {
        case PBE_NPBE:
            Vnm_tprint( 1, "  Nonlinear traditional PBE\n");
            break;
        case PBE_LPBE:
            Vnm_tprint( 1, "  Linearized traditional PBE\n");
            break;
        case PBE_NRPBE:
            Vnm_tprint( 1, "  Nonlinear regularized PBE\n");
            Vnm_tprint( 2, "  ** Sorry, but Nathan broke the nonlinear regularized PBE implementation. **\n");
            Vnm_tprint( 2, "  ** Please let us know if you are interested in using it. **\n");
            VASSERT(0);
            break;
        case PBE_LRPBE:
            Vnm_tprint( 1, "  Linearized regularized PBE\n");
            break;
        case PBE_SMPBE: /* SMPBE Added */
            Vnm_tprint( 1, "  Nonlinear Size-Modified PBE\n");
            break;
        default:
            Vnm_tprint(2, "  Unknown PBE type (%d)!\n", pbeparm->pbetype);
            break;
    }
    if (pbeparm->bcfl == BCFL_ZERO) {
        Vnm_tprint( 1, "  Zero boundary conditions\n");
    } else if (pbeparm->bcfl == BCFL_SDH) {
        Vnm_tprint( 1, "  Single Debye-Huckel sphere boundary \
conditions\n");
    } else if (pbeparm->bcfl == BCFL_MDH) {
        Vnm_tprint( 1, "  Multiple Debye-Huckel sphere boundary \
conditions\n");
    } else if (pbeparm->bcfl == BCFL_FOCUS) {
        Vnm_tprint( 1, "  Boundary conditions from focusing\n");
    } else if (pbeparm->bcfl == BCFL_MAP) {
        Vnm_tprint( 1, "  Boundary conditions from potential map\n");
    } else if (pbeparm->bcfl == BCFL_MEM) {
        Vnm_tprint( 1, "  Membrane potential boundary conditions.\n");
    }
    Vnm_tprint( 1, "  %d ion species (%4.3f M ionic strength):\n",
                pbeparm->nion, ionstr);
    for (i=0; i<pbeparm->nion; i++) {
        Vnm_tprint( 1, "    %4.3f A-radius, %4.3f e-charge, \
%4.3f M concentration\n",
                    pbeparm->ionr[i], pbeparm->ionq[i], pbeparm->ionc[i]);
    }

    if(pbeparm->pbetype == PBE_SMPBE){ /* SMPBE Added */
        Vnm_tprint( 1, "  Lattice spacing: %4.3f A (SMPBE) \n", pbeparm->smvolume);
        Vnm_tprint( 1, "  Relative size parameter: %4.3f  (SMPBE) \n", pbeparm->smsize);
    }

    Vnm_tprint( 1, "  Solute dielectric: %4.3f\n", pbeparm->pdie);
    Vnm_tprint( 1, "  Solvent dielectric: %4.3f\n", pbeparm->sdie);
    switch (pbeparm->srfm) {
        case 0:
            Vnm_tprint( 1, "  Using \"molecular\" surface \
definition; no smoothing\n");
            Vnm_tprint( 1, "  Solvent probe radius: %4.3f A\n",
                        pbeparm->srad);
            break;
        case 1:
            Vnm_tprint( 1, "  Using \"molecular\" surface definition;\
harmonic average smoothing\n");
            Vnm_tprint( 1, "  Solvent probe radius: %4.3f A\n",
                        pbeparm->srad);
            break;
        case 2:
            Vnm_tprint( 1, "  Using spline-based surface definition;\
window = %4.3f\n", pbeparm->swin);
            break;
        default:
            break;
    }
    Vnm_tprint( 1, "  Temperature:  %4.3f K\n", pbeparm->temp);
    if (pbeparm->calcenergy != PCE_NO) Vnm_tprint( 1, "  Electrostatic \
energies will be calculated\n");
    if (pbeparm->calcforce == PCF_TOTAL) Vnm_tprint( 1, "  Net solvent \
forces will be calculated \n");
    if (pbeparm->calcforce == PCF_COMPS) Vnm_tprint( 1, "  All-atom \
solvent forces will be calculated\n");
    for (i=0; i<pbeparm->numwrite; i++) {
        switch (pbeparm->writetype[i]) {
            case VDT_CHARGE:
                Vnm_tprint(1, "  Charge distribution to be written to ");
                break;
            case VDT_POT:
                Vnm_tprint(1, "  Potential to be written to ");
                break;
            case VDT_SMOL:
                Vnm_tprint(1, "  Molecular solvent accessibility \
to be written to ");
                break;
            case VDT_SSPL:
                Vnm_tprint(1, "  Spline-based solvent accessibility \
to be written to ");
                break;
            case VDT_VDW:
                Vnm_tprint(1, "  van der Waals solvent accessibility \
to be written to ");
                break;
            case VDT_IVDW:
                Vnm_tprint(1, "  Ion accessibility to be written to ");
                break;
            case VDT_LAP:
                Vnm_tprint(1, "  Potential Laplacian to be written to ");
                break;
            case VDT_EDENS:
                Vnm_tprint(1, "  Energy density to be written to ");
                break;
            case VDT_NDENS:
                Vnm_tprint(1, "  Ion number density to be written to ");
                break;
            case VDT_QDENS:
                Vnm_tprint(1, "  Ion charge density to be written to ");
                break;
            case VDT_DIELX:
                Vnm_tprint(1, "  X-shifted dielectric map to be written \
to ");
                break;
            case VDT_DIELY:
                Vnm_tprint(1, "  Y-shifted dielectric map to be written \
to ");
                break;
            case VDT_DIELZ:
                Vnm_tprint(1, "  Z-shifted dielectric map to be written \
to ");
                break;
            case VDT_KAPPA:
                Vnm_tprint(1, "  Kappa map to be written to ");
                break;
            case VDT_ATOMPOT:
                Vnm_tprint(1, "  Atom potentials to be written to ");
                break;
            default:
                Vnm_tprint(2, "  Invalid data type for writing!\n");
                break;
        }
        switch (pbeparm->writefmt[i]) {
            case VDF_DX:
                Vnm_tprint(1, "%s.%s\n", pbeparm->writestem[i], "dx");
                break;
            case VDF_GZ:
                Vnm_tprint(1, "%s.%s\n", pbeparm->writestem[i], "dx.gz");
                break;
            case VDF_UHBD:
                Vnm_tprint(1, "%s.%s\n", pbeparm->writestem[i], "grd");
                break;
            case VDF_AVS:
                Vnm_tprint(1, "%s.%s\n", pbeparm->writestem[i], "ucd");
                break;
            case VDF_MCSF:
                Vnm_tprint(1, "%s.%s\n", pbeparm->writestem[i], "mcsf");
                break;
            case VDF_FLAT:
                Vnm_tprint(1, "%s.%s\n", pbeparm->writestem[i], "txt");
                break;
            default:
                Vnm_tprint(2, "  Invalid format for writing!\n");
                break;
        }

    }

}

VPUBLIC void printMGPARM(MGparm *mgparm, double realCenter[3]) {

    switch (mgparm->chgm) {
        case 0:
            Vnm_tprint(1, "  Using linear spline charge discretization.\n");
            break;
        case 1:
            Vnm_tprint(1, "  Using cubic spline charge discretization.\n");
            break;
        default:
            break;
    }
    if (mgparm->type == MCT_PARALLEL) {
        Vnm_tprint( 1, "  Partition overlap fraction = %g\n",
                    mgparm->ofrac);
        Vnm_tprint( 1, "  Processor array = %d x %d x %d\n",
                    mgparm->pdime[0], mgparm->pdime[1], mgparm->pdime[2]);
    }
    Vnm_tprint( 1, "  Grid dimensions: %d x %d x %d\n",
                mgparm->dime[0], mgparm->dime[1], mgparm->dime[2]);
    Vnm_tprint( 1, "  Grid spacings: %4.3f x %4.3f x %4.3f\n",
                mgparm->grid[0], mgparm->grid[1], mgparm->grid[2]);
    Vnm_tprint( 1, "  Grid lengths: %4.3f x %4.3f x %4.3f\n",
                mgparm->glen[0], mgparm->glen[1], mgparm->glen[2]);
    Vnm_tprint( 1, "  Grid center: (%4.3f, %4.3f, %4.3f)\n",
                realCenter[0], realCenter[1], realCenter[2]);
    Vnm_tprint( 1, "  Multigrid levels: %d\n", mgparm->nlev);

}

/**
 * Initialize a multigrid calculation.
 */
VPUBLIC int initMG(int icalc,
                   NOsh *nosh, MGparm *mgparm,
                   PBEparm *pbeparm,
                   double realCenter[3],
                   Vpbe *pbe[NOSH_MAXCALC],
                   Valist *alist[NOSH_MAXMOL],
                   Vgrid *dielXMap[NOSH_MAXMOL],
                   Vgrid *dielYMap[NOSH_MAXMOL],
                   Vgrid *dielZMap[NOSH_MAXMOL],
                   Vgrid *kappaMap[NOSH_MAXMOL],
                   Vgrid *chargeMap[NOSH_MAXMOL],
                   Vpmgp *pmgp[NOSH_MAXCALC],
                   Vpmg *pmg[NOSH_MAXCALC],
                   Vgrid *potMap[NOSH_MAXMOL]
                  ) {

    int j,
        focusFlag,
        iatom;
    size_t bytesTotal,
           highWater;
    double sparm,
           iparm,
           q;
    Vatom *atom = VNULL;
    Vgrid *theDielXMap = VNULL,
          *theDielYMap = VNULL,
          *theDielZMap = VNULL;
    Vgrid *theKappaMap = VNULL,
          *thePotMap = VNULL,
          *theChargeMap = VNULL;
    Valist *myalist = VNULL;

    Vnm_tstart(APBS_TIMER_SETUP, "Setup timer");

    /* Update the grid center */
    for (j=0; j<3; j++) realCenter[j] = mgparm->center[j];

    /* Check for completely-neutral molecule */
    q = 0;
    myalist = alist[pbeparm->molid-1];
    for (iatom=0; iatom<Valist_getNumberAtoms(myalist); iatom++) {
        atom = Valist_getAtom(myalist, iatom);
        q += VSQR(Vatom_getCharge(atom));
    }
    /*  D. Gohara 10/22/09 - disabled
    if (q < (1e-6)) {
        Vnm_tprint(2, "Molecule #%d is uncharged!\n", pbeparm->molid);
        Vnm_tprint(2, "Sum square charge = %g!\n", q);
        return 0;
    }
    */

    /* Set up PBE object */
    Vnm_tprint(0, "Setting up PBE object...\n");
    if (pbeparm->srfm == VSM_SPLINE) {
        sparm = pbeparm->swin;
    } else {
        sparm = pbeparm->srad;
    }
    if (pbeparm->nion > 0) {
        iparm = pbeparm->ionr[0];
    } else {
        iparm = 0.0;
    }
    if (pbeparm->bcfl == BCFL_FOCUS) {
        if (icalc == 0) {
            Vnm_tprint( 2, "Can't focus first calculation!\n");
            return 0;
        }
        focusFlag = 1;
    } else {
        focusFlag = 0;
    }

    // Construct Vpbe object
    pbe[icalc] = Vpbe_ctor(myalist, pbeparm->nion,
                           pbeparm->ionc, pbeparm->ionr, pbeparm->ionq,
                           pbeparm->temp, pbeparm->pdie,
                           pbeparm->sdie, sparm, focusFlag, pbeparm->sdens,
                           pbeparm->zmem, pbeparm->Lmem, pbeparm->mdie,
                           pbeparm->memv);

    /* Set up PDE object */
    Vnm_tprint(0, "Setting up PDE object...\n");
    switch (pbeparm->pbetype) {
        case PBE_NPBE:
            /* TEMPORARY USEAQUA */
            mgparm->nonlintype = NONLIN_NPBE;
            mgparm->method = (mgparm->useAqua == 1) ? VSOL_NewtonAqua : VSOL_Newton;
            pmgp[icalc] = Vpmgp_ctor(mgparm);
            break;
        case PBE_LPBE:
            /* TEMPORARY USEAQUA */
            mgparm->nonlintype = NONLIN_LPBE;
            mgparm->method = (mgparm->useAqua == 1) ? VSOL_CGMGAqua : VSOL_MG;
            pmgp[icalc] = Vpmgp_ctor(mgparm);
            break;
        case PBE_LRPBE:
            Vnm_tprint(2, "Sorry, LRPBE isn't supported with the MG solver!\n");
            return 0;
        case PBE_NRPBE:
            Vnm_tprint(2, "Sorry, NRPBE isn't supported with the MG solver!\n");
            return 0;
        case PBE_SMPBE: /* SMPBE Added */
            mgparm->nonlintype = NONLIN_SMPBE;
            pmgp[icalc] = Vpmgp_ctor(mgparm);

            /* Copy Code */
            pbe[icalc]->smsize = pbeparm->smsize;
            pbe[icalc]->smvolume = pbeparm->smvolume;
            pbe[icalc]->ipkey = pmgp[icalc]->ipkey;

            break;
        default:
            Vnm_tprint(2, "Error!  Unknown PBE type (%d)!\n", pbeparm->pbetype);
            return 0;
    }
    Vnm_tprint(0, "Setting PDE center to local center...\n");
    pmgp[icalc]->bcfl = pbeparm->bcfl;
    pmgp[icalc]->xcent = realCenter[0];
    pmgp[icalc]->ycent = realCenter[1];
    pmgp[icalc]->zcent = realCenter[2];

    if (pbeparm->bcfl == BCFL_FOCUS) {
        if (icalc == 0) {
            Vnm_tprint( 2, "Can't focus first calculation!\n");
            return 0;
        }
        /* Focusing requires the previous calculation in order to setup the
        current run... */
        pmg[icalc] = Vpmg_ctor(pmgp[icalc], pbe[icalc], 1, pmg[icalc-1],
                               mgparm, pbeparm->calcenergy);
        /* ...however, it should be done with the previous calculation now, so
        we should be able to destroy it here. */
        /* Vpmg_dtor(&(pmg[icalc-1])); */
    } else {
        if (icalc>0) Vpmg_dtor(&(pmg[icalc-1]));
        pmg[icalc] = Vpmg_ctor(pmgp[icalc], pbe[icalc], 0, VNULL, mgparm, PCE_NO);
    }
    if (icalc>0) {
        Vpmgp_dtor(&(pmgp[icalc-1]));
        Vpbe_dtor(&(pbe[icalc-1]));
    }
    if (pbeparm->useDielMap) {
        if ((pbeparm->dielMapID-1) < nosh->ndiel) {
            theDielXMap = dielXMap[pbeparm->dielMapID-1];
        } else {
            Vnm_print(2, "Error!  %d is not a valid dielectric map ID!\n",
                      pbeparm->dielMapID);
            return 0;
        }
    }
    if (pbeparm->useDielMap) {
        if ((pbeparm->dielMapID-1) < nosh->ndiel) {
            theDielYMap = dielYMap[pbeparm->dielMapID-1];
        } else {
            Vnm_print(2, "Error!  %d is not a valid dielectric map ID!\n",
                      pbeparm->dielMapID);
            return 0;
        }
    }
    if (pbeparm->useDielMap) {
        if ((pbeparm->dielMapID-1) < nosh->ndiel) {
            theDielZMap = dielZMap[pbeparm->dielMapID-1];
        } else {
            Vnm_print(2, "Error!  %d is not a valid dielectric map ID!\n",
                      pbeparm->dielMapID);
            return 0;
        }
    }
    if (pbeparm->useKappaMap) {
        if ((pbeparm->kappaMapID-1) < nosh->nkappa) {
            theKappaMap = kappaMap[pbeparm->kappaMapID-1];
        } else {
            Vnm_print(2, "Error!  %d is not a valid kappa map ID!\n",
                      pbeparm->kappaMapID);
            return 0;
        }
    }
    if (pbeparm->usePotMap) {
        if ((pbeparm->potMapID-1) < nosh->npot) {
            thePotMap = potMap[pbeparm->potMapID-1];
        } else {
            Vnm_print(2, "Error!  %d is not a valid potential map ID!\n",
                      pbeparm->potMapID);
            return 0;
        }
    }
    if (pbeparm->useChargeMap) {
        if ((pbeparm->chargeMapID-1) < nosh->ncharge) {
            theChargeMap = chargeMap[pbeparm->chargeMapID-1];
        } else {
            Vnm_print(2, "Error!  %d is not a valid charge map ID!\n",
                      pbeparm->chargeMapID);
            return 0;
        }
    }

    if (pbeparm->bcfl == BCFL_MAP && thePotMap == VNULL) {
        Vnm_print(2, "Warning: You specified 'bcfl map' in the input file, but no potential map was found.\n");
        Vnm_print(2, "         You must specify 'usemap pot' statement in the APBS input file!\n");
        Vnm_print(2, "Bailing out ...\n");
        return 0;
    }

    // Initialize calculation coefficients
    if (!Vpmg_fillco(pmg[icalc],
                     pbeparm->srfm, pbeparm->swin, mgparm->chgm,
                     pbeparm->useDielMap, theDielXMap,
                     pbeparm->useDielMap, theDielYMap,
                     pbeparm->useDielMap, theDielZMap,
                     pbeparm->useKappaMap, theKappaMap,
                     pbeparm->usePotMap, thePotMap,
                     pbeparm->useChargeMap, theChargeMap)) {
        Vnm_print(2, "initMG:  problems setting up coefficients (fillco)!\n");
        return 0;
    }

    /* Print a few derived parameters */
#ifndef VAPBSQUIET
    Vnm_tprint(1, "  Debye length:  %g A\n", Vpbe_getDeblen(pbe[icalc]));
#endif

    /* Setup time statistics */
    Vnm_tstop(APBS_TIMER_SETUP, "Setup timer");

    /* Memory statistics */
    bytesTotal = Vmem_bytesTotal();
    highWater = Vmem_highWaterTotal();

#ifndef VAPBSQUIET
    Vnm_tprint( 1, "  Current memory usage:  %4.3f MB total, \
%4.3f MB high water\n", (double)(bytesTotal)/(1024.*1024.),
                (double)(highWater)/(1024.*1024.));
#endif

    return 1;

}

VPUBLIC void killMG(NOsh *nosh, Vpbe *pbe[NOSH_MAXCALC],
                    Vpmgp *pmgp[NOSH_MAXCALC], Vpmg *pmg[NOSH_MAXCALC]) {

        int i;

#ifndef VAPBSQUIET
    Vnm_tprint(1, "Destroying multigrid structures.\n");
#endif

    /*
       There appears to be a relationship (or this is a bug in Linux, can't tell
       at the moment, since Linux is the only OS that seems to be affected)
       between one of the three object types: Vpbe, Vpmg or Vpmgp that requires
       deallocations to be performed in a specific order. This results in a
       bug some of the time when freeing Vpmg objects below. Therefore it
       appears to be important to release the Vpmg structs BEFORE the Vpmgp structs .
    */
    Vpmg_dtor(&(pmg[nosh->ncalc-1]));

    for(i=0;i<nosh->ncalc;i++){
        Vpbe_dtor(&(pbe[i]));
        Vpmgp_dtor(&(pmgp[i]));
    }

}

VPUBLIC int solveMG(NOsh *nosh,
                    Vpmg *pmg,
                    MGparm_CalcType type
                   ) {

    int nx,
        ny,
        nz,
        i;

    if (nosh != VNULL) {
        if (nosh->bogus) return 1;
    }

    Vnm_tstart(APBS_TIMER_SOLVER, "Solver timer");


    if (type != MCT_DUMMY) {
        if (!Vpmg_solve(pmg)) {
            Vnm_print(2, "  Error during PDE solution!\n");
            return 0;
        }
    } else {
        Vnm_tprint( 1,"  Skipping solve for mg-dummy run; zeroing \
solution array\n");
        nx = pmg->pmgp->nx;
        ny = pmg->pmgp->ny;
        nz = pmg->pmgp->nz;
        for (i=0; i<nx*ny*nz; i++) pmg->u[i] = 0.0;
    }
    Vnm_tstop(APBS_TIMER_SOLVER, "Solver timer");

    return 1;

}

VPUBLIC int setPartMG(NOsh *nosh,
                      MGparm *mgparm,
                      Vpmg *pmg
                     ) {

    int j;
    double partMin[3],
           partMax[3];

    if (nosh->bogus) return 1;

    if (mgparm->type == MCT_PARALLEL) {
        for (j=0; j<3; j++) {
            partMin[j] = mgparm->partDisjCenter[j] - 0.5*mgparm->partDisjLength[j];
            partMax[j] = mgparm->partDisjCenter[j] + 0.5*mgparm->partDisjLength[j];
        }
#if 0
        Vnm_tprint(1, "setPartMG (%s, %d):  Disj part center = (%g, %g, %g)\n",
                   __FILE__, __LINE__,
                   mgparm->partDisjCenter[0],
                   mgparm->partDisjCenter[1],
                   mgparm->partDisjCenter[2]
                   );
        Vnm_tprint(1, "setPartMG (%s, %d):  Disj part lower corner = (%g, %g, %g)\n",
                   __FILE__, __LINE__, partMin[0], partMin[1], partMin[2]);
        Vnm_tprint(1, "setPartMG (%s, %d):  Disj part upper corner = (%g, %g, %g)\n",
                   __FILE__, __LINE__,
                   partMax[0], partMax[1], partMax[2]);
#endif
    } else {
        for (j=0; j<3; j++) {
            partMin[j] = mgparm->center[j] - 0.5*mgparm->glen[j];
            partMax[j] = mgparm->center[j] + 0.5*mgparm->glen[j];
        }
    }
    /* Vnm_print(1, "DEBUG (%s, %d):  setPartMG calling setPart with upper corner \
%g %g %g and lower corner %g %g %g\n", __FILE__,  __LINE__,
              partMin[0], partMin[1], partMin[2],
              partMax[0], partMax[1], partMax[2]); */
    Vpmg_setPart(pmg, partMin, partMax, mgparm->partDisjOwnSide);


    return 1;

}

VPUBLIC int energyMG(NOsh *nosh,
                     int icalc,
                     Vpmg *pmg,
                     int *nenergy,
                     double *totEnergy,
                     double *qfEnergy,
                     double *qmEnergy,
                     double *dielEnergy
                    ) {

    Valist *alist;
    Vatom *atom;
    int i,
        extEnergy;
    double tenergy;
    MGparm *mgparm;
    PBEparm *pbeparm;

    mgparm = nosh->calc[icalc]->mgparm;
    pbeparm = nosh->calc[icalc]->pbeparm;

    Vnm_tstart(APBS_TIMER_ENERGY, "Energy timer");
    extEnergy = 1;

    if (pbeparm->calcenergy == PCE_TOTAL) {
        *nenergy = 1;
        /* Some processors don't count */
        if (nosh->bogus == 0) {
            *totEnergy = Vpmg_energy(pmg, extEnergy);
#ifndef VAPBSQUIET
            Vnm_tprint( 1, "  Total electrostatic energy = %1.12E kJ/mol\n",
                        Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*totEnergy));
#endif
        } else *totEnergy = 0;
    } else if (pbeparm->calcenergy == PCE_COMPS) {
        *nenergy = 1;
        *totEnergy = Vpmg_energy(pmg, extEnergy);
        *qfEnergy = Vpmg_qfEnergy(pmg, extEnergy);
        *qmEnergy = Vpmg_qmEnergy(pmg, extEnergy);
        *dielEnergy = Vpmg_dielEnergy(pmg, extEnergy);
#ifndef VAPBSQUIET
        Vnm_tprint( 1, "  Total electrostatic energy = %1.12E \
kJ/mol\n", Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*totEnergy));
        Vnm_tprint( 1, "  Fixed charge energy = %g kJ/mol\n",
                    0.5*Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*qfEnergy));
        Vnm_tprint( 1, "  Mobile charge energy = %g kJ/mol\n",
                    Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*qmEnergy));
        Vnm_tprint( 1, "  Dielectric energy = %g kJ/mol\n",
                    Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*dielEnergy));
        Vnm_tprint( 1, "  Per-atom energies:\n");
#endif
        alist = pmg->pbe->alist;
        for (i=0; i<Valist_getNumberAtoms(alist); i++) {
            atom = Valist_getAtom(alist, i);
            tenergy = Vpmg_qfAtomEnergy(pmg, atom);
#ifndef VAPBSQUIET
            Vnm_tprint( 1, "      Atom %d:  %1.12E kJ/mol\n", i,
                        0.5*Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*tenergy);
#endif
        }
    } else *nenergy = 0;

    Vnm_tstop(APBS_TIMER_ENERGY, "Energy timer");

    return 1;
}

VPUBLIC int forceMG(Vmem *mem,
                    NOsh *nosh,
                    PBEparm *pbeparm,
                    MGparm *mgparm,
                    Vpmg *pmg,
                    int *nforce,
                    AtomForce **atomForce,
                    Valist *alist[NOSH_MAXMOL]
                   ) {

    int j,
        k;
    double qfForce[3],
           dbForce[3],
           ibForce[3];

    Vnm_tstart(APBS_TIMER_FORCE, "Force timer");

#ifndef VAPBSQUIET
    Vnm_tprint( 1,"  Calculating forces...\n");
#endif

    if (pbeparm->calcforce == PCF_TOTAL) {
        *nforce = 1;
        *atomForce = (AtomForce *)Vmem_malloc(mem, 1, sizeof(AtomForce));
        /* Clear out force arrays */
        for (j=0; j<3; j++) {
            (*atomForce)[0].qfForce[j] = 0;
            (*atomForce)[0].ibForce[j] = 0;
            (*atomForce)[0].dbForce[j] = 0;
        }
        for (j=0;j<Valist_getNumberAtoms(alist[pbeparm->molid-1]);j++) {
            if (nosh->bogus == 0) {
                VASSERT(Vpmg_qfForce(pmg, qfForce, j, mgparm->chgm));
                VASSERT(Vpmg_ibForce(pmg, ibForce, j, pbeparm->srfm));
                VASSERT(Vpmg_dbForce(pmg, dbForce, j, pbeparm->srfm));
            } else {
                for (k=0; k<3; k++) {
                    qfForce[k] = 0;
                    ibForce[k] = 0;
                    dbForce[k] = 0;
                }
            }
            for (k=0; k<3; k++) {
                (*atomForce)[0].qfForce[k] += qfForce[k];
                (*atomForce)[0].ibForce[k] += ibForce[k];
                (*atomForce)[0].dbForce[k] += dbForce[k];
            }
        }
#ifndef VAPBSQUIET
        Vnm_tprint( 1, "  Printing net forces for molecule %d (kJ/mol/A)\n",
                    pbeparm->molid);
        Vnm_tprint( 1, "  Legend:\n");
        Vnm_tprint( 1, "    qf  -- fixed charge force\n");
        Vnm_tprint( 1, "    db  -- dielectric boundary force\n");
        Vnm_tprint( 1, "    ib  -- ionic boundary force\n");
        Vnm_tprint( 1, "  qf  %4.3e  %4.3e  %4.3e\n",
                    Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].qfForce[0],
                    Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].qfForce[1],
                    Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].qfForce[2]);
        Vnm_tprint( 1, "  ib  %4.3e  %4.3e  %4.3e\n",
                    Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].ibForce[0],
                    Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].ibForce[1],
                    Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].ibForce[2]);
        Vnm_tprint( 1, "  db  %4.3e  %4.3e  %4.3e\n",
                    Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].dbForce[0],
                    Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].dbForce[1],
                    Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*atomForce)[0].dbForce[2]);
#endif
    } else if (pbeparm->calcforce == PCF_COMPS) {
        *nforce = Valist_getNumberAtoms(alist[pbeparm->molid-1]);
        *atomForce = (AtomForce *)Vmem_malloc(mem, *nforce,
                                              sizeof(AtomForce));
#ifndef VAPBSQUIET
        Vnm_tprint( 1, "  Printing per-atom forces for molecule %d (kJ/mol/A)\n",
                    pbeparm->molid);
        Vnm_tprint( 1, "  Legend:\n");
        Vnm_tprint( 1, "    tot n -- total force for atom n\n");
        Vnm_tprint( 1, "    qf  n -- fixed charge force for atom n\n");
        Vnm_tprint( 1, "    db  n -- dielectric boundary force for atom n\n");
        Vnm_tprint( 1, "    ib  n -- ionic boundary force for atom n\n");
#endif
        for (j=0;j<Valist_getNumberAtoms(alist[pbeparm->molid-1]);j++) {
            if (nosh->bogus == 0) {
                VASSERT(Vpmg_qfForce(pmg, (*atomForce)[j].qfForce, j,
                                     mgparm->chgm));
                VASSERT(Vpmg_ibForce(pmg, (*atomForce)[j].ibForce, j,
                                     pbeparm->srfm));
                VASSERT(Vpmg_dbForce(pmg, (*atomForce)[j].dbForce, j,
                                     pbeparm->srfm));
            } else {
                for (k=0; k<3; k++) {
                    (*atomForce)[j].qfForce[k] = 0;
                    (*atomForce)[j].ibForce[k] = 0;
                    (*atomForce)[j].dbForce[k] = 0;
                }
            }
#ifndef VAPBSQUIET
            Vnm_tprint( 1, "mgF  tot %d  %4.3e  %4.3e  %4.3e\n", j,
                        Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na \
                        *((*atomForce)[j].qfForce[0]+(*atomForce)[j].ibForce[0]+
                          (*atomForce)[j].dbForce[0]),
                        Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na \
                        *((*atomForce)[j].qfForce[1]+(*atomForce)[j].ibForce[1]+
                          (*atomForce)[j].dbForce[1]),
                        Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na \
                        *((*atomForce)[j].qfForce[2]+(*atomForce)[j].ibForce[2]+
                          (*atomForce)[j].dbForce[2]));
            Vnm_tprint( 1, "mgF  qf  %d  %4.3e  %4.3e  %4.3e\n", j,
                        Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na \
                        *(*atomForce)[j].qfForce[0],
                        Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na \
                        *(*atomForce)[j].qfForce[1],
                        Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na \
                        *(*atomForce)[j].qfForce[2]);
            Vnm_tprint( 1, "mgF  ib  %d  %4.3e  %4.3e  %4.3e\n", j,
                        Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na \
                        *(*atomForce)[j].ibForce[0],
                        Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na \
                        *(*atomForce)[j].ibForce[1],
                        Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na \
                        *(*atomForce)[j].ibForce[2]);
            Vnm_tprint( 1, "mgF  db  %d  %4.3e  %4.3e  %4.3e\n", j,
                        Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na \
                        *(*atomForce)[j].dbForce[0],
                        Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na \
                        *(*atomForce)[j].dbForce[1],
                        Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na \
                        *(*atomForce)[j].dbForce[2]);
#endif
        }
    } else *nforce = 0;

    Vnm_tstop(APBS_TIMER_FORCE, "Force timer");

    return 1;
}

VPUBLIC void killEnergy() {

#ifndef VAPBSQUIET
    Vnm_tprint(1, "No energy arrays to destroy.\n");
#endif

}

VPUBLIC void killForce(Vmem *mem, NOsh *nosh, int nforce[NOSH_MAXCALC],
                       AtomForce *atomForce[NOSH_MAXCALC]) {

    int i;

#ifndef VAPBSQUIET
    Vnm_tprint(1, "Destroying force arrays.\n");
#endif

    for (i=0; i<nosh->ncalc; i++) {

        if (nforce[i] > 0) Vmem_free(mem, nforce[i], sizeof(AtomForce),
                                     (void **)&(atomForce[i]));

    }
}

VPUBLIC int writematMG(int rank, NOsh *nosh, PBEparm *pbeparm, Vpmg *pmg) {

    char writematstem[VMAX_ARGLEN];
    char outpath[VMAX_ARGLEN];
    char mxtype[3];
    int strlenmax;

    if (nosh->bogus) return 1;

#ifdef HAVE_MPI_H
    strlenmax = VMAX_ARGLEN-14;
    if (strlen(pbeparm->writematstem) > strlenmax) {
        Vnm_tprint(2, "  Matrix name (%s) too long (%d char max)!\n",
                   pbeparm->writematstem, strlenmax);
        Vnm_tprint(2, "  Not writing matrix!\n");
        return 0;
    }
    sprintf(writematstem, "%s-PE%d", pbeparm->writematstem, rank);
#else
    strlenmax = (int)(VMAX_ARGLEN)-1;
    if ((int)strlen(pbeparm->writematstem) > strlenmax) {
        Vnm_tprint(2, "  Matrix name (%s) too long (%d char max)!\n",
                   pbeparm->writematstem, strlenmax);
        Vnm_tprint(2, "  Not writing matrix!\n");
        return 0;
    }
    if(nosh->ispara == 1){
        sprintf(writematstem, "%s-PE%d", pbeparm->writematstem,nosh->proc_rank);
    }else{
        sprintf(writematstem, "%s", pbeparm->writematstem);
    }
#endif

    if (pbeparm->writemat == 1) {
        strlenmax = VMAX_ARGLEN-5;
        if ((int)strlen(pbeparm->writematstem) > strlenmax) {
            Vnm_tprint(2, "  Matrix name (%s) too long (%d char max)!\n",
                       pbeparm->writematstem, strlenmax);
            Vnm_tprint(2, "  Not writing matrix!\n");
            return 0;
        }
        sprintf(outpath, "%s.%s", writematstem, "mat");
        mxtype[0] = 'R';
        mxtype[1] = 'S';
        mxtype[2] = 'A';
        /* Poisson operator only */
        if (pbeparm->writematflag == 0) {
            Vnm_tprint( 1, "  Writing Poisson operator matrix \
to %s...\n", outpath);

            /* Linearization of Poisson-Boltzmann operator around solution */
        } else if (pbeparm->writematflag == 1) {
            Vnm_tprint( 1, "  Writing linearization of full \
Poisson-Boltzmann operator matrix to %s...\n", outpath);

        } else {
            Vnm_tprint( 2, "  Bogus matrix specification\
(%d)!\n", pbeparm->writematflag);
            return 0;
        }

        Vnm_tprint(0, "  Printing operator...\n");
        //Vpmg_printColComp(pmg, outpath, outpath, mxtype,
        //				  pbeparm->writematflag);
        return 0;

    }

    return 1;
}

VPUBLIC void storeAtomEnergy(Vpmg *pmg, int icalc, double **atomEnergy,
                             int *nenergy){

    Vatom *atom;
    Valist *alist;
    int i;

    alist = pmg->pbe->alist;
    *nenergy = Valist_getNumberAtoms(alist);
    *atomEnergy = (double *)Vmem_malloc(pmg->vmem, *nenergy, sizeof(double));

    for (i=0; i<*nenergy; i++) {
        atom = Valist_getAtom(alist, i);
        (*atomEnergy)[i] = Vpmg_qfAtomEnergy(pmg, atom);
    }
}

VPUBLIC int writedataFlat(
                          NOsh *nosh,
                          Vcom *com,
                          const char *fname,
                          double totEnergy[NOSH_MAXCALC],
                          double qfEnergy[NOSH_MAXCALC],
                          double qmEnergy[NOSH_MAXCALC],
                          double dielEnergy[NOSH_MAXCALC],
                          int nenergy[NOSH_MAXCALC],
                          double *atomEnergy[NOSH_MAXCALC],
                          int nforce[NOSH_MAXCALC],
                          AtomForce *atomForce[NOSH_MAXCALC]) {

    FILE *file;
    time_t now;
    int ielec, icalc, i, j;
    char *timestring = VNULL;
    PBEparm *pbeparm = VNULL;
    MGparm *mgparm = VNULL;
    double conversion, ltenergy, gtenergy, scalar;

    if (nosh->bogus) return 1;

    /* Initialize some variables */

    icalc = 0;

    file = fopen(fname, "w");
    if (file == VNULL) {
        Vnm_print(2, "writedataFlat: Problem opening virtual socket %s\n",
                  fname);
        return 0;
    }

    /* Strip the newline character from the date */

    now = time(VNULL);
    timestring = ctime(&now);
    fprintf(file,"%s\n", timestring);

    for (ielec=0; ielec<nosh->nelec;ielec++) { /* elec loop */

        /* Initialize per-elec pointers */

        mgparm = nosh->calc[icalc]->mgparm;
        pbeparm = nosh->calc[icalc]->pbeparm;

        /* Convert from kT/e to kJ/mol */
        conversion =  Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na;

        fprintf(file,"elec");
        if (Vstring_strcasecmp(nosh->elecname[ielec], "") != 0) {
            fprintf(file," name %s\n", nosh->elecname[ielec]);
        } else fprintf(file, "\n");

        switch (mgparm->type) {
            case MCT_DUMMY:
                fprintf(file,"    mg-dummy\n");
                break;
            case MCT_MANUAL:
                fprintf(file,"    mg-manual\n");
                break;
            case MCT_AUTO:
                fprintf(file,"    mg-auto\n");
                break;
            case MCT_PARALLEL:
                fprintf(file,"    mg-para\n");
                break;
            default:
                break;
        }

        fprintf(file,"    mol %d\n", pbeparm->molid);
        fprintf(file,"    dime %d %d %d\n", mgparm->dime[0], mgparm->dime[1],\
                mgparm->dime[2]);

        switch (pbeparm->pbetype) {
            case PBE_NPBE:
                fprintf(file,"    npbe\n");
                break;
            case PBE_LPBE:
                fprintf(file,"    lpbe\n");
                break;
            default:
                break;
        }

        if (pbeparm->nion > 0) {
            for (i=0; i<pbeparm->nion; i++) {
                fprintf(file,"    ion %4.3f %4.3f %4.3f\n",
                        pbeparm->ionr[i], pbeparm->ionq[i], pbeparm->ionc[i]);
            }
        }

        fprintf(file,"    pdie %4.3f\n", pbeparm->pdie);
        fprintf(file,"    sdie %4.3f\n", pbeparm->sdie);

        switch (pbeparm->srfm) {
            case 0:
                fprintf(file,"    srfm mol\n");
                fprintf(file,"    srad %4.3f\n", pbeparm->srad);
                break;
            case 1:
                fprintf(file,"    srfm smol\n");
                fprintf(file,"    srad %4.3f\n", pbeparm->srad);
                break;
            case 2:
                fprintf(file,"    srfm spl2\n");
                fprintf(file,"    srad %4.3f\n", pbeparm->srad);
                break;
            default:
                break;
        }

        switch (pbeparm->bcfl) {
            case BCFL_ZERO:
                fprintf(file,"    bcfl zero\n");
                break;
            case BCFL_SDH:
                fprintf(file,"    bcfl sdh\n");
                break;
            case BCFL_MDH:
                fprintf(file,"    bcfl mdh\n");
                break;
            case BCFL_FOCUS:
                fprintf(file,"    bcfl focus\n");
                break;
            case BCFL_MAP:
                fprintf(file,"    bcfl map\n");
                break;
            case BCFL_MEM:
                fprintf(file,"    bcfl mem\n");
                break;
            default:
                break;
        }

        fprintf(file,"    temp %4.3f\n", pbeparm->temp);

        for (;icalc<=nosh->elec2calc[ielec];icalc++){ /* calc loop */

            /* Reinitialize per-calc pointers */
            mgparm = nosh->calc[icalc]->mgparm;
            pbeparm = nosh->calc[icalc]->pbeparm;

            fprintf(file,"    calc\n");
            fprintf(file,"        id %i\n", (icalc+1));
            fprintf(file,"        grid %4.3f %4.3f %4.3f\n",
                    mgparm->grid[0], mgparm->grid[1], mgparm->grid[2]);
            fprintf(file,"        glen %4.3f %4.3f %4.3f\n",
                    mgparm->glen[0], mgparm->glen[1], mgparm->glen[2]);

            if (pbeparm->calcenergy == PCE_TOTAL) {
                fprintf(file,"        totEnergy %1.12E kJ/mol\n",
                        (totEnergy[icalc]*conversion));
            } if (pbeparm->calcenergy == PCE_COMPS) {
                    fprintf(file,"        totEnergy %1.12E kJ/mol\n",
                        (totEnergy[icalc]*conversion));
                fprintf(file,"        qfEnergy %1.12E kJ/mol\n",
                        (0.5*qfEnergy[icalc]*conversion));
                fprintf(file,"        qmEnergy %1.12E kJ/mol\n",
                        (qmEnergy[icalc]*conversion));
                fprintf(file,"        dielEnergy %1.12E kJ/mol\n",
                        (dielEnergy[icalc]*conversion));
                for (i=0; i<nenergy[icalc]; i++){
                    fprintf(file,"        atom %i %1.12E kJ/mol\n", i,
                            (0.5*atomEnergy[icalc][i]*conversion));

                }
            }

            if (pbeparm->calcforce == PCF_TOTAL) {
                fprintf(file,"        qfForce %1.12E %1.12E %1.12E kJ/mol/A\n",
                        (atomForce[icalc][0].qfForce[0]*conversion),
                            (atomForce[icalc][0].qfForce[1]*conversion),
                            (atomForce[icalc][0].qfForce[2]*conversion));
                fprintf(file,"        ibForce %1.12E %1.12E %1.12E kJ/mol/A\n",
                        (atomForce[icalc][0].ibForce[0]*conversion),
                            (atomForce[icalc][0].ibForce[1]*conversion),
                            (atomForce[icalc][0].ibForce[2]*conversion));
                fprintf(file,"        dbForce %1.12E %1.12E %1.12E kJ/mol/A\n",
                        (atomForce[icalc][0].dbForce[0]*conversion),
                            (atomForce[icalc][0].dbForce[1]*conversion),
                            (atomForce[icalc][0].dbForce[2]*conversion));
            }
            fprintf(file,"    end\n");
        }

        fprintf(file,"end\n");
    }

/* Handle print energy statements */

for (i=0; i<nosh->nprint; i++) {

    if (nosh->printwhat[i] == NPT_ENERGY) {

        fprintf(file,"print energy");
        fprintf(file," %d", nosh->printcalc[i][0]+1);

        for (j=1; j<nosh->printnarg[i]; j++) {
            if (nosh->printop[i][j-1] == 0) fprintf(file," +");
            else if (nosh->printop[i][j-1] == 1) fprintf(file, " -");
            fprintf(file, " %d", nosh->printcalc[i][j]+1);
        }

        fprintf(file, "\n");
        icalc = nosh->elec2calc[nosh->printcalc[i][0]];

        ltenergy = Vunit_kb * (1e-3) * Vunit_Na * \
            nosh->calc[icalc]->pbeparm->temp * totEnergy[icalc];

        for (j=1; j<nosh->printnarg[i]; j++) {
            icalc = nosh->elec2calc[nosh->printcalc[i][j]];
            /* Add or subtract? */
            if (nosh->printop[i][j-1] == 0) scalar = 1.0;
            else if (nosh->printop[i][j-1] == 1) scalar = -1.0;
            /* Accumulate */
            ltenergy += (scalar * Vunit_kb * (1e-3) * Vunit_Na *
                         nosh->calc[icalc]->pbeparm->temp * totEnergy[icalc]);

            Vcom_reduce(com, &ltenergy, &gtenergy, 1, 2, 0);

        }
        fprintf(file,"    localEnergy %1.12E kJ/mol\n", \
                ltenergy);
        fprintf(file,"    globalEnergy %1.12E kJ/mol\nend\n", \
                gtenergy);
    }
}

fclose(file);

return 1;
}

VPUBLIC int writedataXML(NOsh *nosh, Vcom *com, const char *fname,
                         double totEnergy[NOSH_MAXCALC],
                         double qfEnergy[NOSH_MAXCALC],
                         double qmEnergy[NOSH_MAXCALC],
                         double dielEnergy[NOSH_MAXCALC],
                         int nenergy[NOSH_MAXCALC],
                         double *atomEnergy[NOSH_MAXCALC],
                                     int nforce[NOSH_MAXCALC],
                         AtomForce *atomForce[NOSH_MAXCALC]) {

    FILE *file;
    time_t now;
    int ielec, icalc, i, j;
    char *timestring = VNULL;
    char *c = VNULL;
    PBEparm *pbeparm = VNULL;
    MGparm *mgparm = VNULL;
    double conversion, ltenergy, gtenergy, scalar;

    if (nosh->bogus) return 1;

    /* Initialize some variables */

    icalc = 0;

    file = fopen(fname, "w");
    if (file == VNULL) {
        Vnm_print(2, "writedataXML: Problem opening virtual socket %s\n",
                  fname);
        return 0;
    }

    fprintf(file,"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(file,"<APBS>\n");

    /* Strip the newline character from the date */

    now = time(VNULL);
    timestring = ctime(&now);
    for(c = timestring; *c != '\n'; c++);
    *c = '\0';
    fprintf(file,"    <date>%s</date>\n", timestring);

    for (ielec=0; ielec<nosh->nelec;ielec++){ /* elec loop */

        /* Initialize per-elec pointers */

        mgparm = nosh->calc[icalc]->mgparm;
        pbeparm = nosh->calc[icalc]->pbeparm;

        /* Convert from kT/e to kJ/mol */
        conversion =  Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na;

        fprintf(file,"    <elec>\n");
        if (Vstring_strcasecmp(nosh->elecname[ielec], "") != 0) {
            fprintf(file,"      <name>%s</name>\n", nosh->elecname[ielec]);
        }

        switch (mgparm->type) {
            case MCT_DUMMY:
                fprintf(file,"      <type>mg-dummy</type>\n");
                break;
            case MCT_MANUAL:
                fprintf(file,"      <type>mg-manual</type>\n");
                break;
            case MCT_AUTO:
                fprintf(file,"      <type>mg-auto</type>\n");
                break;
            case MCT_PARALLEL:
                fprintf(file,"      <type>mg-para</type>\n");
                break;
            default:
                break;
        }

        fprintf(file,"      <molid>%d</molid>\n", pbeparm->molid);
        fprintf(file,"      <nx>%d</nx>\n", mgparm->dime[0]);
        fprintf(file,"      <ny>%d</ny>\n", mgparm->dime[1]);
        fprintf(file,"      <nz>%d</nz>\n", mgparm->dime[2]);

        switch (pbeparm->pbetype) {
            case PBE_NPBE:
                fprintf(file,"      <pbe>npbe</pbe>\n");
                break;
            case PBE_LPBE:
                fprintf(file,"      <pbe>lpbe</pbe>\n");
                break;
            default:
                break;
        }

        if (pbeparm->nion > 0) {
            for (i=0; i<pbeparm->nion; i++) {
                fprintf(file, "      <ion>\n");
                fprintf(file,"          <radius>%4.3f A</radius>\n",
                        pbeparm->ionr[i]);
                fprintf(file,"          <charge>%4.3f A</charge>\n",
                        pbeparm->ionq[i]);
                fprintf(file,"          <concentration>%4.3f M</concentration>\n",
                        pbeparm->ionc[i]);
                fprintf(file, "      </ion>\n");

            }
        }

        fprintf(file,"      <pdie>%4.3f</pdie>\n", pbeparm->pdie);
        fprintf(file,"      <sdie>%4.3f</sdie>\n", pbeparm->sdie);

        switch (pbeparm->srfm) {
            case 0:
                fprintf(file,"      <srfm>mol</srfm>\n");
                fprintf(file,"      <srad>%4.3f</srad>\n", pbeparm->srad);
                break;
            case 1:
                fprintf(file,"      <srfm>smol</srfm>\n");
                fprintf(file,"      <srad>%4.3f</srad>\n", pbeparm->srad);
                break;
            case 2:
                fprintf(file,"      <srfm>spl2</srfm>\n");
                break;
            default:
                break;
        }

        switch (pbeparm->bcfl) {
            case BCFL_ZERO:
                fprintf(file,"      <bcfl>zero</bcfl>\n");
                break;
            case BCFL_SDH:
                fprintf(file,"      <bcfl>sdh</bcfl>\n");
                break;
            case BCFL_MDH:
                fprintf(file,"      <bcfl>mdh</bcfl>\n");
                break;
            case BCFL_FOCUS:
                fprintf(file,"      <bcfl>focus</bcfl>\n");
                break;
            case BCFL_MAP:
                fprintf(file,"      <bcfl>map</bcfl>\n");
                break;
            case BCFL_MEM:
                fprintf(file,"      <bcfl>mem</bcfl>\n");
                break;
            default:
                break;
        }

        fprintf(file,"      <temp>%4.3f K</temp>\n", pbeparm->temp);

        for (;icalc<=nosh->elec2calc[ielec];icalc++){ /* calc loop */

            /* Reinitialize per-calc pointers */
            mgparm = nosh->calc[icalc]->mgparm;
            pbeparm = nosh->calc[icalc]->pbeparm;

            fprintf(file,"      <calc>\n");
            fprintf(file,"          <id>%i</id>\n", (icalc+1));
            fprintf(file,"          <hx>%4.3f A</hx>\n", mgparm->grid[0]);
            fprintf(file,"          <hy>%4.3f A</hy>\n", mgparm->grid[1]);
            fprintf(file,"          <hz>%4.3f A</hz>\n", mgparm->grid[2]);
            fprintf(file,"          <xlen>%4.3f A</xlen>\n", mgparm->glen[0]);
            fprintf(file,"          <ylen>%4.3f A</ylen>\n", mgparm->glen[1]);
            fprintf(file,"          <zlen>%4.3f A</zlen>\n", mgparm->glen[2]);

            if (pbeparm->calcenergy == PCE_TOTAL) {
                fprintf(file,"          <totEnergy>%1.12E kJ/mol</totEnergy>\n",
                        (totEnergy[icalc]*conversion));
            } else if (pbeparm->calcenergy == PCE_COMPS) {
                            fprintf(file,"          <totEnergy>%1.12E kJ/mol</totEnergy>\n",
                        (totEnergy[icalc]*conversion));
                fprintf(file,"          <qfEnergy>%1.12E kJ/mol</qfEnergy>\n",
                        (0.5*qfEnergy[icalc]*conversion));
                fprintf(file,"          <qmEnergy>%1.12E kJ/mol</qmEnergy>\n",
                        (qmEnergy[icalc]*conversion));
                fprintf(file,"          <dielEnergy>%1.12E kJ/mol</dielEnergy>\n",
                        (dielEnergy[icalc]*conversion));
                for (i=0; i<nenergy[icalc]; i++){
                    fprintf(file,"          <atom>\n");
                    fprintf(file,"              <id>%i</id>\n", i+1);
                    fprintf(file,"              <energy>%1.12E kJ/mol</energy>\n",
                            (0.5*atomEnergy[icalc][i]*conversion));
                    fprintf(file,"          </atom>\n");
                }
            }


            if (pbeparm->calcforce == PCF_TOTAL) {
                    fprintf(file,"          <qfforce_x>%1.12E</qfforce_x>\n",
                    atomForce[icalc][0].qfForce[0]*conversion);
                fprintf(file,"          <qfforce_y>%1.12E</qfforce_y>\n",
                    atomForce[icalc][0].qfForce[1]*conversion);
                fprintf(file,"          <qfforce_z>%1.12E</qfforce_z>\n",
                    atomForce[icalc][0].qfForce[2]*conversion);
                fprintf(file,"          <ibforce_x>%1.12E</ibforce_x>\n",
                    atomForce[icalc][0].ibForce[0]*conversion);
                fprintf(file,"          <ibforce_y>%1.12E</ibforce_y>\n",
                    atomForce[icalc][0].ibForce[1]*conversion);
                fprintf(file,"          <ibforce_z>%1.12E</ibforce_z>\n",
                    atomForce[icalc][0].ibForce[2]*conversion);
                fprintf(file,"          <dbforce_x>%1.12E</dbforce_x>\n",
                    atomForce[icalc][0].dbForce[0]*conversion);
                fprintf(file,"          <dbforce_y>%1.12E</dbforce_y>\n",
                    atomForce[icalc][0].dbForce[1]*conversion);
                fprintf(file,"          <dbforce_z>%1.12E</dbforce_z>\n",
                    atomForce[icalc][0].dbForce[2]*conversion);
            }

            fprintf(file,"      </calc>\n");
        }

        fprintf(file,"    </elec>\n");
    }

/* Handle print energy statements */

for (i=0; i<nosh->nprint; i++) {

    if (nosh->printwhat[i] == NPT_ENERGY) {

        fprintf(file,"    <printEnergy>\n");
        fprintf(file,"        <equation>%d", nosh->printcalc[i][0]+1);

        for (j=1; j<nosh->printnarg[i]; j++) {
            if (nosh->printop[i][j-1] == 0) fprintf(file," +");
            else if (nosh->printop[i][j-1] == 1) fprintf(file, " -");
            fprintf(file, " %d", nosh->printcalc[i][j] +1);
        }

        fprintf(file, "</equation>\n");
        icalc = nosh->elec2calc[nosh->printcalc[i][0]];

        ltenergy = Vunit_kb * (1e-3) * Vunit_Na * \
            nosh->calc[icalc]->pbeparm->temp * totEnergy[icalc];

        for (j=1; j<nosh->printnarg[i]; j++) {
            icalc = nosh->elec2calc[nosh->printcalc[i][j]];
            /* Add or subtract? */
            if (nosh->printop[i][j-1] == 0) scalar = 1.0;
            else if (nosh->printop[i][j-1] == 1) scalar = -1.0;
            /* Accumulate */
            ltenergy += (scalar * Vunit_kb * (1e-3) * Vunit_Na *
                         nosh->calc[icalc]->pbeparm->temp * totEnergy[icalc]);
        }
        Vcom_reduce(com, &ltenergy, &gtenergy, 1, 2, 0);
        fprintf(file,"        <localEnergy>%1.12E kJ/mol</localEnergy>\n", \
                ltenergy);
        fprintf(file,"        <globalEnergy>%1.12E kJ/mol</globalEnergy>\n", \
                gtenergy);

        fprintf(file,"    </printEnergy>\n");
    }
}

/* Add ending tags and close the file */
fprintf(file,"</APBS>\n");
fclose(file);

return 1;
}

VPUBLIC int writedataMG(int rank,
                        NOsh *nosh,
                        PBEparm *pbeparm,
                        Vpmg *pmg
                       ) {

    char writestem[VMAX_ARGLEN];
    char outpath[VMAX_ARGLEN];
    char title[72];
    int i,
        nx,
        ny,
        nz,
        natoms;
    double hx,
           hy,
           hzed,
           xcent,
           ycent,
           zcent,
           xmin,
           ymin,
           zmin;

    Vgrid *grid;
    Vio *sock;

    if (nosh->bogus) return 1;

    for (i=0; i<pbeparm->numwrite; i++) {

        nx = pmg->pmgp->nx;
        ny = pmg->pmgp->ny;
        nz = pmg->pmgp->nz;
        hx = pmg->pmgp->hx;
        hy = pmg->pmgp->hy;
        hzed = pmg->pmgp->hzed;

        switch (pbeparm->writetype[i]) {

            case VDT_CHARGE:

                Vnm_tprint(1, "  Writing charge distribution to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                VASSERT(Vpmg_fillArray(pmg, pmg->rwork, VDT_CHARGE, 0.0,
                                       pbeparm->pbetype, pbeparm));
                sprintf(title, "CHARGE DISTRIBUTION (e)");
                break;

            case VDT_POT:

                Vnm_tprint(1, "  Writing potential to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                VASSERT(Vpmg_fillArray(pmg, pmg->rwork, VDT_POT, 0.0,
                                       pbeparm->pbetype, pbeparm));
                sprintf(title, "POTENTIAL (kT/e)");
                break;

            case VDT_SMOL:

                Vnm_tprint(1, "  Writing molecular accessibility to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                VASSERT(Vpmg_fillArray(pmg, pmg->rwork, VDT_SMOL,
                                       pbeparm->srad, pbeparm->pbetype, pbeparm));
                sprintf(title,
                        "SOLVENT ACCESSIBILITY -- MOLECULAR (%4.3f PROBE)",
                        pbeparm->srad);
                break;

            case VDT_SSPL:

                Vnm_tprint(1, "  Writing spline-based accessibility to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                VASSERT(Vpmg_fillArray(pmg, pmg->rwork, VDT_SSPL,
                                       pbeparm->swin, pbeparm->pbetype, pbeparm));
                sprintf(title,
                        "SOLVENT ACCESSIBILITY -- SPLINE (%4.3f WINDOW)",
                        pbeparm->swin);
                break;

            case VDT_VDW:

                Vnm_tprint(1, "  Writing van der Waals accessibility to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                VASSERT(Vpmg_fillArray(pmg, pmg->rwork, VDT_VDW, 0.0,
                                       pbeparm->pbetype, pbeparm));
                sprintf(title, "SOLVENT ACCESSIBILITY -- VAN DER WAALS");
                break;

            case VDT_IVDW:

                Vnm_tprint(1, "  Writing ion accessibility to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                VASSERT(Vpmg_fillArray(pmg, pmg->rwork, VDT_IVDW,
                                       pmg->pbe->maxIonRadius, pbeparm->pbetype, pbeparm));
                sprintf(title,
                        "ION ACCESSIBILITY -- SPLINE (%4.3f RADIUS)",
                        pmg->pbe->maxIonRadius);
                break;

            case VDT_LAP:

                Vnm_tprint(1, "  Writing potential Laplacian to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                VASSERT(Vpmg_fillArray(pmg, pmg->rwork, VDT_LAP, 0.0,
                                       pbeparm->pbetype, pbeparm));
                sprintf(title,
                        "POTENTIAL LAPLACIAN (kT/e/A^2)");
                break;

            case VDT_EDENS:

                Vnm_tprint(1, "  Writing energy density to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                VASSERT(Vpmg_fillArray(pmg, pmg->rwork, VDT_EDENS, 0.0,
                                       pbeparm->pbetype, pbeparm));
                sprintf(title, "ENERGY DENSITY (kT/e/A)^2");
                break;

            case VDT_NDENS:

                Vnm_tprint(1, "  Writing number density to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                VASSERT(Vpmg_fillArray(pmg, pmg->rwork, VDT_NDENS, 0.0,
                                       pbeparm->pbetype, pbeparm));
                sprintf(title,
                        "ION NUMBER DENSITY (M)");
                break;

            case VDT_QDENS:

                Vnm_tprint(1, "  Writing charge density to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                VASSERT(Vpmg_fillArray(pmg, pmg->rwork, VDT_QDENS, 0.0,
                                       pbeparm->pbetype, pbeparm));
                sprintf(title,
                        "ION CHARGE DENSITY (e_c * M)");
                break;

            case VDT_DIELX:

                Vnm_tprint(1, "  Writing x-shifted dielectric map to ");
                xcent = pmg->pmgp->xcent + 0.5*hx;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                VASSERT(Vpmg_fillArray(pmg, pmg->rwork, VDT_DIELX, 0.0,
                                       pbeparm->pbetype, pbeparm));
                sprintf(title,
                        "X-SHIFTED DIELECTRIC MAP");
                break;

            case VDT_DIELY:

                Vnm_tprint(1, "  Writing y-shifted dielectric map to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent + 0.5*hy;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                VASSERT(Vpmg_fillArray(pmg, pmg->rwork, VDT_DIELY, 0.0,
                                       pbeparm->pbetype, pbeparm));
                sprintf(title,
                        "Y-SHIFTED DIELECTRIC MAP");
                break;

            case VDT_DIELZ:

                Vnm_tprint(1, "  Writing z-shifted dielectric map to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent + 0.5*hzed;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                VASSERT(Vpmg_fillArray(pmg, pmg->rwork, VDT_DIELZ, 0.0,
                                       pbeparm->pbetype, pbeparm));
                sprintf(title,
                        "Z-SHIFTED DIELECTRIC MAP");
                break;

            case VDT_KAPPA:

                Vnm_tprint(1, "  Writing kappa map to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                VASSERT(Vpmg_fillArray(pmg, pmg->rwork, VDT_KAPPA, 0.0,
                                       pbeparm->pbetype, pbeparm));
                sprintf(title,
                        "KAPPA MAP");
                break;

            case VDT_ATOMPOT:

                Vnm_tprint(1, "  Writing atom potentials to ");
                xcent = pmg->pmgp->xcent;
                ycent = pmg->pmgp->ycent;
                zcent = pmg->pmgp->zcent;
                xmin = xcent - 0.5*(nx-1)*hx;
                ymin = ycent - 0.5*(ny-1)*hy;
                zmin = zcent - 0.5*(nz-1)*hzed;
                VASSERT(Vpmg_fillArray(pmg, pmg->rwork, VDT_ATOMPOT, 0.0,
                                       pbeparm->pbetype, pbeparm));
                sprintf(title,
                        "ATOM POTENTIALS");
                break;
            default:

                Vnm_tprint(2, "Invalid data type for writing!\n");
                return 0;
        }


#ifdef HAVE_MPI_H
        sprintf(writestem, "%s-PE%d", pbeparm->writestem[i], rank);
#else
        if(nosh->ispara){
            sprintf(writestem, "%s-PE%d", pbeparm->writestem[i],nosh->proc_rank);
        }else{
            sprintf(writestem, "%s", pbeparm->writestem[i]);
        }
#endif

        switch (pbeparm->writefmt[i]) {

            case VDF_DX:
                sprintf(outpath, "%s.%s", writestem, "dx");
                Vnm_tprint(1, "%s\n", outpath);
                grid = Vgrid_ctor(nx, ny, nz, hx, hy, hzed, xmin, ymin, zmin,
                                  pmg->rwork);
                Vgrid_writeDX(grid, "FILE", "ASC", VNULL, outpath, title,
                              pmg->pvec);
                Vgrid_dtor(&grid);
                break;

            case VDF_AVS:
                sprintf(outpath, "%s.%s", writestem, "ucd");
                Vnm_tprint(1, "%s\n", outpath);
                Vnm_tprint(2, "Sorry, AVS format isn't supported for \
uniform meshes yet!\n");
                break;

            case VDF_MCSF:
                sprintf(outpath, "%s.%s", writestem, "mcsf");
                Vnm_tprint(1, "%s\n", outpath);
                Vnm_tprint(2, "Sorry, MCSF format isn't supported for \
                           uniform meshes yet!\n");
                break;

            case VDF_UHBD:
                sprintf(outpath, "%s.%s", writestem, "grd");
                Vnm_tprint(1, "%s\n", outpath);
                grid = Vgrid_ctor(nx, ny, nz, hx, hy, hzed, xmin, ymin, zmin,
                                  pmg->rwork);
                Vgrid_writeUHBD(grid, "FILE", "ASC", VNULL, outpath, title,
                                pmg->pvec);
                Vgrid_dtor(&grid);
                break;

            case VDF_GZ:
                sprintf(outpath, "%s.%s", writestem, "dx.gz");
                Vnm_tprint(1, "%s\n", outpath);
                grid = Vgrid_ctor(nx, ny, nz, hx, hy, hzed, xmin, ymin, zmin,
                                  pmg->rwork);
                Vgrid_writeGZ(grid, "FILE", "ASC", VNULL, outpath, title,
                              pmg->pvec);
                Vgrid_dtor(&grid);
                break;
            case VDF_FLAT:
                sprintf(outpath, "%s.%s", writestem, "txt");
                Vnm_tprint(1, "%s\n", outpath);
                Vnm_print(0, "routines:  Opening virtual socket...\n");
                sock = Vio_ctor("FILE","ASC",VNULL,outpath,"w");
                if (sock == VNULL) {
                    Vnm_print(2, "routines:  Problem opening virtual socket %s\n",
                              outpath);
                    return 0;
                }
                if (Vio_connect(sock, 0) < 0) {
                    Vnm_print(2, "routines: Problem connecting virtual socket %s\n",
                              outpath);
                    return 0;
                }
                Vio_printf(sock, "# Data from %s\n", PACKAGE_STRING);
                Vio_printf(sock, "# \n");
                Vio_printf(sock, "# %s\n", title);
                Vio_printf(sock, "# \n");
                natoms = pmg->pbe->alist[pbeparm->molid-1].number;
                for(i=0;i<natoms;i++)
                    Vio_printf(sock, "%12.6e\n", pmg->rwork[i]);
                break;
            default:
                Vnm_tprint(2, "Bogus data format (%d)!\n",
                           pbeparm->writefmt[i]);
                break;
        }

    }

    return 1;
}

VPUBLIC double returnEnergy(Vcom *com,
                            NOsh *nosh,
                            double totEnergy[NOSH_MAXCALC],
                            int iprint
                           ){

    int iarg,
        calcid;
    double ltenergy,
           scalar;

    calcid = nosh->elec2calc[nosh->printcalc[iprint][0]];
    if (nosh->calc[calcid]->pbeparm->calcenergy != PCE_NO) {
        ltenergy = Vunit_kb * (1e-3) * Vunit_Na *
        nosh->calc[calcid]->pbeparm->temp * totEnergy[calcid];
    } else {
        Vnm_tprint( 2, " No energy available in Calculation %d\n", calcid+1);
        return 0.0;
    }
    for (iarg=1; iarg<nosh->printnarg[iprint]; iarg++){
        calcid = nosh->elec2calc[nosh->printcalc[iprint][iarg]];
        /* Add or substract */
        if (nosh->printop[iprint][iarg-1] == 0) scalar = 1.0;
        else if (nosh->printop[iprint][iarg-1] == 1) scalar = -1.0;
        /* Accumulate */
        ltenergy += (scalar * Vunit_kb * (1e-3) * Vunit_Na *
                     nosh->calc[calcid]->pbeparm->temp * totEnergy[calcid]);
    }

    return ltenergy;
}

VPUBLIC int printEnergy(Vcom *com,
                        NOsh *nosh,
                        double totEnergy[NOSH_MAXCALC],
                        int iprint
                       ) {

    int iarg,
        calcid;
    double ltenergy,
           gtenergy,
           scalar;

    Vnm_tprint( 2, "Warning: The 'energy' print keyword is deprecated.\n" \
                   "         Use elecEnergy for electrostatics energy calcs.\n\n");

    if (Vstring_strcasecmp(nosh->elecname[nosh->printcalc[iprint][0]], "") == 0){
        Vnm_tprint( 1, "print energy %d ", nosh->printcalc[iprint][0]+1);
    } else {
        Vnm_tprint( 1, "print energy %d (%s) ", nosh->printcalc[iprint][0]+1,
                    nosh->elecname[nosh->printcalc[iprint][0]]);
    }
    for (iarg=1; iarg<nosh->printnarg[iprint]; iarg++) {
        if (nosh->printop[iprint][iarg-1] == 0)
            Vnm_tprint(1, "+ ");
        else if (nosh->printop[iprint][iarg-1] == 1)
            Vnm_tprint(1, "- ");
        else {
            Vnm_tprint( 2, "Undefined PRINT operation!\n");
            return 0;
        }
        if (Vstring_strcasecmp(nosh->elecname[nosh->printcalc[iprint][iarg]],
                               "") == 0) {
            Vnm_tprint( 1, "%d ", nosh->printcalc[iprint][iarg]+1);
        } else {
            Vnm_tprint( 1, "%d (%s) ", nosh->printcalc[iprint][iarg]+1,
                        nosh->elecname[nosh->printcalc[iprint][iarg]]);
        }
    }
    Vnm_tprint(1, "end\n");
    calcid = nosh->elec2calc[nosh->printcalc[iprint][0]];
    if (nosh->calc[calcid]->pbeparm->calcenergy != PCE_NO) {
        ltenergy = Vunit_kb * (1e-3) * Vunit_Na *
        nosh->calc[calcid]->pbeparm->temp * totEnergy[calcid];
    } else {
        Vnm_tprint( 2, "  Didn't calculate energy in Calculation \
#%d\n", calcid+1);
        return 0;
    }
    for (iarg=1; iarg<nosh->printnarg[iprint]; iarg++) {
        calcid = nosh->elec2calc[nosh->printcalc[iprint][iarg]];
        /* Add or subtract? */
        if (nosh->printop[iprint][iarg-1] == 0) scalar = 1.0;
        else if (nosh->printop[iprint][iarg-1] == 1) scalar = -1.0;
        /* Accumulate */
        ltenergy += (scalar * Vunit_kb * (1e-3) * Vunit_Na *
                     nosh->calc[calcid]->pbeparm->temp * totEnergy[calcid]);
    }

    Vnm_tprint( 1, "  Local net energy (PE %d) = %1.12E kJ/mol\n",
                Vcom_rank(com), ltenergy);
    Vnm_tprint( 0, "printEnergy:  Performing global reduction (sum)\n");
    Vcom_reduce(com, &ltenergy, &gtenergy, 1, 2, 0);
    Vnm_tprint( 1, "  Global net ELEC energy = %1.12E kJ/mol\n", gtenergy);

    return 1;

}

VPUBLIC int printElecEnergy(Vcom *com,
                            NOsh *nosh,
                            double totEnergy[NOSH_MAXCALC],
                            int iprint
                           ) {

    int iarg,
        calcid;
    double ltenergy,
           gtenergy,
           scalar;

    if (Vstring_strcasecmp(nosh->elecname[nosh->printcalc[iprint][0]], "") == 0){
        Vnm_tprint( 1, "\nprint energy %d ", nosh->printcalc[iprint][0]+1);
    } else {
        Vnm_tprint( 1, "\nprint energy %d (%s) ", nosh->printcalc[iprint][0]+1,
                    nosh->elecname[nosh->printcalc[iprint][0]]);
    }
    for (iarg=1; iarg<nosh->printnarg[iprint]; iarg++) {
        if (nosh->printop[iprint][iarg-1] == 0)
            Vnm_tprint(1, "+ ");
        else if (nosh->printop[iprint][iarg-1] == 1)
            Vnm_tprint(1, "- ");
        else {
            Vnm_tprint( 2, "Undefined PRINT operation!\n");
            return 0;
        }
        if (Vstring_strcasecmp(nosh->elecname[nosh->printcalc[iprint][iarg]],
                               "") == 0) {
            Vnm_tprint( 1, "%d ", nosh->printcalc[iprint][iarg]+1);
        } else {
            Vnm_tprint( 1, "%d (%s) ", nosh->printcalc[iprint][iarg]+1,
                        nosh->elecname[nosh->printcalc[iprint][iarg]]);
        }
    }
    Vnm_tprint(1, "end\n");
    calcid = nosh->elec2calc[nosh->printcalc[iprint][0]];
    if (nosh->calc[calcid]->pbeparm->calcenergy != PCE_NO) {
        ltenergy = Vunit_kb * (1e-3) * Vunit_Na *
        nosh->calc[calcid]->pbeparm->temp * totEnergy[calcid];
    } else {
        Vnm_tprint( 2, "  Didn't calculate energy in Calculation \
#%d\n", calcid+1);
        return 0;
    }
    for (iarg=1; iarg<nosh->printnarg[iprint]; iarg++) {
        calcid = nosh->elec2calc[nosh->printcalc[iprint][iarg]];
        /* Add or subtract? */
        if (nosh->printop[iprint][iarg-1] == 0) scalar = 1.0;
        else if (nosh->printop[iprint][iarg-1] == 1) scalar = -1.0;
        /* Accumulate */
        ltenergy += (scalar * Vunit_kb * (1e-3) * Vunit_Na *
                     nosh->calc[calcid]->pbeparm->temp * totEnergy[calcid]);
    }

    Vnm_tprint( 1, "  Local net energy (PE %d) = %1.12E kJ/mol\n",
                Vcom_rank(com), ltenergy);
    Vnm_tprint( 0, "printEnergy:  Performing global reduction (sum)\n");
    Vcom_reduce(com, &ltenergy, &gtenergy, 1, 2, 0);
    Vnm_tprint( 1, "  Global net ELEC energy = %1.12E kJ/mol\n", gtenergy);

    return 1;

}

VPUBLIC int printApolEnergy(NOsh *nosh,
                            int iprint
                           ) {

    int iarg,
        calcid;
    double gtenergy,
           scalar;
    APOLparm *apolparm = VNULL;

    if (Vstring_strcasecmp(nosh->apolname[nosh->printcalc[iprint][0]], "") == 0){
        Vnm_tprint( 1, "\nprint APOL energy %d ", nosh->printcalc[iprint][0]+1);
    } else {
        Vnm_tprint( 1, "\nprint APOL energy %d (%s) ", nosh->printcalc[iprint][0]+1,
                    nosh->apolname[nosh->printcalc[iprint][0]]);
    }
    for (iarg=1; iarg<nosh->printnarg[iprint]; iarg++) {
        if (nosh->printop[iprint][iarg-1] == 0)
            Vnm_tprint(1, "+ ");
        else if (nosh->printop[iprint][iarg-1] == 1)
            Vnm_tprint(1, "- ");
        else {
            Vnm_tprint( 2, "Undefined PRINT operation!\n");
            return 0;
        }
        if (Vstring_strcasecmp(nosh->apolname[nosh->printcalc[iprint][iarg]],
                               "") == 0) {
            Vnm_tprint( 1, "%d ", nosh->printcalc[iprint][iarg]+1);
        } else {
            Vnm_tprint( 1, "%d (%s) ", nosh->printcalc[iprint][iarg]+1,
                        nosh->apolname[nosh->printcalc[iprint][iarg]]);
        }
    }
    Vnm_tprint(1, "end\n");

    calcid = nosh->apol2calc[nosh->printcalc[iprint][0]];
    apolparm = nosh->calc[calcid]->apolparm;

    if (apolparm->calcenergy == ACE_TOTAL) {
        gtenergy = ((apolparm->gamma*apolparm->sasa) + (apolparm->press*apolparm->sav) + (apolparm->wcaEnergy));
    } else {
        Vnm_tprint( 2, "  Didn't calculate energy in Calculation #%d\n", calcid+1);
        return 0;
    }
    for (iarg=1; iarg<nosh->printnarg[iprint]; iarg++) {
        calcid = nosh->apol2calc[nosh->printcalc[iprint][iarg]];
        apolparm = nosh->calc[calcid]->apolparm;

        /* Add or subtract? */
        if (nosh->printop[iprint][iarg-1] == 0) scalar = 1.0;
        else if (nosh->printop[iprint][iarg-1] == 1) scalar = -1.0;
        /* Accumulate */
        gtenergy += (scalar * ((apolparm->gamma*apolparm->sasa) +
                               (apolparm->press*apolparm->sav) +
                               (apolparm->wcaEnergy)));

    }

    Vnm_tprint( 1, "  Global net APOL energy = %1.12E kJ/mol\n", gtenergy);
    return 1;
}

VPUBLIC int printForce(Vcom *com,
                       NOsh *nosh,
                       int nforce[NOSH_MAXCALC],
                       AtomForce *atomForce[NOSH_MAXCALC],
                       int iprint
                      ) {

    int iarg,
        ifr,
        ivc,
        calcid,
        refnforce,
        refcalcforce;
    double temp,
           scalar,
           totforce[3];
    AtomForce *lforce,
              *gforce,
              *aforce;

    Vnm_tprint( 2, "Warning: The 'force' print keyword is deprecated.\n" \
                   "         Use elecForce for electrostatics force calcs.\n\n");

    if (Vstring_strcasecmp(nosh->elecname[nosh->printcalc[iprint][0]], "") == 0){
        Vnm_tprint( 1, "print force %d ", nosh->printcalc[iprint][0]+1);
    } else {
        Vnm_tprint( 1, "print force %d (%s) ", nosh->printcalc[iprint][0]+1,
                    nosh->elecname[nosh->printcalc[iprint][0]]);
    }
    for (iarg=1; iarg<nosh->printnarg[iprint]; iarg++) {
        if (nosh->printop[iprint][iarg-1] == 0)
            Vnm_tprint(1, "+ ");
        else if (nosh->printop[iprint][iarg-1] == 1)
            Vnm_tprint(1, "- ");
        else {
            Vnm_tprint( 2, "Undefined PRINT operation!\n");
            return 0;
        }
        if (Vstring_strcasecmp(nosh->elecname[nosh->printcalc[iprint][iarg]],
                               "") == 0) {
            Vnm_tprint( 1, "%d ", nosh->printcalc[iprint][iarg]+1);
        } else {
            Vnm_tprint( 1, "%d (%s) ", nosh->printcalc[iprint][iarg]+1,
                        nosh->elecname[nosh->printcalc[iprint][iarg]]);
        }
    }
    Vnm_tprint(1, "end\n");

    /* First, go through and make sure we did the same type of force
        * evaluation in each of the requested calculations */
    calcid = nosh->elec2calc[nosh->printcalc[iprint][0]];
    refnforce = nforce[calcid];
    refcalcforce = nosh->calc[calcid]->pbeparm->calcforce;
    if (refcalcforce == PCF_NO) {
        Vnm_tprint( 2, "  Didn't calculate force in calculation \
#%d\n", calcid+1);
        return 0;
    }
    for (iarg=1; iarg<nosh->printnarg[iprint]; iarg++) {
        calcid = nosh->elec2calc[nosh->printcalc[iprint][iarg]-1];
        if (nosh->calc[calcid]->pbeparm->calcforce != refcalcforce) {
            Vnm_tprint(2, "  Inconsistent calcforce declarations in \
calculations %d and %d\n", nosh->elec2calc[nosh->printcalc[iprint][0]]+1,
                       calcid+1);
            return 0;
        }
        if (nforce[calcid] != refnforce) {
            Vnm_tprint(2, "  Inconsistent number of forces evaluated in \
calculations %d and %d\n", nosh->elec2calc[nosh->printcalc[iprint][0]]+1,
                       calcid+1);
            return 0;
        }
    }

    /* Now, allocate space to accumulate the forces */
    lforce = (AtomForce *)Vmem_malloc(VNULL, refnforce, sizeof(AtomForce));
    gforce = (AtomForce *)Vmem_malloc(VNULL, refnforce, sizeof(AtomForce));

    /* Now, accumulate the forces */
    calcid = nosh->elec2calc[nosh->printcalc[iprint][0]];
    aforce = atomForce[calcid];
    temp = nosh->calc[calcid]->pbeparm->temp;

    /* Load up the first calculation */
    if (refcalcforce == PCF_TOTAL) {
        /* Set to total force */
        for (ivc=0; ivc<3; ivc++) {
            lforce[0].qfForce[ivc] =
            Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[0].qfForce[ivc];
            lforce[0].ibForce[ivc] =
                Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[0].ibForce[ivc];
            lforce[0].dbForce[ivc] =
                Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[0].dbForce[ivc];
        }
    } else if (refcalcforce == PCF_COMPS) {
        for (ifr=0; ifr<refnforce; ifr++) {
            for (ivc=0; ivc<3; ivc++) {
                lforce[ifr].qfForce[ivc] =
                Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[ifr].qfForce[ivc];
                lforce[ifr].ibForce[ivc] =
                    Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[ifr].ibForce[ivc];
                lforce[ifr].dbForce[ivc] =
                    Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[ifr].dbForce[ivc];
            }
        }
    }

    /* Load up the rest of the calculations */
    for (iarg=1; iarg<nosh->printnarg[iprint]; iarg++) {
        calcid = nosh->elec2calc[nosh->printcalc[iprint][iarg]];
        temp = nosh->calc[calcid]->pbeparm->temp;
        aforce = atomForce[calcid];
        /* Get operation */
        if (nosh->printop[iprint][iarg-1] == 0) scalar = +1.0;
        else if (nosh->printop[iprint][iarg-1] == 1) scalar = -1.0;
        else scalar = 0.0;
        /* Accumulate */
        if (refcalcforce == PCF_TOTAL) {
            /* Set to total force */
            for (ivc=0; ivc<3; ivc++) {
                lforce[0].qfForce[ivc] +=
                (scalar*Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[0].qfForce[ivc]);
                lforce[0].ibForce[ivc] +=
                    (scalar*Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[0].ibForce[ivc]);
                lforce[0].dbForce[ivc] +=
                    (scalar*Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[0].dbForce[ivc]);
            }
        } else if (refcalcforce == PCF_COMPS) {
            for (ifr=0; ifr<refnforce; ifr++) {
                for (ivc=0; ivc<3; ivc++) {
                    lforce[ifr].qfForce[ivc] +=
                    (scalar*Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[ifr].qfForce[ivc]);
                    lforce[ifr].ibForce[ivc] +=
                        (scalar*Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[ifr].ibForce[ivc]);
                    lforce[ifr].dbForce[ivc] +=
                        (scalar*Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[ifr].dbForce[ivc]);
                }
            }
        }
    }

    Vnm_tprint( 0, "printEnergy:  Performing global reduction (sum)\n");
    for (ifr=0; ifr<refnforce; ifr++) {
        Vcom_reduce(com, lforce[ifr].qfForce, gforce[ifr].qfForce, 3, 2, 0);
        Vcom_reduce(com, lforce[ifr].ibForce, gforce[ifr].ibForce, 3, 2, 0);
        Vcom_reduce(com, lforce[ifr].dbForce, gforce[ifr].dbForce, 3, 2, 0);
    }

#if 0
    if (refcalcforce == PCF_TOTAL) {
        Vnm_tprint( 1, "  Local net fixed charge force = \
(%1.12E, %1.12E, %1.12E) kJ/mol/A\n", lforce[0].qfForce[0],
                    lforce[0].qfForce[1], lforce[0].qfForce[2]);
        Vnm_tprint( 1, "  Local net ionic boundary force = \
(%1.12E, %1.12E, %1.12E) kJ/mol/A\n", lforce[0].ibForce[0],
                    lforce[0].ibForce[1], lforce[0].ibForce[2]);
        Vnm_tprint( 1, "  Local net dielectric boundary force = \
(%1.12E, %1.12E, %1.12E) kJ/mol/A\n", lforce[0].dbForce[0],
                    lforce[0].dbForce[1], lforce[0].dbForce[2]);
    } else if (refcalcforce == PCF_COMPS) {
        for (ifr=0; ifr<refnforce; ifr++) {
            Vnm_tprint( 1, "  Local fixed charge force \
(atom %d) = (%1.12E, %1.12E, %1.12E) kJ/mol/A\n", ifr, lforce[ifr].qfForce[0],
                        lforce[ifr].qfForce[1], lforce[ifr].qfForce[2]);
            Vnm_tprint( 1, "  Local ionic boundary force \
(atom %d) = (%1.12E, %1.12E, %1.12E) kJ/mol/A\n", ifr, lforce[ifr].ibForce[0],
                        lforce[ifr].ibForce[1], lforce[ifr].ibForce[2]);
            Vnm_tprint( 1, "  Local dielectric boundary force \
(atom %d) = (%1.12E, %1.12E, %1.12E) kJ/mol/A\n", ifr, lforce[ifr].dbForce[0],
                        lforce[ifr].dbForce[1], lforce[ifr].dbForce[2]);
        }
    }
#endif

    if (refcalcforce == PCF_TOTAL) {
        Vnm_tprint( 1, "  Printing net forces (kJ/mol/A).\n");
        Vnm_tprint( 1, "  Legend:\n");
        Vnm_tprint( 1, "    tot -- Total force\n");
        Vnm_tprint( 1, "    qf  -- Fixed charge force\n");
        Vnm_tprint( 1, "    db  -- Dielectric boundary force\n");
        Vnm_tprint( 1, "    ib  -- Ionic boundary force\n");

        for (ivc=0; ivc<3; ivc++) {
            totforce[ivc] =
            gforce[0].qfForce[ivc] + gforce[0].ibForce[ivc] \
            + gforce[0].dbForce[ivc];
        }

        Vnm_tprint( 1, "  tot %1.12E  %1.12E  %1.12E\n", totforce[0],
                    totforce[1], totforce[2]);
        Vnm_tprint( 1, "  qf  %1.12E  %1.12E  %1.12E\n", gforce[0].qfForce[0],
                    gforce[0].qfForce[1], gforce[0].qfForce[2]);
        Vnm_tprint( 1, "  ib  %1.12E  %1.12E  %1.12E\n", gforce[0].ibForce[0],
                    gforce[0].ibForce[1], gforce[0].ibForce[2]);
        Vnm_tprint( 1, "  db  %1.12E  %1.12E  %1.12E\n", gforce[0].dbForce[0],
                    gforce[0].dbForce[1], gforce[0].dbForce[2]);

    } else if (refcalcforce == PCF_COMPS) {

        Vnm_tprint( 1, "  Printing per-atom forces (kJ/mol/A).\n");
        Vnm_tprint( 1, "  Legend:\n");
        Vnm_tprint( 1, "    tot n -- Total force for atom n\n");
        Vnm_tprint( 1, "    qf  n -- Fixed charge force for atom n\n");
        Vnm_tprint( 1, "    db  n -- Dielectric boundary force for atom n\n");
        Vnm_tprint( 1, "    ib  n -- Ionic boundary force for atom n\n");
        Vnm_tprint( 1, "    tot all -- Total force for system\n");

        totforce[0] = 0.0;
        totforce[1] = 0.0;
        totforce[2] = 0.0;

        for (ifr=0; ifr<refnforce; ifr++) {
            Vnm_tprint( 1, "  qf  %d  %1.12E  %1.12E  %1.12E\n", ifr,
                        gforce[ifr].qfForce[0], gforce[ifr].qfForce[1],
                        gforce[ifr].qfForce[2]);
            Vnm_tprint( 1, "  ib  %d  %1.12E  %1.12E  %1.12E\n", ifr,
                        gforce[ifr].ibForce[0], gforce[ifr].ibForce[1],
                        gforce[ifr].ibForce[2]);
            Vnm_tprint( 1, "  db  %d  %1.12E  %1.12E  %1.12E\n", ifr,
                        gforce[ifr].dbForce[0], gforce[ifr].dbForce[1],
                        gforce[ifr].dbForce[2]);
            Vnm_tprint( 1, "  tot %d  %1.12E  %1.12E  %1.12E\n", ifr,
                        (gforce[ifr].dbForce[0] \
                         + gforce[ifr].ibForce[0] +
                         gforce[ifr].qfForce[0]),
                        (gforce[ifr].dbForce[1] \
                         + gforce[ifr].ibForce[1] +
                         gforce[ifr].qfForce[1]),
                        (gforce[ifr].dbForce[2] \
                         + gforce[ifr].ibForce[2] +
                         gforce[ifr].qfForce[2]));
            for (ivc=0; ivc<3; ivc++) {
                totforce[ivc] += (gforce[ifr].dbForce[ivc] \
                                + gforce[ifr].ibForce[ivc] \
                                  + gforce[ifr].qfForce[ivc]);
            }
        }
        Vnm_tprint( 1, "  tot all %1.12E  %1.12E  %1.12E\n", totforce[0],
                    totforce[1], totforce[2]);
    }

    Vmem_free(VNULL, refnforce, sizeof(AtomForce), (void **)(&lforce));
    Vmem_free(VNULL, refnforce, sizeof(AtomForce), (void **)(&gforce));

    return 1;

}

VPUBLIC int printElecForce(Vcom *com,
                           NOsh *nosh,
                           int nforce[NOSH_MAXCALC],
                           AtomForce *atomForce[NOSH_MAXCALC],
                           int iprint
                          ) {

    int iarg,
        ifr,
        ivc,
        calcid,
        refnforce,
        refcalcforce;
    double temp,
           scalar,
           totforce[3];
    AtomForce *lforce, *gforce, *aforce;

    if (Vstring_strcasecmp(nosh->elecname[nosh->printcalc[iprint][0]], "") == 0){
        Vnm_tprint( 1, "print force %d ", nosh->printcalc[iprint][0]+1);
    } else {
        Vnm_tprint( 1, "print force %d (%s) ", nosh->printcalc[iprint][0]+1,
                    nosh->elecname[nosh->printcalc[iprint][0]]);
    }
    for (iarg=1; iarg<nosh->printnarg[iprint]; iarg++) {
        if (nosh->printop[iprint][iarg-1] == 0)
            Vnm_tprint(1, "+ ");
        else if (nosh->printop[iprint][iarg-1] == 1)
            Vnm_tprint(1, "- ");
        else {
            Vnm_tprint( 2, "Undefined PRINT operation!\n");
            return 0;
        }
        if (Vstring_strcasecmp(nosh->elecname[nosh->printcalc[iprint][iarg]],
                               "") == 0) {
            Vnm_tprint( 1, "%d ", nosh->printcalc[iprint][iarg]+1);
        } else {
            Vnm_tprint( 1, "%d (%s) ", nosh->printcalc[iprint][iarg]+1,
                        nosh->elecname[nosh->printcalc[iprint][iarg]]);
        }
    }
    Vnm_tprint(1, "end\n");

    /* First, go through and make sure we did the same type of force
        * evaluation in each of the requested calculations */
    calcid = nosh->elec2calc[nosh->printcalc[iprint][0]];
    refnforce = nforce[calcid];
    refcalcforce = nosh->calc[calcid]->pbeparm->calcforce;
    if (refcalcforce == PCF_NO) {
        Vnm_tprint( 2, "  Didn't calculate force in calculation \
#%d\n", calcid+1);
        return 0;
    }
    for (iarg=1; iarg<nosh->printnarg[iprint]; iarg++) {
        calcid = nosh->elec2calc[nosh->printcalc[iprint][iarg]-1];
        if (nosh->calc[calcid]->pbeparm->calcforce != refcalcforce) {
            Vnm_tprint(2, "  Inconsistent calcforce declarations in \
calculations %d and %d\n", nosh->elec2calc[nosh->printcalc[iprint][0]]+1,
                       calcid+1);
            return 0;
        }
        if (nforce[calcid] != refnforce) {
            Vnm_tprint(2, "  Inconsistent number of forces evaluated in \
calculations %d and %d\n", nosh->elec2calc[nosh->printcalc[iprint][0]]+1,
                       calcid+1);
            return 0;
        }
    }

    /* Now, allocate space to accumulate the forces */
    lforce = (AtomForce *)Vmem_malloc(VNULL, refnforce, sizeof(AtomForce));
    gforce = (AtomForce *)Vmem_malloc(VNULL, refnforce, sizeof(AtomForce));

    /* Now, accumulate the forces */
    calcid = nosh->elec2calc[nosh->printcalc[iprint][0]];
    aforce = atomForce[calcid];
    temp = nosh->calc[calcid]->pbeparm->temp;

    /* Load up the first calculation */
    if (refcalcforce == PCF_TOTAL) {
        /* Set to total force */
        for (ivc=0; ivc<3; ivc++) {
            lforce[0].qfForce[ivc] =
            Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[0].qfForce[ivc];
            lforce[0].ibForce[ivc] =
                Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[0].ibForce[ivc];
            lforce[0].dbForce[ivc] =
                Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[0].dbForce[ivc];
        }
    } else if (refcalcforce == PCF_COMPS) {
        for (ifr=0; ifr<refnforce; ifr++) {
            for (ivc=0; ivc<3; ivc++) {
                lforce[ifr].qfForce[ivc] =
                Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[ifr].qfForce[ivc];
                lforce[ifr].ibForce[ivc] =
                    Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[ifr].ibForce[ivc];
                lforce[ifr].dbForce[ivc] =
                    Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[ifr].dbForce[ivc];
            }
        }
    }

    /* Load up the rest of the calculations */
    for (iarg=1; iarg<nosh->printnarg[iprint]; iarg++) {
        calcid = nosh->elec2calc[nosh->printcalc[iprint][iarg]];
        temp = nosh->calc[calcid]->pbeparm->temp;
        aforce = atomForce[calcid];
        /* Get operation */
        if (nosh->printop[iprint][iarg-1] == 0) scalar = +1.0;
        else if (nosh->printop[iprint][iarg-1] == 1) scalar = -1.0;
        else scalar = 0.0;
        /* Accumulate */
        if (refcalcforce == PCF_TOTAL) {
            /* Set to total force */
            for (ivc=0; ivc<3; ivc++) {
                lforce[0].qfForce[ivc] +=
                (scalar*Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[0].qfForce[ivc]);
                lforce[0].ibForce[ivc] +=
                    (scalar*Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[0].ibForce[ivc]);
                lforce[0].dbForce[ivc] +=
                    (scalar*Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[0].dbForce[ivc]);
            }
        } else if (refcalcforce == PCF_COMPS) {
            for (ifr=0; ifr<refnforce; ifr++) {
                for (ivc=0; ivc<3; ivc++) {
                    lforce[ifr].qfForce[ivc] +=
                    (scalar*Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[ifr].qfForce[ivc]);
                    lforce[ifr].ibForce[ivc] +=
                        (scalar*Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[ifr].ibForce[ivc]);
                    lforce[ifr].dbForce[ivc] +=
                        (scalar*Vunit_kb*(1e-3)*Vunit_Na*temp*aforce[ifr].dbForce[ivc]);
                }
            }
        }
    }

    Vnm_tprint( 0, "printEnergy:  Performing global reduction (sum)\n");
    for (ifr=0; ifr<refnforce; ifr++) {
        Vcom_reduce(com, lforce[ifr].qfForce, gforce[ifr].qfForce, 3, 2, 0);
        Vcom_reduce(com, lforce[ifr].ibForce, gforce[ifr].ibForce, 3, 2, 0);
        Vcom_reduce(com, lforce[ifr].dbForce, gforce[ifr].dbForce, 3, 2, 0);
    }

#if 0
    if (refcalcforce == PCF_TOTAL) {
        Vnm_tprint( 1, "  Local net fixed charge force = \
(%1.12E, %1.12E, %1.12E) kJ/mol/A\n", lforce[0].qfForce[0],
                    lforce[0].qfForce[1], lforce[0].qfForce[2]);
        Vnm_tprint( 1, "  Local net ionic boundary force = \
(%1.12E, %1.12E, %1.12E) kJ/mol/A\n", lforce[0].ibForce[0],
                    lforce[0].ibForce[1], lforce[0].ibForce[2]);
        Vnm_tprint( 1, "  Local net dielectric boundary force = \
(%1.12E, %1.12E, %1.12E) kJ/mol/A\n", lforce[0].dbForce[0],
                    lforce[0].dbForce[1], lforce[0].dbForce[2]);
    } else if (refcalcforce == PCF_COMPS) {
        for (ifr=0; ifr<refnforce; ifr++) {
            Vnm_tprint( 1, "  Local fixed charge force \
(atom %d) = (%1.12E, %1.12E, %1.12E) kJ/mol/A\n", ifr, lforce[ifr].qfForce[0],
                        lforce[ifr].qfForce[1], lforce[ifr].qfForce[2]);
            Vnm_tprint( 1, "  Local ionic boundary force \
(atom %d) = (%1.12E, %1.12E, %1.12E) kJ/mol/A\n", ifr, lforce[ifr].ibForce[0],
                        lforce[ifr].ibForce[1], lforce[ifr].ibForce[2]);
            Vnm_tprint( 1, "  Local dielectric boundary force \
(atom %d) = (%1.12E, %1.12E, %1.12E) kJ/mol/A\n", ifr, lforce[ifr].dbForce[0],
                        lforce[ifr].dbForce[1], lforce[ifr].dbForce[2]);
        }
    }
#endif

    if (refcalcforce == PCF_TOTAL) {
        Vnm_tprint( 1, "  Printing net forces (kJ/mol/A).\n");
        Vnm_tprint( 1, "  Legend:\n");
        Vnm_tprint( 1, "    tot -- Total force\n");
        Vnm_tprint( 1, "    qf  -- Fixed charge force\n");
        Vnm_tprint( 1, "    db  -- Dielectric boundary force\n");
        Vnm_tprint( 1, "    ib  -- Ionic boundary force\n");

        for (ivc=0; ivc<3; ivc++) {
            totforce[ivc] =
            gforce[0].qfForce[ivc] + gforce[0].ibForce[ivc] \
            + gforce[0].dbForce[ivc];
        }

        Vnm_tprint( 1, "  tot %1.12E  %1.12E  %1.12E\n", totforce[0],
                    totforce[1], totforce[2]);
        Vnm_tprint( 1, "  qf  %1.12E  %1.12E  %1.12E\n", gforce[0].qfForce[0],
                    gforce[0].qfForce[1], gforce[0].qfForce[2]);
        Vnm_tprint( 1, "  ib  %1.12E  %1.12E  %1.12E\n", gforce[0].ibForce[0],
                    gforce[0].ibForce[1], gforce[0].ibForce[2]);
        Vnm_tprint( 1, "  db  %1.12E  %1.12E  %1.12E\n", gforce[0].dbForce[0],
                    gforce[0].dbForce[1], gforce[0].dbForce[2]);

    } else if (refcalcforce == PCF_COMPS) {

        Vnm_tprint( 1, "  Printing per-atom forces (kJ/mol/A).\n");
        Vnm_tprint( 1, "  Legend:\n");
        Vnm_tprint( 1, "    tot n -- Total force for atom n\n");
        Vnm_tprint( 1, "    qf  n -- Fixed charge force for atom n\n");
        Vnm_tprint( 1, "    db  n -- Dielectric boundary force for atom n\n");
        Vnm_tprint( 1, "    ib  n -- Ionic boundary force for atom n\n");
        Vnm_tprint( 1, "    tot all -- Total force for system\n");

        totforce[0] = 0.0;
        totforce[1] = 0.0;
        totforce[2] = 0.0;

        for (ifr=0; ifr<refnforce; ifr++) {
            Vnm_tprint( 1, "  qf  %d  %1.12E  %1.12E  %1.12E\n", ifr,
                        gforce[ifr].qfForce[0], gforce[ifr].qfForce[1],
                        gforce[ifr].qfForce[2]);
            Vnm_tprint( 1, "  ib  %d  %1.12E  %1.12E  %1.12E\n", ifr,
                        gforce[ifr].ibForce[0], gforce[ifr].ibForce[1],
                        gforce[ifr].ibForce[2]);
            Vnm_tprint( 1, "  db  %d  %1.12E  %1.12E  %1.12E\n", ifr,
                        gforce[ifr].dbForce[0], gforce[ifr].dbForce[1],
                        gforce[ifr].dbForce[2]);
            Vnm_tprint( 1, "  tot %d  %1.12E  %1.12E  %1.12E\n", ifr,
                        (gforce[ifr].dbForce[0] \
                         + gforce[ifr].ibForce[0] +
                         gforce[ifr].qfForce[0]),
                        (gforce[ifr].dbForce[1] \
                         + gforce[ifr].ibForce[1] +
                         gforce[ifr].qfForce[1]),
                        (gforce[ifr].dbForce[2] \
                         + gforce[ifr].ibForce[2] +
                         gforce[ifr].qfForce[2]));
            for (ivc=0; ivc<3; ivc++) {
                totforce[ivc] += (gforce[ifr].dbForce[ivc] \
                                + gforce[ifr].ibForce[ivc] \
                                  + gforce[ifr].qfForce[ivc]);
            }
        }
        Vnm_tprint( 1, "  tot all %1.12E  %1.12E  %1.12E\n", totforce[0],
                    totforce[1], totforce[2]);
    }

    Vmem_free(VNULL, refnforce, sizeof(AtomForce), (void **)(&lforce));
    Vmem_free(VNULL, refnforce, sizeof(AtomForce), (void **)(&gforce));

    return 1;

}

VPUBLIC int printApolForce(Vcom *com,
                           NOsh *nosh,
                           int nforce[NOSH_MAXCALC],
                           AtomForce *atomForce[NOSH_MAXCALC],
                           int iprint
                          ) {

    int iarg,
        ifr,
        ivc,
        calcid,
        refnforce,
        refcalcforce;
    double temp,
           scalar,
           totforce[3];
    AtomForce *lforce,
              *gforce,
              *aforce;

    if (Vstring_strcasecmp(nosh->apolname[nosh->printcalc[iprint][0]], "") == 0){
        Vnm_tprint( 1, "\nprint APOL force %d ", nosh->printcalc[iprint][0]+1);
    } else {
        Vnm_tprint( 1, "\nprint APOL force %d (%s) ", nosh->printcalc[iprint][0]+1,
                    nosh->apolname[nosh->printcalc[iprint][0]]);
    }
    for (iarg=1; iarg<nosh->printnarg[iprint]; iarg++) {
        if (nosh->printop[iprint][iarg-1] == 0)
            Vnm_tprint(1, "+ ");
        else if (nosh->printop[iprint][iarg-1] == 1)
            Vnm_tprint(1, "- ");
        else {
            Vnm_tprint( 2, "Undefined PRINT operation!\n");
            return 0;
        }
        if (Vstring_strcasecmp(nosh->apolname[nosh->printcalc[iprint][iarg]],
                               "") == 0) {
            Vnm_tprint( 1, "%d ", nosh->printcalc[iprint][iarg]+1);
        } else {
            Vnm_tprint( 1, "%d (%s) ", nosh->printcalc[iprint][iarg]+1,
                        nosh->apolname[nosh->printcalc[iprint][iarg]]);
        }
    }
    Vnm_tprint(1, "end\n");

    /* First, go through and make sure we did the same type of force
        * evaluation in each of the requested calculations */
    calcid = nosh->apol2calc[nosh->printcalc[iprint][0]];
    refnforce = nforce[calcid];
    refcalcforce = nosh->calc[calcid]->apolparm->calcforce;
    if (refcalcforce == ACF_NO) {
        Vnm_tprint( 2, "  Didn't calculate force in calculation \
#%d\n", calcid+1);
        return 0;
    }
    for (iarg=1; iarg<nosh->printnarg[iprint]; iarg++) {
        calcid = nosh->apol2calc[nosh->printcalc[iprint][iarg]-1];
        if (nosh->calc[calcid]->apolparm->calcforce != refcalcforce) {
            Vnm_tprint(2, "  Inconsistent calcforce declarations in \
calculations %d and %d\n", nosh->apol2calc[nosh->printcalc[iprint][0]]+1,
                       calcid+1);
            return 0;
        }
        if (nforce[calcid] != refnforce) {
            Vnm_tprint(2, "  Inconsistent number of forces evaluated in \
calculations %d and %d\n", nosh->apol2calc[nosh->printcalc[iprint][0]]+1,
                       calcid+1);
            return 0;
        }
    }

    /* Now, allocate space to accumulate the forces */
    lforce = (AtomForce *)Vmem_malloc(VNULL, refnforce, sizeof(AtomForce));
    gforce = (AtomForce *)Vmem_malloc(VNULL, refnforce, sizeof(AtomForce));

    /* Now, accumulate the forces */
    calcid = nosh->apol2calc[nosh->printcalc[iprint][0]];
    aforce = atomForce[calcid];
    temp = nosh->calc[calcid]->apolparm->temp;

    /* Load up the first calculation */
    if (refcalcforce == ACF_TOTAL) {
        /* Set to total force */
        for (ivc=0; ivc<3; ivc++) {
            lforce[0].sasaForce[ivc] = aforce[0].sasaForce[ivc];
            lforce[0].savForce[ivc] = aforce[0].savForce[ivc];
            lforce[0].wcaForce[ivc] = aforce[0].wcaForce[ivc];
        }
    } else if (refcalcforce == ACF_COMPS) {
        for (ifr=0; ifr<refnforce; ifr++) {
            for (ivc=0; ivc<3; ivc++) {
                lforce[ifr].sasaForce[ivc] = aforce[ifr].sasaForce[ivc];
                lforce[ifr].savForce[ivc] = aforce[ifr].savForce[ivc];
                lforce[ifr].wcaForce[ivc] = aforce[ifr].wcaForce[ivc];
            }
        }
    }

    /* Load up the rest of the calculations */
    for (iarg=1; iarg<nosh->printnarg[iprint]; iarg++) {
        calcid = nosh->apol2calc[nosh->printcalc[iprint][iarg]];
        temp = nosh->calc[calcid]->apolparm->temp;
        aforce = atomForce[calcid];
        /* Get operation */
        if (nosh->printop[iprint][iarg-1] == 0) scalar = +1.0;
        else if (nosh->printop[iprint][iarg-1] == 1) scalar = -1.0;
        else scalar = 0.0;
        /* Accumulate */
        if (refcalcforce == ACF_TOTAL) {
            /* Set to total force */
            for (ivc=0; ivc<3; ivc++) {
                lforce[0].sasaForce[ivc] += aforce[0].sasaForce[ivc];
                lforce[0].savForce[ivc] += aforce[0].savForce[ivc];
                lforce[0].wcaForce[ivc] += aforce[0].wcaForce[ivc];
            }
        } else if (refcalcforce == ACF_COMPS) {
            for (ifr=0; ifr<refnforce; ifr++) {
                for (ivc=0; ivc<3; ivc++) {
                    lforce[ifr].sasaForce[ivc] += aforce[ifr].sasaForce[ivc];
                    lforce[ifr].savForce[ivc] += aforce[ifr].savForce[ivc];
                    lforce[ifr].wcaForce[ivc] += aforce[ifr].wcaForce[ivc];
                }
            }
        }
    }

    Vnm_tprint( 0, "printForce:  Performing global reduction (sum)\n");
    for (ifr=0; ifr<refnforce; ifr++) {
        Vcom_reduce(com, lforce[ifr].sasaForce, gforce[ifr].sasaForce, 3, 2, 0);
        Vcom_reduce(com, lforce[ifr].savForce, gforce[ifr].savForce, 3, 2, 0);
        Vcom_reduce(com, lforce[ifr].wcaForce, gforce[ifr].wcaForce, 3, 2, 0);
    }

    if (refcalcforce == ACF_TOTAL) {
        Vnm_tprint( 1, "  Printing net forces (kJ/mol/A)\n");
        Vnm_tprint( 1, "  Legend:\n");
        Vnm_tprint( 1, "    tot   -- Total force\n");
        Vnm_tprint( 1, "    sasa  -- SASA force\n");
        Vnm_tprint( 1, "    sav   -- SAV force\n");
        Vnm_tprint( 1, "    wca   -- WCA force\n\n");

        for (ivc=0; ivc<3; ivc++) {
            totforce[ivc] =
            gforce[0].sasaForce[ivc] + gforce[0].savForce[ivc] \
            + gforce[0].wcaForce[ivc];
        }

        Vnm_tprint( 1, "  tot %1.12E  %1.12E  %1.12E\n", totforce[0],
                    totforce[1], totforce[2]);
        Vnm_tprint( 1, "  sasa  %1.12E  %1.12E  %1.12E\n", gforce[0].sasaForce[0],
                    gforce[0].sasaForce[1], gforce[0].sasaForce[2]);
        Vnm_tprint( 1, "  sav  %1.12E  %1.12E  %1.12E\n", gforce[0].savForce[0],
                    gforce[0].savForce[1], gforce[0].savForce[2]);
        Vnm_tprint( 1, "  wca  %1.12E  %1.12E  %1.12E\n", gforce[0].wcaForce[0],
                    gforce[0].wcaForce[1], gforce[0].wcaForce[2]);

    } else if (refcalcforce == ACF_COMPS) {

        Vnm_tprint( 1, "  Printing per atom forces (kJ/mol/A)\n");
        Vnm_tprint( 1, "  Legend:\n");
        Vnm_tprint( 1, "    tot   n -- Total force for atom n\n");
        Vnm_tprint( 1, "    sasa  n -- SASA force for atom n\n");
        Vnm_tprint( 1, "    sav   n -- SAV force for atom n\n");
        Vnm_tprint( 1, "    wca   n -- WCA force for atom n\n");
        Vnm_tprint( 1, "    tot all -- Total force for system\n");

        //Vnm_tprint( 1, "    gamma, pressure, bconc are: %f %f %f\n\n",
        //			gamma,press,bconc);

        totforce[0] = 0.0;
        totforce[1] = 0.0;
        totforce[2] = 0.0;

        for (ifr=0; ifr<refnforce; ifr++) {
            Vnm_tprint( 1, "  sasa  %d  %1.12E  %1.12E  %1.12E\n", ifr,
                        gforce[ifr].sasaForce[0], gforce[ifr].sasaForce[1],
                        gforce[ifr].sasaForce[2]);
            Vnm_tprint( 1, "  sav   %d  %1.12E  %1.12E  %1.12E\n", ifr,
                        gforce[ifr].savForce[0], gforce[ifr].savForce[1],
                        gforce[ifr].savForce[2]);
            Vnm_tprint( 1, "  wca   %d  %1.12E  %1.12E  %1.12E\n", ifr,
                        gforce[ifr].wcaForce[0], gforce[ifr].wcaForce[1],
                        gforce[ifr].wcaForce[2]);
            Vnm_tprint( 1, "  tot   %d  %1.12E  %1.12E  %1.12E\n", ifr,
                        (gforce[ifr].wcaForce[0] \
                         + gforce[ifr].savForce[0] +
                         gforce[ifr].sasaForce[0]),
                        (gforce[ifr].wcaForce[1] \
                         + gforce[ifr].savForce[1] +
                         gforce[ifr].sasaForce[1]),
                        (gforce[ifr].wcaForce[2] \
                         + gforce[ifr].savForce[2] +
                         gforce[ifr].sasaForce[2]));
            for (ivc=0; ivc<3; ivc++) {
                totforce[ivc] += (gforce[ifr].wcaForce[ivc] \
                                + gforce[ifr].savForce[ivc] \
                                  + gforce[ifr].sasaForce[ivc]);
            }
        }
        Vnm_tprint( 1, "  tot all  %1.12E  %1.12E  %1.12E\n", totforce[0],
                    totforce[1], totforce[2]);
    }

    Vmem_free(VNULL, refnforce, sizeof(AtomForce), (void **)(&lforce));
    Vmem_free(VNULL, refnforce, sizeof(AtomForce), (void **)(&gforce));

    return 1;
}

#ifdef HAVE_MC_H


VPUBLIC void killFE(NOsh *nosh,
                    Vpbe *pbe[NOSH_MAXCALC],
                    Vfetk *fetk[NOSH_MAXCALC],
                    Gem *gm[NOSH_MAXMOL]
                   ) {

    int i;

#ifndef VAPBSQUIET
    Vnm_tprint(1, "Destroying finite element structures.\n");
#endif

    for(i=0;i<nosh->ncalc;i++){
        Vpbe_dtor(&(pbe[i]));
        Vfetk_dtor(&(fetk[i]));
    }
    for (i=0; i<nosh->nmesh; i++) {
        Gem_dtor(&(gm[i]));
}
}

/**
 * @brief  Initialize FE solver objects
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @bug  THIS FUNCTION IS HARD-CODED TO SOLVE LRPBE
 * @todo  THIS FUNCTION IS HARD-CODED TO SOLVE LRPBE
 */
VPUBLIC Vrc_Codes initFE(int icalc, /**< Index in pb, fetk to initialize (calculation index) */
                         NOsh *nosh, /**< Master parmaeter object */
                         FEMparm *feparm, /**< FE-specific parameters */
                         PBEparm *pbeparm, /**< Generic PBE parameters */
                         Vpbe *pbe[NOSH_MAXCALC],  /**< Array of PBE objects */
                         Valist *alist[NOSH_MAXMOL], /**< Array of atom lists */
                         Vfetk *fetk[NOSH_MAXCALC]  /**< Array of finite element objects */
                        ) {

    int iatom,                  /**< Loop counter */
        imesh,                  /**< Mesh ID */
        i,                      /**< Loop counter */
        j,                      /**< Loop counter */
        theMol,                 /**< Molecule ID */
        focusFlag = 0;          /**< @todo document */
    Vio *sock = VNULL;          /**< I/O socket for reading MCSF mesh data */
    size_t bytesTotal,          /**< Total bytes used by this operation */
           highWater;           /**< High-water memory usage for this operation */
    Vfetk_MeshLoad meshType;    /**< The type of mesh being used (see struct for enum values) */
    double length[3],           /**< @todo document */
           center[3],           /**< @todo document */
           sparm,               /**< @todo document */
           q,                   /**< @todo document */
           iparm = 0.0;         /**< @todo document */
    Vrc_Codes vrc;              /**< Return codes for function calls (see struct for enum value) */
    Valist *myalist;            /**< List of atoms being operated on */
    Vatom *atom = VNULL;        /**< Atom/molecule being operated on */

    Vnm_tstart(27, "Setup timer");

    /* Print some warning messages */
    if (pbeparm->useDielMap)  Vnm_tprint(2, "FEM ignoring dielectric map!\n");
    if (pbeparm->useKappaMap)  Vnm_tprint(2, "FEM ignoring kappa map!\n");
    if (pbeparm->useChargeMap)  Vnm_tprint(2, "FEM ignoring charge map!\n");

    /* Fix mesh center for "GCENT MOL #" types of declarations. */
    Vnm_tprint(0, "Re-centering mesh...\n");
    theMol = pbeparm->molid-1;
    myalist = alist[theMol];
    for (j=0; j<3; j++) {
        if (theMol < nosh->nmol) {
            center[j] = (myalist)->center[j];
        } else{
            Vnm_tprint(2, "ERROR!  Bogus molecule number (%d)!\n",
                       (theMol+1));
            return VRC_FAILURE;
        }
    }

    /* Check for completely-neutral molecule */
    q = 0;
    for (iatom=0; iatom<Valist_getNumberAtoms(myalist); iatom++) {
        atom = Valist_getAtom(myalist, iatom);
        q += VSQR(Vatom_getCharge(atom));
    }
    /* D. Gohara 10/22/09 - disabled
    if (q < (1e-6)) {
        Vnm_tprint(2, "Molecule #%d is uncharged!\n", pbeparm->molid);
        Vnm_tprint(2, "Sum square charge = %g!\n", q);
        return VRC_FAILURE;
    }
    */

    /* Set the femparm pkey value based on the presence of an HB solver */
#ifdef USE_HB
    feparm->pkey = 1;
#else
    feparm->pkey = 0;
#endif

    /* Set up PBE object */
    Vnm_tprint(0, "Setting up PBE object...\n");
    if (pbeparm->srfm == VSM_SPLINE) {
        sparm = pbeparm->swin;
    }
    else {
        sparm = pbeparm->srad;
    }
    if (pbeparm->nion > 0) {
        iparm = pbeparm->ionr[0];
    }

    pbe[icalc] = Vpbe_ctor(myalist, pbeparm->nion,
                           pbeparm->ionc, pbeparm->ionr, pbeparm->ionq,
                           pbeparm->temp, pbeparm->pdie,
                           pbeparm->sdie, sparm, focusFlag, pbeparm->sdens,
                           pbeparm->zmem, pbeparm->Lmem, pbeparm->mdie,
                           pbeparm->memv);

    /* Print a few derived parameters */
    Vnm_tprint(1, "  Debye length:  %g A\n", Vpbe_getDeblen(pbe[icalc]));

    /* Set up FEtk objects */
    Vnm_tprint(0, "Setting up FEtk object...\n");
    fetk[icalc] = Vfetk_ctor(pbe[icalc], pbeparm->pbetype);
    Vfetk_setParameters(fetk[icalc], pbeparm, feparm);

    /* Build mesh - this merely loads an MCSF file from an external source if one is specified or uses the
     * current molecule and sets center/length values based on that molecule if no external source is
     * specified. */
    Vnm_tprint(0, "Setting up mesh...\n");
    sock = VNULL;
    if (feparm->useMesh) {
        imesh = feparm->meshID-1;
        meshType = VML_EXTERNAL;
        for (i=0; i<3; i++) {
            center[i] = 0.0;
            length[i] = 0.0;
        }
        Vnm_print(0, "Using mesh %d (%s) in calculation.\n", imesh+1,
                  nosh->meshpath[imesh]);
        switch (nosh->meshfmt[imesh]) {
            case VDF_DX:
                Vnm_tprint(2, "DX finite element mesh input not supported yet!\n");
                return VRC_FAILURE;
            case VDF_UHBD:
                Vnm_tprint( 2, "UHBD finite element mesh input not supported!\n");
                return VRC_FAILURE;
            case VDF_AVS:
                Vnm_tprint( 2, "AVS finite element mesh input not supported!\n");
                return VRC_FAILURE;
            case VDF_MCSF:
                Vnm_tprint(1, "Reading MCSF-format input finite element mesh from %s.\n",
                           nosh->meshpath[imesh]);
                sock = Vio_ctor("FILE", "ASC", VNULL, nosh->meshpath[imesh], "r");
                if (sock == VNULL) {
                    Vnm_print(2, "Problem opening virtual socket %s!\n",
                              nosh->meshpath[imesh]);
                    return VRC_FAILURE;
                }
                if (Vio_accept(sock, 0) < 0) {
                    Vnm_print(2, "Problem accepting virtual socket %s!\n",
                              nosh->meshpath[imesh]);
                    return VRC_FAILURE;
                }
                break;

            default:
                Vnm_tprint( 2, "Invalid data format (%d)!\n",
                           nosh->meshfmt[imesh]);
                return VRC_FAILURE;
        }
    } else { /* if (feparm->useMesh) */
        meshType = VML_DIRICUBE;
        for (i=0; i<3; i++) {
            center[i] = alist[theMol]->center[i];
            length[i] = feparm->glen[i];
        }
    }

    /* Load the mesh with a particular center and vertex length using the provided input socket */
    vrc = Vfetk_loadMesh(fetk[icalc], center, length, meshType, sock);
    if (vrc == VRC_FAILURE) {
        Vnm_print(2, "Error constructing finite element mesh!\n");
        return VRC_FAILURE;
    }
    //Vnm_redirect(0); // Redirect output to io.mc
    Gem_shapeChk(fetk[icalc]->gm); // Traverse simplices and check shapes using the geometry manager.
    //Vnm_redirect(1);

    /* Uniformly refine the mesh a bit */
    for (j=0; j<2; j++) {
        /* AM_* calls below are from the MC package, mc/src/nam/am.c.  Note that these calls actually are
         * wrappers around Aprx_* functions found in MC as well. */
        /* Mark the mesh for needed refinements */
        AM_markRefine(fetk[icalc]->am, 0, -1, 0, 0.0);
        /* Do actual mesh refinement */
        AM_refine(fetk[icalc]->am, 2, 0, feparm->pkey);
        //Vnm_redirect(0); // Redirect output to io.mc
        Gem_shapeChk(fetk[icalc]->gm); // Traverse simplices and check shapes using the geometry manager.
        //Vnm_redirect(1);
    }

    /* Setup time statistics */
    Vnm_tstop(27, "Setup timer");

    /* Memory statistics */
    bytesTotal = Vmem_bytesTotal();
    highWater = Vmem_highWaterTotal();

#ifndef VAPBSQUIET
    Vnm_tprint( 1, "  Current memory usage:  %4.3f MB total, \
%4.3f MB high water\n", (double)(bytesTotal)/(1024.*1024.),
                (double)(highWater)/(1024.*1024.));
#endif

    return VRC_SUCCESS;
}

VPUBLIC void printFEPARM(int icalc,
                         NOsh *nosh,
                         FEMparm *feparm,
                         Vfetk *fetk[NOSH_MAXCALC]
                        ) {

    Vnm_tprint(1, "  Domain size:  %g A x %g A x %g A\n",
               feparm->glen[0], feparm->glen[1],
               feparm->glen[2]);
    switch(feparm->ekey) {
        case FET_SIMP:
            Vnm_tprint(1, "  Per-simplex error tolerance:  %g\n", feparm->etol);
            break;
        case FET_GLOB:
            Vnm_tprint(1, "  Global error tolerance:  %g\n", feparm->etol);
            break;
        case FET_FRAC:
            Vnm_tprint(1, "  Fraction of simps to refine:  %g\n", feparm->etol);
            break;
        default:
            Vnm_tprint(2, "Invalid ekey (%d)!\n", feparm->ekey);
            VASSERT(0);
            break;
    }
    switch(feparm->akeyPRE) {
        case FRT_UNIF:
            Vnm_tprint(1, "  Uniform pre-solve refinement.\n");
            break;
        case FRT_GEOM:
            Vnm_tprint(1, "  Geometry-based pre-solve refinement.\n");
            break;
        default:
            Vnm_tprint(2, "Invalid akeyPRE (%d)!\n", feparm->akeyPRE);
            VASSERT(0);
            break;
    }
    switch(feparm->akeySOLVE) {
        case FRT_RESI:
            Vnm_tprint(1, "  Residual-based a posteriori refinement.\n");
            break;
        case FRT_DUAL:
            Vnm_tprint(1, "  Dual-based a posteriori refinement.\n");
            break;
        case FRT_LOCA:
            Vnm_tprint(1, "  Local-based a posteriori refinement.\n");
            break;
        default:
            Vnm_tprint(2, "Invalid akeySOLVE (%d)!\n", feparm->akeySOLVE);
            break;
    }
    Vnm_tprint(1, "  Refinement of initial mesh to ~%d vertices\n",
               feparm->targetNum);
    Vnm_tprint(1, "  Geometry-based refinment lower bound:  %g A\n",
               feparm->targetRes);
    Vnm_tprint(1, "  Maximum number of solve-estimate-refine cycles:  %d\n",
               feparm->maxsolve);
    Vnm_tprint(1, "  Maximum number of vertices in mesh:  %d\n",
               feparm->maxvert);

    /* FOLLOWING IS SOLVER-RELATED; BAIL IF NOT SOLVING */
    if (nosh->bogus)  return;
#ifdef USE_HB
    Vnm_tprint(1, "  HB linear solver:  AM_hPcg\n");
#else
    Vnm_tprint(1, "  Non-HB linear solver:  ");
    switch (fetk[icalc]->lkey) {
            case VLT_SLU:
                Vnm_print(1, "SLU direct\n");
                break;
            case VLT_MG:
                Vnm_print(1, "multigrid\n");
                break;
            case VLT_CG:
                Vnm_print(1, "conjugate gradient\n");
                break;
            case VLT_BCG:
                Vnm_print(1, "BiCGStab\n");
                break;
            default:
                Vnm_print(1, "???\n");
                break;
    }
#endif

    Vnm_tprint(1, "  Linear solver tol.:  %g\n", fetk[icalc]->ltol);
    Vnm_tprint(1, "  Linear solver max. iters.:  %d\n", fetk[icalc]->lmax);
    Vnm_tprint(1, "  Linear solver preconditioner:  ");
    switch (fetk[icalc]->lprec) {
        case VPT_IDEN:
            Vnm_print(1, "identity\n");
            break;
        case VPT_DIAG:
            Vnm_print(1, "diagonal\n");
            break;
        case VPT_MG:
            Vnm_print(1, "multigrid\n");
            break;
        default:
            Vnm_print(1, "???\n");
            break;
    }
    Vnm_tprint(1, "  Nonlinear solver:  ");
    switch (fetk[icalc]->nkey) {
        case VNT_NEW:
            Vnm_print(1, "newton\n");
            break;
        case VNT_INC:
            Vnm_print(1, "incremental\n");
            break;
        case VNT_ARC:
            Vnm_print(1, "pseudo-arclength\n");
            break;
        default:
            Vnm_print(1, "??? ");
            break;
    }
    Vnm_tprint(1, "  Nonlinear solver tol.:  %g\n", fetk[icalc]->ntol);
    Vnm_tprint(1, "  Nonlinear solver max. iters.:  %d\n", fetk[icalc]->nmax);
    Vnm_tprint(1, "     Initial guess:  ");
    switch (fetk[icalc]->gues) {
        case VGT_ZERO:
            Vnm_tprint(1, "zero\n");
            break;
        case VGT_DIRI:
            Vnm_tprint(1, "boundary function\n");
            break;
        case VGT_PREV:
            Vnm_tprint(1, "interpolated previous solution\n");
            break;
        default:
            Vnm_tprint(1, "???\n");
            break;
    }

}

VPUBLIC int partFE(int icalc, NOsh *nosh, FEMparm *feparm,
                   Vfetk *fetk[NOSH_MAXCALC]) {

    Vfetk_setAtomColors(fetk[icalc]);
    return 1;
}

VPUBLIC int preRefineFE(int icalc, /* Calculation index */
                        FEMparm *feparm, /* FE-specific parameters */
                        Vfetk *fetk[NOSH_MAXCALC] /* Array of FE solver objects */
                       ) {

    int nverts, /**< Number of vertices in the mesh geometry */
        marked; /**< Essentially a boolean; indicates whether further refinement is required after
                  *  running MC's refinement algorithm. */

    /* Based on the refinement type, alert the user to what we're tryng to refine with. */
    switch(feparm->akeyPRE) {
        case FRT_UNIF:
            Vnm_tprint(1, "  Commencing uniform refinement to %d verts.\n",
                       feparm->targetNum);
            break;
        case FRT_GEOM:
            Vnm_tprint(1, "  Commencing geometry-based refinement to %d \
verts or %g A resolution.\n", feparm->targetNum, feparm->targetRes);
            break;
        case FRT_DUAL:
            Vnm_tprint(2, "What?  You can't do a posteriori error estimation \
before you solve!\n");
            VASSERT(0);
            break;
        case FRT_RESI:
        case FRT_LOCA:
        default:
            VASSERT(0);
            break;
    }

    /**
     * TODO: could this be optimized by moving nverts out of the loop to just
     * above this initial printout? This depends heavily on whether the number
     * of vertices can change during the calculation. - PCE
     */
    Vnm_tprint(1, "  Initial mesh has %d vertices\n",
               Gem_numVV(fetk[icalc]->gm));

    /* As long as we have simplices marked that need to be refined, run MC's
     * AM_markRefine against our data until we hit the error or size tolerance
     * for the refinement. */
    while (1) {
        nverts = Gem_numVV(fetk[icalc]->gm);
        if (nverts > feparm->targetNum) {
            Vnm_tprint(1, "  Hit vertex number limit.\n");
            break;
        }
        marked = AM_markRefine(fetk[icalc]->am, feparm->akeyPRE, -1,
                               feparm->ekey, feparm->etol);
        if (marked == 0) {
            Vnm_tprint(1, "  Marked 0 simps; hit error/size tolerance.\n");
            break;
        }
        Vnm_tprint(1, "    Have %d verts, marked %d.  Refining...\n", nverts,
                   marked);
        AM_refine(fetk[icalc]->am, 0, 0, feparm->pkey);
    }

    nverts = Gem_numVV(fetk[icalc]->gm);
    Vnm_tprint(1, "  Done refining; have %d verts.\n", nverts);

    return 1;
}


/**
 * Call MC's mesh solving equations depending upon the type of PBE we're
 * dealing with.
 */
VPUBLIC int solveFE(int icalc, /**< Calculation index */
                    PBEparm *pbeparm, /**< PBE-specific parameters */
                    FEMparm *feparm, /**< FE-specific parameters */
                    Vfetk *fetk[NOSH_MAXCALC] /**< Array of FE solver objects */
                   ) {

    int lkeyHB = 3,  /**<  AM_hPcg */
        meth = 2,  /**< Coarse-grid solver; 0 = SLU, 1 = MG, 2 = CG, 3 = BCG, 4 = PCG, 5 = PBCG */
        prob = 0,  /**< Primal problem */
        prec = 0;  /** < Preconditioner; 0 = identity. */

    if ((pbeparm->pbetype==PBE_NPBE) ||
        (pbeparm->pbetype == PBE_NRPBE) ||
        (pbeparm->pbetype == PBE_SMPBE)) {

        /* Call MC's nonlinear solver - mc/src/nam/nsolv.c */
        AM_nSolve(
                  fetk[icalc]->am,
                  fetk[icalc]->nkey,
                  fetk[icalc]->nmax,
                  fetk[icalc]->ntol,
                  meth,
                  fetk[icalc]->lmax,
                  fetk[icalc]->ltol,
                  prec,
                  fetk[icalc]->gues,
                  fetk[icalc]->pjac
                  );
    } else if ((pbeparm->pbetype==PBE_LPBE) ||
               (pbeparm->pbetype==PBE_LRPBE)) {
        /* Note: USEHB is a compile time defined macro. The program flow
        is to always take the route using AM_hlSolve when the solver
        is linear. D. Gohara 6/29/06
        */
#ifdef USE_HB
        Vnm_print(2, "SORRY!  DON'T USE HB!!!\n");
        VASSERT(0);

        /* Call MC's hierarchical linear solver - mc/src/nam/lsolv.c */
        AM_hlSolve(fetk[icalc]->am, meth, lkeyHB, fetk[icalc]->lmax,
            fetk[icalc]->ltol, fetk[icalc]->gues, fetk[icalc]->pjac);
#else

        /* Call MC's linear solver - mc/src/nam/lsolv.c */
        AM_lSolve(
                  fetk[icalc]->am,
                  prob,
                  meth,
                  fetk[icalc]->lmax,
                  fetk[icalc]->ltol,
                  prec,
                  fetk[icalc]->gues,
                  fetk[icalc]->pjac
                  );
#endif
    }

    return 1;
}

/**
 * Calculates the electrostatic energies from an FE calculation.
 */
VPUBLIC int energyFE(NOsh *nosh, /**< Object with parsed input file parameters */
                     int icalc, /**< Calculation index */
                     Vfetk *fetk[NOSH_MAXCALC], /**< FE object array */
                     int *nenergy, /**< Set to number of entries in energy arrays */
                     double *totEnergy, /**< Set to total energy (in kT) */
                     double *qfEnergy, /**< Set to charge-potential energy (in kT) */
                     double *qmEnergy, /**< Set to mobile ion energy (in kT) */
                     double *dielEnergy /**< Set to polarization energy (in kT) */
                    ) {

    FEMparm *feparm = nosh->calc[icalc]->femparm; /**< FE-specific parameters */
    PBEparm *pbeparm = nosh->calc[icalc]->pbeparm; /**< PBE-specific parameters */

    *nenergy = 1;
    *totEnergy = 0;

    /**
     * If we're not ignoring this particular NOsh object because it has been
     * rendered invalid, call the Vfetk object's energy calculation function.
     * The flag differences specified have to do with setting specific calculation
     * restrictions (see color variable documentation in function code).
     */
    if (nosh->bogus == 0) {
        if ((pbeparm->pbetype==PBE_NPBE) ||
            (pbeparm->pbetype==PBE_NRPBE) ||
            (pbeparm->pbetype == PBE_SMPBE)) {
            *totEnergy = Vfetk_energy(fetk[icalc], -1, 1); /* Last parameter indicates NPBE */
        } else if ((pbeparm->pbetype==PBE_LPBE) ||
                   (pbeparm->pbetype==PBE_LRPBE)) {
            *totEnergy = Vfetk_energy(fetk[icalc], -1, 0); /* Last parameter indicates LPBE */
        } else {
            VASSERT(0);
        }

#ifndef VAPBSQUIET
        Vnm_tprint(1, "      Total electrostatic energy = %1.12E kJ/mol\n",
                   Vunit_kb*pbeparm->temp*(1e-3)*Vunit_Na*(*totEnergy));
        fflush(stdout);
#endif
    }

    if (pbeparm->calcenergy == PCE_COMPS) {
        Vnm_tprint(2, "Error!  Verbose energy evaluation not available for FEM yet!\n");
        Vnm_tprint(2, "E-mail nathan.baker@pnl.gov if you want this.\n");
        *qfEnergy = 0;
        *qmEnergy = 0;
        *dielEnergy = 0;
    } else {
        *nenergy = 0;
    }
    return 1;
}

/**
 * Estimates the error, marks the mesh, and refines the mesh after solving.
 * @return  1 if successful, 0 otherwise -- note that a 0 will likely imply
 * that either the max number of vertices have been met or no vertices were
 * marked for refinement.  In either case, this should not be treated as a
 * fatal error.
 */
VPUBLIC int postRefineFE(int icalc, /**< Calculation index */
                         FEMparm *feparm, /**< FE-specific parameters */
                         Vfetk *fetk[NOSH_MAXCALC] /**< Array of FE solver objects */
                        ) {

    int nverts, /**< Number of vertices in the molecular geometry */
        marked; /**< Whether vertices are marked for refinement */

    nverts = Gem_numVV(fetk[icalc]->gm);
    if (nverts > feparm->maxvert) {
        Vnm_tprint(1, "    Current number of vertices (%d) exceeds max (%d)!\n",
                   nverts, feparm->maxvert);
        return 0;
    }
    Vnm_tprint(1, "      Mesh currently has %d vertices\n", nverts);

    switch(feparm->akeySOLVE) {
        case FRT_UNIF:
            Vnm_tprint(1, "      Commencing uniform refinement.\n");
            break;
        case FRT_GEOM:
            Vnm_tprint(1, "      Commencing geometry-based refinement.\n");
            break;
        case FRT_RESI:
            Vnm_tprint(1, "      Commencing residual-based refinement.\n");
            break;
        case FRT_DUAL:
            Vnm_tprint(1, "      Commencing dual problem-based refinement.\n");
            break;
        case FRT_LOCA:
            Vnm_tprint(1, "      Commencing local-based refinement.\n.");
            break;
        default:
            Vnm_tprint(2, "      Error -- unknown refinement type (%d)!\n",
                       feparm->akeySOLVE);
            return 0;
            break;
    }

    /* Run MC's refinement algorithm */
    marked = AM_markRefine(fetk[icalc]->am, feparm->akeySOLVE, -1,
                           feparm->ekey, feparm->etol);
    if (marked == 0) {
        Vnm_tprint(1, "      Marked 0 simps; hit error/size tolerance.\n");
        return 0;
    }
    Vnm_tprint(1, "      Have %d verts, marked %d.  Refining...\n", nverts,
               marked);
    AM_refine(fetk[icalc]->am, 0, 0, feparm->pkey);
    nverts = Gem_numVV(fetk[icalc]->gm);
    Vnm_tprint(1, "      Done refining; have %d verts.\n", nverts);
    //Vnm_redirect(0); // Redirect output to io.mc
    Gem_shapeChk(fetk[icalc]->gm); // Traverse simplices and check shapes using the geometry manager.
    //Vnm_redirect(1);

    return 1;
}

/**
 * Write FEM data to file.
 */
VPUBLIC int writedataFE(int rank, /**< Rank of processor (for parallel runs) */
                        NOsh *nosh, /**< NOsh object */
                        PBEparm *pbeparm, /**< PBE-specific parameters */
                        Vfetk *fetk /**< FEtk object (with solution) */
                       ) {

    char writestem[VMAX_ARGLEN];    /**< @todo document */
    char outpath[VMAX_ARGLEN];      /**< @todo document */
    int i,                          /**< Loop counter */
        writeit;                    /**< Flag indicating whether data can be written to output */
    AM *am;                         /**< @todo document */
    Bvec *vec;                      /**< @todo document */

    if (nosh->bogus) return 1;

    am = fetk->am;
    vec = am->w0;

    for (i=0; i<pbeparm->numwrite; i++) {

        writeit = 1;

        switch (pbeparm->writetype[i]) {

            case VDT_CHARGE:

                Vnm_tprint(2, "    Sorry; can't write charge distribution for FEM!\n");
                writeit = 0;
                break;

            case VDT_POT:

                Vnm_tprint(1, "    Writing potential to ");
                Vfetk_fillArray(fetk, vec, VDT_POT);
                break;

            case VDT_SMOL:

                Vnm_tprint(1, "    Writing molecular accessibility to ");
                Vfetk_fillArray(fetk, vec, VDT_SMOL);
                break;

            case VDT_SSPL:

                Vnm_tprint(1, "    Writing spline-based accessibility to ");
                Vfetk_fillArray(fetk, vec, VDT_SSPL);
                break;

            case VDT_VDW:

                Vnm_tprint(1, "    Writing van der Waals accessibility to ");
                Vfetk_fillArray(fetk, vec, VDT_VDW);
                break;

            case VDT_IVDW:

                Vnm_tprint(1, "    Writing ion accessibility to ");
                Vfetk_fillArray(fetk, vec, VDT_IVDW);
                break;

            case VDT_LAP:

                Vnm_tprint(2, "    Sorry; can't write charge distribution for FEM!\n");
                writeit = 0;
                break;

            case VDT_EDENS:

                Vnm_tprint(2, "    Sorry; can't write energy density for FEM!\n");
                writeit = 0;
                break;

            case VDT_NDENS:

                Vnm_tprint(1, "    Writing number density to ");
                Vfetk_fillArray(fetk, vec, VDT_NDENS);
                break;

            case VDT_QDENS:

                Vnm_tprint(1, "    Writing charge density to ");
                Vfetk_fillArray(fetk, vec, VDT_QDENS);
                break;

            case VDT_DIELX:

                Vnm_tprint(2, "    Sorry; can't write x-shifted dielectric map for FEM!\n");
                writeit = 0;
                break;

            case VDT_DIELY:

                Vnm_tprint(2, "    Sorry; can't write y-shifted dielectric map for FEM!\n");
                writeit = 0;
                break;

            case VDT_DIELZ:

                Vnm_tprint(2, "    Sorry; can't write z-shifted dielectric map for FEM!\n");
                writeit = 0;
                break;

            case VDT_KAPPA:

                Vnm_tprint(1, "    Sorry; can't write kappa map for FEM!\n");
                writeit = 0;
                break;

            case VDT_ATOMPOT:

                Vnm_tprint(1, "    Sorry; can't write atom potentials for FEM!\n");
                writeit = 0;
                break;

            default:

                Vnm_tprint(2, "Invalid data type for writing!\n");
                writeit = 0;
                return 0;
        }

        if (!writeit) return 0;


#ifdef HAVE_MPI_H
        sprintf(writestem, "%s-PE%d", pbeparm->writestem[i], rank);
#else
        if(nosh->ispara){
            sprintf(writestem, "%s-PE%d", pbeparm->writestem[i],nosh->proc_rank);
        }else{
            sprintf(writestem, "%s", pbeparm->writestem[i]);
        }
#endif

        switch (pbeparm->writefmt[i]) {

            case VDF_DX:
                sprintf(outpath, "%s.%s", writestem, "dx");
                Vnm_tprint(1, "%s\n", outpath);
                Vfetk_write(fetk, "FILE", "ASC", VNULL, outpath, vec, VDF_DX);
                break;

            case VDF_AVS:
                sprintf(outpath, "%s.%s", writestem, "ucd");
                Vnm_tprint(1, "%s\n", outpath);
                Vfetk_write(fetk, "FILE", "ASC", VNULL, outpath, vec, VDF_AVS);
                break;

            case VDF_UHBD:
                Vnm_tprint(2, "UHBD format not supported for FEM!\n");
                break;

            case VDF_MCSF:
                Vnm_tprint(2, "MCSF format not supported yet!\n");
                break;

            default:
                Vnm_tprint(2, "Bogus data format (%d)!\n",
                           pbeparm->writefmt[i]);
                break;
        }

    }

    return 1;
}
#endif /* ifdef HAVE_MCX_H */

VPUBLIC int initAPOL(NOsh *nosh, /**< Input parameter object */
                     Vmem *mem, /**< Memory manager */
                     Vparam *param, /**< Atom parameters */
                     APOLparm *apolparm, /**< Apolar calculation parameters */
                     int *nforce, /**< Number of force calculations */
                     AtomForce **atomForce, /**< Atom force storage object */
                     Valist *alist /**< Atom list */
                    ) {
    int i,          /**< @todo document */
        natoms,     /**< Number of atoms */
        len,        /**< Used to capture length of loops to prevent multiple calls in counters */
        inhash[3],  /**< @todo document */
        rc = 0;     /**< @todo document */

    time_t ts;      /**< Temporary timing variable for debugging (PCE) */
    Vclist *clist = VNULL;  /**< @todo document */
    Vacc *acc = VNULL;      /**< @todo document */
    Vatom *atom = VNULL;    /**< @todo document */
    Vparam_AtomData *atomData = VNULL;  /**< @todo document */

    double sasa,        /**< @todo document */
           sav,         /**< @todo document */
           nhash[3],    /**< @todo document */
           sradPad,     /**< @todo document */
           x,           /**< @todo document */
           y,           /**< @todo document */
           z,           /**< @todo document */
           atomRadius,  /**< @todo document */
           srad,        /**< @todo document */
           *atomsasa,   /**< @todo document */
           *atomwcaEnergy,  /**< @todo document */
           energy = 0.0,    /**< WCA energy per atom */
           dist,        /**< @todo document */
           charge,      /**< @todo document */
           xmin,        /**< @todo document */
           xmax,        /**< @todo document */
           ymin,        /**< @todo document */
           ymax,        /**< @todo document */
           zmin,        /**< @todo document */
           zmax,        /**< @todo document */
           disp[3],     /**< @todo document */
           center[3],   /**< @todo document */
           soluteXlen,  /**< @todo document */
           soluteYlen,  /**< @todo document */
           soluteZlen;  /**< @todo document */

    atomsasa = (double *)Vmem_malloc(VNULL, Valist_getNumberAtoms(alist), sizeof(double));
    atomwcaEnergy = (double *)Vmem_malloc(VNULL, Valist_getNumberAtoms(alist), sizeof(double));

    /* Determine solute length and charge*/
    atom = Valist_getAtom(alist, 0);
    xmin = Vatom_getPosition(atom)[0];
    xmax = Vatom_getPosition(atom)[0];
    ymin = Vatom_getPosition(atom)[1];
    ymax = Vatom_getPosition(atom)[1];
    zmin = Vatom_getPosition(atom)[2];
    zmax = Vatom_getPosition(atom)[2];
    charge = 0;
    natoms = Valist_getNumberAtoms(alist);

    for (i=0; i < natoms; i++) {
        atom = Valist_getAtom(alist, i);
        atomRadius = Vatom_getRadius(atom);
        x = Vatom_getPosition(atom)[0];
        y = Vatom_getPosition(atom)[1];
        z = Vatom_getPosition(atom)[2];
        if ((x+atomRadius) > xmax) xmax = x + atomRadius;
        if ((x-atomRadius) < xmin) xmin = x - atomRadius;
        if ((y+atomRadius) > ymax) ymax = y + atomRadius;
        if ((y-atomRadius) < ymin) ymin = y - atomRadius;
        if ((z+atomRadius) > zmax) zmax = z + atomRadius;
        if ((z-atomRadius) < zmin) zmin = z - atomRadius;
        disp[0] = (x - center[0]);
        disp[1] = (y - center[1]);
        disp[2] = (z - center[2]);
        dist = (disp[0]*disp[0]) + (disp[1]*disp[1]) + (disp[2]*disp[2]);
        dist = VSQRT(dist) + atomRadius;
        charge += Vatom_getCharge(Valist_getAtom(alist, i));
    }
    soluteXlen = xmax - xmin;
    soluteYlen = ymax - ymin;
    soluteZlen = zmax - zmin;

    /* Set up the hash table for the cell list */
    Vnm_print(0, "APOL: Setting up hash table and accessibility object...\n");
    nhash[0] = soluteXlen/0.5;
    nhash[1] = soluteYlen/0.5;
    nhash[2] = soluteZlen/0.5;
    for (i=0; i<3; i++) inhash[i] = (int)(nhash[i]);

    for (i=0;i<3;i++){
        if (inhash[i] < 3) inhash[i] = 3;
        if (inhash[i] > MAX_HASH_DIM) inhash[i] = MAX_HASH_DIM;
    }

    /* Pad the radius by 2x the maximum displacement value */
    srad = apolparm->srad;
    sradPad = srad + (2*apolparm->dpos);
    clist = Vclist_ctor(alist, sradPad , inhash, CLIST_AUTO_DOMAIN,
                                    VNULL, VNULL);
    acc = Vacc_ctor(alist, clist, apolparm->sdens);

    /* Get WAT (water) LJ parameters from Vparam object */
    if (param == VNULL && (apolparm->bconc != 0.0)) {
        Vnm_tprint(2, "initAPOL:  Got NULL Vparam object!\n");
        Vnm_tprint(2, "initAPOL:  You are performing an apolar calculation with the van der Waals integral term,\n");
        Vnm_tprint(2, "initAPOL:  this term requires van der Waals parameters which are not available from the \n");
        Vnm_tprint(2, "initAPOL:  PQR file. Therefore, you need to supply a parameter file with the parm keyword,\n");
        Vnm_tprint(2, "initAPOL:  for example,\n");
        Vnm_tprint(2, "initAPOL:    read parm flat amber94.dat end\n");
        Vnm_tprint(2, "initAPOL:  where the relevant parameter files can be found in apbs/tools/conversion/param/vparam.\n");
        return VRC_FAILURE;
    }

    if (apolparm->bconc != 0.0){
        atomData = Vparam_getAtomData(param, "WAT", "OW");
        if (atomData == VNULL) atomData = Vparam_getAtomData(param, "WAT", "O");
        if (atomData == VNULL) {
            Vnm_tprint(2, "initAPOL:  Couldn't find parameters for WAT OW or WAT O!\n");
            Vnm_tprint(2, "initAPOL:  These parameters must be present in your file\n");
            Vnm_tprint(2, "initAPOL:  for apolar calculations.\n");
            return VRC_FAILURE;
        }
        apolparm->watepsilon = atomData->epsilon;
        apolparm->watsigma = atomData->radius;
        apolparm->setwat = 1;
    }

    /* Calculate Energy and Forces */
    if(apolparm->calcforce) {
        rc = forceAPOL(acc, mem, apolparm, nforce, atomForce, alist, clist);
        if(rc == VRC_FAILURE) {
            Vnm_print(2, "Error in apolar force calculation!\n");
            return VRC_FAILURE;
        }
    }

    /* Get the SAV and SAS */
    sasa = 0.0;
    sav = 0.0;

    if (apolparm->calcenergy) {
        len = Valist_getNumberAtoms(alist);

        if (VABS(apolparm->gamma) > VSMALL) {
            /* Total Solvent Accessible Surface Area (SASA) */
            apolparm->sasa = Vacc_totalSASA(acc, srad);
            /* SASA for each atom */
            for (i = 0; i < len; i++) {
                atom = Valist_getAtom(alist, i);
                atomsasa[i] = Vacc_atomSASA(acc, srad, atom);
            }
        } else {
            /* Total Solvent Accessible Surface Area (SASA) set to zero */
            apolparm->sasa = 0.0;
            /* SASA for each atom set to zero*/
            for (i = 0; i < len; i++) {
                atomsasa[i] = 0.0;
            }
        }

        /* Inflated van der Waals accessibility */
        if (VABS(apolparm->press) > VSMALL){
            apolparm->sav = Vacc_totalSAV(acc, clist, apolparm, srad);
        } else {
            apolparm->sav = 0.0;
        }

        /* wcaEnergy integral code */
        if (VABS(apolparm->bconc) > VSMALL) {
            /* wcaEnergy for each atom */
            for (i = 0; i < len; i++) {
                rc = Vacc_wcaEnergyAtom(acc, apolparm, alist, clist, i, &energy);
                if (rc == 0)  {
                    Vnm_print(2, "Error in apolar energy calculation!\n");
                    return 0;
                }
                atomwcaEnergy[i] = energy;
            }
            /* Total WCA Energy */
            rc = Vacc_wcaEnergy(acc, apolparm, alist, clist);
            if (rc == 0) {
                Vnm_print(2, "Error in apolar energy calculation!\n");
                return 0;
            }
        } else {
            apolparm->wcaEnergy = 0.0;
        }
        energyAPOL(apolparm, apolparm->sasa, apolparm->sav, atomsasa, atomwcaEnergy, Valist_getNumberAtoms(alist));
    }

    Vmem_free(VNULL, Valist_getNumberAtoms(alist), sizeof(double), (void **)&(atomsasa));
    Vmem_free(VNULL, Valist_getNumberAtoms(alist), sizeof(double), (void **)&(atomwcaEnergy));
    Vclist_dtor(&clist);
    Vacc_dtor(&acc);

    return VRC_SUCCESS;
}

VPUBLIC int energyAPOL(APOLparm *apolparm,
                       double sasa,
                       double sav,
                       double atomsasa[],
                       double atomwcaEnergy[],
                       int numatoms
                      ){

    double energy = 0.0;
    int i = 0;

#ifndef VAPBSQUIET
    Vnm_print(1,"\nSolvent Accessible Surface Area (SASA) for each atom:\n");
    for (i = 0; i < numatoms; i++) {
        Vnm_print(1,"  SASA for atom %i: %1.12E\n", i, atomsasa[i]);
    }

    Vnm_print(1,"\nTotal solvent accessible surface area: %g A^2\n",sasa);
#endif

    switch(apolparm->calcenergy){
        case ACE_NO:
            break;
        case ACE_COMPS:
            Vnm_print(1,"energyAPOL: Cannot calculate component energy, skipping.\n");
            break;
        case ACE_TOTAL:
            energy = (apolparm->gamma*sasa) + (apolparm->press*sav)
                        + (apolparm->wcaEnergy);
#ifndef VAPBSQUIET
            Vnm_print(1,"\nSurface tension*area energies (gamma * SASA) for each atom:\n");
            for (i = 0; i < numatoms; i++) {
                Vnm_print(1,"  Surface tension*area energy for atom %i: %1.12E\n", i, apolparm->gamma*atomsasa[i]);
            }

            Vnm_print(1,"\nTotal surface tension energy: %g kJ/mol\n", apolparm->gamma*sasa);
            Vnm_print(1,"\nTotal solvent accessible volume: %g A^3\n", sav);
            Vnm_print(1,"\nTotal pressure*volume energy: %g kJ/mol\n", apolparm->press*sav);
            Vnm_print(1,"\nWCA dispersion Energies for each atom:\n");
            for (i = 0; i < numatoms; i++) {
                Vnm_print(1,"  WCA energy for atom %i: %1.12E\n", i, atomwcaEnergy[i]);
            }

            Vnm_print(1,"\nTotal WCA energy: %g kJ/mol\n", (apolparm->wcaEnergy));
            Vnm_print(1,"\nTotal non-polar energy = %1.12E kJ/mol\n", energy);
#endif
            break;
        default:
            Vnm_print(2,"energyAPOL: Error in energyAPOL. Unknown option.\n");
            break;
    }

    return VRC_SUCCESS;
}

VPUBLIC int forceAPOL(Vacc *acc,
                      Vmem *mem,
                      APOLparm *apolparm,
                      int *nforce,
                      AtomForce **atomForce,
                      Valist *alist,
                      Vclist *clist
                     ) {
                         time_t ts, ts_main, ts_sub;
    int i,
        j,
        natom;

    double srad, /* Probe radius */
           xF,
           yF,
           zF,	/* Individual forces */
           press,
           gamma,
           offset,
           bconc,
           dSASA[3],
           dSAV[3],
           force[3],
           *apos;

    Vatom *atom = VNULL;
    ts_main = clock();

    srad = apolparm->srad;
    press = apolparm->press;
    gamma = apolparm->gamma;
    offset = apolparm->dpos;
    bconc = apolparm->bconc;

    natom = Valist_getNumberAtoms(alist);

    /* Check to see if we need to build the surface */
    Vnm_print(0, "forceAPOL: Trying atom surf...\n");
    ts = clock();
    if (acc->surf == VNULL) {
        acc->surf = (VaccSurf**)Vmem_malloc(acc->mem, natom, sizeof(VaccSurf *));
        for (i=0; i<natom; i++) {
            atom = Valist_getAtom(acc->alist, i);
            //apos = Vatom_getPosition(atom); // apos never referenced? - Peter
            /* NOTE:  RIGHT NOW WE DO THIS FOR THE ENTIRE MOLECULE WHICH IS
             * INCREDIBLY INEFFICIENT, PARTICULARLY DURING FOCUSING!!! */
            acc->surf[i] = Vacc_atomSurf(acc, atom, acc->refSphere, srad);
        }
    }
    Vnm_print(0, "forceAPOL: atom surf: Time elapsed: %f\n", ((double)clock() - ts) / CLOCKS_PER_SEC);

    if(apolparm->calcforce == ACF_TOTAL){
        Vnm_print(0, "forceAPOL: calcforce == ACF_TOTAL\n");
        ts = clock();

        *nforce = 1;
        if(*atomForce != VNULL){
            Vmem_free(mem,*nforce,sizeof(AtomForce), (void **)atomForce);
        }

            *atomForce = (AtomForce *)Vmem_malloc(mem, *nforce,
                                                  sizeof(AtomForce));

        /* Clear out force arrays */
        for (j=0; j<3; j++) {
            (*atomForce)[0].sasaForce[j] = 0.0;
            (*atomForce)[0].savForce[j] = 0.0;
            (*atomForce)[0].wcaForce[j] = 0.0;
        }

        // problem block
        for (i=0; i<natom; i++) {
            atom = Valist_getAtom(alist, i);

            for(j=0;j<3;j++){
                dSASA[j] = 0.0;
                dSAV[j] = 0.0;
                force[j] = 0.0;
            }

            if(VABS(gamma) > VSMALL) {
                Vacc_atomdSASA(acc, offset, srad, atom, dSASA);
            }
            if(VABS(press) > VSMALL) {
                Vacc_atomdSAV(acc, srad, atom, dSAV);
            }
            if(VABS(bconc) > VSMALL) {
                Vacc_wcaForceAtom(acc, apolparm, clist, atom, force);
            }

            for(j=0;j<3;j++){
                (*atomForce)[0].sasaForce[j] += dSASA[j];
                (*atomForce)[0].savForce[j] += dSAV[j];
                (*atomForce)[0].wcaForce[j] += force[j];
            }
        }
        // end block

        Vnm_tprint( 1, "  Printing net forces (kJ/mol/A)\n");
        Vnm_tprint( 1, "  Legend:\n");
        Vnm_tprint( 1, "    sasa  -- SASA force\n");
        Vnm_tprint( 1, "    sav   -- SAV force\n");
        Vnm_tprint( 1, "    wca   -- WCA force\n\n");

        Vnm_tprint( 1, "  sasa  %4.3e %4.3e %4.3e\n",
                    (*atomForce)[0].sasaForce[0],
                    (*atomForce)[0].sasaForce[1],
                    (*atomForce)[0].sasaForce[2]);
        Vnm_tprint( 1, "  sav   %4.3e %4.3e %4.3e\n",
                    (*atomForce)[0].savForce[0],
                    (*atomForce)[0].savForce[1],
                    (*atomForce)[0].savForce[2]);
        Vnm_tprint( 1, "  wca   %4.3e %4.3e %4.3e\n",
                    (*atomForce)[0].wcaForce[0],
                    (*atomForce)[0].wcaForce[1],
                    (*atomForce)[0].wcaForce[2]);

        Vnm_print(0, "forceAPOL: calcforce == ACF_TOTAL: %f\n", ((double)clock() - ts) / CLOCKS_PER_SEC);
    } else if (apolparm->calcforce == ACF_COMPS ){
        *nforce = Valist_getNumberAtoms(alist);
        if(*atomForce == VNULL){
            *atomForce = (AtomForce *)Vmem_malloc(mem, *nforce,
                                                  sizeof(AtomForce));
        }else{
            Vmem_free(mem,*nforce,sizeof(AtomForce), (void **)atomForce);
            *atomForce = (AtomForce *)Vmem_malloc(mem, *nforce,
                                                  sizeof(AtomForce));
        }

#ifndef VAPBSQUIET
        Vnm_tprint( 1, "  Printing per atom forces (kJ/mol/A)\n");
        Vnm_tprint( 1, "  Legend:\n");
        Vnm_tprint( 1, "    tot  n -- Total force for atom n\n");
        Vnm_tprint( 1, "    sasa n -- SASA force for atom n\n");
        Vnm_tprint( 1, "    sav  n -- SAV force for atom n\n");
        Vnm_tprint( 1, "    wca  n -- WCA force for atom n\n\n");

        Vnm_tprint( 1, "    gamma    %f\n" \
                       "    pressure %f\n" \
                       "    bconc    %f \n\n",
                            gamma,press,bconc);
#endif

        for (i=0; i<natom; i++) {
            atom = Valist_getAtom(alist, i);

            for(j=0;j<3;j++){
                dSASA[j] = 0.0;
                dSAV[j] = 0.0;
                force[j] = 0.0;
            }

            /* Clear out force arrays */
            for (j=0; j<3; j++) {
                (*atomForce)[i].sasaForce[j] = 0.0;
                (*atomForce)[i].savForce[j] = 0.0;
                (*atomForce)[i].wcaForce[j] = 0.0;
            }

            if(VABS(gamma) > VSMALL) Vacc_atomdSASA(acc, offset, srad, atom, dSASA);
            if(VABS(press) > VSMALL) Vacc_atomdSAV(acc, srad, atom, dSAV);
            if(VABS(bconc) > VSMALL) Vacc_wcaForceAtom(acc,apolparm,clist,atom,force);

            xF = -((gamma*dSASA[0]) + (press*dSAV[0]) + (bconc*force[0]));
            yF = -((gamma*dSASA[1]) + (press*dSAV[1]) + (bconc*force[1]));
            zF = -((gamma*dSASA[2]) + (press*dSAV[2]) + (bconc*force[2]));

            for(j=0;j<3;j++){
                (*atomForce)[i].sasaForce[j] += dSASA[j];
                (*atomForce)[i].savForce[j] += dSAV[j];
                (*atomForce)[i].wcaForce[j] += force[j];
            }

#ifndef VAPBSQUIET
            Vnm_print( 1, "  tot  %i %4.3e %4.3e %4.3e\n",
                        i,
                        xF,
                        yF,
                        zF);
            Vnm_print( 1, "  sasa %i %4.3e %4.3e %4.3e\n",
                        i,
                        (*atomForce)[i].sasaForce[0],
                        (*atomForce)[i].sasaForce[1],
                        (*atomForce)[i].sasaForce[2]);
            Vnm_print( 1, "  sav  %i %4.3e %4.3e %4.3e\n",
                        i,
                        (*atomForce)[i].savForce[0],
                        (*atomForce)[i].savForce[1],
                        (*atomForce)[i].savForce[2]);
            Vnm_print( 1, "  wca  %i %4.3e %4.3e %4.3e\n",
                        i,
                        (*atomForce)[i].wcaForce[0],
                        (*atomForce)[i].wcaForce[1],
                        (*atomForce)[i].wcaForce[2]);
#endif

        }
    } else *nforce = 0;

#ifndef VAPBSQUIET
    Vnm_print(1,"\n");
#endif

    Vnm_print(0, "forceAPOL: Time elapsed: %f\n", ((double)clock() - ts_main) / CLOCKS_PER_SEC);
    return VRC_SUCCESS;
}


