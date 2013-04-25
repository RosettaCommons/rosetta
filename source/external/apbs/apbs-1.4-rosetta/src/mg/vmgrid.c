/**
 *  @file    vmgrid.c
 *  @author  Nathan Baker
 *  @brief   Class Vmgrid methods
 *  @ingroup Vmgrid
 *  @version $Id: vmgrid.c 1615 2010-10-20 19:16:35Z sobolevnrm $
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
 * Copyright (c) 2010, Pacific Northwest National Laboratory.  Portions Copyright (c) 2002-2010, Washington University in St. Louis.  Portions Copyright (c) 2002-2010, Nathan A. Baker.  Portions Copyright (c) 1999-2002, The Regents of the University of California. Portions Copyright (c) 1995, Michael Holst.
 * All rights reserved.
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

#include "vmgrid.h"

VEMBED(rcsid="$Id: vmgrid.c 1615 2010-10-20 19:16:35Z sobolevnrm $")

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vmgrid_ctor
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC Vmgrid* Vmgrid_ctor() {

    Vmgrid *thee = VNULL;

    thee = Vmem_malloc(VNULL, 1, sizeof(Vmgrid));
    VASSERT(thee != VNULL);
    VASSERT(Vmgrid_ctor2(thee));

    return thee;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vmgrid_ctor2
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vmgrid_ctor2(Vmgrid *thee) {

    int i;

    if (thee == VNULL) return 0;

    thee->ngrids = 0;
    for (i=0; i<VMGRIDMAX; i++) thee->grids[i] = VNULL;

    return 1;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vmgrid_dtor
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vmgrid_dtor(Vmgrid **thee) {

    if ((*thee) != VNULL) {
        Vmgrid_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vmgrid), (void **)thee);
        (*thee) = VNULL;
    }
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vmgrid_dtor2
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC void Vmgrid_dtor2(Vmgrid *thee) { ; }

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vmgrid_value
// Author:   Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vmgrid_value(Vmgrid *thee, double pt[3], double *value) {

    int i, rc;
    double tvalue;

    VASSERT(thee != VNULL);

    for (i=0; i<thee->ngrids; i++) {
        rc = Vgrid_value(thee->grids[i], pt, &tvalue);
        if (rc) {
            *value = tvalue;
            return 1;
        }
    }

    Vnm_print(2, "Vmgrid_value:  Point (%g, %g, %g) not found in \
hiearchy!\n", pt[0], pt[1], pt[2]);

    return 0;
}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vmgrid_curvature
//
//   Notes:  cflag=0 ==> Reduced Maximal Curvature
//           cflag=1 ==> Mean Curvature (Laplace)
//           cflag=2 ==> Gauss Curvature
//           cflag=3 ==> True Maximal Curvature
//
// Authors:  Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vmgrid_curvature(Vmgrid *thee, double pt[3], int cflag,
  double *value) {

    int i, rc;
    double tvalue;

    VASSERT(thee != VNULL);

    for (i=0; i<thee->ngrids; i++) {
        rc = Vgrid_curvature(thee->grids[i], pt, cflag, &tvalue);
        if (rc) {
            *value = tvalue;
            return 1;
        }
    }

    Vnm_print(2, "Vmgrid_curvature:  Point (%g, %g, %g) not found in \
hiearchy!\n", pt[0], pt[1], pt[2]);

    return 0;


}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vmgrid_gradient
//
// Authors:  Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vmgrid_gradient(Vmgrid *thee, double pt[3], double grad[3]) {

    int i, j, rc;
    double tgrad[3];

    VASSERT(thee != VNULL);

    for (i=0; i<thee->ngrids; i++) {
        rc = Vgrid_gradient(thee->grids[i], pt, tgrad);
        if (rc) {
            for (j=0; j<3; j++) grad[j] = tgrad[j];
            return 1;
        }
    }

    Vnm_print(2, "Vmgrid_gradient:  Point (%g, %g, %g) not found in \
hiearchy!\n", pt[0], pt[1], pt[2]);

    return 0;


}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vmgrid_addGrid
//
// Authors:  Nathan Baker
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vmgrid_addGrid(Vmgrid *thee, Vgrid *grid) {

    int i, j, rc;
    double tgrad[3];

    VASSERT(thee != VNULL);

    if (grid == VNULL) {
        Vnm_print(2, "Vmgrid_addGrid:  Not adding VNULL grid!\n");
        return 0;
    }

    if (thee->ngrids >= VMGRIDMAX) {
        Vnm_print(2, "Vmgrid_addGrid:  Too many grids in hierarchy (max = \
%d)!\n", VMGRIDMAX);
        Vnm_print(2, "Vmgrid_addGrid:  Not adding grid!\n");
        return 0;
    }

    thee->grids[thee->ngrids] = grid;
    (thee->ngrids)++;

    return 1;

}
