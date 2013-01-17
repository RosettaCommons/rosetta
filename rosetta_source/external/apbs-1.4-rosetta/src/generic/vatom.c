/**
 *  @file    vatom.c
 *  @ingroup Vatom
 *  @author  Nathan Baker
 *  @brief   Class Vatom methods
 *  @version $Id$
 *  @attention
 *  @verbatim
 *
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
 *
 * @endverbatim
 */

#include "vatom.h"

VEMBED(rcsid="$Id$")

#if !defined(VINLINE_VATOM)

VPUBLIC double *Vatom_getPosition(Vatom *thee) {

   VASSERT(thee != VNULL);
   return thee->position;

}

VPUBLIC double Vatom_getPartID(Vatom *thee) {

   VASSERT(thee != VNULL);
   return thee->partID;

}

VPUBLIC void Vatom_setPartID(Vatom *thee, int partID) {

   VASSERT(thee != VNULL);
   thee->partID = (double)partID;

}

VPUBLIC double Vatom_getAtomID(Vatom *thee) {

   VASSERT(thee != VNULL);
   return thee->id;

}

VPUBLIC void Vatom_setAtomID(Vatom *thee, int atomID) {

   VASSERT(thee != VNULL);
   thee->id = atomID;

}

VPUBLIC void Vatom_setRadius(Vatom *thee, double radius) {

   VASSERT(thee != VNULL);
   thee->radius = radius;

}

VPUBLIC double Vatom_getRadius(Vatom *thee) {

   VASSERT(thee != VNULL);
   return thee->radius;

}

VPUBLIC void Vatom_setCharge(Vatom *thee, double charge) {

   VASSERT(thee != VNULL);
   thee->charge = charge;

}

VPUBLIC double Vatom_getCharge(Vatom *thee) {

   VASSERT(thee != VNULL);
   return thee->charge;

}

VPUBLIC void Vatom_setEpsilon(Vatom *thee, double epsilon) {

   VASSERT(thee != VNULL);
   thee->epsilon = epsilon;
}

VPUBLIC double Vatom_getEpsilon(Vatom *thee) {

   VASSERT(thee != VNULL);
   return thee->epsilon;
}

VPUBLIC unsigned long int Vatom_memChk(Vatom *thee) { return sizeof(Vatom); }

#endif /* if !defined(VINLINE_VATOM) */

VPUBLIC Vatom* Vatom_ctor() {

    /* Set up the structure */
    Vatom *thee = VNULL;
    thee = (Vatom *)Vmem_malloc( VNULL, 1, sizeof(Vatom) );
    VASSERT( thee != VNULL);
    VASSERT( Vatom_ctor2(thee));

    return thee;
}

VPUBLIC int Vatom_ctor2(Vatom *thee) {
    thee->partID = -1;
    return 1;
}

VPUBLIC void Vatom_dtor(Vatom **thee) {
    if ((*thee) != VNULL) {
        Vatom_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vatom), (void **)thee);
        (*thee) = VNULL;
    }
}

VPUBLIC void Vatom_dtor2(Vatom *thee) { ; }

VPUBLIC void Vatom_setPosition(Vatom *thee, double position[3]) {

   VASSERT(thee != VNULL);
   (thee->position)[0] = position[0];
   (thee->position)[1] = position[1];
   (thee->position)[2] = position[2];

}

VPUBLIC void Vatom_copyTo(Vatom *thee, Vatom *dest) {

    VASSERT(thee != VNULL);
    VASSERT(dest != VNULL);

    memcpy(dest, thee, sizeof(Vatom));

}

VPUBLIC void Vatom_copyFrom(Vatom *thee, Vatom *src) {

    Vatom_copyTo(src, thee);

}

VPUBLIC void Vatom_setResName(Vatom *thee, char resName[VMAX_RECLEN]) {

    VASSERT(thee != VNULL);
    strcpy(thee->resName, resName);

}

VPUBLIC void Vatom_getResName(Vatom *thee, char resName[VMAX_RECLEN]) {


    VASSERT(thee != VNULL);
    strcpy(resName,thee->resName);

}

VPUBLIC void Vatom_setAtomName(Vatom *thee, char atomName[VMAX_RECLEN]) {

    VASSERT(thee != VNULL);
    strcpy(thee->atomName, atomName);

}

VPUBLIC void Vatom_getAtomName(Vatom *thee, char atomName[VMAX_RECLEN]) {

    VASSERT(thee != VNULL);
    strcpy(atomName,thee->atomName);

}

#if defined(WITH_TINKER)

VPUBLIC void Vatom_setDipole(Vatom *thee, double dipole[3]) {

   VASSERT(thee != VNULL);
   (thee->dipole)[0] = dipole[0];
   (thee->dipole)[1] = dipole[1];
   (thee->dipole)[2] = dipole[2];

}

VPUBLIC void Vatom_setQuadrupole(Vatom *thee, double quadrupole[9]) {

   int i;
   VASSERT(thee != VNULL);
   for (i=0; i<9; i++)  (thee->quadrupole)[i] = quadrupole[i];
}

VPUBLIC void Vatom_setInducedDipole(Vatom *thee, double dipole[3]) {

   VASSERT(thee != VNULL);
   (thee->inducedDipole)[0] = dipole[0];
   (thee->inducedDipole)[1] = dipole[1];
   (thee->inducedDipole)[2] = dipole[2];
}

VPUBLIC void Vatom_setNLInducedDipole(Vatom *thee, double dipole[3]) {

   VASSERT(thee != VNULL);
   (thee->nlInducedDipole)[0] = dipole[0];
   (thee->nlInducedDipole)[1] = dipole[1];
   (thee->nlInducedDipole)[2] = dipole[2];

}

VPUBLIC double *Vatom_getDipole(Vatom *thee) {

   VASSERT(thee != VNULL);
   return thee->dipole;

}

VPUBLIC double *Vatom_getQuadrupole(Vatom *thee) {

   VASSERT(thee != VNULL);
   return thee->quadrupole;

}

VPUBLIC double *Vatom_getInducedDipole(Vatom *thee) {

   VASSERT(thee != VNULL);
   return thee->inducedDipole;

}

VPUBLIC double *Vatom_getNLInducedDipole(Vatom *thee) {

   VASSERT(thee != VNULL);
   return thee->nlInducedDipole;

}

#endif /* if defined(WITH_TINKER) */
