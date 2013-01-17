/**
 *  @file    pbeparm.c
 *  @ingroup PBEparm
 *  @author  Nathan Baker
 *  @brief   Class PBEparm methods
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

#include "pbeparm.h"

VEMBED(rcsid="$Id$")

#if !defined(VINLINE_MGPARM)

#endif /* if !defined(VINLINE_MGPARM) */

VPUBLIC double PBEparm_getIonCharge(PBEparm *thee, int i) {
    VASSERT(thee != VNULL);
    VASSERT(i < thee->nion);
    return thee->ionq[i];
}

VPUBLIC double PBEparm_getIonConc(PBEparm *thee, int i) {
    VASSERT(thee != VNULL);
    VASSERT(i < thee->nion);
    return thee->ionc[i];
}

VPUBLIC double PBEparm_getIonRadius(PBEparm *thee, int i) {
    VASSERT(thee != VNULL);
    VASSERT(i < thee->nion);
    return thee->ionr[i];
}

/*----------------------------------------------------------------------*/
/* Added by Michael Grabe                                               */
/*----------------------------------------------------------------------*/

VPUBLIC double PBEparm_getzmem(PBEparm *thee) {
    VASSERT(thee != VNULL);
    return thee->zmem;
}
VPUBLIC double PBEparm_getLmem(PBEparm *thee) {
    VASSERT(thee != VNULL);
    return thee->Lmem;
}
VPUBLIC double PBEparm_getmembraneDiel(PBEparm *thee) {
    VASSERT(thee != VNULL);
    return thee->mdie;
}
VPUBLIC double PBEparm_getmemv(PBEparm *thee) {
    VASSERT(thee != VNULL);
    return thee->memv;
}

VPUBLIC PBEparm* PBEparm_ctor() {

    /* Set up the structure */
    PBEparm *thee = VNULL;
    thee = (PBEparm*)Vmem_malloc(VNULL, 1, sizeof(PBEparm));
    VASSERT( thee != VNULL);
    VASSERT( PBEparm_ctor2(thee) );

    return thee;
}

VPUBLIC int PBEparm_ctor2(PBEparm *thee) {

    int i;

    if (thee == VNULL) return 0;

    thee->parsed = 0;

    thee->setmolid = 0;
    thee->setpbetype = 0;
    thee->setbcfl = 0;
    thee->setnion = 0;
    for (i=0; i<MAXION; i++){
        thee->setion[i] = 0;
        thee->ionq[i] = 0.0;
        thee->ionc[i] = 0.0;
        thee->ionr[i] = 0.0;
    }
    thee->setpdie = 0;
    thee->setsdie = 0;
    thee->setsrfm = 0;
    thee->setsrad = 0;
    thee->setswin = 0;
    thee->settemp = 0;
    thee->setcalcenergy = 0;
    thee->setcalcforce = 0;
    thee->setsdens = 0;
    thee->numwrite = 0;
    thee->setwritemat = 0;
    thee->nion = 0;
    thee->sdens = 0;
    thee->swin = 0;
    thee->srad = 1.4;
    thee->useDielMap = 0;
    thee->useKappaMap = 0;
    thee->usePotMap = 0;
    thee->useChargeMap = 0;

    /*----------------------------------------------*/
    /* Added by Michael Grabe                       */
    /*----------------------------------------------*/

    thee->setzmem = 0;
    thee->setLmem = 0;
    thee->setmdie = 0;
    thee->setmemv = 0;

    /*----------------------------------------------*/

    thee->smsize = 0.0;
    thee->smvolume = 0.0;

    thee->setsmsize = 0;
    thee->setsmvolume = 0;

    return 1;
}

VPUBLIC void PBEparm_dtor(PBEparm **thee) {
    if ((*thee) != VNULL) {
        PBEparm_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(PBEparm), (void **)thee);
        (*thee) = VNULL;
    }
}

VPUBLIC void PBEparm_dtor2(PBEparm *thee) { ; }

VPUBLIC int PBEparm_check(PBEparm *thee) {

    int i;

    /* Check to see if we were even filled... */
    if (!thee->parsed) {
        Vnm_print(2, "PBEparm_check:  not filled!\n");
        return 0;
    }

    if (!thee->setmolid) {
        Vnm_print(2, "PBEparm_check:  MOL not set!\n");
        return 0;
    }
    if (!thee->setpbetype) {
        Vnm_print(2, "PBEparm_check:  LPBE/NPBE/LRPBE/NRPBE/SMPBE not set!\n");
        return 0;
    }
    if (!thee->setbcfl) {
        Vnm_print(2, "PBEparm_check:  BCFL not set!\n");
        return 0;
    }
    if (!thee->setnion) {
        thee->setnion = 1;
        thee->nion = 0;
    }
    for (i=0; i<thee->nion; i++) {
        if (!thee->setion[i]) {
            Vnm_print(2, "PBEparm_check: ION #%d not set!\n",i);
            return 0;
        }
    }
    if (!thee->setpdie) {
        Vnm_print(2, "PBEparm_check: PDIE not set!\n");
        return 0;
    }
    if (((thee->srfm==VSM_MOL) || (thee->srfm==VSM_MOLSMOOTH)) \
      && (!thee->setsdens) && (thee->srad > VSMALL)) {
        Vnm_print(2, "PBEparm_check: SDENS not set!\n");
        return 0;
    }
    if (!thee->setsdie) {
        Vnm_print(2, "PBEparm_check: SDIE not set!\n");
        return 0;
    }
    if (!thee->setsrfm) {
        Vnm_print(2, "PBEparm_check: SRFM not set!\n");
        return 0;
    }
    if (((thee->srfm==VSM_MOL) || (thee->srfm==VSM_MOLSMOOTH)) \
      && (!thee->setsrad)) {
        Vnm_print(2, "PBEparm_check: SRAD not set!\n");
        return 0;
    }
    if ((thee->srfm==VSM_SPLINE) && (!thee->setswin)) {
        Vnm_print(2, "PBEparm_check: SWIN not set!\n");
        return 0;
    }
    if ((thee->srfm==VSM_SPLINE3) && (!thee->setswin)) {
        Vnm_print(2, "PBEparm_check: SWIN not set!\n");
        return 0;
    }
    if ((thee->srfm==VSM_SPLINE4) && (!thee->setswin)) {
        Vnm_print(2, "PBEparm_check: SWIN not set!\n");
        return 0;
    }
    if (!thee->settemp) {
        Vnm_print(2, "PBEparm_check: TEMP not set!\n");
        return 0;
    }
    if (!thee->setcalcenergy) thee->calcenergy = PCE_NO;
    if (!thee->setcalcforce) thee->calcforce = PCF_NO;
    if (!thee->setwritemat) thee->writemat = 0;

    /*--------------------------------------------------------*/
    /* Added by Michael Grabe                                 */
    /*--------------------------------------------------------*/

    if ((!thee->setzmem) && (thee->bcfl == 3)){
        Vnm_print(2, "PBEparm_check: ZMEM not set!\n");
        return 0;
    }
    if ((!thee->setLmem) && (thee->bcfl == 3)){
        Vnm_print(2, "PBEparm_check: LMEM not set!\n");
        return 0;
    }
    if ((!thee->setmdie) && (thee->bcfl == 3)){
        Vnm_print(2, "PBEparm_check: MDIE not set!\n");
        return 0;
    }
    if ((!thee->setmemv) && (thee->bcfl == 3)){
        Vnm_print(2, "PBEparm_check: MEMV not set!\n");
        return 0;
    }

    /*--------------------------------------------------------*/

    return 1;
}

VPUBLIC void PBEparm_copy(PBEparm *thee, PBEparm *parm) {

    int i, j;

    VASSERT(thee != VNULL);
    VASSERT(parm != VNULL);

    thee->molid = parm->molid;
    thee->setmolid = parm->setmolid;
    thee->useDielMap = parm->useDielMap;
    thee->dielMapID = parm->dielMapID;
    thee->useKappaMap = parm->useKappaMap;
    thee->kappaMapID = parm->kappaMapID;
    thee->usePotMap = parm->usePotMap;
    thee->potMapID = parm->potMapID;
    thee->useChargeMap = parm->useChargeMap;
    thee->chargeMapID = parm->chargeMapID;
    thee->pbetype = parm->pbetype;
    thee->setpbetype = parm->setpbetype;
    thee->bcfl = parm->bcfl;
    thee->setbcfl = parm->setbcfl;
    thee->nion = parm->nion;
    thee->setnion = parm->setnion;
    for (i=0; i<MAXION; i++) {
        thee->ionq[i] = parm->ionq[i];
        thee->ionc[i] = parm->ionc[i];
        thee->ionr[i] = parm->ionr[i];
        thee->setion[i] = parm->setion[i];
    };
    thee->pdie = parm->pdie;
    thee->setpdie = parm->setpdie;
    thee->sdens = parm->sdens;
    thee->setsdens = parm->setsdens;
    thee->sdie = parm->sdie;
    thee->setsdie = parm->setsdie;
    thee->srfm = parm->srfm;
    thee->setsrfm = parm->setsrfm;
    thee->srad = parm->srad;
    thee->setsrad = parm->setsrad;
    thee->swin = parm->swin;
    thee->setswin = parm->setswin;
    thee->temp = parm->temp;
    thee->settemp = parm->settemp;
    thee->calcenergy = parm->calcenergy;
    thee->setcalcenergy = parm->setcalcenergy;
    thee->calcforce = parm->calcforce;
    thee->setcalcforce = parm->setcalcforce;

    /*----------------------------------------------------*/
    /* Added by Michael Grabe                             */
    /*----------------------------------------------------*/

    thee->zmem = parm->zmem;
    thee->setzmem = parm->setzmem;
    thee->Lmem = parm->Lmem;
    thee->setLmem = parm->setLmem;
    thee->mdie = parm->mdie;
    thee->setmdie = parm->setmdie;
    thee->memv = parm->memv;
    thee->setmemv = parm->setmemv;

    /*----------------------------------------------------*/

    thee->numwrite = parm->numwrite;
    for (i=0; i<PBEPARM_MAXWRITE; i++) {
        thee->writetype[i] = parm->writetype[i];
        thee->writefmt[i] = parm->writefmt[i];
        for (j=0; j<VMAX_ARGLEN; j++)
          thee->writestem[i][j] = parm->writestem[i][j];
    }
    thee->writemat = parm->writemat;
    thee->setwritemat = parm->setwritemat;
    for (i=0; i<VMAX_ARGLEN; i++) thee->writematstem[i] = parm->writematstem[i];
    thee->writematflag = parm->writematflag;

    thee->smsize = parm->smsize;
    thee->smvolume = parm->smvolume;

    thee->setsmsize = parm->setsmsize;
    thee->setsmvolume = parm->setsmvolume;

    thee->parsed = parm->parsed;

}

VPRIVATE int PBEparm_parseLPBE(PBEparm *thee, Vio *sock) {
    Vnm_print(0, "NOsh: parsed lpbe\n");
    thee->pbetype = PBE_LPBE;
    thee->setpbetype = 1;
    return 1;
}

VPRIVATE int PBEparm_parseNPBE(PBEparm *thee, Vio *sock) {
    Vnm_print(0, "NOsh: parsed npbe\n");
    thee->pbetype = PBE_NPBE;
    thee->setpbetype = 1;
    return 1;
}

VPRIVATE int PBEparm_parseMOL(PBEparm *thee, Vio *sock) {
    int ti;
    char tok[VMAX_BUFSIZE];

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%d", &ti) == 0) {
        Vnm_print(2, "NOsh:  Read non-int (%s) while parsing MOL \
keyword!\n", tok);
        return -1;
    }
    thee->molid = ti;
    thee->setmolid = 1;
    return 1;

    VERROR1:
        Vnm_print(2, "parsePBE:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int PBEparm_parseLRPBE(PBEparm *thee, Vio *sock) {
    Vnm_print(0, "NOsh: parsed lrpbe\n");
    thee->pbetype = PBE_LRPBE;
    thee->setpbetype = 1;
    return 1;
}

VPRIVATE int PBEparm_parseNRPBE(PBEparm *thee, Vio *sock) {
    Vnm_print(0, "NOsh: parsed nrpbe\n");
    thee->pbetype = PBE_NRPBE;
    thee->setpbetype = 1;
    return 1;
}

VPRIVATE int PBEparm_parseSMPBE(PBEparm *thee, Vio *sock) {

    int i;

    char type[VMAX_BUFSIZE]; /* vol or size (keywords) */
    char value[VMAX_BUFSIZE]; /* floating point value */

    char setVol = 1;
    char setSize = 1;
    char keyValuePairs = 2;

    double size, volume;

    for(i=0;i<keyValuePairs;i++){

        /* The line two tokens at a time */
        VJMPERR1(Vio_scanf(sock, "%s", type) == 1);
        VJMPERR1(Vio_scanf(sock, "%s", value) == 1);

        if(!strcmp(type,"vol")){
            if ((setVol = sscanf(value, "%lf", &volume)) == 0){
                Vnm_print(2,"NOsh:  Read non-float (%s) while parsing smpbe keyword!\n", value);
                return VRC_FAILURE;
            }
        }else if(!strcmp(type,"size")){
            if ((setSize = sscanf(value, "%lf", &size)) == 0){
                Vnm_print(2,"NOsh:  Read non-float (%s) while parsing smpbe keyword!\n", value);
                return VRC_FAILURE;
            }
        }else{
            Vnm_print(2,"NOsh:  Read non-float (%s) while parsing smpbe keyword!\n", value);
            return VRC_FAILURE;
        }
    }

    /* If either the volume or size isn't set, throw an error */
    if((setVol == 0) || (setSize == 0)){
        Vnm_print(2,"NOsh:  Error while parsing smpbe keywords! Only size or vol was specified.\n");
        return VRC_FAILURE;
    }

    Vnm_print(0, "NOsh: parsed smpbe\n");
    thee->pbetype = PBE_SMPBE;
    thee->setpbetype = 1;

    thee->smsize = size;
    thee->setsmsize = 1;

    thee->smvolume = volume;
    thee->setsmvolume = 1;

    return VRC_SUCCESS;

VERROR1:
    Vnm_print(2, "parsePBE:  ran out of tokens!\n");
    return VRC_FAILURE;

}

VPRIVATE int PBEparm_parseBCFL(PBEparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    int ti;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);

    /* We can either parse int flag... */
    if (sscanf(tok, "%d", &ti) == 1) {

        thee->bcfl = (Vbcfl)ti;
        thee->setbcfl = 1;
        /* Warn that this usage is deprecated */
        Vnm_print(2, "parsePBE:  Warning -- parsed deprecated \"bcfl %d\" \
statement\n", ti);
        Vnm_print(2, "parsePBE:  Please use \"bcfl ");
        switch (thee->bcfl) {
            case BCFL_ZERO:
                Vnm_print(2, "zero");
                break;
            case BCFL_SDH:
                Vnm_print(2, "sdh");
                break;
            case BCFL_MDH:
                Vnm_print(2, "mdh");
                break;
            case BCFL_FOCUS:
                Vnm_print(2, "focus");
                break;
            case BCFL_MEM:
                Vnm_print(2, "mem");
                break;
            case BCFL_MAP:
                Vnm_print(2, "map");
                break;
            default:
                Vnm_print(2, "UKNOWN");
                break;
        }
        Vnm_print(2, "\" instead.\n");
        return 1;

    /* ...or the word */
    } else {

        if (Vstring_strcasecmp(tok, "zero") == 0) {
            thee->bcfl = BCFL_ZERO;
            thee->setbcfl = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "sdh") == 0) {
            thee->bcfl = BCFL_SDH;
            thee->setbcfl = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "mdh") == 0) {
            thee->bcfl = BCFL_MDH;
            thee->setbcfl = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "focus") == 0) {
            thee->bcfl = BCFL_FOCUS;
            thee->setbcfl = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "mem") == 0) {
            thee->bcfl = BCFL_MEM;
            thee->setbcfl = 1;
            return 1;
        } else if (Vstring_strcasecmp(tok, "map") == 0) {
            thee->bcfl = BCFL_MAP;
            thee->setbcfl = 1;
            return 1;
        } else {
            Vnm_print(2, "NOsh:  parsed unknown BCFL parameter (%s)!\n",
              tok);
            return -1;
        }
    }
    return 0;

    VERROR1:
        Vnm_print(2, "parsePBE:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int PBEparm_parseION(PBEparm *thee, Vio *sock) {

    int i;
    int meth = 0;

    char tok[VMAX_BUFSIZE];
    char value[VMAX_BUFSIZE];

    double tf;
    double charge, conc, radius;

    int setCharge = 0;
    int setConc = 0;
    int setRadius = 0;
    int keyValuePairs = 3;

    /* Get the initial token for the ION statement */
    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);

    /* Scan the token once to determine the type (old style or new keyValue pair) */
    meth = sscanf(tok, "%lf", &tf);
    /* If tok is a non-zero float value, we are using the old method */
    if(meth != 0){

        Vnm_print(2, "NOsh:  Deprecated use of ION keyword! Use key-value pairs\n", tok);

        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing ION keyword!\n", tok);
            return VRC_FAILURE;
        }
        thee->ionq[thee->nion] = tf;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing ION keyword!\n", tok);
            return VRC_FAILURE;
        }
        thee->ionc[thee->nion] = tf;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%lf", &tf) == 0) {
            Vnm_print(2, "NOsh:  Read non-float (%s) while parsing ION keyword!\n", tok);
            return VRC_FAILURE;
        }
        thee->ionr[thee->nion] = tf;

    }else{

        /* Three key-value pairs (charge, radius and conc) */
        for(i=0;i<keyValuePairs;i++){

            /* Now scan for the value (float) to be used with the key token parsed
             * above the if-else statement */
            VJMPERR1(Vio_scanf(sock, "%s", value) == 1);
            if(!strcmp(tok,"charge")){
                setCharge = sscanf(value, "%lf", &charge);
                if (setCharge == 0){
                    Vnm_print(2,"NOsh:  Read non-float (%s) while parsing ION %s keyword!\n", value, tok);
                    return VRC_FAILURE;
                }
                thee->ionq[thee->nion] = charge;
            }else if(!strcmp(tok,"radius")){
                setRadius = sscanf(value, "%lf", &radius);
                if (setRadius == 0){
                    Vnm_print(2,"NOsh:  Read non-float (%s) while parsing ION %s keyword!\n", value, tok);
                    return VRC_FAILURE;
                }
                thee->ionr[thee->nion] = radius;
            }else if(!strcmp(tok,"conc")){
                setConc = sscanf(value, "%lf", &conc);
                if (setConc == 0){
                    Vnm_print(2,"NOsh:  Read non-float (%s) while parsing ION %s keyword!\n", value, tok);
                    return VRC_FAILURE;
                }
                thee->ionc[thee->nion] = conc;
            }else{
                Vnm_print(2,"NOsh:  Illegal or missing key-value pair for ION keyword!\n");
                return VRC_FAILURE;
            }

            /* If all three values haven't be set (setValue = 0) then read the next token */
            if((setCharge != 1) || (setConc != 1) || (setRadius != 1)){
                VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
            }

        } /* end for */
    } /* end if */

    /* Finally set the setion, nion and setnion flags and return success */
    thee->setion[thee->nion] = 1;
    (thee->nion)++;
    thee->setnion = 1;
    return VRC_SUCCESS;

    VERROR1:
        Vnm_print(2, "parsePBE:  ran out of tokens!\n");
        return VRC_FAILURE;
}

VPRIVATE int PBEparm_parsePDIE(PBEparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    double tf;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing PDIE \
keyword!\n", tok);
        return -1;
    }
    thee->pdie = tf;
    thee->setpdie = 1;
    return 1;

    VERROR1:
        Vnm_print(2, "parsePBE:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int PBEparm_parseSDENS(PBEparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    double tf;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing SDENS \
keyword!\n", tok);
        return -1;
    }
    thee->sdens = tf;
    thee->setsdens = 1;
    return 1;

    VERROR1:
        Vnm_print(2, "parsePBE:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int PBEparm_parseSDIE(PBEparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    double tf;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing SDIE \
keyword!\n", tok);
        return -1;
    }
    thee->sdie = tf;
    thee->setsdie = 1;
    return 1;

    VERROR1:
        Vnm_print(2, "parsePBE:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int PBEparm_parseSRFM(PBEparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    int ti;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);

    /* Parse old-style int arg */
    if (sscanf(tok, "%d", &ti) == 1) {
        thee->srfm = (Vsurf_Meth)ti;
        thee->setsrfm = 1;

        Vnm_print(2, "parsePBE:  Warning -- parsed deprecated \"srfm %d\" \
statement.\n", ti);
        Vnm_print(2, "parsePBE:  Please use \"srfm ");
        switch (thee->srfm) {
            case VSM_MOL:
                Vnm_print(2, "mol");
                break;
            case VSM_MOLSMOOTH:
                Vnm_print(2, "smol");
                break;
            case VSM_SPLINE:
                Vnm_print(2, "spl2");
                break;
            case VSM_SPLINE3:
                Vnm_print(2, "spl3");
                break;
            case VSM_SPLINE4:
                Vnm_print(2, "spl4");
                break;
            default:
                Vnm_print(2, "UNKNOWN");
                break;
        }
        Vnm_print(2, "\" instead.\n");
        return 1;

    /* Parse newer text-based args */
    } else if (Vstring_strcasecmp(tok, "mol") == 0) {
        thee->srfm = VSM_MOL;
        thee->setsrfm = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "smol") == 0) {
        thee->srfm = VSM_MOLSMOOTH;
        thee->setsrfm = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "spl2") == 0) {
        thee->srfm = VSM_SPLINE;
        thee->setsrfm = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "spl3") == 0) {
        thee->srfm = VSM_SPLINE3;
        thee->setsrfm = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "spl4") == 0) {
        thee->srfm = VSM_SPLINE4;
        thee->setsrfm = 1;
        return 1;
    } else {
        Vnm_print(2, "NOsh:  Unrecongnized keyword (%s) when parsing \
srfm!\n", tok);
        return -1;
    }

    return 0;

    VERROR1:
        Vnm_print(2, "parsePBE:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int PBEparm_parseSRAD(PBEparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    double tf;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing SRAD \
keyword!\n", tok);
        return -1;
    }
    thee->srad = tf;
    thee->setsrad = 1;
    return 1;

    VERROR1:
        Vnm_print(2, "parsePBE:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int PBEparm_parseSWIN(PBEparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    double tf;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing SWIN \
keyword!\n", tok);
        return -1;
    }
    thee->swin = tf;
    thee->setswin = 1;
    return 1;

    VERROR1:
        Vnm_print(2, "parsePBE:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int PBEparm_parseTEMP(PBEparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    double tf;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing TEMP \
keyword!\n", tok);
        return -1;
    }
    thee->temp = tf;
    thee->settemp = 1;
    return 1;

    VERROR1:
        Vnm_print(2, "parsePBE:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int PBEparm_parseUSEMAP(PBEparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    int ti;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    Vnm_print(0, "PBEparm_parseToken:  Read %s...\n", tok);
    if (Vstring_strcasecmp(tok, "diel") == 0) {
        thee->useDielMap = 1;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0) {
            Vnm_print(2, "NOsh:  Read non-int (%s) while parsing \
USEMAP DIEL keyword!\n", tok);
            return -1;
        }
        thee->dielMapID = ti;
        return 1;
    } else if (Vstring_strcasecmp(tok, "kappa") == 0) {
        thee->useKappaMap = 1;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0) {
            Vnm_print(2, "NOsh:  Read non-int (%s) while parsing \
USEMAP KAPPA keyword!\n", tok);
            return -1;
        }
        thee->kappaMapID = ti;
        return 1;
    } else if (Vstring_strcasecmp(tok, "pot") == 0) {
        thee->usePotMap = 1;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0) {
            Vnm_print(2, "NOsh:  Read non-int (%s) while parsing \
                      USEMAP POT keyword!\n", tok);
            return -1;
        }
        thee->potMapID = ti;
        return 1;
    } else if (Vstring_strcasecmp(tok, "charge") == 0) {
        thee->useChargeMap = 1;
        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        if (sscanf(tok, "%d", &ti) == 0) {
            Vnm_print(2, "NOsh:  Read non-int (%s) while parsing \
USEMAP CHARGE keyword!\n", tok);
            return -1;
        }
        thee->chargeMapID = ti;
        return 1;
    } else {
        Vnm_print(2, "NOsh:  Read undefined keyword (%s) while parsing \
USEMAP statement!\n", tok);
        return -1;
    }
    return 0;

    VERROR1:
        Vnm_print(2, "parsePBE:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int PBEparm_parseCALCENERGY(PBEparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    int ti;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    /* Parse number */
    if (sscanf(tok, "%d", &ti) == 1) {
        thee->calcenergy = (PBEparm_calcEnergy)ti;
        thee->setcalcenergy = 1;

        Vnm_print(2, "parsePBE:  Warning -- parsed deprecated \"calcenergy \
%d\" statement.\n", ti);
        Vnm_print(2, "parsePBE:  Please use \"calcenergy ");
        switch (thee->calcenergy) {
            case PCE_NO:
                Vnm_print(2, "no");
                break;
            case PCE_TOTAL:
                Vnm_print(2, "total");
                break;
            case PCE_COMPS:
                Vnm_print(2, "comps");
                break;
            default:
                Vnm_print(2, "UNKNOWN");
                break;
        }
        Vnm_print(2, "\" instead.\n");
        return 1;
    } else if (Vstring_strcasecmp(tok, "no") == 0) {
        thee->calcenergy = PCE_NO;
        thee->setcalcenergy = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "total") == 0) {
        thee->calcenergy = PCE_TOTAL;
        thee->setcalcenergy = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "comps") == 0) {
        thee->calcenergy = PCE_COMPS;
        thee->setcalcenergy = 1;
        return 1;
    } else {
        Vnm_print(2, "NOsh:  Unrecognized parameter (%s) while parsing \
calcenergy!\n", tok);
        return -1;
    }
    return 0;

    VERROR1:
        Vnm_print(2, "parsePBE:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int PBEparm_parseCALCFORCE(PBEparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    int ti;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    /* Parse number */
    if (sscanf(tok, "%d", &ti) == 1) {
        thee->calcforce = (PBEparm_calcForce)ti;
        thee->setcalcforce = 1;

        Vnm_print(2, "parsePBE:  Warning -- parsed deprecated \"calcforce \
%d\" statement.\n", ti);
        Vnm_print(2, "parsePBE:  Please use \"calcforce ");
        switch (thee->calcenergy) {
            case PCF_NO:
                Vnm_print(2, "no");
                break;
            case PCF_TOTAL:
                Vnm_print(2, "total");
                break;
            case PCF_COMPS:
                Vnm_print(2, "comps");
                break;
            default:
                Vnm_print(2, "UNKNOWN");
                break;
        }
        Vnm_print(2, "\" instead.\n");
        return 1;
    } else if (Vstring_strcasecmp(tok, "no") == 0) {
        thee->calcforce = PCF_NO;
        thee->setcalcforce = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "total") == 0) {
        thee->calcforce = PCF_TOTAL;
        thee->setcalcforce = 1;
        return 1;
    } else if (Vstring_strcasecmp(tok, "comps") == 0) {
        thee->calcforce = PCF_COMPS;
        thee->setcalcforce = 1;
        return 1;
    } else {
        Vnm_print(2, "NOsh:  Unrecognized parameter (%s) while parsing \
calcforce!\n", tok);
        return -1;
    }
    return 0;

    VERROR1:
        Vnm_print(2, "parsePBE:  ran out of tokens!\n");
        return -1;
}

/*----------------------------------------------------------*/
/* Added by Michael Grabe                                   */
/*----------------------------------------------------------*/

VPRIVATE int PBEparm_parseZMEM(PBEparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    double tf;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing ZMEM \
                  keyword!\n", tok);
        return -1;
    }
    thee->zmem = tf;
    thee->setzmem = 1;
    return 1;

VERROR1:
    Vnm_print(2, "parsePBE:  ran out of tokens!\n");
    return -1;
}


VPRIVATE int PBEparm_parseLMEM(PBEparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    double tf;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing LMEM \
                  keyword!\n", tok);
        return -1;
    }
    thee->Lmem = tf;
    thee->setLmem = 1;
    return 1;

VERROR1:
    Vnm_print(2, "parsePBE:  ran out of tokens!\n");
    return -1;
}

VPRIVATE int PBEparm_parseMDIE(PBEparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    double tf;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing MDIE \
                  keyword!\n", tok);
        return -1;
    }
    thee->mdie = tf;
    thee->setmdie  = 1;
    return 1;

VERROR1:
    Vnm_print(2, "parsePBE:  ran out of tokens!\n");
    return -1;
}

VPRIVATE int PBEparm_parseMEMV(PBEparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    double tf;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing MEMV \
                  keyword!\n", tok);
        return -1;
    }
    thee->memv = tf;
    thee->setmemv = 1;
    return 1;

VERROR1:
    Vnm_print(2, "parsePBE:  ran out of tokens!\n");
    return -1;
}

/*----------------------------------------------------------*/

VPRIVATE int PBEparm_parseWRITE(PBEparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE], str[VMAX_BUFSIZE]="", strnew[VMAX_BUFSIZE]="";
    Vdata_Type writetype;
    Vdata_Format writefmt;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (Vstring_strcasecmp(tok, "pot") == 0) {
        writetype = VDT_POT;
    } else if (Vstring_strcasecmp(tok, "atompot") == 0) {
        writetype = VDT_ATOMPOT;
    }  else if (Vstring_strcasecmp(tok, "charge") == 0) {
        writetype = VDT_CHARGE;
    } else if (Vstring_strcasecmp(tok, "smol") == 0) {
        writetype = VDT_SMOL;
    } else if (Vstring_strcasecmp(tok, "dielx") == 0) {
        writetype = VDT_DIELX;
    } else if (Vstring_strcasecmp(tok, "diely") == 0) {
        writetype = VDT_DIELY;
    } else if (Vstring_strcasecmp(tok, "dielz") == 0) {
        writetype = VDT_DIELZ;
    } else if (Vstring_strcasecmp(tok, "kappa") == 0) {
        writetype = VDT_KAPPA;
    } else if (Vstring_strcasecmp(tok, "sspl") == 0) {
        writetype = VDT_SSPL;
    } else if (Vstring_strcasecmp(tok, "vdw") == 0) {
        writetype = VDT_VDW;
    } else if (Vstring_strcasecmp(tok, "ivdw") == 0) {
        writetype = VDT_IVDW;
    } else if (Vstring_strcasecmp(tok, "lap") == 0) {
        writetype = VDT_LAP;
    } else if (Vstring_strcasecmp(tok, "edens") == 0) {
        writetype = VDT_EDENS;
    } else if (Vstring_strcasecmp(tok, "ndens") == 0) {
        writetype = VDT_NDENS;
    } else if (Vstring_strcasecmp(tok, "qdens") == 0) {
        writetype = VDT_QDENS;
    } else {
        Vnm_print(2, "PBEparm_parse:  Invalid data type (%s) to write!\n",
           tok);
        return -1;
    }
    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (Vstring_strcasecmp(tok, "dx") == 0) {
        writefmt = VDF_DX;
    } else if (Vstring_strcasecmp(tok, "uhbd") == 0) {
        writefmt = VDF_UHBD;
    } else if (Vstring_strcasecmp(tok, "avs") == 0) {
        writefmt = VDF_AVS;
    } else if (Vstring_strcasecmp(tok, "gz") == 0) {
        writefmt = VDF_GZ;
    } else if (Vstring_strcasecmp(tok, "flat") == 0) {
        writefmt = VDF_FLAT;
    } else {
        Vnm_print(2, "PBEparm_parse:  Invalid data format (%s) to write!\n",
           tok);
        return -1;
    }
    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (tok[0]=='"') {
        strcpy(strnew, "");
        while (tok[strlen(tok)-1] != '"') {
            strcat(str, tok);
            strcat(str, " ");
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        }
        strcat(str, tok);
        strncpy(strnew, str+1, strlen(str)-2);
        strcpy(tok, strnew);
    }
    if (thee->numwrite < (PBEPARM_MAXWRITE-1)) {
        strncpy(thee->writestem[thee->numwrite], tok, VMAX_ARGLEN);
        thee->writetype[thee->numwrite] = writetype;
        thee->writefmt[thee->numwrite] = writefmt;
        (thee->numwrite)++;
    } else {
        Vnm_print(2, "PBEparm_parse:  You have exceeded the maximum number of write statements!\n");
        Vnm_print(2, "PBEparm_parse:  Ignoring additional write statements!\n");
    }
    return 1;

    VERROR1:
        Vnm_print(2, "parsePBE:  ran out of tokens!\n");
        return -1;
}

VPRIVATE int PBEparm_parseWRITEMAT(PBEparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE], str[VMAX_BUFSIZE]="", strnew[VMAX_BUFSIZE]="";

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (Vstring_strcasecmp(tok, "poisson") == 0) {
        thee->writematflag = 0;
    } else if (Vstring_strcasecmp(tok, "full") == 0) {
        thee->writematflag = 1;
    } else {
        Vnm_print(2, "NOsh:  Invalid format (%s) while parsing \
WRITEMAT keyword!\n", tok);
        return -1;
    }

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (tok[0]=='"') {
        strcpy(strnew, "");
        while (tok[strlen(tok)-1] != '"') {
            strcat(str, tok);
            strcat(str, " ");
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
        }
        strcat(str, tok);
        strncpy(strnew, str+1, strlen(str)-2);
        strcpy(tok, strnew);
    }
    strncpy(thee->writematstem, tok, VMAX_ARGLEN);
    thee->setwritemat = 1;
    thee->writemat = 1;
    return 1;

    VERROR1:
        Vnm_print(2, "parsePBE:  ran out of tokens!\n");
        return -1;

}

VPUBLIC int PBEparm_parseToken(PBEparm *thee, char tok[VMAX_BUFSIZE],
  Vio *sock) {

    if (thee == VNULL) {
        Vnm_print(2, "parsePBE:  got NULL thee!\n");
        return -1;
    }
    if (sock == VNULL) {
        Vnm_print(2, "parsePBE:  got NULL socket!\n");
        return -1;
    }

    Vnm_print(0, "PBEparm_parseToken:  trying %s...\n", tok);

    if (Vstring_strcasecmp(tok, "mol") == 0) {
        return PBEparm_parseMOL(thee, sock);
    } else if (Vstring_strcasecmp(tok, "lpbe") == 0) {
        return PBEparm_parseLPBE(thee, sock);
    } else if (Vstring_strcasecmp(tok, "npbe") == 0) {
        return PBEparm_parseNPBE(thee, sock);
    } else if (Vstring_strcasecmp(tok, "lrpbe") == 0) {
        return PBEparm_parseLRPBE(thee, sock);
    } else if (Vstring_strcasecmp(tok, "nrpbe") == 0) {
        return PBEparm_parseNRPBE(thee, sock);
    } else if (Vstring_strcasecmp(tok, "smpbe") == 0) {
        return PBEparm_parseSMPBE(thee, sock);
    } else if (Vstring_strcasecmp(tok, "bcfl") == 0) {
        return PBEparm_parseBCFL(thee, sock);
    } else if (Vstring_strcasecmp(tok, "ion") == 0) {
        return PBEparm_parseION(thee, sock);
    } else if (Vstring_strcasecmp(tok, "pdie") == 0) {
        return PBEparm_parsePDIE(thee, sock);
    } else if (Vstring_strcasecmp(tok, "sdens") == 0) {
        return PBEparm_parseSDENS(thee, sock);
    } else if (Vstring_strcasecmp(tok, "sdie") == 0) {
        return PBEparm_parseSDIE(thee, sock);
    } else if (Vstring_strcasecmp(tok, "srfm") == 0) {
        return PBEparm_parseSRFM(thee, sock);
    } else if (Vstring_strcasecmp(tok, "srad") == 0) {
        return PBEparm_parseSRAD(thee, sock);
    } else if (Vstring_strcasecmp(tok, "swin") == 0) {
        return PBEparm_parseSWIN(thee, sock);
    } else if (Vstring_strcasecmp(tok, "temp") == 0) {
        return PBEparm_parseTEMP(thee, sock);
    } else if (Vstring_strcasecmp(tok, "usemap") == 0) {
        return PBEparm_parseUSEMAP(thee, sock);
    } else if (Vstring_strcasecmp(tok, "calcenergy") == 0) {
        return PBEparm_parseCALCENERGY(thee, sock);
    } else if (Vstring_strcasecmp(tok, "calcforce") == 0) {
        return PBEparm_parseCALCFORCE(thee, sock);
    } else if (Vstring_strcasecmp(tok, "write") == 0) {
        return PBEparm_parseWRITE(thee, sock);
    } else if (Vstring_strcasecmp(tok, "writemat") == 0) {
        return PBEparm_parseWRITEMAT(thee, sock);

    /*----------------------------------------------------------*/
    /* Added by Michael Grabe                                   */
    /*----------------------------------------------------------*/

    } else if (Vstring_strcasecmp(tok, "zmem") == 0) {
        return PBEparm_parseZMEM(thee, sock);
    } else if (Vstring_strcasecmp(tok, "Lmem") == 0) {
        return PBEparm_parseLMEM(thee, sock);
    } else if (Vstring_strcasecmp(tok, "mdie") == 0) {
        return PBEparm_parseMDIE(thee, sock);
    } else if (Vstring_strcasecmp(tok, "memv") == 0) {
        return PBEparm_parseMEMV(thee, sock);
    }

    /*----------------------------------------------------------*/

    return 0;

}
