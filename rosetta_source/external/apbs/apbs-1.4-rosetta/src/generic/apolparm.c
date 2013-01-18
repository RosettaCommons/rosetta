/**
 *  @file    apolparm.c
 *  @ingroup APOLparm
 *  @author  David Gohara
 *  @brief   Class APOLparm methods
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

#include "apolparm.h"

VEMBED(rcsid="$Id$")

#if !defined(VINLINE_MGPARM)

#endif /* if !defined(VINLINE_MGPARM) */

VPUBLIC APOLparm* APOLparm_ctor() {

    /* Set up the structure */
    APOLparm *thee = VNULL;
    thee = (APOLparm*)Vmem_malloc(VNULL, 1, sizeof(APOLparm));
    VASSERT( thee != VNULL);
    VASSERT( APOLparm_ctor2(thee) == VRC_SUCCESS );

    return thee;
}

VPUBLIC Vrc_Codes APOLparm_ctor2(APOLparm *thee) {

    int i;

    if (thee == VNULL) return VRC_FAILURE;

    thee->parsed = 0;

    thee->setgrid = 0;
    thee->setmolid = 0;
    thee->setbconc = 0;
    thee->setsdens = 0;
    thee->setdpos = 0;
    thee->setpress = 0;
    thee->setsrfm = 0;
    thee->setsrad = 0;
    thee->setswin = 0;

    thee->settemp = 0;
    thee->setgamma = 0;

    thee->setwat = 0;

    thee->sav = 0.0;
    thee->sasa = 0.0;
    thee->wcaEnergy = 0.0;

    for(i=0;i<3;i++) thee->totForce[i] = 0.0;

    return VRC_SUCCESS;
}

VPUBLIC void APOLparm_copy(
                          APOLparm *thee,
                          APOLparm *source
                          ) {

    int i;

    thee->parsed = source->parsed;

    for (i=0; i<3; i++) thee->grid[i] = source->grid[i];
    thee->setgrid = source->setgrid;

    thee->molid = source->molid;
    thee->setmolid = source->setmolid;

    thee->bconc = source->bconc ;
    thee->setbconc= source->setbconc ;

    thee->sdens = source->sdens ;
    thee->setsdens= source->setsdens ;

    thee->dpos = source->dpos ;
    thee->setdpos= source->setdpos ;

    thee->press = source->press ;
    thee->setpress = source->setpress ;

    thee->srfm = source->srfm ;
    thee->setsrfm = source->setsrfm ;

    thee->srad = source->srad ;
    thee->setsrad = source->setsrad ;

    thee->swin = source->swin ;
    thee->setswin = source->setswin ;

    thee->temp = source->temp ;
    thee->settemp = source->settemp ;

    thee->gamma = source->gamma ;
    thee->setgamma = source->setgamma ;

    thee->calcenergy = source->calcenergy ;
    thee->setcalcenergy = source->setcalcenergy ;

    thee->calcforce = source->calcforce ;
    thee->setcalcforce = source->setcalcforce ;

    thee->setwat = source->setwat ;

    thee->sav = source->sav;
    thee->sasa = source->sasa;
    thee->wcaEnergy = source->wcaEnergy;

    for(i=0;i<3;i++) thee->totForce[i] = source->totForce[i];

    return;
}

VPUBLIC void APOLparm_dtor(APOLparm **thee) {
    if ((*thee) != VNULL) {
        APOLparm_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(APOLparm), (void **)thee);
        (*thee) = VNULL;
    }

    return;
}

VPUBLIC void APOLparm_dtor2(APOLparm *thee) { ; }

VPUBLIC Vrc_Codes APOLparm_check(APOLparm *thee) {


    Vrc_Codes rc;
    rc = VRC_SUCCESS;

    if (!thee->parsed) {
        Vnm_print(2, "APOLparm_check:  not filled!\n");
        return VRC_FAILURE;
    }
    if (!thee->setgrid) {
        Vnm_print(2, "APOLparm_check:  grid not set!\n");
        rc = VRC_FAILURE;
    }
    if (!thee->setmolid) {
        Vnm_print(2, "APOLparm_check:  molid not set!\n");
        rc = VRC_FAILURE;
    }
    if (!thee->setbconc) {
        Vnm_print(2, "APOLparm_check:  bconc not set!\n");
        rc = VRC_FAILURE;
    }
    if (!thee->setsdens) {
        Vnm_print(2, "APOLparm_check:  sdens not set!\n");
        rc = VRC_FAILURE;
    }
    if (!thee->setdpos) {
        Vnm_print(2, "APOLparm_check:  dpos not set!\n");
        rc = VRC_FAILURE;
    }
    if (!thee->setpress) {
        Vnm_print(2, "APOLparm_check:  press not set!\n");
        rc = VRC_FAILURE;
    }
    if (!thee->setsrfm) {
        Vnm_print(2, "APOLparm_check:  srfm not set!\n");
        rc = VRC_FAILURE;
    }
    if (!thee->setsrad) {
        Vnm_print(2, "APOLparm_check:  srad not set!\n");
        rc = VRC_FAILURE;
    }
    if (!thee->setswin) {
        Vnm_print(2, "APOLparm_check:  swin not set!\n");
        rc = VRC_FAILURE;
    }
    if (!thee->settemp) {
        Vnm_print(2, "APOLparm_check:  temp not set!\n");
        rc = VRC_FAILURE;
    }
    if (!thee->setgamma) {
        Vnm_print(2, "APOLparm_check:  gamma not set!\n");
        rc = VRC_FAILURE;
    }
    return rc;

}

VPRIVATE Vrc_Codes APOLparm_parseGRID(APOLparm *thee, Vio *sock) {

    char tok[VMAX_BUFSIZE];
    double tf;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GRID \
keyword!\n", tok);
        return VRC_WARNING;
    } else thee->grid[0] = tf;
    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GRID \
keyword!\n", tok);
        return VRC_WARNING;
    } else thee->grid[1] = tf;
    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GRID \
keyword!\n", tok);
        return VRC_WARNING;
    } else thee->grid[2] = tf;
    thee->setgrid = 1;
    return VRC_SUCCESS;

VERROR1:
        Vnm_print(2, "parseAPOL:  ran out of tokens!\n");
    return VRC_WARNING;
}

VPRIVATE Vrc_Codes APOLparm_parseMOL(APOLparm *thee, Vio *sock) {
    int ti;
    char tok[VMAX_BUFSIZE];

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%d", &ti) == 0) {
        Vnm_print(2, "NOsh:  Read non-int (%s) while parsing MOL \
keyword!\n", tok);
        return VRC_WARNING;
    }
    thee->molid = ti;
    thee->setmolid = 1;
    return VRC_SUCCESS;

VERROR1:
        Vnm_print(2, "parseAPOL:  ran out of tokens!\n");
    return VRC_WARNING;
}

VPRIVATE Vrc_Codes APOLparm_parseSRFM(APOLparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);

    if (Vstring_strcasecmp(tok, "sacc") == 0) {
        thee->srfm = VSM_MOL;
        thee->setsrfm = 1;
        return VRC_SUCCESS;
    } else {
        Vnm_print(2, "parseAPOL: Unrecongnized keyword (%s) when parsing srfm!\n", tok);
        Vnm_print(2, "parseAPOL: Accepted values for srfm = sacc\n");
        return VRC_WARNING;
    }

    return VRC_FAILURE;

VERROR1:
        Vnm_print(2, "parseAPOL:  ran out of tokens!\n");
    return VRC_WARNING;
}

VPRIVATE Vrc_Codes APOLparm_parseSRAD(APOLparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    double tf;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing SRAD \
keyword!\n", tok);
        return VRC_WARNING;
    }
    thee->srad = tf;
    thee->setsrad = 1;
    return VRC_SUCCESS;

VERROR1:
        Vnm_print(2, "parseAPOL:  ran out of tokens!\n");
    return VRC_WARNING;
}

VPRIVATE Vrc_Codes APOLparm_parseSWIN(APOLparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    double tf;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing SWIN \
keyword!\n", tok);
        return VRC_WARNING;
    }
    thee->swin = tf;
    thee->setswin = 1;
    return VRC_SUCCESS;

VERROR1:
        Vnm_print(2, "parseAPOL:  ran out of tokens!\n");
    return VRC_WARNING;
}

VPRIVATE Vrc_Codes APOLparm_parseTEMP(APOLparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    double tf;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing TEMP \
keyword!\n", tok);
        return VRC_WARNING;
    }
    thee->temp = tf;
    thee->settemp = 1;
    return VRC_SUCCESS;

VERROR1:
        Vnm_print(2, "parseAPOL:  ran out of tokens!\n");
    return VRC_WARNING;
}

VPRIVATE Vrc_Codes APOLparm_parseGAMMA(APOLparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    double tf;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing GAMMA \
keyword!\n", tok);
        return VRC_WARNING;
    }
    thee->gamma = tf;
    thee->setgamma = 1;
    return VRC_SUCCESS;

VERROR1:
        Vnm_print(2, "parseAPOL:  ran out of tokens!\n");
    return VRC_WARNING;
}

VPRIVATE Vrc_Codes APOLparm_parseCALCENERGY(APOLparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    int ti;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    /* Parse number */
    if (sscanf(tok, "%d", &ti) == 1) {
        thee->calcenergy = (APOLparm_calcEnergy)ti;
        thee->setcalcenergy = 1;

        Vnm_print(2, "parseAPOL:  Warning -- parsed deprecated \"calcenergy \
%d\" statement.\n", ti);
        Vnm_print(2, "parseAPOL:  Please use \"calcenergy ");
        switch (thee->calcenergy) {
            case ACE_NO:
                Vnm_print(2, "no");
                break;
            case ACE_TOTAL:
                Vnm_print(2, "total");
                break;
            case ACE_COMPS:
                Vnm_print(2, "comps");
                break;
            default:
                Vnm_print(2, "UNKNOWN");
                break;
        }
        Vnm_print(2, "\" instead.\n");
        return VRC_SUCCESS;
    } else if (Vstring_strcasecmp(tok, "no") == 0) {
        thee->calcenergy = ACE_NO;
        thee->setcalcenergy = 1;
        return VRC_SUCCESS;
    } else if (Vstring_strcasecmp(tok, "total") == 0) {
        thee->calcenergy = ACE_TOTAL;
        thee->setcalcenergy = 1;
        return VRC_SUCCESS;
    } else if (Vstring_strcasecmp(tok, "comps") == 0) {
        thee->calcenergy = ACE_COMPS;
        thee->setcalcenergy = 1;
        return VRC_SUCCESS;
    } else {
        Vnm_print(2, "NOsh:  Unrecognized parameter (%s) while parsing \
calcenergy!\n", tok);
        return VRC_WARNING;
    }
    return VRC_FAILURE;

VERROR1:
        Vnm_print(2, "parseAPOL:  ran out of tokens!\n");
    return VRC_WARNING;
}

VPRIVATE Vrc_Codes APOLparm_parseCALCFORCE(APOLparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    int ti;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    /* Parse number */
    if (sscanf(tok, "%d", &ti) == 1) {
        thee->calcforce = (APOLparm_calcForce)ti;
        thee->setcalcforce = 1;

        Vnm_print(2, "parseAPOL:  Warning -- parsed deprecated \"calcforce \
%d\" statement.\n", ti);
        Vnm_print(2, "parseAPOL:  Please use \"calcforce ");
        switch (thee->calcenergy) {
            case ACF_NO:
                Vnm_print(2, "no");
                break;
            case ACF_TOTAL:
                Vnm_print(2, "total");
                break;
            case ACF_COMPS:
                Vnm_print(2, "comps");
                break;
            default:
                Vnm_print(2, "UNKNOWN");
                break;
        }
        Vnm_print(2, "\" instead.\n");
        return VRC_SUCCESS;
    } else if (Vstring_strcasecmp(tok, "no") == 0) {
        thee->calcforce = ACF_NO;
        thee->setcalcforce = 1;
        return VRC_SUCCESS;
    } else if (Vstring_strcasecmp(tok, "total") == 0) {
        thee->calcforce = ACF_TOTAL;
        thee->setcalcforce = 1;
        return VRC_SUCCESS;
    } else if (Vstring_strcasecmp(tok, "comps") == 0) {
        thee->calcforce = ACF_COMPS;
        thee->setcalcforce = 1;
        return VRC_SUCCESS;
    } else {
        Vnm_print(2, "NOsh:  Unrecognized parameter (%s) while parsing \
calcforce!\n", tok);
        return VRC_WARNING;
    }
    return VRC_FAILURE;

VERROR1:
        Vnm_print(2, "parseAPOL:  ran out of tokens!\n");
    return VRC_WARNING;
}

VPRIVATE Vrc_Codes APOLparm_parseBCONC(APOLparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    double tf;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing BCONC \
keyword!\n", tok);
        return VRC_WARNING;
    }
    thee->bconc = tf;
    thee->setbconc = 1;
    return VRC_SUCCESS;

VERROR1:
        Vnm_print(2, "parseAPOL:  ran out of tokens!\n");
    return VRC_WARNING;
}

VPRIVATE Vrc_Codes APOLparm_parseSDENS(APOLparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    double tf;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing SDENS \
keyword!\n", tok);
        return VRC_WARNING;
    }
    thee->sdens = tf;
    thee->setsdens = 1;
    return VRC_SUCCESS;

VERROR1:
        Vnm_print(2, "parseAPOL:  ran out of tokens!\n");
    return VRC_WARNING;
}

VPRIVATE Vrc_Codes APOLparm_parseDPOS(APOLparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    double tf;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing SDENS \
keyword!\n", tok);
        return VRC_WARNING;
    }
    thee->dpos = tf;
    thee->setdpos = 1;

    if(thee->dpos < 0.001){
        Vnm_print(1,"\nWARNING WARNING WARNING WARNING WARNING\n");
        Vnm_print(1,"NOsh: dpos is set to a very small value.\n");
        Vnm_print(1,"NOsh: If you are not using a PQR file, you can \
safely ignore this message.\n");
        Vnm_print(1,"NOsh: Otherwise please choose a value greater than \
or equal to 0.001.\n\n");
    }

    return VRC_SUCCESS;

VERROR1:
        Vnm_print(2, "parseAPOL:  ran out of tokens!\n");
    return VRC_WARNING;
}

VPRIVATE Vrc_Codes APOLparm_parsePRESS(APOLparm *thee, Vio *sock) {
    char tok[VMAX_BUFSIZE];
    double tf;

    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &tf) == 0) {
        Vnm_print(2, "NOsh:  Read non-float (%s) while parsing PRESS \
keyword!\n", tok);
        return VRC_WARNING;
    }
    thee->press = tf;
    thee->setpress = 1;
    return VRC_SUCCESS;

VERROR1:
        Vnm_print(2, "parseAPOL:  ran out of tokens!\n");
    return VRC_WARNING;
}

VPUBLIC Vrc_Codes APOLparm_parseToken(APOLparm *thee, char tok[VMAX_BUFSIZE],
  Vio *sock) {

    if (thee == VNULL) {
        Vnm_print(2, "parseAPOL:  got NULL thee!\n");
        return VRC_WARNING;
    }

    if (sock == VNULL) {
        Vnm_print(2, "parseAPOL:  got NULL socket!\n");
        return VRC_WARNING;
    }

    if (Vstring_strcasecmp(tok, "mol") == 0) {
        return APOLparm_parseMOL(thee, sock);
    } else if (Vstring_strcasecmp(tok, "grid") == 0) {
        return APOLparm_parseGRID(thee, sock);
    } else if (Vstring_strcasecmp(tok, "dime") == 0) {
        Vnm_print(2, "APOLparm_parseToken:  The DIME and GLEN keywords for APOLAR have been replaced with GRID.\n");
        Vnm_print(2, "APOLparm_parseToken:  Please see the APBS User Guide for more information.\n");
        return VRC_WARNING;
    } else if (Vstring_strcasecmp(tok, "glen") == 0) {
        Vnm_print(2, "APOLparm_parseToken:  The DIME and GLEN keywords for APOLAR have been replaced with GRID.\n");
        Vnm_print(2, "APOLparm_parseToken:  Please see the APBS User Guide for more information.\n");
        return VRC_WARNING;
    } else if (Vstring_strcasecmp(tok, "bconc") == 0) {
        return APOLparm_parseBCONC(thee, sock);
    } else if (Vstring_strcasecmp(tok, "sdens") == 0) {
        return APOLparm_parseSDENS(thee, sock);
    } else if (Vstring_strcasecmp(tok, "dpos") == 0) {
        return APOLparm_parseDPOS(thee, sock);
    }  else if (Vstring_strcasecmp(tok, "srfm") == 0) {
        return APOLparm_parseSRFM(thee, sock);
    } else if (Vstring_strcasecmp(tok, "srad") == 0) {
        return APOLparm_parseSRAD(thee, sock);
    } else if (Vstring_strcasecmp(tok, "swin") == 0) {
        return APOLparm_parseSWIN(thee, sock);
    } else if (Vstring_strcasecmp(tok, "temp") == 0) {
        return APOLparm_parseTEMP(thee, sock);
    } else if (Vstring_strcasecmp(tok, "gamma") == 0) {
        return APOLparm_parseGAMMA(thee, sock);
    } else if (Vstring_strcasecmp(tok, "press") == 0) {
        return APOLparm_parsePRESS(thee, sock);
    } else if (Vstring_strcasecmp(tok, "calcenergy") == 0) {
        return APOLparm_parseCALCENERGY(thee, sock);
    } else if (Vstring_strcasecmp(tok, "calcforce") == 0) {
        return APOLparm_parseCALCFORCE(thee, sock);
    } else {
        Vnm_print(2, "parseAPOL:  Unrecognized keyword (%s)!\n", tok);
        return VRC_WARNING;
    }


    return VRC_FAILURE;

}
