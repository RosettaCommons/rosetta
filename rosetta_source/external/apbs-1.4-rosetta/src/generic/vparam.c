/**
 *  @file    vparam.c
 *  @ingroup Vparam
 *  @author  Nathan Baker
 *  @brief   Class Vparam methods
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

#include "vparam.h"

VEMBED(rcsid="$Id$")

/**
 * @brief  Whitespace characters for socket reads
 * @ingroup  Vparam
 */
VPRIVATE char *MCwhiteChars = " =,;\t\n\r";

/**
 * @brief  Comment characters for socket reads
 * @ingroup  Vparam
 */
VPRIVATE char *MCcommChars  = "#%";

/**
 * @brief  Whitespace characters for XML socket reads
 * @ingroup Vparam
 */
VPRIVATE char *MCxmlwhiteChars = " =,;\t\n\r<>";

/**
 * @brief  Read a single line of the flat file database
 * @author  Nathan Baker
 * @ingroup  Vparam
 * @param  sock  Socket ready for reading
 * @param  atom  Atom to hold parsed data
 * @returns 1 if successful, 0 otherwise
 */
VPRIVATE int readFlatFileLine(Vio *sock, Vparam_AtomData *atom);

/**
 * @brief  Read atom information from an XML file
 * @author  Todd Dolinsky
 * @ingroup  Vparam
 * @param  sock  Socket ready for reading
 * @param  atom  Atom to hold parsed data
 * @returns 1 if successful, 0 otherwise
 */
VPRIVATE int readXMLFileAtom(Vio *sock, Vparam_AtomData *atom);


#if !defined(VINLINE_VPARAM)

VPUBLIC unsigned long int Vparam_memChk(Vparam *thee) {
    if (thee == VNULL) return 0;
    return Vmem_bytes(thee->vmem);
}

#endif /* if !defined(VINLINE_VPARAM) */

VPUBLIC Vparam_AtomData* Vparam_AtomData_ctor() {

    Vparam_AtomData *thee = VNULL;

    /* Set up the structure */
    thee = (Vparam_AtomData*)Vmem_malloc(VNULL, 1, sizeof(Vparam_AtomData) );
    VASSERT(thee != VNULL);
    VASSERT(Vparam_AtomData_ctor2(thee));

    return thee;
}

VPUBLIC int Vparam_AtomData_ctor2(Vparam_AtomData *thee) { return 1; }

VPUBLIC void Vparam_AtomData_dtor(Vparam_AtomData **thee) {

    if ((*thee) != VNULL) {
        Vparam_AtomData_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vparam_AtomData), (void **)thee);
        (*thee) = VNULL;
    }

}

VPUBLIC void Vparam_AtomData_dtor2(Vparam_AtomData *thee) { ; }

VPUBLIC Vparam_ResData* Vparam_ResData_ctor(Vmem *mem) {

    Vparam_ResData *thee = VNULL;

    /* Set up the structure */
    thee = (Vparam_ResData*)Vmem_malloc(mem, 1, sizeof(Vparam_ResData) );
    VASSERT(thee != VNULL);
    VASSERT(Vparam_ResData_ctor2(thee, mem));

    return thee;
}

VPUBLIC int Vparam_ResData_ctor2(Vparam_ResData *thee, Vmem *mem) {

    if (thee == VNULL) {
        Vnm_print(2, "Vparam_ResData_ctor2:  Got VNULL thee!\n");
        return 0;
    }
    thee->vmem = mem;
    thee->nAtomData = 0;
    thee->atomData = VNULL;

    return 1;
}

VPUBLIC void Vparam_ResData_dtor(Vparam_ResData **thee) {

    if ((*thee) != VNULL) {
        Vparam_ResData_dtor2(*thee);
        Vmem_free((*thee)->vmem, 1, sizeof(Vparam_ResData), (void **)thee);
        (*thee) = VNULL;
    }

}

VPUBLIC void Vparam_ResData_dtor2(Vparam_ResData *thee) {

    if (thee == VNULL) return;
    if (thee->nAtomData > 0) {
        Vmem_free(thee->vmem, thee->nAtomData, sizeof(Vparam_AtomData),
          (void **)&(thee->atomData));
    }
    thee->nAtomData = 0;
    thee->atomData = VNULL;
}

VPUBLIC Vparam* Vparam_ctor() {

    Vparam *thee = VNULL;

    /* Set up the structure */
    thee = (Vparam*)Vmem_malloc(VNULL, 1, sizeof(Vparam) );
    VASSERT(thee != VNULL);
    VASSERT(Vparam_ctor2(thee));

    return thee;
}

VPUBLIC int Vparam_ctor2(Vparam *thee) {

    if (thee == VNULL) {
        Vnm_print(2, "Vparam_ctor2: got VNULL thee!\n");
        return 0;
    }

    thee->vmem = VNULL;
    thee->vmem = Vmem_ctor("APBS:VPARAM");
    if (thee->vmem == VNULL) {
        Vnm_print(2, "Vparam_ctor2:  failed to init Vmem!\n");
        return 0;
    }

    thee->nResData = 0;
    thee->resData = VNULL;

    return 1;
}

VPUBLIC void Vparam_dtor(Vparam **thee) {

    if ((*thee) != VNULL) {
        Vparam_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Vparam), (void **)thee);
        (*thee) = VNULL;
    }

}

VPUBLIC void Vparam_dtor2(Vparam *thee) {

    int i;

    if (thee == VNULL) return;

    /* Destroy the residue data */
    for (i=0; i<thee->nResData; i++) Vparam_ResData_dtor2(&(thee->resData[i]));
    if (thee->nResData > 0) Vmem_free(thee->vmem, thee->nResData,
      sizeof(Vparam_ResData), (void **)&(thee->resData));
    thee->nResData = 0;
    thee->resData = VNULL;

    if (thee->vmem != VNULL) Vmem_dtor(&(thee->vmem));
    thee->vmem = VNULL;

}

VPUBLIC Vparam_ResData* Vparam_getResData(Vparam *thee,
  char resName[VMAX_ARGLEN]) {

    int i;
    Vparam_ResData *res = VNULL;

    VASSERT(thee != VNULL);

    if ((thee->nResData == 0) || (thee->resData == VNULL)) {
        res = VNULL;
        return res;
    }

    /* Look for the matching residue */
    for (i=0; i<thee->nResData; i++) {
        res = &(thee->resData[i]);
        if (Vstring_strcasecmp(resName, res->name) == 0) return res;

    }

    /* Didn't find a matching residue */
    res = VNULL;
    Vnm_print(2, "Vparam_getResData:  unable to find res=%s\n", resName);
    return res;
}

VPUBLIC Vparam_AtomData* Vparam_getAtomData(Vparam *thee,
  char resName[VMAX_ARGLEN], char atomName[VMAX_ARGLEN]) {

    int i;
    Vparam_ResData *res = VNULL;
    Vparam_AtomData *atom = VNULL;

    VASSERT(thee != VNULL);

    if ((thee->nResData == 0) || (thee->resData == VNULL)) {
        atom = VNULL;
        return atom;
    }

    /* Look for the matching residue */
    res = Vparam_getResData(thee, resName);
    if (res == VNULL) {
        atom = VNULL;
        Vnm_print(2, "Vparam_getAtomData:  Unable to find residue %s!\n", resName);
        return atom;
    }
    for (i=0; i<res->nAtomData; i++) {
        atom = &(res->atomData[i]);
        if (atom == VNULL) {
            Vnm_print(2, "Vparam_getAtomData:  got NULL atom!\n");
            return VNULL;
        }
        if (Vstring_strcasecmp(atomName, atom->atomName) == 0) {
            return atom;
        }
    }

    /* Didn't find a matching atom/residue */
    atom = VNULL;
    Vnm_print(2, "Vparam_getAtomData:  unable to find atom '%s', res '%s'\n",
      atomName, resName);
    return atom;
}

VPUBLIC int Vparam_readXMLFile(Vparam *thee, const char *iodev,
  const char *iofmt, const char *thost, const char *fname) {

    int i, ires, natoms, nalloc, ralloc;
    Vparam_AtomData *atoms = VNULL;
    Vparam_AtomData *tatoms = VNULL;
    Vparam_AtomData *atom = VNULL;
    Vparam_ResData *res = VNULL;
    Vparam_ResData *residues = VNULL;
    Vparam_ResData *tresidues = VNULL;
    Vio *sock = VNULL;
    char currResName[VMAX_ARGLEN];
    char tok[VMAX_ARGLEN];
    char endtag[VMAX_ARGLEN];

    VASSERT(thee != VNULL);

    /* Setup communication */
    sock = Vio_ctor(iodev,iofmt,thost,fname,"r");
    if (sock == VNULL) {
        Vnm_print(2, "Vparam_readXMLFile: Problem opening virtual socket %s\n",
          fname);
        return 0;
    }
    if (Vio_accept(sock, 0) < 0) {
        Vnm_print(2, "Vparam_readXMLFile: Problem accepting virtual socket %s\n",
          fname);
        return 0;
    }
    Vio_setWhiteChars(sock, MCxmlwhiteChars);
    Vio_setCommChars(sock, MCcommChars);

    /* Clear existing parameters */
    if (thee->nResData > 0) {
        Vnm_print(2, "WARNING -- CLEARING PARAMETER DATABASE!\n");
        for (i=0; i<thee->nResData; i++) {
            Vparam_ResData_dtor2(&(thee->resData[i]));
        }
        Vmem_free(thee->vmem, thee->nResData,
          sizeof(Vparam_ResData), (void **)&(thee->resData));
    }

    strcpy(endtag,"/");

    /* Set up temporary residue list */

    ralloc = 50;
    residues = (Vparam_ResData*)Vmem_malloc(thee->vmem, ralloc, sizeof(Vparam_ResData));

    /* Read until we run out of entries, allocating space as needed */
    while (1) {

        VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);

        /* The first token should be the start tag */

        if (Vstring_strcasecmp(endtag, "/") == 0) strcat(endtag, tok);

        if (Vstring_strcasecmp(tok, "residue") == 0) {
          if (thee->nResData >= ralloc) {
                tresidues = (Vparam_ResData*)Vmem_malloc(thee->vmem, 2*ralloc, sizeof(Vparam_ResData));
                VASSERT(tresidues != VNULL);
                for (i=0; i<thee->nResData; i++) {
                    Vparam_ResData_copyTo(&(residues[i]), &(tresidues[i]));
                }
                Vmem_free(thee->vmem, ralloc, sizeof(Vparam_ResData),
                          (void **)&(residues));
                residues = tresidues;
                tresidues = VNULL;
                ralloc = 2*ralloc;
            }

          /* Initial space for this residue's atoms */
          nalloc = 20;
          natoms = 0;
          atoms = (Vparam_AtomData*)Vmem_malloc(thee->vmem, nalloc, sizeof(Vparam_AtomData));

        } else if (Vstring_strcasecmp(tok, "name") == 0) {
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);  /* value */
            strcpy(currResName, tok);
            VJMPERR1(Vio_scanf(sock, "%s", tok) == 1); /* </name> */
        } else if (Vstring_strcasecmp(tok, "atom") == 0) {
            if (natoms >= nalloc) {
                tatoms = (Vparam_AtomData*)Vmem_malloc(thee->vmem, 2*nalloc, sizeof(Vparam_AtomData));
                VASSERT(tatoms != VNULL);
                for (i=0; i<natoms; i++) {
                    Vparam_AtomData_copyTo(&(atoms[i]), &(tatoms[i]));
                }
                Vmem_free(thee->vmem, nalloc, sizeof(Vparam_AtomData),
                          (void **)&(atoms));
                atoms = tatoms;
                tatoms = VNULL;
                nalloc = 2*nalloc;
            }
            atom = &(atoms[natoms]);
            if (!readXMLFileAtom(sock, atom)) break;
            natoms++;

        } else if (Vstring_strcasecmp(tok, "/residue") == 0) {

          res = &(residues[thee->nResData]);
          Vparam_ResData_ctor2(res, thee->vmem);
          res->atomData = (Vparam_AtomData*)Vmem_malloc(thee->vmem, natoms,
                                      sizeof(Vparam_AtomData));
          res->nAtomData = natoms;
          strcpy(res->name, currResName);
          for (i=0; i<natoms; i++) {
              strcpy(atoms[i].resName, currResName);
              Vparam_AtomData_copyTo(&(atoms[i]), &(res->atomData[i]));
          }
          Vmem_free(thee->vmem, nalloc, sizeof(Vparam_AtomData), (void **)&(atoms));
          (thee->nResData)++;

        } else if (Vstring_strcasecmp(tok, endtag) == 0) break;
    }

    /* Initialize and copy the residues into the Vparam object */

    thee->resData = (Vparam_ResData*)Vmem_malloc(thee->vmem, thee->nResData,
                                sizeof(Vparam_ResData));
    for (ires=0; ires<thee->nResData; ires++) {
        Vparam_ResData_copyTo(&(residues[ires]), &(thee->resData[ires]));
    }

    /* Destroy temporary atom space */
    Vmem_free(thee->vmem, ralloc, sizeof(Vparam_ResData), (void **)&(residues));

    /* Shut down communication */
    Vio_acceptFree(sock);
    Vio_dtor(&sock);

    return 1;

VERROR1:
    Vnm_print(2, "Vparam_readXMLFile: Got unexpected EOF reading parameter file!\n");
    return 0;

}

VPUBLIC int Vparam_readFlatFile(Vparam *thee, const char *iodev,
  const char *iofmt, const char *thost, const char *fname) {

    int i, iatom, jatom, ires, natoms, nalloc;
    Vparam_AtomData *atoms = VNULL;
    Vparam_AtomData *tatoms = VNULL;
    Vparam_AtomData *atom = VNULL;
    Vparam_ResData *res = VNULL;
    Vio *sock = VNULL;
    char currResName[VMAX_ARGLEN];

    VASSERT(thee != VNULL);

    /* Setup communication */
    sock = Vio_ctor(iodev,iofmt,thost,fname,"r");
    if (sock == VNULL) {
        Vnm_print(2, "Vparam_readFlatFile: Problem opening virtual socket %s\n",
          fname);
        return 0;
    }
    if (Vio_accept(sock, 0) < 0) {
        Vnm_print(2, "Vparam_readFlatFile: Problem accepting virtual socket %s\n",
          fname);
        return 0;
    }
    Vio_setWhiteChars(sock, MCwhiteChars);
    Vio_setCommChars(sock, MCcommChars);

    /* Clear existing parameters */
    if (thee->nResData > 0) {
        Vnm_print(2, "WARNING -- CLEARING PARAMETER DATABASE!\n");
        for (i=0; i<thee->nResData; i++) {
            Vparam_ResData_dtor2(&(thee->resData[i]));
        }
        Vmem_free(thee->vmem, thee->nResData,
          sizeof(Vparam_ResData), (void **)&(thee->resData));
    }

    /* Initial space for atoms */
    nalloc = 200;
    natoms = 0;
    atoms = (Vparam_AtomData*)Vmem_malloc(thee->vmem, nalloc, sizeof(Vparam_AtomData));

    /* Read until we run out of entries, allocating space as needed */
    while (1) {
        if (natoms >= nalloc) {
            tatoms = (Vparam_AtomData*)Vmem_malloc(thee->vmem, 2*nalloc, sizeof(Vparam_AtomData));
            VASSERT(tatoms != VNULL);
            for (i=0; i<natoms; i++) {
                Vparam_AtomData_copyTo(&(atoms[i]), &(tatoms[i]));
            }
            Vmem_free(thee->vmem, nalloc, sizeof(Vparam_AtomData),
              (void **)&(atoms));
            atoms = tatoms;
            tatoms = VNULL;
            nalloc = 2*nalloc;
        }
        atom = &(atoms[natoms]);
        if (!readFlatFileLine(sock, atom)) break;
        natoms++;
    }
    if (natoms == 0) return 0;

    /* Count the number of residues */
    thee->nResData = 1;
    strcpy(currResName, atoms[0].resName);
    for (i=1; i<natoms; i++) {
        if (Vstring_strcasecmp(atoms[i].resName, currResName) != 0) {
            strcpy(currResName, atoms[i].resName);
            (thee->nResData)++;
        }
    }

    /* Create the residues */
    thee->resData = (Vparam_ResData*)Vmem_malloc(thee->vmem, thee->nResData,
      sizeof(Vparam_ResData));
    VASSERT(thee->resData != VNULL);
    for (i=0; i<(thee->nResData); i++) {
        res = &(thee->resData[i]);
        Vparam_ResData_ctor2(res, thee->vmem);
    }

    /* Count the number of atoms per residue */
    ires = 0;
    res = &(thee->resData[ires]);
    res->nAtomData = 1;
    strcpy(res->name, atoms[0].resName);
    for (i=1; i<natoms; i++) {
        if (Vstring_strcasecmp(atoms[i].resName, res->name) != 0) {
            (ires)++;
            res = &(thee->resData[ires]);
            res->nAtomData = 1;
            strcpy(res->name, atoms[i].resName);
        } else (res->nAtomData)++;
    }

    /* Allocate per-residue space for atoms */
    for (ires=0; ires<thee->nResData; ires++) {
        res = &(thee->resData[ires]);
        res->atomData = (Vparam_AtomData*)Vmem_malloc(thee->vmem, res->nAtomData,
          sizeof(Vparam_AtomData));
    }

    /* Copy atoms into residues */
    iatom = 0;
    Vparam_AtomData_copyTo(&(atoms[0]), &(res->atomData[iatom]));
    for (ires=0; ires<thee->nResData; ires++) {
        res = &(thee->resData[ires]);
        for (jatom=0; jatom<res->nAtomData; jatom++) {
            Vparam_AtomData_copyTo(&(atoms[iatom]), &(res->atomData[jatom]));
            iatom++;
        }
    }


    /* Shut down communication */
    Vio_acceptFree(sock);
    Vio_dtor(&sock);

    /* Destroy temporary atom space */
    Vmem_free(thee->vmem, nalloc, sizeof(Vparam_AtomData), (void **)&(atoms));

    return 1;

}

VEXTERNC void Vparam_AtomData_copyTo(Vparam_AtomData *thee,
  Vparam_AtomData *dest) {

    VASSERT(thee != VNULL);
    VASSERT(dest != VNULL);

    strcpy(dest->atomName, thee->atomName);
    strcpy(dest->resName, thee->resName);
    dest->charge = thee->charge;
    dest->radius = thee->radius;
    dest->epsilon = thee->epsilon;

}

VEXTERNC void Vparam_ResData_copyTo(Vparam_ResData *thee,
  Vparam_ResData *dest) {

    int i;

    VASSERT(thee != VNULL);
    VASSERT(dest != VNULL);

    strcpy(dest->name, thee->name);
    dest->vmem = thee->vmem;
    dest->nAtomData = thee->nAtomData;

    dest->atomData = (Vparam_AtomData*)Vmem_malloc(thee->vmem, dest->nAtomData,
                                 sizeof(Vparam_AtomData));

    for (i=0; i<dest->nAtomData; i++) {
      Vparam_AtomData_copyTo(&(thee->atomData[i]), &(dest->atomData[i]));
    }
    Vmem_free(thee->vmem, thee->nAtomData, sizeof(Vparam_AtomData),
              (void **)&(thee->atomData));
}

VEXTERNC void Vparam_AtomData_copyFrom(Vparam_AtomData *thee,
  Vparam_AtomData *src) {  Vparam_AtomData_copyTo(src, thee); }

VPRIVATE int readXMLFileAtom(Vio *sock, Vparam_AtomData *atom) {

    double dtmp;
    char tok[VMAX_BUFSIZE];
    int chgflag, radflag, nameflag;

    VASSERT(atom != VNULL);

    if (Vio_scanf(sock, "%s", tok) != 1) return 0;

    chgflag = 0;
    radflag = 0;
    nameflag = 0;

    while (1)
      {
          if (Vstring_strcasecmp(tok, "name") == 0) {
              VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
              if (strlen(tok) > VMAX_ARGLEN) {
                  Vnm_print(2, "Vparam_readXMLFileAtom:  string (%s) too long \
(%d)!\n", tok, strlen(tok));
                  return 0;
              }
              nameflag = 1;
              strcpy(atom->atomName, tok);
          } else if (Vstring_strcasecmp(tok, "charge") == 0) {
              VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
              if (sscanf(tok, "%lf", &dtmp) != 1) {
                  Vnm_print(2, "Vparam_readXMLFileAtom:  Unexpected token (%s) while \
parsing charge!\n", tok);
                  return 0;
              }
              chgflag = 1;
              atom->charge = dtmp;
          }  else if (Vstring_strcasecmp(tok, "radius") == 0) {
              VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
              if (sscanf(tok, "%lf", &dtmp) != 1) {
                  Vnm_print(2, "Vparam_readXMLFileAtom:  Unexpected token (%s) while \
parsing radius!\n", tok);
                  return 0;
              }
              radflag = 1;
              atom->radius = dtmp;
          }  else if (Vstring_strcasecmp(tok, "epsilon") == 0) {
              VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
              if (sscanf(tok, "%lf", &dtmp) != 1) {
                  Vnm_print(2, "Vparam_readXMLFileAtom:  Unexpected token (%s) while \
parsing epsilon!\n", tok);
                  return 0;
              }
              atom->epsilon = dtmp;
          } else if ((Vstring_strcasecmp(tok, "/atom") == 0) ||
                     (Vstring_strcasecmp(tok, "atom") == 0)){
                if (chgflag && radflag && nameflag) return 1;
                else if (!chgflag) {
                  Vnm_print(2, "Vparam_readXMLFileAtom: Reached end of atom without \
setting the charge!\n");
                  return 0;
                } else if (!radflag) {
                  Vnm_print(2, "Vparam_readXMLFileAtom: Reached end of atom without \
setting the radius!\n");
                  return 0;
                } else if (!nameflag) {
                  Vnm_print(2, "Vparam_readXMLFileAtom: Reached end of atom without \
setting the name!\n");
                  return 0;
                }
          }
          VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
      }

    /* If we get here something wrong has happened */

    VJMPERR1(1);

VERROR1:
    Vnm_print(2, "Vparam_readXMLFileAtom: Got unexpected EOF reading parameter file!\n");
    return 0;

}

VPRIVATE int readFlatFileLine(Vio *sock, Vparam_AtomData *atom) {

    double dtmp;
    char tok[VMAX_BUFSIZE];

    VASSERT(atom != VNULL);

    if (Vio_scanf(sock, "%s", tok) != 1) return 0;
    if (strlen(tok) > VMAX_ARGLEN) {
        Vnm_print(2, "Vparam_readFlatFile:  string (%s) too long (%d)!\n",
          tok, strlen(tok));
        return 0;
    }
    strcpy(atom->resName, tok);
    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (strlen(tok) > VMAX_ARGLEN) {
        Vnm_print(2, "Vparam_readFlatFile:  string (%s) too long (%d)!\n",
          tok, strlen(tok));
        return 0;
    }
    strcpy(atom->atomName, tok);
    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &dtmp) != 1) {
        Vnm_print(2, "Vparam_readFlatFile:  Unexpected token (%s) while \
parsing charge!\n", tok);
        return 0;
    }
    atom->charge = dtmp;
    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &dtmp) != 1) {
        Vnm_print(2, "Vparam_readFlatFile:  Unexpected token (%s) while \
parsing radius!\n", tok);
        return 0;
    }
    atom->radius = dtmp;
    VJMPERR1(Vio_scanf(sock, "%s", tok) == 1);
    if (sscanf(tok, "%lf", &dtmp) != 1) {
        Vnm_print(2, "Vparam_readFlatFile:  Unexpected token (%s) while \
parsing radius!\n", tok);
        return 0;
    }
    atom->epsilon = dtmp;

    return 1;

VERROR1:
    Vnm_print(2, "Vparam_readFlatFile: Got unexpected EOF reading parameter file!\n");
    return 0;
}
