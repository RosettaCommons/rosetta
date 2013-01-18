/**
 *  @file    valist.c
 *  @author  Nathan Baker
 *  @brief   Class Valist methods
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

#include "valist.h"

VEMBED(rcsid="$Id$")

VPRIVATE char *Valist_whiteChars = " \t\r\n";
VPRIVATE char *Valist_commChars  = "#%";
VPRIVATE char *Valist_xmlwhiteChars = " \t\r\n<>";

#if !defined(VINLINE_VATOM)

VPUBLIC double Valist_getCenterX(Valist *thee) {

    if (thee == NULL) {
        Vnm_print(2, "Valist_getCenterX:  Found null pointer when getting the center of X coordinate!\n");
        VASSERT(0);
    }
    return thee->center[0];

}

VPUBLIC double Valist_getCenterY(Valist *thee) {

    if (thee == NULL) {
        Vnm_print(2, "Valist_getCenterY:  Found null pointer when getting the center of Y coordinate!\n");
        VASSERT(0);
    }
    return thee->center[1];

}
VPUBLIC double Valist_getCenterZ(Valist *thee) {

    if (thee == NULL) {
        Vnm_print(2, "Valist_getCenterZ:  Found null pointer when getting the center of Z coordinate!\n");
        VASSERT(0);
    }
    return thee->center[2];

}

VPUBLIC Vatom* Valist_getAtomList(Valist *thee) {

    if (thee == NULL) {
        Vnm_print(2, "Valist_getAtomList:  Found null pointer when getting the atom list!\n");
        VASSERT(0);
    }
    return thee->atoms;

}

VPUBLIC int Valist_getNumberAtoms(Valist *thee) {

    if (thee == NULL) {
        Vnm_print(2, "Valist_getNumberAtoms:  Found null pointer when getting the number of atoms!\n");
        VASSERT(0);
    }
    return thee->number;

}

VPUBLIC Vatom* Valist_getAtom(Valist *thee, int i) {

    if (thee == NULL) {
        Vnm_print(2, "Valist_getAtom:  Found null pointer when getting atoms!\n");
        VASSERT(0);
    }
    if (i >= thee->number) {
        Vnm_print(2, "Valist_getAtom:  Requested atom number (%d) outside of atom list range (%d)!\n", i, thee->number);
        VASSERT(0);
    }
    return &(thee->atoms[i]);

}

VPUBLIC unsigned long int Valist_memChk(Valist *thee) {

    if (thee == NULL) return 0;
    return Vmem_bytes(thee->vmem);

}

#endif /* if !defined(VINLINE_VATOM) */

VPUBLIC Valist* Valist_ctor() {

    /* Set up the structure */
    Valist *thee = VNULL;
    thee = (Valist*)Vmem_malloc(VNULL, 1, sizeof(Valist));
    if ( thee == VNULL) {
        Vnm_print(2, "Valist_ctor:  Got NULL pointer when constructing the atom list object!\n");
        VASSERT(0);
    }
    if ( Valist_ctor2(thee) != VRC_SUCCESS) {
        Vnm_print(2, "Valist_ctor:   Error in constructing the atom list object!\n");
        VASSERT(0);
    }

    return thee;
}

VPUBLIC Vrc_Codes Valist_ctor2(Valist *thee) {

    thee->atoms = VNULL;
    thee->number = 0;

    /* Initialize the memory management object */
    thee->vmem = Vmem_ctor("APBS:VALIST");

    return VRC_SUCCESS;

}

VPUBLIC void Valist_dtor(Valist **thee)
{
    if ((*thee) != VNULL) {
        Valist_dtor2(*thee);
        Vmem_free(VNULL, 1, sizeof(Valist), (void **)thee);
        (*thee) = VNULL;
    }
}

VPUBLIC void Valist_dtor2(Valist *thee) {

    Vmem_free(thee->vmem, thee->number, sizeof(Vatom), (void **)&(thee->atoms));
    thee->atoms = VNULL;
    thee->number = 0;

    Vmem_dtor(&(thee->vmem));
}

/* Read serial number from PDB ATOM/HETATM field */
VPRIVATE Vrc_Codes Valist_readPDBSerial(Valist *thee, Vio *sock, int *serial) {

    char tok[VMAX_BUFSIZE];
    int ti = 0;

    if (Vio_scanf(sock, "%s", tok) != 1) {
        Vnm_print(2, "Valist_readPDB:  Ran out of tokens while parsing serial!\n");
        return VRC_FAILURE;
    }
    if (sscanf(tok, "%d", &ti) != 1) {
        Vnm_print(2, "Valist_readPDB:  Unable to parse serial token (%s) as int!\n",
                tok);
        return VRC_FAILURE;
    }
    *serial = ti;

    return VRC_SUCCESS;
}

/* Read atom name from PDB ATOM/HETATM field */
VPRIVATE Vrc_Codes Valist_readPDBAtomName(Valist *thee, Vio *sock,
        char atomName[VMAX_ARGLEN]) {

    char tok[VMAX_BUFSIZE];

    if (Vio_scanf(sock, "%s", tok) != 1) {
        Vnm_print(2, "Valist_readPDB:  Ran out of tokens while parsing atom name!\n");
        return VRC_FAILURE;
    }
    if (strlen(tok) < VMAX_ARGLEN) strcpy(atomName, tok);
    else {
        Vnm_print(2, "Valist_readPDB:  Atom name (%s) too long!\n", tok);
        return VRC_FAILURE;
    }
    return VRC_SUCCESS;
}

/* Read residue name from PDB ATOM/HETATM field */
VPRIVATE Vrc_Codes Valist_readPDBResidueName(Valist *thee, Vio *sock,
        char resName[VMAX_ARGLEN]) {

    char tok[VMAX_BUFSIZE];

    if (Vio_scanf(sock, "%s", tok) != 1) {
        Vnm_print(2, "Valist_readPDB:  Ran out of tokens while parsing residue name!\n");
        return VRC_FAILURE;
    }
    if (strlen(tok) < VMAX_ARGLEN) strcpy(resName, tok);
    else {
        Vnm_print(2, "Valist_readPDB:  Residue name (%s) too long!\n", tok);
        return VRC_FAILURE;
    }
    return VRC_SUCCESS;
}

/* Read residue number from PDB ATOM/HETATM field */
VPRIVATE Vrc_Codes Valist_readPDBResidueNumber(
        Valist *thee, Vio *sock, int *resSeq) {

    char tok[VMAX_BUFSIZE];
    char *resstring;
    int ti = 0;

    if (Vio_scanf(sock, "%s", tok) != 1) {
        Vnm_print(2, "Valist_readPDB:  Ran out of tokens while parsing resSeq!\n");
        return VRC_FAILURE;
    }
    if (sscanf(tok, "%d", &ti) != 1) {

        /* One of three things can happen here:
            1)  There is a chainID in the line:    THR A   1
            2)  The chainID is merged with resSeq: THR A1001
            3)  An actual error:                   THR foo

        */

        if (strlen(tok) == 1) {
            /* Case 1: Chain ID Present
                       Read the next field and hope its a float */

            if (Vio_scanf(sock, "%s", tok) != 1) {
                Vnm_print(2, "Valist_readPDB:  Ran out of tokens while parsing resSeq!\n");
                return VRC_FAILURE;
            }
            if (sscanf(tok, "%d", &ti) != 1) {
                Vnm_print(2, "Valist_readPDB:  Unable to parse resSeq token (%s) as int!\n",
                tok);
                return VRC_FAILURE;
            }

        } else {
            /* Case 2: Chain ID, merged string.
                       Move pointer forward past the chainID and check
            */
            //strcpy(resstring, tok);
            resstring = tok;
            resstring++;

            if (sscanf(resstring, "%d", &ti) != 1) {
                /* Case 3:  More than one non-numeral char is present. Error.*/
                Vnm_print(2, "Valist_readPDB:  Unable to parse resSeq token (%s) as int!\n",
                resstring);
                return VRC_FAILURE;
            }
        }
    }
    *resSeq = ti;

    return VRC_SUCCESS;
}

/* Read atom coordinate from PDB ATOM/HETATM field */
VPRIVATE Vrc_Codes Valist_readPDBAtomCoord(Valist *thee, Vio *sock, double *coord) {

    char tok[VMAX_BUFSIZE];
    double tf = 0;

    if (Vio_scanf(sock, "%s", tok) != 1) {
        Vnm_print(2, "Valist_readPDB:  Ran out of tokens while parsing atom coordinate!\n");
        return VRC_FAILURE;
    }
    if (sscanf(tok, "%lf", &tf) != 1) {
        return VRC_FAILURE;
    }
    *coord = tf;

    return VRC_SUCCESS;
}

/* Read charge and radius from PQR ATOM/HETATM field */
VPRIVATE Vrc_Codes Valist_readPDBChargeRadius(Valist *thee, Vio *sock,
        double *charge, double *radius) {

    char tok[VMAX_BUFSIZE];
    double tf = 0;

    if (Vio_scanf(sock, "%s", tok) != 1) {
        Vnm_print(2, "Valist_readPQR:  Ran out of tokens while parsing charge!\n");
        return VRC_FAILURE;
    }
    if (sscanf(tok, "%lf", &tf) != 1) {
        return VRC_FAILURE;
    }
    *charge = tf;

    if (Vio_scanf(sock, "%s", tok) != 1) {
        Vnm_print(2, "Valist_readPQR:  Ran out of tokens while parsing radius!\n");
        return VRC_FAILURE;
    }
    if (sscanf(tok, "%lf", &tf) != 1) {
        return VRC_FAILURE;
    }
    *radius = tf;

    return VRC_SUCCESS;
}

/* Read ATOM/HETATM field of PDB through the X/Y/Z fields */
VPRIVATE Vrc_Codes Valist_readPDB_throughXYZ(
        Valist *thee,
        Vio *sock, /* Socket ready for reading */
        int *serial, /* Set to atom number */
        char atomName[VMAX_ARGLEN], /* Set to atom name */
        char resName[VMAX_ARGLEN], /* Set to residue name */
        int *resSeq, /* Set to residue number */
        double *x, /* Set to x-coordinate */
        double *y, /* Set to y-coordinate */
        double *z  /* Set to z-coordinate */
        ) {


    int i, njunk, gotit;

    /* Grab serial */
    if (Valist_readPDBSerial(thee, sock, serial) == VRC_FAILURE) {
        Vnm_print(2, "Valist_readPDB:  Error while parsing serial!\n");
    }

    /* Grab atom name */
    if (Valist_readPDBAtomName(thee, sock, atomName) == VRC_FAILURE) {
        Vnm_print(2, "Valist_readPDB:  Error while parsing atom name!\n");
        return VRC_FAILURE;
    }

    /* Grab residue name */
    if (Valist_readPDBResidueName(thee, sock, resName) == VRC_FAILURE) {
        Vnm_print(2, "Valist_readPDB:  Error while parsing residue name!\n");
        return VRC_FAILURE;
    }


    /* Grab residue number */
    if (Valist_readPDBResidueNumber(thee, sock, resSeq) == VRC_FAILURE) {
        Vnm_print(2, "Valist_readPDB:  Error while parsing residue number!\n");
        return VRC_FAILURE;
    }


    /* Read tokens until we find one that can be parsed as an atom
     * x-coordinate.  We will allow njunk=1 intervening field that
     * cannot be parsed as a coordinate */
    njunk = 1;
    gotit = 0;
    for (i=0; i<(njunk+1); i++) {
        if (Valist_readPDBAtomCoord(thee, sock, x) == VRC_SUCCESS) {
            gotit = 1;
            break;
        }
    }
    if (!gotit) {
        Vnm_print(2, "Valist_readPDB:  Can't find x!\n");
        return VRC_FAILURE;
    }
    /* Read y-coordinate */
    if (Valist_readPDBAtomCoord(thee, sock, y) == VRC_FAILURE) {
        Vnm_print(2, "Valist_readPDB:  Can't find y!\n");
        return VRC_FAILURE;
    }
    /* Read z-coordinate */
    if (Valist_readPDBAtomCoord(thee, sock, z) == VRC_FAILURE) {
        Vnm_print(2, "Valist_readPDB:  Can't find z!\n");
        return VRC_FAILURE;
    }

#if 0 /* Set to 1 if you want to debug */
    Vnm_print(1, "Valist_readPDB:  serial = %d\n", *serial);
    Vnm_print(1, "Valist_readPDB:  atomName = %s\n", atomName);
    Vnm_print(1, "Valist_readPDB:  resName = %s\n", resName);
    Vnm_print(1, "Valist_readPDB:  resSeq = %d\n", *resSeq);
    Vnm_print(1, "Valist_readPDB:  pos = (%g, %g, %g)\n",
              *x, *y, *z);
#endif

    return VRC_SUCCESS;
}

/* Get a the next available atom storage location, increasing the storage
 * space if necessary.  Return VNULL if something goes wrong. */
VPRIVATE Vatom* Valist_getAtomStorage(
        Valist *thee,
        Vatom **plist, /* Pointer to existing list of atoms */
        int *pnlist, /* Size of existing list, may be changed */
        int *pnatoms /* Existing number of atoms in list; incremented
                       before exit */
        ) {

    Vatom *oldList, *newList, *theList;
    Vatom *oldAtom, *newAtom;
    int iatom, inext, oldLength, newLength, natoms;

    newList = VNULL;

    /* See if we need more space */
    if (*pnatoms >= *pnlist) {

        /* Double the storage space */
        oldLength = *pnlist;
        newLength = 2*oldLength;
        newList = (Vatom*)Vmem_malloc(thee->vmem, newLength, sizeof(Vatom));
        oldList = *plist;

        /* Check the allocation */
        if (newList == VNULL) {
            Vnm_print(2, "Valist_readPDB:  failed to allocate space for %d (Vatom)s!\n", newLength);
            return VNULL;
        }

        /* Copy the atoms over */
        natoms = *pnatoms;
        for (iatom=0; iatom<natoms; iatom++) {
            oldAtom = &(oldList[iatom]);
            newAtom = &(newList[iatom]);
            Vatom_copyTo(oldAtom, newAtom);
            Vatom_dtor2(oldAtom);
        }

        /* Free the old list */
        Vmem_free(thee->vmem, oldLength, sizeof(Vatom), (void **)plist);

        /* Copy new list to plist */
        *plist = newList;
        *pnlist = newLength;
    }

    theList = *plist;
    inext = *pnatoms;

    /* Get the next available spot and increment counters */
    newAtom = &(theList[inext]);
    *pnatoms = inext + 1;

    return newAtom;
}

VPRIVATE Vrc_Codes Valist_setAtomArray(Valist *thee,
        Vatom **plist, /* Pointer to list of atoms to store */
        int nlist, /* Length of list */
        int natoms /* Number of real atom entries in list */
        ) {

    Vatom *list, *newAtom, *oldAtom;
    int i;

    list = *plist;

    /* Allocate necessary space */
    thee->number = 0;
    thee->atoms = (Vatom*)Vmem_malloc(thee->vmem, natoms, sizeof(Vatom));
    if (thee->atoms == VNULL) {
        Vnm_print(2, "Valist_readPDB:  Unable to allocate space for %d (Vatom)s!\n",
                natoms);
        return VRC_FAILURE;
    }
    thee->number = natoms;

    /* Copy over data */
    for (i=0; i<thee->number; i++) {
        newAtom = &(thee->atoms[i]);
        oldAtom = &(list[i]);
        Vatom_copyTo(oldAtom, newAtom);
        Vatom_dtor2(oldAtom);
    }

    /* Free old array */
    Vmem_free(thee->vmem, nlist, sizeof(Vatom), (void **)plist);

    return VRC_SUCCESS;
}

VPUBLIC Vrc_Codes Valist_readPDB(Valist *thee, Vparam *param, Vio *sock) {

    /* WE DO NOT DIRECTLY CONFORM TO PDB STANDARDS -- TO ALLOW LARGER FILES, WE
     * REQUIRE ALL FIELDS TO BE WHITESPACE DELIMITED */

    Vatom *atoms = VNULL;
    Vatom *nextAtom = VNULL;
    Vparam_AtomData *atomData = VNULL;

    char tok[VMAX_BUFSIZE];
    char atomName[VMAX_ARGLEN], resName[VMAX_ARGLEN];

    int nlist, natoms, serial, resSeq;

    double x, y, z, charge, radius, epsilon;
    double pos[3];

    if (thee == VNULL) {
        Vnm_print(2, "Valist_readPDB:  Got NULL pointer when reading PDB file!\n");
        VASSERT(0);
    }
    thee->number = 0;

    Vio_setWhiteChars(sock, Valist_whiteChars);
    Vio_setCommChars(sock, Valist_commChars);

    /* Allocate some initial space for the atoms */
    nlist = 200;
    atoms = (Vatom*)Vmem_malloc(thee->vmem, nlist, sizeof(Vatom));

    natoms = 0;
    /* Read until we run out of lines */
    while (Vio_scanf(sock, "%s", tok) == 1) {

        /* Parse only ATOM/HETATOM fields */
        if ((Vstring_strcasecmp(tok, "ATOM") == 0) ||
            (Vstring_strcasecmp(tok, "HETATM") == 0)) {

            /* Read ATOM/HETATM field of PDB through the X/Y/Z fields */
            if (Valist_readPDB_throughXYZ(thee, sock, &serial, atomName,
                        resName, &resSeq, &x, &y, &z) == VRC_FAILURE) {
                Vnm_print(2, "Valist_readPDB:  Error parsing atom %d!\n",
                          serial);
                return VRC_FAILURE;
            }

            /* Try to find the parameters. */
            atomData = Vparam_getAtomData(param, resName, atomName);
            if (atomData == VNULL) {
                Vnm_print(2, "Valist_readPDB:  Couldn't find parameters for \
atom = %s, residue = %s\n", atomName, resName);
                return VRC_FAILURE;
            }
            charge = atomData->charge;
            radius = atomData->radius;
            epsilon = atomData->epsilon;

            /* Get pointer to next available atom position */
            nextAtom = Valist_getAtomStorage(thee, &atoms, &nlist, &natoms);
            if (nextAtom == VNULL) {
                Vnm_print(2, "Valist_readPDB:  Error in allocating spacing for atoms!\n");
                return VRC_FAILURE;
            }

            /* Store the information */
            pos[0] = x; pos[1] = y; pos[2] = z;
            Vatom_setPosition(nextAtom, pos);
            Vatom_setCharge(nextAtom, charge);
            Vatom_setRadius(nextAtom, radius);
            Vatom_setEpsilon(nextAtom, epsilon);
            Vatom_setAtomID(nextAtom, natoms-1);
            Vatom_setResName(nextAtom, resName);
            Vatom_setAtomName(nextAtom, atomName);

        } /* if ATOM or HETATM */
    } /* while we haven't run out of tokens */

    Vnm_print(0, "Valist_readPDB: Counted %d atoms\n", natoms);
    fflush(stdout);

    /* Store atoms internally */
    if (Valist_setAtomArray(thee, &atoms, nlist, natoms) == VRC_FAILURE) {
        Vnm_print(2, "Valist_readPDB:  unable to store atoms!\n");
        return VRC_FAILURE;
    }

    return Valist_getStatistics(thee);


}

VPUBLIC Vrc_Codes Valist_readPQR(Valist *thee, Vparam *params, Vio *sock) {

    /* WE DO NOT DIRECTLY CONFORM TO PDB STANDARDS -- TO ALLOW LARGER FILES, WE
     * REQUIRE ALL FIELDS TO BE WHITESPACE DELIMITED */


    Vatom *atoms = VNULL;
    Vatom *nextAtom = VNULL;
    Vparam_AtomData *atomData = VNULL;

    char tok[VMAX_BUFSIZE];
    char atomName[VMAX_ARGLEN], resName[VMAX_ARGLEN];
    char chs[VMAX_BUFSIZE];
    char c, *ch;

    int use_params = 0;
    int nlist, natoms, serial, resSeq;

    double x, y, z, charge, radius, epsilon;
    double pos[3];

    epsilon = 0.0;
    c = 'a';

    if (thee == VNULL) {
        Vnm_print(2, "Valist_readPQR:  Got NULL pointer when reading PQR file!\n");
        VASSERT(0);
    }
    thee->number = 0;

    Vio_setWhiteChars(sock, Valist_whiteChars);
    Vio_setCommChars(sock, Valist_commChars);

    /* Allocate some initial space for the atoms */
    nlist = 200;
    atoms = (Vatom*)Vmem_malloc(thee->vmem, nlist, sizeof(Vatom));

    /* Check if we are using a parameter file or not */
    if(params != VNULL) use_params = 1;

    natoms = 0;
    /* Read until we run out of lines */
    while (Vio_scanf(sock, "%s", tok) == 1) {

        /* Parse only ATOM/HETATOM fields */
        if ((Vstring_strcasecmp(tok, "ATOM") == 0) ||
            (Vstring_strcasecmp(tok, "HETATM") == 0)) {

            /* Read ATOM/HETATM field of PDB through the X/Y/Z fields */
            if (Valist_readPDB_throughXYZ(thee, sock, &serial, atomName,
                        resName, &resSeq, &x, &y, &z) == VRC_FAILURE) {
                Vnm_print(2, "Valist_readPQR:  Error parsing atom %d!\n",serial);
                Vnm_print(2, "Please double check this atom in the pqr file, e.g., make sure there are no concatenated fields.\n");
                return VRC_FAILURE;
            }

            /* Read Q/R fields */
            if (Valist_readPDBChargeRadius(thee, sock, &charge, &radius) == VRC_FAILURE) {
                Vnm_print(2, "Valist_readPQR:  Error parsing atom %d!\n",
                          serial);
                Vnm_print(2, "Please double check this atom in the pqr file, e.g., make sure there are no concatenated fields.\n");
                return VRC_FAILURE;
            }

            if(use_params){
                /* Try to find the parameters. */
                atomData = Vparam_getAtomData(params, resName, atomName);
                if (atomData == VNULL) {
                    Vnm_print(2, "Valist_readPDB:  Couldn't find parameters for \
atom = %s, residue = %s\n", atomName, resName);
                    return VRC_FAILURE;
                }
                charge = atomData->charge;
                radius = atomData->radius;
                epsilon = atomData->epsilon;
            }

            /* Get pointer to next available atom position */
            nextAtom = Valist_getAtomStorage(thee, &atoms, &nlist, &natoms);
            if (nextAtom == VNULL) {
                Vnm_print(2, "Valist_readPQR:  Error in allocating spacing for atoms!\n");
                return VRC_FAILURE;
            }

            /* Store the information */
            pos[0] = x; pos[1] = y; pos[2] = z;
            Vatom_setPosition(nextAtom, pos);
            Vatom_setCharge(nextAtom, charge);
            Vatom_setRadius(nextAtom, radius);
            Vatom_setEpsilon(nextAtom, epsilon);
            Vatom_setAtomID(nextAtom, natoms-1);
            Vatom_setResName(nextAtom, resName);
            Vatom_setAtomName(nextAtom, atomName);

        } /* if ATOM or HETATM */
        //if the line doesn't start with ATOM or HETATM, then throw it away
        else {
            //if we didn't see ATOM or HETATM then skip to the next line
            while (Vio_getc(sock) != '\n');
        }
    } /* while we haven't run out of tokens */

    Vnm_print(0, "Valist_readPQR: Counted %d atoms\n", natoms);
    fflush(stdout);

    /* Store atoms internally */
    if (Valist_setAtomArray(thee, &atoms, nlist, natoms) == VRC_FAILURE) {
        Vnm_print(2, "Valist_readPDB:  unable to store atoms!\n");
        return VRC_FAILURE;
    }

    return Valist_getStatistics(thee);


}

VPUBLIC Vrc_Codes Valist_readXML(Valist *thee, Vparam *params, Vio *sock) {

    Vatom *atoms = VNULL;
    Vatom *nextAtom = VNULL;

    char tok[VMAX_BUFSIZE];
    char endtag[VMAX_BUFSIZE];

    int nlist, natoms;
    int xset, yset, zset, chgset, radset;

    double x, y, z, charge, radius, dtmp;
    double pos[3];

    if (thee == VNULL) {
        Vnm_print(2, "Valist_readXML:  Got NULL pointer when reading XML file!\n");
        VASSERT(0);
    }
    thee->number = 0;

    Vio_setWhiteChars(sock, Valist_xmlwhiteChars);
    Vio_setCommChars(sock, Valist_commChars);

    /* Allocate some initial space for the atoms */
    nlist = 200;
    atoms = (Vatom*)Vmem_malloc(thee->vmem, nlist, sizeof(Vatom));

    /* Initialize some variables */
    natoms = 0;
    xset = 0;
    yset = 0;
    zset = 0;
    chgset = 0;
    radset = 0;
    strcpy(endtag,"/");

    if(params == VNULL){
        Vnm_print(1,"\nValist_readXML: Warning Warning Warning Warning Warning\n");
        Vnm_print(1,"Valist_readXML: The use of XML input files with parameter\n");
        Vnm_print(1,"Valist_readXML: files is currently not supported.\n");
        Vnm_print(1,"Valist_readXML: Warning Warning Warning Warning Warning\n\n");
    }

    /* Read until we run out of lines */
    while (Vio_scanf(sock, "%s", tok) == 1) {

        /* The first tag taken is the start tag - save it to detect end */
        if (Vstring_strcasecmp(endtag, "/") == 0) strcat(endtag, tok);

        if (Vstring_strcasecmp(tok, "x") == 0) {
            Vio_scanf(sock, "%s", tok);
            if (sscanf(tok, "%lf", &dtmp) != 1) {
                Vnm_print(2, "Valist_readXML:  Unexpected token (%s) while \
reading x!\n", tok);
                  return VRC_FAILURE;
              }
            x = dtmp;
            xset = 1;
        } else if (Vstring_strcasecmp(tok, "y") == 0) {
            Vio_scanf(sock, "%s", tok);
            if (sscanf(tok, "%lf", &dtmp) != 1) {
                Vnm_print(2, "Valist_readXML:  Unexpected token (%s) while \
reading y!\n", tok);
                  return VRC_FAILURE;
              }
            y = dtmp;
            yset = 1;
        } else if (Vstring_strcasecmp(tok, "z") == 0) {
            Vio_scanf(sock, "%s", tok);
            if (sscanf(tok, "%lf", &dtmp) != 1) {
                Vnm_print(2, "Valist_readXML:  Unexpected token (%s) while \
reading z!\n", tok);
                  return VRC_FAILURE;
              }
            z = dtmp;
            zset = 1;
        } else if (Vstring_strcasecmp(tok, "charge") == 0) {
            Vio_scanf(sock, "%s", tok);
            if (sscanf(tok, "%lf", &dtmp) != 1) {
                Vnm_print(2, "Valist_readXML:  Unexpected token (%s) while \
reading charge!\n", tok);
                  return VRC_FAILURE;
              }
            charge = dtmp;
            chgset = 1;
        } else if (Vstring_strcasecmp(tok, "radius") == 0) {
            Vio_scanf(sock, "%s", tok);
            if (sscanf(tok, "%lf", &dtmp) != 1) {
                Vnm_print(2, "Valist_readXML:  Unexpected token (%s) while \
reading radius!\n", tok);
                  return VRC_FAILURE;
              }
            radius = dtmp;
            radset = 1;
        } else if (Vstring_strcasecmp(tok, "/atom") == 0) {

          /* Get pointer to next available atom position */
            nextAtom = Valist_getAtomStorage(thee, &atoms, &nlist, &natoms);
            if (nextAtom == VNULL) {
                Vnm_print(2, "Valist_readXML:  Error in allocating spacing for atoms!\n");
                return VRC_FAILURE;
            }

            if (xset && yset && zset && chgset && radset){

                /* Store the information */
                pos[0] = x; pos[1] = y; pos[2] = z;
                Vatom_setPosition(nextAtom, pos);
                Vatom_setCharge(nextAtom, charge);
                Vatom_setRadius(nextAtom, radius);
                Vatom_setAtomID(nextAtom, natoms-1);

                /* Reset the necessary flags */
                xset = 0;
                yset = 0;
                zset = 0;
                chgset = 0;
                radset = 0;
            } else {
                Vnm_print(2,  "Valist_readXML:  Missing field(s) in atom tag:\n");
                if (!xset) Vnm_print(2,"\tx value not set!\n");
                if (!yset) Vnm_print(2,"\ty value not set!\n");
                if (!zset) Vnm_print(2,"\tz value not set!\n");
                if (!chgset) Vnm_print(2,"\tcharge value not set!\n");
                if (!radset) Vnm_print(2,"\tradius value not set!\n");
                return VRC_FAILURE;
            }
        } else if (Vstring_strcasecmp(tok, endtag) == 0) break;
    }

    Vnm_print(0, "Valist_readXML: Counted %d atoms\n", natoms);
    fflush(stdout);

    /* Store atoms internally */
    if (Valist_setAtomArray(thee, &atoms, nlist, natoms) == VRC_FAILURE) {
        Vnm_print(2, "Valist_readXML:  unable to store atoms!\n");
        return VRC_FAILURE;
    }

    return Valist_getStatistics(thee);

}

/* Load up Valist with various statistics */
VPUBLIC Vrc_Codes Valist_getStatistics(Valist *thee) {

    Vatom *atom;
    int i, j;

    if (thee == VNULL) {
        Vnm_print(2, "Valist_getStatistics:  Got NULL pointer when loading up Valist with various statistics!\n");
        VASSERT(0);
    }

    thee->center[0] = 0.;
    thee->center[1] = 0.;
    thee->center[2] = 0.;
    thee->maxrad = 0.;
    thee->charge = 0.;

    if (thee->number == 0) return VRC_FAILURE;

    /* Reset stat variables */
    atom = &(thee->atoms[0]);
    for (i=0; i<3; i++) {
        thee->maxcrd[i] = thee->mincrd[i] = atom->position[i];
    }
    thee->maxrad = atom->radius;
    thee->charge = 0.0;

    for (i=0; i<thee->number; i++) {

        atom = &(thee->atoms[i]);
        for (j=0; j<3; j++) {
            if (atom->position[j] < thee->mincrd[j])
              thee->mincrd[j] = atom->position[j];
            if (atom->position[j] > thee->maxcrd[j])
              thee->maxcrd[j] = atom->position[j];
        }
        if (atom->radius > thee->maxrad) thee->maxrad = atom->radius;
        thee->charge = thee->charge + atom->charge;
    }

    thee->center[0] = 0.5*(thee->maxcrd[0] + thee->mincrd[0]);
    thee->center[1] = 0.5*(thee->maxcrd[1] + thee->mincrd[1]);
    thee->center[2] = 0.5*(thee->maxcrd[2] + thee->mincrd[2]);

    Vnm_print(0, "Valist_getStatistics:  Max atom coordinate:  (%g, %g, %g)\n",
              thee->maxcrd[0], thee->maxcrd[1], thee->maxcrd[2]);
    Vnm_print(0, "Valist_getStatistics:  Min atom coordinate:  (%g, %g, %g)\n",
              thee->mincrd[0], thee->mincrd[1], thee->mincrd[2]);
    Vnm_print(0, "Valist_getStatistics:  Molecule center:  (%g, %g, %g)\n",
              thee->center[0], thee->center[1], thee->center[2]);

    return VRC_SUCCESS;
}
