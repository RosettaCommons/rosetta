/**
 * @defgroup  Frontend  High-level front-end routines
 */

/**
 *  @file    routines.h
 *  @author  Nathan Baker
 *  @brief   Header file for front end auxiliary routines
 *  @ingroup  Frontend
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

#ifndef _APBSROUTINES_H_
#define _APBSROUTINES_H_

#include "apbs.h"

#ifdef HAVE_MC_H
#    include "mc/mc.h"
#    include "apbs/vfetk.h"
#endif
#ifdef HAVE_MCX_H
#    include "mcx/mcx.h"
#endif

/**
 * @brief  Return code for APBS during failure
 * @ingroup  Frontend */
#define APBSRC 13

/**
 * @brief  Structure to hold atomic forces
 * @ingroup  Frontend
 * @author  Nathan Baker */
struct AtomForce {
   double ibForce[3];  /**< Ion-boundary force */
   double qfForce[3];  /**< Charge-field force */
   double dbForce[3];  /**< Dielectric boundary force */
   double sasaForce[3];  /**< SASA force (coupled to gamma) */
   double savForce[3];  /**< SAV force (coupled to press) */
   double wcaForce[3];  /**< WCA integral force (coupled to bconc) */
};

/**
 * @brief  Define AtomForce type
 * @ingroup  Frontend */
typedef struct AtomForce AtomForce;

/**
 * @brief  Loads and returns parameter object
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @returns  Pointer to parameter object or NULL */
VEXTERNC Vparam* loadParameter(
                               NOsh *nosh  /**< Pointer to NOsh object with input
                               file information */
                               );

/**
 * @brief  Load the molecules given in NOsh into atom lists
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @returns  1 if successful, 0 otherwise */
VEXTERNC int loadMolecules(
                           NOsh *nosh, /**< NOsh object with input file information */
                           Vparam *param,  /**< NULL (if PQR files only) or pointer
                           to parameter object */
                           Valist *alist[NOSH_MAXMOL]  /**< List of atom list objects
                           (to be populated) */
                           );

/**
 * @brief  Destroy the loaded molecules
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  nosh  NOsh object with input file information
 * @param  alist  List of atom list objects */
VEXTERNC void killMolecules(NOsh *nosh, Valist *alist[NOSH_MAXMOL]);

/**
 * @brief  Load the dielectric maps given in NOsh into grid objects
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  nosh  NOsh object with input file information
 * @param  dielXMap  List of x-shifted dielectric maps
 * @param  dielYMap  List of y-shifted dielectric maps
 * @param  dielZMap  List of x-shifted dielectric maps
 * @returns  1 if successful, 0 otherwise */
VEXTERNC int loadDielMaps(NOsh *nosh,
                          Vgrid *dielXMap[NOSH_MAXMOL],
                          Vgrid *dielYMap[NOSH_MAXMOL],
                          Vgrid *dielZMap[NOSH_MAXMOL]
                         );

/**
 * @brief  Destroy the loaded dielectric
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  nosh  NOsh object with input file information
 * @param  dielXMap  List of x-shifted dielectric maps
 * @param  dielYMap  List of y-shifted dielectric maps
 * @param  dielZMap  List of x-shifted dielectric maps */
VEXTERNC void killDielMaps(NOsh *nosh, Vgrid *dielXMap[NOSH_MAXMOL],
  Vgrid *dielYMap[NOSH_MAXMOL], Vgrid *dielZMap[NOSH_MAXMOL]);

/**
 * @brief  Load the kappa maps given in NOsh into grid objects
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  nosh  NOsh object with input file information
 * @param  kappa  List of kappa maps
 * @returns  1 if successful, 0 otherwise */
VEXTERNC int loadKappaMaps(NOsh *nosh, Vgrid *kappa[NOSH_MAXMOL]);

/**
 * @brief  Destroy the loaded kappa maps
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  nosh  NOsh object with input file information
 * @param  kappa  List of kappa maps */
VEXTERNC void killKappaMaps(NOsh *nosh, Vgrid *kappa[NOSH_MAXMOL]);

/**
 * @brief  Load the potential maps given in NOsh into grid objects
 * @ingroup  Frontend
 * @author  David Gohara
 * @param  nosh  NOsh object with input file information
 * @param  pot  List of potential maps
 * @returns  1 if successful, 0 otherwise */
VEXTERNC int loadPotMaps(NOsh *nosh, Vgrid *pot[NOSH_MAXMOL]);

/**
 * @brief  Destroy the loaded potential maps
 * @ingroup  Frontend
 * @author  David Gohara
 * @param  nosh  NOsh object with input file information
 * @param  pot  List of potential maps */
VEXTERNC void killPotMaps(NOsh *nosh, Vgrid *pot[NOSH_MAXMOL]);

/**
 * @brief  Load the charge maps given in NOsh into grid objects
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  nosh  NOsh object with input file information
 * @param  charge  List of kappa maps
 * @returns  1 if successful, 0 otherwise */
VEXTERNC int loadChargeMaps(NOsh *nosh, Vgrid *charge[NOSH_MAXMOL]);

/**
 * @brief  Destroy the loaded charge maps
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  nosh  NOsh object with input file information
 * @param  charge  List of charge maps */
VEXTERNC void killChargeMaps(NOsh *nosh, Vgrid *charge[NOSH_MAXMOL]);

/**
 * @brief  Print out generic PBE params loaded from input
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  pbeparm  PBEparm object */
VEXTERNC void printPBEPARM(PBEparm *pbeparm);

/**
 * @brief  Print out MG-specific params loaded from input
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  realCenter  Center of mesh for actual calculation
 * @param  mgparm  MGparm object */
VEXTERNC void printMGPARM(MGparm *mgparm, double realCenter[3]);

/**
 * @brief  Initialize an MG calculation
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @return  1 if succesful, 0 otherwise */
VEXTERNC int initMG(
                    int icalc,  /**< Index of calculation in pmg/pmpg arrays */
                    NOsh *nosh,  /**< Object with parsed input file parameters */
                    MGparm *mgparm,  /**< Object with MG-specific parameters */
                    PBEparm *pbeparm,  /**< Object with generic PBE parameters  */
                    double realCenter[3],  /**< The actual center of the current mesh */
                    Vpbe *pbe[NOSH_MAXCALC],  /**< Array of Vpbe objects (one for each calc) */
                    Valist *alist[NOSH_MAXMOL],  /**< Array of atom lists */
                    Vgrid *dielXMap[NOSH_MAXMOL],  /**< Array of x-shifted dielectric maps */
                    Vgrid *dielYMap[NOSH_MAXMOL],  /**< Array of y-shifted dielectric maps */
                    Vgrid *dielZMap[NOSH_MAXMOL],  /**< Array of z-shifted dielectric maps */
                    Vgrid *kappaMap[NOSH_MAXMOL],  /**< Array of kappa maps  */
                    Vgrid *chargeMap[NOSH_MAXMOL],  /**< Array of charge maps */
                    Vpmgp *pmgp[NOSH_MAXCALC],  /**< Array of MG parameter objects (one for each calc) */
                    Vpmg *pmg[NOSH_MAXCALC],  /**< Array of MG objects (one for each calc) */
                    Vgrid *potMap[NOSH_MAXMOL]  /**< Array of potential maps  */
                    );

/**
 * @brief  Kill structures initialized during an MG calculation
 * @ingroup  Frontend
 * @author  Nathan Baker
 */
VEXTERNC void killMG(
                     NOsh *nosh,  /** Object with parsed input file parameters */
                     Vpbe *pbe[NOSH_MAXCALC],  /** Array of Vpbe objects for each calc */
                     Vpmgp *pmgp[NOSH_MAXCALC],  /** Array of MG parameter objects for each calc */
                     Vpmg *pmg[NOSH_MAXCALC]  /** Array of MG objects for each calc */
);

/**
 * @brief  Solve the PBE with MG
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param nosh  Object with parsed input file parameters
 * @param pmg  MG objects for this calculation
 * @param type  Type of MG calculation
 * @return  1 if successful, 0 otherwise */
VEXTERNC int solveMG(NOsh *nosh, Vpmg *pmg, MGparm_CalcType type);

/**
 * @brief  Set MG partitions for calculating observables and performing I/O
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param nosh  Object with parsed input file parameters
 * @param mgparm  MG parameters from input file
 * @param pmg  MG object
 * @return  1 if successful, 0 otherwise */
VEXTERNC int setPartMG(NOsh *nosh, MGparm *mgparm, Vpmg *pmg);

/**
 * @brief  Calculate electrostatic energies from MG solution
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param nosh  Object with parsed input file parameters
 * @param icalc  Index of calculation
 * @param pmg  MG object
 * @param nenergy  Set to number of entries in energy arrays
 * @param totEnergy  Set to total energy (in kT)
 * @param qfEnergy  Set to charge-potential energy (in kT)
 * @param qmEnergy  Set to mobile ion energy (in kT)
 * @param dielEnergy  Set to polarization energy (in kT)
 * @return  1 if successful, 0 otherwise */
VEXTERNC int energyMG(NOsh* nosh, int icalc, Vpmg *pmg,
  int *nenergy, double *totEnergy, double *qfEnergy, double *qmEnergy,
  double *dielEnergy);

/**
 * @brief  Kill arrays allocated for energies
 * @ingroup  Frontend
 * @author  Nathan Baker */
VEXTERNC void killEnergy();

/**
 * @brief  Calculate forces from MG solution
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  mem  Memory management object
 * @param  nosh  Parameters from input file
 * @param  pbeparm  Generic PBE parameters
 * @param  mgparm  MG-specific parmaeters
 * @param pmg  MG object
 * @param nforce  Set to number of forces in arrays
 * @param  atomForce  List of atom forces
 * @param  alist  List of atom lists
 * @return  1 if successful, 0 otherwise */
VEXTERNC int forceMG(Vmem *mem, NOsh *nosh, PBEparm *pbeparm,  MGparm *mgparm,
  Vpmg *pmg, int *nforce, AtomForce **atomForce, Valist *alist[NOSH_MAXMOL]);

/**
 * @brief  Free memory from MG force calculation
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  mem  Memory management object
 * @param  nosh  Parameters from input file
 * @param  nforce  Number of forces in arrays
 * @param  atomForce  List of atom forces */
VEXTERNC void killForce(Vmem *mem, NOsh *nosh, int nforce[NOSH_MAXCALC],
  AtomForce *atomForce[NOSH_MAXCALC]);

/**
 * @brief Store energy in arrays for future use
 * @ingroup Frontend
 * @author Todd Dolinsky */
VEXTERNC void storeAtomEnergy(
        Vpmg *pmg, /**< MG object */
        int icalc, /**< Calculation number */
        double **atomEnergy, /**< Pointer to storage array of doubles */
        int *nenergy /**< Stores number of atoms per calc */
        );

/**
 * @brief Write out information to a flat file
 * @ingroup Frontend
 * @author  Todd Dolinsky
 * @param nosh  Parameters from input file
 * @param com   The communications object
 * @param fname The target XML file name
 * @param totEnergy An array with per-calc total energies (in kT)
 * @param qfEnergy  An array with per-calc charge-potential energies (in kT)
 * @param qmEnergy  An array with per-calc mobile energies (in kT)
 * @param dielEnergy  An array with per-calc polarization energies (in kT)
 * @param nenergy  An array containing the number of atoms per-calc
 * @param atomEnergy An array containing per-atom energies (in KT) per calc
 * @param nforce  An array containing the number of forces calculated per-calc
 * @param atomForce An array containing per-atom forces per calc
 * @return 1 if successful, 0 otherwise */
VEXTERNC int writedataFlat(NOsh *nosh, Vcom *com, const char *fname,
  double totEnergy[NOSH_MAXCALC], double qfEnergy[NOSH_MAXCALC],
  double qmEnergy[NOSH_MAXCALC], double dielEnergy[NOSH_MAXCALC],
  int nenergy[NOSH_MAXCALC], double *atomEnergy[NOSH_MAXCALC],
  int nforce[NOSH_MAXCALC], AtomForce *atomForce[NOSH_MAXCALC]);

/**
 * @brief Write out information to an XML file
 * @ingroup Frontend
 * @author  Todd Dolinsky
 * @param nosh  Parameters from input file
 * @param com   The communications object
 * @param fname The target XML file name
 * @param totEnergy An array with per-calc total energies (in kT)
 * @param qfEnergy  An array with per-calc charge-potential energies (in kT)
 * @param qmEnergy  An array with per-calc mobile energies (in kT)
 * @param dielEnergy  An array with per-calc polarization energies (in kT)
 * @param nenergy  An array containing the number of atoms per-calc
 * @param atomEnergy An array containing per-atom energies (in KT) per calc
 * @param nforce  An array containing the number of forces calculated per-calc
 * @param atomForce An array containing per-atom forces per calc
 * @return 1 if successful, 0 otherwise */
VEXTERNC int writedataXML(NOsh *nosh, Vcom *com, const char *fname,
  double totEnergy[NOSH_MAXCALC], double qfEnergy[NOSH_MAXCALC],
  double qmEnergy[NOSH_MAXCALC], double dielEnergy[NOSH_MAXCALC],
  int nenergy[NOSH_MAXCALC], double *atomEnergy[NOSH_MAXCALC],
  int nforce[NOSH_MAXCALC], AtomForce *atomForce[NOSH_MAXCALC]);

/**
 * @brief  Write out observables from MG calculation to file
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  rank  Processor rank (if parallel calculation)
 * @param  nosh  Parameters from input file
 * @param  pbeparm  Generic PBE parameters
 * @param pmg  MG object
 * @return  1 if successful, 0 otherwise */
VEXTERNC int writedataMG(int rank, NOsh *nosh, PBEparm *pbeparm, Vpmg *pmg);

/**
 * @brief  Write out operator matrix from MG calculation to file
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  rank  Processor rank (if parallel calculation)
 * @param  nosh  Parameters from input file
 * @param  pbeparm  Generic PBE parameters
 * @param pmg  MG object
 * @return  1 if successful, 0 otherwise */
VEXTERNC int writematMG(int rank, NOsh *nosh, PBEparm *pbeparm, Vpmg *pmg);

/**
* @brief  Access net local energy
 * @ingroup  Frontend
 * @author  Justin Xiang
 * @param  com  Communications object
 * @param  nosh  Parameters from input file
 * @param  totEnergy  Array of energies from different calculations
 * @param  iprint  Index of energy statement to print
 * @return  Net local energy */
VEXTERNC double returnEnergy(Vcom *com, NOsh *nosh, double totEnergy[NOSH_MAXCALC], int iprint);

/**
* @brief  Combine and pretty-print energy data (deprecated...see printElecEnergy)
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @return  1 if successful, 0 otherwise */
VEXTERNC int printEnergy(
                         Vcom *com, /** Communications object */
                         NOsh *nosh, /** Parameters from input file */
                         double totEnergy[NOSH_MAXCALC], /** Array of energies
                         from different calculations */
                         int iprint /** Index of energy statement to print */
                         );

/**
* @brief  Combine and pretty-print energy data
* @ingroup  Frontend
* @author  David Gohara
* @return  1 if successful, 0 otherwise */
VEXTERNC int printElecEnergy(
                         Vcom *com, /** Communications object */
                         NOsh *nosh, /** Parameters from input file */
                         double totEnergy[NOSH_MAXCALC], /** Array of energies
                             from different calculations */
                         int iprint /** Index of energy statement to print */
                         );

/**
* @brief  Combine and pretty-print energy data
* @ingroup  Frontend
* @author  David Gohara
* @return  1 if successful, 0 otherwise */
VEXTERNC int printApolEnergy(
                         NOsh *nosh,  /**< Parameters from input file */
                         int iprint  /**< Index of energy statement to print */
                         );

/**
 * @brief  Combine and pretty-print force data (deprecated...see printElecForce)
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @return  1 if successful, 0 otherwise */
VEXTERNC int printForce(
                        Vcom *com, /** Communications object */
                        NOsh *nosh, /** Parameters from input file */
                        int nforce[NOSH_MAXCALC], /** Number of forces calculated */
                        AtomForce *atomForce[NOSH_MAXCALC], /** Array of force structures */
                        int i /** Index of force statement to print */
                        );

/**
* @brief  Combine and pretty-print force data
* @ingroup  Frontend
* @author  David Gohara
* @return  1 if successful, 0 otherwise */
VEXTERNC int printElecForce(
                            Vcom *com, /** Communications object */
                            NOsh *nosh, /** Parameters from input file */
                            int nforce[NOSH_MAXCALC], /** Number of forces calculated */
                            AtomForce *atomForce[NOSH_MAXCALC], /** Array of force structures */
                            int i /** Index of force statement to print */
                            );

/**
* @brief  Combine and pretty-print force data
* @ingroup  Frontend
* @author  David Gohara
* @return  1 if successful, 0 otherwise */
VEXTERNC int printApolForce(
                            Vcom *com, /** Communications object */
                            NOsh *nosh, /** Parameters from input file */
                            int nforce[NOSH_MAXCALC], /** Number of forces calculated */
                            AtomForce *atomForce[NOSH_MAXCALC], /** Array of force structures */
                            int i /** Index of force statement to print */
                            );

/**
 * @brief  Wrapper to start MALOC Vio layer
 * @ingroup  Frontend
 * @author  Nathan Baker and Robert Konecny */
VEXTERNC void startVio();

/**
 * @brief  Calculate non-polar energies
 * @ingroup  Frontend
 * @author  David Gohara
 * @return  1 if successful, 0 otherwise */
VEXTERNC int energyAPOL(
                        APOLparm *apolparm, /** APOLparm object */
                        double sasa, /** Solvent accessible surface area */
                        double sav, /** Solvent accessible volume */
                        double atomsasa[], /** Array for SASA per atom **/
                        double atomwcaEnergy[], /** Array for WCA energy per atom **/
                        int numatoms /** Number of atoms (or size of the above arrays) **/
                        );

/**
 * @brief  Calculate non-polar forces
 * @ingroup  Frontend
 * @author  David Gohara
 * @return  1 if successful, 0 otherwise */
VEXTERNC int forceAPOL(
                       Vacc *acc,  /**< Accessiblity object */
                       Vmem *mem,  /**< Memory manager */
                       APOLparm *apolparm,  /**< Apolar calculation parameter
                       object */
                       int *nforce,  /**< Number of atomic forces to calculate
                       statements for */
                       AtomForce **atomForce,  /**< Object for storing atom
                       forces */
                       Valist *alist,  /**< Atom list */
                       Vclist *clist  /**< Cell list for accessibility object */
                       );

/**
 * @brief  Upperlevel routine to the non-polar energy and force routines
 * @ingroup  Frontend
 * @author  David Gohara
 * @return  1 if successful, 0 otherwise */
VEXTERNC int initAPOL(
                      NOsh *nosh,  /**< Input parameter object */
                      Vmem *mem,  /**< Memory manager */
                      Vparam *param,  /**< Atom parameters */
                      APOLparm *apolparm,  /**< Apolar calculation parameters */
                      int *nforce,  /**< Number of force calculations */
                      AtomForce **atomForce,  /**< Atom force storage object */
                      Valist *alist  /**< Atom list */
                      );


#ifdef HAVE_MC_H
#include "apbs/vfetk.h"

/**
 * @brief  Print out FE-specific params loaded from input
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  icalc  Calculation index
 * @param  nosh  Master parameter object
 * @param feparm  FE-specific parameters
 * @param fetk  Array of FE solver objects  */
VEXTERNC void printFEPARM(int icalc, NOsh *nosh, FEMparm *feparm,
  Vfetk *fetk[NOSH_MAXCALC]);

/**
 * @brief  Calculate electrostatic energies from FE solution
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param nosh  Object with parsed input file parameters
 * @param icalc  Index of calculation
 * @param fetk  FE object  array
 * @param nenergy  Set to number of entries in energy arrays
 * @param totEnergy  Set to total energy (in kT)
 * @param qfEnergy  Set to charge-potential energy (in kT)
 * @param qmEnergy  Set to mobile ion energy (in kT)
 * @param dielEnergy  Set to polarization energy (in kT)
 * @bug  "calcenergy 2" does not work
 * @return  1 if successful, 0 otherwise */
VEXTERNC int energyFE(NOsh* nosh, int icalc, Vfetk *fetk[NOSH_MAXCALC],
  int *nenergy, double *totEnergy, double *qfEnergy, double *qmEnergy,
  double *dielEnergy);

/**
 * @brief  Initialize FE solver objects
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @bug  THIS FUNCTION IS HARD-CODED TO SOLVE LRPBE
 * @todo  THIS FUNCTION IS HARD-CODED TO SOLVE LRPBE
 */
VEXTERNC Vrc_Codes initFE(
                    int icalc, /** Index in pb, fetk to initialize (calculation index) */
                    NOsh *nosh,  /** Master parmaeter object */
                    FEMparm *feparm,  /** FE-specific parameters */
                    PBEparm *pbeparm,  /** Generic PBE parameters */
                    Vpbe *pbe[NOSH_MAXCALC],  /** Array of PBE objects */
                    Valist *alist[NOSH_MAXMOL],  /** Array of atom lists */
                    Vfetk *fetk[NOSH_MAXCALC]  /** Array of finite element objects */
);

/**
 * @brief  Kill structures initialized during an FE calculation
 * @ingroup  Frontend
 * @author  Nathan Baker
 */
VEXTERNC void killFE(
                     NOsh *nosh,  /** Object with parsed input file parameters */
                     Vpbe *pbe[NOSH_MAXCALC],  /** Array of Vpbe objects for each calc */
                     Vfetk *fetk[NOSH_MAXCALC],  /** Array of FEtk objects for each calc */
                     Gem *gem[NOSH_MAXMOL]  /** Array of geometry manager objects for each calc */
);

/**
 * @brief  Pre-refine mesh before solve
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  i  Calculation index
 * @param  nosh  Master parameter object
 * @param feparm  FE-specific parameters
 * @param fetk  Array of FE solver objects
 * @return  1 if successful, 0 otherwise */
VEXTERNC int preRefineFE(int i,
                        FEMparm *feparm,
                        Vfetk *fetk[NOSH_MAXCALC]
);

/**
 * @brief  Partition mesh (if applicable)
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  i  Calculation index
 * @param  nosh  Master parameter object
 * @param feparm  FE-specific parameters
 * @param fetk  Array of FE solver objects
 * @return  1 if successful, 0 otherwise */
VEXTERNC int partFE(int i, NOsh *nosh, FEMparm *feparm,
  Vfetk *fetk[NOSH_MAXCALC]);

/**
 * @brief  Solve-estimate-refine
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  i  Calculation index
 * @param feparm  FE-specific parameters
 * @param pbeparm  Generic PBE parameters
 * @param fetk  Array of FE solver objects
 * @return  1 if successful, 0 otherwise */
VEXTERNC int solveFE(int i,
                     PBEparm *pbeparm,
                     FEMparm *feparm,
                     Vfetk *fetk[NOSH_MAXCALC]
);

/**
 * @brief  Estimate error, mark mesh, and refine mesh after solve
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  icalc  Calculation index
 * @param feparm  FE-specific parameters
 * @param fetk  Array of FE solver objects
 * @return  1 if successful, 0 otherwise -- note that a 0 will likely imply
 * that either the max number of vertices have been met or no vertices were
 * marked for refinement.  In either case, this should not be treated as a
 * fatal error.  */
VEXTERNC int postRefineFE(int icalc,
                          FEMparm *feparm,
                          Vfetk *fetk[NOSH_MAXCALC]
);

/**
 * @brief  Write FEM data to files
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @param  rank  Rank of processor (for parallel runs)
 * @param  nosh  NOsh object
 * @param  pbeparm  PBEparm object
 * @param  fetk  FEtk object (with solution)
 * @return  1 if successful, 0 otherwise */
VEXTERNC int writedataFE(int rank, NOsh *nosh, PBEparm *pbeparm, Vfetk *fetk);

/**
 * @brief  Load the meshes given in NOsh into geometry objects
 * @ingroup  Frontend
 * @author  Nathan Baker
 * @returns  Error code on success/failure */
VEXTERNC Vrc_Codes loadMeshes(
                            NOsh *nosh, /**< NOsh object with input file information */
                            Gem *gm[NOSH_MAXMOL]  /**< List of geometry objects
                            (to be populated) */
                           );

/**
 * @brief  Destroy the loaded meshes
 * @ingroup  Frontend
 * @author  Nathan Baker */
VEXTERNC void killMeshes(
                            NOsh *nosh, /**< NOsh object with input file information */
                            Gem *alist[NOSH_MAXMOL]  /**< Populated list of geometry objects to be destroyed */
                            );
#endif

#endif
