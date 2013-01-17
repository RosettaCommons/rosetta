/** @defgroup Vpmg Vpmg class
 *  @brief  A wrapper for Mike Holst's PMG multigrid code.
 *  @note   Many of the routines and macros are borrowed from the main.c driver
 *          (written by Mike Holst) provided with the PMG code.
 */

/**
 *  @file     vpmg.h
 *  @ingroup  Vpmg
 *  @brief    Contains declarations for class Vpmg
 *  @version  $Id$
 *  @author   Nathan A. Baker
 *
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

#ifndef _VPMG_H_
#define _VPMG_H_

#include "apbscfg.h"

#include "maloc/maloc.h"

#include "generic/vhal.h"
#include "generic/vacc.h"
#include "generic/vcap.h"
#include "generic/vpbe.h"
#include "generic/mgparm.h"
#include "generic/pbeparm.h"
#include "generic/vmatrix.h"
#include "pmgc/mgdrvd.h"
#include "pmgc/newdrvd.h"
#include "pmgc/mgsubd.h"
#include "pmgc/mikpckd.h"
#include "pmgc/matvecd.h"
#include "mg/vpmgp.h"
#include "mg/vgrid.h"

/** @def VPMGMAXPART The maximum number of partitions the mesh can be divided into
 *  @ingroup Vpmg
 */
#define VPMGMAXPART 2000

/**
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @brief   Contains public data members for Vpmg class/module
 *
 *  Many of the routines and macros are borrowed from the main.c driver
 *  (written by Mike Holst) provided with the PMG code.
 *
 */
struct sVpmg {

  Vmem *vmem;  /**< Memory management object for this class */
  Vpmgp *pmgp;  /**< Parameters */
  Vpbe *pbe;  /**< Information about the PBE system */

#ifdef BURY_FORTRAN
  Vpde *pde;            /**< @todo doc */
  Vmgdriver *mgdriver;  /**< @todo doc */
#endif

  double *epsx;  /**< X-shifted dielectric map */
  double *epsy;  /**< Y-shifted dielectric map */
  double *epsz;  /**< Y-shifted dielectric map */
  double *kappa;  /**< Ion accessibility map (0 <= kappa(x) <= 1) */
  double *pot;  /**< Potential map */
  double *charge;  /**< Charge map */

  int *iparm;  /**< Passing int parameters to FORTRAN */
  double *rparm;  /**< Passing real parameters to FORTRAN */
  int *iwork;  /**< Work array */
  double *rwork;  /**< Work array */
  double *a1cf;  /**< Operator coefficient values (a11) -- this array can be
                  * overwritten */
  double *a2cf;  /**< Operator coefficient values (a22) -- this array can be
                   overwritten */
  double *a3cf;  /**< Operator coefficient values (a33) -- this array can be
                   overwritten */
  double *ccf;  /**< Helmholtz term -- this array can be overwritten */
  double *fcf;  /**< Right-hand side -- this array can be overwritten */
  double *tcf;  /**< True solution */
  double *u;  /**< Solution */
  double *xf;  /**< Mesh point x coordinates */
  double *yf;  /**< Mesh point y coordinates */
  double *zf;  /**< Mesh point z coordinates */
  double *gxcf;  /**< Boundary conditions for x faces */
  double *gycf;  /**< Boundary conditions for y faces */
  double *gzcf;  /**< Boundary conditions for z faces */
  double *pvec;  /**< Partition mask array */
  double extDiEnergy;  /**< Stores contributions to the dielectric energy from
                        * regions outside the problem domain */
  double extQmEnergy;  /**< Stores contributions to the mobile ion energy from
                        * regions outside the problem domain */
  double extQfEnergy;  /**< Stores contributions to the fixed charge energy
                        * from regions outside the problem domain */
  double extNpEnergy;  /**< Stores contributions to the apolar energy from
                        * regions outside the problem domain */
  Vsurf_Meth surfMeth;  /**< Surface definition method */
  double splineWin;  /**< Spline window parm for surf defs */
  Vchrg_Meth chargeMeth;  /**< Charge discretization method */
  Vchrg_Src chargeSrc;  /**< Charge source */

  int filled;  /**< Indicates whether Vpmg_fillco has been called */

  int useDielXMap;  /**< Indicates whether Vpmg_fillco was called with an
                      external x-shifted dielectric map */
  Vgrid *dielXMap;  /**< External x-shifted dielectric map */
  int useDielYMap;  /**< Indicates whether Vpmg_fillco was called with an
                     * external y-shifted dielectric map */
  Vgrid *dielYMap;  /**< External y-shifted dielectric map */
  int useDielZMap;  /**< Indicates whether Vpmg_fillco was called with an
                     * external z-shifted dielectric map */
  Vgrid *dielZMap;  /**< External z-shifted dielectric map */
  int useKappaMap;  /**< Indicates whether Vpmg_fillco was called with an
                     * external kappa map */
  Vgrid *kappaMap;  /**< External kappa map */
  int usePotMap;    /**< Indicates whether Vpmg_fillco was called with an
                       * external potential map */
  Vgrid *potMap;    /**< External potential map */

  int useChargeMap;  /**< Indicates whether Vpmg_fillco was called with an
                      * external charge distribution map */
  Vgrid *chargeMap;  /**< External charge distribution map */
};

/**
 *  @ingroup Vpmg
 *  @brief   Declaration of the Vpmg class as the Vpmg structure
 */
typedef struct sVpmg Vpmg;

/* /////////////////////////////////////////////////////////////////////////
/// Inlineable methods
//////////////////////////////////////////////////////////////////////////// */
#if !defined(VINLINE_VPMG)

    /** @brief   Return the memory used by this structure (and its contents)
     *           in bytes
     *  @ingroup Vpmg
     *  @author  Nathan Baker
     *  @return  The memory used by this structure and its contents in bytes
     */
    VEXTERNC unsigned long int Vpmg_memChk(
            Vpmg *thee  /**< Object for memory check */
            );

#else /* if defined(VINLINE_VPMG) */

#   define Vpmg_memChk(thee) (Vmem_bytes((thee)->vmem))

#endif /* if !defined(VINLINE_VPMG) */

/* /////////////////////////////////////////////////////////////////////////
/// Non-inlineable methods
//////////////////////////////////////////////////////////////////////////// */
/** @brief   Constructor for the Vpmg class (allocates new memory)
 *  @author  Nathan Baker
 *  @ingroup Vpmg
 *  @returns Pointer to newly allocated Vpmg object
 */
VEXTERNC Vpmg* Vpmg_ctor(
        Vpmgp *parms,  /**< PMG parameter object */
        Vpbe *pbe,  /**< PBE-specific variables */
        int focusFlag,  /**< 1 for focusing, 0 otherwise */
        Vpmg *pmgOLD,  /**< Old Vpmg object to use for boundary conditions */
        MGparm *mgparm,  /**< MGparm parameter object for boundary conditions */
        PBEparm_calcEnergy energyFlag  /**< What types of energies to calculate */
        );

/**
 * @brief  FORTRAN stub constructor for the Vpmg class (uses
 *         previously-allocated memory)
 *  @author  Nathan Baker
 *  @ingroup Vpmg
 *  @returns 1 if successful, 0 otherwise
 */
VEXTERNC int Vpmg_ctor2(
        Vpmg *thee,  /**< Memory location for object */
        Vpmgp *parms,  /**< PMG parameter object */
        Vpbe *pbe,  /**< PBE-specific variables */
        int focusFlag,  /**< 1 for focusing, 0 otherwise */
        Vpmg *pmgOLD,  /**< Old Vpmg object to use for boundary conditions (can
                         be VNULL if focusFlag = 0) */
        MGparm *mgparm,  /**< MGparm parameter object for boundary
                          * conditions (can be VNULL if focusFlag = 0) */
        PBEparm_calcEnergy energyFlag  /**< What types of energies to
                                        * calculate (ignored if focusFlag
                                        * = 0) */
        );

/** @brief   Object destructor
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 */
VEXTERNC void Vpmg_dtor(
        Vpmg **thee  /**< Pointer to memory location of object to be
                      * destroyed */
        );

/** @brief   FORTRAN stub object destructor
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 */
VEXTERNC void Vpmg_dtor2(
        Vpmg *thee  /**< Pointer to object to be destroyed */
        );

/** @brief  Fill the coefficient arrays prior to solving the equation
 *  @ingroup  Vpmg
 *  @author  Nathan Baker
 *  @returns  1 if successful, 0 otherwise
 *  @bug useDielMap could only be passed once, not three times, to this
 *       function - why not just once? that's what the call in routines.c
 *       ends up doing - just passing useDielMap three times. - P. Ellis 11/3/11
 */
VEXTERNC int Vpmg_fillco(
        Vpmg *thee,  /**< Vpmg object */
        Vsurf_Meth surfMeth,  /**< Surface discretization method */
        double splineWin,  /**< Spline window (in A) for surfMeth =
                            * VSM_SPLINE */
        Vchrg_Meth chargeMeth,  /**< Charge discretization method */
        int useDielXMap,  /**< Boolean to use dielectric map argument */
        Vgrid *dielXMap,  /**< External dielectric map */
        int useDielYMap,  /**< Boolean to use dielectric map argument */
        Vgrid *dielYMap,  /**< External dielectric map */
        int useDielZMap,  /**< Boolean to use dielectric map argument */
        Vgrid *dielZMap,  /**< External dielectric map */
        int useKappaMap,  /**< Boolean to use kappa map argument */
        Vgrid *kappaMap,  /**< External kappa map */
        int usePotMap,  /**< Boolean to use potential map argument */
        Vgrid *potMap,  /**< External potential map */
        int useChargeMap,  /**< Boolean to use charge map argument */
        Vgrid *chargeMap  /**< External charge map */
        );

/** @brief   Solve the PBE using PMG
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @returns  1 if successful, 0 otherwise
 */
VEXTERNC int Vpmg_solve(
        Vpmg *thee  /**< Vpmg object */
        );

/** @brief   Solve Poisson's equation with a homogeneous Laplacian operator
 *           using the solvent dielectric constant.  This solution is
 *           performed by a sine wave decomposition.
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @returns  1 if successful, 0 otherwise
 *  @note    This function is really only for testing purposes as the
 *           PMG multigrid solver can solve the homogeneous system much more
 *           quickly.  Perhaps we should implement an FFT version at some
 *           point...
 */
VEXTERNC int Vpmg_solveLaplace(
        Vpmg *thee  /**< Vpmg object */
        );

/** @brief   Get the total electrostatic energy.
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    The value of this observable may be modified by setting
 *           restrictions on the subdomain over which it is calculated.  Such
 *           limits can be set via Vpmg_setPart and are generally useful for
 *           parallel runs.
 *  @returns The electrostatic energy in units of k_B T.
 */
VEXTERNC double Vpmg_energy(
        Vpmg *thee,  /**< Vpmg object */
        int extFlag  /**< If this was a focused calculation, include (1 -- for
                      * serial calculations) or ignore (0 -- for parallel
                      * calculations) energy contributions from outside the
                      * focusing domain */
        );

/** @brief   Get the "fixed charge" contribution to the electrostatic energy
 *
 *           Using the solution at the finest mesh level, get the
 *           electrostatic energy due to the interaction of the fixed charges
 *           with the potential: \f[ G = \sum_i q_i u(r_i) \f]
 *           and return the result in units of k_B T.  Clearly, no
 *           self-interaction terms are removed.  A factor a 1/2 has to be
 *           included to convert this to a real energy.
 *
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    The value of this observable may be modified by setting
 *           restrictions on the subdomain over which it is calculated.  Such
 *           limits can be set via Vpmg_setPart and are generally useful for
 *           parallel runs.
 *  @returns The fixed charge electrostatic energy in units of k_B T.
 */
VEXTERNC double Vpmg_qfEnergy(
        Vpmg *thee,  /**< Vpmg object */
        int extFlag  /**< If this was a focused calculation, include (1 -- for
                      * serial calculations) or ignore (0 -- for parallel
                      * calculations) energy contributions from outside the
                      * focusing domain */
        );

/** @brief   Get the per-atom "fixed charge" contribution to the electrostatic
 *           energy
 *
 *           Using the solution at the finest mesh level, get the
 *           electrostatic energy due to the interaction of the fixed charges
 *           with the potential: \f[ G = q u(r), \f] where q$ is the
 *           charge and r is the location of the atom of interest.  The
 *           result is returned in units of k_B T.  Clearly, no
 *           self-interaction terms are removed.  A factor a 1/2 has to be
 *           included to convert this to a real energy.
 *
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    The value of this observable may be modified by setting
 *           restrictions on the subdomain over which it is calculated.  Such
 *           limits can be set via Vpmg_setPart and are generally useful for
 *           parallel runs.
 *  @returns The fixed charge electrostatic energy in units of k_B T.
 */
VEXTERNC double Vpmg_qfAtomEnergy(
        Vpmg *thee,  /**< The Vpmg object */
        Vatom *atom  /**< The atom for energy calculations */
        );

/** @brief Get the "mobile charge" contribution to the electrostatic energy.
 *
 *           Using the solution at the finest mesh level, get the
 *           electrostatic energy due to the interaction of the mobile charges
 *           with the potential:
 *              \f[ G = \frac{1}{4 I_s} \sum_i c_i q_i^2 \int
 *              \kappa^2(x) e^{-q_i u(x)} dx \f]
 *           for the NPBE and
 *              \f[ G = \frac{1}{2} \int \overline{\kappa}^2(x) u^2(x) dx \f]
 *           for the LPBE.  Here i denotes the counterion species,
 *           I_s is the bulk ionic strength, kappa^2(x)
 *           is the modified Debye-Huckel parameter, c_i is the
 *           concentration of species i, q_i is the charge of
 *           species i, and u(x) is the dimensionless electrostatic
 *           potential.  The energy is scaled to units of k_b T.
 *
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    The value of this observable may be modified by setting
 *           restrictions on the subdomain over which it is calculated.  Such
 *           limits can be set via Vpmg_setPart and are generally useful for
 *           parallel runs.
 *  @returns The mobile charge electrostatic energy in units of k_B T.
 */
VEXTERNC double Vpmg_qmEnergy(
        Vpmg *thee,  /**< Vpmg object */
        int extFlag  /**< If this was a focused calculation, include (1 -- for
                      * serial calculations) or ignore (0 -- for parallel
                      * calculations) energy contributions from outside the
                      * focusing domain */
        );


/** @brief Get the "polarization" contribution to the electrostatic energy.
 *
 *           Using the solution at the finest mesh level, get the
 *           electrostatic energy due to the interaction of the mobile charges
 *           with the potential:
 *              \f[ G = \frac{1}{2} \int \epsilon (\nabla u)^2 dx \f]
 *           where epsilon is the dielectric parameter and u(x) is
 *           the dimensionless electrostatic potential.  The energy is scaled
 *           to units of k_b T.
 *
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    The value of this observable may be modified by setting
 *           restrictions on the subdomain over which it is calculated.  Such
 *           limits can be set via Vpmg_setPart and are generally useful for
 *           parallel runs.
 *  @returns The polarization electrostatic energy in units of k_B T.
 */
VEXTERNC double Vpmg_dielEnergy(
        Vpmg *thee,  /**< Vpmg object */
        int extFlag  /**< If this was a focused calculation, include (1 -- for
                      * serial calculations) or ignore (0 -- for parallel
                      * calculations) energy contributions from outside the
                      * focusing domain */
        );


/** @brief Get the integral of the gradient of the dielectric function
 *
 *           Using the dielectric map at the finest mesh level, calculate the
 *           integral of the norm of the dielectric function gradient
 *           routines of Im et al (see Vpmg_dbForce for reference):
 *              \f[ \int \| \nabla \epsilon \| dx \f]
 *           where epsilon is the dielectric parameter.
 *           The integral is returned in units of A^2.
 *
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *           restrictions on the subdomain over which it is calculated.  Such
 *           limits can be set via Vpmg_setPart and are generally useful for
 *           parallel runs.
 *  @returns The integral in units of A^2.
 */
VEXTERNC double Vpmg_dielGradNorm(
        Vpmg *thee  /**< Vpmg object */
        );

/** @brief    Calculate the total force on the specified atom in units of
 *            k_B T/AA
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    \li Using the force evaluation methods of Im et al (Roux group),
 *             Comput Phys Commun, 111, 59--75 (1998).  However, this gives the
 *             whole (self-interactions included) force -- reaction field
 *             forces will have to be calculated at higher level.
 *           \li No contributions are made from higher levels of focusing.
 * @returns  1 if successful, 0 otherwise
 */
VEXTERNC int Vpmg_force(
        Vpmg *thee,  /**< Vpmg object */
        double *force, /**< 3*sizeof(double) space to hold the force in units
                         of k_B T/AA */
        int atomID,  /**< Valist ID of desired atom */
        Vsurf_Meth srfm,  /**< Surface discretization method */
        Vchrg_Meth chgm  /**< Charge discretization method */
        );

/** @brief    Calculate the "charge-field" force on the specified atom in units
 *           of k_B T/AA
 * @ingroup  Vpmg
 * @author   Nathan Baker
 * @note     \li Using the force evaluation methods of Im et al (Roux group),
 *             Comput Phys Commun, 111, 59--75 (1998).  However, this gives the
 *             whole (self-interactions included) force -- reaction field
 *             forces will have to be calculated at higher level.
 *           \li No contributions are made from higher levels of focusing.
 * @returns  1 if sucessful, 0 otherwise
 */
VEXTERNC int Vpmg_qfForce(
        Vpmg *thee,  /**< Vpmg object */
        double *force, /**< 3*sizeof(double) space to hold the force in units
                         of k_B T/A */
        int atomID,  /**< Valist ID of desired atom */
        Vchrg_Meth chgm  /**< Charge discretization method */
        );

/** @brief   Calculate the dielectric boundary forces on the
 *           specified atom in units of k_B T/AA
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    \li Using the force evaluation methods of Im et al (Roux group),
 *             Comput Phys Commun, 111, 59--75 (1998).  However, this gives the
 *             whole (self-interactions included) force -- reaction field
 *             forces will have to be calculated at higher level.
 *           \li No contributions are made from higher levels of focusing.
 * @returns  1 if successful, 0 otherwise
 */
VEXTERNC int Vpmg_dbForce(
        Vpmg *thee,  /**< Vpmg object */
        double *dbForce, /**< 3*sizeof(double) space to hold the dielectric
                           boundary force in units of k_B T/AA */
        int atomID,  /**< Valist ID of desired atom */
        Vsurf_Meth srfm  /**< Surface discretization method */
        );

/** @brief   Calculate the osmotic pressure on the specified atom in units of
 *           k_B T/AA
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @note    \li Using the force evaluation methods of Im et al (Roux group),
 *             Comput Phys Commun, 111, 59--75 (1998).  However, this gives the
 *             whole (self-interactions included) force -- reaction field
 *             forces will have to be calculated at higher level.
 *           \li No contributions are made from higher levels of focusing.
 *  @returns  1 if successful, 0 otherwise
 */
VEXTERNC int Vpmg_ibForce(
        Vpmg *thee,  /**< Vpmg object */
        double *force, /**< 3*sizeof(double) space to hold the
                           boundary force in units of k_B T/AA */
        int atomID,  /**< Valist ID of desired atom */
        Vsurf_Meth srfm  /**< Surface discretization method */
        );

/** @brief   Set partition information which restricts the calculation of
 *           observables to a (rectangular) subset of the problem domain
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 */
VEXTERNC void Vpmg_setPart(
        Vpmg *thee,  /**< Vpmg object */
        double lowerCorner[3],  /**< Partition lower corner */
        double upperCorner[3],  /**< Partition upper corner */
        int bflags[6]  /**< Booleans indicating whether a particular processor
                         is on the boundary with another partition.  0 if the
                         face is not bounded (next to) another partition, and
                         1 otherwise. */
        );

/** @brief  Remove partition restrictions
 *  @ingroup  Vpmg
 *  @author  Nathan Baker
 */
VEXTERNC void Vpmg_unsetPart(
        Vpmg *thee  /**< Vpmg object */
        );

/** @brief  Fill the specified array with accessibility values
 *  @ingroup  Vpmg
 *  @author  Nathan Baker
 *  @returns  1 if successful, 0 otherwise
 */
VEXTERNC int Vpmg_fillArray(
        Vpmg *thee,  /**< Vpmg object */
        double *vec,  /**< A nx*ny*nz*sizeof(double) array to contain the
                        values to be written */
        Vdata_Type type,  /**< What to write */
        double parm,  /**< Parameter for data type definition (if needed) */
        Vhal_PBEType pbetype, /**< Parameter for PBE type (if needed) */
        PBEparm * pbeparm /**< Pass in the PBE parameters (if needed) */
        );

/** @brief   Computes the field at an atomic center using a stencil based
 *           on the first derivative of a 5th order B-spline
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */
VPUBLIC void Vpmg_fieldSpline4(
             Vpmg *thee,     /**< Vpmg object */
             int atomID,     /**< Atom index */
             double field[3] /**< The (returned) electric field */
             );

/** @brief   Computes the permanent multipole electrostatic hydration
 *           energy (the polarization component of the hydration energy
 *           currently computed in TINKER).
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 *  @returns The permanent multipole electrostatic hydration energy
 */
VEXTERNC double Vpmg_qfPermanentMultipoleEnergy(
             Vpmg *thee,     /**< Vpmg object */
             int atomID      /**< Atom index */
             );

/** @brief   Computes the q-Phi Force for permanent multipoles based on
 *           5th order B-splines
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */
VEXTERNC void Vpmg_qfPermanentMultipoleForce(
             Vpmg *thee,      /**< Vpmg object */
             int atomID,      /**< Atom index */
             double force[3], /**< (returned) force */
             double torque[3] /**< (returned) torque */
             );

/** @brief   Compute the ionic boundary force for permanent multipoles.
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */
VEXTERNC void Vpmg_ibPermanentMultipoleForce(
             Vpmg *thee,      /**< Vpmg object */
             int atomID,      /**< Atom index */
             double force[3]  /**< (returned) force */
             );

/** @brief   Compute the dielectric boundary force for permanent multipoles.
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */
VEXTERNC void Vpmg_dbPermanentMultipoleForce(
             Vpmg *thee,      /**< Vpmg object */
             int atomID,      /**< Atom index */
             double force[3]  /**< (returned) force */
             );

/** @brief   q-Phi direct polarization force between permanent multipoles and
 *           induced dipoles, which are induced by the sum of the permanent
 *           intramolecular field and the permanent reaction field.
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */
VEXTERNC void Vpmg_qfDirectPolForce(
             Vpmg *thee,      /**< Vpmg object */
             Vgrid *perm,     /**< Permanent multipole potential */
             Vgrid *induced,  /**< Induced dipole potential */
             int atomID,      /**< Atom index */
             double force[3], /**< (returned) force */
             double torque[3] /**< (returned) torque */
             );

/** @brief   q-Phi direct polarization force between permanent multipoles and
 *           non-local induced dipoles based on 5th Order B-Splines.
 *           Keep in mind that the "non-local" induced dipooles are just
 *           a mathematical quantity that result from differentiation of
 *           the AMOEBA polarization energy.
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */
VEXTERNC void Vpmg_qfNLDirectPolForce(
             Vpmg *thee,      /**< Vpmg object */
             Vgrid *perm,     /**< Permanent multipole potential */
             Vgrid *nlInduced,/**< Non-local induced dipole potential */
             int atomID,      /**< Atom index */
             double force[3], /**< (returned) force */
             double torque[3] /**< (returned) torque */
             );

/** @brief   Ionic boundary direct polarization force between permanent
 *           multipoles and induced dipoles, which are induced by the
 *           sum of the permanent intramolecular field and the permanent
 *           reaction field.
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */
VEXTERNC void Vpmg_ibDirectPolForce(
             Vpmg *thee,      /**< Vpmg object */
             Vgrid *perm,     /**< Permanent multipole potential */
             Vgrid *induced,  /**< Induced dipole potential */
             int atomID,      /**< Atom index */
             double force[3]  /**< (returned) force */
             );

/** @brief   Ionic boundary direct polarization force between permanent
 *           multipoles and non-local induced dipoles based on 5th order
 *           Keep in mind that the "non-local" induced dipooles are just
 *           a mathematical quantity that result from differentiation of
 *           the AMOEBA polarization energy.
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */
VEXTERNC void Vpmg_ibNLDirectPolForce(
             Vpmg *thee,      /**< Vpmg object */
             Vgrid *perm,     /**< Permanent multipole potential */
             Vgrid *nlInduced,/**< Induced dipole potential */
             int atomID,      /**< Atom index */
             double force[3]  /**< (returned) force */
             );

/** @brief   Dielectric boundary direct polarization force between permanent
 *           multipoles and induced dipoles, which are induced by the
 *           sum of the permanent intramolecular field and the permanent
 *           reaction field.
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */
VEXTERNC void Vpmg_dbDirectPolForce(
             Vpmg *thee,      /**< Vpmg object */
             Vgrid *perm,     /**< Permanent multipole potential */
             Vgrid *induced,  /**< Induced dipole potential */
             int atomID,      /**< Atom index */
             double force[3]  /**< (returned) force */
             );

/** @brief   Dielectric bounday direct polarization force between
 *           permanent multipoles and non-local induced dipoles.
 *           Keep in mind that the "non-local" induced dipooles are just
 *           a mathematical quantity that result from differentiation of
 *           the AMOEBA polarization energy.
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */
VEXTERNC void Vpmg_dbNLDirectPolForce(
             Vpmg *thee,      /**< Vpmg object */
             Vgrid *perm,     /**< Permanent multipole potential */
             Vgrid *nlInduced,/**< Non-local induced dipole potential */
             int atomID,      /**< Atom index */
             double force[3]  /**< (returned) force */
             );

/** @brief   Mutual polarization force for induced dipoles based on 5th
 *           order B-Splines. This force arises due to self-consistent
 *           convergence of the solute induced dipoles and reaction field.
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */
VEXTERNC void Vpmg_qfMutualPolForce(
             Vpmg *thee,      /**< Vpmg object */
             Vgrid *induced,  /**< Induced dipole potential */
             Vgrid *nlInduced,/**< Non-local induced dipole potential */
             int atomID,      /**< Atom index */
             double force[3]  /**< (returned) force */
             );

/** @brief   Ionic boundary mutual polarization force for induced dipoles
 *           based on 5th order B-Splines. This force arises due to
 *           self-consistent convergence of the solute induced dipoles
 *           and reaction field.
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */
VEXTERNC void Vpmg_ibMutualPolForce(
             Vpmg *thee,      /**< Vpmg object */
             Vgrid *induced,  /**< Induced dipole potential */
             Vgrid *nlInduced,/**< Non-local induced dipole potential */
             int atomID,      /**< Atom index */
             double force[3]  /**< (returned) force */
             );

/** @brief   Dielectric boundary mutual polarization force for induced dipoles
 *           based on 5th order B-Splines. This force arises due to
 *           self-consistent convergence of the solute induced dipoles
 *           and reaction field.
 *  @ingroup Vpmg
 *  @author  Michael Schnieders
 */
VEXTERNC void Vpmg_dbMutualPolForce(
             Vpmg *thee,      /**< Vpmg object */
             Vgrid *induced,  /**< Induced dipole potential */
             Vgrid *nlInduced,/**< Non-local induced dipole potential */
             int atomID,      /**< Atom index */
             double force[3]  /**< (returned) force */
             );

/** @brief   Print out a column-compressed sparse matrix in Harwell-Boeing
 *           format.
 *  @ingroup Vpmg
 *  @author  Nathan Baker
 *  @bug  Can this path variable be replaced with a Vio socket?
 */
VEXTERNC void Vpmg_printColComp(
        Vpmg *thee,  /**<  Vpmg object */
        char path[72],  /**< The file to which the matrix is to be written */
        char title[72],  /**< The title of the matrix */
        char mxtype[3],   /**< The type of REAL-valued matrix, a 3-character
                            string of the form "R_A" where the '_' can be one
                            of:
                            \li S:  symmetric matrix
                            \li U:  unsymmetric matrix
                            \li H:  Hermitian matrix
                            \li Z:  skew-symmetric matrix
                            \li R:  rectangular matrix */
        int flag  /**< The operator to compress:
                    \li 0:  Poisson operator
                    \li 1:  Linearization of the full Poisson-Boltzmann
                            operator around the current solution */
        );



/** @brief   Build a column-compressed matrix in Harwell-Boeing format
 *  @ingroup Vpmg
 *  @author  Tucker Beck [C Translation]
 *           Nathan Baker [Original] (mostly ripped off from Harwell-Boeing
 *                                    format documentation)Michael Schnieders)
 */
VPRIVATE void bcolcomp(
        int    *iparm,  ///< @todo Document
        double *rparm,  ///< @todo Document
        int    *iwork,  ///< @todo Document
        double *rwork,  ///< @todo Document
        double *values, ///< @todo Document
        int    *rowind, ///< @todo Document
        int    *colptr, ///< @todo Document
        int    *flag    /**< Operation selection parameter
                         *     0 = Use Poisson operator only
                         *     1 = Use linearization of full operation around
                         *         current solution.
                         */
        );



/** @brief   Build a column-compressed matrix in Harwell-Boeing format
 *  @ingroup Vpmg
 *  @author  Tucker Beck [C Translation]
 *           Nathan Baker [Original] (mostly ripped off from Harwell-Boeing
 *                                    format documentation)Michael Schnieders)
 */
VPRIVATE void bcolcomp2(
        int    *iparm,  ///< @todo Document
        double *rparm,  ///< @todo Document
        int    *nx,     ///< @todo Document
        int    *ny,     ///< @todo Document
        int    *nz,     ///< @todo Document
        int    *iz,     ///< @todo Document
        int    *ipc,    ///< @todo Document
        double *rpc,    ///< @todo Document
        double *ac,     ///< @todo Document
        double *cc,     ///< @todo Document
        double *values, ///< @todo Document
        int    *rowind, ///< @todo Document
        int    *colptr, ///< @todo Document
        int    *flag    /**< Operation selection parameter
                         *     0 = Use Poisson operator only
                         *     1 = Use linearization of full operation around
                         *         current solution.
                         */
        );



/** @brief   Build a column-compressed matrix in Harwell-Boeing format
 *  @ingroup Vpmg
 *  @author  Tucker Beck [C Translation]
 *           Nathan Baker [Original] (mostly ripped off from Harwell-Boeing
 *                                    format documentation)Michael Schnieders)
 */
VPRIVATE void bcolcomp3(
        int    *nx,     ///< @todo Document
        int    *ny,     ///< @todo Document
        int    *nz,     ///< @todo Document
        int    *ipc,    ///< @todo Document
        double *rpc,    ///< @todo Document
        double *ac,     ///< @todo Document
        double *cc,     ///< @todo Document
        double *values, ///< @todo Document
        int    *rowind, ///< @todo Document
        int    *colptr, ///< @todo Document
        int    *flag    ///< @todo Document
        );



/** @brief   Build a column-compressed matrix in Harwell-Boeing format
 *  @ingroup Vpmg
 *  @author  Tucker Beck [C Translation]
 *           Nathan Baker [Original] (mostly ripped off from Harwell-Boeing
 *                                    format documentation)Michael Schnieders)
 */
VPRIVATE void bcolcomp4(
        int    *nx,     ///< @todo Document
        int    *ny,     ///< @todo Document
        int    *nz,     ///< @todo Document
        int    *ipc,    ///< @todo Document
        double *rpc,    ///< @todo Document
        double *oC,     ///< @todo Document
        double *cc,     ///< @todo Document
        double *oE,     ///< @todo Document
        double *oN,     ///< @todo Document
        double *uC,     ///< @todo Document
        double *values, ///< @todo Document
        int    *rowind, ///< @todo Document
        int    *colptr, ///< @todo Document
        int    *flag    ///< @todo Document
        );



/** @brief   Print a column-compressed matrix in Harwell-Boeing format
 *  @ingroup Vpmg
 *  @author  Tucker Beck [C Translation]
 *           Nathan Baker [Original] (mostly ripped off from Harwell-Boeing
 *                                    format documentation)Michael Schnieders)
 */
VPRIVATE void pcolcomp(
        int    *nrow,   ///< @todo Document
        int    *ncol,   ///< @todo Document
        int    *nnzero, ///< @todo Document
        double *values, ///< @todo Document
        int    *rowind, ///< @todo Document
        int    *colptr, ///< @todo Document
        char   *path,   ///< @todo Document
        char   *title,  ///< @todo Document
        char   *mxtype  ///< @todo Document
        );



/* ///////////////////////////////////////////////////////////////////////////
// Internal routines
/////////////////////////////////////////////////////////////////////////// */

/**
 * @brief  Evaluate a cubic B-spline
 * @author  Nathan Baker
 * @return  Cubic B-spline value
 */
VPRIVATE double bspline2(
        double x  /** Position */
        );

/**
 * @brief  Evaluate a cubic B-spline derivative
 * @author  Nathan Baker
 * @return  Cubic B-spline derivative
 */
VPRIVATE double dbspline2(
        double x  /** Position */
        );

/**
 * @brief   Return 2.5 plus difference of i - f
 * @author  Michael Schnieders
 * @return  (2.5+((double)(i)-(f)))
 */
VPRIVATE double VFCHI4(
        int i,
        double f
        );

/**
 * @brief   Evaluate a 5th Order B-Spline (4th order polynomial)
 * @author: Michael Schnieders
 * @return  5th Order B-Spline
 */
VPRIVATE double bspline4(
         double x /** Position */
         );

/**
 * @brief   Evaluate a 5th Order B-Spline derivative (4th order polynomial)
 * @author: Michael Schnieders
 * @return  5th Order B-Spline derivative
 */
VPRIVATE double dbspline4(
         double x /** Position */
         );

/**
 * @brief   Evaluate the 2nd derivative of a 5th Order B-Spline
 * @author: Michael Schnieders
 * @return  2nd derivative of a 5th Order B-Spline
 */
VPRIVATE double d2bspline4(
         double x /** Position */
         );

/**
 * @brief   Evaluate the 3rd derivative of a 5th Order B-Spline
 * @author: Michael Schnieders
 * @return  3rd derivative of a 5th Order B-Spline
 */
VPRIVATE double d3bspline4(
         double x /** Position */
         );

/**
 * @brief  Determines energy from polarizeable charge and interaction with
 *         fixed charges according to Rocchia et al.
 * @author  Nathan Baker
 * @return  Energy in kT
 */
VPRIVATE double Vpmg_polarizEnergy(
         Vpmg *thee,
         int extFlag  /** If 1, add external energy contributions to
                       result */
         );
/**
 * @brief  Calculates charge-potential energy using summation over delta
 *         function positions (i.e. something like an Linf norm)
 * @author  Nathan Baker
 * @return  Energy in kT
 */
VPRIVATE double Vpmg_qfEnergyPoint(
        Vpmg *thee,
        int extFlag  /** If 1, add external energy contributions to
                       result */
        );

/**
 * @brief  Calculates charge-potential energy as integral over a volume
 * @author  Nathan Baker
 * @return  Energy in kT
 */
VPRIVATE double Vpmg_qfEnergyVolume(
        Vpmg *thee,
        int extFlag  /** If 1, add external energy contributions to
                       result */
        );

/**
* @brief Selects a spline based surface method from either VSM_SPLINE,
 *        VSM_SPLINE5 or VSM_SPLINE7
 * @author David Gohara
 */
VPRIVATE void Vpmg_splineSelect(
        int srfm,		/** Surface method, currently VSM_SPLINE,
        VSM_SPLINE5, or VSM_SPLINE7 */
        Vacc *acc,		/** Accessibility object */
        double *gpos,	/** Position array -> array[3] */
        double win,		/** Spline window */
        double infrad,	/** Inflation radius */
        Vatom *atom,	/** Atom object */
        double *force	/** Force array -> array[3] */
        );

/**
 * @brief  For focusing, fill in the boundaries of the new mesh based on the
 * potential values in the old mesh
 * @author  Nathan Baker
 */
VPRIVATE void focusFillBound(
        Vpmg *thee,  /** New PMG object (the one just created) */
        Vpmg *pmg  /** Old PMG object */
        );

/**
 * @brief  Increment all boundary points by
 *         pre1*(charge/d)*(exp(-xkappa*(d-size))/(1+xkappa*size) to add the
 *         effect of the Debye-Huckel potential due to a single charge
 * @author  Nathan Baker
 */
VPRIVATE void bcfl1(
        double size,  /** Size of the ion */
        double *apos,  /** Position of the ion */
        double charge,  /** Charge of the ion */
        double xkappa,  /** Exponential screening factor */
        double pre1,  /** Unit- and dielectric-dependent prefactor */
        double *gxcf,  /** Set to x-boundary values */
        double *gycf,  /** Set to y-boundary values */
        double *gzcf,  /** Set to z-boundary values */
        double *xf,  /** Boundary point x-coordinates */
        double *yf,  /** Boundary point y-coordinates */
        double *zf,  /** Boundary point z-coordinates */
        int nx,  /** Number of grid points in x-direction */
        int ny,  /** Number of grid points in y-direction */
        int nz /** Number of grid points in y-direction */
        );

/**
 * @brief  Increment all boundary points to include the Debye-Huckel
 *         potential due to a single multipole site. (truncated at quadrupole)
 * @author Michael Schnieders
 */
VPRIVATE void bcfl2(
        double size,  /** Size of the ion */
        double *apos,  /** Position of the ion */
        double charge,  /** Charge of the ion */
        double *dipole, /** Dipole of the ion */
        double *quad,   /** Traceless Quadrupole of the ion */
        double xkappa,  /** Exponential screening factor */
        double eps_p,   /** Solute dielectric */
        double eps_w,   /** Solvent dielectric */
        double T,       /** Temperature */
        double *gxcf,  /** Set to x-boundary values */
        double *gycf,  /** Set to y-boundary values */
        double *gzcf,  /** Set to z-boundary values */
        double *xf,  /** Boundary point x-coordinates */
        double *yf,  /** Boundary point y-coordinates */
        double *zf,  /** Boundary point z-coordinates */
        int nx,  /** Number of grid points in x-direction */
        int ny,  /** Number of grid points in y-direction */
        int nz /** Number of grid points in y-direction */
        );

/**
 * @brief  This routine serves bcfl2. It returns (in tsr) the contraction
 *         independent portion of the Debye-Huckel potential tensor
 *         for a spherical ion with a central charge, dipole and quadrupole.
 *         See the code for an in depth description.
 *
 * @author Michael Schnieders
 */
VPRIVATE void multipolebc(
        double r,      /** Distance to the boundary */
        double kappa,  /** Exponential screening factor */
        double eps_p,  /** Solute dielectric */
        double eps_w,  /** Solvent dielectric */
        double rad,    /** Radius of the sphere */
        double tsr[3]  /** Contraction-independent portion of each tensor */
        );

/**
 * @brief  Calculate
 *         pre1*(charge/d)*(exp(-xkappa*(d-size))/(1+xkappa*size) due to a
 *         specific ion at a specific point
 * @author  Nathan Baker
 * @returns  Value of above function in arbitrary units (dependent on
 *           pre-factor)
 */
VPRIVATE double bcfl1sp(
        double size,  /** Atom size */
        double *apos,  /** Atom position */
        double charge,  /** Atom charge */
        double xkappa,  /** Exponential screening factor */
        double pre1,  /** Unit- and dielectric-dependent prefactor */
        double *pos  /** Function evaluation position */
        );

/**
 * @brief  Fill boundary condition arrays
 * @author  Nathan Baker
 */
VPRIVATE void bcCalc(
        Vpmg *thee
        );

/**
 * @brief  Top-level driver to fill all operator coefficient arrays
 * @author  Nathan Baker
 */
VPRIVATE void fillcoCoef(
        Vpmg *thee
        );

/**
 * @brief  Fill operator coefficient arrays from pre-calculated maps
 * @author  Nathan Baker
 */
VPRIVATE void fillcoCoefMap(
        Vpmg *thee
        );

/**
 * @brief  Fill operator coefficient arrays from a molecular surface
 *         calculation
 * @author  Nathan Baker
 */
VPRIVATE void fillcoCoefMol(
        Vpmg *thee
        );

/**
 * @brief  Fill ion (nonlinear) operator coefficient array from a molecular
 * surface calculation
 * @author  Nathan Baker
 */
VPRIVATE void fillcoCoefMolIon(
        Vpmg *thee
        );

/**
 * @brief  Fill differential operator coefficient arrays from a molecular
 *         surface calculation
 * @author  Nathan Baker
 */
VPRIVATE void fillcoCoefMolDiel(
        Vpmg *thee
        );

/**
 * @brief  Fill differential operator coefficient arrays from a molecular
 *         surface calculation without smoothing
 * @author  Nathan Baker
 */
VPRIVATE void fillcoCoefMolDielNoSmooth(
        Vpmg *thee
        );

/**
 * @brief  Fill differential operator coefficient arrays from a molecular
 *         surface calculation with smoothing.
 *
 *         Molecular surface, dielectric smoothing following an implementation
 *         of Bruccoleri, et al.  J Comput Chem 18 268-276 (1997).
 *
 *         This algorithm uses a 9 point harmonic smoothing technique - the point
 *         in question and all grid points 1/sqrt(2) grid spacings away.
 *
 * @note   This uses thee->a1cf, thee->a2cf, thee->a3cf as temporary storage.
 * @author  Todd Dolinsky
 */
VPRIVATE void fillcoCoefMolDielSmooth(
        Vpmg *thee
        );

/**
 * @brief  Fill operator coefficient arrays from a spline-based surface
 *         calculation
 * @author  Nathan Baker
 */
VPRIVATE void fillcoCoefSpline(
        Vpmg *thee
        );

/**
* @brief  Fill operator coefficient arrays from a 5th order polynomial
*         based surface calculation
* @author  Michael Schnieders
*/
VPRIVATE void fillcoCoefSpline3(
        Vpmg *thee
        );

/**
 * @brief  Fill operator coefficient arrays from a 7th order polynomial
 *         based surface calculation
 * @author  Michael Schnieders
 */
VPRIVATE void fillcoCoefSpline4(
        Vpmg *thee
        );

/**
 * @brief  Top-level driver to fill source term charge array
 * @returns  Success/failure status
 * @author  Nathan Baker
 */
VPRIVATE Vrc_Codes fillcoCharge(
        Vpmg *thee
        );

/**
 * @brief  Fill source term charge array from a pre-calculated map
 * @returns  Success/failure status
 * @author  Nathan Baker
 */
VPRIVATE Vrc_Codes fillcoChargeMap(
        Vpmg *thee
        );

/**
 * @brief  Fill source term charge array from linear interpolation
 * @author  Nathan Baker
 */
VPRIVATE void fillcoChargeSpline1(
        Vpmg *thee
        );

/**
 * @brief  Fill source term charge array from cubic spline interpolation
 * @author  Nathan Baker
 */
VPRIVATE void fillcoChargeSpline2(
        Vpmg *thee
        );

/**
 * @brief  Fill source term charge array for the use of permanent multipoles
 * @author  Michael Schnieders
 */
VPRIVATE void fillcoPermanentMultipole(
        Vpmg *thee
        );

/**
 * @brief  Fill source term charge array for use of induced dipoles
 * @author  Michael Schnieders
 */
VPRIVATE void fillcoInducedDipole(
        Vpmg *thee
        );

/**
 * @brief  Fill source term charge array for non-local induced
 * dipoles
 * @author  Michael Schnieders
 */
VPRIVATE void fillcoNLInducedDipole(
        Vpmg *thee
        );

/**
 * @brief  For focusing, set external energy data members in new Vpmg object
 *         based on energy calculations on old Vpmg object from regions
 *         outside the indicated partition.
 * @author  Nathan Baker, Todd Dolinsky
 */
VPRIVATE void extEnergy(
        Vpmg *thee,  /** Newly created PMG manager */
        Vpmg *pmgOLD,  /** Old PMG manager, source of energies */
        PBEparm_calcEnergy extFlag,  /** Energy calculation flag */
        double partMin[3],  /** Partition lower corner */
        double partMax[3],  /** Partition upper corner */
        int bflags[6]  /** Which boundaries to include the calculation */
        );

/**
 * @brief  Charge-field force due to a linear spline charge function
 * @author  Nathan Baker
 */
VPRIVATE void qfForceSpline1(
        Vpmg *thee,
        double *force,  /** Set to force */
        int atomID  /** Valist atom ID */
        );

/**
 * @brief  Charge-field force due to a cubic spline charge function
 * @author  Nathan Baker
 */
VPRIVATE void qfForceSpline2(
        Vpmg *thee,
        double *force,  /** Set to force */
        int atomID  /** Valist atom ID */
        );

/**
* @brief  Charge-field force due to a quintic spline charge function
* @author  Michael Schnieders
*/
VPRIVATE void qfForceSpline4(
                             Vpmg *thee,
                             double *force,  /** Set to force */
                             int atomID  /** Valist atom ID */
                             );


/**
 * @brief  Calculate the solution to Poisson's equation with a simple
 *         Laplacian operator and zero-valued Dirichlet boundary conditions.
 *         Store the solution in thee->u.
 * @author  Nathan Baker
 * @note  Vpmg_fillco must be called first
 */
VPRIVATE void zlapSolve(
        Vpmg *thee,
        double **solution,  /** Solution term vector */
        double **source,  /** Source term vector */
        double **work1  /** Work vector */
        );

/**
 * @brief  Mark the grid points inside a sphere with a particular value.  This
 *         marks by resetting the the grid points inside the sphere to the
 *         specified value.
 * @author  Nathan Baker
 */
VPRIVATE void markSphere(
        double rtot,  /** Sphere radius */
        double *tpos,  /** Sphere position */
        int nx,  /** Number of grid points */
        int ny,  /** Number of grid points */
        int nz,  /** Number of grid points */
        double hx,  /** Grid spacing */
        double hy,  /** Grid spacing */
        double hzed,  /** Grid spacing */
        double xmin,  /** Grid lower corner */
        double ymin,  /** Grid lower corner */
        double zmin,  /** Grid lower corner */
        double *array,  /** Grid values */
        double markVal  /** Value to mark with */
        );

/**
 * @brief Vpmg_qmEnergy for SMPBE
 * @author Vincent Chu
 */
VPRIVATE double Vpmg_qmEnergySMPBE(Vpmg *thee, int extFlag);
VPRIVATE double Vpmg_qmEnergyNONLIN(Vpmg *thee, int extFlag);



// Additional macros and definitions.  May not be needed

// Added by Vincent Chu 9/13/06 for SMPB
#define VCUB(x)            ((x)*(x)*(x))
#define VLOG(x)            (log(x))

#define IJK(i,j,k)  (((k)*(nx)*(ny))+((j)*(nx))+(i))
#define IJKx(j,k,i) (((i)*(ny)*(nz))+((k)*(ny))+(j))
#define IJKy(i,k,j) (((j)*(nx)*(nz))+((k)*(nx))+(i))
#define IJKz(i,j,k) (((k)*(nx)*(ny))+((j)*(nx))+(i))
#define VFCHI(iint,iflt) (1.5+((double)(iint)-(iflt)))


#endif    /* ifndef _VPMG_H_ */

