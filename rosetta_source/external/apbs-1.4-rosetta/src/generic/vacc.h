/** @defgroup Vacc Vacc class
 *  @brief    Solvent- and ion-accessibility oracle
 */

/**
 *  @file     vacc.h
 *  @ingroup  Vacc
 *  @brief    Contains declarations for class Vacc
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

#ifndef _VACC_H_
#define _VACC_H_

#include "apbscfg.h"

#include "maloc/maloc.h"
#if defined(HAVE_MC_H)
#include "mc/mc.h"
#endif

#include "generic/vhal.h"
#include "generic/valist.h"
#include "generic/vclist.h"
#include "generic/vatom.h"
#include "generic/vunit.h"
#include "generic/apolparm.h"

/**
 * @ingroup  Vacc
 * @author  Nathan Baker
 * @brief  Surface object list of per-atom surface points
 */
struct sVaccSurf {
    Vmem *mem;  /**< Memory object */
    double *xpts;  /**< Array of point x-locations */
    double *ypts;  /**< Array of point y-locations */
    double *zpts;  /**< Array of point z-locations */
    char *bpts;  /**< Array of booleans indicating whether a point is (1) or is
                 * not (0) part of the surface */
    double area;  /**< Area spanned by these points */
    int npts;  /**< Length of thee->xpts, ypts, zpts arrays */
    double probe_radius;  /**< Probe radius (A) with which this surface was
                     * constructed */
};

/**
 *  @ingroup Vacc
 *  @brief   Declaration of the VaccSurf class as the VaccSurf structure
 */
typedef struct sVaccSurf VaccSurf;

/**
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @brief   Oracle for solvent- and ion-accessibility around a biomolecule
 */
struct sVacc {

  Vmem *mem;  /**< Memory management object for this class */
  Valist *alist;  /**< Valist structure for list of atoms */
  Vclist *clist;  /**< Vclist structure for atom cell list */
  int *atomFlags;  /**< Array of boolean flags of length
                    * Valist_getNumberAtoms(thee->alist) to prevent
                    * double-counting atoms during calculations */
  VaccSurf *refSphere;  /**< Reference sphere for SASA calculations */
  VaccSurf **surf;  /**< Array of surface points for each atom; is not
                    * initialized until needed (test against VNULL to
                    * determine initialization state) */
  Vset acc;  /**< An integer array (to be treated as bitfields) of Vset type
              * with length equal to the number of vertices in the mesh */
  double surf_density;  /**< Minimum solvent accessible surface point density
                         * (in pts/A^2) */

};

/**
 *  @ingroup Vacc
 *  @brief   Declaration of the Vacc class as the Vacc structure
 */
typedef struct sVacc Vacc;

#if !defined(VINLINE_VACC)

    /** @brief   Get number of bytes in this object and its members
     *  @ingroup Vacc
     *  @author  Nathan Baker
     *  @returns Number of bytes allocated for object
     */
    VEXTERNC unsigned long int Vacc_memChk(
            Vacc *thee /**< Object for memory check */
            );

#else /* if defined(VINLINE_VACC) */

#   define Vacc_memChk(thee) (Vmem_bytes((thee)->mem))

#endif /* if !defined(VINLINE_VACC) */

/**
 * @brief  Allocate and construct the surface object; do not assign surface
 *         points to positions
 * @ingroup  Vacc
 * @author  Nathan Baker
 * @returns  Newly allocated and constructed surface object
 */
VEXTERNC VaccSurf* VaccSurf_ctor(
        Vmem *mem,  /**< Memory manager (can be VNULL) */
        double probe_radius,  /**< Probe radius (in A) for this surface */
        int nsphere  /**< Number of points in sphere */
        );

/**
 * @brief  Construct the surface object using previously allocated memory; do
 *         not assign surface points to positions
 * @ingroup  Vacc
 * @author  Nathan Baker
 * @returns  1 if successful, 0 otherwise
 */
VEXTERNC int VaccSurf_ctor2(
        VaccSurf *thee,  /**< Allocated memory */
        Vmem *mem,  /**< Memory manager (can be VNULL) */
        double probe_radius,  /**< Probe radius (in A) for this surface */
        int nsphere  /**< Number of points in sphere */
        );

/**
 * @brief  Destroy the surface object and free its memory
 * @ingroup  Vacc
 * @author  Nathan Baker
 */
VEXTERNC void VaccSurf_dtor(
        VaccSurf **thee  /**< Object to be destroyed */
        );

/**
 * @brief  Destroy the surface object
 * @ingroup  Vacc
 * @author  Nathan Baker
 */
VEXTERNC void VaccSurf_dtor2(
        VaccSurf *thee  /**< Object to be destroyed */
        );

/**
 * @brief Set up an array of points for a reference sphere of unit radius
 *
 * Generates approximately npts # of points (actual number stored in
 * thee->npts) somewhat uniformly distributed across a sphere of unit radius
 * centered at the origin.
 *
 * @note  This routine was shamelessly ripped off from sphere.f from UHBD as
 *        developed by Michael K. Gilson.
 *
 * @ingroup Vacc
 * @author  Nathan Baker (original FORTRAN code by Mike Gilson)
 * @return  Reference sphere surface object
 */
VEXTERNC VaccSurf* VaccSurf_refSphere(
        Vmem *mem,  /**< Memory object */
        int npts /**< Requested number of points on sphere */
        );

/**
 * @brief  Set up an array of points corresponding to the SAS due to a
 *         particular atom.
 * @ingroup  Vacc
 * @author  Nathan Baker
 * @return  Atom sphere surface object
 */
VEXTERNC VaccSurf* Vacc_atomSurf(
        Vacc *thee,  /**< Accessibility object for molecule */
        Vatom *atom,  /**< Atom for which the surface should be constructed */
        VaccSurf *ref,  /**< Reference sphere which sets the resolution for the
                         * surface. @see VaccSurf_refSphere */
        double probe_radius  /**< Probe radius (in A) */
        );

/** @brief   Construct the accessibility object
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @returns Newly allocated Vacc object */
VEXTERNC Vacc* Vacc_ctor(
        Valist *alist,  /**< Molecule for accessibility queries */
        Vclist *clist,  /**< Pre-constructed cell list for looking up atoms
                        * near specific positions */
        double surf_density  /**< Minimum per-atom solvent accessible surface
                             * point density (in pts/A^2)*/
        );

/** @brief   FORTRAN stub to construct the accessibility object
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @returns 1 if successful, 0 otherwise */
VEXTERNC int Vacc_ctor2(
        Vacc *thee, /**< Memory for Vacc objet */
        Valist *alist, /**< Molecule for accessibility queries */
        Vclist *clist, /**< Pre-constructed cell list for looking up atoms
                        * near specific positions */
        double surf_density  /**< Minimum per-atom solvent accessible surface
                             * point density (in pts/A^2)*/
        );

/** @brief   Destroy object
 *  @ingroup Vacc
 *  @author  Nathan Baker
 */
VEXTERNC void Vacc_dtor(
        Vacc **thee /**< Pointer to memory location of object */
        );

/** @brief   FORTRAN stub to destroy object
 *  @ingroup Vacc
 *  @author  Nathan Baker
 */
VEXTERNC void Vacc_dtor2(
        Vacc *thee /**< Pointer to object */
        );

/** @brief   Report van der Waals accessibility
 *
 *  Determines if a point is within the union of the atomic spheres (with
 *  radii equal to their van der Waals radii).
 *
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @returns Characteristic function value between 1.0 (accessible) and 0.0
 *          (inaccessible)
 */
VEXTERNC double Vacc_vdwAcc(
        Vacc *thee,  /**< Accessibility object */
        double center[VAPBS_DIM] /**< Probe center coordinates */
        );

/** @brief   Report inflated van der Waals accessibility
 *
 *  Determines if a point is within the union of the spheres centered at the
 *  atomic centers with radii equal to the sum of the atomic van der Waals
 *  radius and the probe radius.
 *
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @returns Characteristic function value between 1.0 (accessible) and 0.0
 *          (inaccessible)
 */
VEXTERNC double Vacc_ivdwAcc(
        Vacc *thee, /**< Accessibility object */
        double center[VAPBS_DIM], /**< Probe center coordinates */
        double radius /**< Probe radius (&Aring;) */
        );

/** @brief   Report molecular accessibility
 *
 * Determine accessibility of a probe (of radius radius) at a given point,
 * given a collection of atomic spheres.  Uses molecular (Connolly) surface
 * definition.
 *
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @returns Characteristic function value between 1.0 (accessible) and 0.0
 *          (inaccessible)
 *  @bug     This routine has a slight bug which can generate very small
 *           internal regions of high dielectric (thanks to John Mongan and
 *           Jess Swanson for finding this)
 */
VEXTERNC double Vacc_molAcc(
        Vacc *thee, /**< Accessibility object */
        double center[VAPBS_DIM], /**< Probe center coordinates */
        double radius /**< Probe radius (in &Aring;) */
        );

/** @brief   Report molecular accessibility quickly
 *
 *  Given a point which is INSIDE the collection of inflated van der Waals
 *  spheres, but OUTSIDE the collection of non-inflated van der Waals spheres,
 *  determine accessibility of a probe (of radius radius) at a given point,
 *  given a collection of atomic spheres.  Uses molecular (Connolly) surface
 *  definition.
 *
 *  @note    THIS ASSUMES YOU HAVE TESTED THAT THIS POINT IS DEFINITELY INSIDE
 *           THE INFLATED AND NON-INFLATED VAN DER WAALS SURFACES!
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @returns Characteristic function value between 1.0 (accessible) and 0.0
 *          (inaccessible)
 *  @bug     This routine has a slight bug which can generate very small
 *           internal regions of high dielectric (thanks to John Mongan and
 *           Jess Swanson for finding this)
 */
VEXTERNC double Vacc_fastMolAcc(
        Vacc *thee,  /**< Accessibility object */
        double center[VAPBS_DIM],  /**< Probe center coordinates */
        double radius /**< Probe radius (in &Aring;) */
        );

/** @brief   Report spline-based accessibility
 *
 *  Determine accessibility at a given point, given a collection of atomic
 *  spheres.  Uses Benoit Roux (Im et al, Comp Phys Comm, 111, 59--75, 1998)
 *  definition suitable for force evalation; basically a cubic spline.
 *
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @returns Characteristic function value between 1.0 (accessible) and 0.0
 *          (inaccessible)
 */
VEXTERNC double Vacc_splineAcc(
        Vacc *thee, /**< Accessibility object */
        double center[VAPBS_DIM], /**< Probe center coordinates */
        double win, /**< Spline window (&Aring;) */
        double infrad  /**< Inflation radius (&Aring;) for ion access. */
        );

/** @brief   Report gradient of spline-based accessibility.
 *
 *  @ingroup Vacc
 *  @author  Nathan Baker
 */
VEXTERNC void Vacc_splineAccGrad(
        Vacc *thee, /**< Accessibility object */
        double center[VAPBS_DIM], /**< Probe center coordinates */
        double win, /**< Spline window (&Aring;) */
        double infrad,  /**< Inflation radius (&Aring;) for ion access. */
        double *grad /**< 3-vector set to gradient of accessibility */
        );

/** @brief   Report spline-based accessibility for a given atom
 *
 *  Determine accessibility at a given point for a given atomic
 *  spheres.  Uses Benoit Roux (Im et al, Comp Phys Comm, 111, 59--75, 1998)
 *  definition suitable for force evalation; basically a cubic spline.
 *
 *  @ingroup Vacc
 *  @author  Nathan Baker
 *  @returns Characteristic function value between 1.0 (accessible) and 0.0
 *          (inaccessible)
 */
VEXTERNC double Vacc_splineAccAtom(
        Vacc *thee, /**< Accessibility object */
        double center[VAPBS_DIM], /**< Probe center coordinates */
        double win, /**< Spline window (&Aring;) */
        double infrad, /**< Inflation radius (&Aring;) for ion access. */
        Vatom *atom /**< Atom */
        );

/** @brief   Report gradient of spline-based accessibility with respect to a
 *           particular atom (see Vpmg_splineAccAtom)
 *
 *  Determine accessibility at a given point, given a collection of atomic
 *  spheres.  Uses Benoit Roux (Im et al, Comp Phys Comm, 111, 59--75, 1998)
 *  definition suitable for force evalation; basically a cubic spline.
 *
 *  @ingroup Vacc
 *  @author  Nathan Baker
 */
VEXTERNC void Vacc_splineAccGradAtomUnnorm(
        Vacc *thee, /**< Accessibility object */
        double center[VAPBS_DIM], /**< Probe center coordinates */
        double win, /**< Spline window (&Aring;) */
        double infrad, /**< Inflation radius (&Aring;) for ion access. */
        Vatom *atom, /**< Atom */
        double *force /**< VAPBS_DIM-vector set to gradient of accessibility */
        );

/** @brief   Report gradient of spline-based accessibility with respect to a
 *           particular atom normalized by the accessibility value due to that
 *           atom at that point (see Vpmg_splineAccAtom)
 *
 *  Determine accessibility at a given point, given a collection of atomic
 *  spheres.  Uses Benoit Roux (Im et al, Comp Phys Comm, 111, 59--75, 1998)
 *  definition suitable for force evalation; basically a cubic spline.
 *
 *  @ingroup Vacc
 *  @author  Nathan Baker
 */
VEXTERNC void Vacc_splineAccGradAtomNorm(
        Vacc *thee, /**< Accessibility object */
        double center[VAPBS_DIM], /**< Probe center coordinates */
        double win, /**< Spline window (&Aring;) */
        double infrad, /**< Inflation radius (&Aring;) for ion access. */
        Vatom *atom, /**< Atom */
        double *force /**< VAPBS_DIM-vector set to gradient of accessibility */
        );

/** @brief   Report gradient of spline-based accessibility with respect to a
 *           particular atom normalized by a 4th order accessibility value due
 *           to that atom at that point (see Vpmg_splineAccAtom)
 *
 *  @ingroup Vacc
 *  @author  Michael Schnieders
 */
VEXTERNC void Vacc_splineAccGradAtomNorm4(
        Vacc *thee, /**< Accessibility object */
        double center[VAPBS_DIM], /**< Probe center coordinates */
        double win, /**< Spline window (&Aring;) */
        double infrad, /**< Inflation radius (&Aring;) for ion access. */
        Vatom *atom, /**< Atom */
        double *force /**< VAPBS_DIM-vector set to gradient of accessibility */
        );

/** @brief   Report gradient of spline-based accessibility with respect to a
*           particular atom normalized by a 3rd order accessibility value due
*           to that atom at that point (see Vpmg_splineAccAtom)
*
*  @ingroup Vacc
*  @author  Michael Schnieders
*/
VEXTERNC void Vacc_splineAccGradAtomNorm3(
        Vacc *thee, /**< Accessibility object */
        double center[VAPBS_DIM], /**< Probe center coordinates */
        double win, /**< Spline window (&Aring;) */
        double infrad, /**< Inflation radius (&Aring;) for ion access. */
        Vatom *atom, /**< Atom */
        double *force /**< VAPBS_DIM-vector set to gradient of accessibility */
        );


/**
 * @brief  Build the solvent accessible surface (SAS) and calculate the
 *         solvent accessible surface area
 * @ingroup Vacc
 * @note  Similar to UHBD FORTRAN routine by Brock Luty
 *        (returns UHBD's asas2)
 * @author  Nathan Baker (original FORTRAN routine by Brock Luty)
 * @return  Total solvent accessible area (A^2)
 */
VEXTERNC double Vacc_SASA(
        Vacc *thee,  /**< Accessibility object */
        double radius  /**< Probe molecule radius (&Aring;) */
        );

/**
 * @brief  Return the total solvent accessible surface area (SASA)
 * @ingroup  Vacc
 * @note  Alias for Vacc_SASA
 * @author  Nathan Baker
 * @return  Total solvent accessible area (A^2)
 */
VEXTERNC double Vacc_totalSASA(
        Vacc *thee,  /**< Accessibility object */
        double radius  /**< Probe molecule radius (&Aring;) */
        );

/**
 * @brief  Return the atomic solvent accessible surface area (SASA)
 * @ingroup  Vacc
 * @note  Alias for Vacc_SASA
 * @author  Nathan Baker
 * @return  Atomic solvent accessible area (A^2)
 */
VEXTERNC double Vacc_atomSASA(
        Vacc *thee,  /**< Accessibility object */
        double radius,  /**< Probe molecule radius (&Aring;) */
        Vatom *atom  /**< Atom of interest */
        );

/**
 * @brief  Get the set of points for this atom's solvent-accessible surface
 * @ingroup  Vacc
 * @author  Nathan Baker
 * @return  Pointer to VaccSurf object for this atom
 */
VEXTERNC VaccSurf* Vacc_atomSASPoints(
        Vacc *thee,  /**< Accessibility object */
        double radius,  /**< Probe molecule radius (&Aring;) */
        Vatom *atom  /**< Atom of interest */
        );

/**
* @brief  Get the derivatve of solvent accessible volume
 * @ingroup  Vacc
 * @author  Jason Wagoner, Nathan Baker
 */
VEXTERNC void Vacc_atomdSAV(
                            Vacc *thee, /**< Acessibility object */
                            double radius, /**< Probe radius (&Aring;) */
                            Vatom *atom, /**< Atom of interest */
                            double *dSA /**< Array holding answers of calc */
                            );

/**
* @brief  Get the derivatve of solvent accessible area
 * @ingroup  Vacc
 * @author  Jason Wagoner, David Gohara, Nathan Baker
 */
VEXTERNC void Vacc_atomdSASA(
                            Vacc *thee, /**< Acessibility object */
                            double dpos, /**< Atom position offset */
                            double radius, /**< Probe radius (&Aring;) */
                            Vatom *atom, /**< Atom of interest */
                            double *dSA /**< Array holding answers of calc */
                            );

/**
* @brief  Testing purposes only
 * @ingroup  Vacc
 * @author  David Gohara, Nathan Baker
 */
VEXTERNC void Vacc_totalAtomdSASA(
                             Vacc *thee, /**< Acessibility object */
                             double dpos, /**< Atom position offset */
                             double radius, /**< Probe radius (&Aring;) */
                             Vatom *atom, /**< Atom of interest */
                             double *dSA /**< Array holding answers of calc */
                             );

/**
* @brief  Total solvent accessible volume
 * @ingroup  Vacc
 * @author  David Gohara, Nathan Baker
 */
VEXTERNC void Vacc_totalAtomdSAV(
                                 Vacc *thee, /**< Acessibility object */
                                 double dpos, /**< Atom position offset */
                                 double radius, /**< Probe radius (&Aring;) */
                                 Vatom *atom, /**< Atom of interest */
                                 double *dSA, /**< Array holding answers of calc */
                                 Vclist *clist /**< clist for this calculation */
                                 );

/**
 * @brief  Return the total solvent accessible volume (SAV)
 * @ingroup  Vacc
 * @note  Alias for Vacc_SAV
 * @author  David Gohara
 * @return  Total solvent accessible volume (A^3)
 */
VEXTERNC double Vacc_totalSAV(
    Vacc *thee,  /**< Accessibility object */
    Vclist *clist, /**< Clist for acc object */
    APOLparm *apolparm,  /**< Apolar parameters -- could be VNULL if none required for this calculation.
        * If VNULL, then default settings are used */
    double radius  /**< Probe molecule radius (&Aring;) */
    );

/**
 * @brief  Return the WCA integral energy
 * @ingroup  Vacc
 * @author  David Gohara
 * @return Success flag
 */
VEXTERNC int Vacc_wcaEnergy(
                           Vacc *thee,  /**< Accessibility object */
                             APOLparm *apolparm,  /**< Apolar calculation parameters */
                             Valist *alist, /**< Alist for acc object */
                             Vclist *clist /**< Clist for acc object */
                             );
/**
 * @brief  Return the WCA integral force
 * @ingroup  Vacc
 * @author  David Gohara
 * @return WCA energy (kJ/mol/A)
 */
VEXTERNC int Vacc_wcaForceAtom(Vacc *thee, /**< Accessibility object */
                              APOLparm *apolparm,  /**< Apolar calculation parameters */
                              Vclist *clist, /**< Clist for acc object */
                              Vatom *atom, /**< Current atom */
                              double *force /**< Force for atom */
                           );

/**	@brief	Calculate the WCA energy for an atom
    @ingroup Vacc
    @author	Dave Gohara and Nathan Baker
    @return	Success flag
    */
VEXTERNC int Vacc_wcaEnergyAtom(
    Vacc *thee,	/**< Accessibility object */
    APOLparm *apolparm,	/**< Apolar calculation parameters */
    Valist *alist,	/**< Atom list */
    Vclist *clist,	/**< Cell list associated with Vacc object */
    int iatom,	/**< Index for atom of interest */
    double *value	/**< Set to energy value */
    );

#endif    /* ifndef _VACC_H_ */
