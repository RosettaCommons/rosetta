/** @defgroup Vfetk Vfetk class
 *  @brief    FEtk master class (interface between FEtk and APBS)
 */

/**
 *  @file     vfetk.h
 *  @ingroup  Vfetk
 *  @brief    Contains declarations for class Vfetk
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

#ifndef _VFETK_H_
#define _VFETK_H_

#include "apbscfg.h"

#include "maloc/maloc.h"
#include "mc/mc.h"

#include "generic/vhal.h"
#include "generic/vatom.h"
// #include "generic/valist.h"
#include "generic/vpbe.h"
#include "generic/vunit.h"
#include "generic/vgreen.h"
#include "generic/vcap.h"
#include "generic/pbeparm.h"
#include "generic/femparm.h"
#include "fem/vcsm.h"

/**
 * @brief  Linear solver type
 * @ingroup Vfetk
 * @note  Do not change these values; they correspond to settings in FEtk
 */
enum eVfetk_LsolvType {
    VLT_SLU=0,  /**< SuperLU direct solve */
    VLT_MG=1,  /**< Multigrid */
    VLT_CG=2,  /**< Conjugate gradient */
    VLT_BCG=3  /**< BiCGStab */
};

/**
 * @brief  Declare FEMparm_LsolvType type
 * @ingroup  Vfetk
 */
typedef enum eVfetk_LsolvType Vfetk_LsolvType;


/**
 * @brief	Mesh loading operation
 * @ingroup	Vfetk
 */
enum eVfetk_MeshLoad {
    VML_DIRICUBE,  /**< Dirichlet cube */
    VML_NEUMCUBE,  /**< Neumann cube */
    VML_EXTERNAL  /**< External mesh (from socket) */
};

/**
 * @brief  Declare FEMparm_GuessType type
 * @ingroup  Vfetk
 */
typedef enum eVfetk_MeshLoad Vfetk_MeshLoad;

/**
 * @brief  Non-linear solver type
 * @ingroup Vfetk
 * @note  Do not change these values; they correspond to settings in FEtk
 */
enum eVfetk_NsolvType {
    VNT_NEW=0,  /**< Newton solver */
    VNT_INC=1,  /**< Incremental */
    VNT_ARC=2  /**< Psuedo-arclength */
};

/**
 * @brief  Declare FEMparm_NsolvType type
 * @ingroup  Vfetk
 */
typedef enum eVfetk_NsolvType Vfetk_NsolvType;

/**
 * @brief  Initial guess type
 * @ingroup Vfetk
 * @note  Do not change these values; they correspond to settings in FEtk
 */
enum eVfetk_GuessType {
    VGT_ZERO=0,  /**< Zero initial guess */
    VGT_DIRI=1,  /**< Dirichlet boundary condition initial guess */
    VGT_PREV=2  /**< Previous level initial guess */
};

/**
 * @brief  Declare FEMparm_GuessType type
 * @ingroup  Vfetk
 */
typedef enum eVfetk_GuessType Vfetk_GuessType;

/**
 * @brief  Preconditioner type
 * @ingroup Vfetk
 * @note  Do not change these values; they correspond to settings in FEtk
 */
enum eVfetk_PrecType {
    VPT_IDEN=0,  /**< Identity matrix */
    VPT_DIAG=1,  /**< Diagonal scaling */
    VPT_MG=2  /**< Multigrid */
};

/**
 * @brief  Declare FEMparm_GuessType type
 * @ingroup  Vfetk
 */
typedef enum eVfetk_PrecType Vfetk_PrecType;

/**
 *  @ingroup Vfetk
 *  @author  Nathan Baker
 *  @brief   Contains public data members for Vfetk class/module
 *
 *  Many of the routines and macros are borrowed from the main.c driver
 *  (written by Mike Holst) provided with the PMG code.
 *
 */
struct sVfetk {

  Vmem *vmem;  /**< Memory management object */
  Gem *gm;  /**< Grid manager (container class for master vertex
             * and simplex lists as well as prolongation operator for updating
             * after refinement).  */
  AM *am;  /**< Multilevel algebra manager. */
  Aprx *aprx;  /**< Approximation manager. */
  PDE *pde;  /**< FEtk PDE object */
  Vpbe *pbe;  /**< Poisson-Boltzmann object */
  Vcsm *csm;  /**< Charge-simplex map */
  Vfetk_LsolvType lkey;  /**< Linear solver method */
  int lmax;  /**< Maximum number of linear solver iterations */
  double ltol;  /**< Residual tolerance for linear solver */
  Vfetk_NsolvType nkey;  /**< Nonlinear solver method */
  int nmax;  /**< Maximum number of nonlinear solver iterations */
  double ntol;  /**< Residual tolerance for nonlinear solver */
  Vfetk_GuessType gues;  /**< Initial guess method */
  Vfetk_PrecType lprec;  /**< Linear preconditioner */
  int pjac;  /**< Flag to print the jacobians (usually set this to -1,
              * please) */
  PBEparm *pbeparm;  /**<  Generic PB parameters */
  FEMparm *feparm;  /**<  FEM-specific parameters */
  Vhal_PBEType type;  /**< Version of PBE to solve */
  int level;  /**< Refinement level (starts at 0) */

};

/** @typedef Vfetk
 *  @ingroup Vfetk
 *  @brief   Declaration of the Vfetk class as the Vfetk structure */
typedef struct sVfetk Vfetk;

/**
 * @brief  Vfetk LocalVar subclass
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @brief  Contains variables used when solving the PDE with FEtk
 */
struct sVfetk_LocalVar {
    double nvec[VAPBS_DIM];  /**< Normal vector for a simplex face */
    double vx[4][VAPBS_DIM];  /**< Vertex coordinates */
    double xq[VAPBS_DIM];  /**< Quadrature pt */
    double U[MAXV];  /**< Solution value */
    double dU[MAXV][VAPBS_DIM];  /**< Solution gradient */
    double W;  /**< Coulomb regularization term scalar value */
    double dW[VAPBS_DIM];  /**< Coulomb regularization term gradient */
    double d2W;  /**< Coulomb regularization term Laplacia */
    int sType;  /**< Simplex type */
    int fType;  /**< Face type */
    double diel;  /**< Dielectric value */
    double ionacc;  /**< Ion accessibility value */
    double A;  /**< Second-order differential term */
    double F;  /**< RHS characteristic function value */
    double B;  /**< Entire ionic strength term */
    double DB;  /**< Entire ionic strength term derivative */
    double jumpDiel;  /**< Dielectric value on one side of a simplex face */
    Vfetk *fetk;  /**< Pointer to the VFETK object */
    Vgreen *green;  /**< Pointer to a Green's function object */
    int initGreen;  /**< Boolean to designate whether Green's function
                     * has been initialized */
    SS *simp;  /**< Pointer to the latest simplex object; set in initElement()
                *  and delta() */
    VV *verts[4];  /**< Pointer to the latest vertices; set in initElement */
    int nverts;  /**< number of vertices in the simplex */
    double ionConc[MAXION];  /**< Counterion species' concentrations */
    double ionQ[MAXION];  /**< Counterion species' valencies */
    double ionRadii[MAXION];  /**< Counterion species' radii */
    double zkappa2; /**< Ionic strength parameters */
    double zks2; /**< Ionic strength parameters */
    double ionstr; /**< Ionic strength parameters (M) */
    int nion;  /**<  Number of ion species */
    double Fu_v;  /**< Store Fu_v value */
    double DFu_wv;  /**< Store DFu_wv value */
    double delta;  /**< Store delta value */
    double u_D;  /**< Store Dirichlet value */
    double u_T;  /**< Store true value */
};

/**
 *  @ingroup Vfetk
 *  @brief   Declaration of the Vfetk_LocalVar subclass as the Vfetk_LocalVar
 *           structure */
typedef struct sVfetk_LocalVar Vfetk_LocalVar;

#if !defined(VINLINE_VFETK)

    /** @brief   Get a pointer to the Gem (grid manager) object
     *  @ingroup Vfetk
     *  @author  Nathan Baker
     *  @return  Pointer to the Gem (grid manager) object
     */
    VEXTERNC Gem* Vfetk_getGem(
            Vfetk *thee /**< Vfetk object */
            );

    /** @brief   Get a pointer to the AM (algebra manager) object
     *  @ingroup Vfetk
     *  @author  Nathan Baker
     *  @return  Pointer to the AM (algebra manager) object
     */
    VEXTERNC AM* Vfetk_getAM(
            Vfetk *thee /**< The Vfetk object */
            );

    /** @brief   Get a pointer to the Vpbe (PBE manager) object
     *  @ingroup Vfetk
     *  @author  Nathan Baker
     *  @return  Pointer to the Vpbe (PBE manager) object
     */
    VEXTERNC Vpbe* Vfetk_getVpbe(
            Vfetk *thee /**< The Vfetk object */
            );

    /** @brief   Get a pointer to the Vcsm (charge-simplex map) object
     *  @ingroup Vfetk
     *  @author  Nathan Baker
     *  @return  Pointer to the Vcsm (charge-simplex map) object
     */
    VEXTERNC Vcsm* Vfetk_getVcsm(
            Vfetk *thee /**< The Vfetk object */
            );

    /** @brief   Get the partition information for a particular atom
     *  @ingroup Vfetk
     *  @author  Nathan Baker
     *  @note    Friend function of Vatom
     *  @returns Partition ID
     */
    VEXTERNC int Vfetk_getAtomColor(
            Vfetk *thee, /**< The Vfetk object */
            int iatom  /**< Valist atom index */
            );

#else /* if defined(VINLINE_VFETK) */
#   define Vfetk_getGem(thee) ((thee)->gm)
#   define Vfetk_getAM(thee) ((thee)->am)
#   define Vfetk_getVpbe(thee) ((thee)->pbe)
#   define Vfetk_getVcsm(thee) ((thee)->csm)
#   define Vfetk_getAtomColor(thee, iatom) (Vatom_getPartID(Valist_getAtom(Vpbe_getValist(thee->pbe), iatom)))
#endif /* if !defined(VINLINE_VFETK) */

/* ///////////////////////////////////////////////////////////////////////////
// Class Vfetk: Non-Inlineable methods (vfetk.c)
/////////////////////////////////////////////////////////////////////////// */

/**
 * @brief  Constructor for Vfetk object
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @return  Pointer to newly allocated Vfetk object
 * @note  This sets up the Gem, AM, and Aprx FEtk objects but does not create
 *         a mesh.  The easiest way to create a mesh is to then call
 *         Vfetk_genCube
 */
VEXTERNC Vfetk* Vfetk_ctor(
        Vpbe *pbe, /**< Vpbe (PBE manager object) */
        Vhal_PBEType type /**< Version of PBE to solve */
        );

/**
 * @brief  FORTRAN stub constructor for Vfetk object
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @return  1 if successful, 0 otherwise
 * @note  This sets up the Gem, AM, and Aprx FEtk objects but does not create
 *         a mesh.  The easiest way to create a mesh is to then call
 *         Vfetk_genCube
 */
VEXTERNC int Vfetk_ctor2(
        Vfetk *thee, /**< Vfetk object memory */
        Vpbe *pbe, /**< PBE manager object */
        Vhal_PBEType type /**< Version of PBE to solve */
        );

/**
 * @brief   Object destructor
 * @ingroup Vfetk
 * @author  Nathan Baker
 */
VEXTERNC void Vfetk_dtor(
        Vfetk **thee /**< Pointer to memory location of Vfetk object */
        );

/**
 * @brief   FORTRAN stub object destructor
 * @ingroup Vfetk
 * @author  Nathan Baker
 */
VEXTERNC void Vfetk_dtor2(
        Vfetk *thee /**< Pointer to Vfetk object to be destroyed */
        );

/**
 * @brief   Create an array containing the solution (electrostatic potential
 *          in units of \f$k_B T/e\f$) at the finest mesh level.
 * @ingroup Vfetk
 * @author  Nathan Baker and Michael Holst
 * @note    The user is responsible for destroying the newly created array
 * @return  Newly created array of length "length" (see above); the user is
 *           responsible for destruction
 */
VEXTERNC double* Vfetk_getSolution(
        Vfetk *thee, /**< Vfetk object with solution */
        int *length /**< Ste to length of the newly created solution array */
        );

/**
 * @brief  Set the parameter objects
 * @ingroup  Vfetk
 * @author  Nathan Baker
 */
VEXTERNC void Vfetk_setParameters(
        Vfetk *thee, /**< The Vfetk object */
        PBEparm *pbeparm, /**< Parameters for solution of the PBE */
        FEMparm *feparm /**< FEM-speecific solution parameters */
        );

/**
 * @brief   Return the total electrostatic energy
 *
 *  Using the solution at the finest mesh level, get the electrostatic energy
 *  using the free energy functional for the Poisson-Boltzmann equation
 *  without removing any self-interaction terms (i.e., removing the reference
 *  state of isolated charges present in an infinite dielectric continuum with
 *  the same relative permittivity as the interior of the protein) and return
 *  the result in units of \f$k_B T\f$.  The argument color allows the user to
 *  control the partition on which this energy is calculated; if (color == -1)
 *  no restrictions are used.  The solution is obtained from the finest level
 *  of the passed AM object, but atomic data from the Vfetk object is used to
 *  calculate the energy.
 *
 * @ingroup Vfetk
 * @author  Nathan Baker
 * @return Total electrostatic energy in units of \f$k_B T\f$.
 */
VEXTERNC double  Vfetk_energy(
        Vfetk *thee, /**< THe Vfetk object */
        int color, /**< Partition restriction for energy calculation; if
                     non-negative, energy calculation is restricted to the
                     specified partition (indexed by simplex and atom colors
                     */
        int nonlin /**< If 1, the NPBE energy functional is used; otherwise,
                     the LPBE energy functional is used. If -2, SMPBE is used. */
        );

/**
 * @brief   Get the "mobile charge" and "polarization" contributions to the
 *          electrostatic energy.
 *
 * Using the solution at the finest mesh level, get the
 * electrostatic energy due to the interaction of the mobile charges
 * with the potential and polarization of the dielectric medium:
 *   \f[ G = \frac{1}{4 I_s} \sum_i c_i q_i^2 \int
 *   \overline{\kappa}^2(x) e^{-q_i u(x)} dx + \frac{1}{2} \int
 *   \epsilon ( \nabla u )^2 dx \f]
 * for the NPBE and
 *   \f[ G = \frac{1}{2} \int \overline{\kappa}^2(x) u^2(x) dx +
 *   \frac{1}{2} \int \epsilon ( \nabla u )^2 dx \f]
 * for the LPBE.  Here \f$i\f$ denotes the counterion species,
 * \f$I_s\f$ is the bulk ionic strength, \f$\overline{\kappa}^2(x)\f$
 * is the modified Debye-Huckel parameter, \f$c_i\f$ is the
 * concentration of species \f$i\f$, \f$q_i\f$ is the charge of
 * species \f$i\f$, \f$\epsilon\f$ is the dielectric function, and
 * \f$u(x)\f$ is the dimensionless electrostatic potential.  The
 * energy is scaled to units of \f$k_b T\f$.
 *
 * @ingroup Vfetk
 * @author  Nathan Baker
 * @param thee  Vfetk object
 * @param color Partition restriction for energy evaluation, only used if
 *               non-negative
 * @return The "mobile charge" and "polarization" contributions to the
 *          electrostatic energy in units of \f$k_B T\f$.
 */
VEXTERNC double  Vfetk_dqmEnergy(
        Vfetk *thee, /**< The Vfetk object */
        int color  /**< Partition restriction for energy calculation; if
                     non-negative, energy calculation is restricted to the
                     specified partition (indexed by simplex and atom colors
                     */
        );

/**
 * @brief   Get the "fixed charge" contribution to the electrostatic energy
 *
 *          Using the solution at the finest mesh level, get the
 *          electrostatic energy due to the interaction of the fixed charges
 *          with the potential: \f[ G = \sum_i q_i u(r_i) \f]
 *          and return the result in units of \f$k_B T\f$.  Clearly, no
 *          self-interaction terms are removed.  A factor a 1/2 has to be
 *          included to convert this to a real energy.
 *
 * @ingroup Vfetk
 * @author  Nathan Baker
 * @param   thee   Vfetk object
 * @param   color Partition restriction for energy evaluation, only used if
 *              non-negative
 *  @returns The fixed charge electrostatic energy in units of \f$k_B T\f$.
 */
VEXTERNC double Vfetk_qfEnergy(
        Vfetk *thee, /**< The Vfetk object */
        int color /**< Partition restriction for energy evaluation, only used
                    if non-negative */
        );

/**
 * @brief   Return the memory used by this structure (and its contents)
 *          in bytes
 * @ingroup Vfetk
 * @author  Nathan Baker
 * @return  The memory used by this structure and its contents in bytes
 */
VEXTERNC unsigned long int Vfetk_memChk(
        Vfetk *thee /**< THe Vfetk object */
        );

/**
 * @brief   Transfer color (partition ID) information frmo a partitioned mesh
 *          to the atoms.
 *
 *          Transfer color information from partitioned mesh to the atoms.
 *          In the case that a charge is shared between two partitions, the
 *          partition color of the first simplex is selected.  Due to the
 *          arbitrary nature of this selection, THIS METHOD SHOULD ONLY BE
 *          USED IMMEDIATELY AFTER PARTITIONING!!!
 * @warning This function should only be used immediately after mesh
 *          partitioning
 * @ingroup Vfetk
 * @author  Nathan Baker
 * @note    This is a friend function of Vcsm
 */
VEXTERNC void Vfetk_setAtomColors(
        Vfetk *thee /**< THe Vfetk object */
        );

/**
 * @brief   Writes a Bmat to disk in Harwell-Boeing sparse matrix format.
 *
 * @ingroup Vfetk
 * @author  Stephen Bond
 * @note    This is a friend function of Bmat
 * @bug     Hardwired to only handle the single block symmetric case.
 */
VEXTERNC void Bmat_printHB(
        Bmat *thee, /**< The matrix to write */
        char *fname /**< Filename for output */
        );

/**
 * @brief	Construct a rectangular mesh (in the current Vfetk object)
 * @ingroup Vfetk
 * @author	Nathan Baker
 */
VEXTERNC Vrc_Codes Vfetk_genCube(
                                 Vfetk *thee,  /**< Vfetk object */
                                 double center[3],  /**< Center for mesh */
                                 double length[3],  /**< Mesh lengths */
                                 Vfetk_MeshLoad meshType  /**< Mesh boundary conditions */
                                 );

/**
 * @brief	Loads a mesh into the Vfetk (and associated) object(s).
 * @ingroup	Vfetk
 * @author	Nathan Baker
 */
VEXTERNC Vrc_Codes Vfetk_loadMesh(
                                  Vfetk *thee, /**< Vfetk object to load into */
                                  double center[3],  /**< Center for mesh (if constructed) */
                                  double length[3],  /**< Mesh lengths (if constructed) */
                                  Vfetk_MeshLoad meshType,  /**< Type of mesh to load */
                                  Vio *sock  /**< Socket for external mesh data (NULL otherwise) */
                                  );

/**
 * @brief  Constructs the FEtk PDE object
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @return  Newly-allocated PDE object
 * @bug  Not thread-safe */
VEXTERNC PDE* Vfetk_PDE_ctor(
        Vfetk *fetk /**< The Vfetk object */
        );

/**
 * @brief  Intializes the FEtk PDE object
 * @ingroup  Vfetk
 * @author  Nathan Baker (with code by Mike Holst)
 * @return  1 if successful, 0 otherwise
 * @bug  Not thread-safe */
VEXTERNC int Vfetk_PDE_ctor2(
        PDE *thee, /**< The newly-allocated PDE object */
        Vfetk *fetk /**< The parent Vfetk object */
        );

/**
 * @brief  Destroys FEtk PDE object
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @note  Thread-safe
 */
VEXTERNC void Vfetk_PDE_dtor(
        PDE **thee /**< Pointer to PDE object memory */
        );

/**
 * @brief  FORTRAN stub:  destroys FEtk PDE object
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @note  Thread-safe
 */
VEXTERNC void Vfetk_PDE_dtor2(
        PDE *thee /**< PDE object memory */
        );

/**
 * @brief  Do once-per-assembly initialization
 * @ingroup  Vfetk
 * @author  Nathan Baker and Mike Holst
 * @note  Thread-safe */
VEXTERNC void Vfetk_PDE_initAssemble(
        PDE *thee, /**< PDE object */
        int ip[], /**< Integer parameter array (not used) */
        double rp[] /**< Double parameter array (not used) */
        );

/**
 * @brief  Do once-per-element initialization
 * @ingroup  Vfetk
 * @author  Nathan Baker and Mike Holst
 * @todo  Jump term is not implemented
 * @bug This function is not thread-safe */
VEXTERNC void Vfetk_PDE_initElement(
        PDE *thee,  /**< PDE object */
        int elementType,  /**< Material type (not used) */
        int chart,  /**< Chart in which the vertex coordinates are provided,
                      used here as a bitfield to store molecular accessibility
                      */
        double tvx[][VAPBS_DIM],  /**< Vertex coordinates */
        void *data /**< Simplex pointer (hack) */
        );

/**
 * @brief  Do once-per-face initialization
 * @ingroup  Vfetk
 * @author  Nathan Baker and Mike Holst
 * @bug This function is not thread-safe */
VEXTERNC void Vfetk_PDE_initFace(
        PDE *thee, /**< THe PDE object */
        int faceType, /**< Simplex face type (interior or various boundary
                        types) */
        int chart, /**< Chart in which the vertex coordinates are provided,
                     used here as a bitfield for molecular accessibility */
        double tnvec[] /**< Coordinates of outward normal vector for face */
        );

/**
 * @brief  Do once-per-point initialization
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @bug This function is not thread-safe
 * @bug This function uses pre-defined boudnary definitions for the molecular
 *      surface. */
VEXTERNC void Vfetk_PDE_initPoint(
        PDE *thee,  /**< The PDE object */
        int pointType, /**< The type of point -- interior or various faces */
        int chart,  /**< The chart in which the point coordinates are provided,
                      used here as bitfield for molecular accessibility */
        double txq[],  /**< Point coordinates */
        double tU[],  /**< Solution value at point */
        double tdU[][VAPBS_DIM] /**< Solution derivative at point */
        );

/**
 * @brief  Evaluate strong form of PBE.  For interior points, this is:
 *  \f[ -\nabla \cdot \epsilon \nabla u + b(u) - f \f]
 *  where \f$b(u)\f$ is the (possibly nonlinear) mobile ion term and \f$f\f$ is
 *  the source charge distribution term (for PBE) or the induced surface charge
 *  distribution (for RPBE).  For an interior-boundary (simplex face) point,
 *  this is:
 *  \f[ [\epsilon(x) \nabla u(x) \cdot n(x)]_{x=0^+} - [\epsilon(x) \nabla u(x)
 *  \cdot n(x)]_{x=0^-} \f]
 *  where \f$n(x)\f$ is the normal to the simplex face and the term represents
 *  the jump in dielectric displacement across the face.  There is no
 *  outer-boundary contribution for this problem.
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @bug This function is not thread-safe
 * @bug This function is not implemented (sets error to zero)
 */
VEXTERNC void Vfetk_PDE_Fu(
        PDE *thee, /**< The PDE object */
        int key, /**< Type of point (0 = interior, 1 = boundary, 2 = interior
                   boundary */
        double F[] /**< Set to value of residual */
        );

/**
 * @brief  This is the weak form of the PBE; i.e. the strong form integrated
 * with a test function to give:
 * \f[ \int_\Omega \left[ \epsilon \nabla u \cdot \nabla v + b(u) v - f v
 * \right] dx \f]
 * where \f$b(u)\f$ denotes the mobile ion term.
 * @ingroup  Vfetk
 * @author  Nathan Baker and Mike Holst
 * @return  Integrand value
 * @bug  This function is not thread-safe */
VEXTERNC double Vfetk_PDE_Fu_v(
        PDE *thee, /**< The PDE object */
        int key, /**< Integrand to evaluate (0 = interior weak form, 1 =
                   boundary weak form */
        double V[],  /**< Test function at current point */
        double dV[][VAPBS_DIM] /**< Test function derivative at current point */
        );

/**
 * @brief  This is the linearization of the weak form of the PBE; e.g., for use
 * in a Newton iteration. This is the functional linearization of the strong
 * form integrated with a test function to give:
 * \f[ \int_\Omega \left[ \epsilon \nabla w \cdot \nabla v + b'(u) w v - f v
 * \right] dx \f]
 * where \f$b'(u)\f$ denotes the functional derivation of the mobile ion term.
 * @ingroup  Vfetk
 * @author  Nathan Baker and Mike Holst
 * @return  Integrand value
 * @bug  This function is not thread-safe */
VEXTERNC double Vfetk_PDE_DFu_wv(
        PDE *thee, /**< The PDE object */
        int key, /**< Integrand to evaluate (0 = interior weak form, 1 =
                   boundary weak form) */
        double W[], /**< Trial function value at current point */
        double dW[][VAPBS_DIM], /**< Trial function gradient at current point */
        double V[], /**< Test function value at current point */
        double dV[][VAPBS_DIM] /**< Test function gradient */
        );

/**
 * @brief  Evaluate a (discretized) delta function source term at the given
 * point
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @bug  This function is not thread-safe */
VEXTERNC void Vfetk_PDE_delta(
        PDE *thee, /**< PDE object */
        int type, /**< Vertex type */
        int chart, /**< Chart for point coordinates */
        double txq[], /**< Point coordinates */
        void *user, /**< Vertex object pointer */
        double F[] /**< Set to delta function value */
        );

/**
 * @brief  Evaluate the Dirichlet boundary condition at the given point
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @bug  This function is hard-coded to call only multiple-sphere
 * Debye-H&uuml; functions.
 * @bug  This function is not thread-safe. */
VEXTERNC void Vfetk_PDE_u_D(
        PDE *thee, /**< PDE object */
        int type, /**< Vertex boundary type */
        int chart, /**< Chart for point coordinates */
        double txq[], /**< Point coordinates */
        double F[] /**< Set to boundary values */
        );

/**
 * @brief  Evaluate the "true solution" at the given point for comparison with
 * the numerical solution
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @note  This function only returns zero.
 * @bug  This function is not thread-safe. */
VEXTERNC void Vfetk_PDE_u_T(
        PDE *thee, /**< PDE object */
        int type, /**< Point type */
        int chart, /**< Chart for point coordinates */
        double txq[], /**< Point coordinates */
        double F[] /**< Set to value at point */
        );

/**
 * @brief  Define the way manifold edges are bisected
 * @ingroup  Vfetk
 * @author  Nathan Baker and Mike Holst
 * @note  This function is thread-safe. */
VEXTERNC void Vfetk_PDE_bisectEdge(
        int dim, /**< Intrinsic dimension of manifold */
        int dimII, /**< Embedding dimension of manifold */
        int edgeType,  /**< Type of edge being refined */
        int chart[], /**< Chart for edge vertices, used here as accessibility
                       bitfields */
        double vx[][VAPBS_DIM] /**< Edge vertex coordindates */
        );

/**
 * @brief  Map a boundary point to some pre-defined shape
 * @ingroup  Vfetk
 * @author  Nathan Baker and Mike Holst
 * @note  This function is thread-safe and is a no-op */
VEXTERNC void Vfetk_PDE_mapBoundary(
        int dim, /**< Intrinsic dimension of manifold */
        int dimII, /**< Embedding dimension of manifold */
        int vertexType,  /**< Type of vertex */
        int chart, /**< Chart for vertex coordinates */
        double vx[VAPBS_DIM] /**< Vertex coordinates */
        );

/**
 * @brief  User-defined error estimator -- in our case, a geometry-based
 * refinement method; forcing simplex refinement at the dielectric boundary and
 * (for non-regularized PBE) the charges.
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @return 1 if mark simplex for refinement, 0 otherwise
 * @bug  This function is not thread-safe */
VEXTERNC int Vfetk_PDE_markSimplex(
        int dim, /**< Intrinsic manifold dimension */
        int dimII, /**< Embedding manifold dimension */
        int simplexType,  /**< Type of simplex being refined */
        int faceType[VAPBS_NVS], /**< Types of faces in simplex */
        int vertexType[VAPBS_NVS], /**< Types of vertices in simplex */
        int chart[], /**< Charts for vertex coordinates */
        double vx[][VAPBS_DIM], /**< Vertex coordinates */
        void *simplex /**< Simplex pointer */
        );

/**
 * @brief  Unify the chart for different coordinate systems -- a no-op for us.
 * @ingroup  Vfetk
 * @author Nathan Baker
 * @note  Thread-safe; a no-op */
VEXTERNC void Vfetk_PDE_oneChart(
        int dim, /**< Intrinsic manifold dimension */
        int dimII, /**< Embedding manifold dimension */
        int objType, /**< ??? */
        int chart[], /**< Charts of vertices' coordinates */
        double vx[][VAPBS_DIM], /**< Vertices' coordinates */
        int dimV /**< Number of vertices */
        );

/**
 * @brief  Energy functional.  This returns the energy (less delta function
 * terms) in the form:
 *   \f[ c^{-1}/2 \int (\epsilon (\nabla u)^2 + \kappa^2 (cosh u - 1)) dx \f]
 * for a 1:1 electrolyte where \f$c\f$ is the output from Vpbe_getZmagic.
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @returns  Energy value (in kT)
 * @bug  This function is not thread-safe. */
VEXTERNC double Vfetk_PDE_Ju(
        PDE *thee, /**< The PDE object */
        int key /**< What to evluate:  interior (0) or boundary (1)? */
        );

/**
 * @brief  External hook to simplex subdivision routines in Gem.  Called each
 * time a simplex is subdivided (we use it to update the charge-simplex map)
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @bug  This function is not thread-safe.
 */
VEXTERNC void Vfetk_externalUpdateFunction(
        SS **simps, /**< List of parent (simps[0]) and children (remainder)
                      simplices */
        int num /**< Number of simplices in list */
        );


/**
 * @brief Initialize the bases for the trial or the test space, for a
 * particular component of the system, at all quadrature points on the master
 * simplex element.
 * @ingroup Vfetk
 * @author Mike Holst
 * @note
 *   @verbatim
 *   The basis ordering is important.  For a fixed quadrature
 *   point iq, you must follow the following ordering in p[iq][],
 *   based on how you specify the degrees of freedom in dof[]:
 *
 *   <v_0 vDF_0>,      <v_1 vDF_0>,      ..., <v_{nv} vDF_0>
 *   <v_0 vDF_1>,      <v_1 vDF_1>,      ..., <v_{nv} vDF_1>
 *                           ...
 *   <v_0 vDF_{nvDF}>, <v_0 vDF_{nvDF}>, ..., <v_{nv} vDF_{nvDF}>
 *
 *   <e_0 eDF_0>,      <e_1 eDF_0>,      ..., <e_{ne} eDF_0>
 *   <e_0 eDF_1>,      <e_1 eDF_1>,      ..., <e_{ne} eDF_1>
 *                           ...
 *   <e_0 eDF_{neDF}>, <e_1 eDF_{neDF}>, ..., <e_{ne} eDF_{neDF}>
 *
 *   <f_0 fDF_0>,      <f_1 fDF_0>,      ..., <f_{nf} fDF_0>
 *   <f_0 fDF_1>,      <f_1 fDF_1>,      ..., <f_{nf} fDF_1>
 *                           ...
 *   <f_0 fDF_{nfDF}>, <f_1 fDF_{nfDF}>, ..., <f_{nf} fDF_{nfDF}>
 *
 *   <s_0 sDF_0>,      <s_1 sDF_0>,      ..., <s_{ns} sDF_0>
 *   <s_0 sDF_1>,      <s_1 sDF_1>,      ..., <s_{ns} sDF_1>
 *                           ...
 *   <s_0 sDF_{nsDF}>, <s_1 sDF_{nsDF}>, ..., <s_{ns} sDF_{nsDF}>
 *
 *   For example, linear elements in R^3, with one degree of freedom at each *
 *   vertex, would use the following ordering:
 *
 *     <v_0 vDF_0>, <v_1 vDF_0>, <v_2 vDF_0>, <v_3 vDF_0>
 *
 *   Quadratic elements in R^2, with one degree of freedom at each vertex and
 *   edge, would use the following ordering:
 *
 *     <v_0 vDF_0>, <v_1 vDF_0>, <v_2 vDF_0>
 *     <e_0 eDF_0>, <e_1 eDF_0>, <e_2 eDF_0>
 *
 *   You can use different trial and test spaces for each component of the
 *   elliptic system, thereby allowing for the use of Petrov-Galerkin methods.
 *   You MUST then tag the bilinear form symmetry entries as nonsymmetric in
 *   your PDE constructor to reflect that DF(u)(w,v) will be different from
 *   DF(u)(v,w), even if your form acts symmetrically when the same basis is
 *   used for w and v.
 *
 *   You can also use different trial spaces for each component of the elliptic
 *   system, and different test spaces for each component of the elliptic
 *   system.  This allows you to e.g.  use a basis which is vertex-based for
 *   one component, and a basis which is edge-based for another.  This is
 *   useful in fluid mechanics, eletromagnetics, or simply to play around with
 *   different elements.
 *
 *   This function is called by MC to build new master elements whenever it
 *   reads in a new mesh.  Therefore, this function does not have to be all
 *   that fast, and e.g.  could involve symbolic computation.
 *   @endverbatim
 */
VEXTERNC int Vfetk_PDE_simplexBasisInit(
        int key, /**< Basis type to evaluate (0 = trial, 1 = test, 2 = trialB,
                   3 = testB) */
        int dim, /**< Spatial dimension */
        int comp, /**< Which component of elliptic system to produce basis
                    for?  */
        int *ndof,  /**< Set to the number of degrees of freedom */
        int dof[] /**< Set to degree of freedom per v/e/f/s */
        );

/**
 * @brief Evaluate the bases for the trial or test space, for a particular
 * component of the system, at all quadrature points on the master simplex
 * element.
 * @ingroup Vfetk
 * @author Mike Holst
 */
VEXTERNC void Vfetk_PDE_simplexBasisForm(
        int key, /**< Basis type to evaluate (0 = trial, 1 = test, 2 = trialB,
                   3 = testB) */
        int dim, /**< Spatial dimension */
        int comp /**< Which component of elliptic system to produce basis for */,
        int pdkey, /**< Basis partial differential equation evaluation key:
                     \li 0 = evaluate basis(x,y,z)
                     \li 1 = evaluate basis_x(x,y,z)
                     \li 2 = evaluate basis_y(x,y,z)
                     \li 3 = evaluate basis_z(x,y,z)
                     \li 4 = evaluate basis_xx(x,y,z)
                     \li 5 = evaluate basis_yy(x,y,z)
                     \li 6 = evaluate basis_zz(x,y,z)
                     \li 7 = etc... */
        double xq[], /**< Set to quad pt coordinate */
        double basis[] /**< Set to all basis functions evaluated at all
                         quadrature pts */
        );

/**
 * @brief  Read in mesh and initialize associated internal structures
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @note  @see Vfetk_genCube */
VEXTERNC void Vfetk_readMesh(
        Vfetk *thee, /**< THe Vfetk object */
        int skey, /**< The sock format key (0 = MCSF simplex format) */
        Vio *sock /**< Socket object ready for reading */
        );

/**
 * @brief  Debugging routine to print out local variables used by PDE object
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @bug  This function is not thread-safe */
VEXTERNC void Vfetk_dumpLocalVar();

/**
 * @brief  Fill an array with the specified data
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @note  This function is thread-safe
 * @bug  Several values of type are not implemented
 * @return  1 if successful, 0 otherwise */
VEXTERNC int Vfetk_fillArray(
        Vfetk *thee, /**< The Vfetk object with the data */
        Bvec *vec, /**< The vector to hold the data */
        Vdata_Type type /**< THe type of data to write */
        );

/**
 * @brief  Write out data
 * @ingroup  Vfetk
 * @author  Nathan Baker
 * @param  thee  Vfetk object
 * @param  vec  FEtk Bvec vector to use
 * @param  format  Format for data
 * @param iodev  Output device type (FILE/BUFF/UNIX/INET)
 * @param iofmt  Output device format (ASCII/XDR)
 * @param thost  Output hostname (for sockets)
 * @param fname  Output FILE/BUFF/UNIX/INET name
 * @note  This function is thread-safe
 * @bug  Some values of format are not implemented
 * @return  1 if successful, 0 otherwise */
VEXTERNC int Vfetk_write(
        Vfetk *thee, /**< The Vfetk object */
        const char *iodev, /**< Output device type (FILE = file, BUFF = buffer,
                             UNIX = unix pipe, INET = network socket) */
        const char *iofmt, /**< Output device format (ASCII = ascii/plaintext,
                             XDR = xdr) */
        const char *thost, /**< Output hostname for sockets */
        const char *fname, /**< Output filename for other */
        Bvec *vec, /**< Data vector */
        Vdata_Format format /**< Data format */
        );

/**
 * @brief  Load a Gem geometry manager object into Vfetk
 * @ingroup Vfetk
 * @author  Nathan Baker
 */
VEXTERNC Vrc_Codes Vfetk_loadGem(
                                 Vfetk *thee, /**< Destination */
                                 Gem *gm /**< Geometry manager source */
                                 );


#endif /* ifndef _VFETK_H_ */
