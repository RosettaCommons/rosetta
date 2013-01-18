c
c $Id: wrapper_inc.f 554 2012-01-10 02:45:26Z rok $
c


c-----------------------------------------------------------------------
c
c                 iAPBS variables definition
c
c Similar type of varibles are bunched together (integer/integer and
c real/real) to mimize fortran/C passing incompatability.
c
c For detailed description of individual variables please see
c APBS manual.
c
c-----------------------------------------------------------------------
c      INTEGER(kind=4), allocatable, dimension(:) :: ARRINFO, MAXPROCS
c      INTEGER(kind=4), allocatable, dimension(:) :: ARRERR
c      INTEGER(kind=4) NPROC,MPIINFO4,MPICOMM4,ME4,MPISELF,MPIERRS(1)
c      INTEGER(kind=4) ME,NP,MPIDP,MEROOT,IONE,IERR

      integer :: MAXION, NATOMS
      parameter (MAXION = 10, NATOMS = 1000)

c-----------------------------------------------------------------------
c     int ispara;          /**< 1 => is a parallel calculation,
c                           *  0 => is not */
      integer ispara

c-----------------------------------------------------------------------
      double precision r_param(9)

c-----------------------------------------------------------------------
      integer :: i_param(25)


c-----------------------------------------------------------------------
c   int dime[3];               /**< Grid dimensions */
c   int pdime[3];              /**< Grid of processors to be used in
c                               * calculation */
c
      integer :: dime(3), pdime(3)

c-----------------------------------------------------------------------
c   double grid[3];            /**< Grid spacings */
c   double glen[3];            /**< Grid side lengths. */
c   double center[3];          /**< Grid center. If ispart = 0, then this is
c                               * only meaningful if cmeth = 0.  However, if
c                               * ispart = 1 and cmeth = 0, then this is the
c                               * center of the non-disjoint (overlapping)
c                               * partition.  If ispart = 1 and cmeth = 1, then
c                               * this is the vector that must be added to the
c                               * center of the molecule to give the center of
c                               * the non-disjoint partition.  */c
c   double cglen[3];           /**< Coarse grid side lengths */
c   double fglen[3];           /**< Fine grid side lengths */
c   double ccenter[3];         /**< Coarse grid center.  */
c   double fcenter[3];         /**< Fine grid center.  */
c   double ofrac;              /**< Overlap fraction between procs */
c
      double precision grid(3), glen(3), center(3), cglen(3), fglen(3)
      double precision ccenter(3), fcenter(3), ofrac

c-----------------------------------------------------------------------
c mobile ion definition
c
c   double ionq[MAXION]; /**< Counterion charges (in e) */
c   double ionc[MAXION]; /**< Counterion concentrations (in M) */
c   double ionr[MAXION]; /**< Counterion radii (in A) */
c
      double precision ionq(MAXION), ionc(MAXION), ionrr(MAXION)

c-----------------------------------------------------------------------
c atom
c    double position[3];     /**< Atomic position */
c    double radius;          /**< Atomic radius   */
c    double charge;          /**< Atomic charge   */
c these are defined inside of CHARMM

c-----------------------------------------------------------------------
c radii from wmain get saved in a_radius
c      double precision a_radius(NATOMS)

