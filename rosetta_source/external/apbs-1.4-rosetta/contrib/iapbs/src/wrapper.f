      program wrapper
c
c apbs interface
c
c $Id: wrapper.f 556 2012-01-10 03:03:33Z rok $
c
c simple reference Fortran code
c this shows how to call iapbs library from a fortran application
c
c reads in a file and performs single APBS calculation
c
      implicit none
      integer rc, apbsdrv, natom, i, j, loop
      character*80 rcsid, finput, pqr
      data rcsid /'$Id: wrapper.f 556 2012-01-10 03:03:33Z rok $'/

c      integer MAXAIM
c      parameter (MAXAIM = 150000)

      integer ierr
      include "wrapper_inc.f"

c disabled for now, rok 2012.12.4
c#ifdef HAVE_MPI_H
c      include "mpif.h"
c#endif
      double precision x(NATOMS), y(NATOMS), z(NATOMS)
      double precision radius(NATOMS), charge(NATOMS)

      double precision esenergy(1), npenergy(1)
      double precision apbsdx(NATOMS), apbsdy(NATOMS), apbsdz(NATOMS)
      double precision apbsqfx(NATOMS), apbsqfy(NATOMS), apbsqfz(NATOMS)
      double precision apbsdbx(NATOMS), apbsdby(NATOMS), apbsdbz(NATOMS)
      double precision apbsnpx(NATOMS), apbsnpy(NATOMS), apbsnpz(NATOMS)
      double precision apbsibx(NATOMS), apbsiby(NATOMS), apbsibz(NATOMS)
      double precision apbsnp(3)

      integer apbs_debug
c local
      double precision  pdie, sdie, srad, swin, temp, gamma, sdens
      double precision smvolume, smsize
      integer nonlin, bcfl, nion, srfm, calcenergy, calcforce
      integer calc_type, nlev, cmeth, ccmeth, fcmeth, chgm
      integer wpot, wchg, wsmol, wkappa, wdiel, rchg, rkappa
      integer watompot, rpot, rdiel
      integer calcnpenergy, calcnpforce

      NAMELIST /apbs/ dime, pdime, cglen, fglen, grid, 
     + nonlin, bcfl, nion, pdie, sdie, srfm, chgm, srad, swin, 
     + temp, gamma, sdens, calc_type, nlev,
     + cmeth, ccmeth, fcmeth, ionq, ionc, ionrr, 
     + calcenergy, calcforce, calcnpenergy, calcnpforce, apbs_debug, 
     + wpot, wchg, wsmol, ispara, pqr, loop, smvolume, smsize,
     + wkappa, wdiel, rchg, rkappa, rdiel, watompot, rpot

      integer dummyi
      character dummyc
      double precision dummyr, maxx, minx, maxy, miny, maxz, minz

c#ifdef HAVE_MPI_H
c      call MPI_INIT(ierr)
c#endif

c     initialization of input data
      do i=1, 3
       grid(i) = 0.
       dime(i) = 0
       glen(i) = 0.
       center(i) = 0.
       pdime(i) = 0
       cglen(i) = 0.
       fglen(i) = 0.
       ccenter(i) = 0.
       fcenter(i) = 0.
      end do

c     defaults
      grid(1) = 0.5
      grid(2) = 0.5
      grid(3) = 0.5

      pdie = 2.0
      sdie = 78.4
      srad = 1.4
      swin = 0.3
      temp = 298.15
      sdens = 10.0
      gamma = 0.105
      smvolume = 10.0
      smsize = 1000.0
      smvolume = 10.0
      smsize = 1000.0

      calc_type = 0
      nlev = 4
      cmeth = 1 
      ccmeth = 1 
      fcmeth = 1 
      chgm = 1
      nonlin = 0
      bcfl = 1
      srfm = 2
      calcenergy = 1
      calcforce = 0
      calcnpenergy = 1
      calcnpforce = 0
      wpot = 0
      wchg = 0
      wsmol = 0
      wkappa = 0
      wdiel= 0
      watompot = 0
      rchg = 0
      rkappa = 0
      rdiel = 0
      rpot = 0

      nion = 0
      do i=1, MAXION
         ionq(i) = 0.
         ionc(i) = 0.
         ionrr(i) = 0.
      end do

      ofrac = 0.3

      apbs_debug = 1
      loop = 1

c     read in APBS parameters (from a file specified as cmd line param)
      call getarg(1, finput)
      write ( *, "('Reading parameter file ', a)") finput
      open(1, file=finput, status="old")
      read(1, nml=apbs)
      close(1)

c     read in PQR data
      print *, 'Reading PQR file ...'
      open(1, file=pqr, status="old")
      natom = 0
      do i = 1, 1000
         read(1,*, end=120) dummyc, dummyi, dummyc, dummyc, dummyi,
     +    x(i), y(i), z(i),
     +    charge(i), radius(i)
         natom = natom + 1
      end do
 120  close(1)

      maxx = x(1)
      minx = x(1)
      maxy = y(1)
      miny = y(1)
      maxz = z(1)
      minz = z(1)
      do i = 1, natom
         if(maxx < x(i)+radius(i)) maxx = x(i)+radius(i)
         if(minx > x(i)-radius(i)) minx = x(i)-radius(i)
         if(maxy < y(i)+radius(i)) maxy = y(i)+radius(i)
         if(miny > y(i)-radius(i)) miny = y(i)-radius(i)
         if(maxz < z(i)+radius(i)) maxz = z(i)+radius(i)
         if(minz > z(i)-radius(i)) minz = z(i)-radius(i)
      end do         

      write(6,'(a, 3f8.3)') 'Mol. dimensions: ', maxx-minx, maxy-miny,
     +     maxz-minz

c if we are doing mg-auto calculate recommended grid values
c including dime, if not specified

      if ((calc_type == 0 .OR. calc_type == 1) .AND. dime(1) == 0) then
         cglen(1) = 1.7 * (maxx-minx)
         cglen(2) = 1.7 * (maxy-miny)
         cglen(3) = 1.7 * (maxz-minz)
         fglen(1) = 20.0 + (maxx-minx)
         fglen(2) = 20.0 + (maxy-miny)
         fglen(3) = 20.0 + (maxz-minz)

         do i = 1, 3
            if (fglen(i) > cglen(i)) cglen(i) = fglen(i)
         end do

         if (dime(1) == 0 ) then
            print *, 'Grid dime not specified, calculating ...'
            do i = 1, 3
               dime(i) = 
     +              32*(int((int(fglen(i)/grid(i)+0.5)-1)/32 + 0.5))+ 1
               if (dime(i) < 33) dime(i) = 33
            end do
         end if


         print *, 'Grid values: '
         write(*, '(a, 3f8.3)') 'fglen: ', fglen(1), fglen(2), fglen(3)
         write(*, '(a, 3f8.3)') 'cglen: ', cglen(1), cglen(2), cglen(3)
         write(*, '(a, 3i4)') 'dime: ', dime(1), dime(2), dime(3)
         write(*, '(a, 3f8.3)') 'grid: ', grid(1), grid(2), grid(3)

         write(*, '(a, f10.3)') 'Required memory (in MB): ', 
     +        dime(1)*dime(2)*dime(3)*200.0/1024/1024

      end if

c     print molecule data
      if (apbs_debug > 2) then
         print *, 'PQR file (x, y, z, charge, radius):'
         do i = 1, natom
            write(*, '(i4, 5f8.3)') i, x(i), y(i), z(i), 
     +           charge(i), radius(i)
         end do
      end if

      i_param(1) = calc_type
      i_param(2) = nlev 
      i_param(3) = cmeth
      i_param(4) = ccmeth
      i_param(5) = fcmeth
      i_param(6) = chgm
      i_param(7) = nonlin
      i_param(8) = bcfl
      i_param(9) = srfm
      i_param(10) = calcenergy
      i_param(11) = calcforce
      i_param(12) = wpot
      i_param(13) = wchg
      i_param(14) = wsmol
      i_param(15) = wkappa
      i_param(16) = wdiel
      i_param(17) = watompot
      i_param(18) = rpot
      i_param(19) = 0
      i_param(20) = calcnpforce
      i_param(21) = calcnpenergy
      i_param(22) = nion
      i_param(23) = rchg
      i_param(24) = rkappa
      i_param(25) = rdiel

      r_param(1) = pdie
      r_param(2) = sdie
      r_param(3) = srad
      r_param(4) = swin
      r_param(5) = temp
      r_param(6) = sdens
      r_param(7) = gamma
      r_param(8) = smvolume
      r_param(9) = smsize

      if (apbs_debug > 2) then
         print *, 'i_param:'
         write(*, '(a, i4)') 'calc_type', i_param(1)
         write(*, '(a, i4)') 'nlev', i_param(2)
         write(*, '(a, i4)') 'cmeth', i_param(3)
         write(*, '(a, i4)') 'ccmeth', i_param(4)
         write(*, '(a, i4)') 'fcmeth', i_param(5)
         write(*, '(a, i4)') 'chgm', i_param(6)
         write(*, '(a, i4)') 'nonlin', i_param(7)
         write(*, '(a, i4)') 'bcfl', i_param(8)
         write(*, '(a, i4)') 'srfm', i_param(9)
         write(*, '(a, i4)') 'calcenergy', i_param(10)
         write(*, '(a, i4)') 'calcforce', i_param(11)
         write(*, '(a, i4)') 'wpot', i_param(12)
         write(*, '(a, i4)') 'wchg', i_param(13)
         write(*, '(a, i4)') 'wsmol', i_param(14)
         write(*, '(a, i4)') 'wkappa', i_param(15)
         write(*, '(a, i4)') 'wdiel', i_param(16)
         write(*, '(a, i4)') 'watompot', i_param(17)
         write(*, '(a, i4)') 'rpot', i_param(18)
         write(*, '(a, i4)') '0', i_param(19)
         write(*, '(a, i4)') 'calcnpforce', i_param(20)
         write(*, '(a, i4)') 'calcnpenergy', i_param(21)
         write(*, '(a, i4)') 'nion', i_param(22)
         write(*, '(a, i4)') 'rchg', i_param(23)
         write(*, '(a, i4)') 'rkappa', i_param(24)
         write(*, '(a, i4)') 'rdiel', i_param(25)

         print *, 'r_param:'
         write(*, '(a, f8.3)') 'pdie', r_param(1)
      end if

c more intialization
      do i = 1, natom
         apbsdx = 0.0
         apbsdy = 0.0
         apbsdz = 0.0
         apbsqfx = 0.0
         apbsqfy = 0.0
         apbsqfz = 0.0
         apbsibx = 0.0
         apbsiby = 0.0
         apbsibz = 0.0
         apbsdbx = 0.0
         apbsdby = 0.0
         apbsdbz = 0.0
         apbsnpx = 0.0
         apbsnpy = 0.0
         apbsnpz = 0.0
      end do
      esenergy = 0.0
      npenergy = 0.0

c OK, now we are ready to call the apbs_driver and start the show
      do j = 1, loop
      rc = apbsdrv(natom,x,y,z,radius,charge,r_param,i_param,grid,dime,
     +     pdime, glen, center, cglen, fglen,
     +     ccenter, fcenter, ofrac, apbs_debug,
     +     ionq, ionc, ionrr,
     +     esenergy,npenergy,
     +     apbsdx, apbsdy, apbsdz,
     +     apbsqfx, apbsqfy, apbsqfz,
     +     apbsibx, apbsiby, apbsibz,
     +     apbsnpx, apbsnpy, apbsnpz,
     +     apbsdbx, apbsdby, apbsdbz)

      print *, "main.f: apbs return code: ", rc

      write(*, '(a, f14.8)'), 'esenergy (kJ/mol): ', esenergy(1)
      write(*, '(a, f14.8)'), 'npenergy (kJ/mol): ', npenergy(1)
      write(*, '(a, f14.8)'), 'esenergy (kcal/mol): ', esenergy(1)/4.184
      write(*, '(a, f14.8)'), 'npenergy (kcal/mol): ', npenergy(1)/4.184

      end do

      apbsnp(1) = 0.0
      apbsnp(2) = 0.0
      apbsnp(3) = 0.0

      if (apbs_debug > 2 .and. calcforce == 1) then
         do i = 1, natom
            apbsnp(1) = apbsnp(1) + apbsnpx(i)
            apbsnp(2) = apbsnp(2) + apbsnpy(i)
            apbsnp(3) = apbsnp(3) + apbsnpz(i)
         end do
         print *, "main: Total force on molecule kJ/(mol/A)"
         write(*, '("qf  ", 3f14.8)') apbsqfx(1), apbsqfy(1), apbsqfz(1)
         write(*, '("ib  ", 3f14.8)') apbsibx(1), apbsiby(1), apbsibz(1)
         write(*, '("db  ", 3f14.8)') apbsdbx(1), apbsdby(1), apbsdbz(1)
         write(*, '("pol ", 3f14.8)') apbsdx(1), apbsdy(1), apbsdz(1)
         write(*, '("apol", 3f14.8)') apbsnp(1), apbsnp(2), apbsnp(3)
         write(*, '("tot ", 3f14.8)') apbsdx(1) + apbsnp(1), 
     +        apbsdy(1) + apbsnp(2), apbsdz(1) + apbsnp(3)
      end if


      if (apbs_debug > 2 .and. calcforce == 2) then
         print *, "main: Comps force on atoms kJ/(mol/A)"
         do i = 1, natom
            write(*, '(a, i5, 3f10.4)') 'qf  ', i, apbsqfx(i), 
     +           apbsqfy(i), apbsqfz(i)
            write(*, '(a, i5, 3f10.4)') 'ib  ', i, apbsibx(i), 
     +           apbsiby(i), apbsibz(i)
            write(*, '(a, i5, 3f10.4)') 'db  ', i, apbsdbx(i), 
     +           apbsdby(i), apbsdbz(i)
            write(*, '(a, i5, 3f10.4)') 'ptot', i, apbsdx(i), 
     +           apbsdy(i), apbsdz(i)
            write(*, '(a, i5, 3f10.4)') 'np  ', i, apbsnpx(i), 
     +           apbsnpy(i), apbsnpz(i)
         end do
      end if


c#ifdef HAVE_MPI_H
c      call MPI_FINALIZE(ierr)
c#endif

      end
