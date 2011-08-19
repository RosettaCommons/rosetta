c	Adjust.f		Version 1 6/17/2005	Patrice Koehl

c	This file contains a set of routines that adjust the number of
c	spheres included in the calculation to be always at least 4,
c	such that the regular triangulation can be computed, and the
c	alpha shape derived

c	Copyright (C) 2005 Patrice Koehl

c	This library is free software; you can redistribute it and/or
c	modify it under the terms of the GNU Lesser General Public
c	License as published by the Free Software Foundation; either
c	version 2.1 of the License, or (at your option) any later version.

c	This library is distributed in the hope that it will be useful,
c	but WITHOUT ANY WARRANTY; without even the implied warranty of
c	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
c	Lesser General Public License for more details.

c	You should have received a copy of the GNU Lesser General Public
c	License along with this library; if not, write to the Free Software
c	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

#include "defines.h"

c	adjust_nsphere.f		Version 1 6/17/2005	Patrice Koehl

c	This subroutine adds spheres, if needed

	subroutine adjust_nsphere

	integer	npointmax
	parameter (npointmax = MAX_POINT)

	integer i,j,k
	integer	nvertex,npoint
	integer	sign(3,3)

	integer*1	redinfo(npointmax)

c	Information on the vertices

	real*8	Rmax
	real*8	Dmax(3)
	real*8	radius(npointmax)
	real*8	coord(3*npointmax)
	real*8	coord4(npointmax)

	data ((sign(i,j),j=1,3),i=1,3) /1,1,1,1,-1,-1,-1,1,1/

	common /xyz_vertex/ coord,radius,coord4
	common /vertex_zone/ npoint,nvertex,redinfo

	save


c	Do nothing if we already have at least 4 balls

	if(npoint.ge.4) return

c	Get bounding box of the current set of spheres, as well as max
c	radius

	do 100 i = 1,3
		Dmax(i) = coord(i)
100	continue
	Rmax = radius(1)

	do 300 i = 2,npoint
		do 200 j = 1,3
			if(Dmax(j).lt.coord(3*(i-1)+j)) then
				Dmax(j) = coord(3*(i-1)+j)
			endif
200		continue
		if(Rmax.lt.radius(i)) Rmax = radius(i)
300	continue

c	Now add point(s) with center at Dmax + 3*Rmax, and radius Rmax/20

	do 500 i = npoint+1,4
		j = i - npoint
		do 400 k = 1,3
			coord(3*(i-1)+k)=sign(j,k)*(Dmax(k)+2*Rmax)
400		continue
		radius(i) = Rmax/20
500	continue

	npoint = 4

	return
	end
c	readjust_nsphere.f		Version 1 6/17/2005	Patrice Koehl

c	This subroutine removes the artificial spheres, if needed

	subroutine readjust_nsphere(nsphere,nred,listred)

	integer	npointmax
	parameter (npointmax = MAX_POINT)

	integer i,j
	integer	npoints,nvertex,nsphere,nred

	integer	listred(nred)

c	Information on the vertices

	integer*1 redinfo(npointmax)

	common  /vertex_zone/   npoints,nvertex,redinfo

	save

c	Do nothing if we already have at least 4 balls

	if(nsphere.ge.4) return

	do 100 i = nsphere+5,8
		redinfo(i) = 1
100	continue

	npoints = nsphere
	nvertex = npoints + 4

	j = 0
	do 200 i = 1,nred
		if(listred(i).le.nsphere) then
			j = j + 1
			listred(j) = listred(i)
		endif
200	continue

	return
	end
