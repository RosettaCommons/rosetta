c	Alfcx.f		Version 1 6/17/2005	Patrice Koehl

c	This file contains a suite of routine for generating the alpha
c	shape based on a given weighted Delaunay triangulation, and
c	a fixed value of Alpha (0 when we are interested in a union of balls)

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

c	voids.f		Version 2 3/30/2007	Patrice Koehl

c	This subroutine builds the alpha complex based on the weighted
c	Delaunay triangulation

	subroutine cavballs( x, y, z, r, ncav )

	integer	npointmax,ntetra_max,nlink_max,nredmax

	parameter	(npointmax   = MAX_POINT)
	parameter	(ntetra_max  = MAX_TETRA)
	parameter	(nlink_max   = MAX_LINK)
	parameter	(nredmax     = MAX_RED)

	integer i,idx,jdx,count, ia,ib,ic,id

	integer ntetra,npoints,nvertex
	integer tetra(4,ntetra_max)
	integer tetra_neighbour(4,ntetra_max)

	integer*1 tetra_info(ntetra_max)
	integer*1 tetra_nindex(ntetra_max)

	integer   depth(ntetra_max)
	integer*1 isvoid(ntetra_max)
	integer*1 vertex_info(npointmax)
	
	integer	ncav
	real*8	x(MAX_TETRA)
	real*8  y(MAX_TETRA)
	real*8  z(MAX_TETRA)
	real*8  r(MAX_TETRA)
	
	real*8	radius(npointmax),weight(npointmax)
	real*8	coord(3*npointmax)

	real*8	a(3),b(3),c(3),d(3),wa,wb,wc,wd,c_abcd(3),c_rad

	common  /tetra_zone/    ntetra,tetra,tetra_neighbour
	common  /tetra_stat/    tetra_info,tetra_nindex
	common  /voids/			depth,isvoid
	common  /xyz_vertex/		coord,radius,weight
	common  /vertex_zone/   npoints,nvertex,vertex_info
	
	ncav = 0
		do 300 idx = 1,ntetra
			if( btest(tetra_info(idx),1) ) then
				ia = tetra(1,idx)
				ib = tetra(2,idx)
				ic = tetra(3,idx)
				id = tetra(4,idx)
	    		do 400 i = 1,3
					a(i) = coord(3*(ia-1)+i)
					b(i) = coord(3*(ib-1)+i)
					c(i) = coord(3*(ic-1)+i)
					d(i) = coord(3*(id-1)+i)
400	continue
				wa = a(1)**2 + a(2)**2 + a(3)**2 - 2*radius(ia)**2
				wb = b(1)**2 + b(2)**2 + b(3)**2 - 2*radius(ib)**2
				wc = c(1)**2 + c(2)**2 + c(3)**2 - 2*radius(ic)**2
				wd = d(1)**2 + d(2)**2 + d(3)**2 - 2*radius(id)**2
				
! 				write(*,*) "tetra:",ia-4,ib-4,ic-4,id-4
				
				call center4( a, b, c, d, wa/2, wb/2, wc/2, wd/2, c_abcd )
				
				c_rad = 0.0
				c_rad = c_rad + (c_abcd(1)-a(1))**2
				c_rad = c_rad + (c_abcd(2)-a(2))**2
				c_rad = c_rad + (c_abcd(3)-a(3))**2
				c_rad = c_rad**(0.5) - radius(ia)
				
! 				write(*,*) c_abcd, c_rad
				if( c_rad > 0 ) then
					if( c_rad < 3 ) then
						ncav = ncav + 1
						x(ncav) = c_abcd(1)
						y(ncav) = c_abcd(2)
						z(ncav) = c_abcd(3)
						r(ncav) = c_rad
					endif
				endif
		endif						
300   continue
		
	return
	end
