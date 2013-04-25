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

	subroutine markvoids

	parameter	(npointmax   = MAX_POINT)
	parameter	(ntetra_max  = MAX_TETRA)
	parameter	(nlink_max   = MAX_LINK)
	parameter	(nredmax     = MAX_RED)

	integer i,idx,jdx,count

	integer ntetra
	integer tetra(4,ntetra_max)
	integer tetra_neighbour(4,ntetra_max)

	integer*1 tetra_info(ntetra_max)
	integer*1 tetra_nindex(ntetra_max)

	integer   depth(ntetra_max)
	integer*1 isvoid(ntetra_max)

	common  /tetra_zone/    ntetra,tetra,tetra_neighbour
	common  /tetra_stat/    tetra_info,tetra_nindex
	common  /voids/			depth,isvoid
	
	do 100 idx = 1,ntetra
		depth(idx) = idx 
		! do 110 itrig = 1,4
		! 	jdx = tetra_neighbour(itrig,idx) 
		! 	if(jdx.eq.0)                      depth(idx) = 0
		! 	if(.not.btest(tetra_info(jdx),1)) depth(idx) = 0
110   continue
		if(.not.btest(tetra_info(idx),1)) depth(idx) = 0
		if(btest(tetra_info(idx),7))      depth(idx) = 1234567890
100	continue
	
c	write(*,*) "COMPUTING VOIDS"
	
c  loop over trig of each tetra and merge if trig not in alpha complex
	do 200 idx = 1,ntetra
		
c  skip dead tetra
		if(.not.btest(tetra_info(idx),1)) goto 200
		
c  if tetra in alpha complex, skip
		if(btest(tetra_info(idx),7)) goto 200
		
c  loop over trigs		
		do 300 itrig = 1,4

c  if trig is in alpha complex, skip it
			if(btest(tetra_info(idx),2+itrig)) goto 300

c  if neighbor tetra in alpha complex, no merge
			jdx = tetra_neighbour(itrig,idx)
			if(btest(tetra_info(jdx),7)) goto 300

c  merge tetra on either side of trig by taking lower depth
			! write(*,*) "merging"
			depth(idx) = min( depth(idx), depth(jdx) )
			depth(jdx) = min( depth(idx), depth(jdx) )			

300   continue
200   continue

		count = 0
	   do 400 idx = 1,ntetra
			isvoid(idx) = 0
			if(depth(idx).eq.0         ) goto 400
			if(depth(idx).eq.1234567890) goto 400
			! write(*,*) "FOUND VOID"
			isvoid(idx) = 1
			count = count + 1
400   continue

		! write(*,*) "found",count,"voids"

		do 500 idx = 1,ntetra
		write(*,*) idx,tetra(1,idx),tetra(2,idx),tetra(3,idx),
	1	 tetra(4,idx),btest(tetra_info(idx),1),btest(tetra_info(idx),7),
	2   depth(idx),isvoid(idx)
500   continue
		
	return
	end
