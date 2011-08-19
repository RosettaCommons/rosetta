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

c	Alfcx.f		Version 2 3/30/2007	Patrice Koehl

c	This subroutine builds the alpha complex based on the weighted
c	Delaunay triangulation

	subroutine alfcx(alpha,nred,listred)

	integer	npointmax,ntetra_max,nlink_max,nredmax

	parameter	(npointmax   = MAX_POINT)
	parameter	(ntetra_max  = MAX_TETRA)
	parameter	(nlink_max   = MAX_LINK)
	parameter	(nredmax     = MAX_RED)

	integer i,j,k,l,m,n
	integer	ia,ib,i1,i2
	integer	ntetra,ntrig,nedge
	integer	npoints,nvertex,nred
	integer idx,iflag
	integer irad,iattach
	integer	itrig,jtrig,iedge
	integer	trig1,trig2,trig_in,trig_out,triga,trigb
	integer	ncheck,icheck,jtetra,itetra,ktetra
	integer	npass,ipair
	integer	i_out

	integer*1 ival

	integer	other3(3,4)
	integer face_info(2,6),face_pos(2,6)
	integer pair(2,6)
	integer	listred(nredmax)
	integer	listcheck(nlink_max)

	integer*1 tetra_info(ntetra_max)
	integer*1 tetra_nindex(ntetra_max)
	integer*1 tetra_mask(ntetra_max)
	integer*1 tetra_edge(ntetra_max)

c	Information on the tetrahedra of the regular
c	triangulation

	integer tetra(4,ntetra_max)
	integer tetra_neighbour(4,ntetra_max)

c	Information on the vertices

	integer*1 vertex_info(npointmax)

	real*8	wa,wb,wc,wd,we
	real*8	alpha,scale,eps
	real*8	a(3),b(3),c(3),d(3),e(3)
	real*8	radius(npointmax),weight(npointmax)
	real*8	coord(3*npointmax)

	logical testa,testb

	common  /xyz_vertex/	coord,radius,weight
	common  /tetra_zone/    ntetra,tetra,tetra_neighbour
	common  /tetra_stat/    tetra_info,tetra_nindex
	common  /alp_zone/	tetra_edge
	common  /vertex_zone/   npoints,nvertex,vertex_info
	common  /gmp_info/	scale,eps
	common  /workspace/	tetra_mask

	data other3 /2,3,4,1,3,4,1,2,4,1,2,3/
	data face_info/1,2,1,3,1,4,2,3,2,4,3,4/
	data face_pos/2,1,3,1,4,1,3,2,4,2,4,3/
	data pair/3,4,2,4,2,3,1,4,1,3,1,2/

	save

	ival = 0
	do 100 i = 1,ntetra
		tetra_mask(i) = 0
		call mvbits(ival,0,5,tetra_info(i),3)
100	continue

	call set_alf_gmp

c	Loop over all tetrahedra:

	do 300 idx = 1,ntetra

c		"Dead" tetrahedron are ignored

		if(.not.btest(tetra_info(idx),1)) goto 300

		i = tetra(1,idx)
		j = tetra(2,idx)
		k = tetra(3,idx)
		l = tetra(4,idx)

		call alf_tetra(i,j,k,l,iflag,eps,scale,alpha)

		if(iflag.eq.1) then
			tetra_info(idx)=ibset(tetra_info(idx),7)
		endif

300	continue

c	Now loop over all triangles: each triangle is defined implicitly
c	as the interface between two tetrahedra i and j with i < j

	do 500 idx = 1,ntetra

c		"Dead" tetrahedron are ignored

	     if(.not.btest(tetra_info(idx),1)) goto 500

	     do 400 itrig = 1,4

		jtetra = tetra_neighbour(itrig,idx)
		ival = ibits(tetra_nindex(idx),2*(itrig-1),2)
		jtrig = ival + 1

		if(jtetra.eq.0.or.jtetra.gt.idx) then

c			We are checking the triangle defined
c			by itetra and jtetra
c			If one of them belongs to the alpha complex,
c			the triangle belongs to the alpha complex

			if(btest(tetra_info(idx),7)) then
				tetra_info(idx) = ibset(
     1				tetra_info(idx),2+itrig)
				if(jtetra.ne.0) then
					tetra_info(jtetra) = ibset(
     1					tetra_info(jtetra),2+jtrig)
				endif
				goto 400
			endif

			if(jtetra.ne.0) then
				if(btest(tetra_info(jtetra),7)) then
					tetra_info(idx) = ibset(
     1					tetra_info(idx),2+itrig)
					tetra_info(jtetra) = ibset(
     1					tetra_info(jtetra),2+jtrig)
					goto 400
				endif
			endif

c			If we are here, it means that the two
c			attached tetrahedra do not belong to the
c			alpha complex: need to check the triangle
c			itself

c			Define the 3 vertices of the triangle, as well as the 2
c			remaining vertices of the two tetrahedra attached to the
c			triangle

			i = tetra(other3(1,itrig),idx)
			j = tetra(other3(2,itrig),idx)
			k = tetra(other3(3,itrig),idx)

			l = tetra(itrig,idx)
			if(jtetra.ne.0) then
				m = tetra(jtrig,jtetra)
			else
				m = 0
			endif

			call alf_trig(i,j,k,l,m,irad,iattach,eps,
     1				scale,alpha)

			if(iattach.eq.0.and.irad.eq.1) then
				l = 1
				tetra_info(idx) = 
     1				ibset(tetra_info(idx),2+itrig)
				if(jtetra.ne.0) then
				     tetra_info(jtetra) = 
     1				     ibset(tetra_info(jtetra),2+jtrig)
				endif
			endif

		endif

400	    continue

500	continue

c	Now loop over all edges: each edge is defined implicitly
c	by the tetrahedra to which it belongs

	do 600 idx = 1,ntetra
		tetra_mask(idx) = 0
		tetra_edge(idx) = 0
600	continue

	do 1100 itetra = 1,ntetra

		if(.not.btest(tetra_info(itetra),1)) goto 1100

	    do 1000 iedge = 1,6

		if(btest(tetra_mask(itetra),iedge-1)) goto 1000

c		For each edge, check triangles attached to the edge
c		if at least one of these triangles is in alpha complex,
c		then the edge is in the alpha complex
c		We then put the two vertices directly in the alpha complex
c		Otherwise, build list of triangles to check

c		itetra is one tetrahedron (a,b,c,d) containing the edge

c		iedge is the edge number in the tetrahedron itetra, with:
c		iedge = 1		(c,d)
c		iedge = 2		(b,d)
c		iedge = 3		(b,c)
c		iedge = 4		(a,d)
c		iedge = 5		(a,c)
c		iedge = 6		(a,b)
c		
c		Define indices of the edge

		i = tetra(pair(1,iedge),itetra)
		j = tetra(pair(2,iedge),itetra)

c		trig1 and trig2 are the two faces of itetra that share
c		iedge
c		i1 and i2 are the positions of the third vertices of
c		trig1 and trig2

		trig1 = face_info(1,iedge)
		i1    = face_pos(1,iedge)
		trig2 = face_info(2,iedge)
		i2    = face_pos(2,iedge)

		ia = tetra(i1,itetra)
		ib = tetra(i2,itetra)

		icheck = 0
		if(btest(tetra_info(itetra),2+trig1)) then
			tetra_edge(itetra)=ibset(tetra_edge(itetra),
     1			iedge-1)
			vertex_info(i) = ibset(vertex_info(i),7)
			vertex_info(j) = ibset(vertex_info(j),7)
			goto 1000
		else
			icheck = icheck + 1
			listcheck(icheck) = ia
		endif
		if(btest(tetra_info(itetra),2+trig2)) then
			tetra_edge(itetra)=ibset(tetra_edge(itetra),
     1			iedge-1)
			vertex_info(i) = ibset(vertex_info(i),7)
			vertex_info(j) = ibset(vertex_info(j),7)
			goto 1000
		else
			icheck = icheck + 1
			listcheck(icheck) = ib
		endif

c		Now we look at the star of the edge:

		ktetra = itetra
		npass = 1
		trig_out = trig1
		jtetra = tetra_neighbour(trig_out,ktetra)

700		continue

c		Leave this side of the star if we hit the convex hull

		if(jtetra.eq.0) goto 800

c		Leave the loop completely if we have described the
c		full cycle

		if(jtetra.eq.itetra) goto 900

c		Identify the position of iedge in tetrahedron jtetra

		if(i.eq.tetra(1,jtetra)) then
			if(j.eq.tetra(2,jtetra)) then
				ipair = 6
			elseif(j.eq.tetra(3,jtetra)) then
				ipair = 5
			else
				ipair = 4
			endif
		elseif(i.eq.tetra(2,jtetra)) then
			if(j.eq.tetra(3,jtetra)) then
				ipair = 3
			else
				ipair = 2
			endif
		else
			ipair = 1
		endif

		tetra_mask(jtetra) = ibset(tetra_mask(jtetra),ipair-1)

c		Find out the face we "went in"

		ival = ibits(tetra_nindex(ktetra),2*(trig_out-1),2)
		trig_in = ival + 1

c		We know the two faces of jtetra that share iedge:

		triga = face_info(1,ipair)
		i1    = face_pos(1,ipair)
		trigb = face_info(2,ipair)
		i2    = face_pos(2,ipair)

		trig_out = triga
		i_out = i1
		if(trig_in.eq.triga) then
			i_out = i2
			trig_out = trigb
		endif

c		Now check if trig_out is already in the alpha complex.
c		if it is, iedge is in; otherwise, will need an attach 
c	 	test...

		if(btest(tetra_info(jtetra),2+trig_out)) then
			tetra_edge(itetra)=ibset(tetra_edge(itetra),
     1			iedge-1)
			vertex_info(i) = ibset(vertex_info(i),7)
			vertex_info(j) = ibset(vertex_info(j),7)
			goto 1000
		endif

		ktetra = jtetra
		jtetra = tetra_neighbour(trig_out,ktetra)
		if(jtetra.eq.itetra) goto 900

		icheck = icheck + 1
		listcheck(icheck) = tetra(i_out,ktetra)

		goto 700

800		continue

		if(npass.eq.2) goto 900
		npass = npass + 1
		ktetra = itetra
		trig_out = trig2
		jtetra = tetra_neighbour(trig_out,ktetra)
		goto 700

900		continue

c		If we got here, it means that none of the triangles
c		in the star of the edge belongs to the alpha complex:
c		this is a singular edge.
c		We check if the edge is attached, and if alpha is
c		smaller than the radius of the sphere orthogonal
c		to the two balls corresponding to the edge

		call alf_edge(i,j,icheck,listcheck,irad,iattach,
     1		eps,scale,alpha)

		if(iattach.eq.0.and.irad.eq.1) then
			tetra_edge(itetra)=ibset(tetra_edge(itetra),
     1			iedge-1)
			vertex_info(i) = ibset(vertex_info(i),7)
			vertex_info(j) = ibset(vertex_info(j),7)
			goto 1000
		endif

c		Edge is not in alpha complex: now check if the two vertices
c		could be attached to each other: 

		call vertex_attach(i,j,testa,testb,eps)

		if(testa) vertex_info(i) = ibset(vertex_info(i),6)
		if(testb) vertex_info(j) = ibset(vertex_info(j),6)

1000	   continue

1100	continue

c	Now loop over vertices

	nred = 0
	do 1200 i = 1,nvertex

		if(.not.btest(vertex_info(i),0)) goto 1200

c		if vertex was already set in alpha-complex, nothing to do...

		if(btest(vertex_info(i),7)) goto 1150

c		Vertex is in alpha complex...unless it is attached!

		if(.not.btest(vertex_info(i),6)) then
			vertex_info(i) = ibset(vertex_info(i),7)
		else
			nred = nred + 1
			listred(nred) = i
		endif

1150		continue

1200	continue

c	We are done!

	return
	end
