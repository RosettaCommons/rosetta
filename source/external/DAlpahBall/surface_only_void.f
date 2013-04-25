c	Surface_only.f		Version 1 9/1/2008	Patrice Koehl

c	This file contains a suite of routines for computed the (weighted) surface
c	area of a union of balls. The calculation allows also computing
c	derivatives

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

c	surface_only.f		Version 2 3/30/2007	Patrice Koehl

c	This subroutine builds the alpha complex based on the weighted
c	Delaunay triangulation

	subroutine surface_only_void(coef,WSurf,Surf,ballwsurf)

	integer	npointmax,ntetra_max,ntrig_max,nedge_max,nlink_max
	integer	ncortot

	parameter	(npointmax   = MAX_POINT)
	parameter	(ntetra_max  = MAX_TETRA)
	parameter	(nlink_max   = MAX_LINK)
	parameter	(ncortot     = 3*npointmax)

!	integer i,j,k,l,m,n
	integer	ia,ib,ic,id
	integer	i1,i2,i3,i4
	integer	ntetra,ntrig,nedge
	integer	npoints,nvertex,nred
	integer idx,iflag
	integer irad,iattach
	integer	itrig,jtrig,iedge
	integer	trig1,trig2,trig_in,trig_out,triga,trigb
	integer	ncheck,icheck,option,isave
	integer jtetra,ktetra,npass,ipair,i_out,sum

	integer*1 ival,it1,it2
	
	real*8  ballwsurfa,ballwsurfb,ballwsurfc,ballwsurfd

	integer	other3(3,4)
	integer face_info(2,6),face_pos(2,6)
	integer pair(2,6)
	integer	listcheck(nlink_max),tag(nlink_max)

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

	real*8	ang1,ang2,ang3,ang4,ang5,ang6
	real*8	ra,rb,rc,rd,ra2,rb2,rc2,rd2
	real*8	rab,rac,rad,rbc,rbd,rcd
	real*8	oab,oac,oad,obc,obd,ocd
	real*8	rab2,rac2,rad2,rbc2,rbd2,rcd2
	real*8	wa,wb,wc,wd,we,val
	real*8	alpha,scale,eps
	real*8	coefa,coefb,coefc,coefd,coefval
	real*8	surfa,surfb,surfc,surfd
	real*8	pi,twopi,precision
	real*8	a(3),b(3),c(3),d(3),e(3)
	real*8	radius(npointmax),weight(npointmax)
	real*8	coord(3*npointmax)

c	Results:

	real*8	WSurf,Surf
	real*8	coef(npointmax),ballwsurf(npointmax)

	integer   depth(ntetra_max)
	integer*1 isvoid(ntetra_max)

	common  /xyz_vertex/	 coord,radius,weight
	common  /tetra_zone/  ntetra,tetra,tetra_neighbour
	common  /tetra_stat/  tetra_info,tetra_nindex
	common  /alp_zone/	 tetra_edge
	common  /vertex_zone/ npoints,nvertex,vertex_info
	common  /gmp_info/	 scale,eps
	common  /workspace/	 tetra_mask
	common  /constants/	 pi,twopi,precision
	common  /voids/       depth,isvoid

	data other3 /2,3,4,1,3,4,1,2,4,1,2,3/
	data face_info/1,2,1,3,1,4,2,3,2,4,3,4/
	data face_pos/2,1,3,1,4,1,3,2,4,2,4,3/
	data pair/3,4,2,4,2,3,1,4,1,3,1,2/
	data isave /0/

	save
	do 10 i = 1,nvertex
		ballwsurf(i) = 0
10    continue

	do 100 idx = 1,ntetra
	
		if(isvoid(idx).eq.0) goto 100
		
		! write(*,*) "void calc",idx

		ia = tetra(1,idx)
		ib = tetra(2,idx)
		ic = tetra(3,idx)
		id = tetra(4,idx)

		ra = radius(ia)
		rb = radius(ib)
		rc = radius(ic)
		rd = radius(id)

		ra2 = ra*ra
		rb2 = rb*rb
		rc2 = rc*rc
		rd2 = rd*rd

		wa = 0.5*weight(ia)
		wb = 0.5*weight(ib)
		wc = 0.5*weight(ic)
		wd = 0.5*weight(id)

		do 200 m = 1,3
			a(m) = coord(3*(ia-1)+m)
			b(m) = coord(3*(ib-1)+m)
			c(m) = coord(3*(ic-1)+m)
			d(m) = coord(3*(id-1)+m)
200	continue

	    call distance2(coord,ia,ib,rab2,ncortot)
	    call distance2(coord,ia,ic,rac2,ncortot)
	    call distance2(coord,ia,id,rad2,ncortot)
	    call distance2(coord,ib,ic,rbc2,ncortot)
	    call distance2(coord,ib,id,rbd2,ncortot)
	    call distance2(coord,ic,id,rcd2,ncortot)

	    rab = sqrt(rab2)
	    rac = sqrt(rac2)
	    rad = sqrt(rad2)
	    rbc = sqrt(rbc2)
	    rbd = sqrt(rbd2)
	    rcd = sqrt(rcd2)

       oab = rab-ra-rb
       oac = rac-ra-rc
       obc = rbc-rb-rc
       oad = rad-ra-rd
       obd = rbd-rb-rd
       ocd = rcd-rc-rd

		call tetra_6dihed(a,b,c,d,ang1,ang2,ang3,ang4,ang5,ang6)
		write(*,*) "tetra",idx,"balls",ia,ib,ic,id
		write(*,*) "angles",ang1,ang2,ang3,ang4,ang5,ang6

c	First check all vertices of the tetrahedron

		ballwsurfa = 2*pi*radius(ia)*radius(ia)*(ang1+ang2+ang4-0.5)
		ballwsurfb = 2*pi*radius(ib)*radius(ib)*(ang1+ang3+ang5-0.5)
		ballwsurfc = 2*pi*radius(ic)*radius(ic)*(ang2+ang3+ang6-0.5)
		ballwsurfd = 2*pi*radius(id)*radius(id)*(ang4+ang5+ang6-0.5)

		write(*,*) "oneway"
		write(*,*) "ballsurf(ia)",ballwsurfa,(ang1+ang2+ang4-0.5)/2
		write(*,*) "ballsurf(ib)",ballwsurfb,(ang1+ang3+ang5-0.5)/2
		write(*,*) "ballsurf(ic)",ballwsurfc,(ang2+ang3+ang6-0.5)/2
		write(*,*) "ballsurf(id)",ballwsurfd,(ang4+ang5+ang6-0.5)/2

c	Now check all edges
c		iedge = 1		(c,d)
c		iedge = 2		(b,d)
c		iedge = 3		(b,c)
c		iedge = 4		(a,d)
c		iedge = 5		(a,c)
c		iedge = 6		(a,b)

		if(oab.lt.0) then
			call twosphere_surf(a,b,ra,ra2,rb,rb2,rab,rab2,surfa,surfb)
			ballwsurfa = ballwsurfa - ang1*surfa
			ballwsurfb = ballwsurfb - ang1*surfb
			write(*,*) "edge ab",ang1,surfa,surfb,rab-ra-rb
		endif
		if(oac.lt.0) then
			call twosphere_surf(a,c,ra,ra2,rc,rc2,rac,rac2,surfa,surfc)
			ballwsurfa = ballwsurfa - ang2*surfa
			ballwsurfc = ballwsurfc - ang2*surfc
			write(*,*) "edge ac",ang2,surfa,surfc,rac-ra-rc
		endif
		if(obc.lt.0) then
			call twosphere_surf(b,c,rb,rb2,rc,rc2,rbc,rbc2,surfb,surfc)
			ballwsurfb = ballwsurfb - ang3*surfb
			ballwsurfc = ballwsurfc - ang3*surfc
			write(*,*) "edge bc",ang3,surfb,surfc,rbc-rb-rc
		endif
		if(oad.lt.0) then
			call twosphere_surf(a,d,ra,ra2,rd,rd2,rad,rad2,surfa,surfd)
			ballwsurfa = ballwsurfa - ang4*surfa
			ballwsurfd = ballwsurfd - ang4*surfd
			write(*,*) "edge ad",ang4,surfa,surfd,rad-ra-rd
		endif
		if(obd.lt.0) then
			call twosphere_surf(b,d,rb,rb2,rd,rd2,rbd,rbd2,surfb,surfd)
			ballwsurfb = ballwsurfb - ang5*surfb
			ballwsurfd = ballwsurfd - ang5*surfd
			write(*,*) "edge bd",ang5,surfb,surfd,rbd-rb-rd
		endif
		if(ocd.lt.0) then
			call twosphere_surf(c,d,rc,rc2,rd,rd2,rcd,rcd2,surfc,surfd)
			ballwsurfc = ballwsurfc - ang6*surfc
			ballwsurfd = ballwsurfd - ang6*surfd
			write(*,*) "edge cd",ang6,surfc,surfd,rcd-rc-rd
		endif
	
		write(*,*) "twoway"
		write(*,*) "ballsurf(ia)",ballwsurfa
		write(*,*) "ballsurf(ib)",ballwsurfb
		write(*,*) "ballsurf(ic)",ballwsurfc
		write(*,*) "ballsurf(id)",ballwsurfd

c	Finally check faces

		if(btest(tetra_info(idx),3)) then
			surfa = 0
			call threesphere_surf(b,c,d,rb,
	1					rc,rd,rb2,rc2,rd2,wb,wc,wd,rbc,
	2					rbd,rcd,rbc2,rbd2,rcd2,surfb,
	3					surfc,surfd)
			ballwsurfb = ballwsurfb+0.5*surfb
			ballwsurfc = ballwsurfc+0.5*surfc
			ballwsurfd = ballwsurfd+0.5*surfd
			write(*,*) "trig1",surfa,surfb,surfc,surfd
		endif
		if(btest(tetra_info(idx),4)) then
			surfb = 0
			call threesphere_surf(a,c,d,ra,
	1					rc,rd,ra2,rc2,rd2,wa,wc,wd,rac,
	2					rad,rcd,rac2,rad2,rcd2,surfa,
	3					surfc,surfd)
			ballwsurfa = ballwsurfa+0.5*surfa
			ballwsurfc = ballwsurfc+0.5*surfc
			ballwsurfd = ballwsurfd+0.5*surfd						
			write(*,*) "trig2",surfa,surfb,surfc,surfd
		endif
		if(btest(tetra_info(idx),5)) then
			surfc = 0
			call threesphere_surf(a,b,d,ra,
	1					rb,rd,ra2,rb2,rd2,wa,wb,wd,rab,
	2					rad,rbd,rab2,rad2,rbd2,surfa,
	3					surfb,surfd)
			ballwsurfa = ballwsurfa+0.5*surfa
			ballwsurfb = ballwsurfb+0.5*surfb
			ballwsurfd = ballwsurfd+0.5*surfd						
			write(*,*) "trig3",surfa,surfb,surfc,surfd
		endif
		if(btest(tetra_info(idx),6)) then
			surfd = 0
			call threesphere_surf(a,b,c,ra,
	1					rb,rc,ra2,rb2,rc2,wa,wb,wc,rab,
	2					rac,rbc,rab2,rac2,rbc2,surfa,
	3					surfb,surfc)
			ballwsurfa = ballwsurfa+0.5*surfa						
			ballwsurfb = ballwsurfb+0.5*surfb
			ballwsurfc = ballwsurfc+0.5*surfc
			write(*,*) "trig4",surfa,surfb,surfc,surfd
		endif

		write(*,*) "threeway"
		write(*,*) "ballsurf(ia)",ballwsurfa
		write(*,*) "ballsurf(ib)",ballwsurfb
		write(*,*) "ballsurf(ic)",ballwsurfc
		write(*,*) "ballsurf(id)",ballwsurfd

		if(ballwsurfa.gt.0) ballwsurf(ia) = ballwsurf(ia) + ballwsurfa
		if(ballwsurfb.gt.0) ballwsurf(ib) = ballwsurf(ib) + ballwsurfb
		if(ballwsurfc.gt.0) ballwsurf(ic) = ballwsurf(ic) + ballwsurfc
		if(ballwsurfd.gt.0) ballwsurf(id) = ballwsurf(id) + ballwsurfd
		
		write(*,*)

100	   continue

	WSurf = 0
	Surf = 0
	do 1800 i = 1,nvertex
		WSurf = Wsurf + coef(i)*ballwsurf(i)
		Surf  = Surf  + ballwsurf(i)
		ballwsurf(i) = ballwsurf(i)*coef(i)
1800	continue

	return
	end
