c	Volume_deriv_dist.f		Version 1 9/1/2008	Patrice Koehl
c
c	This file contains a suite of routines for computed the (weighted) surface
c	area and volume of a union of balls. The calculation allows also computing
c	derivatives with respect to internal distances
c
c	Copyright (C) 2005 Patrice Koehl
c
c	This library is free software; you can redistribute it and/or
c	modify it under the terms of the GNU Lesser General Public
c	License as published by the Free Software Foundation; either
c	version 2.1 of the License, or (at your option) any later version.
c
c	This library is distributed in the hope that it will be useful,
c	but WITHOUT ANY WARRANTY; without even the implied warranty of
c	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
c	Lesser General Public License for more details.
c
c	You should have received a copy of the GNU Lesser General Public
c	License along with this library; if not, write to the Free Software
c	Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
c
#include "defines.h"
c
c	volume_deriv_dist.f		Version 1 9/1/2008	Patrice Koehl
c
c	This subroutine builds the alpha complex based on the weighted
c	Delaunay triangulation
c
	subroutine volume_deriv_dist(coef,WSurf,WVol,Surf,Vol,
     1          ballwsurf,ballwvol,dwsurf,dwvol,dwsurf_dist,dwvol_dist)
c
	integer	npointmax,ntetra_max,ntrig_max,nedge_max,nlink_max
	integer	ncortot,npairmax
c
	parameter	(npointmax   = MAX_POINT)
	parameter	(ntetra_max  = MAX_TETRA)
	parameter	(nlink_max   = MAX_LINK)
	parameter	(ncortot     = 3*npointmax)
        parameter       (npairmax    = MAX_PAIR)
c
	integer i,j,k,l,m,n
        integer iab,iac,iad,ibc,ibd,icd
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
c
        integer ipos
c
	integer*1 ival,it1,it2
c
	integer	other3(3,4)
	integer face_info(2,6),face_pos(2,6)
	integer pair(2,6)
	integer	listcheck(nlink_max),tag(nlink_max)
c
	integer*1 tetra_info(ntetra_max)
	integer*1 tetra_nindex(ntetra_max)
	integer*1 tetra_mask(ntetra_max)
	integer*1 tetra_edge(ntetra_max)
c
c	Information on the tetrahedra of the regular
c	triangulation
c
	integer tetra(4,ntetra_max)
	integer tetra_neighbour(4,ntetra_max)
c
c	Information on the vertices
c
	integer*1 vertex_info(npointmax)
c
	real*8	ra,rb,rc,rd,ra2,rb2,rc2,rd2
	real*8	rab,rac,rad,rbc,rbd,rcd
	real*8	rab2,rac2,rad2,rbc2,rbd2,rcd2
	real*8	wa,wb,wc,wd,we,val
	real*8	alpha,scale,eps
	real*8	coefa,coefb,coefc,coefd,coefval
	real*8	surfa,surfb,surfc,surfd
	real*8	vola,volb,volc,vold
	real*8	pi,twopi,precision
        real*8  sh_abc,sh_acb,sh_bca,sh_abd,sh_adb,sh_bda
        real*8  sh_acd,sh_adc,sh_cda,sh_bcd,sh_bdc,sh_cdb
        real*8  eps_abc,eps_abd,eps_acd,eps_bcd
c
	real*8	radius(npointmax),weight(npointmax)
	real*8	coord(3*npointmax)
	real*8	a(3),b(3),c(3),d(3),e(3),u_ab(3)
        real*8  pabc(3),pacb(3),pabd(3),padb(3),pacd(3),padc(3)
        real*8  pbcd(3),pbdc(3)
	real*8	dsurfa(6),dsurfb(6),dsurfc(6),dsurfd(6)
	real*8	dsurfa2,dsurfb2,dsurfc2,dsurfd2
	real*8	dsurfa3(3),dsurfb3(3),dsurfc3(3),dsurfd3(3)
	real*8	dvola(6),dvolb(6),dvolc(6),dvold(6)
	real*8	dvola2,dvolb2,dvolc2,dvold2
	real*8	dvola3(3),dvolb3(3),dvolc3(3),dvold3(3)
        real*8  der_shabc(3),der_shacb(3),der_shbca(3),der_shabd(3)
        real*8  der_shadb(3),der_shbda(3),der_shacd(3),der_shadc(3)
        real*8  der_shcda(3),der_shbcd(3),der_shbdc(3),der_shcdb(3)
        real*8  der_abc(6),der_acb(6),der_bca(6),der_abd(6),der_adb(6)
        real*8  der_bda(6),der_acd(6),der_adc(6),der_cda(6),der_bcd(6)
        real*8  der_bdc(6),der_cdb(6)
c
c	Results:
c
	real*8	WSurf,WVol,Surf,Vol
	real*8	coef(npointmax),ballwsurf(npointmax)
        real*8  ballwvol(npointmax)
	real*8	dwsurf(3,npointmax)
	real*8	dwvol(3,npointmax)
        real*8  dwsurf_dist(npairmax),dwvol_dist(npairmax)
c
	common  /xyz_vertex/	coord,radius,weight
	common  /tetra_zone/    ntetra,tetra,tetra_neighbour
	common  /tetra_stat/    tetra_info,tetra_nindex
	common  /alp_zone/	tetra_edge
	common  /vertex_zone/   npoints,nvertex,vertex_info
	common  /gmp_info/	scale,eps
	common  /workspace/	tetra_mask
	common  /constants/	pi,twopi,precision
c
	data other3 /2,3,4,1,3,4,1,2,4,1,2,3/
	data face_info/1,2,1,3,1,4,2,3,2,4,3,4/
	data face_pos/2,1,3,1,4,1,3,2,4,2,4,3/
	data pair/3,4,2,4,2,3,1,4,1,3,1,2/
	data isave /0/
c
	save
c
	if(isave.eq.0) then
		pi = acos(-1.d0)
		twopi = 2.d0 * pi
		isave = 1
		precision = 1.0d-10
	endif
c
        do 50 i = 1,(nvertex*(nvertex-1))/2
                dwsurf_dist(i) = 0
                dwvol_dist(i) = 0
50      continue
c
	Wsurf = 0
	Surf  = 0
c
	do 200 i = 1,nvertex
		ballwsurf(i) = 0
		ballwvol(i) = 0
		do 100 j = 1,3
			dwsurf(j,i) = 0
			dwvol(j,i) = 0
100		continue
200	continue
c
	do 300 i = 1,ntetra
		tetra_mask(i) = 0
300	continue
c
c	loop over tetra, triangles and edges
c
	do 1600 idx = 1,ntetra
c
	    if(.not.btest(tetra_info(idx),1)) goto 1600
c
	    ia = tetra(1,idx)
	    ib = tetra(2,idx)
	    ic = tetra(3,idx)
	    id = tetra(4,idx)
c
            iab = ipos(ia,ib)
            iac = ipos(ia,ic)
            iad = ipos(ia,id)
            ibc = ipos(ib,ic)
            ibd = ipos(ib,id)
            icd = ipos(ic,id)
c
	    do 400 i = 1,3
		a(i) = coord(3*(ia-1)+i)	   
		b(i) = coord(3*(ib-1)+i)	   
		c(i) = coord(3*(ic-1)+i)	   
		d(i) = coord(3*(id-1)+i)	   
400	    continue
c
	    ra = radius(ia)
	    rb = radius(ib)
	    rc = radius(ic)
	    rd = radius(id)
c
	    ra2 = ra*ra
	    rb2 = rb*rb
	    rc2 = rc*rc
	    rd2 = rd*rd
c
	    call distance2(coord,ia,ib,rab2,ncortot)
	    call distance2(coord,ia,ic,rac2,ncortot)
	    call distance2(coord,ia,id,rad2,ncortot)
	    call distance2(coord,ib,ic,rbc2,ncortot)
	    call distance2(coord,ib,id,rbd2,ncortot)
	    call distance2(coord,ic,id,rcd2,ncortot)
c
	    rab = sqrt(rab2)
	    rac = sqrt(rac2)
	    rad = sqrt(rad2)
	    rbc = sqrt(rbc2)
	    rbd = sqrt(rbd2)
	    rcd = sqrt(rcd2)
c
	    wa = 0.5d0*weight(ia)
	    wb = 0.5d0*weight(ib)
	    wc = 0.5d0*weight(ic)
	    wd = 0.5d0*weight(id)
c
	    coefa = coef(ia)
	    coefb = coef(ib)
	    coefc = coef(ic)
	    coefd = coef(id)
c
	    if(btest(tetra_info(idx),7)) then
c
                call threesphere_info(a,b,c,ra,rb,rc,ra2,rb2,rc2,
     1          wa,wb,wc,rab,rac,rbc,rab2,rac2,rbc2,pabc,pacb,eps_abc,
     2          sh_abc,sh_acb,sh_bca,der_shabc,der_shacb,der_shbca)
c
                call threesphere_info(a,b,d,ra,rb,rd,ra2,rb2,rd2,
     1          wa,wb,wd,rab,rad,rbd,rab2,rad2,rbd2,pabd,padb,eps_abd,
     2          sh_abd,sh_adb,sh_bda,der_shabd,der_shadb,der_shbda)
c
                call threesphere_info(a,c,d,ra,rc,rd,ra2,rc2,rd2,
     1          wa,wc,wd,rac,rad,rcd,rac2,rad2,rcd2,pacd,padc,eps_acd,
     2          sh_acd,sh_adc,sh_cda,der_shacd,der_shadc,der_shcda)
c
                call threesphere_info(b,c,d,rb,rc,rd,rb2,rc2,rd2,
     1          wb,wc,wd,rbc,rbd,rcd,rbc2,rbd2,rcd2,pbcd,pbdc,eps_bcd,
     2          sh_bcd,sh_bdc,sh_cdb,der_shbcd,der_shbdc,der_shcdb)
c
	    	if(btest(tetra_info(idx),0)) then
c
                        call make_shder(1,der_shabc,der_shacb,
     1                  der_shbca,der_shabd,der_shadb,der_shbda,
     2                  der_shacd,der_shadc,der_shcda,der_shbcd,
     3                  der_shbdc,der_shcdb,
     6                  der_abc,der_acb,der_bca,der_abd,
     7                  der_adb,der_bda,der_acd,der_adc,
     7                  der_cda,der_bcd,der_bdc,der_cdb)
c
			call foursphere_dvol_dist(a,b,c,d,ra,rb,rc,rd,
     1			ra2,rb2,rc2,rd2,rab,rac,rad,rbc,rbd,rcd,
     2			rab2,rac2,rad2,rbc2,rbd2,rcd2,wa,wb,wc,wd,
     3                  eps_abc,eps_abd,eps_acd,eps_bcd,
     4                  sh_abc,sh_acb,sh_bca,sh_abd,sh_adb,sh_bda,
     5                  sh_acd,sh_adc,sh_cda,sh_bcd,sh_bdc,sh_cdb,
     6                  der_abc,der_acb,der_bca,der_abd,
     7                  der_adb,der_bda,der_acd,der_adc,
     7                  der_cda,der_bcd,der_bdc,der_cdb,
     8                  pacb,pabd,padc,pbcd,
     9			surfa,surfb,surfc,surfd,vola,volb,volc,vold,
     *                  dsurfa,dsurfb,dsurfc,dsurfd,dvola,dvolb,dvolc,
     *                  dvold)
c
		else
c
                        call make_shder(2,der_shabc,der_shacb,
     1                  der_shbca,der_shabd,der_shadb,der_shbda,
     2                  der_shacd,der_shadc,der_shcda,der_shbcd,
     3                  der_shbdc,der_shcdb,
     4                  der_abc,der_acb,der_bca,der_abd,
     5                  der_adb,der_bda,der_acd,der_adc,
     6                  der_cda,der_bcd,der_bdc,der_cdb)
c
			call foursphere_dvol_dist(a,b,d,c,ra,rb,rd,rc,
     1			ra2,rb2,rd2,rc2,rab,rad,rac,rbd,rbc,rcd,
     2			rab2,rad2,rac2,rbd2,rbc2,rcd2,wa,wb,wd,wc,
     3                  eps_abd,eps_abc,eps_acd,eps_bcd,
     4                  sh_abd,sh_adb,sh_bda,sh_abc,sh_acb,sh_bca,
     5                  sh_adc,sh_acd,sh_cda,sh_bdc,sh_bcd,sh_cdb,
     6                  der_abd,der_adb,der_bda,der_abc,
     7                  der_acb,der_bca,der_adc,der_acd,
     7                  der_cda,der_bdc,der_bcd,der_cdb,
     8                  padb,pabc,pacd,pbdc,
     9			surfa,surfb,surfd,surfc,vola,volb,vold,volc,
     *                  dsurfa,dsurfb,dsurfd,dsurfc,dvola,dvolb,dvold,
     *                  dvolc)
c
                        call permute_ccw(dsurfa)
                        call permute_ccw(dsurfb)
                        call permute_ccw(dsurfc)
                        call permute_ccw(dsurfd)
                        call permute_ccw(dvola)
                        call permute_ccw(dvolb)
                        call permute_ccw(dvolc)
                        call permute_ccw(dvold)
		endif
c
		ballwsurf(ia) = ballwsurf(ia) - surfa
		ballwsurf(ib) = ballwsurf(ib) - surfb
		ballwsurf(ic) = ballwsurf(ic) - surfc
		ballwsurf(id) = ballwsurf(id) - surfd
c
		ballwvol(ia) = ballwvol(ia) - vola
		ballwvol(ib) = ballwvol(ib) - volb
		ballwvol(ic) = ballwvol(ic) - volc
		ballwvol(id) = ballwvol(id) - vold
c
                dwsurf_dist(iab) = dwsurf_dist(iab) - coefa*dsurfa(1)
     1          -coefb*dsurfb(1) - coefc*dsurfc(1)-coefd*dsurfd(1)
                dwsurf_dist(iac) = dwsurf_dist(iac) - coefa*dsurfa(2)
     1          -coefb*dsurfb(2) - coefc*dsurfc(2)-coefd*dsurfd(2)
                dwsurf_dist(iad) = dwsurf_dist(iad) - coefa*dsurfa(3)
     1          -coefb*dsurfb(3) - coefc*dsurfc(3)-coefd*dsurfd(3)
                dwsurf_dist(ibc) = dwsurf_dist(ibc) - coefa*dsurfa(4)
     1          -coefb*dsurfb(4) - coefc*dsurfc(4)-coefd*dsurfd(4)
                dwsurf_dist(ibd) = dwsurf_dist(ibd) - coefa*dsurfa(5)
     1          -coefb*dsurfb(5) - coefc*dsurfc(5)-coefd*dsurfd(5)
                dwsurf_dist(icd) = dwsurf_dist(icd) - coefa*dsurfa(6)
     1          -coefb*dsurfb(6) - coefc*dsurfc(6)-coefd*dsurfd(6)
c
                dwvol_dist(iab) = dwvol_dist(iab) - coefa*dvola(1)
     1          -coefb*dvolb(1) - coefc*dvolc(1)-coefd*dvold(1)
                dwvol_dist(iac) = dwvol_dist(iac) - coefa*dvola(2)
     1          -coefb*dvolb(2) - coefc*dvolc(2)-coefd*dvold(2)
                dwvol_dist(iad) = dwvol_dist(iad) - coefa*dvola(3)
     1          -coefb*dvolb(3) - coefc*dvolc(3)-coefd*dvold(3)
                dwvol_dist(ibc) = dwvol_dist(ibc) - coefa*dvola(4)
     1          -coefb*dvolb(4) - coefc*dvolc(4)-coefd*dvold(4)
                dwvol_dist(ibd) = dwvol_dist(ibd) - coefa*dvola(5)
     1          -coefb*dvolb(5) - coefc*dvolc(5)-coefd*dvold(5)
                dwvol_dist(icd) = dwvol_dist(icd) - coefa*dvola(6)
     1          -coefb*dvolb(6) - coefc*dvolc(6)-coefd*dvold(6)
c
	    endif
c
c	    Check triangles
c
	    do 700 itrig = 1,4
c
		jtetra = tetra_neighbour(itrig,idx)
		ival = ibits(tetra_nindex(idx),2*(itrig-1),2)
		jtrig = ival + 1
c
		if(jtetra.eq.0.or.jtetra.gt.idx) then
c
			if(btest(tetra_info(idx),2+itrig)) then
				call mvbits(tetra_info(idx),7,1,it1,0)
				if(jtetra.ne.0) then
			    		call mvbits(tetra_info(jtetra),
     1					7,1,it2,0)
				else
					it2 = 0
				endif
				coefval = 2 - it1 - it2
				if(coefval.eq.0) goto 700
c
				coefval = 0.5d0*coefval
c
				if(itrig.eq.1) then
                                     surfa = 0
                                     vola = 0
				     call threesphere_dvol_dist(b,c,
     1				     d,rb,rc,rd,rb2,rc2,rd2,wb,wc,wd,
     2				     rbc,rbd,rcd,rbc2,rbd2,rcd2,
     3                               surfb,surfc,surfd,volb,volc,vold,
     4                               dsurfb3,dsurfc3,dsurfd3,dvolb3,
     5                               dvolc3,dvold3)
                                     dwsurf_dist(ibc)=dwsurf_dist(ibc)+
     1                               coefval*(coefb*dsurfb3(1)+coefc*
     2                               dsurfc3(1)+coefd*dsurfd3(1))
                                     dwsurf_dist(ibd)=dwsurf_dist(ibd)+
     1                               coefval*(coefb*dsurfb3(2)+coefc*
     2                               dsurfc3(2)+coefd*dsurfd3(2))
                                     dwsurf_dist(icd)=dwsurf_dist(icd)+
     1                               coefval*(coefb*dsurfb3(3)+coefc*
     2                               dsurfc3(3)+coefd*dsurfd3(3))
                                     dwvol_dist(ibc)=dwvol_dist(ibc)+
     1                               coefval*(coefb*dvolb3(1)+coefc*
     2                               dvolc3(1)+coefd*dvold3(1))
                                     dwvol_dist(ibd)=dwvol_dist(ibd)+
     1                               coefval*(coefb*dvolb3(2)+coefc*
     2                               dvolc3(2)+coefd*dvold3(2))
                                     dwvol_dist(icd)=dwvol_dist(icd)+
     1                               coefval*(coefb*dvolb3(3)+coefc*
     2                               dvolc3(3)+coefd*dvold3(3))
				elseif(itrig.eq.2) then
				     surfb = 0
                                     volb = 0
				     call threesphere_dvol_dist(a,c,
     1				     d,ra,rc,rd,ra2,rc2,rd2,wa,wc,wd,
     2				     rac,rad,rcd,rac2,rad2,rcd2,
     3                               surfa,surfc,surfd,vola,volc,vold,
     4                               dsurfa3,dsurfc3,dsurfd3,dvola3,
     5                               dvolc3,dvold3)
                                     dwsurf_dist(iac)=dwsurf_dist(iac)+
     1                               coefval*(coefa*dsurfa3(1)+coefc*
     2                               dsurfc3(1)+coefd*dsurfd3(1))
                                     dwsurf_dist(iad)=dwsurf_dist(iad)+
     1                               coefval*(coefa*dsurfa3(2)+coefc*
     2                               dsurfc3(2)+coefd*dsurfd3(2))
                                     dwsurf_dist(icd)=dwsurf_dist(icd)+
     1                               coefval*(coefa*dsurfa3(3)+coefc*
     2                               dsurfc3(3)+coefd*dsurfd3(3))
                                     dwvol_dist(iac)=dwvol_dist(iac)+
     1                               coefval*(coefa*dvola3(1)+coefc*
     2                               dvolc3(1)+coefd*dvold3(1))
                                     dwvol_dist(iad)=dwvol_dist(iad)+
     1                               coefval*(coefa*dvola3(2)+coefc*
     2                               dvolc3(2)+coefd*dvold3(2))
                                     dwvol_dist(icd)=dwvol_dist(icd)+
     1                               coefval*(coefa*dvola3(3)+coefc*
     2                               dvolc3(3)+coefd*dvold3(3))
				elseif(itrig.eq.3) then
				     surfc = 0
                                     volc  = 0
				     call threesphere_dvol_dist(a,b,
     1				     d,ra,rb,rd,ra2,rb2,rd2,wa,wb,wd,
     2				     rab,rad,rbd,rab2,rad2,rbd2,
     3                               surfa,surfb,surfd,vola,volb,vold,
     4                               dsurfa3,dsurfb3,dsurfd3,dvola3,
     5                               dvolb3,dvold3)
                                     dwsurf_dist(iab)=dwsurf_dist(iab)+
     1                               coefval*(coefa*dsurfa3(1)+coefb*
     2                               dsurfb3(1)+coefd*dsurfd3(1))
                                     dwsurf_dist(iad)=dwsurf_dist(iad)+
     1                               coefval*(coefa*dsurfa3(2)+coefb*
     2                               dsurfb3(2)+coefd*dsurfd3(2))
                                     dwsurf_dist(ibd)=dwsurf_dist(ibd)+
     1                               coefval*(coefa*dsurfa3(3)+coefb*
     2                               dsurfb3(3)+coefd*dsurfd3(3))
                                     dwvol_dist(iab)=dwvol_dist(iab)+
     1                               coefval*(coefa*dvola3(1)+coefb*
     2                               dvolb3(1)+coefd*dvold3(1))
                                     dwvol_dist(iad)=dwvol_dist(iad)+
     1                               coefval*(coefa*dvola3(2)+coefb*
     2                               dvolb3(2)+coefd*dvold3(2))
                                     dwvol_dist(ibd)=dwvol_dist(ibd)+
     1                               coefval*(coefa*dvola3(3)+coefb*
     2                               dvolb3(3)+coefd*dvold3(3))
				elseif(itrig.eq.4) then
				     surfd = 0
                                     vold  = 0
				     call threesphere_dvol_dist(a,b,
     1				     c,ra,rb,rc,ra2,rb2,rc2,wa,wb,wc,
     2				     rab,rac,rbc,rab2,rac2,rbc2,
     3                               surfa,surfb,surfc,vola,volb,volc,
     4                               dsurfa3,dsurfb3,dsurfc3,dvola3,
     5                               dvolb3,dvolc3)
                                     dwsurf_dist(iab)=dwsurf_dist(iab)+
     1                               coefval*(coefa*dsurfa3(1)+coefb*
     2                               dsurfb3(1)+coefc*dsurfc3(1))
                                     dwsurf_dist(iac)=dwsurf_dist(iac)+
     1                               coefval*(coefa*dsurfa3(2)+coefb*
     2                               dsurfb3(2)+coefc*dsurfc3(2))
                                     dwsurf_dist(ibc)=dwsurf_dist(ibc)+
     1                               coefval*(coefa*dsurfa3(3)+coefb*
     2                               dsurfb3(3)+coefc*dsurfc3(3))
                                     dwvol_dist(iab)=dwvol_dist(iab)+
     1                               coefval*(coefa*dvola3(1)+coefb*
     2                               dvolb3(1)+coefc*dvolc3(1))
                                     dwvol_dist(iac)=dwvol_dist(iac)+
     1                               coefval*(coefa*dvola3(2)+coefb*
     2                               dvolb3(2)+coefc*dvolc3(2))
                                     dwvol_dist(ibc)=dwvol_dist(ibc)+
     1                               coefval*(coefa*dvola3(3)+coefb*
     2                               dvolb3(3)+coefc*dvolc3(3))
				endif
c
				ballwsurf(ia) = ballwsurf(ia) + 
     1				coefval*surfa
				ballwsurf(ib) = ballwsurf(ib) + 
     1				coefval*surfb
				ballwsurf(ic) = ballwsurf(ic) + 
     1				coefval*surfc
				ballwsurf(id) = ballwsurf(id) + 
     1				coefval*surfd
				ballwvol(ia) = ballwvol(ia) + 
     1				coefval*vola
				ballwvol(ib) = ballwvol(ib) + 
     1				coefval*volb
				ballwvol(ic) = ballwvol(ic) + 
     1				coefval*volc
				ballwvol(id) = ballwvol(id) + 
     1				coefval*vold
c
			endif
		endif
c
700	    continue
c
c	    Now check edges
c					
	    do 1500 iedge = 1,6
c
		if(btest(tetra_mask(idx),iedge-1)) goto 1500
c
		if(.not.btest(tetra_edge(idx),iedge-1)) goto 1500
c
		surfa = 0
		surfb = 0
		surfc = 0
		surfd = 0
		vola = 0
		volb = 0
		volc = 0
		vold = 0
c
		if(iedge.eq.1) then
			call twosphere_dvol_dist(c,d,rc,rc2,rd,rd2,
     1			rcd,rcd2,surfc,surfd,volc,vold,dsurfc2,dsurfd2,
     2                  dvolc2,dvold2)
                        dwsurf_dist(icd)=dwsurf_dist(icd) -coefc*
     1                  dsurfc2 - coefd*dsurfd2
                        dwvol_dist(icd)=dwvol_dist(icd) -coefc*
     1                  dvolc2 - coefd*dvold2
		elseif(iedge.eq.2) then
			call twosphere_dvol_dist(b,d,rb,rb2,rd,rd2,
     1			rbd,rbd2,surfb,surfd,volb,vold,dsurfb2,dsurfd2,
     2                  dvolb2,dvold2)
                        dwsurf_dist(ibd)=dwsurf_dist(ibd) -coefb*
     1                  dsurfb2 - coefd*dsurfd2
                        dwvol_dist(ibd)=dwvol_dist(ibd) -coefb*
     1                  dvolb2 - coefd*dvold2
		elseif(iedge.eq.3) then
			call twosphere_dvol_dist(b,c,rb,rb2,rc,rc2,
     1			rbc,rbc2,surfb,surfc,volb,volc,dsurfb2,dsurfc2,
     2                  dvolb2,dvolc2)
                        dwsurf_dist(ibc)=dwsurf_dist(ibc) -coefb*
     1                  dsurfb2 - coefc*dsurfc2
                        dwvol_dist(ibc)=dwvol_dist(ibc) -coefb*
     1                  dvolb2 - coefc*dvolc2
		elseif(iedge.eq.4) then
			call twosphere_dvol_dist(a,d,ra,ra2,rd,rd2,
     1			rad,rad2,surfa,surfd,vola,vold,dsurfa2,dsurfd2,
     2                  dvola2,dvold2)
                        dwsurf_dist(iad)=dwsurf_dist(iad) -coefa*
     1                  dsurfa2 - coefd*dsurfd2
                        dwvol_dist(iad)=dwvol_dist(iad) -coefa*
     1                  dvola2 - coefd*dvold2
		elseif(iedge.eq.5) then
			call twosphere_dvol_dist(a,c,ra,ra2,rc,rc2,
     1			rac,rac2,surfa,surfc,vola,volc,dsurfa2,dsurfc2,
     2                  dvola2,dvolc2)
                        dwsurf_dist(iac)=dwsurf_dist(iac) -coefa*
     1                  dsurfa2 - coefc*dsurfc2
                        dwvol_dist(iac)=dwvol_dist(iac) -coefa*
     1                  dvola2 - coefc*dvolc2
		elseif(iedge.eq.6) then
			call twosphere_dvol_dist(a,b,ra,ra2,rb,rb2,
     1			rab,rab2,surfa,surfb,vola,volb,dsurfa2,dsurfb2,
     2                  dvola2,dvolb2)
                        dwsurf_dist(iab)=dwsurf_dist(iab) -coefa*
     1                  dsurfa2 - coefb*dsurfb2
                        dwvol_dist(iab)=dwvol_dist(iab) -coefa*
     1                  dvola2 - coefb*dvolb2
		endif
c
		ballwsurf(ia) = ballwsurf(ia) - surfa
		ballwsurf(ib) = ballwsurf(ib) - surfb
		ballwsurf(ic) = ballwsurf(ic) - surfc
		ballwsurf(id) = ballwsurf(id) - surfd
		ballwvol(ia) = ballwvol(ia) - vola
		ballwvol(ib) = ballwvol(ib) - volb
		ballwvol(ic) = ballwvol(ic) - volc
		ballwvol(id) = ballwvol(id) - vold
c
		icheck = 0
c
c		iedge is the edge number in the tetrahedron idx, with:
c		iedge = 1		(c,d)
c		iedge = 2		(b,d)
c		iedge = 3		(b,c)
c		iedge = 4		(a,d)
c		iedge = 5		(a,c)
c		iedge = 6		(a,b)
c		
c		Define indices of the edge
c
		i = tetra(pair(1,iedge),idx)
		j = tetra(pair(2,iedge),idx)
c
c		trig1 and trig2 are the two faces of idx that share
c		iedge
c		i1 and i2 are the positions of the third vertices of
c		trig1 and trig2
c
		trig1 = face_info(1,iedge)
		i1    = face_pos(1,iedge)
		trig2 = face_info(2,iedge)
		i2    = face_pos(2,iedge)
c
		i3 = tetra(i1,idx)
		i4 = tetra(i2,idx)
c
		call mvbits(tetra_info(idx),7,1,it1,0)
c
		if(btest(tetra_info(idx),2+trig2)) then
			icheck = icheck + 1
			listcheck(icheck) = i4
			jtetra = tetra_neighbour(trig2,idx)
			if(jtetra.ne.0) then
			    call mvbits(tetra_info(jtetra),7,1,it2,0)
			    tag(icheck) = it1 + it2
			else
			    tag(icheck) = it1
			endif
		endif
c
		if(btest(tetra_info(idx),2+trig1)) then
			icheck = icheck + 1
			listcheck(icheck) = i3
			jtetra = tetra_neighbour(trig1,idx)
			if(jtetra.ne.0) then
			    call mvbits(tetra_info(jtetra),7,1,it2,0)
			    tag(icheck) = it1 + it2
			else
			    tag(icheck) = it1
			endif
		endif
c
c		Now we look at the star of the edge:
c
		ktetra = idx
		npass = 1
		trig_out = trig1
		jtetra = tetra_neighbour(trig_out,ktetra)
c
800		continue
c
c
c		Leave this side of the star if we hit the convex hull
c
		if(jtetra.eq.0) goto 900
		call mvbits(tetra_info(jtetra),7,1,it1,0)
c
c		Leave the loop completely if we have described the
c		full cycle
c
		if(jtetra.eq.idx) goto 1000
c
c		Identify the position of iedge in tetrahedron jtetra
c
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
c
		tetra_mask(jtetra) = ibset(tetra_mask(jtetra),ipair-1)
c
c		Find out the face we "went in"
c
		ival = ibits(tetra_nindex(ktetra),2*(trig_out-1),2)
		trig_in = ival + 1
c
c		We know the two faces of jtetra that share iedge:
c
		triga = face_info(1,ipair)
		i1    = face_pos(1,ipair)
		trigb = face_info(2,ipair)
		i2    = face_pos(2,ipair)
c
		trig_out = triga
		i_out = i1
		if(trig_in.eq.triga) then
			i_out = i2
			trig_out = trigb
		endif
c
c
		ktetra = jtetra
		jtetra = tetra_neighbour(trig_out,ktetra)
		if(jtetra.eq.idx) goto 1000
c
		if(btest(tetra_info(ktetra),2+trig_out)) then
			icheck = icheck + 1
			listcheck(icheck) = tetra(i_out,ktetra)
			if(jtetra.ne.0) then
			    call mvbits(tetra_info(jtetra),7,1,it2,0)
			    tag(icheck) = it1 + it2
			else
			    tag(icheck) = it1
			endif
		endif
c
		goto 800
c
900		continue
c
		if(npass.eq.2) goto 1000
		npass = npass + 1
		ktetra = idx
		trig_out = trig2
		jtetra = tetra_neighbour(trig_out,ktetra)
		goto 800
c
1000		continue
c
c		Now we have the list of triangles in the alpha complex
c		that are attached to the edge
c
c		If the edge is fully covered, skip!
c
		sum = 0
		do 1100 l = 1,icheck
			sum = sum + tag(l)
1100		continue
c
		if(sum.eq.2*icheck) goto 1500
c
1500	    continue
c
1600	continue
c
c	Now loop over vertices
c
	nred = 0
	do 1700 i = 1,nvertex
c
		if(.not.btest(vertex_info(i),0)) goto 1700
c
c		if vertex is not in alpha-complex, nothing to do...
c
		if(.not.btest(vertex_info(i),7)) goto 1700
c
c		Vertex is in alpha complex...
c
		ra = radius(i)
		ballwsurf(i) = ballwsurf(i) + 4*pi*ra*ra
		ballwvol(i) = ballwvol(i) + 4*pi*ra*ra*ra/3
c		
1700	continue
c
c	Compute total surface (weighted, and unweighted):
c
	do 1800 i = 1,nvertex
		WSurf        = Wsurf + coef(i)*ballwsurf(i)
		WVol         = WVol + coef(i)*ballwvol(i)
		Surf         = Surf  + ballwsurf(i)
		Vol          = Vol  + ballwvol(i)
		ballwsurf(i) = ballwsurf(i)*coef(i)
		ballwvol(i)  = ballwvol(i)*coef(i)
1800	continue
c
        do 2200 i = 1,nvertex-1
c
                do 2100 j = i+1,nvertex
c
                        iab = ipos(i,j)
c
                        if(dwsurf_dist(iab).ne.0) then
c
                                call distance2(coord,i,j,rab2,ncortot)
                                rab = dsqrt(rab2)
c
                                do 1900 k = 1,3
                                        u_ab(k)=(coord(3*(i-1)+k)-
     1                                  coord(3*(j-1)+k))/rab
1900                            continue
c
                                do 2000 k = 1,3
                                        coefa=dwsurf_dist(iab)*u_ab(k)
                                        dwsurf(k,i) = dwsurf(k,i)+coefa
                                        dwsurf(k,j) = dwsurf(k,j)-coefa
                                        coefa=dwvol_dist(iab)*u_ab(k)
                                        dwvol(k,i) = dwvol(k,i)+coefa
                                        dwvol(k,j) = dwvol(k,j)-coefa
2000                            continue
c
                        endif
c
2100            continue
c
2200    continue
c
c	We are done...
c
	return
	end
c
c       permute_ccw.f
c
        subroutine permute_ccw(array)
c
        real*8  val,array(6)
c
        val = array(2)
        array(2) = array(3)
        array(3) = val
        val = array(4)
        array(4) = array(5)
        array(5) = val
c
        return
        end
c
        function ipos(ia,ib)
c
        integer ia,ib,ipos,i,j
c
        if(ia.le.ib) then
                i = ia
                j = ib
        else
                i = ib
                j = ia
        endif
c
        ipos= (j*(j-1))/2 + i - 1
c
        return
        end
