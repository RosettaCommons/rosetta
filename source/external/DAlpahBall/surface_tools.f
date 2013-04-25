c	surface_tools.f		Version 1 : 4/17/2007

c	Copyright (C) 2007 Patrice Koehl

c	This file contains a series of subroutines used to compute the
c	volume area of a union of balls, and optionally the derivatives
c	of the volume with respect to the coordinates of the centers of the ball.
c	As a side result, it also provides the surface area (and its derivatives)
c	of the union of balls.

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

c	twosphere_surf.f

c	This subroutine calculates the surface of the
c	intersection of two spheres; it is only called when the
c	intersection exists
c	
c	Copyright (C) 2007 Patrice Koehl

	subroutine twosphere_surf(a,b,ra,ra2,rb,rb2,rab,rab2,
     1			surfa,surfb)

c	Input:
c			a,b	: position of the centers of the 2 spheres
c			rab	: distance between the centers of the 2 spheres
c			rab2	: distance between the centers of the 2 spheres
c				  (squared)
c			ra,rb	: radii of sphere A and B, respectively
c			ra2 and rb2 are the squared of the quantities
c			above)
c	Output
c			surfa	: partial contribution of A to the total
c				  surface of the intersection
c			surfb	: partial contribution of B to the total
c				  surface of the intersection


	integer	i,option

	real*8	ra,rb,surfa,surfb
	real*8	vala,valb,lamda
	real*8	ra2,rb2,rab,rab2,ha,hb
	real*8	pi,twopi,precision

	real*8	a(3),b(3),c(3)

	common /constants/ pi,twopi,precision

	save

c	Get "center" of the two spheres

	call center2(a,b,ra2,rb2,rab2,c,lamda)

	valb = lamda*rab
	vala = rab-valb

c	Get height of the cap of sphere A occluded by sphere B

	ha = ra - vala

c	same for sphere B ...

	hb = rb - valb

c	Get surfaces of intersection

	surfa = twopi*ra*ha
	surfb = twopi*rb*hb

	return
	end

c	twosphere_dsurf_coord.f

c	This subroutine calculates the surface of the
c	intersection of two spheres; it is only called when the
c	intersection exists
c	
c	Copyright (C) 2007 Patrice Koehl

c	It also provides the derivative of the surface with
c	respect to the coordinates of the center of the atoms

	subroutine twosphere_dsurf_coord(a,b,ra,ra2,rb,rb2,rab,rab2,
     1			surfa,surfb,dsurfa,dsurfb)

c	Input:
c			a,b	: position of the centers of the 2 spheres
c			rab	: distance between the centers of the 2 spheres
c			rab2	: distance between the centers of the 2 spheres
c				  (squared)
c			ra,rb	: radii of sphere A and B, respectively
c			ra2 and rb2 are the squared of the quantities
c			above)
c	Output
c			surfa	: partial contribution of A to the total
c				  surface of the intersection
c			surfb	: partial contribution of B to the total
c				  surface of the intersection
c                       dsurfa,dsurfb: derivatives of surfa and surfb
c                                      with respect to the coordinates of A and B

	integer	i

	real*8	ra,rb,surfa,surfb
	real*8	vala,valb,lamda
	real*8	ra2,rb2,rab,rab2,ha,hb
	real*8	dera,derb,der1,der2,coef1,coef2
	real*8	pi,twopi,precision

	real*8	a(3),b(3),c(3)
	real*8  u_ab(3),dsurfa(3,2),dsurfb(3,2)

	common /constants/ pi,twopi,precision

	save

c	Get "center" of the two spheres

	call center2(a,b,ra2,rb2,rab2,c,lamda)

	valb = lamda*rab
	vala = rab-valb

c	Get height of the cap of sphere A occluded by sphere B

	ha = ra - vala

c	same for sphere B ...

	hb = rb - valb

c	Get surfaces of intersection

	surfa = twopi*ra*ha
	surfb = twopi*rb*hb

c       Compute derivatives with respect to atom centers

	do 50 i = 1,3
		u_ab(i) = (a(i)-b(i))/rab
50	continue

	dera = -lamda
	derb = lamda-1

	der1 = twopi*ra*dera
	der2 = twopi*rb*derb

	do 100 i = 1,3
		coef1 = der1*u_ab(i)
		coef2 = der2*u_ab(i)
		dsurfa(i,1) = coef1
		dsurfa(i,2) = -coef1
		dsurfb(i,1) = coef2
		dsurfb(i,2) = -coef2
100	continue

	return
	end
c	twosphere_dsurf_dist.f

c	This subroutine calculates the surface of the
c	intersection of two spheres; it is only called when the
c	intersection exists
c	
c	Copyright (C) 2007 Patrice Koehl

c	It also provides the derivative of the surface with
c	respect to the internal distance Rab

	subroutine twosphere_dsurf_dist(a,b,ra,ra2,rb,rb2,rab,rab2,
     1			surfa,surfb,dsurfa,dsurfb)

c	Input:
c			a,b	: position of the centers of the 2 spheres
c			rab	: distance between the centers of the 2 spheres
c			rab2	: distance between the centers of the 2 spheres
c				  (squared)
c			ra,rb	: radii of sphere A and B, respectively
c			ra2 and rb2 are the squared of the quantities
c			above)
c	Output
c			surfa	: partial contribution of A to the total
c				  surface of the intersection
c			surfb	: partial contribution of B to the total
c				  surface of the intersection
c                       dsurfa,dsurfb: derivatives of surfa, surfb with respect
c                                       to Rab

	integer	i

	real*8	ra,rb,surfa,surfb
	real*8	vala,valb,lamda
	real*8	ra2,rb2,rab,rab2,ha,hb
	real*8	dera,derb
	real*8	pi,twopi,precision
        real*8  dsurfa,dsurfb

	real*8	a(3),b(3),c(3)

	common /constants/ pi,twopi,precision

	save

c	Get "center" of the two spheres

	call center2(a,b,ra2,rb2,rab2,c,lamda)

	valb = lamda*rab
	vala = rab-valb

c	Get height of the cap of sphere A occluded by sphere B

	ha = ra - vala

c	same for sphere B ...

	hb = rb - valb

c	Get surfaces of intersection

	surfa = twopi*ra*ha
	surfb = twopi*rb*hb

c       Compute derivatives with respect to Rab

	dera = -lamda
	derb = lamda-1

	dsurfa = twopi*ra*dera
	dsurfb = twopi*rb*derb

	return
	end

c	threesphere_surf.f

c	Copyright (C) 2007 Patrice Koehl

c	This subroutine calculates the surface of the intersection 
c	of three spheres; it is only called when the intersection exists

	subroutine threesphere_surf(a,b,c,ra,rb,rc,ra2,rb2,rc2,
     1		wa,wb,wc,rab,rac,rbc,rab2,rac2,rbc2,surfa,surfb,
     2		surfc)

c	Input:
c			ra,rb,rc  : radii of sphere A, B and C, respectively
c			ra2,rb2,rc2: radii squared
c			wa,wb,wc   : "weights" of A, B and C
c			(for example, wa = 0.5*(xa**2+ya**2+za**2-ra**2)
c			rab,rab2: distance between the centers of sphere A and B
c			rac,rac2: distance between the centers of sphere A and C
c			rbc,rbc2: distance between the centers of sphere B and C
c	Output
c			surfa,surfb,surfc : contribution of A, B and C to
c			the total surface of the intersection of A,B,C


	integer	i,option

	real*8	surfa,surfb,surfc
	real*8	ra,rb,rc,rab,rac,rbc,rab2,rac2,rbc2
	real*8	ra2,rb2,rc2,wa,wb,wc
	real*8	a1,a2,a3,eps
	real*8	seg_ang_abc,seg_ang_acb
	real*8	seg_ang_bca
	real*8	ang_dih_abc,ang_dih_cab
	real*8	ang_dih_bac
	real*8	val1,val2,val3,l1,l2,l3
	real*8	val1b,val2b,val3b
	real*8	pi,twopi,precision

	real*8	a(3),b(3),c(3),center(3),n(3)
        real*8  c_ab(3),c_ac(3),c_bc(3)
	real*8	pabc(3),pacb(3)
        real*8  cosine(3),sine(3)

	common /constants/ pi,twopi,precision

	save

	call center2(a,b,ra2,rb2,rab2,c_ab,l1)
	call center2(a,c,ra2,rc2,rac2,c_ac,l2)
	call center2(b,c,rb2,rc2,rbc2,c_bc,l3)

	val1 = l1*rab
	val2 = l2*rac
	val3 = l3*rbc

	val1b = rab-val1
	val2b = rac-val2
	val3b = rbc-val3

	call center3(a,b,c,wa,wb,wc,center)

	call triangle_dual(a,b,c,center,eps,ra2,pabc,pacb,n)

	option = 1
	call tetra3_noder(a,b,c,pabc,rab,rac,rbc,
     1	ra,rb,rc,seg_ang_abc,seg_ang_acb,
     2	seg_ang_bca,ang_dih_abc,ang_dih_bac,ang_dih_cab,
     2	cosine,sine,option)

	a1 = ra*(1-2*ang_dih_abc)
	a2 = 2*seg_ang_abc*val1b
	a3 = 2*seg_ang_acb*val2b

	surfa = twopi*ra*(a1 - a2 - a3)

	a1 = rb*(1-2*ang_dih_bac)
	a2 = 2*seg_ang_abc*val1
	a3 = 2*seg_ang_bca*val3b

	surfb = twopi*rb*(a1 - a2 - a3)

	a1 = rc*(1-2*ang_dih_cab)
	a2 = 2*seg_ang_acb*val2
	a3 = 2*seg_ang_bca*val3

	surfc = twopi*rc*(a1 - a2 - a3)

	return
	end

c	threesphere_dsurf_coord.f	Version 1 08/25/2000	Patrice Koehl

c	This subroutine calculates the surface of the intersection 
c	of three spheres; it is only called when it exists
c	The computation follows the method described in VOLBL by Herbert
c	Edelsbrunner

c	This subroutine supposes that each pair of spheres intersect (i.e.
c	a three sphere intersection is possible)

	subroutine threesphere_dsurf_coord(a,b,c,ra,rb,rc,ra2,rb2,rc2,
     1		wa,wb,wc,rab,rac,rbc,rab2,rac2,rbc2,surfa,surfb,
     2		surfc,dsurfa,dsurfb,dsurfc)


c	Input:
c			ra,rb,rc  : radii of sphere A, B and C, respectively
c			ra2,rb2,rc2: radii squared
c			rab,rab2: distance between the centers of sphere A and B
c			rac,rac2: distance between the centers of sphere A and C
c			rbc,rbc2: distance between the centers of sphere B and C
c	Output
c			surfa,surfb,surfc : contribution of A, B and C to
c			the total surface of the intersection of A,B,C
c			dsurfa,dsurfb,dsurfc : derivatives of surfa, surfb, surfc
c                                              with respect to the coordinates
c                                              A, B and C

	integer	i

	real*8	surfa,surfb,surfc
	real*8	ra,rb,rc,rab,rac,rbc,rab2,rac2,rbc2
	real*8	ra2,rb2,rc2,wa,wb,wc
	real*8	a1,a2,a3,s2,c1,c2,eps
	real*8	seg_ang_abc,seg_ang_acb
	real*8	seg_ang_bca
	real*8	ang_abc,ang_acb,ang_bca
	real*8	cos_abc,cos_acb,cos_bca
	real*8	sin_abc,sin_acb,sin_bca
	real*8	ang_dih_abc,ang_dih_cab
	real*8	ang_dih_bac
	real*8	val1,val2,val3,l1,l2,l3
	real*8	val1b,val2b,val3b
	real*8	der_ab,der_ac,der_bc,diff_ab,diff_ac,diff_bc
	real*8	dsurfa_ab,dsurfa_ac,dsurfa_bc
	real*8	dsurfb_ab,dsurfb_ac,dsurfb_bc
	real*8	dsurfc_ab,dsurfc_ac,dsurfc_bc
	real*8	coef_ab,coef_ac,coef_bc

	real*8	a(3),b(3),c(3),center(3)
        real*8  c_ab(3),c_ac(3),c_bc(3)
	real*8	pabc(3),pacb(3),n(3)
	real*8	der1(3),der2(3),der3(3),der4(3),der5(3),der6(3)
	real*8	dsurfa(3,3),dsurfb(3,3),dsurfc(3,3)
	real*8	u_ab(3),u_ac(3),u_bc(3)
	real*8	cosine(3),sine(3)

	real*8	pi,twopi,precision

	common /constants/ pi,twopi,precision

	save

	call center2(a,b,ra2,rb2,rab2,c_ab,l1)
	call center2(a,c,ra2,rc2,rac2,c_ac,l2)
	call center2(b,c,rb2,rc2,rbc2,c_bc,l3)

	val1 = l1*rab
	val2 = l2*rac
	val3 = l3*rbc

	val1b = rab-val1
	val2b = rac-val2
	val3b = rbc-val3

	call center3(a,b,c,wa,wb,wc,center)

	call triangle_dual(a,b,c,center,eps,ra2,pabc,pacb,n)

	call tetra3(a,b,c,pabc,rab,rac,rbc,rab2,rac2,
     1	rbc2,ra,rb,rc,ra2,rb2,rc2,seg_ang_abc,seg_ang_acb,
     2	seg_ang_bca,ang_dih_abc,ang_dih_bac,ang_dih_cab,
     2	der1,der2,der3,der4,der5,der6,cosine,sine)

	a1 = ra*(1-2*ang_dih_abc)
	a2 = 2*seg_ang_abc*val1b
	a3 = 2*seg_ang_acb*val2b

	surfa = twopi*ra*(a1 - a2 - a3)

	a1 = rb*(1-2*ang_dih_bac)
	a2 = 2*seg_ang_abc*val1
	a3 = 2*seg_ang_bca*val3b

	surfb = twopi*rb*(a1 - a2 - a3)

	a1 = rc*(1-2*ang_dih_cab)
	a2 = 2*seg_ang_acb*val2
	a3 = 2*seg_ang_bca*val3

	surfc = twopi*rc*(a1 - a2 - a3)

	do 50 i = 1,3
		u_ab(i) = (a(i)-b(i))/rab
		u_ac(i) = (a(i)-c(i))/rac
		u_bc(i) = (b(i)-c(i))/rbc
50	continue

	ang_abc = twopi*seg_ang_abc
	ang_acb = twopi*seg_ang_acb
	ang_bca = twopi*seg_ang_bca

	dsurfa_ab = -2*ra*(ra*der4(1) + der1(1)*val1b +
     1		ang_abc*l1 +der2(1)*val2b)
	dsurfa_ac = - 2*ra*(ra*der4(2) +der1(2)*val1b +
     1		der2(2)*val2b + ang_acb*l2)
	dsurfa_bc = -2*ra*(ra*der4(3) +der1(3)*val1b +
     1			der2(3)*val2b)

	do 100 i = 1,3
		diff_ab = dsurfa_ab*u_ab(i)
		diff_ac = dsurfa_ac*u_ac(i)
		diff_bc = dsurfa_bc*u_bc(i)
		dsurfa(i,1) = diff_ab + diff_ac
		dsurfa(i,2) = -diff_ab + diff_bc
		dsurfa(i,3) = -diff_ac - diff_bc
100	continue

	dsurfb_ab = -2*rb*(rb*der5(1) + der1(1)*val1 +
     1		ang_abc*(1-l1) +der3(1)*val3b)
	dsurfb_ac = - 2*rb*(rb*der5(2) +der1(2)*val1 +
     1		der3(2)*val3b )
	dsurfb_bc = -2*rb*(rb*der5(3) +der1(3)*val1 +
     1			der3(3)*val3b+ang_bca*l3)

	do 200 i = 1,3
		diff_ab = dsurfb_ab*u_ab(i)
		diff_ac = dsurfb_ac*u_ac(i)
		diff_bc = dsurfb_bc*u_bc(i)
		dsurfb(i,1) = diff_ab + diff_ac
		dsurfb(i,2) = -diff_ab + diff_bc
		dsurfb(i,3) = -diff_ac - diff_bc
200	continue

	dsurfc_ab = -2*rc*(rc*der6(1) + der2(1)*val2 +
     1		der3(1)*val3)
	dsurfc_ac = - 2*rc*(rc*der6(2) +der2(2)*val2 +
     1		ang_acb*(1-l2) + der3(2)*val3 )
	dsurfc_bc = -2*rc*(rc*der6(3) +der2(3)*val2 +
     1			der3(3)*val3+ang_bca*(1-l3))

	do 300 i = 1,3
		diff_ab = dsurfc_ab*u_ab(i)
		diff_ac = dsurfc_ac*u_ac(i)
		diff_bc = dsurfc_bc*u_bc(i)
		dsurfc(i,1) = diff_ab + diff_ac
		dsurfc(i,2) = -diff_ab + diff_bc
		dsurfc(i,3) = -diff_ac - diff_bc
300	continue

	return
	end
c	threesphere_dsurf_dist.f	Version 1 9/1/2008	Patrice Koehl

c	This subroutine calculates the surface of the intersection 
c	of three spheres; it is only called when it exists
c	The computation follows the method described in VOLBL by Herbert
c	Edelsbrunner

c	This subroutine supposes that each pair of spheres intersect (i.e.
c	a three sphere intersection is possible)

	subroutine threesphere_dsurf_dist(a,b,c,ra,rb,rc,ra2,rb2,rc2,
     1		wa,wb,wc,rab,rac,rbc,rab2,rac2,rbc2,surfa,surfb,
     2		surfc,dsurfa,dsurfb,dsurfc)

c	Input:
c			ra,rb,rc  : radii of sphere A, B and C, respectively
c			ra2,rb2,rc2: radii squared
c			rab,rab2: distance between the centers of sphere A and B
c			rac,rac2: distance between the centers of sphere A and C
c			rbc,rbc2: distance between the centers of sphere B and C
c	Output
c			surfa,surfb,surfc : contribution of A, B and C to
c			the total surface of the intersection of A,B,C
c			dsurfa,dsurfb,dsurfc : derivatives of surfa, surfb, surfc
c                                              with respect to Rab, Rac, Rbc

	integer	i

	real*8	surfa,surfb,surfc
	real*8	ra,rb,rc,rab,rac,rbc,rab2,rac2,rbc2
	real*8	ra2,rb2,rc2,wa,wb,wc
	real*8	a1,a2,a3,s2,c1,c2,eps
	real*8	seg_ang_abc,seg_ang_acb
	real*8	seg_ang_bca
	real*8	ang_abc,ang_acb,ang_bca
	real*8	cos_abc,cos_acb,cos_bca
	real*8	sin_abc,sin_acb,sin_bca
	real*8	ang_dih_abc,ang_dih_cab
	real*8	ang_dih_bac
	real*8	val1,val2,val3,l1,l2,l3
	real*8	val1b,val2b,val3b
	real*8	der_ab,der_ac,der_bc,diff_ab,diff_ac,diff_bc
	real*8	dsurfa_ab,dsurfa_ac,dsurfa_bc
	real*8	dsurfb_ab,dsurfb_ac,dsurfb_bc
	real*8	dsurfc_ab,dsurfc_ac,dsurfc_bc
	real*8	coef_ab,coef_ac,coef_bc

	real*8	a(3),b(3),c(3),center(3)
        real*8  c_ab(3),c_ac(3),c_bc(3)
	real*8	pabc(3),pacb(3),n(3)
	real*8	der1(3),der2(3),der3(3),der4(3),der5(3),der6(3)
	real*8	dsurfa(3),dsurfb(3),dsurfc(3)
	real*8	cosine(3),sine(3)

	real*8	pi,twopi,precision

	common /constants/ pi,twopi,precision

	save

	call center2(a,b,ra2,rb2,rab2,c_ab,l1)
	call center2(a,c,ra2,rc2,rac2,c_ac,l2)
	call center2(b,c,rb2,rc2,rbc2,c_bc,l3)

	val1 = l1*rab
	val2 = l2*rac
	val3 = l3*rbc

	val1b = rab-val1
	val2b = rac-val2
	val3b = rbc-val3

	call center3(a,b,c,wa,wb,wc,center)

	call triangle_dual(a,b,c,center,eps,ra2,pabc,pacb,n)

	call tetra3(a,b,c,pabc,rab,rac,rbc,rab2,rac2,
     1	rbc2,ra,rb,rc,ra2,rb2,rc2,seg_ang_abc,seg_ang_acb,
     2	seg_ang_bca,ang_dih_abc,ang_dih_bac,ang_dih_cab,
     2	der1,der2,der3,der4,der5,der6,cosine,sine)

	a1 = ra*(1-2*ang_dih_abc)
	a2 = 2*seg_ang_abc*val1b
	a3 = 2*seg_ang_acb*val2b

	surfa = twopi*ra*(a1 - a2 - a3)

	a1 = rb*(1-2*ang_dih_bac)
	a2 = 2*seg_ang_abc*val1
	a3 = 2*seg_ang_bca*val3b

	surfb = twopi*rb*(a1 - a2 - a3)

	a1 = rc*(1-2*ang_dih_cab)
	a2 = 2*seg_ang_acb*val2
	a3 = 2*seg_ang_bca*val3

	surfc = twopi*rc*(a1 - a2 - a3)

	ang_abc = twopi*seg_ang_abc
	ang_acb = twopi*seg_ang_acb
	ang_bca = twopi*seg_ang_bca

	dsurfa_ab = -2*ra*(ra*der4(1) + der1(1)*val1b +
     1		ang_abc*l1 +der2(1)*val2b)
	dsurfa_ac = - 2*ra*(ra*der4(2) +der1(2)*val1b +
     1		der2(2)*val2b + ang_acb*l2)
	dsurfa_bc = -2*ra*(ra*der4(3) +der1(3)*val1b +
     1			der2(3)*val2b)

        dsurfa(1) = dsurfa_ab
        dsurfa(2) = dsurfa_ac
        dsurfa(3) = dsurfa_bc

	dsurfb_ab = -2*rb*(rb*der5(1) + der1(1)*val1 +
     1		ang_abc*(1-l1) +der3(1)*val3b)
	dsurfb_ac = - 2*rb*(rb*der5(2) +der1(2)*val1 +
     1		der3(2)*val3b )
	dsurfb_bc = -2*rb*(rb*der5(3) +der1(3)*val1 +
     1			der3(3)*val3b+ang_bca*l3)

        dsurfb(1) = dsurfb_ab
        dsurfb(2) = dsurfb_ac
        dsurfb(3) = dsurfb_bc

	dsurfc_ab = -2*rc*(rc*der6(1) + der2(1)*val2 +
     1		der3(1)*val3)
	dsurfc_ac = - 2*rc*(rc*der6(2) +der2(2)*val2 +
     1		ang_acb*(1-l2) + der3(2)*val3 )
	dsurfc_bc = -2*rc*(rc*der6(3) +der2(3)*val2 +
     1			der3(3)*val3+ang_bca*(1-l3))

        dsurfc(1) = dsurfc_ab
        dsurfc(2) = dsurfc_ac
        dsurfc(3) = dsurfc_bc

	return
	end
c	Foursphere_surf.f

c	This subroutine calculates the volume and surface area of the
c	intersection of four spheres; this intersection is supposed to exist

c	This routine assumes that the 4 points (a,b,c,d) are in ccw order

c	Copyright (C) 2002 Patrice Koehl

	subroutine foursphere_surf(a,b,c,d,ra,rb,rc,rd,
     1			ra2,rb2,rc2,rd2,rab,rac,rad,rbc,rbd,rcd,
     2			rab2,rac2,rad2,rbc2,rbd2,rcd2,wa,wb,wc,wd,
     3			surfa,surfb,surfc,surfd)

	integer	i,option

	real*8	pi,twopi,precision
	real*8	ra,rb,rc,rd,ra2,rb2,rc2,rd2
	real*8	rab,rac,rad,rbc,rbd,rcd
	real*8	rab2,rac2,rad2,rbc2,rbd2,rcd2
	real*8	surfa,surfb,surfc,surfd
	real*8	wa,wb,wc,wd
	real*8	val_ab,val_ac,val_ad,val_bc,val_bd,val_cd
	real*8	val2_ab,val2_ac,val2_ad,val2_bc,val2_bd,val2_cd
	real*8	ang1,ang2,ang3,ang4,ang5,ang6
	real*8	l_ab,l_ac,l_ad,l_bc,l_bd,l_cd

	real*8	a(3),b(3),c(3),d(3)
	real*8	c_ab(3),c_ac(3),c_ad(3),c_bc(3),c_bd(3),c_cd(3)
	real*8	c_abcd(3)
	real*8	pacb(3),pabd(3),padc(3),pbcd(3)
	real*8	cosine(6),sine(6)

	common /constants/ pi,twopi,precision

	save

	call tetra_6dihed(a,b,c,d,ang1,ang2,ang4,ang3,ang5,ang6)

	call center2(a,b,ra2,rb2,rab2,c_ab,l_ab)
	call center2(a,c,ra2,rc2,rac2,c_ac,l_ac)
	call center2(a,d,ra2,rd2,rad2,c_ad,l_ad)
	call center2(b,c,rb2,rc2,rbc2,c_bc,l_bc)
	call center2(b,d,rb2,rd2,rbd2,c_bd,l_bd)
	call center2(c,d,rc2,rd2,rcd2,c_cd,l_cd)

	val_ab = l_ab*rab
	val_ac = l_ac*rac
	val_ad = l_ad*rad
	val_bc = l_bc*rbc
	val_bd = l_bd*rbd
	val_cd = l_cd*rcd

	val2_ab = rab - val_ab
	val2_ac = rac - val_ac
	val2_ad = rad - val_ad
	val2_bc = rbc - val_bc
	val2_bd = rbd - val_bd
	val2_cd = rcd - val_cd

	surfa = -0.5d0*ra + ang1*val2_ab + ang2*val2_ac +
     1		ang3*val2_ad
	surfa = twopi*ra*surfa

	surfb = -0.5d0*rb + ang1*val_ab + ang5*val2_bd +
     1		ang4*val2_bc
	surfb = twopi*rb*surfb

	surfc = -0.5d0*rc + ang2*val_ac + ang4*val_bc +
     1		ang6*val2_cd
	surfc = twopi*rc*surfc

	surfd = -0.5d0*rd + ang3*val_ad + ang6*val_cd +
     1		ang5*val_bd
	surfd = twopi*rd*surfd

	return
	end

c	Foursphere_dsurf_coord.f	Version 1 08/28/2000	Patrice Koehl

c	This subroutine calculates the surface area of the
c	intersection of four spheres; this intersection is supposed to exist
c	The computation follows the method described in VOLBL by Herbert
c	Edelsbrunner

c       It also computes the derivatives of the surface area with respect
c       to the coordinates of A, B, C and D

	subroutine foursphere_dsurf_coord(a,b,c,d,ra,rb,rc,rd,
     1			ra2,rb2,rc2,rd2,rab,rac,rad,rbc,rbd,rcd,
     2			rab2,rac2,rad2,rbc2,rbd2,rcd2,wa,wb,wc,wd,
     7			surfa,surfb,surfc,surfd,
     8			dsurfa,dsurfb,dsurfc,dsurfd)

	integer	i,option

	real*8	pi,twopi,precision
	real*8	ra,rb,rc,rd,ra2,rb2,rc2,rd2
	real*8	rab,rac,rad,rbc,rbd,rcd
	real*8	rab2,rac2,rad2,rbc2,rbd2,rcd2
	real*8	surfa,surfb,surfc,surfd
	real*8	wa,wb,wc,wd
	real*8	val_ab,val_ac,val_ad,val_bc,val_bd,val_cd
	real*8	val2_ab,val2_ac,val2_ad,val2_bc,val2_bd,val2_cd
	real*8	ang1,ang2,ang3,ang4,ang5,ang6
	real*8	l_ab,l_ac,l_ad,l_bc,l_bd,l_cd
	real*8	dsurfa_ab,dsurfa_ac,dsurfa_ad,dsurfa_bc
	real*8	dsurfa_bd,dsurfa_cd
	real*8	dsurfb_ab,dsurfb_ac,dsurfb_ad,dsurfb_bc
	real*8	dsurfb_bd,dsurfb_cd
	real*8	dsurfc_ab,dsurfc_ac,dsurfc_ad,dsurfc_bc
	real*8	dsurfc_bd,dsurfc_cd
	real*8	dsurfd_ab,dsurfd_ac,dsurfd_ad,dsurfd_bc
	real*8	dsurfd_bd,dsurfd_cd
	real*8	der_ab,der_ac,der_ad,der_bc,der_bd,der_cd
	real*8	diff_ab,diff_ac,diff_ad,diff_bc,diff_bd,diff_cd

	real*8	a(3),b(3),c(3),d(3)
	real*8	c_ab(3),c_ac(3),c_ad(3),c_bc(3),c_bd(3),c_cd(3)
	real*8	c_abcd(3)
	real*8  der1(6),der2(6),der3(6),der4(6),der5(6),der6(6)
	real*8	dsurfa(3,4),dsurfb(3,4),dsurfc(3,4),dsurfd(3,4)
	real*8	cosine(6),sine(6)

	real*8	u_ab(3),u_ac(3),u_ad(3),u_bc(3),u_bd(3),u_cd(3)


	common /constants/ pi,twopi,precision

	save

	do 30 i = 1,3
		u_ab(i) = (a(i)-b(i))/rab
		u_ac(i) = (a(i)-c(i))/rac
		u_ad(i) = (a(i)-d(i))/rad
		u_bc(i) = (b(i)-c(i))/rbc
		u_bd(i) = (b(i)-d(i))/rbd
		u_cd(i) = (c(i)-d(i))/rcd
30	continue

	call tetra6(a,b,c,d,rab,rac,rad,rbc,rbd,rcd,rab2,rac2,
     1		rad2,rbc2,rbd2,rcd2,ang1,ang2,ang3,
     2		ang4,ang5,ang6,der1,der2,der3,der4,der5,der6,
     3		cosine,sine)

	call center2(a,b,ra2,rb2,rab2,c_ab,l_ab)
	call center2(a,c,ra2,rc2,rac2,c_ac,l_ac)
	call center2(a,d,ra2,rd2,rad2,c_ad,l_ad)
	call center2(b,c,rb2,rc2,rbc2,c_bc,l_bc)
	call center2(b,d,rb2,rd2,rbd2,c_bd,l_bd)
	call center2(c,d,rc2,rd2,rcd2,c_cd,l_cd)

	val_ab = l_ab*rab
	val_ac = l_ac*rac
	val_ad = l_ad*rad
	val_bc = l_bc*rbc
	val_bd = l_bd*rbd
	val_cd = l_cd*rcd

	val2_ab = rab - val_ab
	val2_ac = rac - val_ac
	val2_ad = rad - val_ad
	val2_bc = rbc - val_bc
	val2_bd = rbd - val_bd
	val2_cd = rcd - val_cd

	surfa = -0.5d0*ra + ang1*val2_ab + ang2*val2_ac +
     1		ang3*val2_ad

	surfa = twopi*ra*surfa

	dsurfa_ab = ra*(der1(1)*val2_ab+twopi*ang1*l_ab+
     1		der2(1)*val2_ac + der3(1)*val2_ad)
	dsurfa_ac = ra*(der1(2)*val2_ab+der2(2)*val2_ac +
     1		twopi*ang2*l_ac + der3(2)*val2_ad)
	dsurfa_ad = ra*(der1(3)*val2_ab + der2(3)*val2_ac +
     1		der3(3)*val2_ad + twopi*ang3*l_ad)
	dsurfa_bc = ra*(der1(4)*val2_ab + der2(4)*val2_ac +
     1		 der3(4)*val2_ad)
	dsurfa_bd = ra*(der1(5)*val2_ab + der2(5)*val2_ac +
     1		 der3(5)*val2_ad)
	dsurfa_cd = ra*(der1(6)*val2_ab + der2(6)*val2_ac +
     1		 der3(6)*val2_ad)

	do 100 i = 1,3
		diff_ab = dsurfa_ab*u_ab(i)
		diff_ac = dsurfa_ac*u_ac(i)
		diff_ad = dsurfa_ad*u_ad(i)
		diff_bc = dsurfa_bc*u_bc(i)
		diff_bd = dsurfa_bd*u_bd(i)
		diff_cd = dsurfa_cd*u_cd(i)
		dsurfa(i,1) = diff_ab+diff_ac+ diff_ad
		dsurfa(i,2) = -diff_ab+diff_bc+diff_bd
		dsurfa(i,3) = -diff_ac-diff_bc+diff_cd
		dsurfa(i,4) = -diff_ad-diff_bd-diff_cd
100	continue

	surfb = -0.5d0*rb + ang1*val_ab + ang5*val2_bd +
     1		ang4*val2_bc
	surfb = twopi*rb*surfb

	dsurfb_ab = rb*(der1(1)*val_ab+twopi*ang1*(1-l_ab)+
     1		der5(1)*val2_bd +der4(1)*val2_bc)
	dsurfb_ac = rb*(der1(2)*val_ab+der5(2)*val2_bd +
     1		der4(2)*val2_bc)
	dsurfb_ad = rb*(der1(3)*val_ab + der5(3)*val2_bd +
     1		der4(3)*val2_bc)
	dsurfb_bc = rb*(der1(4)*val_ab + der5(4)*val2_bd +
     1		 der4(4)*val2_bc+twopi*ang4*l_bc)
	dsurfb_bd = rb*(der1(5)*val_ab +der5(5)*val2_bd +
     1		 twopi*ang5*l_bd + der4(5)*val2_bc)
	dsurfb_cd = rb*(der1(6)*val_ab + der5(6)*val2_bd +
     1		 der4(6)*val2_bc)

	do 200 i = 1,3
		diff_ab = dsurfb_ab*u_ab(i)
		diff_ac = dsurfb_ac*u_ac(i)
		diff_ad = dsurfb_ad*u_ad(i)
		diff_bc = dsurfb_bc*u_bc(i)
		diff_bd = dsurfb_bd*u_bd(i)
		diff_cd = dsurfb_cd*u_cd(i)
		dsurfb(i,1) = diff_ab+diff_ac+ diff_ad
		dsurfb(i,2) = -diff_ab+diff_bc+diff_bd
		dsurfb(i,3) = -diff_ac-diff_bc+diff_cd
		dsurfb(i,4) = -diff_ad-diff_bd-diff_cd
200	continue

	surfc = -0.5d0*rc + ang2*val_ac + ang4*val_bc +
     1		ang6*val2_cd
	surfc = twopi*rc*surfc

	dsurfc_ab = rc*(der2(1)*val_ac+
     1		der4(1)*val_bc + der6(1)*val2_cd)
	dsurfc_ac = rc*(der2(2)*val_ac+twopi*ang2*(1-l_ac)+
     1		der4(2)*val_bc + der6(2)*val2_cd)
	dsurfc_ad = rc*(der2(3)*val_ac + der4(3)*val_bc +
     1		der6(3)*val2_cd)
	dsurfc_bc = rc*(der2(4)*val_ac + der4(4)*val_bc +
     1		 twopi*ang4*(1-l_bc)+der6(4)*val2_cd)
	dsurfc_bd = rc*(der2(5)*val_ac + der4(5)*val_bc +
     1		 der6(5)*val2_cd)
	dsurfc_cd = rc*(der2(6)*val_ac +der4(6)*val_bc +
     1		 der6(6)*val2_cd + twopi*ang6*l_cd)

	do 300 i = 1,3
		diff_ab = dsurfc_ab*u_ab(i)
		diff_ac = dsurfc_ac*u_ac(i)
		diff_ad = dsurfc_ad*u_ad(i)
		diff_bc = dsurfc_bc*u_bc(i)
		diff_bd = dsurfc_bd*u_bd(i)
		diff_cd = dsurfc_cd*u_cd(i)
		dsurfc(i,1) = diff_ab+diff_ac+ diff_ad
		dsurfc(i,2) = -diff_ab+diff_bc+diff_bd
		dsurfc(i,3) = -diff_ac-diff_bc+diff_cd
		dsurfc(i,4) = -diff_ad-diff_bd-diff_cd
300	continue

	surfd = -0.5d0*rd + ang3*val_ad + ang6*val_cd +
     1		ang5*val_bd
	surfd = twopi*rd*surfd

	dsurfd_ab = rd*(der3(1)*val_ad+
     1		der6(1)*val_cd + der5(1)*val_bd)
	dsurfd_ac = rd*(der3(2)*val_ad+
     1		der6(2)*val_cd + der5(2)*val_bd)
	dsurfd_ad = rd*(der3(3)*val_ad + twopi*ang3*(1-l_ad) +
     1		der6(3)*val_cd + der5(3)*val_bd)
	dsurfd_bc = rd*(der3(4)*val_ad + der6(4)*val_cd +
     1		 der5(4)*val_bd)
	dsurfd_bd = rd*(der3(5)*val_ad + der6(5)*val_cd +
     1		 der5(5)*val_bd+ twopi*ang5*(1-l_bd))
	dsurfd_cd = rd*(der3(6)*val_ad + der6(6)*val_cd +
     1		 twopi*ang6*(1-l_cd) + der5(6)*val_bd)

	do 400 i = 1,3
		diff_ab = dsurfd_ab*u_ab(i)
		diff_ac = dsurfd_ac*u_ac(i)
		diff_ad = dsurfd_ad*u_ad(i)
		diff_bc = dsurfd_bc*u_bc(i)
		diff_bd = dsurfd_bd*u_bd(i)
		diff_cd = dsurfd_cd*u_cd(i)
		dsurfd(i,1) = diff_ab+diff_ac+ diff_ad
		dsurfd(i,2) = -diff_ab+diff_bc+diff_bd
		dsurfd(i,3) = -diff_ac-diff_bc+diff_cd
		dsurfd(i,4) = -diff_ad-diff_bd-diff_cd
400	continue

        return
        end

c	Foursphere_dsurf_dist.f	Version 1 08/28/2000	Patrice Koehl

c	This subroutine calculates the surface area of the
c	intersection of four spheres; this intersection is supposed to exist
c	The computation follows the method described in VOLBL by Herbert
c	Edelsbrunner

c       It also computes the derivatives of the surface area with respect
c       to the coordinates of A, B, C and D

	subroutine foursphere_dsurf_dist(a,b,c,d,ra,rb,rc,rd,
     1			ra2,rb2,rc2,rd2,rab,rac,rad,rbc,rbd,rcd,
     2			rab2,rac2,rad2,rbc2,rbd2,rcd2,wa,wb,wc,wd,
     7			surfa,surfb,surfc,surfd,
     8			dsurfa,dsurfb,dsurfc,dsurfd)

	integer	i,option

	real*8	pi,twopi,precision
	real*8	ra,rb,rc,rd,ra2,rb2,rc2,rd2
	real*8	rab,rac,rad,rbc,rbd,rcd
	real*8	rab2,rac2,rad2,rbc2,rbd2,rcd2
	real*8	surfa,surfb,surfc,surfd
	real*8	wa,wb,wc,wd
	real*8	val_ab,val_ac,val_ad,val_bc,val_bd,val_cd
	real*8	val2_ab,val2_ac,val2_ad,val2_bc,val2_bd,val2_cd
	real*8	ang1,ang2,ang3,ang4,ang5,ang6
	real*8	l_ab,l_ac,l_ad,l_bc,l_bd,l_cd
	real*8	dsurfa_ab,dsurfa_ac,dsurfa_ad,dsurfa_bc
	real*8	dsurfa_bd,dsurfa_cd
	real*8	dsurfb_ab,dsurfb_ac,dsurfb_ad,dsurfb_bc
	real*8	dsurfb_bd,dsurfb_cd
	real*8	dsurfc_ab,dsurfc_ac,dsurfc_ad,dsurfc_bc
	real*8	dsurfc_bd,dsurfc_cd
	real*8	dsurfd_ab,dsurfd_ac,dsurfd_ad,dsurfd_bc
	real*8	dsurfd_bd,dsurfd_cd
	real*8	der_ab,der_ac,der_ad,der_bc,der_bd,der_cd

	real*8	a(3),b(3),c(3),d(3)
	real*8	c_ab(3),c_ac(3),c_ad(3),c_bc(3),c_bd(3),c_cd(3)
	real*8	c_abcd(3)
	real*8  der1(6),der2(6),der3(6),der4(6),der5(6),der6(6)
	real*8	dsurfa(6),dsurfb(6),dsurfc(6),dsurfd(6)
	real*8	cosine(6),sine(6)


	common /constants/ pi,twopi,precision

	save

	call tetra6(a,b,c,d,rab,rac,rad,rbc,rbd,rcd,rab2,rac2,
     1		rad2,rbc2,rbd2,rcd2,ang1,ang2,ang3,
     2		ang4,ang5,ang6,der1,der2,der3,der4,der5,der6,
     3		cosine,sine)

	call center2(a,b,ra2,rb2,rab2,c_ab,l_ab)
	call center2(a,c,ra2,rc2,rac2,c_ac,l_ac)
	call center2(a,d,ra2,rd2,rad2,c_ad,l_ad)
	call center2(b,c,rb2,rc2,rbc2,c_bc,l_bc)
	call center2(b,d,rb2,rd2,rbd2,c_bd,l_bd)
	call center2(c,d,rc2,rd2,rcd2,c_cd,l_cd)

	val_ab = l_ab*rab
	val_ac = l_ac*rac
	val_ad = l_ad*rad
	val_bc = l_bc*rbc
	val_bd = l_bd*rbd
	val_cd = l_cd*rcd

	val2_ab = rab - val_ab
	val2_ac = rac - val_ac
	val2_ad = rad - val_ad
	val2_bc = rbc - val_bc
	val2_bd = rbd - val_bd
	val2_cd = rcd - val_cd

	surfa = -0.5d0*ra + ang1*val2_ab + ang2*val2_ac +
     1		ang3*val2_ad

	surfa = twopi*ra*surfa

	dsurfa_ab = ra*(der1(1)*val2_ab+twopi*ang1*l_ab+
     1		der2(1)*val2_ac + der3(1)*val2_ad)
	dsurfa_ac = ra*(der1(2)*val2_ab+der2(2)*val2_ac +
     1		twopi*ang2*l_ac + der3(2)*val2_ad)
	dsurfa_ad = ra*(der1(3)*val2_ab + der2(3)*val2_ac +
     1		der3(3)*val2_ad + twopi*ang3*l_ad)
	dsurfa_bc = ra*(der1(4)*val2_ab + der2(4)*val2_ac +
     1		 der3(4)*val2_ad)
	dsurfa_bd = ra*(der1(5)*val2_ab + der2(5)*val2_ac +
     1		 der3(5)*val2_ad)
	dsurfa_cd = ra*(der1(6)*val2_ab + der2(6)*val2_ac +
     1		 der3(6)*val2_ad)

        dsurfa(1) = dsurfa_ab
        dsurfa(2) = dsurfa_ac
        dsurfa(3) = dsurfa_ad
        dsurfa(4) = dsurfa_bc
        dsurfa(5) = dsurfa_bd
        dsurfa(6) = dsurfa_cd

	surfb = -0.5d0*rb + ang1*val_ab + ang5*val2_bd +
     1		ang4*val2_bc
	surfb = twopi*rb*surfb

	dsurfb_ab = rb*(der1(1)*val_ab+twopi*ang1*(1-l_ab)+
     1		der5(1)*val2_bd +der4(1)*val2_bc)
	dsurfb_ac = rb*(der1(2)*val_ab+der5(2)*val2_bd +
     1		der4(2)*val2_bc)
	dsurfb_ad = rb*(der1(3)*val_ab + der5(3)*val2_bd +
     1		der4(3)*val2_bc)
	dsurfb_bc = rb*(der1(4)*val_ab + der5(4)*val2_bd +
     1		 der4(4)*val2_bc+twopi*ang4*l_bc)
	dsurfb_bd = rb*(der1(5)*val_ab +der5(5)*val2_bd +
     1		 twopi*ang5*l_bd + der4(5)*val2_bc)
	dsurfb_cd = rb*(der1(6)*val_ab + der5(6)*val2_bd +
     1		 der4(6)*val2_bc)

        dsurfb(1) = dsurfb_ab
        dsurfb(2) = dsurfb_ac
        dsurfb(3) = dsurfb_ad
        dsurfb(4) = dsurfb_bc
        dsurfb(5) = dsurfb_bd
        dsurfb(6) = dsurfb_cd

	surfc = -0.5d0*rc + ang2*val_ac + ang4*val_bc +
     1		ang6*val2_cd
	surfc = twopi*rc*surfc

	dsurfc_ab = rc*(der2(1)*val_ac+
     1		der4(1)*val_bc + der6(1)*val2_cd)
	dsurfc_ac = rc*(der2(2)*val_ac+twopi*ang2*(1-l_ac)+
     1		der4(2)*val_bc + der6(2)*val2_cd)
	dsurfc_ad = rc*(der2(3)*val_ac + der4(3)*val_bc +
     1		der6(3)*val2_cd)
	dsurfc_bc = rc*(der2(4)*val_ac + der4(4)*val_bc +
     1		 twopi*ang4*(1-l_bc)+der6(4)*val2_cd)
	dsurfc_bd = rc*(der2(5)*val_ac + der4(5)*val_bc +
     1		 der6(5)*val2_cd)
	dsurfc_cd = rc*(der2(6)*val_ac +der4(6)*val_bc +
     1		 der6(6)*val2_cd + twopi*ang6*l_cd)

        dsurfc(1) = dsurfc_ab
        dsurfc(2) = dsurfc_ac
        dsurfc(3) = dsurfc_ad
        dsurfc(4) = dsurfc_bc
        dsurfc(5) = dsurfc_bd
        dsurfc(6) = dsurfc_cd

	surfd = -0.5d0*rd + ang3*val_ad + ang6*val_cd +
     1		ang5*val_bd
	surfd = twopi*rd*surfd

	dsurfd_ab = rd*(der3(1)*val_ad+
     1		der6(1)*val_cd + der5(1)*val_bd)
	dsurfd_ac = rd*(der3(2)*val_ad+
     1		der6(2)*val_cd + der5(2)*val_bd)
	dsurfd_ad = rd*(der3(3)*val_ad + twopi*ang3*(1-l_ad) +
     1		der6(3)*val_cd + der5(3)*val_bd)
	dsurfd_bc = rd*(der3(4)*val_ad + der6(4)*val_cd +
     1		 der5(4)*val_bd)
	dsurfd_bd = rd*(der3(5)*val_ad + der6(5)*val_cd +
     1		 der5(5)*val_bd+ twopi*ang5*(1-l_bd))
	dsurfd_cd = rd*(der3(6)*val_ad + der6(6)*val_cd +
     1		 twopi*ang6*(1-l_cd) + der5(6)*val_bd)

        dsurfd(1) = dsurfd_ab
        dsurfd(2) = dsurfd_ac
        dsurfd(3) = dsurfd_ad
        dsurfd(4) = dsurfd_bc
        dsurfd(5) = dsurfd_bd
        dsurfd(6) = dsurfd_cd

        return
        end
