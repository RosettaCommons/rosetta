c	volume_tools.f		Version 1 : 4/17/2007
c
c	Copyright (C) 2007 Patrice Koehl
c
c	This file contains a series of subroutines used to compute the
c	volume area of a union of balls, and optionally the derivatives
c	of the volume with respect to the coordinates of the centers of the ball.
c	As a side result, it also provides the surface area (and its derivatives)
c	of the union of balls.
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
c	twosphere_vol.f
c
c	This subroutine calculates the volume and surface of the
c	intersection of two spheres; it is only called when the
c	intersection exists
c	
c	Copyright (C) 2007 Patrice Koehl
c
	subroutine twosphere_vol(a,b,ra,ra2,rb,rb2,rab,rab2,
     1			surfa,surfb,vola,volb)
c
c	Input:
c			a,b	: position of the centers of the 2 spheres
c			rab	: distance between the centers of the 2 spheres
c			rab2	: distance between the centers of the 2 spheres
c				  (squared)
c			ra,rb	: radii of sphere A and B, respectively
c			ra2 and rb2 are the squared of the quantities
c			above)
c			option	: 0 for surf and vol only, 1 if the derivatives
c				  are also computed
c	Output
c			surfa	: partial contribution of A to the total
c				  surface of the intersection
c			surfb	: partial contribution of B to the total
c				  surface of the intersection
c			vola	: partial contribution of A to the total
c				  volume of the intersection
c			volb	: partial contribution of B to the total
c				  volume of the intersection
c
	integer	i
c
	real*8	ra,rb,surfa,surfb,vola,volb
	real*8	vala,valb,lamda
	real*8	ra2,rb2,rab,rab2,ha,hb,sa,ca,sb,cb
	real*8	dera,derb,der1,der2,coef1,coef2
	real*8	coefa,coefb,Aab
	real*8	pi,twopi,precision
c
	real*8	a(3),b(3),c(3)
c
	common /constants/ pi,twopi,precision
c
	save
c
c	Get "center" of the two spheres
c
	call center2(a,b,ra2,rb2,rab2,c,lamda)
c
	valb = lamda*rab
	vala = rab-valb
c
c	Get height of the cap of sphere A occluded by sphere B
c
	ha = ra - vala
c
c	same for sphere B ...
c
	hb = rb - valb
c
c	Get surfaces of intersection
c
	surfa = twopi*ra*ha
	surfb = twopi*rb*hb
c
c	Now get volume
c
	Aab = pi*(ra2-vala*vala)
c
	sa = ra*surfa
	ca = vala*Aab
c
	vola = (sa-ca)/3
c
	sb = rb*surfb
	cb = valb*Aab
c
	volb = (sb-cb)/3
c
	return
	end
c
c	twosphere_dvol_coord.f
c
c	This subroutine calculates the volume and surface of the
c	intersection of two spheres; it is only called when the
c	intersection exists
c	
c	Copyright (C) 2007 Patrice Koehl
c
c	It also provides the derivative of the surface and volume with
c	respect to the coordinates of the center of the atoms
c
	subroutine twosphere_dvol_coord(a,b,ra,ra2,rb,rb2,rab,rab2,
     1			surfa,surfb,vola,volb,dsurfa,dsurfb,
     2			dvola,dvolb)
c
c	Input:
c			a,b	: position of the centers of the 2 spheres
c			rab	: distance between the centers of the 2 spheres
c			rab2	: distance between the centers of the 2 spheres
c				  (squared)
c			ra,rb	: radii of sphere A and B, respectively
c			ra2 and rb2 are the squared of the quantities
c			above)
c			option	: 0 for surf and vol only, 1 if the derivatives
c				  are also computed
c	Output
c			surfa	: partial contribution of A to the total
c				  surface of the intersection
c			surfb	: partial contribution of B to the total
c				  surface of the intersection
c			vola	: partial contribution of A to the total
c				  volume of the intersection
c			volb	: partial contribution of B to the total
c				  volume of the intersection
c                       dsurfa,dsurfb,dvola,dvolb: derivatives of surfa,
c                                                  surfb,vola and volb
c                                                  with respect to the
c                                                  coordinates of A and B
	integer	i
c
	real*8	ra,rb,surfa,surfb,vola,volb
	real*8	vala,valb,lamda
	real*8	ra2,rb2,rab,rab2,ha,hb,sa,ca,sb,cb
	real*8	dera,derb,der1,der2,coef1,coef2
	real*8	coefa,coefb,Aab
	real*8	pi,twopi,precision
c
	real*8	a(3),b(3),c(3)
	real*8  u_ab(3),dsurfa(3,2),dsurfb(3,2)
	real*8  dvola(3,2),dvolb(3,2)
c
	common /constants/ pi,twopi,precision
c
	save
c
c	Get "center" of the two spheres
c
	call center2(a,b,ra2,rb2,rab2,c,lamda)
c
	valb = lamda*rab
	vala = rab-valb
c
c	Get height of the cap of sphere A occluded by sphere B
c
	ha = ra - vala
c
c	same for sphere B ...
c
	hb = rb - valb
c
c	Get surfaces of intersection
c
	surfa = twopi*ra*ha
	surfb = twopi*rb*hb
c
c	Now get volume
c
	Aab = pi*(ra2-vala*vala)
c
	sa = ra*surfa
	ca = vala*Aab
c
	vola = (sa-ca)/3
c
	sb = rb*surfb
	cb = valb*Aab
c
	volb = (sb-cb)/3
c
c       Compute derivatives
c
	do 50 i = 1,3
		u_ab(i) = (a(i)-b(i))/rab
50	continue
c
	dera = -lamda
	derb = lamda-1
c
	der1 = twopi*ra*dera
	der2 = twopi*rb*derb
c
	do 100 i = 1,3
		coef1 = der1*u_ab(i)
		coef2 = der2*u_ab(i)
		dsurfa(i,1) = coef1
		dsurfa(i,2) = -coef1
		dsurfb(i,1) = coef2
		dsurfb(i,2) = -coef2
100	continue
c
	coefa = Aab*lamda
	coefb = Aab - coefa
c
	do 150 i = 1,3
		coef1 = -coefa*u_ab(i)
		coef2 = -coefb*u_ab(i)
		dvola(i,1) = coef1
		dvola(i,2) = -coef1
		dvolb(i,1) = coef2
		dvolb(i,2) = -coef2
150     continue
c
	return
        end
c
c	twosphere_dvol_dist.f
c
c	This subroutine calculates the volume and surface of the
c	intersection of two spheres; it is only called when the
c	intersection exists
c	
c	Copyright (C) 2007 Patrice Koehl
c
c	It also provides the derivative of the surface and volume with
c	respect to the distance Rab
c
	subroutine twosphere_dvol_dist(a,b,ra,ra2,rb,rb2,rab,rab2,
     1			surfa,surfb,vola,volb,dsurfa,dsurfb,
     2			dvola,dvolb)
c
c	Input:
c			a,b	: position of the centers of the 2 spheres
c			rab	: distance between the centers of the 2 spheres
c			rab2	: distance between the centers of the 2 spheres
c				  (squared)
c			ra,rb	: radii of sphere A and B, respectively
c			ra2 and rb2 are the squared of the quantities
c			above)
c			option	: 0 for surf and vol only, 1 if the derivatives
c				  are also computed
c	Output
c			surfa	: partial contribution of A to the total
c				  surface of the intersection
c			surfb	: partial contribution of B to the total
c				  surface of the intersection
c			vola	: partial contribution of A to the total
c				  volume of the intersection
c			volb	: partial contribution of B to the total
c				  volume of the intersection
c                       dsurfa,dsurfb,dvola,dvolb: derivatives of surfa,
c                                                  surfb,vola and volb
c                                                  with respect to Rab
c
	integer	i
c
	real*8	ra,rb,surfa,surfb,vola,volb
	real*8	vala,valb,lamda
	real*8	ra2,rb2,rab,rab2,ha,hb,sa,ca,sb,cb
	real*8	dera,derb,der1,der2,coef1,coef2
	real*8	coefa,coefb,Aab
	real*8	pi,twopi,precision
c
	real*8  dsurfa,dsurfb
	real*8  dvola,dvolb
c
	real*8	a(3),b(3),c(3)
c
	common /constants/ pi,twopi,precision
c
	save
c
c	Get "center" of the two spheres
c
	call center2(a,b,ra2,rb2,rab2,c,lamda)
c
	valb = lamda*rab
	vala = rab-valb
c
c	Get height of the cap of sphere A occluded by sphere B
c
	ha = ra - vala
c
c	same for sphere B ...
c
	hb = rb - valb
c
c	Get surfaces of intersection
c
	surfa = twopi*ra*ha
	surfb = twopi*rb*hb
c
c	Now get volume
c
	Aab = pi*(ra2-vala*vala)
c
	sa = ra*surfa
	ca = vala*Aab
c
	vola = (sa-ca)/3
c
	sb = rb*surfb
	cb = valb*Aab
c
	volb = (sb-cb)/3
c
c       Compute derivatives
c
	dera = -lamda
	derb = lamda-1
c
	dsurfa = twopi*ra*dera
	dsurfb = twopi*rb*derb
c
	dvola = -Aab*lamda
	dvolb = -Aab + Aab*lamda
c
	return
	end
c
c	threesphere_vol.f
c
c	Copyright (C) 2007 Patrice Koehl
c
c	This subroutine calculates the volume and surface of the intersection 
c	of three spheres; it is only called when the intersection exists
c
	subroutine threesphere_vol(a,b,c,ra,rb,rc,ra2,rb2,rc2,
     1		wa,wb,wc,rab,rac,rbc,rab2,rac2,rbc2,
     2		surfa,surfb,surfc,vola,volb,volc)
c
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
c			vola,volb,volc : contribution of A, B and C to
c			the total volume of the intersection of A,B,C
c
	integer	i,option2
c
	real*8	surfa,surfb,surfc
	real*8	vola,volb,volc
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
	real*8  sh_abc,sh_acb,sh_bca
	real*8	s_abc,s_acb,s_bca
	real*8	val1,val2,val3,l1,l2,l3
	real*8	val1b,val2b,val3b
	real*8	rho_ab2,rho_ac2,rho_bc2
	real*8  coef_ab,coef_ac,coef_bc
	real*8	coef,coef1,coef2,coef3
	real*8	pi,twopi,precision
	real*8	a(3),b(3),c(3),center(3),n(3)
	real*8	c_ab(3),c_ac(3),c_bc(3)
	real*8	pabc(3),pacb(3)
	real*8	cosine(3),sine(3)
c
	common /constants/ pi,twopi,precision
c
	save
c
        call center2(a,b,ra2,rb2,rab2,c_ab,l1)
        call center2(a,c,ra2,rc2,rac2,c_ac,l2)
        call center2(b,c,rb2,rc2,rbc2,c_bc,l3)
c
        val1 = l1*rab
        val2 = l2*rac
        val3 = l3*rbc
c
        val1b = rab-val1
        val2b = rac-val2
        val3b = rbc-val3
c
	call center3(a,b,c,wa,wb,wc,center)
c
	call triangle_dual(a,b,c,center,eps,ra2,pabc,pacb,n)
c
	option2 = 1
	call tetra3_noder(a,b,c,pabc,rab,rac,rbc,
     1	ra,rb,rc,seg_ang_abc,seg_ang_acb,
     2	seg_ang_bca,ang_dih_abc,ang_dih_bac,ang_dih_cab,
     2	cosine,sine,option2)
c
	a1 = ra*(1-2*ang_dih_abc)
	a2 = 2*seg_ang_abc*val1b
	a3 = 2*seg_ang_acb*val2b
c
	surfa = twopi*ra*(a1 - a2 - a3)
c
	a1 = rb*(1-2*ang_dih_bac)
	a2 = 2*seg_ang_abc*val1
	a3 = 2*seg_ang_bca*val3b
c
	surfb = twopi*rb*(a1 - a2 - a3)
c
	a1 = rc*(1-2*ang_dih_cab)
	a2 = 2*seg_ang_acb*val2
	a3 = 2*seg_ang_bca*val3
c
	surfc = twopi*rc*(a1 - a2 - a3)
c
	ang_abc = twopi*seg_ang_abc
	ang_acb = twopi*seg_ang_acb
	ang_bca = twopi*seg_ang_bca
c
	cos_abc = cosine(1)
	sin_abc = sine(1)
	cos_acb = cosine(2)
	sin_acb = sine(2)
	cos_bca = cosine(3)
	sin_bca = sine(3)
c
	rho_ab2 = ra2 - val1b*val1b
	rho_ac2 = ra2 - val2b*val2b
	rho_bc2 = rb2 - val3b*val3b
c
	s_abc = rho_ab2*(ang_abc - sin_abc*cos_abc)
	s_acb = rho_ac2*(ang_acb - sin_acb*cos_acb)
	s_bca = rho_bc2*(ang_bca - sin_bca*cos_bca)
c
	s2 = ra*surfa
	c1 = val1b*s_abc
	c2 = val2b*s_acb
c
	vola = (s2 - c1 - c2)/3
c
	s2 = rb*surfb
	c1 = val1*s_abc
	c2 = val3b*s_bca
c
	volb = (s2 - c1 - c2)/3
c
	s2 = rc*surfc
	c1 = val2*s_acb
	c2 = val3*s_bca
c
	volc = (s2 - c1 - c2)/3
c
	return
	end
c
c	threesphere_dvol_coord.f	Version 1 08/25/2000	Patrice Koehl
c
c	This subroutine calculates the volume and surface of the intersection 
c	of three spheres; it is only called when it exists
c	The computation follows the method described in VOLBL by Herbert
c	Edelsbrunner
c
c	This subroutine supposes that each pair of spheres intersect (i.e.
c	a three sphere intersection is possible)
c
	subroutine threesphere_dvol_coord(a,b,c,ra,rb,rc,ra2,rb2,rc2,
     1		wa,wb,wc,rab,rac,rbc,rab2,rac2,rbc2,
     3		surfa,surfb,surfc,vola,volb,volc,dsurfa,dsurfb,dsurfc,
     4		dvola,dvolb,dvolc)
c
c	Input:
c			ra,rb,rc  : radii of sphere A, B and C, respectively
c			ra2,rb2,rc2: radii squared
c			rab,rab2: distance between the centers of sphere A and B
c			rac,rac2: distance between the centers of sphere A and C
c			rbc,rbc2: distance between the centers of sphere B and C
c	Output
c			surfa,surfb,surfc : contribution of A, B and C to
c			the total surface of the intersection of A,B,C
c			vola,volb,volc : contribution of A, B and C to
c			the total surface of the intersection of A,B,C
c			seg_ang_abc,seg_ang_acb,seg_ang_bca: quantities
c			required to compute derivative of surface area
c
	integer	i
c
	real*8	pi,twopi,precision
c
	real*8	surfa,surfb,surfc
	real*8	vola,volb,volc
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
	real*8	s_abc,s_acb,s_bca
	real*8	val1,val2,val3,l1,l2,l3
	real*8	val1b,val2b,val3b
	real*8	rho_ab2,rho_ac2,rho_bc2
	real*8	sin_2abc,sin_2acb,sin_2bca
	real*8	rss_abc,rss_acb,rss_bca
	real*8	der_ab,der_ac,der_bc,diff_ab,diff_ac,diff_bc
	real*8	dsurfa_ab,dsurfa_ac,dsurfa_bc
	real*8	dsurfb_ab,dsurfb_ac,dsurfb_bc
	real*8	dsurfc_ab,dsurfc_ac,dsurfc_bc
	real*8	coef_ab,coef_ac,coef_bc
	real*8	a(3),b(3),c(3),center(3),n(3)
	real*8	c_ab(3),c_ac(3),c_bc(3)
	real*8	pabc(3),pacb(3)
	real*8	der1(3),der2(3),der3(3),der4(3),der5(3),der6(3)
	real*8	dsurfa(3,3),dsurfb(3,3),dsurfc(3,3)
	real*8	dvola(3,3),dvolb(3,3),dvolc(3,3)
	real*8	u_ab(3),u_ac(3),u_bc(3)
	real*8	cosine(3),sine(3)
c
	common /constants/ pi,twopi,precision
c
	save
c
        call center2(a,b,ra2,rb2,rab2,c_ab,l1)
        call center2(a,c,ra2,rc2,rac2,c_ac,l2)
        call center2(b,c,rb2,rc2,rbc2,c_bc,l3)
c
        val1 = l1*rab
        val2 = l2*rac
        val3 = l3*rbc
c
        val1b = rab-val1
        val2b = rac-val2
        val3b = rbc-val3
c
	call center3(a,b,c,wa,wb,wc,center)
c
	call triangle_dual(a,b,c,center,eps,ra2,pabc,pacb,n)
c
	call tetra3(a,b,c,pabc,rab,rac,rbc,rab2,rac2,
     1	rbc2,ra,rb,rc,ra2,rb2,rc2,seg_ang_abc,seg_ang_acb,
     2	seg_ang_bca,ang_dih_abc,ang_dih_bac,ang_dih_cab,
     2	der1,der2,der3,der4,der5,der6,cosine,sine)
c
	a1 = ra*(1-2*ang_dih_abc)
	a2 = 2*seg_ang_abc*val1b
	a3 = 2*seg_ang_acb*val2b
c
	surfa = twopi*ra*(a1 - a2 - a3)
c
	a1 = rb*(1-2*ang_dih_bac)
	a2 = 2*seg_ang_abc*val1
	a3 = 2*seg_ang_bca*val3b
c
	surfb = twopi*rb*(a1 - a2 - a3)
c
	a1 = rc*(1-2*ang_dih_cab)
	a2 = 2*seg_ang_acb*val2
	a3 = 2*seg_ang_bca*val3
c
	surfc = twopi*rc*(a1 - a2 - a3)
c
	do 50 i = 1,3
		u_ab(i) = (a(i)-b(i))/rab
		u_ac(i) = (a(i)-c(i))/rac
		u_bc(i) = (b(i)-c(i))/rbc
50	continue
c
	ang_abc = twopi*seg_ang_abc
	ang_acb = twopi*seg_ang_acb
	ang_bca = twopi*seg_ang_bca
c
	dsurfa_ab = -2*ra*(ra*der4(1) + der1(1)*val1b +
     1		ang_abc*l1 +der2(1)*val2b)
	dsurfa_ac = - 2*ra*(ra*der4(2) +der1(2)*val1b +
     1		der2(2)*val2b + ang_acb*l2)
	dsurfa_bc = -2*ra*(ra*der4(3) +der1(3)*val1b +
     1			der2(3)*val2b)
c
	do 100 i = 1,3
		diff_ab = dsurfa_ab*u_ab(i)
		diff_ac = dsurfa_ac*u_ac(i)
		diff_bc = dsurfa_bc*u_bc(i)
		dsurfa(i,1) = diff_ab + diff_ac
		dsurfa(i,2) = -diff_ab + diff_bc
		dsurfa(i,3) = -diff_ac - diff_bc
100	continue
c
	dsurfb_ab = -2*rb*(rb*der5(1) + der1(1)*val1 +
     1		ang_abc*(1-l1) +der3(1)*val3b)
	dsurfb_ac = - 2*rb*(rb*der5(2) +der1(2)*val1 +
     1		der3(2)*val3b )
	dsurfb_bc = -2*rb*(rb*der5(3) +der1(3)*val1 +
     1			der3(3)*val3b+ang_bca*l3)
c
	do 200 i = 1,3
		diff_ab = dsurfb_ab*u_ab(i)
		diff_ac = dsurfb_ac*u_ac(i)
		diff_bc = dsurfb_bc*u_bc(i)
		dsurfb(i,1) = diff_ab + diff_ac
		dsurfb(i,2) = -diff_ab + diff_bc
		dsurfb(i,3) = -diff_ac - diff_bc
200	continue
c
	dsurfc_ab = -2*rc*(rc*der6(1) + der2(1)*val2 +
     1		der3(1)*val3)
	dsurfc_ac = - 2*rc*(rc*der6(2) +der2(2)*val2 +
     1		ang_acb*(1-l2) + der3(2)*val3 )
	dsurfc_bc = -2*rc*(rc*der6(3) +der2(3)*val2 +
     1			der3(3)*val3+ang_bca*(1-l3))
c
	do 300 i = 1,3
		diff_ab = dsurfc_ab*u_ab(i)
		diff_ac = dsurfc_ac*u_ac(i)
		diff_bc = dsurfc_bc*u_bc(i)
		dsurfc(i,1) = diff_ab + diff_ac
		dsurfc(i,2) = -diff_ab + diff_bc
		dsurfc(i,3) = -diff_ac - diff_bc
300	continue
c
	cos_abc = cosine(1)
	sin_abc = sine(1)
	cos_acb = cosine(2)
	sin_acb = sine(2)
	cos_bca = cosine(3)
	sin_bca = sine(3)
c
	rho_ab2 = ra2 - val1b*val1b
	rho_ac2 = ra2 - val2b*val2b
	rho_bc2 = rb2 - val3b*val3b
c
	s_abc = rho_ab2*(ang_abc - sin_abc*cos_abc)
	s_acb = rho_ac2*(ang_acb - sin_acb*cos_acb)
	s_bca = rho_bc2*(ang_bca - sin_bca*cos_bca)
c
	s2 = ra*surfa
	c1 = val1b*s_abc
	c2 = val2b*s_acb
c
	vola = (s2 - c1 - c2)/3
c
	s2 = rb*surfb
	c1 = val1*s_abc
	c2 = val3b*s_bca
c
	volb = (s2 - c1 - c2)/3
c
	s2 = rc*surfc
	c1 = val2*s_acb
	c2 = val3*s_bca
c
	volc = (s2 - c1 - c2)/3
c
	sin_2abc = 2*sin_abc*cos_abc
	sin_2acb = 2*sin_acb*cos_acb
	sin_2bca = 2*sin_bca*cos_bca
c
	rss_abc  = 2*rho_ab2*sin_abc*sin_abc
	rss_acb  = 2*rho_ac2*sin_acb*sin_acb
	rss_bca  = 2*rho_bc2*sin_bca*sin_bca
c
	der_ab = (ra*dsurfa_ab-s_abc*l1 -val1b*(-2*val1b*
     1	l1*ang_abc + val1b*l1*sin_2abc+ rss_abc*der1(1))
     2	- val2b*rss_acb*der2(1))/3
c
	der_ac = (ra*dsurfa_ac-s_acb*l2-val1b*rss_abc*der1(2)
     1	- val2b*(rss_acb *der2(2)
     2	-2*val2b*l2*ang_acb+val2b*l2*sin_2acb))/3
c
	der_bc = (ra*dsurfa_bc-val1b*rss_abc*der1(3)
     2	- val2b*rss_acb*der2(3))/3
c
	do 400 i = 1,3
		coef_ab = der_ab*u_ab(i)
		coef_ac = der_ac*u_ac(i)
		coef_bc = der_bc*u_bc(i)
		dvola(i,1) = coef_ab + coef_ac
		dvola(i,2) = -coef_ab + coef_bc
		dvola(i,3) = -coef_ac - coef_bc
400	continue
c
	der_ab = (rb*dsurfb_ab-s_abc*(1-l1) -val1*(-2*val1b*
     1	l1*ang_abc + val1b*l1*sin_2abc	+ rss_abc*der1(1))
     2	- val3b*rss_bca*der3(1))/3
c
	der_ac = (rb*dsurfb_ac-val1*rss_abc*der1(2)
     1	- val3b*rss_bca*der3(2))/3
c
	der_bc = (rb*dsurfb_bc-s_bca*l3-val1*rss_abc*der1(3)
     1	- val3b*(rss_bca*der3(3)
     2  -2*val3b*l3*ang_bca+val3b*l3*sin_2bca))/3
c
	do 500 i = 1,3
		coef_ab = der_ab*u_ab(i)
		coef_ac = der_ac*u_ac(i)
		coef_bc = der_bc*u_bc(i)
		dvolb(i,1) = coef_ab + coef_ac
		dvolb(i,2) = -coef_ab + coef_bc
		dvolb(i,3) = -coef_ac - coef_bc
500	continue
c
	der_ab = (rc*dsurfc_ab-val2*rss_acb*der2(1)
     2	- val3*rss_bca*der3(1))/3
c
	der_ac = (rc*dsurfc_ac-s_acb*(1-l2) -val2*(-2*val2b*
     1	l2*ang_acb + val2b*l2*sin_2acb+ rss_acb*der2(2))
     3	- val3*rss_bca*der3(2))/3
	der_bc = (rc*dsurfc_bc-s_bca*(1-l3)-val2*rss_acb*der2(3)
     2	- val3*(rss_bca*der3(3)
     3	-2*val3b*l3*ang_bca+val3b*l3*sin_2bca))/3
c
	do 600 i = 1,3
		coef_ab = der_ab*u_ab(i)
		coef_ac = der_ac*u_ac(i)
		coef_bc = der_bc*u_bc(i)
		dvolc(i,1) = coef_ab + coef_ac
		dvolc(i,2) = -coef_ab + coef_bc
		dvolc(i,3) = -coef_ac - coef_bc
600	continue
c
	return
	end
c
c	threesphere_dvol_dist.f	Version 1 08/25/2000	Patrice Koehl
c
c	This subroutine calculates the volume and surface of the intersection 
c	of three spheres; it is only called when it exists
c	The computation follows the method described in VOLBL by Herbert
c	Edelsbrunner
c
c	This subroutine supposes that each pair of spheres intersect (i.e.
c	a three sphere intersection is possible)
c
	subroutine threesphere_dvol_dist(a,b,c,ra,rb,rc,ra2,rb2,rc2,
     1		wa,wb,wc,rab,rac,rbc,rab2,rac2,rbc2,
     3		surfa,surfb,surfc,vola,volb,volc,dsurfa,dsurfb,dsurfc,
     4		dvola,dvolb,dvolc)
c
c	Input:
c			ra,rb,rc  : radii of sphere A, B and C, respectively
c			ra2,rb2,rc2: radii squared
c			rab,rab2: distance between the centers of sphere A and B
c			rac,rac2: distance between the centers of sphere A and C
c			rbc,rbc2: distance between the centers of sphere B and C
c	Output
c			surfa,surfb,surfc : contribution of A, B and C to
c			the total surface of the intersection of A,B,C
c			vola,volb,volc : contribution of A, B and C to
c			the total surface of the intersection of A,B,C
c			seg_ang_abc,seg_ang_acb,seg_ang_bca: quantities
c			required to compute derivative of surface area
c
	integer	i
c
	real*8	pi,twopi,precision
c
	real*8	surfa,surfb,surfc
	real*8	vola,volb,volc
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
	real*8	s_abc,s_acb,s_bca
	real*8	sh_abc,sh_acb,sh_bca
	real*8	val1,val2,val3,l1,l2,l3
	real*8	val1b,val2b,val3b
	real*8	rho_ab2,rho_ac2,rho_bc2
	real*8	sin_2abc,sin_2acb,sin_2bca
	real*8	rss_abc,rss_acb,rss_bca
	real*8	der_ab,der_ac,der_bc,diff_ab,diff_ac,diff_bc
	real*8	dsurfa_ab,dsurfa_ac,dsurfa_bc
	real*8	dsurfb_ab,dsurfb_ac,dsurfb_bc
	real*8	dsurfc_ab,dsurfc_ac,dsurfc_bc
	real*8	coef_ab,coef_ac,coef_bc
	real*8	a(3),b(3),c(3),center(3),n(3)
	real*8	c_ab(3),c_ac(3),c_bc(3)
	real*8	pabc(3),pacb(3)
	real*8	der1(3),der2(3),der3(3),der4(3),der5(3),der6(3)
	real*8	dsurfa(3),dsurfb(3),dsurfc(3)
	real*8	dvola(3),dvolb(3),dvolc(3)
	real*8	cosine(3),sine(3)
c
	common /constants/ pi,twopi,precision
c
	save
c
        call center2(a,b,ra2,rb2,rab2,c_ab,l1)
        call center2(a,c,ra2,rc2,rac2,c_ac,l2)
        call center2(b,c,rb2,rc2,rbc2,c_bc,l3)
c
        val1 = l1*rab
        val2 = l2*rac
        val3 = l3*rbc
c
        val1b = rab-val1
        val2b = rac-val2
        val3b = rbc-val3
c
	call center3(a,b,c,wa,wb,wc,center)
c
	call triangle_dual(a,b,c,center,eps,ra2,pabc,pacb,n)
c
	call tetra3(a,b,c,pabc,rab,rac,rbc,rab2,rac2,
     1	rbc2,ra,rb,rc,ra2,rb2,rc2,seg_ang_abc,seg_ang_acb,
     2	seg_ang_bca,ang_dih_abc,ang_dih_bac,ang_dih_cab,
     2	der1,der2,der3,der4,der5,der6,cosine,sine)
c
	a1 = ra*(1-2*ang_dih_abc)
	a2 = 2*seg_ang_abc*val1b
	a3 = 2*seg_ang_acb*val2b
c
	surfa = twopi*ra*(a1 - a2 - a3)
c
	a1 = rb*(1-2*ang_dih_bac)
	a2 = 2*seg_ang_abc*val1
	a3 = 2*seg_ang_bca*val3b
c
	surfb = twopi*rb*(a1 - a2 - a3)
c
	a1 = rc*(1-2*ang_dih_cab)
	a2 = 2*seg_ang_acb*val2
	a3 = 2*seg_ang_bca*val3
c
	surfc = twopi*rc*(a1 - a2 - a3)
c
	ang_abc = twopi*seg_ang_abc
	ang_acb = twopi*seg_ang_acb
	ang_bca = twopi*seg_ang_bca
c
	dsurfa_ab = -2*ra*(ra*der4(1) + der1(1)*val1b +
     1		ang_abc*l1 +der2(1)*val2b)
	dsurfa_ac = - 2*ra*(ra*der4(2) +der1(2)*val1b +
     1		der2(2)*val2b + ang_acb*l2)
	dsurfa_bc = -2*ra*(ra*der4(3) +der1(3)*val1b +
     1			der2(3)*val2b)
c
        dsurfa(1) = dsurfa_ab
        dsurfa(2) = dsurfa_ac
        dsurfa(3) = dsurfa_bc
c
	dsurfb_ab = -2*rb*(rb*der5(1) + der1(1)*val1 +
     1		ang_abc*(1-l1) +der3(1)*val3b)
	dsurfb_ac = - 2*rb*(rb*der5(2) +der1(2)*val1 +
     1		der3(2)*val3b )
	dsurfb_bc = -2*rb*(rb*der5(3) +der1(3)*val1 +
     1			der3(3)*val3b+ang_bca*l3)
c
        dsurfb(1) = dsurfb_ab
        dsurfb(2) = dsurfb_ac
        dsurfb(3) = dsurfb_bc
c
	dsurfc_ab = -2*rc*(rc*der6(1) + der2(1)*val2 +
     1		der3(1)*val3)
	dsurfc_ac = - 2*rc*(rc*der6(2) +der2(2)*val2 +
     1		ang_acb*(1-l2) + der3(2)*val3 )
	dsurfc_bc = -2*rc*(rc*der6(3) +der2(3)*val2 +
     1			der3(3)*val3+ang_bca*(1-l3))
c
        dsurfc(1) = dsurfc_ab
        dsurfc(2) = dsurfc_ac
        dsurfc(3) = dsurfc_bc
c
	cos_abc = cosine(1)
	sin_abc = sine(1)
	cos_acb = cosine(2)
	sin_acb = sine(2)
	cos_bca = cosine(3)
	sin_bca = sine(3)
c
	rho_ab2 = ra2 - val1b*val1b
	rho_ac2 = ra2 - val2b*val2b
	rho_bc2 = rb2 - val3b*val3b
c
	s_abc = rho_ab2*(ang_abc - sin_abc*cos_abc)
	s_acb = rho_ac2*(ang_acb - sin_acb*cos_acb)
	s_bca = rho_bc2*(ang_bca - sin_bca*cos_bca)
c
	s2 = ra*surfa
	c1 = val1b*s_abc
	c2 = val2b*s_acb
c
	vola = (s2 - c1 - c2)/3
c
	s2 = rb*surfb
	c1 = val1*s_abc
	c2 = val3b*s_bca
c
	volb = (s2 - c1 - c2)/3
c
	s2 = rc*surfc
	c1 = val2*s_acb
	c2 = val3*s_bca
c
	volc = (s2 - c1 - c2)/3
c
	sin_2abc = 2*sin_abc*cos_abc
	sin_2acb = 2*sin_acb*cos_acb
	sin_2bca = 2*sin_bca*cos_bca
c
	rss_abc  = 2*rho_ab2*sin_abc*sin_abc
	rss_acb  = 2*rho_ac2*sin_acb*sin_acb
	rss_bca  = 2*rho_bc2*sin_bca*sin_bca
c
	der_ab = (ra*dsurfa_ab-s_abc*l1 -val1b*(-2*val1b*
     1	l1*ang_abc + val1b*l1*sin_2abc+ rss_abc*der1(1))
     2	- val2b*rss_acb*der2(1))/3
c
	der_ac = (ra*dsurfa_ac-s_acb*l2-val1b*rss_abc*der1(2)
     1	- val2b*(rss_acb *der2(2)
     2	-2*val2b*l2*ang_acb+val2b*l2*sin_2acb))/3
c
	der_bc = (ra*dsurfa_bc-val1b*rss_abc*der1(3)
     2	- val2b*rss_acb*der2(3))/3
c
        dvola(1) = der_ab
        dvola(2) = der_ac
        dvola(3) = der_bc
c
	der_ab = (rb*dsurfb_ab-s_abc*(1-l1) -val1*(-2*val1b*
     1	l1*ang_abc + val1b*l1*sin_2abc	+ rss_abc*der1(1))
     2	- val3b*rss_bca*der3(1))/3
c
	der_ac = (rb*dsurfb_ac-val1*rss_abc*der1(2)
     1	- val3b*rss_bca*der3(2))/3
c
	der_bc = (rb*dsurfb_bc-s_bca*l3-val1*rss_abc*der1(3)
     1	- val3b*(rss_bca*der3(3)
     2  -2*val3b*l3*ang_bca+val3b*l3*sin_2bca))/3
c
        dvolb(1) = der_ab
        dvolb(2) = der_ac
        dvolb(3) = der_bc
c
	der_ab = (rc*dsurfc_ab-val2*rss_acb*der2(1)
     2	- val3*rss_bca*der3(1))/3
c
	der_ac = (rc*dsurfc_ac-s_acb*(1-l2) -val2*(-2*val2b*
     1	l2*ang_acb + val2b*l2*sin_2acb+ rss_acb*der2(2))
     3	- val3*rss_bca*der3(2))/3
	der_bc = (rc*dsurfc_bc-s_bca*(1-l3)-val2*rss_acb*der2(3)
     2	- val3*(rss_bca*der3(3)
     3	-2*val3b*l3*ang_bca+val3b*l3*sin_2bca))/3
c
        dvolc(1) = der_ab
        dvolc(2) = der_ac
        dvolc(3) = der_bc
c
	return
	end
c
c	Foursphere_vol.f
c
c	This subroutine calculates the volume and surface area of the
c	intersection of four spheres; this intersection is supposed to exist
c
c	This routine assumes that the 4 points (a,b,c,d) are in ccw order
c
c	Copyright (C) 2002 Patrice Koehl
c
	subroutine foursphere_vol(a,b,c,d,ra,rb,rc,rd,
     1	ra2,rb2,rc2,rd2,rab,rac,rad,rbc,rbd,rcd,
     2  rab2,rac2,rad2,rbc2,rbd2,rcd2,wa,wb,wc,wd,
     3  shabc,shacb,shbca,shabd,shadb,shbda,shacd,shadc,
     4  shcda,shbcd,shbdc,shcdb,pacb,pabd,padc,pbcd,eps1,eps3,eps5,eps7,
     5  surfa,surfb,surfc,surfd,vola,volb,volc,vold)
c
	integer	i,option
c
	real*8	pi,twopi,precision
	real*8	ra,rb,rc,rd,ra2,rb2,rc2,rd2
	real*8	rab,rac,rad,rbc,rbd,rcd
	real*8	rab2,rac2,rad2,rbc2,rbd2,rcd2
	real*8	surfa,surfb,surfc,surfd
	real*8	vola,volb,volc,vold
	real*8	wa,wb,wc,wd
	real*8	val_ab,val_ac,val_ad,val_bc,val_bd,val_cd
	real*8	val2_ab,val2_ac,val2_ad,val2_bc,val2_bd,val2_cd
	real*8	dab2,dac2,dad2,dbc2,dbd2,dcd2
	real*8	cap_ab,cap_ac,cap_ad,cap_bc,cap_bd,cap_cd
	real*8	ang1,ang2,ang3,ang4,ang5,ang6
	real*8	l_ab,l_ac,l_ad,l_bc,l_bd,l_cd
	real*8	dist1,dist3,dist5,dist7
	real*8	h1,h3,h5,h7,eps1,eps3,eps5,eps7
	real*8	s1,t1,t2
        real*8  shabc,shacb,shbca,shabd,shadb,shbda,shacd,shadc,shcda
        real*8  shbcd,shbdc,shcdb
c
	real*8	a(3),b(3),c(3),d(3)
	real*8	c_ab(3),c_ac(3),c_ad(3),c_bc(3),c_bd(3),c_cd(3)
	real*8	c_abcd(3)
	real*8	pacb(3),pabd(3),padc(3),pbcd(3)
c
	common /constants/ pi,twopi,precision
c
	save
c
        call tetra_6dihed(a,b,c,d,ang1,ang2,ang4,ang3,ang5,ang6)
c
        call center2(a,b,ra2,rb2,rab2,c_ab,l_ab)
        call center2(a,c,ra2,rc2,rac2,c_ac,l_ac)
        call center2(a,d,ra2,rd2,rad2,c_ad,l_ad)
        call center2(b,c,rb2,rc2,rbc2,c_bc,l_bc)
        call center2(b,d,rb2,rd2,rbd2,c_bd,l_bd)
        call center2(c,d,rc2,rd2,rcd2,c_cd,l_cd)
c
        val_ab = l_ab*rab
        val_ac = l_ac*rac
        val_ad = l_ad*rad
        val_bc = l_bc*rbc
        val_bd = l_bd*rbd
        val_cd = l_cd*rcd
c
        val2_ab = rab - val_ab
        val2_ac = rac - val_ac
        val2_ad = rad - val_ad
        val2_bc = rbc - val_bc
        val2_bd = rbd - val_bd
        val2_cd = rcd - val_cd
c
	surfa = -0.5d0*ra + ang1*val2_ab + ang2*val2_ac +
     1		ang3*val2_ad
	surfa = twopi*ra*surfa
c
	surfb = -0.5d0*rb + ang1*val_ab + ang5*val2_bd +
     1		ang4*val2_bc
	surfb = twopi*rb*surfb
c
	surfc = -0.5d0*rc + ang2*val_ac + ang4*val_bc +
     1		ang6*val2_cd
	surfc = twopi*rc*surfc
c
	surfd = -0.5d0*rd + ang3*val_ad + ang6*val_cd +
     1		ang5*val_bd
	surfd = twopi*rd*surfd
c
c       Now compute volume
c
	call center4(a,b,c,d,wa,wb,wc,wd,c_abcd)
c
        dab2 = ra2 - val2_ab*val2_ab
        dac2 = ra2 - val2_ac*val2_ac
        dad2 = ra2 - val2_ad*val2_ad
        dbc2 = rb2 - val2_bc*val2_bc
        dbd2 = rb2 - val2_bd*val2_bd
        dcd2 = rc2 - val2_cd*val2_cd
c
	dist1 = 0
	dist3 = 0
	dist5 = 0
	dist7 = 0
	do 50 i = 1,3
		dist1 = dist1 + (pacb(i)-c_abcd(i))**2
		dist3 = dist3 + (pabd(i)-c_abcd(i))**2
		dist5 = dist5 + (padc(i)-c_abcd(i))**2
		dist7 = dist7 + (pbcd(i)-c_abcd(i))**2
50	continue
	dist1 = sqrt(dist1)
	dist3 = sqrt(dist3)
	dist5 = sqrt(dist5)
	dist7 = sqrt(dist7)
c
	h1 = dist1-eps1
	h3 = dist3-eps3
	h5 = dist5-eps5
	h7 = dist7-eps7
c
	s1 = -twopi*dab2*ang1
	t1 = shabc*h1
	t2 = shabd*h3
c
	cap_ab = s1 -t1 -t2
c
	s1 = -twopi*dac2*ang2
	t1 = shacd*h5
	t2 = shacb*h1
c
	cap_ac = s1 -t1 -t2
c
	s1 = -twopi*dad2*ang3
	t1 = shadb*h3
	t2 = shadc*h5
c
	cap_ad = s1 - t1 -t2
c
	s1 = -twopi*dbc2*ang4
	t1 = shbca*h1
	t2 = shbcd*h7
c
	cap_bc = s1 - t1 -t2
c
	s1 = -twopi*dbd2*ang5
	t1 = shbdc*h7
	t2 = shbda*h3
c
	cap_bd = s1 - t1 -t2
c
	s1 = -twopi*dcd2*ang6
	t1 = shcda*h5
	t2 = shcdb*h7
c
	cap_cd = s1 - t1 -t2
c
	vola = 2*ra*surfa-val2_ab*cap_ab-val2_ac*cap_ac-val2_ad*cap_ad
	vola = vola/6
c
	volb = 2*rb*surfb - val_ab*cap_ab-val2_bd*cap_bd-val2_bc*cap_bc
	volb = volb /6
c
	volc = 2*rc*surfc -val_ac*cap_ac-val_bc*cap_bc-val2_cd*cap_cd
	volc = volc/6
c
	vold = 2*rd*surfd-val_ad*cap_ad-val_bd*cap_bd-val_cd*cap_cd
	vold = vold/6
c
        return
        end
c
c	Foursphere_dvol_coord.f	Version 1 08/28/2000	Patrice Koehl
c
c	This subroutine calculates the volume and surface area of the
c	intersection of four spheres; this intersection is supposed to exist
c	The computation follows the method described in VOLBL by Herbert
c	Edelsbrunner
c
	subroutine foursphere_dvol_coord(a,b,c,d,ra,rb,rc,rd,
     1			ra2,rb2,rc2,rd2,rab,rac,rad,rbc,rbd,rcd,
     2			rab2,rac2,rad2,rbc2,rbd2,rcd2,wa,wb,wc,wd,
     3			eps1,eps3,eps5,eps7,
     4			shabc,shacb,shbca,shabd,shadb,shbda,
     5			shacd,shadc,shcda,shbcd,shbdc,shcdb,
     4			der_abc,der_acb,der_bca,der_abd,der_adb,der_bda,
     5			der_acd,der_adc,der_cda,der_bcd,der_bdc,der_cdb,
     6			pacb,pabd,padc,pbcd,
     7			surfa,surfb,surfc,surfd,vola,volb,volc,vold,
     8			dsurfa,dsurfb,dsurfc,dsurfd,
     9			dvola,dvolb,dvolc,dvold)
c
	integer	i
c
	real*8	pi,twopi,precision
c
	real*8	ra,rb,rc,rd,ra2,rb2,rc2,rd2
	real*8	rab,rac,rad,rbc,rbd,rcd
	real*8	rab2,rac2,rad2,rbc2,rbd2,rcd2
	real*8	surfa,surfb,surfc,surfd
	real*8	wa,wb,wc,wd
	real*8	val_ab,val_ac,val_ad,val_bc,val_bd,val_cd
	real*8	val2_ab,val2_ac,val2_ad,val2_bc,val2_bd,val2_cd
	real*8	ang1,ang2,ang3,ang4,ang5,ang6
	real*8	l_ab,l_ac,l_ad,l_bc,l_bd,l_cd
	real*8	shabc,shabd,shacd,shbcd,shacb,shadb,shadc,shbdc
	real*8	shbca,shbda,shcda,shcdb
	real*8	dist1,dist3,dist5,dist7
	real*8	h1,h3,h5,h7
	real*8	eps1,eps3,eps5,eps7
	real*8	dab2,dac2,dad2,dbc2,dbd2,dcd2
	real*8	s1,t1,t2
	real*8	vola,volb,volc,vold
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
	real*8	cap_ab,cap_ac,cap_ad,cap_bc,cap_bd,cap_cd
	real*8	dcap_ab_ab,dcap_ab_ac,dcap_ab_ad,dcap_ab_bc,dcap_ab_bd
	real*8	dcap_ab_cd
	real*8	dcap_ac_ab,dcap_ac_ac,dcap_ac_ad,dcap_ac_bc,dcap_ac_bd
	real*8	dcap_ac_cd
	real*8	dcap_ad_ab,dcap_ad_ac,dcap_ad_ad,dcap_ad_bc,dcap_ad_bd
	real*8	dcap_ad_cd
	real*8	dcap_bc_ab,dcap_bc_ac,dcap_bc_ad,dcap_bc_bc,dcap_bc_bd
	real*8	dcap_bc_cd
	real*8	dcap_bd_ab,dcap_bd_ac,dcap_bd_ad,dcap_bd_bc,dcap_bd_bd
	real*8	dcap_bd_cd
	real*8	dcap_cd_ab,dcap_cd_ac,dcap_cd_ad,dcap_cd_bc,dcap_cd_bd
	real*8	dcap_cd_cd
c
	real*8	a(3),b(3),c(3),d(3)
	real*8	c_ab(3),c_ac(3),c_ad(3),c_bc(3),c_bd(3),c_cd(3)
	real*8	c_abcd(3)
	real*8	pacb(3),pabd(3),padc(3),pbcd(3)
	real*8  der1(6),der2(6),der3(6),der4(6),der5(6),der6(6)
	real*8	dsurfa(3,4),dsurfb(3,4),dsurfc(3,4),dsurfd(3,4)
	real*8	der_abc(6),der_abd(6),der_acd(6),der_bcd(6)
	real*8	der_acb(6),der_adb(6),der_adc(6),der_bdc(6)
	real*8	der_bca(6),der_bda(6),der_cda(6),der_cdb(6)
	real*8	der_h1(6),der_h3(6),der_h5(6),der_h7(6)
	real*8	cosine(6),sine(6)
c
	real*8	dvola(3,4),dvolb(3,4),dvolc(3,4),dvold(3,4)
	real*8	u_ab(3),u_ac(3),u_ad(3),u_bc(3),u_bd(3),u_cd(3)
c
	common /constants/ pi,twopi,precision
c
	save
c
	do 30 i = 1,3
		u_ab(i) = (a(i)-b(i))/rab
		u_ac(i) = (a(i)-c(i))/rac
		u_ad(i) = (a(i)-d(i))/rad
		u_bc(i) = (b(i)-c(i))/rbc
		u_bd(i) = (b(i)-d(i))/rbd
		u_cd(i) = (c(i)-d(i))/rcd
30	continue
c
	call tetra6(a,b,c,d,rab,rac,rad,rbc,rbd,rcd,rab2,rac2,
     1		rad2,rbc2,rbd2,rcd2,ang1,ang2,ang3,
     2		ang4,ang5,ang6,der1,der2,der3,der4,der5,der6,
     3		cosine,sine)
c
	call center2(a,b,ra2,rb2,rab2,c_ab,l_ab)
	call center2(a,c,ra2,rc2,rac2,c_ac,l_ac)
	call center2(a,d,ra2,rd2,rad2,c_ad,l_ad)
	call center2(b,c,rb2,rc2,rbc2,c_bc,l_bc)
	call center2(b,d,rb2,rd2,rbd2,c_bd,l_bd)
	call center2(c,d,rc2,rd2,rcd2,c_cd,l_cd)
c
	val_ab = l_ab*rab
	val_ac = l_ac*rac
	val_ad = l_ad*rad
	val_bc = l_bc*rbc
	val_bd = l_bd*rbd
	val_cd = l_cd*rcd
c
	val2_ab = rab - val_ab
	val2_ac = rac - val_ac
	val2_ad = rad - val_ad
	val2_bc = rbc - val_bc
	val2_bd = rbd - val_bd
	val2_cd = rcd - val_cd
c
	surfa = -0.5d0*ra + ang1*val2_ab + ang2*val2_ac +
     1		ang3*val2_ad
c
	surfa = twopi*ra*surfa
c
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
c
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
c
	surfb = -0.5d0*rb + ang1*val_ab + ang5*val2_bd +
     1		ang4*val2_bc
	surfb = twopi*rb*surfb
c
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
c
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
c
	surfc = -0.5d0*rc + ang2*val_ac + ang4*val_bc +
     1		ang6*val2_cd
	surfc = twopi*rc*surfc
c
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
c
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
c
	surfd = -0.5d0*rd + ang3*val_ad + ang6*val_cd +
     1		ang5*val_bd
	surfd = twopi*rd*surfd
c
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
c
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
c
	call center4(a,b,c,d,wa,wb,wc,wd,c_abcd)
c
	dab2 = ra2 - val2_ab*val2_ab
	dac2 = ra2 - val2_ac*val2_ac
	dad2 = ra2 - val2_ad*val2_ad
	dbc2 = rb2 - val2_bc*val2_bc
	dbd2 = rb2 - val2_bd*val2_bd
	dcd2 = rc2 - val2_cd*val2_cd
c
	dist1 = 0
	dist3 = 0
	dist5 = 0
	dist7 = 0
	do 500 i = 1,3
		dist1 = dist1 + (pacb(i)-c_abcd(i))**2
		dist3 = dist3 + (pabd(i)-c_abcd(i))**2
		dist5 = dist5 + (padc(i)-c_abcd(i))**2
		dist7 = dist7 + (pbcd(i)-c_abcd(i))**2
500	continue
	dist1 = sqrt(dist1)
	dist3 = sqrt(dist3)
	dist5 = sqrt(dist5)
	dist7 = sqrt(dist7)
c
	h1 = dist1 - eps1
	h3 = dist3 - eps3
	h5 = dist5 - eps5
	h7 = dist7 - eps7
c
	call deriv_h(h1,h3,h5,h7,shabc,shabd,shbca,shbcd,shadb,
     1		     shadc,der_abc,der_abd,der_bca,
     2		     der_bcd,der_adb,der_adc,der1,
     3		     cosine(1),sine(1),der_h1,der_h3,der_h5,der_h7)
c
	s1 = -twopi*dab2*ang1
	t1 = shabc*h1
	t2 = shabd*h3
c
	cap_ab = s1 -t1 -t2
	dcap_ab_ab = 2*twopi*val2_ab*l_ab*ang1 -dab2*der1(1)
     1		-der_abc(1)*h1 - shabc*der_h1(1)
     2		-der_abd(1)*h3 - shabd*der_h3(1)
	dcap_ab_ac =-dab2*der1(2)-der_abc(2)*h1-shabc*der_h1(2)
     1		-der_abd(2)*h3 - shabd*der_h3(2)
	dcap_ab_ad =-dab2*der1(3)-der_abc(3)*h1-shabc*der_h1(3)
     1		-der_abd(3)*h3 - shabd*der_h3(3)
	dcap_ab_bc =-dab2*der1(4)-der_abc(4)*h1-shabc*der_h1(4)
     1		-der_abd(4)*h3 - shabd*der_h3(4)
	dcap_ab_bd =-dab2*der1(5)-der_abc(5)*h1-shabc*der_h1(5)
     1		-der_abd(5)*h3 - shabd*der_h3(5)
	dcap_ab_cd =-dab2*der1(6)-der_abc(6)*h1-shabc*der_h1(6)
     1		-der_abd(6)*h3 - shabd*der_h3(6)
c
	s1 = -twopi*dac2*ang2
	t1 = shacd*h5
	t2 = shacb*h1
c
	cap_ac = s1 -t1 -t2
	dcap_ac_ab =-dac2*der2(1)-der_acd(1)*h5-shacd*der_h5(1)
     1		-der_acb(1)*h1 - shacb*der_h1(1)
	dcap_ac_ac =-dac2*der2(2)-der_acd(2)*h5-shacd*der_h5(2)
     1		-der_acb(2)*h1 - shacb*der_h1(2)
     2		+2*twopi*val2_ac*l_ac*ang2
	dcap_ac_ad =-dac2*der2(3)-der_acd(3)*h5-shacd*der_h5(3)
     1		-der_acb(3)*h1 - shacb*der_h1(3)
	dcap_ac_bc =-dac2*der2(4)-der_acd(4)*h5-shacd*der_h5(4)
     1		-der_acb(4)*h1 - shacb*der_h1(4)
	dcap_ac_bd =-dac2*der2(5)-der_acd(5)*h5-shacd*der_h5(5)
     1		-der_acb(5)*h1 - shacb*der_h1(5)
	dcap_ac_cd =-dac2*der2(6)-der_acd(6)*h5-shacd*der_h5(6)
     1		-der_acb(6)*h1 - shacb*der_h1(6)
c
	s1 = -twopi*dad2*ang3
	t1 = shadb*h3
	t2 = shadc*h5
c
	cap_ad = s1 - t1 -t2
c
	dcap_ad_ab = -dad2*der3(1)-der_adb(1)*h3-shadb*der_h3(1)
     1			-der_adc(1)*h5 - shadc*der_h5(1)
	dcap_ad_ac = -dad2*der3(2)-der_adb(2)*h3-shadb*der_h3(2)
     1			-der_adc(2)*h5 - shadc*der_h5(2)
	dcap_ad_ad = -dad2*der3(3)+2*val2_ad*l_ad*twopi*ang3
     1	-der_adb(3)*h3-shadb*der_h3(3) -der_adc(3)*h5 - shadc*der_h5(3)
	dcap_ad_bc = -dad2*der3(4)-der_adb(4)*h3-shadb*der_h3(4)
     1			-der_adc(4)*h5 - shadc*der_h5(4)
	dcap_ad_bd = -dad2*der3(5)-der_adb(5)*h3-shadb*der_h3(5)
     1			-der_adc(5)*h5 - shadc*der_h5(5)
	dcap_ad_cd = -dad2*der3(6)-der_adb(6)*h3-shadb*der_h3(6)
     1			-der_adc(6)*h5 - shadc*der_h5(6)
c
	s1 = -twopi*dbc2*ang4
	t1 = shbca*h1
	t2 = shbcd*h7
c
	cap_bc = s1 - t1 -t2
c
	dcap_bc_ab = -dbc2*der4(1)-der_bca(1)*h1-shbca*der_h1(1)
     1			-der_bcd(1)*h7 - shbcd*der_h7(1)
	dcap_bc_ac = -dbc2*der4(2)-der_bca(2)*h1-shbca*der_h1(2)
     1			-der_bcd(2)*h7 - shbcd*der_h7(2)
	dcap_bc_ad = -dbc2*der4(3)-der_bca(3)*h1-shbca*der_h1(3)
     1			-der_bcd(3)*h7 - shbcd*der_h7(3)
	dcap_bc_bc = -dbc2*der4(4)-der_bca(4)*h1-shbca*der_h1(4)
     1			-der_bcd(4)*h7 - shbcd*der_h7(4)
     2			+2*twopi*val2_bc*l_bc*ang4
	dcap_bc_bd = -dbc2*der4(5)-der_bca(5)*h1-shbca*der_h1(5)
     1			-der_bcd(5)*h7 - shbcd*der_h7(5)
	dcap_bc_cd = -dbc2*der4(6)-der_bca(6)*h1-shbca*der_h1(6)
     1			-der_bcd(6)*h7 - shbcd*der_h7(6)
c
	s1 = -twopi*dbd2*ang5
	t1 = shbdc*h7
	t2 = shbda*h3
c
	cap_bd = s1 - t1 -t2
c
	dcap_bd_ab = -dbd2*der5(1)-der_bdc(1)*h7-shbdc*der_h7(1)
     1			-der_bda(1)*h3 - shbda*der_h3(1)
	dcap_bd_ac = -dbd2*der5(2)-der_bdc(2)*h7-shbdc*der_h7(2)
     1			-der_bda(2)*h3 - shbda*der_h3(2)
	dcap_bd_ad = -dbd2*der5(3)-der_bdc(3)*h7-shbdc*der_h7(3)
     1			-der_bda(3)*h3 - shbda*der_h3(3)
	dcap_bd_bc = -dbd2*der5(4)-der_bdc(4)*h7-shbdc*der_h7(4)
     1			-der_bda(4)*h3 - shbda*der_h3(4)
	dcap_bd_bd = -dbd2*der5(5)-der_bdc(5)*h7-shbdc*der_h7(5)
     1			-der_bda(5)*h3 - shbda*der_h3(5)
     2			+2*twopi*val2_bd*l_bd*ang5
	dcap_bd_cd = -dbd2*der5(6)-der_bdc(6)*h7-shbdc*der_h7(6)
     1			-der_bda(6)*h3 - shbda*der_h3(6)
c
	s1 = -twopi*dcd2*ang6
	t1 = shcda*h5
	t2 = shcdb*h7
c
	cap_cd = s1 - t1 -t2
c
	dcap_cd_ab = -dcd2*der6(1)-der_cda(1)*h5-shcda*der_h5(1)
     1			-der_cdb(1)*h7 - shcdb*der_h7(1)
	dcap_cd_ac = -dcd2*der6(2)-der_cda(2)*h5-shcda*der_h5(2)
     1			-der_cdb(2)*h7 - shcdb*der_h7(2)
	dcap_cd_ad = -dcd2*der6(3)-der_cda(3)*h5-shcda*der_h5(3)
     1			-der_cdb(3)*h7 - shcdb*der_h7(3)
	dcap_cd_bc = -dcd2*der6(4)-der_cda(4)*h5-shcda*der_h5(4)
     1			-der_cdb(4)*h7 - shcdb*der_h7(4)
	dcap_cd_bd = -dcd2*der6(5)-der_cda(5)*h5-shcda*der_h5(5)
     1			-der_cdb(5)*h7 - shcdb*der_h7(5)
	dcap_cd_cd = -dcd2*der6(6)-der_cda(6)*h5-shcda*der_h5(6)
     1			-der_cdb(6)*h7 - shcdb*der_h7(6)
     2			+2*twopi*val2_cd*l_cd*ang6
c
	vola = 2*ra*surfa-val2_ab*cap_ab-val2_ac*cap_ac-val2_ad*cap_ad
c
	vola = vola/6
c
	der_ab = (2*ra*dsurfa_ab-l_ab*cap_ab-val2_ab*dcap_ab_ab
     1	-val2_ac*dcap_ac_ab -val2_ad*dcap_ad_ab)/6
	der_ac = (2*ra*dsurfa_ac-val2_ab*dcap_ab_ac-l_ac*cap_ac
     1	-val2_ac*dcap_ac_ac - val2_ad*dcap_ad_ac)/6
	der_ad = (2*ra*dsurfa_ad-val2_ab*dcap_ab_ad-val2_ac*dcap_ac_ad
     1	-l_ad*cap_ad - val2_ad*dcap_ad_ad)/6
	der_bc = (2*ra*dsurfa_bc-val2_ab*dcap_ab_bc-val2_ac*dcap_ac_bc
     1	-val2_ad*dcap_ad_bc)/6
	der_bd = (2*ra*dsurfa_bd-val2_ab*dcap_ab_bd-val2_ac*dcap_ac_bd
     1	-val2_ad*dcap_ad_bd)/6
	der_cd = (2*ra*dsurfa_cd-val2_ab*dcap_ab_cd-val2_ac*dcap_ac_cd
     1	-val2_ad*dcap_ad_cd)/6
c
	do 600 i = 1,3
		diff_ab = der_ab*u_ab(i)
		diff_ac = der_ac*u_ac(i)
		diff_ad = der_ad*u_ad(i)
		diff_bc = der_bc*u_bc(i)
		diff_bd = der_bd*u_bd(i)
		diff_cd = der_cd*u_cd(i)
		dvola(i,1) = diff_ab+diff_ac+ diff_ad
		dvola(i,2) = -diff_ab+diff_bc+diff_bd
		dvola(i,3) = -diff_ac-diff_bc+diff_cd
		dvola(i,4) = -diff_ad-diff_bd-diff_cd
600	continue
c
	volb = 2*rb*surfb - val_ab*cap_ab-val2_bd*cap_bd-val2_bc*cap_bc
	volb = volb /6
c
	der_ab = (2*rb*dsurfb_ab-val_ab*dcap_ab_ab-val2_bd*dcap_bd_ab
     1		-val2_bc*dcap_bc_ab-(1-l_ab)*cap_ab)/6
	der_ac = (2*rb*dsurfb_ac-val_ab*dcap_ab_ac-val2_bd*dcap_bd_ac
     1		-val2_bc*dcap_bc_ac)/6
	der_ad = (2*rb*dsurfb_ad-val_ab*dcap_ab_ad-val2_bd*dcap_bd_ad
     1		-val2_bc*dcap_bc_ad)/6
	der_bc = (2*rb*dsurfb_bc-val_ab*dcap_ab_bc-val2_bd*dcap_bd_bc
     1		-val2_bc*dcap_bc_bc - l_bc*cap_bc)/6
	der_bd = (2*rb*dsurfb_bd-val_ab*dcap_ab_bd-val2_bd*dcap_bd_bd
     1		-val2_bc*dcap_bc_bd - l_bd*cap_bd)/6
	der_cd = (2*rb*dsurfb_cd-val_ab*dcap_ab_cd-val2_bd*dcap_bd_cd
     1		-val2_bc*dcap_bc_cd)/6
c
	do 700 i = 1,3
		diff_ab = der_ab*u_ab(i)
		diff_ac = der_ac*u_ac(i)
		diff_ad = der_ad*u_ad(i)
		diff_bc = der_bc*u_bc(i)
		diff_bd = der_bd*u_bd(i)
		diff_cd = der_cd*u_cd(i)
		dvolb(i,1) = diff_ab+diff_ac+ diff_ad
		dvolb(i,2) = -diff_ab+diff_bc+diff_bd
		dvolb(i,3) = -diff_ac-diff_bc+diff_cd
		dvolb(i,4) = -diff_ad-diff_bd-diff_cd
700	continue
c
	volc = 2*rc*surfc -val_ac*cap_ac-val_bc*cap_bc-val2_cd*cap_cd
	volc = volc/6
c
	der_ab = (2*rc*dsurfc_ab-val_ac*dcap_ac_ab-val_bc*dcap_bc_ab
     1		-val2_cd*dcap_cd_ab)/6
	der_ac = (2*rc*dsurfc_ac-val_ac*dcap_ac_ac-val_bc*dcap_bc_ac
     1		-val2_cd*dcap_cd_ac-(1-l_ac)*cap_ac)/6
	der_ad = (2*rc*dsurfc_ad-val_ac*dcap_ac_ad-val_bc*dcap_bc_ad
     1		-val2_cd*dcap_cd_ad)/6
	der_bc = (2*rc*dsurfc_bc-val_ac*dcap_ac_bc-val_bc*dcap_bc_bc
     1		-val2_cd*dcap_cd_bc-(1-l_bc)*cap_bc)/6
	der_bd = (2*rc*dsurfc_bd-val_ac*dcap_ac_bd-val_bc*dcap_bc_bd
     1		-val2_cd*dcap_cd_bd)/6
	der_cd = (2*rc*dsurfc_cd-val_ac*dcap_ac_cd-val_bc*dcap_bc_cd
     1		-val2_cd*dcap_cd_cd-l_cd*cap_cd)/6
c
	do 800 i = 1,3
		diff_ab = der_ab*u_ab(i)
		diff_ac = der_ac*u_ac(i)
		diff_ad = der_ad*u_ad(i)
		diff_bc = der_bc*u_bc(i)
		diff_bd = der_bd*u_bd(i)
		diff_cd = der_cd*u_cd(i)
		dvolc(i,1) = diff_ab+diff_ac+ diff_ad
		dvolc(i,2) = -diff_ab+diff_bc+diff_bd
		dvolc(i,3) = -diff_ac-diff_bc+diff_cd
		dvolc(i,4) = -diff_ad-diff_bd-diff_cd
800	continue
c
	vold = 2*rd*surfd-val_ad*cap_ad-val_bd*cap_bd-val_cd*cap_cd
	vold = vold/6
c
	der_ab = (2*rd*dsurfd_ab-val_ad*dcap_ad_ab-val_bd*dcap_bd_ab
     1		-val_cd*dcap_cd_ab)/6
	der_ac = (2*rd*dsurfd_ac-val_ad*dcap_ad_ac-val_bd*dcap_bd_ac
     1		-val_cd*dcap_cd_ac)/6
	der_ad = (2*rd*dsurfd_ad-val_ad*dcap_ad_ad-val_bd*dcap_bd_ad
     1		-val_cd*dcap_cd_ad-(1-l_ad)*cap_ad)/6
	der_bc = (2*rd*dsurfd_bc-val_ad*dcap_ad_bc-val_bd*dcap_bd_bc
     1		-val_cd*dcap_cd_bc)/6
	der_bd = (2*rd*dsurfd_bd-val_ad*dcap_ad_bd-val_bd*dcap_bd_bd
     1		-val_cd*dcap_cd_bd-(1-l_bd)*cap_bd)/6
	der_cd = (2*rd*dsurfd_cd-val_ad*dcap_ad_cd-val_bd*dcap_bd_cd
     1		-val_cd*dcap_cd_cd-(1-l_cd)*cap_cd)/6
c
	do 900 i = 1,3
		diff_ab = der_ab*u_ab(i)
		diff_ac = der_ac*u_ac(i)
		diff_ad = der_ad*u_ad(i)
		diff_bc = der_bc*u_bc(i)
		diff_bd = der_bd*u_bd(i)
		diff_cd = der_cd*u_cd(i)
		dvold(i,1) = diff_ab+diff_ac+ diff_ad
		dvold(i,2) = -diff_ab+diff_bc+diff_bd
		dvold(i,3) = -diff_ac-diff_bc+diff_cd
		dvold(i,4) = -diff_ad-diff_bd-diff_cd
900	continue
c
	return
	end
c
c	Foursphere_dvol_dist.f	Version 1 08/28/2000	Patrice Koehl
c
c	This subroutine calculates the volume and surface area of the
c	intersection of four spheres; this intersection is supposed to exist
c	The computation follows the method described in VOLBL by Herbert
c	Edelsbrunner
c
	subroutine foursphere_dvol_dist(a,b,c,d,ra,rb,rc,rd,
     1			ra2,rb2,rc2,rd2,rab,rac,rad,rbc,rbd,rcd,
     2			rab2,rac2,rad2,rbc2,rbd2,rcd2,wa,wb,wc,wd,
     3			eps1,eps3,eps5,eps7,
     4			shabc,shacb,shbca,shabd,shadb,shbda,
     5			shacd,shadc,shcda,shbcd,shbdc,shcdb,
     4			der_abc,der_acb,der_bca,der_abd,der_adb,der_bda,
     5			der_acd,der_adc,der_cda,der_bcd,der_bdc,der_cdb,
     6			pacb,pabd,padc,pbcd,
     7			surfa,surfb,surfc,surfd,vola,volb,volc,vold,
     8			dsurfa,dsurfb,dsurfc,dsurfd,
     9			dvola,dvolb,dvolc,dvold)
c
	integer	i
c
	real*8	pi,twopi,precision
c
	real*8	ra,rb,rc,rd,ra2,rb2,rc2,rd2
	real*8	rab,rac,rad,rbc,rbd,rcd
	real*8	rab2,rac2,rad2,rbc2,rbd2,rcd2
	real*8	surfa,surfb,surfc,surfd
	real*8	wa,wb,wc,wd
	real*8	val_ab,val_ac,val_ad,val_bc,val_bd,val_cd
	real*8	val2_ab,val2_ac,val2_ad,val2_bc,val2_bd,val2_cd
	real*8	ang1,ang2,ang3,ang4,ang5,ang6
	real*8	l_ab,l_ac,l_ad,l_bc,l_bd,l_cd
	real*8	shabc,shabd,shacd,shbcd,shacb,shadb,shadc,shbdc
	real*8	shbca,shbda,shcda,shcdb
	real*8	dist1,dist3,dist5,dist7
	real*8	h1,h3,h5,h7
	real*8	eps1,eps3,eps5,eps7
	real*8	dab2,dac2,dad2,dbc2,dbd2,dcd2
	real*8	s1,t1,t2
	real*8	vola,volb,volc,vold
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
	real*8	cap_ab,cap_ac,cap_ad,cap_bc,cap_bd,cap_cd
	real*8	dcap_ab_ab,dcap_ab_ac,dcap_ab_ad,dcap_ab_bc,dcap_ab_bd
	real*8	dcap_ab_cd
	real*8	dcap_ac_ab,dcap_ac_ac,dcap_ac_ad,dcap_ac_bc,dcap_ac_bd
	real*8	dcap_ac_cd
	real*8	dcap_ad_ab,dcap_ad_ac,dcap_ad_ad,dcap_ad_bc,dcap_ad_bd
	real*8	dcap_ad_cd
	real*8	dcap_bc_ab,dcap_bc_ac,dcap_bc_ad,dcap_bc_bc,dcap_bc_bd
	real*8	dcap_bc_cd
	real*8	dcap_bd_ab,dcap_bd_ac,dcap_bd_ad,dcap_bd_bc,dcap_bd_bd
	real*8	dcap_bd_cd
	real*8	dcap_cd_ab,dcap_cd_ac,dcap_cd_ad,dcap_cd_bc,dcap_cd_bd
	real*8	dcap_cd_cd
c
	real*8	a(3),b(3),c(3),d(3)
	real*8	c_ab(3),c_ac(3),c_ad(3),c_bc(3),c_bd(3),c_cd(3)
c	real*8	c_abc(3),c_abd(3),c_acd(3),c_bcd(3)
	real*8	c_abcd(3)
	real*8	pacb(3),pabd(3),padc(3),pbcd(3)
	real*8  der1(6),der2(6),der3(6),der4(6),der5(6),der6(6)
	real*8	der_abc(6),der_abd(6),der_acd(6),der_bcd(6)
	real*8	der_acb(6),der_adb(6),der_adc(6),der_bdc(6)
	real*8	der_bca(6),der_bda(6),der_cda(6),der_cdb(6)
	real*8	der_h1(6),der_h3(6),der_h5(6),der_h7(6)
	real*8	cosine(6),sine(6)
c
	real*8	dsurfa(6),dsurfb(6),dsurfc(6),dsurfd(6)
	real*8	dvola(6),dvolb(6),dvolc(6),dvold(6)
c
	common /constants/ pi,twopi,precision
c
	save
c
	call tetra6(a,b,c,d,rab,rac,rad,rbc,rbd,rcd,rab2,rac2,
     1		rad2,rbc2,rbd2,rcd2,ang1,ang2,ang3,
     2		ang4,ang5,ang6,der1,der2,der3,der4,der5,der6,
     3		cosine,sine)
c
	call center2(a,b,ra2,rb2,rab2,c_ab,l_ab)
	call center2(a,c,ra2,rc2,rac2,c_ac,l_ac)
	call center2(a,d,ra2,rd2,rad2,c_ad,l_ad)
	call center2(b,c,rb2,rc2,rbc2,c_bc,l_bc)
	call center2(b,d,rb2,rd2,rbd2,c_bd,l_bd)
	call center2(c,d,rc2,rd2,rcd2,c_cd,l_cd)
c
	val_ab = l_ab*rab
	val_ac = l_ac*rac
	val_ad = l_ad*rad
	val_bc = l_bc*rbc
	val_bd = l_bd*rbd
	val_cd = l_cd*rcd
c
	val2_ab = rab - val_ab
	val2_ac = rac - val_ac
	val2_ad = rad - val_ad
	val2_bc = rbc - val_bc
	val2_bd = rbd - val_bd
	val2_cd = rcd - val_cd
c
	surfa = -0.5d0*ra + ang1*val2_ab + ang2*val2_ac +
     1		ang3*val2_ad
c
	surfa = twopi*ra*surfa
c
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
c
        dsurfa(1) = dsurfa_ab
        dsurfa(2) = dsurfa_ac
        dsurfa(3) = dsurfa_ad
        dsurfa(4) = dsurfa_bc
        dsurfa(5) = dsurfa_bd
        dsurfa(6) = dsurfa_cd
c
	surfb = -0.5d0*rb + ang1*val_ab + ang5*val2_bd +
     1		ang4*val2_bc
	surfb = twopi*rb*surfb
c
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
c
        dsurfb(1) = dsurfb_ab
        dsurfb(2) = dsurfb_ac
        dsurfb(3) = dsurfb_ad
        dsurfb(4) = dsurfb_bc
        dsurfb(5) = dsurfb_bd
        dsurfb(6) = dsurfb_cd
c
	surfc = -0.5d0*rc + ang2*val_ac + ang4*val_bc +
     1		ang6*val2_cd
	surfc = twopi*rc*surfc
c
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
c
        dsurfc(1) = dsurfc_ab
        dsurfc(2) = dsurfc_ac
        dsurfc(3) = dsurfc_ad
        dsurfc(4) = dsurfc_bc
        dsurfc(5) = dsurfc_bd
        dsurfc(6) = dsurfc_cd
c
	surfd = -0.5d0*rd + ang3*val_ad + ang6*val_cd +
     1		ang5*val_bd
	surfd = twopi*rd*surfd
c
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
c
        dsurfd(1) = dsurfd_ab
        dsurfd(2) = dsurfd_ac
        dsurfd(3) = dsurfd_ad
        dsurfd(4) = dsurfd_bc
        dsurfd(5) = dsurfd_bd
        dsurfd(6) = dsurfd_cd
c
	call center4(a,b,c,d,wa,wb,wc,wd,c_abcd)
c
	dab2 = ra2 - val2_ab*val2_ab
	dac2 = ra2 - val2_ac*val2_ac
	dad2 = ra2 - val2_ad*val2_ad
	dbc2 = rb2 - val2_bc*val2_bc
	dbd2 = rb2 - val2_bd*val2_bd
	dcd2 = rc2 - val2_cd*val2_cd
c
	dist1 = 0
	dist3 = 0
	dist5 = 0
	dist7 = 0
	do 500 i = 1,3
		dist1 = dist1 + (pacb(i)-c_abcd(i))**2
		dist3 = dist3 + (pabd(i)-c_abcd(i))**2
		dist5 = dist5 + (padc(i)-c_abcd(i))**2
		dist7 = dist7 + (pbcd(i)-c_abcd(i))**2
500	continue
	dist1 = sqrt(dist1)
	dist3 = sqrt(dist3)
	dist5 = sqrt(dist5)
	dist7 = sqrt(dist7)
c
	h1 = dist1 - eps1
	h3 = dist3 - eps3
	h5 = dist5 - eps5
	h7 = dist7 - eps7
c
c
	call deriv_h(h1,h3,h5,h7,shabc,shabd,shbca,shbcd,shadb,
     1		     shadc,der_abc,der_abd,der_bca,
     2		     der_bcd,der_adb,der_adc,der1,
     3		     cosine(1),sine(1),der_h1,der_h3,der_h5,der_h7)
c
	s1 = -twopi*dab2*ang1
	t1 = shabc*h1
	t2 = shabd*h3
c
	cap_ab = s1 -t1 -t2
	dcap_ab_ab = 2*twopi*val2_ab*l_ab*ang1 -dab2*der1(1)
     1		-der_abc(1)*h1 - shabc*der_h1(1)
     2		-der_abd(1)*h3 - shabd*der_h3(1)
	dcap_ab_ac =-dab2*der1(2)-der_abc(2)*h1-shabc*der_h1(2)
     1		-der_abd(2)*h3 - shabd*der_h3(2)
	dcap_ab_ad =-dab2*der1(3)-der_abc(3)*h1-shabc*der_h1(3)
     1		-der_abd(3)*h3 - shabd*der_h3(3)
	dcap_ab_bc =-dab2*der1(4)-der_abc(4)*h1-shabc*der_h1(4)
     1		-der_abd(4)*h3 - shabd*der_h3(4)
	dcap_ab_bd =-dab2*der1(5)-der_abc(5)*h1-shabc*der_h1(5)
     1		-der_abd(5)*h3 - shabd*der_h3(5)
	dcap_ab_cd =-dab2*der1(6)-der_abc(6)*h1-shabc*der_h1(6)
     1		-der_abd(6)*h3 - shabd*der_h3(6)
c
	s1 = -twopi*dac2*ang2
	t1 = shacd*h5
	t2 = shacb*h1
c
	cap_ac = s1 -t1 -t2
	dcap_ac_ab =-dac2*der2(1)-der_acd(1)*h5-shacd*der_h5(1)
     1		-der_acb(1)*h1 - shacb*der_h1(1)
	dcap_ac_ac =-dac2*der2(2)-der_acd(2)*h5-shacd*der_h5(2)
     1		-der_acb(2)*h1 - shacb*der_h1(2)
     2		+2*twopi*val2_ac*l_ac*ang2
	dcap_ac_ad =-dac2*der2(3)-der_acd(3)*h5-shacd*der_h5(3)
     1		-der_acb(3)*h1 - shacb*der_h1(3)
	dcap_ac_bc =-dac2*der2(4)-der_acd(4)*h5-shacd*der_h5(4)
     1		-der_acb(4)*h1 - shacb*der_h1(4)
	dcap_ac_bd =-dac2*der2(5)-der_acd(5)*h5-shacd*der_h5(5)
     1		-der_acb(5)*h1 - shacb*der_h1(5)
	dcap_ac_cd =-dac2*der2(6)-der_acd(6)*h5-shacd*der_h5(6)
     1		-der_acb(6)*h1 - shacb*der_h1(6)
c
	s1 = -twopi*dad2*ang3
	t1 = shadb*h3
	t2 = shadc*h5
c
	cap_ad = s1 - t1 -t2
c
	dcap_ad_ab = -dad2*der3(1)-der_adb(1)*h3-shadb*der_h3(1)
     1			-der_adc(1)*h5 - shadc*der_h5(1)
	dcap_ad_ac = -dad2*der3(2)-der_adb(2)*h3-shadb*der_h3(2)
     1			-der_adc(2)*h5 - shadc*der_h5(2)
	dcap_ad_ad = -dad2*der3(3)+2*val2_ad*l_ad*twopi*ang3
     1	-der_adb(3)*h3-shadb*der_h3(3) -der_adc(3)*h5 - shadc*der_h5(3)
	dcap_ad_bc = -dad2*der3(4)-der_adb(4)*h3-shadb*der_h3(4)
     1			-der_adc(4)*h5 - shadc*der_h5(4)
	dcap_ad_bd = -dad2*der3(5)-der_adb(5)*h3-shadb*der_h3(5)
     1			-der_adc(5)*h5 - shadc*der_h5(5)
	dcap_ad_cd = -dad2*der3(6)-der_adb(6)*h3-shadb*der_h3(6)
     1			-der_adc(6)*h5 - shadc*der_h5(6)
c
	s1 = -twopi*dbc2*ang4
	t1 = shbca*h1
	t2 = shbcd*h7
c
	cap_bc = s1 - t1 -t2
c
	dcap_bc_ab = -dbc2*der4(1)-der_bca(1)*h1-shbca*der_h1(1)
     1			-der_bcd(1)*h7 - shbcd*der_h7(1)
	dcap_bc_ac = -dbc2*der4(2)-der_bca(2)*h1-shbca*der_h1(2)
     1			-der_bcd(2)*h7 - shbcd*der_h7(2)
	dcap_bc_ad = -dbc2*der4(3)-der_bca(3)*h1-shbca*der_h1(3)
     1			-der_bcd(3)*h7 - shbcd*der_h7(3)
	dcap_bc_bc = -dbc2*der4(4)-der_bca(4)*h1-shbca*der_h1(4)
     1			-der_bcd(4)*h7 - shbcd*der_h7(4)
     2			+2*twopi*val2_bc*l_bc*ang4
	dcap_bc_bd = -dbc2*der4(5)-der_bca(5)*h1-shbca*der_h1(5)
     1			-der_bcd(5)*h7 - shbcd*der_h7(5)
	dcap_bc_cd = -dbc2*der4(6)-der_bca(6)*h1-shbca*der_h1(6)
     1			-der_bcd(6)*h7 - shbcd*der_h7(6)
c
	s1 = -twopi*dbd2*ang5
	t1 = shbdc*h7
	t2 = shbda*h3
c
	cap_bd = s1 - t1 -t2
c
	dcap_bd_ab = -dbd2*der5(1)-der_bdc(1)*h7-shbdc*der_h7(1)
     1			-der_bda(1)*h3 - shbda*der_h3(1)
	dcap_bd_ac = -dbd2*der5(2)-der_bdc(2)*h7-shbdc*der_h7(2)
     1			-der_bda(2)*h3 - shbda*der_h3(2)
	dcap_bd_ad = -dbd2*der5(3)-der_bdc(3)*h7-shbdc*der_h7(3)
     1			-der_bda(3)*h3 - shbda*der_h3(3)
	dcap_bd_bc = -dbd2*der5(4)-der_bdc(4)*h7-shbdc*der_h7(4)
     1			-der_bda(4)*h3 - shbda*der_h3(4)
	dcap_bd_bd = -dbd2*der5(5)-der_bdc(5)*h7-shbdc*der_h7(5)
     1			-der_bda(5)*h3 - shbda*der_h3(5)
     2			+2*twopi*val2_bd*l_bd*ang5
	dcap_bd_cd = -dbd2*der5(6)-der_bdc(6)*h7-shbdc*der_h7(6)
     1			-der_bda(6)*h3 - shbda*der_h3(6)
c
	s1 = -twopi*dcd2*ang6
	t1 = shcda*h5
	t2 = shcdb*h7
c
	cap_cd = s1 - t1 -t2
c
	dcap_cd_ab = -dcd2*der6(1)-der_cda(1)*h5-shcda*der_h5(1)
     1			-der_cdb(1)*h7 - shcdb*der_h7(1)
	dcap_cd_ac = -dcd2*der6(2)-der_cda(2)*h5-shcda*der_h5(2)
     1			-der_cdb(2)*h7 - shcdb*der_h7(2)
	dcap_cd_ad = -dcd2*der6(3)-der_cda(3)*h5-shcda*der_h5(3)
     1			-der_cdb(3)*h7 - shcdb*der_h7(3)
	dcap_cd_bc = -dcd2*der6(4)-der_cda(4)*h5-shcda*der_h5(4)
     1			-der_cdb(4)*h7 - shcdb*der_h7(4)
	dcap_cd_bd = -dcd2*der6(5)-der_cda(5)*h5-shcda*der_h5(5)
     1			-der_cdb(5)*h7 - shcdb*der_h7(5)
	dcap_cd_cd = -dcd2*der6(6)-der_cda(6)*h5-shcda*der_h5(6)
     1			-der_cdb(6)*h7 - shcdb*der_h7(6)
     2			+2*twopi*val2_cd*l_cd*ang6
c
	vola = 2*ra*surfa-val2_ab*cap_ab-val2_ac*cap_ac-val2_ad*cap_ad
c
	vola = vola/6
c
	der_ab = (2*ra*dsurfa_ab-l_ab*cap_ab-val2_ab*dcap_ab_ab
     1	-val2_ac*dcap_ac_ab -val2_ad*dcap_ad_ab)/6
	der_ac = (2*ra*dsurfa_ac-val2_ab*dcap_ab_ac-l_ac*cap_ac
     1	-val2_ac*dcap_ac_ac - val2_ad*dcap_ad_ac)/6
	der_ad = (2*ra*dsurfa_ad-val2_ab*dcap_ab_ad-val2_ac*dcap_ac_ad
     1	-l_ad*cap_ad - val2_ad*dcap_ad_ad)/6
	der_bc = (2*ra*dsurfa_bc-val2_ab*dcap_ab_bc-val2_ac*dcap_ac_bc
     1	-val2_ad*dcap_ad_bc)/6
	der_bd = (2*ra*dsurfa_bd-val2_ab*dcap_ab_bd-val2_ac*dcap_ac_bd
     1	-val2_ad*dcap_ad_bd)/6
	der_cd = (2*ra*dsurfa_cd-val2_ab*dcap_ab_cd-val2_ac*dcap_ac_cd
     1	-val2_ad*dcap_ad_cd)/6
c
        dvola(1) = der_ab
        dvola(2) = der_ac
        dvola(3) = der_ad
        dvola(4) = der_bc
        dvola(5) = der_bd
        dvola(6) = der_cd
c
	volb = 2*rb*surfb - val_ab*cap_ab-val2_bd*cap_bd-val2_bc*cap_bc
	volb = volb /6
c
	der_ab = (2*rb*dsurfb_ab-val_ab*dcap_ab_ab-val2_bd*dcap_bd_ab
     1		-val2_bc*dcap_bc_ab-(1-l_ab)*cap_ab)/6
	der_ac = (2*rb*dsurfb_ac-val_ab*dcap_ab_ac-val2_bd*dcap_bd_ac
     1		-val2_bc*dcap_bc_ac)/6
	der_ad = (2*rb*dsurfb_ad-val_ab*dcap_ab_ad-val2_bd*dcap_bd_ad
     1		-val2_bc*dcap_bc_ad)/6
	der_bc = (2*rb*dsurfb_bc-val_ab*dcap_ab_bc-val2_bd*dcap_bd_bc
     1		-val2_bc*dcap_bc_bc - l_bc*cap_bc)/6
	der_bd = (2*rb*dsurfb_bd-val_ab*dcap_ab_bd-val2_bd*dcap_bd_bd
     1		-val2_bc*dcap_bc_bd - l_bd*cap_bd)/6
	der_cd = (2*rb*dsurfb_cd-val_ab*dcap_ab_cd-val2_bd*dcap_bd_cd
     1		-val2_bc*dcap_bc_cd)/6
c
        dvolb(1) = der_ab
        dvolb(2) = der_ac
        dvolb(3) = der_ad
        dvolb(4) = der_bc
        dvolb(5) = der_bd
        dvolb(6) = der_cd
c
	volc = 2*rc*surfc -val_ac*cap_ac-val_bc*cap_bc-val2_cd*cap_cd
	volc = volc/6
c
	der_ab = (2*rc*dsurfc_ab-val_ac*dcap_ac_ab-val_bc*dcap_bc_ab
     1		-val2_cd*dcap_cd_ab)/6
	der_ac = (2*rc*dsurfc_ac-val_ac*dcap_ac_ac-val_bc*dcap_bc_ac
     1		-val2_cd*dcap_cd_ac-(1-l_ac)*cap_ac)/6
	der_ad = (2*rc*dsurfc_ad-val_ac*dcap_ac_ad-val_bc*dcap_bc_ad
     1		-val2_cd*dcap_cd_ad)/6
	der_bc = (2*rc*dsurfc_bc-val_ac*dcap_ac_bc-val_bc*dcap_bc_bc
     1		-val2_cd*dcap_cd_bc-(1-l_bc)*cap_bc)/6
	der_bd = (2*rc*dsurfc_bd-val_ac*dcap_ac_bd-val_bc*dcap_bc_bd
     1		-val2_cd*dcap_cd_bd)/6
	der_cd = (2*rc*dsurfc_cd-val_ac*dcap_ac_cd-val_bc*dcap_bc_cd
     1		-val2_cd*dcap_cd_cd-l_cd*cap_cd)/6
c
        dvolc(1) = der_ab
        dvolc(2) = der_ac
        dvolc(3) = der_ad
        dvolc(4) = der_bc
        dvolc(5) = der_bd
        dvolc(6) = der_cd
c
	vold = 2*rd*surfd-val_ad*cap_ad-val_bd*cap_bd-val_cd*cap_cd
	vold = vold/6
c
	der_ab = (2*rd*dsurfd_ab-val_ad*dcap_ad_ab-val_bd*dcap_bd_ab
     1		-val_cd*dcap_cd_ab)/6
	der_ac = (2*rd*dsurfd_ac-val_ad*dcap_ad_ac-val_bd*dcap_bd_ac
     1		-val_cd*dcap_cd_ac)/6
	der_ad = (2*rd*dsurfd_ad-val_ad*dcap_ad_ad-val_bd*dcap_bd_ad
     1		-val_cd*dcap_cd_ad-(1-l_ad)*cap_ad)/6
	der_bc = (2*rd*dsurfd_bc-val_ad*dcap_ad_bc-val_bd*dcap_bd_bc
     1		-val_cd*dcap_cd_bc)/6
	der_bd = (2*rd*dsurfd_bd-val_ad*dcap_ad_bd-val_bd*dcap_bd_bd
     1		-val_cd*dcap_cd_bd-(1-l_bd)*cap_bd)/6
	der_cd = (2*rd*dsurfd_cd-val_ad*dcap_ad_cd-val_bd*dcap_bd_cd
     1		-val_cd*dcap_cd_cd-(1-l_cd)*cap_cd)/6
c
        dvold(1) = der_ab
        dvold(2) = der_ac
        dvold(3) = der_ad
        dvold(4) = der_bc
        dvold(5) = der_bd
        dvold(6) = der_cd
c
	return
	end
c
c	make_shder.f		Version 1 10/25/2002	Patrice Koehl
c
c	This subroutine prepares all derivatives of the segment heights
c	for volume derivative calculation
c
	subroutine make_shder(itype,abc,acb,bca,abd,adb,bda,
     1          acd,adc,cda,bcd,bdc,cdb,
     2		der_abc,der_acb,der_bca,der_abd,der_adb,der_bda,
     3		der_acd,der_adc,der_cda,der_bcd,der_bdc,der_cdb)
c
	integer	i,itype
	integer	trig1,trig2,trig3,trig4
	integer	ipos(6)
c
        real*8  abc(3),acb(3),bca(3),abd(3),adb(3),bda(3)
        real*8  acd(3),adc(3),cda(3),bcd(3),bdc(3),cdb(3)
	real*8  der_abc(6),der_acb(6),der_bca(6),der_abd(6)
	real*8  der_adb(6),der_bda(6)
	real*8  der_acd(6),der_adc(6),der_cda(6),der_bcd(6)
	real*8  der_bdc(6),der_cdb(6)
c		
c	The tetrahedron is (abcd). The four triangles trig1, trig2, trig3
c	and trig4 are the triangles (abc), (abd), (acd) and (bcd).
c
	do 50 i = 1,6
		ipos(i) = i
50	continue
c
	if(itype.eq.2) then
		ipos(2) = 3
		ipos(3) = 2
		ipos(4) = 5
		ipos(5) = 4
	endif
c
	do 100 i = 1,6
		der_abc(i) = 0
		der_acb(i) = 0
		der_abd(i) = 0
		der_adb(i) = 0
		der_acd(i) = 0
		der_adc(i) = 0
		der_bcd(i) = 0
		der_bdc(i) = 0
		der_bca(i) = 0
		der_bda(i) = 0
		der_cda(i) = 0
		der_cdb(i) = 0
100	continue
c
c	We first consider trig1=(abc). 
c	trig_der(i,trig1) for i = 1,3 is : der_shabc/drab, der_sh_abc/drac,
c					   der_shabc/drbc
c	trig_der(i,trig1) for i = 4,6 is : der_shacb/drab, der_sh_acb/drac,
c					   der_shacb/drbc
c	trig_der(i,trig1) for i = 7,9 is : der_shbca/drab, der_sh_bca/drac,
c					   der_shbca/drbc
c	All these values are now stored in the arrays der_***, with some
c	rearrangements (since the der_ are given with respect to the 6 distances
c	ab, ac, ad, bc, bd and cd
c
	der_abc(ipos(1)) = abc(1)
	der_abc(ipos(2)) = abc(2)
	der_abc(ipos(4)) = abc(3)
c
	der_acb(ipos(1)) = acb(1)
	der_acb(ipos(2)) = acb(2)
	der_acb(ipos(4)) = acb(3)
c
	der_bca(ipos(1)) = bca(1)
	der_bca(ipos(2)) = bca(2)
	der_bca(ipos(4)) = bca(3)
c
c	We now consider trig2=(ipos(abd)). 
c	trig_der(ipos(i,trig2)) for i = 1,3 is : der_shabd/drab, der_sh_abd/drad,
c					   der_shabd/drbd
c	trig_der(ipos(i,trig2)) for i = 4,6 is : der_shadb/drab, der_sh_adb/drad,
c					   der_shadb/drbd
c	trig_der(ipos(i,trig2)) for i = 7,9 is : der_shbda/drab, der_sh_bda/drad,
c					   der_shbda/drbd
c
	der_abd(ipos(1)) = abd(1)
	der_abd(ipos(3)) = abd(2)
	der_abd(ipos(5)) = abd(3)
c
	der_adb(ipos(1)) = adb(1)
	der_adb(ipos(3)) = adb(2)
	der_adb(ipos(5)) = adb(3)
c
	der_bda(ipos(1)) = bda(1)
	der_bda(ipos(3)) = bda(2)
	der_bda(ipos(5)) = bda(3)
c
c	We now consider trig3=(ipos(acd)). 
c	trig_der(ipos(i,trig3)) for i = 1,3 is : der_shacd/drac, der_sh_acd/drad,
c					   der_shacd/drcd
c	trig_der(ipos(i,trig3)) for i = 4,6 is : der_shadc/drac, der_sh_adc/drad,
c					   der_shadc/drcd
c	trig_der(ipos(i,trig3)) for i = 7,9 is : der_shcda/drac, der_sh_cda/drad,
c					   der_shcda/drcd
c
	der_acd(ipos(2)) = acd(1)
	der_acd(ipos(3)) = acd(2)
	der_acd(ipos(6)) = acd(3)
c
	der_adc(ipos(2)) = adc(1)
	der_adc(ipos(3)) = adc(2)
	der_adc(ipos(6)) = adc(3)
c
	der_cda(ipos(2)) = cda(1)
	der_cda(ipos(3)) = cda(2)
	der_cda(ipos(6)) = cda(3)
c
c	We now consider trig4=(ipos(bcd)). 
c	trig_der(ipos(i,trig4)) for i = 1,3 is : der_shbcd/drbc, der_sh_bcd/drbd,
c					   der_shbcd/drcd
c	trig_der(ipos(i,trig4)) for i = 4,6 is : der_shbdc/drbc, der_sh_bdc/drbd,
c					   der_shbdc/drcd
c	trig_der(ipos(i,trig4)) for i = 7,9 is : der_shcdb/drbc, der_sh_cdb/drbd,
c					   der_shcdb/drcd
c
	der_bcd(ipos(4)) = bcd(1)
	der_bcd(ipos(5)) = bcd(2)
	der_bcd(ipos(6)) = bcd(3)
c
	der_bdc(ipos(4)) = bdc(1)
	der_bdc(ipos(5)) = bdc(2)
	der_bdc(ipos(6)) = bdc(3)
c
	der_cdb(ipos(4)) = cdb(1)
	der_cdb(ipos(5)) = cdb(2)
	der_cdb(ipos(6)) = cdb(3)
c
	return
	end
c
c	Deriv_h		Version 1 10/25/2002	Patrice Koehl
c
c	This subroutine evaluates the derivatives of four parameters
c	needed for computing the volume of the intersection of 4 balls
c
	subroutine deriv_h(h1,h3,h5,h7,sh_abc,sh_abd,sh_bca,
     1			sh_bcd,sh_adb,sh_adc,
     2			der_abc,der_abd,der_bca,der_bcd,der_adb,
     3			der_adc,der_ang1,
     4			cos_ang1,sin_ang1,der_h1,der_h3,der_h5,der_h7)
c
	integer i
c
	real*8	pi,twopi,precision
c
	real*8	a,b,c,d,e,f,coef1,coef2,coef3,den
	real*8	h1,h3,h5,h7,cos_ang1,sin_ang1,cotg_ang1
	real*8	sh_abc,sh_abd,sh_bca,sh_bcd,sh_adb,sh_adc
	real*8	der_abc(6),der_abd(6),der_bca(6),der_bcd(6),der_adb(6)
	real*8	der_adc(6),der_ang1(6)
	real*8	der_h1(6),der_h3(6),der_h5(6),der_h7(6)
c
	common /constants/ pi,twopi,precision
c
	save
c
c	To get the derivatives of h1 (= eps1-dist(pabc-pabcd)),
c	h3 (=eps3-dist(dabd-pabcd)), h5 (=eps5-dist(dacd-pabcd))
c	and h7 (=eps7-dist(dbcd-pabcd))
c	we use the four relations:
c
c	sh_abc**2 + h1**2 = sh_abd**2 + h3**2	(1)
c
c	sh_bca**2 + h1**2 = sh_bcd**2 + h7**2	(2)
c
c	sh_adb**2 + h3**2 = sh_adc**2 + h5**2   (3)
c
c	tan(ang1) = (sh_abd*h1 - sh_abc*h3)/(sh_abc*sh_abd-h1*h3) (4)
c
c	We first derive equations (1) and (4). This leads to a 2x2
c	linear system, whose solution are the der_h1 and der_h3 :
c
	cotg_ang1 = - cos_ang1/sin_ang1
c
	a = h3 + sh_abd*cotg_ang1
	b = h1 + sh_abc*cotg_ang1
	d = h1
	e = -h3
	den = a*e-d*b
	coef1 = sh_abd-h3*cotg_ang1
	coef2 = sh_abc-h1*cotg_ang1
	coef3 = -(sh_abc*h3 +h1*sh_abd)/(sin_ang1*sin_ang1)
c
	do 100 i = 1,6
		c = coef1*der_abc(i) + coef2*der_abd(i)+
     1			coef3*der_ang1(i)
		f = -sh_abc*der_abc(i)+sh_abd*der_abd(i)
		der_h1(i) = (c*e-f*b)/den
		der_h3(i) = (a*f-d*c)/den
100	continue
c
c	We use equation 3 to get der_h5:
c
	do 200 i = 1,6
		der_h5(i) = (h3*der_h3(i)+sh_adb*der_adb(i)
     1				-sh_adc*der_adc(i))/h5
200	continue
c
c	We use equation 2 to get der_h7:
c
	do 300 i = 1,6
		der_h7(i) = (h1*der_h1(i)+sh_bca*der_bca(i)
     1				-sh_bcd*der_bcd(i))/h7
300	continue
c
	return
	end
c	threesphere_info.f	Version 1 08/25/2000	Patrice Koehl
c
c	This subroutine pre-computes some quantities specific to 3-sphere
c       intersections
c
c
	subroutine threesphere_info(a,b,c,ra,rb,rc,ra2,rb2,rc2,
     1		wa,wb,wc,rab,rac,rbc,rab2,rac2,rbc2,pabc,pacb,
     2		eps,sh_abc,sh_acb,sh_bca,
     3		der_shabc,der_shacb,der_shbca)
c
c
c	Input:
c			ra,rb,rc  : radii of sphere A, B and C, respectively
c			ra2,rb2,rc2: radii squared
c			rab,rab2: distance between the centers of sphere A and B
c			rac,rac2: distance between the centers of sphere A and C
c			rbc,rbc2: distance between the centers of sphere B and C
c	Output
c
	integer	i
c
	real*8	pi,twopi,precision
c
	real*8	ra,rb,rc,rab,rac,rbc,rab2,rac2,rbc2
	real*8	ra2,rb2,rc2,wa,wb,wc
	real*8	eps
	real*8	seg_ang_abc,seg_ang_acb
	real*8	seg_ang_bca
	real*8	cos_abc,cos_acb,cos_bca
	real*8	sin_abc,sin_acb,sin_bca
	real*8	ang_dih_abc,ang_dih_cab
	real*8	ang_dih_bac
	real*8	sh_abc,sh_acb,sh_bca
	real*8	val1,val2,val3,l1,l2,l3
	real*8	val1b,val2b,val3b
	real*8	a(3),b(3),c(3),center(3)
	real*8	c_ab(3),c_ac(3),c_bc(3)
	real*8	pabc(3),pacb(3),n(3)
	real*8	der_shabc(3),der_shacb(3),der_shbca(3)
	real*8	der1(3),der2(3),der3(3),der4(3),der5(3),der6(3)
	real*8	cosine(3),sine(3)
c
	common /constants/ pi,twopi,precision
c
	save
c
	call center2(a,b,ra2,rb2,rab2,c_ab,l1)
	call center2(a,c,ra2,rc2,rac2,c_ac,l2)
	call center2(b,c,rb2,rc2,rbc2,c_bc,l3)
c
	val1 = l1*rab
	val2 = l2*rac
	val3 = l3*rbc
c
	val1b = rab-val1
	val2b = rac-val2
	val3b = rbc-val3
c
	call center3(a,b,c,wa,wb,wc,center)
c
	call triangle_dual(a,b,c,center,eps,ra2,pabc,pacb,n)
c
	call tetra3(a,b,c,pabc,rab,rac,rbc,rab2,rac2,
     1	rbc2,ra,rb,rc,ra2,rb2,rc2,seg_ang_abc,seg_ang_acb,
     2	seg_ang_bca,ang_dih_abc,ang_dih_bac,ang_dih_cab,
     2	der1,der2,der3,der4,der5,der6,cosine,sine)
c
	cos_abc = cosine(1)
	sin_abc = sine(1)
	cos_acb = cosine(2)
	sin_acb = sine(2)
	cos_bca = cosine(3)
	sin_bca = sine(3)
c
	sh_abc = eps*cos_abc/sin_abc
	sh_acb = eps*cos_acb/sin_acb
	sh_bca = eps*cos_bca/sin_bca
c
	der_shabc(1) = - val1b*l1*sin_abc*cos_abc/eps - eps*der1(1)
	der_shabc(2) = - eps*der1(2)
	der_shabc(3) = - eps*der1(3)
c
	der_shacb(1) = - eps*der2(1)
	der_shacb(2) = - val2b*l2*sin_acb*cos_acb/eps - eps*der2(2)
	der_shacb(3) = - eps*der2(3)
c
	der_shbca(1) = - eps*der3(1)
	der_shbca(2) = - eps*der3(2)
	der_shbca(3) = - val3b*l3*sin_bca*cos_bca/eps - eps*der3(3)
c
	return
	end
c
