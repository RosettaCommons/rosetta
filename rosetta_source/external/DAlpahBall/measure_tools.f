c	measure_tools.f		Version 1 : 9/1/2008

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


c	center2.f

c	Copyright (C) 2002 Patrice Koehl

c	This subroutine finds the orthogonal center C of two spheres A and B
c	(if radii of A and B are equal, C is just the midpoint of the centers
c	of A and B)

	subroutine center2(a,b,ra2,rb2,rab2,c,lamda)

	integer	i

	real*8	ra2,rb2,rab2,lamda,uml

	real*8	a(3),b(3),c(3)

	save

	lamda = 0.5d0 - (ra2-rb2)/(2*rab2)
	uml   = 1-lamda

	do 100 i = 1,3
		c(i) = lamda*a(i) + uml*b(i)
100	continue

	return
	end

c	Center3.f

c	Copyright (C) 2002 Patrice Koehl

c	This subroutine finds the orthogonal center Y of three spheres A and B
c	and C

	subroutine center3(a,b,c,i0,j0,k0,y)

	real*8	i0,j0,k0,a1,a2,a3,a4
	real*8	dx,dy,dz,d0
	real*8	val1,val2,val3
	real*8	a(3),b(3),c(3),y(3)

	save

	a1=(b(2)-a(2))*(c(3)-a(3)) - (c(2)-a(2))*(b(3)-a(3))
	a2=(b(3)-a(3))*(c(1)-a(1)) - (c(3)-a(3))*(b(1)-a(1))

	val1=b(1)*c(2)-c(1)*b(2)
	val2=a(1)*c(2)-c(1)*a(2)
	val3=a(1)*b(2)-b(1)*a(2)

	a3 = val1-val2+val3
	a4 = a(3)*val1 - b(3)*val2 + c(3)*val3

	d0 = -a1*a1 - a2*a2 - a3*a3

	val1 = i0*(b(3)-c(3)) - j0*(a(3)-c(3)) + k0*(a(3)-b(3))
	val2 = i0*(b(2)-c(2)) - j0*(a(2)-c(2)) + k0*(a(2)-b(2))
	val3 = i0*(b(1)-c(1)) - j0*(a(1)-c(1)) + k0*(a(1)-b(1))

	dx = -a4*a1 + a2*val1 -a3*val2
	dy = -a1*val1 - a4*a2 + a3*val3
	dz = a1*val2 -a2*val3 -a4*a3

	y(1) = dx/d0
	y(2) = dy/d0
	y(3) = dz/d0

	return
	end

c	Center4.f

c	Copyright (C) 2002 Patrice Koehl

c	This subroutine finds the orthogonal center Y of four spheres A, B,
c	C and D


	subroutine center4(a,b,c,d,SA,SB,SC,SD,y)

	integer i

	real*8  SA,SB,SC,SD,A1,B1,C1
	real*8  S(3),a(3),b(3),c(3),d(3),y(3),detval(4)
	real*8  ad(3),bd(3),cd(3)

	save

	S(1) = SA - SD
	S(2) = SB - SD
	S(3) = SC - SD

	do 200 i = 1,3
		ad(i) = a(i) - d(i)
		bd(i) = b(i) - d(i)
		cd(i) = c(i) - d(i)
200     continue

	A1 = bd(2)*cd(3) - cd(2)*bd(3)
	B1 = ad(2)*cd(3) - ad(3)*cd(2)
	C1 = ad(2)*bd(3) - ad(3)*bd(2)

	detval(1) = ad(1)*A1 - bd(1)*B1 + cd(1)*C1
	detval(2) = S(1)*A1 - S(2)*B1 + S(3)*C1

	A1 = S(2)*cd(1) - bd(1)*S(3)
	B1 = S(1)*cd(1) - S(3)*ad(1)
	C1 = S(1)*bd(1) - S(2)*ad(1)

	detval(3) = -ad(3)*A1 + bd(3)*B1 - cd(3)*C1
	detval(4) =  ad(2)*A1 - bd(2)*B1 + cd(2)*C1

	if(detval(1).ne.0) then
		
		y(1) = detval(2)/detval(1)
		y(2) = detval(3)/detval(1)
		y(3) = detval(4)/detval(1)

! 		write(*,*) "============ measure tools ============"
! 		write(*,*) "cen1",a
! 		write(*,*) "cen2",b
! 		write(*,*) "cen3",c
! 		write(*,*) "cen4",d
! 		write(*,*) "w",SA,SB,SC,SD
! 		write(*,*) "rtn:",y(1), y(2), y(3)
! 		write(*,*) "======================================="
		
	else

c	       here we are in the special case that the four spheres
c	       are co-planar, and the center3 of each subset of three
c	       spheres are equal. We take this center3 as the "center4"

		call center3(a,b,c,SA,SB,SC,y)

	endif

	return
	end
c	det4_1.f

c	Copyright (C) 2002 Patrice Koehl

c	This function computes a 4 by 4 determinant, in which the last column
c	only contains the number 1

	function det4_1(mat4)

	integer	i,j

	real*8	det4_1
	real*8	val1,val2,val3
	real*8 mat4(4,4), mat2(3,3)

	do 200 i = 1,3
		do 100 j = 1,3
			mat2(i,j) = mat4(i,j)-mat4(4,j)
100		continue
200	continue

	val1 = mat2(2,2)*mat2(3,3)-mat2(3,2)*mat2(2,3)
	val2 = mat2(1,2)*mat2(3,3)-mat2(3,2)*mat2(1,3)
	val3 = mat2(1,2)*mat2(2,3)-mat2(2,2)*mat2(1,3)

	det4_1 = -mat2(1,1)*val1 +mat2(2,1)*val2 -mat2(3,1)*val3

	return
	end

c	triangle_dual.f

c	Copyright (C) 2002 Patrice Koehl

c	This subroutine computes the dual of a triangle (a,b,c), i.e.
c	the point P_ABC at which the 3 spheres centered at a, b and c
c	with radii ra, rb and rc intersect, and (a,b,c,P_ABC) has
c	positive orientation

	subroutine triangle_dual(a,b,c,center,eps,ra2,d1,d2,n)

	integer	i

	real*8	ra2
	real*8	s3,eps
	real*8	a(3),b(3),c(3),d1(3),d2(3)
	real*8	u1(3),u2(3),n(3),center(3)

	save

	call diffvect(b,a,u1)
	call diffvect(c,a,u2)

	call crossvect(u1,u2,d1)
	call unitvector(d1,n)
	call diffvect(center,a,u1)

	call dotvect(u1,u1,s3)

	eps = sqrt(ra2-s3)

	do 100 i = 1,3
		d1(i) = center(i) + eps*n(i)
		d2(i) = center(i) - eps*n(i)
100	continue

	return
	end

c	tetra_volume.f

c	Copyright (C) 2002 Patrice Koehl

c	This subroutine computes the volume of a tetrahedron

	function tetra_volume(a,b,c,d)

	integer	i

	real*8	det,det4_1
	real*8	tetra_volume
	real*8	a(3),b(3),c(3),d(3)
	real*8	mat4(4,4)

	save

	do 100 i = 1,4
		mat4(i,4) = 1.d0
100	continue

	do 200 i = 1,3
		mat4(1,i) = a(i)
		mat4(2,i) = b(i)
		mat4(3,i) = c(i)
		mat4(4,i) = d(i)
200	continue

	det = det4_1(mat4)

	tetra_volume = abs(det)/6

	return
	end

c	Angle_dihed.f

c	Copyright (C) 2002 Patrice Koehl

	subroutine angle_dihed(a,b,c,d,ang,cos_ang)

	real*8	ang,cos_ang
	real*8	pi,twopi,precision

	real*8  a(3),b(3),c(3),d(3),u1(3),u2(3),m(3),n1(3),n2(3)

	common /constants/ pi,twopi,precision

	save

	call diffvect(c,a,u1)
	call diffvect(c,b,u2)

	call crossvect(u1,u2,m)
	call unitvector(m,n1)

	call diffvect(d,a,u1)
	call diffvect(d,b,u2)

	call crossvect(u1,u2,m)
	call unitvector(m,n2)

	call dotvect(n1,n2,cos_ang)

	ang = acos(cos_ang)/twopi

	return
	end

c	Distance.f

c	Copyright (C) 2002 Patrice Koehl

c	This subroutine computes the square of the distance between two
c	sphere centers

	subroutine distance(coord,n1,n2,dist,ncortot)

	integer	i,n1,n2,ncortot

	real*8	dist,val
	real*8	coord(ncortot)

	save

	dist = 0
	do 100 i = 1,3
		val = coord(3*(n1-1)+i)-coord(3*(n2-1)+i)
		dist = dist + val*val
100	continue
	dist = sqrt(dist)
	
	return
	end

c	Distance2.f

c	Copyright (C) 2002 Patrice Koehl

c	This subroutine computes the square of the distance between two
c	sphere centers

	subroutine distance2(coord,n1,n2,dist,ncortot)

	integer	i,n1,n2,ncortot

	real*8	dist,val
	real*8	coord(ncortot)

	save

	dist = 0
	do 100 i = 1,3
		val = coord(3*(n1-1)+i)-coord(3*(n2-1)+i)
		dist = dist + val*val
100	continue

	return
	end

C	tetra3_noder.f

c	Copyright (C) 2002 Patrice Koehl

c	This subroutine computes the volume of a tetrahedron, its
c	six dihedral angles, as well as the derivatives of the six
c	dihedral angles with respect to 3 distances (the 3 other
c	distances are considered fixed)

	subroutine tetra3_noder(a,b,c,p,rab,rac,rbc,ra,rb,rc,ang1,
     1			ang2,ang3,ang4,ang5,ang6,cosine,sine,option)

	integer	option

	real*8	tetra_volume

	real*8	pi,twopi,precision
	real*8	ra,rb,rc
	real*8	rab,rac,rbc
	real*8	ang1,ang2,ang3,ang4,ang5,ang6
	real*8	cos_ang1,cos_ang2,cos_ang3
	real*8	cos_ang4,cos_ang5,cos_ang6
	real*8	sin_ang1,sin_ang2,sin_ang3
	real*8	vol
	real*8	val1,val2
	real*8	sum
	real*8	s1_2,s1
	real*8	s2_2,s2
	real*8	s3_2,s3
	real*8	s4_2,s4,sum_s_2
	real*8	c1,c2,c3,c4
	real*8	a(3),b(3),c(3),p(3)
	real*8	cosine(3),sine(3)

	common /constants/pi,twopi,precision

	save

c	Volume of the tetrahedron, based on all edge lengths
c	(in fact, this is the square of the volume, multiplied
c	by 288)
c	The derivative are really (derivative/vol) and are
c	computed from the fact that

c	vol2 = 2*ra2*(4*rb2*rc2 - val_bc*val_bc)
c     1		+ val_ab*(-2*rc2*val_ab - val_bc*val_ac)
c     2		- val_ac*(val_ab*val_bc + 2*rb2*val_ac)

c	where vol2 = 288*vol^2

	vol = tetra_volume(a,b,c,p)

c	Surfaces s1,s2,s3,s4 of the four faces of the tetrahedron,

c	(We use the fact that for a triangle T with side lengths
c	a,b,c, then

c	P (=perimeter) = (a+b+c)/2
c	Surf**2 = p*(p-a)*(p-b)*(p-c)

c	The four triangles considered are:

c		Triangle	Surface

c		T1 : BCP	s1
c		T2 : ACP	s2
c		T3 : ABP	s3
c		T4 : ABC	s4

	sum = (rb + rc + rbc)/2
	s1_2 = sum*(sum-rb)*(sum-rc)*(sum-rbc)
	s1   = sqrt(s1_2)

	sum = (ra + rc + rac)/2
	s2_2 = sum*(sum-ra)*(sum-rc)*(sum-rac)
	s2   = sqrt(s2_2)

	sum = (ra + rb + rab)/2
	s3_2 = sum*(sum-ra)*(sum-rb)*(sum-rab)
	s3   = sqrt(s3_2)

	sum = (rab + rac + rbc)/2
	s4_2 = sum*(sum-rab)*(sum-rac)*(sum-rbc)
	s4 = sqrt(s4_2)

	sum_s_2 = s1_2 + s2_2 + s3_2 + s4_2

c	Get all six dihedral angles

c	ang1 = angle_dihed(a,b,c,p) = angle_dihed(T3,T4)
c	ang2 = angle_dihed(a,c,b,p) = angle_dihed(T2,T4)
c	ang3 = angle_dihed(b,c,a,p) = angle_dihed(T1,T4)
c	ang4 = angle_dihed(a,p,b,c) = angle_dihed(T2,T3)
c	ang5 = angle_dihed(b,p,a,c) = angle_dihed(T1,T3)
c	ang6 = angle_dihed(c,p,a,b) = angle_dihed(T1,T2)

	call angle_dihed(b,p,a,c,ang5,cos_ang5)
	call angle_dihed(c,p,a,b,ang6,cos_ang6)

c	Get 4 other dihedral angle using cosine rule

	val1 = 2*s1*s3*cos_ang5
	val2 = 2*s1*s2*cos_ang6

	c1 = sum_s_2 - 2*s1_2
	c2 = sum_s_2 - 2*s2_2 - val1
	c3 = sum_s_2 - 2*s3_2 - val2
	c4 = sum_s_2 - 2*s4_2 - val1 - val2

	cos_ang1 = (c1+c2-c3-c4)/(4*s3*s4)
	cos_ang2 = (c1-c2+c3-c4)/(4*s2*s4)
	cos_ang3 = (-c1+c2+c3+c4)/(4*s1*s4)
	cos_ang4 = c4/(2*s2*s3)

	if(option.eq.1) then
		ang1 = acos(cos_ang1)/twopi
		ang2 = acos(cos_ang2)/twopi
		ang3 = acos(cos_ang3)/twopi
		ang4 = acos(cos_ang4)/twopi
	endif

	sin_ang1 = 1.5d0*rab*vol/(s3*s4)
	sin_ang2 = 1.5d0*rac*vol/(s2*s4)
	sin_ang3 = 1.5d0*rbc*vol/(s1*s4)
c	sin_ang4 = 1.5d0*ra*vol/(s2*s3)
c	sin_ang5 = 1.5d0*rb*vol/(s1*s3)
c	sin_ang6 = 1.5d0*rc*vol/(s1*s2)

	cosine(1) = cos_ang1
	cosine(2) = cos_ang2
	cosine(3) = cos_ang3
	sine(1) = sin_ang1
	sine(2) = sin_ang2
	sine(3) = sin_ang3

	return
	end
c	tetra6.f

c	Copyright (C) 2002 Patrice Koehl

c	This subroutine computes the volume of a tetrahedron, its
c	six dihedral angles, as well as the derivatives of the six
c	dihedral angles with respect to 6 distances (the 6 edge lengths)

	subroutine tetra6(a,b,c,d,rab,rac,rad,rbc,rbd,rcd,rab2,rac2,
     1			rad2,rbc2,rbd2,rcd2,ang1,ang2,ang3,
     2			ang4,ang5,ang6,der1,der2,der3,der4,der5,der6,
     3			cosine,sine)

	integer	i,flag(6)

	real*8	tetra_volume

	real*8	pi,twopi,precision
	real*8	rad,rbd,rcd,rad2,rbd2,rcd2
	real*8	rab,rac,rbc,rab2,rac2,rbc2
	real*8	ang1,ang2,ang3,ang4,ang5,ang6
	real*8	cos_ang1,cos_ang2,cos_ang3
	real*8	cos_ang4,cos_ang5,cos_ang6
	real*8	sin_ang1,sin_ang2,sin_ang3
	real*8	sin_ang4,sin_ang5,sin_ang6
	real*8	vol,vol2
	real*8	val1,val2,val3,val4,val5,val6
	real*8	val,val_ab,val_ac,val_ad,val_bc,val_bd,val_cd
	real*8	val_bc2,val_bd2,val_ab2
	real*8	sum,deter
	real*8	s1_2,s1
	real*8	s2_2,s2
	real*8	s3_2,s3
	real*8	s4_2,s4,sum_s_2
	real*8	c1,c2,c3,c4
	real*8	a(3),b(3),c(3),d(3)
	real*8	a1(6),a2(6),a3(6),a4(6),b1(6),b2(6),b3(6),b4(6)
	real*8	der1(6),der2(6),der3(6),der4(6),der5(6),der6(6)
	real*8	dvol(6),ds1(6),ds2(6),ds3(6),ds4(6)
	real*8	cosine(6),sine(6)

	logical testder

	common /constants/ pi,twopi,precision

	save

	val_ab = rab2-rad2-rbd2
	val_ac = rac2-rad2-rcd2
	val_bc = rbc2-rbd2-rcd2

c	Volume of the tetrahedron, based on all edge lengths
c	(in fact, this is the square of the volume, multiplied
c	by 288)

c	288*V**2 = | 	0	1	1	1	1	|
c		   |	1	0	rad2	rbd2	rcd2	|
c		   |	1	rad2	0	rab2	rac2	|
c		   |	1	rbd2	rab2	0	rbc2	|
c		   |	1	rcd2	rac2	rbc2	0	|

c	vol2 = 2*rad2*(4*rbd2*rcd2 - val_bc*val_bc)
c     1		+ val_ab*(-2*rcd2*val_ab - val_bc*val_ac)
c     2		- val_ac*(val_ab*val_bc + 2*rbd2*val_ac)

c	Notice: 
c	The volume of a tetrahedron can also be computed from:

c	6*V = | a(1)	a(2)	a(3)	1|
c	      | b(1)	b(2)	b(3)	1|
c	      | c(1)	c(2)	c(3)	1|
c	      | d(1)	d(2)	d(3)	1|

c	where a, b, c and d are the 4 vertices of the tetrahedron

c	The derivative are really (derivative/vol)

	vol = tetra_volume(a,b,c,d)
	vol2 = 144*vol*vol

	val = -rab*(2*rcd2*val_ab + val_bc*val_ac)
	dvol(1) = val/vol2

	val = - rac*(val_ab*val_bc + 2*rbd2*val_ac)
	dvol(2) = val/vol2

	val = - rbc*(val_ab*val_ac + 2*rad2*val_bc)
	dvol(4) = val/vol2

c	By changing the computation of the determinant, 
c	vol2 is also:
c	vol2 = 2*rcd2(4*rac2*rbc2-val_ab_c**2)
c		+ val_ad_c*(-2*rbc2*val_ad_c-val_ab_c*val_bd_c)
c		- val_bd_c*(val_ad_c*val_ab_c + 2*rac2*val_bd_c)

c	with 
c		val_ad_c = rad2 - rac2-rcd2
c		val_ab_c = rab2 - rac2 -rbc2
c		val_bd_c = rbd2 - rbc2 -rcd2

c	This formula provides a better start to compute dvol_ad and
c	dvol_bd:

	val_ad = rad2 - rac2-rcd2
	val_ab2 = rab2 - rac2 - rbc2
	val_bd = rbd2 - rbc2 - rcd2

	val = -rad*(2*rbc2*val_ad +val_ab2*val_bd)
	dvol(3) = val/vol2

	val = -rbd*(2*rac2*val_bd + val_ad*val_ab2)
	dvol(5) = val/vol2

c	By re-modifying the computation of Vol, we obtain:

	val_cd = rcd2-rad2-rac2
	val_bd2 = rbd2 - rab2-rad2
	val_bc2 = rbc2 - rab2 -rac2

	val = -rcd*(2*rab2*val_cd + val_bd2*val_bc2)
	dvol(6) = val/vol2

c	Surfaces s1,s2,s3,s4 of the four faces of the tetrahedron,

c	(We use the fact that for a triangle T with side lengths
c	a,b,c, then

c	P (=perimeter) = (a+b+c)/2
c	Surf**2 = p*(p-a)*(p-b)*(p-c)

c	The four triangles considered are:

c		Triangle	Surface

c		T1 : BCP	s1
c		T2 : ACP	s2
c		T3 : ABP	s3
c		T4 : ABC	s4

c	and their derivatives with respect to distances (in fact,
c	Der(S)/S )

	sum = (rbd + rcd + rbc)/2
	s1_2 = sum*(sum-rbd)*(sum-rcd)*(sum-rbc)
	s1   = sqrt(s1_2)
	val1 = 1/sum
	val2 = 1/(sum-rbd)
	val3 = 1/(sum-rcd)
	val4 = 1/(sum-rbc)
	ds1(1) = 0
	ds1(2) = 0
	ds1(3) = 0
	ds1(4) = 0.25d0*(val1+val2+val3-val4)
	ds1(5) = 0.25d0*(val1-val2+val3+val4)
	ds1(6) = 0.25d0*(val1+val2-val3+val4)

	sum = (rad + rcd + rac)/2
	s2_2 = sum*(sum-rad)*(sum-rcd)*(sum-rac)
	s2   = sqrt(s2_2)
	val1 = 1/sum
	val2 = 1/(sum-rad)
	val3 = 1/(sum-rcd)
	val4 = 1/(sum-rac)
	ds2(1) = 0
	ds2(2) = 0.25d0*(val1+val2+val3-val4)
	ds2(3) = 0.25d0*(val1-val2+val3+val4)
	ds2(4) = 0
	ds2(5) = 0
	ds2(6) = 0.25d0*(val1+val2-val3+val4)

	sum = (rad + rbd + rab)/2
	s3_2 = sum*(sum-rad)*(sum-rbd)*(sum-rab)
	s3   = sqrt(s3_2)
	val1 = 1/sum
	val2 = 1/(sum-rad)
	val3 = 1/(sum-rbd)
	val4 = 1/(sum-rab)
	ds3(1) = 0.25d0*(val1+val2+val3-val4)
	ds3(2) = 0
	ds3(3) = 0.25d0*(val1-val2+val3+val4)
	ds3(4) = 0
	ds3(5) = 0.25d0*(val1+val2-val3+val4)
	ds3(6) = 0

	sum = (rab + rac + rbc)/2
	s4_2 = sum*(sum-rab)*(sum-rac)*(sum-rbc)
	s4 = sqrt(s4_2)
	val1 = 1/sum
	val2 = 1/(sum-rab)
	val3 = 1/(sum-rac)
	val4 = 1/(sum-rbc)
	ds4(1) = 0.25d0*(val1-val2+val3+val4)
	ds4(2) = 0.25d0*(val1+val2-val3+val4)
	ds4(3) = 0
	ds4(4) = 0.25d0*(val1+val2+val3-val4)
	ds4(5) = 0
	ds4(6) = 0

	sum_s_2 = s1_2 + s2_2 + s3_2 + s4_2

c	Get all six dihedral angles

c	ang1 = angle_dihed(a,b,c,d) = angle_dihed(T3,T4)
c	ang2 = angle_dihed(a,c,b,d) = angle_dihed(T2,T4)
c	ang3 = angle_dihed(b,c,a,d) = angle_dihed(T1,T4)
c	ang4 = angle_dihed(a,d,b,c) = angle_dihed(T2,T3)
c	ang5 = angle_dihed(b,d,a,c) = angle_dihed(T1,T3)
c	ang6 = angle_dihed(c,d,a,b) = angle_dihed(T1,T2)

	call angle_dihed(b,d,a,c,ang5,cos_ang5)
	call angle_dihed(c,d,a,b,ang6,cos_ang6)

c	Get 4 other dihedral angle using cosine rule

	val1 = 2*s1*s3*cos_ang5
	val2 = 2*s1*s2*cos_ang6

	c1 = sum_s_2 - 2*s1_2
	c2 = sum_s_2 - 2*s2_2 - val1
	c3 = sum_s_2 - 2*s3_2 - val2
	c4 = sum_s_2 - 2*s4_2 - val1 - val2

	cos_ang1 = (c1+c2-c3-c4)/(4*s3*s4)
	cos_ang2 = (c1-c2+c3-c4)/(4*s2*s4)
	cos_ang3 = (-c1+c2+c3+c4)/(4*s1*s4)
	cos_ang4 = c4/(2*s2*s3)

	ang1 = acos(cos_ang1)/twopi
	ang2 = acos(cos_ang2)/twopi
	ang3 = acos(cos_ang3)/twopi
	ang4 = acos(cos_ang4)/twopi

	sin_ang1 = 1.5d0*rab*vol/(s3*s4)
	sin_ang2 = 1.5d0*rac*vol/(s2*s4)
	sin_ang3 = 1.5d0*rbc*vol/(s1*s4)
	sin_ang4 = 1.5d0*rad*vol/(s2*s3)
	sin_ang5 = 1.5d0*rbd*vol/(s1*s3)
	sin_ang6 = 1.5d0*rcd*vol/(s1*s2)

	cosine(1) = cos_ang1
	cosine(2) = cos_ang2
	cosine(3) = cos_ang3
	cosine(4) = cos_ang4
	cosine(5) = cos_ang5
	cosine(6) = cos_ang6

	sine(1) = sin_ang1
	sine(2) = sin_ang2
	sine(3) = sin_ang3
	sine(4) = sin_ang4
	sine(5) = sin_ang5
	sine(6) = sin_ang6

c	Now get derivatives

	testder = .true.

	if(abs(cos_ang1).gt.precision) then
		val = sin_ang1/cos_ang1
		do 100 i = 1,6
			der1(i) = val*(dvol(i)-ds3(i)-ds4(i))
100		continue
		der1(1) = der1(1) + val/rab
		flag(1) = 1
	else
		flag(1) = 0
		testder = .false.
	endif

	if(abs(cos_ang2).gt.precision) then
		val = sin_ang2/cos_ang2
		do 300 i = 1,6
			der2(i) = val*(dvol(i)-ds2(i)-ds4(i))
300		continue
		der2(2) = der2(2) + val/rac
		flag(2) = 1
	else
		flag(2) = 0
		testder = .false.
	endif

	if(abs(cos_ang3).gt.precision) then
		val = sin_ang3/cos_ang3
		do 500 i = 1,6
			der3(i)=val*(dvol(i) -ds1(i)-ds4(i))
500		continue
		der3(4) = der3(4) +val/rbc
		flag(3) = 1
	else
		flag(3) = 0
		testder = .false.
	endif

	if(abs(cos_ang4).gt.precision) then
		val = sin_ang4/cos_ang4
		do 700 i = 1,6
			der4(i) = val*(dvol(i)-ds2(i)-ds3(i))
700		continue
		der4(3) = der4(3) +val/rad
		flag(4) = 1
	else
		flag(4) = 0
		testder = .false.
	endif

	if(abs(cos_ang5).gt.precision) then
		val = sin_ang5/cos_ang5
		do 900 i = 1,6
			der5(i) = val*(dvol(i)-ds1(i)-ds3(i))
900		continue
		der5(5) = der5(5) + val/rbd
		flag(5) = 1
	else
		flag(5) = 0
		testder = .false.
	endif

	if(abs(cos_ang6).gt.precision) then
		val = sin_ang6/cos_ang6
		do 1100 i = 1,6
			der6(i) = val*(dvol(i)-ds1(i)-ds2(i))
1100		continue
		der6(6) = der6(6) +val/rcd
		flag(6) = 1
	else
		flag(6) = 0
		testder = .false.
	endif

	if(testder) goto 3400

c	Use cosine rules in case one of the dihedral angle is pi/2

	do 1300 i = 1,6
		sum_s_2 = -s1_2*ds1(i)-s2_2*ds2(i)-
     1			s3_2*ds3(i)-s4_2*ds4(i)
		val1 = s3*s4*(ds3(i)+ds4(i))*cos_ang1
		val2 = s2*s4*(ds2(i)+ds4(i))*cos_ang2
		val3 = s1*s4*(ds1(i)+ds4(i))*cos_ang3
		val4 = s2*s3*(ds2(i)+ds3(i))*cos_ang4
		val5 = s1*s3*(ds1(i)+ds3(i))*cos_ang5
		val6 = s1*s2*(ds1(i)+ds2(i))*cos_ang6
		a1(i) = sum_s_2 + 2*s1_2*ds1(i)
     1		+ val1+val2+val4
		a2(i) = sum_s_2 + 2*s2_2*ds2(i)
     1		+ val1+val3+val5
		a3(i) = sum_s_2 + 2*s3_2*ds3(i)
     1		+ val2+val3+val6
		a4(i) = sum_s_2 + 2*s4_2*ds4(i)
     1		+ val4+val5+val6
1300	continue

c	Since a tetrahedron contains a maximum of 3 right dihedral
c	angles, we only have to consider ten pairs:

	val1 = 1.5d0*rab*vol
	val2 = 1.5d0*rac*vol
	val3 = 1.5d0*rbc*vol
	val4 = 1.5d0*rad*vol
	val5 = 1.5d0*rbd*vol
	val6 = 1.5d0*rcd*vol

	if(flag(1).ne.0.and.flag(2).ne.0) then

		deter = 2
		do 1400 i = 1,6
			b1(i) = a1(i) - val1*der1(i) - val2*der2(i)
			b2(i) = a2(i) - val1*der1(i)
			b3(i) = a3(i) - val2*der2(i)
			b4(i) = a4(i)
1400		continue

		do 1500 i = 1,6
			if(flag(3).eq.0) der3(i) = (b1(i)+b2(i)
     1					+b3(i)-b4(i))/(deter*val3)
			if(flag(4).eq.0) der4(i) = 2*b1(i)/
     1					(deter*val4)
			if(flag(5).eq.0) der5(i) = (-b1(i)+b2(i)
     1					-b3(i)+b4(i))/(deter*val5)
			if(flag(6).eq.0) der6(i) = (-b1(i)-b2(i)
     1					+b3(i)+b4(i))/(deter*val6)
1500		continue

	elseif(flag(1).ne.0.and.flag(3).ne.0) then

		deter = 2
		do 1600 i = 1,6
			b1(i) = a1(i) - val1*der1(i) 
			b2(i) = a2(i) - val1*der1(i) -val3*der3(i)
			b3(i) = a3(i) - val3*der3(i)
			b4(i) = a4(i)
1600		continue

		do 1700 i = 1,6
			if(flag(2).eq.0) der2(i) = (b1(i)+b2(i)
     1					+b3(i)-b4(i))/(deter*val2)
			if(flag(4).eq.0) der4(i) = (b1(i)-b2(i)
     1					-b3(i)+b4(i))/(deter*val4)
			if(flag(5).eq.0) der5(i) = 2*b2(i)/
     1					(deter*val5)
			if(flag(6).eq.0) der6(i) = (-b1(i)-b2(i)
     1					+b3(i)+b4(i))/(deter*val6)
1700		continue

	elseif(flag(1).ne.0.and.flag(4).ne.0) then

		deter = -2
		do 1800 i = 1,6
			b1(i) = a1(i) - val1*der1(i)  -val4*der4(i)
			b2(i) = a2(i) - val1*der1(i)
			b3(i) = a3(i) 
			b4(i) = a4(i) -val4*der4(i)
1800		continue

		do 1900 i = 1,6
			if(flag(2).eq.0) der2(i) = -2*b1(i)/
     1					(deter*val2)
			if(flag(3).eq.0) der3(i) = (b1(i)-b2(i)
     1					-b3(i)+b4(i))/(deter*val3)
			if(flag(5).eq.0) der5(i) = (-b1(i)-b2(i)
     1					+b3(i)-b4(i))/(deter*val5)
			if(flag(6).eq.0) der6(i) = (b1(i)+b2(i)
     1					-b3(i)-b4(i))/(deter*val6)
1900		continue

	elseif(flag(1).ne.0.and.flag(5).ne.0) then

		deter = -2
		do 2000 i = 1,6
			b1(i) = a1(i) - val1*der1(i) 
			b2(i) = a2(i) - val1*der1(i)  -val5*der5(i)
			b3(i) = a3(i) 
			b4(i) = a4(i) -val5*der5(i)
2000		continue

		do 2100 i = 1,6
			if(flag(2).eq.0) der2(i) = (-b1(i)+b2(i)
     1					-b3(i)+b4(i))/(deter*val2)
			if(flag(3).eq.0) der3(i) = -2*b2(i)/
     1					(deter*val3)
			if(flag(4).eq.0) der4(i) = (-b1(i)-b2(i)
     1					+b3(i)-b4(i))/(deter*val4)
			if(flag(6).eq.0) der6(i) = (b1(i)+b2(i)-
     1					b3(i)-b4(i))/(deter*val6)
2100		continue

	elseif(flag(2).ne.0.and.flag(3).ne.0) then

		deter = 2
		do 2200 i = 1,6
			b1(i) = a1(i) - val2*der2(i)
			b2(i) = a2(i) - val3*der3(i)
			b3(i) = a3(i) - val2*der2(i) -val3*der3(i)
			b4(i) = a4(i) 
2200		continue

		do 2300 i = 1,6
			if(flag(1).eq.0) der1(i) = (b1(i)+b2(i)
     1					+b3(i)-b4(i))/(deter*val1)
			if(flag(4).eq.0) der4(i) = (b1(i)-b2(i)
     1					-b3(i)+b4(i))/(deter*val4)
			if(flag(5).eq.0) der5(i) = (-b1(i)+b2(i)
     1					-b3(i)+b4(i))/(deter*val5)
			if(flag(6).eq.0) der6(i) = 2*b3(i)/
     1					(deter*val6)
2300		continue

	elseif(flag(2).ne.0.and.flag(4).ne.0) then

		deter = -2
		do 2400 i = 1,6
			b1(i) = a1(i) - val2*der2(i) - val4*der4(i)
			b2(i) = a2(i) 
			b3(i) = a3(i) - val2*der2(i)
			b4(i) = a4(i) -val4*der4(i)
2400		continue

		do 2500 i = 1,6
			if(flag(1).eq.0) der1(i) = -2*b1(i)/
     1					(deter*val1)
			if(flag(3).eq.0) der3(i) = (b1(i)-b2(i)
     1					-b3(i)+b4(i))/(deter*val3)
			if(flag(5).eq.0) der5(i) = (b1(i)-b2(i)
     1					+b3(i)-b4(i))/(deter*val5)
			if(flag(6).eq.0) der6(i) = (-b1(i)+b2(i)
     1					-b3(i)-b4(i))/(deter*val6)
2500		continue

	elseif(flag(2).ne.0.and.flag(6).ne.0) then

		deter = 2
		do 2600 i = 1,6
			b1(i) = a1(i) - val2*der2(i)
			b2(i) = a2(i) 
			b3(i) = a3(i) - val2*der2(i)  -val6*der6(i)
			b4(i) = a4(i) -val6*der6(i)
2600		continue

		do 2700 i = 1,6
			if(flag(1).eq.0) der1(i) = (b1(i)+b2(i)
     1					-b3(i)-b4(i))/(deter*val1)
			if(flag(3).eq.0) der3(i) = 2*b3(i)/
     1					(deter*val3)
			if(flag(4).eq.0) der4(i) = (b1(i)-b2(i)
     1					+b3(i)+b4(i))/(deter*val4)
			if(flag(5).eq.0) der5(i) = (-b1(i)+b2(i)
     1					-b3(i)+b4(i))/(deter*val5)
2700		continue

	elseif(flag(3).ne.0.and.flag(5).ne.0) then

		deter = 2
		do 2800 i = 1,6
			b1(i) = a1(i)
			b2(i) = a2(i) -val3*der3(i) - val5*der5(i)
			b3(i) = a3(i) - val3*der3(i)
			b4(i) = a4(i) -val5*der5(i)
2800		continue

		do 2900 i = 1,6
			if(flag(1).eq.0) der1(i) = 2*b2(i)/
     1					(deter*val1)
			if(flag(2).eq.0) der2(i) = (b1(i)-b2(i)
     1					+b3(i)-b4(i))/(deter*val2)
			if(flag(4).eq.0) der4(i) = (b1(i)-b2(i)
     1					-b3(i)+b4(i))/(deter*val4)
			if(flag(6).eq.0) der6(i) = (-b1(i)+b2(i)
     1					+b3(i)+b4(i))/(deter*val6)
2900		continue

	elseif(flag(3).ne.0.and.flag(6).ne.0) then

		deter = 2
		do 3000 i = 1,6
			b1(i) = a1(i)
			b2(i) = a2(i) -val3*der3(i)
			b3(i) = a3(i) - val3*der3(i) - val6*der6(i)
			b4(i) = a4(i) -val6*der6(i)
3000		continue

		do 3100 i = 1,6
			if(flag(1).eq.0) der1(i) = (b1(i)+b2(i)
     1					-b3(i)-b4(i))/(deter*val1)
			if(flag(2).eq.0) der2(i) = 2*b3(i)/
     1					(deter*val2)
			if(flag(4).eq.0) der4(i) = (b1(i)-b2(i)
     1					-b3(i)+b4(i))/(deter*val4)
			if(flag(5).eq.0) der5(i) = (-b1(i)+b2(i)+
     1					b3(i)+b4(i))/(deter*val5)
3100		continue

	elseif(flag(4).ne.0.and.flag(5).ne.0) then

		deter = -2
		do 3200 i = 1,6
			b1(i) = a1(i) -val4*der4(i)
			b2(i) = a2(i) -val5*der5(i)
			b3(i) = a3(i) 
			b4(i) = a4(i) -val4*der4(i) - val5*der5(i)
3200		continue

		do 3300 i = 1,6
			if(flag(1).eq.0) der1(i) = (-b1(i)-b2(i)
     1					+b3(i)-b4(i))/(deter*val1)
			if(flag(2).eq.0) der2(i) = (-b1(i)+b2(i)
     1					-b3(i)+b4(i))/(deter*val2)
			if(flag(3).eq.0) der3(i) = (b1(i)-b2(i)
     1					-b3(i)+b4(i))/(deter*val3)
			if(flag(6).eq.0) der6(i) = -2*b4(i)/
     1					(deter*val6)
3300		continue

	endif

c	For sake of consistency with calling program, invert ang3 and ang4

3400	continue

	val = ang3
	ang3 = ang4
	ang4 = val

	do 3500 i = 1,6
		val = der3(i)
		der3(i) = der4(i)
		der4(i) = val
3500	continue

	return
	end

c	tetra_6dihed.f

c	Copyright (C) 2002 Patrice Koehl

c	Get all six dihedral angles

c	ang1 = angle_dihed(a,b,c,d) = angle_dihed(T3,T4)
c	ang2 = angle_dihed(a,c,b,d) = angle_dihed(T2,T4)
c	ang3 = angle_dihed(b,c,a,d) = angle_dihed(T1,T4)
c	ang4 = angle_dihed(a,d,b,c) = angle_dihed(T2,T3)
c	ang5 = angle_dihed(b,d,a,c) = angle_dihed(T1,T3)
c	ang6 = angle_dihed(c,d,a,b) = angle_dihed(T1,T2)

	subroutine tetra_6dihed(a,b,c,d,ang1,ang2,ang3,ang4,ang5,ang6)

	real*8	ang1,ang2,ang3,ang4,ang5,ang6
	real*8	cos_ang
	real*8	pi,twopi,precision

	real*8  a(3),b(3),c(3),d(3)
	real*8	u_ab(3),u_ac(3),u_ad(3),u_bc(3),u_bd(3)

	real*8	u_abc(3),u_abd(3),u_acd(3),u_bcd(3)
	real*8	n_abc(3),n_abd(3),n_acd(3),n_bcd(3)

	common /constants/ pi,twopi,precision

	save

	call diffvect(a,b,u_ab)
	call diffvect(a,c,u_ac)
	call diffvect(a,d,u_ad)
	call diffvect(b,c,u_bc)
	call diffvect(b,d,u_bd)

	call crossvect(u_ab,u_ac,u_abc)
	call crossvect(u_ab,u_ad,u_abd)
	call crossvect(u_ac,u_ad,u_acd)
	call crossvect(u_bc,u_bd,u_bcd)

	call unitvector(u_abc,n_abc)
	call unitvector(u_abd,n_abd)
	call unitvector(u_acd,n_acd)
	call unitvector(u_bcd,n_bcd)

	call dotvect(n_abc,n_abd,cos_ang)
	ang1 = acos(cos_ang)/twopi

	call dotvect(n_abc,n_acd,cos_ang)
	ang2 = acos(-cos_ang)/twopi

	call dotvect(n_abc,n_bcd,cos_ang)
	ang3 = acos(cos_ang)/twopi

	call dotvect(n_abd,n_acd,cos_ang)
	ang4 = acos(cos_ang)/twopi

	call dotvect(n_abd,n_bcd,cos_ang)
	ang5 = acos(-cos_ang)/twopi

	call dotvect(n_acd,n_bcd,cos_ang)
	ang6 = acos(cos_ang)/twopi

	! write(*,*) "infunc",ang1,ang2,ang3,ang4,ang5,ang6

	return
	end
c	tetra3.f

c	Copyright (C) 2002 Patrice Koehl

c	This subroutine computes the volume of a tetrahedron, its
c	six dihedral angles, as well as the derivatives of the six
c	dihedral angles with respect to 3 distances (the 3 other
c	distances are considered fixed)

	subroutine tetra3(a,b,c,p,rab,rac,rbc,rab2,rac2,rbc2,
     1			ra,rb,rc,ra2,rb2,rc2,ang1,ang2,ang3,
     2			ang4,ang5,ang6,der1,der2,der3,der4,der5,
     3			der6,cosine,sine)

	integer	i,flag(6)

	real*8	tetra_volume
        real*8  pi,twopi,precision

	real*8	ra,rb,rc,ra2,rb2,rc2
	real*8	rab,rac,rbc,rab2,rac2,rbc2
	real*8	ang1,ang2,ang3,ang4,ang5,ang6
	real*8	cos_ang1,cos_ang2,cos_ang3
	real*8	cos_ang4,cos_ang5,cos_ang6
	real*8	sin_ang1,sin_ang2,sin_ang3
	real*8	sin_ang4,sin_ang5,sin_ang6
	real*8	vol,vol2
	real*8	val,val_ab,val_ac,val_bc
	real*8	val1,val2,val3,val4,val5,val6
	real*8	sum,deter
	real*8	s1_2,s1
	real*8	s2_2,s2
	real*8	s3_2,s3
	real*8	s4_2,s4,sum_s_2
	real*8	c1,c2,c3,c4
	real*8	a(3),b(3),c(3),p(3)
	real*8	a1(3),a2(3),a3(3),a4(3)
	real*8	b1(3),b2(3),b3(3),b4(3)
	real*8	der1(3),der2(3),der3(3),der4(3),der5(3),der6(3)
	real*8	dvol(3),ds1(3),ds2(3),ds3(3),ds4(3)
	real*8	cosine(3),sine(3)

	logical	testder

	common /constants/pi,twopi,precision

	save

	val_ab = rab2-ra2-rb2
	val_ac = rac2-ra2-rc2
	val_bc = rbc2-rb2-rc2

c	Volume of the tetrahedron, based on all edge lengths
c	(in fact, this is the square of the volume, multiplied
c	by 288)
c	The derivative are really (derivative/vol) and are
c	computed from the fact that

c	vol2 = 2*ra2*(4*rb2*rc2 - val_bc*val_bc)
c     1		+ val_ab*(-2*rc2*val_ab - val_bc*val_ac)
c     2		- val_ac*(val_ab*val_bc + 2*rb2*val_ac)

c	where vol2 = 288*vol^2

	vol = tetra_volume(a,b,c,p)
	vol2 = vol*vol*144

	val = -rab*(2*rc2*val_ab + val_bc*val_ac)
	dvol(1) = val/vol2

	val = - rac*(val_ab*val_bc + 2*rb2*val_ac)
	dvol(2) = val/vol2

	val = - rbc*(val_ab*val_ac + 2*ra2*val_bc)
	dvol(3) = val/vol2

c	Surfaces s1,s2,s3,s4 of the four faces of the tetrahedron,

c	(We use the fact that for a triangle T with side lengths
c	a,b,c, then

c	P (=perimeter) = (a+b+c)/2
c	Surf**2 = p*(p-a)*(p-b)*(p-c)

c	The four triangles considered are:

c		Triangle	Surface

c		T1 : BCP	s1
c		T2 : ACP	s2
c		T3 : ABP	s3
c		T4 : ABC	s4

c	and their derivatives with respect to distances (in fact,
c	Der(S)/S )

	sum = (rb + rc + rbc)/2
	s1_2 = sum*(sum-rb)*(sum-rc)*(sum-rbc)
	s1   = sqrt(s1_2)
	ds1(1) = 0
	ds1(2) = 0
	ds1(3) = 0.25d0*(1/sum + 1/(sum-rb) + 1/(sum-rc) -1/(sum-rbc))

	sum = (ra + rc + rac)/2
	s2_2 = sum*(sum-ra)*(sum-rc)*(sum-rac)
	s2   = sqrt(s2_2)
	ds2(1) = 0
	ds2(2) = 0.25d0*(1/sum+1/(sum-ra)+1/(sum-rc)-1/(sum-rac))
	ds2(3) = 0

	sum = (ra + rb + rab)/2
	s3_2 = sum*(sum-ra)*(sum-rb)*(sum-rab)
	s3   = sqrt(s3_2)
	ds3(1) = 0.25d0*(1/sum+1/(sum-ra) + 1/(sum-rb) - 1/(sum-rab))
	ds3(2) = 0
	ds3(3) = 0

	sum = (rab + rac + rbc)/2
	s4_2 = sum*(sum-rab)*(sum-rac)*(sum-rbc)
	s4 = sqrt(s4_2)
	ds4(1) = 0.25d0*(1/sum-1/(sum-rab) + 1/(sum-rac) + 1/(sum-rbc))
	ds4(2) = 0.25d0*(1/sum+1/(sum-rab) - 1/(sum-rac) + 1/(sum-rbc))
	ds4(3) = 0.25d0*(1/sum+1/(sum-rab) + 1/(sum-rac) - 1/(sum-rbc))

	sum_s_2 = s1_2 + s2_2 + s3_2 + s4_2

c	Get all six dihedral angles

c	ang1 = angle_dihed(a,b,c,p) = angle_dihed(T3,T4)
c	ang2 = angle_dihed(a,c,b,p) = angle_dihed(T2,T4)
c	ang3 = angle_dihed(b,c,a,p) = angle_dihed(T1,T4)
c	ang4 = angle_dihed(a,p,b,c) = angle_dihed(T2,T3)
c	ang5 = angle_dihed(b,p,a,c) = angle_dihed(T1,T3)
c	ang6 = angle_dihed(c,p,a,b) = angle_dihed(T1,T2)

	call angle_dihed(b,p,a,c,ang5,cos_ang5)
	call angle_dihed(c,p,a,b,ang6,cos_ang6)

c	Get 4 other dihedral angle using cosine rule

	val1 = 2*s1*s3*cos_ang5
	val2 = 2*s1*s2*cos_ang6

	c1 = sum_s_2 - 2*s1_2
	c2 = sum_s_2 - 2*s2_2 - val1
	c3 = sum_s_2 - 2*s3_2 - val2
	c4 = sum_s_2 - 2*s4_2 - val1 - val2

	cos_ang1 = (c1+c2-c3-c4)/(4*s3*s4)
	cos_ang2 = (c1-c2+c3-c4)/(4*s2*s4)
	cos_ang3 = (-c1+c2+c3+c4)/(4*s1*s4)
	cos_ang4 = c4/(2*s2*s3)

	ang1 = acos(cos_ang1)/twopi
	ang2 = acos(cos_ang2)/twopi
	ang3 = acos(cos_ang3)/twopi
	ang4 = acos(cos_ang4)/twopi

	sin_ang1 = 1.5d0*rab*vol/(s3*s4)
	sin_ang2 = 1.5d0*rac*vol/(s2*s4)
	sin_ang3 = 1.5d0*rbc*vol/(s1*s4)
	sin_ang4 = 1.5d0*ra*vol/(s2*s3)
	sin_ang5 = 1.5d0*rb*vol/(s1*s3)
	sin_ang6 = 1.5d0*rc*vol/(s1*s2)

	cosine(1) = cos_ang1
	cosine(2) = cos_ang2
	cosine(3) = cos_ang3
	sine(1) = sin_ang1
	sine(2) = sin_ang2
	sine(3) = sin_ang3

c	Now get derivatives

	testder = .true.

	if(abs(cos_ang1).gt.precision) then
		val = sin_ang1/cos_ang1
		do 100 i = 1,3
			der1(i) = val*(dvol(i) -ds3(i) -ds4(i))
100		continue
		der1(1) = der1(1) +val/rab
		flag(1) = 1
	else
		flag(1) = 0
		testder = .false.
	endif

	if(abs(cos_ang2).gt.precision) then
		val = sin_ang2/cos_ang2
		do 300 i = 1,3
			der2(i) = val*(dvol(i)-ds2(i)-ds4(i))
300		continue
		der2(2) = der2(2) + val/rac
		flag(2) = 1
	else
		flag(2) = 0
		testder = .false.
	endif

	if(abs(cos_ang3).gt.precision) then
		val = sin_ang3/cos_ang3
		do 500 i = 1,3
			der3(i)=val*(dvol(i) -ds1(i)-ds4(i))
500		continue
		der3(3) = der3(3) +val/rbc
		flag(3) = 1
	else
		flag(3) = 0
		testder = .false.
	endif

	if(abs(cos_ang4).gt.precision) then
		val = sin_ang4/cos_ang4
		do 700 i = 1,3
			der4(i) = val*(dvol(i)-ds2(i)-ds3(i))
700		continue
		flag(4) = 1
	else
		flag(4) = 0
		testder = .false.
	endif

	if(abs(cos_ang5).gt.precision) then
		val = sin_ang5/cos_ang5
		do 900 i = 1,3
			der5(i) = val*(dvol(i)-ds1(i)-ds3(i))
900		continue
		flag(5) = 1
	else
		flag(5) = 0
		testder = .false.
	endif

	if(abs(cos_ang6).gt.precision) then
		val = sin_ang6/cos_ang6
		do 1100 i = 1,3
			der6(i) = val*(dvol(i)-ds1(i)-ds2(i))
1100		continue
		flag(6) = 1
	else
		flag(6) = 0
		testder = .false.
	endif

	if(testder) return

c	Use cosine rules in case one of the dihedral angle is pi/2

	do 1300 i = 1,3
		sum_s_2 = -s1_2*ds1(i)-s2_2*ds2(i)-s3_2*ds3(i)-
     1		s4_2*ds4(i)
		val1 = s3*s4*(ds3(i)+ds4(i))*cos_ang1
		val2 = s2*s4*(ds2(i)+ds4(i))*cos_ang2
		val3 = s1*s4*(ds1(i)+ds4(i))*cos_ang3
		val4 = s2*s3*(ds2(i)+ds3(i))*cos_ang4
		val5 = s1*s3*(ds1(i)+ds3(i))*cos_ang5
		val6 = s1*s2*(ds1(i)+ds2(i))*cos_ang6
		a1(i) = sum_s_2 + 2*s1_2*ds1(i)
     1		+ val1+val2+val4
		a2(i) = sum_s_2 + 2*s2_2*ds2(i)
     1		+ val1+val3+val5
		a3(i) = sum_s_2 + 2*s3_2*ds3(i)
     1		+ val2+val3+val6
		a4(i) = sum_s_2 + 2*s4_2*ds4(i)
     1		+ val4+val5+val6
1300	continue

c	Since a tetrahedron contains a maximum of 3 right dihedral
c	angles, we only have to consider ten pairs:

	val1 = 1.5d0*rab*vol
	val2 = 1.5d0*rac*vol
	val3 = 1.5d0*rbc*vol
	val4 = 1.5d0*ra*vol
	val5 = 1.5d0*rb*vol
	val6 = 1.5d0*rc*vol

	if(flag(1).ne.0.and.flag(2).ne.0) then

		deter = 2
		do 1400 i = 1,3
			b1(i) = a1(i) - val1*der1(i) - val2*der2(i)
			b2(i) = a2(i) - val1*der1(i)
			b3(i) = a3(i) - val2*der2(i)
			b4(i) = a4(i)
1400		continue

		do 1500 i = 1,3
			if(flag(3).eq.0) der3(i) = (b1(i)+b2(i)+
     1					b3(i)-b4(i))/(deter*val3)
			if(flag(4).eq.0) der4(i) = 2*b1(i)/
     1					(deter*val4)
			if(flag(5).eq.0) der5(i) = (-b1(i)+b2(i)
     1					-b3(i)+b4(i))/(deter*val5)
			if(flag(6).eq.0) der6(i) = (-b1(i)-b2(i)
     1					+b3(i)+b4(i))/(deter*val6)
1500		continue

	elseif(flag(1).ne.0.and.flag(3).ne.0) then

		deter = 2
		do 1600 i = 1,3
			b1(i) = a1(i) - val1*der1(i) 
			b2(i) = a2(i) - val1*der1(i) -val3*der3(i)
			b3(i) = a3(i) - val3*der3(i)
			b4(i) = a4(i)
1600		continue

		do 1700 i = 1,3
			if(flag(2).eq.0) der2(i) = (b1(i)+b2(i)+
     1					b3(i)-b4(i))/(deter*val2)
			if(flag(4).eq.0) der4(i) = (b1(i)-b2(i)-
     1					b3(i)+b4(i))/(deter*val4)
			if(flag(5).eq.0) der5(i) = 2*b2(i)/
     1					(deter*val5)
			if(flag(6).eq.0) der6(i) = (-b1(i)-b2(i)
     1					+b3(i)+b4(i))/(deter*val6)
1700		continue

	elseif(flag(1).ne.0.and.flag(4).ne.0) then

		deter = -2
		do 1800 i = 1,3
			b1(i) = a1(i) - val1*der1(i)  -val4*der4(i)
			b2(i) = a2(i) - val1*der1(i)
			b3(i) = a3(i) 
			b4(i) = a4(i) -val4*der4(i)
1800		continue

		do 1900 i = 1,3
			if(flag(2).eq.0) der2(i) = -2*b1(i)/
     1					(deter*val2)
			if(flag(3).eq.0) der3(i) = (b1(i)-b2(i)
     1					-b3(i)+b4(i))/(deter*val3)
			if(flag(5).eq.0) der5(i) = (-b1(i)-b2(i)
     1					+b3(i)-b4(i))/(deter*val5)
			if(flag(6).eq.0) der6(i) = (b1(i)+b2(i)
     1					-b3(i)-b4(i))/(deter*val6)
1900		continue

	elseif(flag(1).ne.0.and.flag(5).ne.0) then

		deter = -2
		do 2000 i = 1,3
			b1(i) = a1(i) - val1*der1(i) 
			b2(i) = a2(i) - val1*der1(i)  -val5*der5(i)
			b3(i) = a3(i) 
			b4(i) = a4(i) -val5*der5(i)
2000		continue

		do 2100 i = 1,3
			if(flag(2).eq.0) der2(i) = (-b1(i)+b2(i)
     1					-b3(i)+b4(i))/(deter*val2)
			if(flag(3).eq.0) der3(i) = -2*b2(i)/
     1					(deter*val3)
			if(flag(4).eq.0) der4(i) = (-b1(i)-b2(i)
     1					+b3(i)-b4(i))/(deter*val4)
			if(flag(6).eq.0) der6(i) = (b1(i)+b2(i)
     1					-b3(i)-b4(i))/(deter*val6)
2100		continue

	elseif(flag(2).ne.0.and.flag(3).ne.0) then

		deter = 2
		do 2200 i = 1,3
			b1(i) = a1(i) - val2*der2(i)
			b2(i) = a2(i) - val3*der3(i)
			b3(i) = a3(i) - val2*der2(i) -val3*der3(i)
			b4(i) = a4(i) 
2200		continue

		do 2300 i = 1,3
			if(flag(1).eq.0) der1(i) = (b1(i)+b2(i)
     1					+b3(i)-b4(i))/(deter*val1)
			if(flag(4).eq.0) der4(i) = (b1(i)-b2(i)
     1					-b3(i)+b4(i))/(deter*val4)
			if(flag(5).eq.0) der5(i) = (-b1(i)+b2(i)
     1					-b3(i)+b4(i))/(deter*val5)
			if(flag(6).eq.0) der6(i) = 2*b3(i)/
     1					(deter*val6)
2300		continue

	elseif(flag(2).ne.0.and.flag(4).ne.0) then

		deter = -2
		do 2400 i = 1,3
			b1(i) = a1(i) - val2*der2(i) - val4*der4(i)
			b2(i) = a2(i) 
			b3(i) = a3(i) - val2*der2(i)
			b4(i) = a4(i) -val4*der4(i)
2400		continue

		do 2500 i = 1,3
			if(flag(1).eq.0) der1(i) = -2*b1(i)/
     1					(deter*val1)
			if(flag(3).eq.0) der3(i) = (b1(i)-b2(i)
     1					-b3(i)+b4(i))/(deter*val3)
			if(flag(5).eq.0) der5(i) = (b1(i)-b2(i)
     1					+b3(i)-b4(i))/(deter*val5)
			if(flag(6).eq.0) der6(i) = (-b1(i)+b2(i)
     1					-b3(i)-b4(i))/(deter*val6)
2500		continue

	elseif(flag(2).ne.0.and.flag(6).ne.0) then

		deter = 2
		do 2600 i = 1,3
			b1(i) = a1(i) - val2*der2(i)
			b2(i) = a2(i) 
			b3(i) = a3(i) - val2*der2(i)  -val6*der6(i)
			b4(i) = a4(i) -val6*der6(i)
2600		continue

		do 2700 i = 1,3
			if(flag(1).eq.0) der1(i) = (b1(i)+b2(i)
     1					-b3(i)-b4(i))/(deter*val1)
			if(flag(3).eq.0) der3(i) = 2*b3(i)/
     1					(deter*val3)
			if(flag(4).eq.0) der4(i) = (b1(i)-b2(i)
     1					+b3(i)+b4(i))/(deter*val4)
			if(flag(5).eq.0) der5(i) = (-b1(i)+b2(i)
     1					-b3(i)+b4(i))/(deter*val5)
2700		continue

	elseif(flag(3).ne.0.and.flag(5).ne.0) then

		deter = 2
		do 2800 i = 1,3
			b1(i) = a1(i)
			b2(i) = a2(i) -val3*der3(i) - val5*der5(i)
			b3(i) = a3(i) - val3*der3(i)
			b4(i) = a4(i) -val5*der5(i)
2800		continue

		do 2900 i = 1,3
			if(flag(1).eq.0) der1(i) = 2*b2(i)/
     1					(deter*val1)
			if(flag(2).eq.0) der2(i) = (b1(i)-b2(i)
     1					+b3(i)-b4(i))/(deter*val2)
			if(flag(4).eq.0) der4(i) = (b1(i)-b2(i)
     1					-b3(i)+b4(i))/(deter*val4)
			if(flag(6).eq.0) der6(i) = (-b1(i)+b2(i)
     1					+b3(i)+b4(i))/(deter*val6)
2900		continue

	elseif(flag(3).ne.0.and.flag(6).ne.0) then

		deter = 2
		do 3000 i = 1,3
			b1(i) = a1(i)
			b2(i) = a2(i) -val3*der3(i)
			b3(i) = a3(i) - val3*der3(i) - val6*der6(i)
			b4(i) = a4(i) -val6*der6(i)
3000		continue

		do 3100 i = 1,3
			if(flag(1).eq.0) der1(i) = (b1(i)+b2(i)
     1					-b3(i)-b4(i))/(deter*val1)
			if(flag(2).eq.0) der2(i) = 2*b3(i)/
     1					(deter*val2)
			if(flag(4).eq.0) der4(i) = (b1(i)-b2(i)
     1					-b3(i)+b4(i))/(deter*val4)
			if(flag(5).eq.0) der5(i) = (-b1(i)+b2(i)
     1					+b3(i)+b4(i))/(deter*val5)
3100		continue

	elseif(flag(4).ne.0.and.flag(5).ne.0) then

		deter = -2
		do 3200 i = 1,3
			b1(i) = a1(i) -val4*der4(i)
			b2(i) = a2(i) -val5*der5(i)
			b3(i) = a3(i) 
			b4(i) = a4(i) -val4*der4(i) - val5*der5(i)
3200		continue

		do 3300 i = 1,3
			if(flag(1).eq.0) der1(i) = (-b1(i)-b2(i)
     1					+b3(i)-b4(i))/(deter*val1)
			if(flag(2).eq.0) der2(i) = (-b1(i)+b2(i)
     1					-b3(i)+b4(i))/(deter*val2)
			if(flag(3).eq.0) der3(i) = (b1(i)-b2(i)
     1					-b3(i)+b4(i))/(deter*val3)
			if(flag(6).eq.0) der6(i) = -2*b4(i)/
     1					(deter*val6)
3300		continue

	endif

	return
	end
