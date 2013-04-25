c	alf_tetra.f	Version 1 11/24/2000	Patrice Koehl
c
c	This subroutine computes the radius R of the sphere orthogonal
c	to the four spheres that define a tetrahedron [A,B,C,D]
c
c	Since we are only interested at how R compares to Alpha, we do not
c	output R, rather the result of the comparison
c
c	Computation is first done in floating point; however if the
c	radius R is found to close to Alpha, we switch
c	to multiple precision integer arithmetics.
c	The package GMP is used for multiple precision (with a C wrapper)
c
#include "defines.h"
c
	subroutine alf_tetra(ia,ib,ic,id,iflag,eps,SCALE,ALPHA)
c
c	Input:
c		ia,ib,ic,id	: index of the 4 vertices of the tetrahedron
c		eps 		: cutoff value for floating point
c				  filter; if value below, switch
c				  to GMP
c		SCALE		: factor used to convert floating points
c				  to multi precision integers
c		ALPHA		: value of alpha for the alpha shape
c				  (usually 0 in Biology)
c
c	Output:
c		iflag		: 1 if tetrahedron belongs to alpha complex,
c				  0 otherwise
c
	integer	npointmax
	parameter	(npointmax=MAX_POINT)
c
	integer	i,j,k,ierr
c
	integer	ia,ib,ic,id
	integer ires,iflag
c
	real*8	Dabc,Dabd,Dacd,Dbcd
	real*8	D1,D2,D3,D4,Det
	real*8	wa,wb,wc,wd
	real*8	num,den
	real*8  eps,ALPHA,SCALE
	real*8	test
	real*8	a_xyz(3),b_xyz(3),c_xyz(3),d_xyz(3)
	real*8	a(4),b(4),c(4),d(4)
	real*8	Sab(3),Sac(3),Sad(3),Sbc(3),Sbd(3),Scd(3)
	real*8	Dab(3),Dac(3),Dad(3),Dbc(3),Dbd(3),Dcd(3)
	real*8	Sa(3),Sb(3),Sc(3),Sd(3)
	real*8	Sam1(3),Sbm1(3),Scm1(3),Sdm1(3)
	real*8	Deter(3)
	real*8  coord(3*npointmax),weight(npointmax),radius(npointmax)
c
	save
c
	common /xyz_vertex/	coord,radius,weight
c
c	first create vectors in 4D space, by adding the weights as the
c	fourth coordinate
c
	do 200 i = 1,3
		a(i) = coord(3*(ia-1)+i)
		b(i) = coord(3*(ib-1)+i)
		c(i) = coord(3*(ic-1)+i)
		d(i) = coord(3*(id-1)+i)
200	continue
	a(4) = weight(ia)
	b(4) = weight(ib)
	c(4) = weight(ic)
	d(4) = weight(id)
c
c	Perform computation in floating points; if a problem occurs,
c	switch to GMP
c
c	1. Computes all Minors Smn(i+j-2)= M(m,n,i,j) = Det | m(i)  m(j) |
c						            | n(i)  n(j) |
c	for all i in [1,2] and all j in [i+1,3]
c
	do 400 i = 1,2
		do 300 j = i+1,3
			k = i+j-2
			Sab(k) = a(i)*b(j)-a(j)*b(i)
			Sac(k) = a(i)*c(j)-a(j)*c(i)
			Sad(k) = a(i)*d(j)-a(j)*d(i)
			Sbc(k) = b(i)*c(j)-b(j)*c(i)
			Sbd(k) = b(i)*d(j)-b(j)*d(i)
			Scd(k) = c(i)*d(j)-c(j)*d(i)
300		continue
400	continue
c
c	Now compute all Minors 
c		Sq(i+j-2) = M(m,n,p,i,j,0) = Det | m(i) m(j) 1 |
c		       			         | n(i) n(j) 1 |
c						 | p(i) p(j) 1 |
c
c	and all Minors
c		Det(i+j-2) = M(m,n,p,q,i,j,4,0) = Det | m(i) m(j) m(4) 1 |
c						      | n(i) n(j) n(4) 1 |
c						      | p(i) p(j) p(4) 1 |
c						      | q(i) q(j) q(4) 1 |
c
c	m,n,p,q are the four vertices of the tetrahedron, i and j correspond
c	to two of the coordinates of the vertices, and m(4) refers to the
c	"weight" of vertices m
c
	do 500 i = 1,3
		Sa(i) = Scd(i) - Sbd(i) + Sbc(i)
		Sb(i) = Scd(i) - Sad(i) + Sac(i)
		Sc(i) = Sbd(i) - Sad(i) + Sab(i)
		Sd(i) = Sbc(i) - Sac(i) + Sab(i)
		Sam1(i) = -Sa(i)
		Sbm1(i) = -Sb(i)
		Scm1(i) = -Sc(i)
		Sdm1(i) = -Sd(i)
500	continue
c
	do 600 i = 1,3
		Deter(i) = a(4)*Sa(i)-b(4)*Sb(i)+c(4)*Sc(i)-d(4)*Sd(i)
600	continue
c
c	Now compute the determinant needed to compute the radius of the
c	sphere orthogonal to the four balls that define the tetrahedron :
c
c		D1 = Minor(a,b,c,d,4,2,3,0)
c		D2 = Minor(a,b,c,d,1,3,4,0)
c		D3 = Minor(a,b,c,d,1,2,4,0)
c		D4 = Minor(a,b,c,d,1,2,3,0)
c
	D1 = Deter(3)
	D2 = Deter(2)
	D3 = Deter(1)
	D4 = a(1)*Sa(3)-b(1)*Sb(3)+c(1)*Sc(3)-d(1)*Sd(3)
c
c	Now compute all minors:
c		Dmnp = Minor(m,n,p,1,2,3) = Det | m(1) m(2) m(3) |
c						| n(1) n(2) n(3) |
c						| p(1) p(2) p(3) |
c
	Dabc = a(1)*Sbc(3)-b(1)*Sac(3) + c(1)*Sab(3)
	Dabd = a(1)*Sbd(3)-b(1)*Sad(3) + d(1)*Sab(3)
	Dacd = a(1)*Scd(3)-c(1)*Sad(3) + d(1)*Sac(3)
	Dbcd = b(1)*Scd(3)-c(1)*Sbd(3) + d(1)*Sbc(3)
c
c	We also need :
c		Det = Det | m(1) m(2) m(3) m(4) |
c			  | n(1) n(2) n(3) n(4) |
c			  | p(1) p(2) p(3) p(4) |
c			  | q(1) q(2) q(3) q(4) |
c
	Det = -a(4)*Dbcd + b(4)*Dacd -c(4)*Dabd + d(4)*Dabc
c
c	The radius of the circumsphere of the weighted tetrahedron is then:
c
	num = (D1*D1 + D2*D2 + D3*D3 + 4*D4*Det)
	den = (4*D4*D4)
c
c	If this radius is too close to the value of ALPHA, we switch to
c	GMP
c
	test = alpha*den - num
	if(abs(test).lt.eps) then
		call tetra_radius_gmp(ia,ib,ic,id,ires,SCALE,ALPHA)
		test = ires
	endif
c
c	The spectrum for a tetrahedron is [R_t Infinity[. If ALPHA is in
c	that interval, the tetrahedron is part of the alpha shape, otherwise
c	it is discarded
c	If tetrahedron is part of the alpha shape, then the 4 triangles,
c	the 6 edges and the four vertices are also part of the alpha
c	complex
c
	iflag = 0
	if(test.gt.0) iflag = 1
c
	return
	end
c
c	alf_trig.f		Version 1 3/29/07	Patrice Koehl
c
c	This subroutine checks if a triangle belongs to the alpha complex.
c	It computes the radius of the sphere orthogonal to the three
c	balls that define the triangle; if this	radius is smaller than
c	alpha, the triangle belongs to the alpha complex.
c	We also checked if the triangle is "attached", i.e. if
c	the fourth vertex of any of the tetrahedron attached to the
c	triangle is "hidden" by the triangle
c
	subroutine alf_trig(ia,ib,ic,id,ie,irad,iattach,eps,scale,alpha)
c
c	Input:
c		ia,ib,ic	: indices of the three vertices of the triangle
c		id,ie		: indices of the two vertices such that
c			  	  (ia,ib,ic,id) and (ia,ib,ic,ie) are the
c			  	  two tetrahedra that share the triangle
c		eps 		: cutoff value for floating point
c			  	  filter; if value below, switch
c			  	  to GMP
c		SCALE		: factor used to convert floating points
c			  	  to multi precision integers
c		ALPHA		: value of alpha for the alpha shape
c			  	(usually 0 for measures of molecule)
c	Output:
c		irad		: flag; set to 1 if radius(trig) < alpha
c		iattach		: flag; set to 1 if trig is attached
c
	integer	npointmax
	parameter	(npointmax=MAX_POINT)
c
	integer	i,j
	integer	ia,ib,ic,id,ie
	integer	irad,iattach
	integer memory
c
	real*8	wa,wb,wc,wd,we
	real*8	eps,scale,alpha
c
	real*8	Dabc
	real*8	a_xyz(3),b_xyz(3),c_xyz(3),d_xyz(3),e_xyz(3)
	real*8	a(4),b(4),c(4),d(4),e(4)
	real*8	Sab(3,4),Sac(3,4),Sbc(3,4)
	real*8	S(3,4),T(2,3)
	real*8	coord(3*npointmax),radius(npointmax),weight(npointmax)
c
	logical	attach,testr
c
	common /xyz_vertex/	coord,radius,weight
c
	save
c
	iattach = 0
	irad = 0
c
	do 100 i = 1,3
		a(i) = coord(3*(ia-1)+i)
		b(i) = coord(3*(ib-1)+i)
		c(i) = coord(3*(ic-1)+i)
		d(i) = coord(3*(id-1)+i)
		if(ie.ne.0) e(i) = coord(3*(ie-1)+i)
100	continue
	a(4) = weight(ia)
	b(4) = weight(ib)
	c(4) = weight(ic)
	d(4) = weight(id)
	if(ie.ne.0) e(4) = weight(ie)
c
c	Perform computation in floating points; if a problem occurs,
c	switch to GMP
c
c	1. Computes all Minors Smn(i,j)= M(m,n,i,j)   = Det | m(i)  m(j) |
c						            | n(i)  n(j) |
c	m,n are two vertices of the triangle, i and j correspond
c	to two of the coordinates of the vertices
c
c	for all i in [1,3] and all j in [i+1,4]
c
	do 400 i = 1,3
		do 300 j = i+1,4
			Sab(i,j) = a(i)*b(j)-a(j)*b(i)
			Sac(i,j) = a(i)*c(j)-a(j)*c(i)
			Sbc(i,j) = b(i)*c(j)-b(j)*c(i)
300		continue
400	continue
c
c	Now compute all Minors 
c		S(i,j) = M(a,b,c,i,j,0)    = Det | a(i) a(j) 1 |
c		       			         | b(i) b(j) 1 |
c						 | c(i) c(j) 1 |
c
c	a,b,c are the 3 vertices of the triangle, i and j correspond
c	to two of the coordinates of the vertices
c
c	for all i in [1,3] and all j in [i+1,4]
c
	do 600 i = 1,3
		do 500 j = i+1,4
			S(i,j)=Sbc(i,j)-Sac(i,j)+Sab(i,j)
500		continue
600	continue
c
c	Now compute all Minors
c		T(i,j) = M(a,b,c,i,j,4)    = Det | a(i) a(j) a(4) |
c		       			         | b(i) b(j) b(4) |
c						 | c(i) c(j) c(4) |
c
c	for all i in [1,2] and all j in [i+1,3]
c
	do 800 i = 1,2
		do 700 j = i+1,3
			T(i,j)=a(4)*Sbc(i,j)-b(4)*Sac(i,j)+c(4)*Sab(i,j)
700		continue
800	continue
c
c	Finally,  need Dabc = M(a,b,c,1,2,3) Det | a(1) a(2) a(3) |
c		       			         | b(2) b(2) b(3) |
c						 | c(3) c(2) c(3) |
c
	Dabc = a(1)*Sbc(2,3)-b(1)*Sac(2,3) + c(1)*Sab(2,3)
c
c	First check if a,b,c attached to d:
c
	memory = 0
	call triangle_attach(ia,ib,ic,id,S,T,Dabc,d,attach,eps,memory)
c
c	If attached, we can stop there, the triangle will not be part of the
c	alpha complex
c
	if(attach) then
		iattach = 1
		return
	endif
c
c	If e exists, check if a,b,c attached to e:
c
	if(ie.ne.0) then
		call triangle_attach(ia,ib,ic,ie,S,T,Dabc,e,attach,eps,memory)
c
c		If attached, we can stop there, the triangle will not be part of the
c		alpha complex
c
		if(attach) then
			iattach = 1
			return
		endif
	endif
c
c	Now check if alpha is bigger than the radius of the sphere orthogonal
c	to the three balls at A, B, C:
c
	call triangle_radius(ia,ib,ic,S,T,Dabc,testr,alpha,eps,scale,
     1		memory)
c
	if(testr) irad = 1
c
	return
	end
c	alf_edge.f		Version 1 3/29/07	Patrice Koehl
c
c	This subroutine checks if an edge belongs to the alpha complex.
c	It computes the radius of the sphere orthogonal to the two
c	balls that define the edge; if this radius is smaller than
c	alpha, the edge belongs to the alpha complex.
c	We also checked if the edge is "attached", i.e. if
c	the third vertex of any of the triangle attached to the
c	edge is "hidden" by the edge
c
	subroutine alf_edge(ia,ib,ncheck,listcheck,irad,iattach,
     1			eps,scale,alpha)
c
c	Input:
c		ia,ib		: indices of the two vertices of the edge
c		ncheck		: number of triangles in the star of the edge
c		listcheck	: list of vertices to check
c		eps 		: cutoff value for floating point
c			  	  filter; if value below, switch
c			  	  to GMP
c		SCALE		: factor used to convert floating points
c			  	  to multi precision integers
c		ALPHA		: value of alpha for the alpha shape
c			  	(usually 0 for measures of molecule)
c	Output:
c		irad		: flag; set to 1 if radius(edge) < alpha
c		iattach		: flag; set to 1 if edge is attached
c
	integer	npointmax
	parameter	(npointmax=MAX_POINT)
c
	integer	i,j,k
	integer	ia,ib,ic
	integer	ncheck
	integer memory,irad
	integer iattach
c
	integer	listcheck(ncheck)
c
	real*8	eps,scale,alpha
	real*8	Dab(4),Sab(3),Tab(3)
	real*8	a(4),b(4),c(4)
	real*8	coord(3*npointmax),radius(npointmax),weight(npointmax)
c
	logical attach,rad
c
	common /xyz_vertex/ coord,radius,weight
c
	save
c
	iattach = 1
	irad = 0
c
c	0. define coordinates
c
	do 100 i = 1,3
		a(i) = coord(3*(ia-1)+i)
		b(i) = coord(3*(ib-1)+i)
100	continue
	a(4) = weight(ia)
	b(4) = weight(ib)
c
c	1. Compute all Minors Dab(i) = M(a,b,i,0) = Det | a(i) 1 |
c							| b(i) 1 |
c
c	for all i in [1,4]
c
	do 200 i = 1,4
		Dab(i) = a(i) - b(i)
200	continue
c
c	2. Computes all Minors Sab(i,j)= M(a,b,i,j)   = Det | a(i)  a(j) |
c						            | b(i)  b(j) |
c
	do 400 i = 1,2
		do 300 j = i+1,3
			k=i+j-2
			Sab(k) = a(i)*b(j)-b(i)*a(j)
300		continue
400	continue
c
c	3. Computes all Minors Tab(i)= M(a,b,i,4)   = Det | a(i)  a(4) |
c					                  | b(i)  b(4) |
c
	do 500 i = 1,3
		Tab(i) = a(i)*b(4) - b(i)*a(4)
500	continue
c
c	First check attachment
c
	memory = 0
	do 700 i = 1,ncheck
c
		ic = listcheck(i)
c
		do 600 j = 1,3
			c(j) = coord(3*(ic-1)+j)
600		continue
		c(4) = weight(ic)
c
		call edge_attach(ia,ib,a,b,Dab,Sab,Tab,ic,c,
     1			attach,eps,memory)
c
		if(attach) return
c
700	continue
c
	iattach = 0
c
c	Edge is not attached; check radius
c
	call edge_radius(ia,ib,a,b,Dab,Sab,Tab,rad,alpha,eps,scale,
     1		memory)
c
	if(rad) irad = 1
c
	return
	end
c
c	edge_radius.f		Version 1 11/25/2000	Patrice Koehl
c
	subroutine edge_radius(ia,ib,a,b,Dab,Sab,Tab,testr,alpha,eps,
     1		scale,memory)
c
c	This subroutine computes the radius of the smallest circumsphere to
c	an edge, and compares it to alpha.
c	For that, it needs:
c
c	Input:
c		a,b	: coordinate of the two vertices defining the edge
c		Dab	: minor(a,b,i,0) for all i=1,2,3,4
c		Sab	: minor(a,b,i,j) for i = 1,2 and j =i+1,3
c		Tab	: minor(a,b,i,4) for i = 1,2,3
c		alpha	: value of alpha considered
c		eps	: precision: if a floating point test lead to a
c			  value below this precision, computation
c			  switches to GMP
c	Ouput:
c		testr	: flag that defines if radius smaller than alpha
c
	integer	i,val,memory,ia,ib
c
	real*8	d0,d1,d2,d3,d4
	real*8	alpha,eps,scale,num,den,rho2
	real*8	r_11,r_22,r_33,r_14,r_313,r_212,diff
	real*8	a(3),b(3)
	real*8	Sab(3),Dab(4),Tab(3)
	real*8	res(0:3,1:4)
c
	logical testr
c
	save
c
	testr = .false.
c
c	Formulas have been derived by projection on 4D space,
c	which requires some precaution when some coordinates are
c	equal.
c
	res(0,4) = Dab(4)
c
	if(a(1).ne.b(1)) then
		do 100 i = 1,3
			res(0,i) = Dab(i)
			res(i,4) = Tab(i)
100		continue
		res(1,2) = Sab(1)
		res(1,3) = Sab(2)
		res(2,3) = Sab(3)
	elseif(a(2).ne.b(2)) then
		res(0,1) = Dab(2)
		res(0,2) = Dab(3)
		res(0,3) = Dab(1)
		res(1,2) = Sab(3)
		res(1,3) = -Sab(1)
		res(2,3) = -Sab(2)
		res(1,4) = Tab(2)
		res(2,4) = Tab(3)
		res(3,4) = Tab(1)
	elseif(a(3).ne.b(3)) then
		res(0,1) = Dab(3)
		res(0,2) = Dab(1)
		res(0,3) = Dab(2)
		res(1,2) = -Sab(2)
		res(1,3) = -Sab(3)
		res(2,3) = Sab(1)
		res(1,4) = Tab(3)
		res(2,4) = Tab(1)
		res(3,4) = Tab(2)
	else
		write(6,*) 'Problem in hidden1: edges defined from',
     1		' a single point'
		stop
	endif
c
	r_11 = res(0,1)*res(0,1)
	r_22 = res(0,2)*res(0,2)
	r_33 = res(0,3)*res(0,3)
	r_14 = res(0,1)*res(0,4)
	r_313 = res(0,3)*res(1,3)
	r_212 = res(0,2)*res(1,2)
	diff = res(0,3)*res(1,2)-res(0,2)*res(1,3)
c
c	First compute radius of circumsphere
c
	d0 = -2*res(0,1)*(r_11+r_22+r_33)
	d1 = res(0,1)*(2*(r_313 + r_212)-r_14)
	d2 = -2*res(1,2)*(r_11+r_33) - res(0,2)*(r_14-2*r_313)
	d3 = -2*res(1,3)*(r_11+r_22) -res(0,3)*(r_14-2*r_212)
	d4 = 2*res(0,1)*(res(0,1)*res(1,4)+res(0,2)*res(2,4)+
     1		res(0,3)*res(3,4)) +4*(res(2,3)*diff
     2		    - res(0,1)*(res(1,2)*res(1,2)+res(1,3)*res(1,3)))
c
	num = d1*d1 + d2*d2 + d3*d3 - d0*d4
	den = d0*d0
c
c	For efficiency purpose, I assume that this routine is only used to compute
c	the dual complex (i.e. alpha=0), and therefore I do not consider the denominator as
c	it is always positive)
c
c	rho2 = num/den
	rho2 = num
c
c	if(ia.eq.92.and.ib.eq.98) then
c		write(6,*) 'ia,ib,,num:',ia,ib,num
c	endif
c	if(ia.eq.98.and.ib.eq.92) then
c		write(6,*) 'ia,ib,,num:',ia,ib,num
c	endif
	if(abs(alpha*den-rho2).lt.eps) then
		call edge_radius_gmp(ia,ib,val,scale,alpha,memory)
		memory = 1
		if(val.eq.1) testr = .true.
		return
	endif
c
	if(alpha.gt.rho2) testr = .true.
c
	return
	end
c
c	edge_attach.f		Version 1 6/17/2005	Patrice Koehl
c
	subroutine edge_attach(ia,ib,a,b,Dab,Sab,Tab,ic,c,testa,
     1		eps,memory)
c
c	This subroutine check if an edge ab of a tetrahedron is "attached"
c	to a given vertex c
c	For that, it needs:
c
c	Input:
c		ia,ib:	: indices of the two vertices of the edge
c		Dab	: minor(a,b,i,0) for all i=1,2,3,4
c		Sab	: minor(a,b,i,j) for i = 1,2 and j =i+1,3
c		Tab	: minor(a,b,i,0) for all i=1,2,3
c		eps	: precision: if a floating point test lead to a
c			  value below this precision, computation
c			  switches to LIA
c	Ouput:
c		testa	: flag that defines if edge is attached or not
c
	integer	i,j,k
	integer	ia,ib,ic
	integer	val,memory
c
	real*8	eps,dtest
	real*8	r_11,r_22,r_33,diff,d0,d5
	real*8	Sab(3),Dab(4),Tab(3),Sc(3),Tc(3)
	real*8	a(4),b(4),c(4)
	real*8	res(0:3,1:3),res2_c(3,4)
c
	logical testa
c
	save
c
	testa = .false.
c
c	Need to compute:
c	Sc	: minor(a,b,c,i,j,0) for i=1,2 and j = i+1,3
c	Tc	: minor(a,b,c,i,4,0) for i = 1,2,3
c
	Do 200 i = 1,2
		do 100 j = i+1,3
			k=i+j-2
			Sc(k)=c(i)*Dab(j)-c(j)*Dab(i)+Sab(k)
100		continue
200	continue
c
	do 300 i = 1,3
		Tc(i) = c(i)*Dab(4)-c(4)*Dab(i) + Tab(i)
300	continue
c
c	Formulas have been derived by projection on 4D space,
c	which requires some precaution when some coordinates are
c	equal.
c
	if(a(1).ne.b(1)) then
		do 400 i = 1,3
			res(0,i) = Dab(i)
			res2_c(i,4) = Tc(i)
400		continue
		res(1,2) = Sab(1)
		res(1,3) = Sab(2)
		res(2,3) = Sab(3)
		res2_c(1,2) = Sc(1)
		res2_c(1,3) = Sc(2)
		res2_c(2,3) = Sc(3)
	elseif(a(2).ne.b(2)) then
		res(0,1) = Dab(2)
		res(0,2) = Dab(3)
		res(0,3) = Dab(1)
		res(1,2) = Sab(3)
		res(1,3) = -Sab(1)
		res(2,3) = -Sab(2)
		res2_c(1,2) = Sc(3)
		res2_c(1,3) = -Sc(1)
		res2_c(2,3) = -Sc(2)
		res2_c(1,4) = Tc(2)
		res2_c(2,4) = Tc(3)
		res2_c(3,4) = Tc(1)
	elseif(a(3).ne.b(3)) then
		res(0,1) = Dab(3)
		res(0,2) = Dab(1)
		res(0,3) = Dab(2)
		res(1,2) = -Sab(2)
		res(1,3) = -Sab(3)
		res(2,3) = Sab(1)
		res2_c(1,2) = -Sc(2)
		res2_c(1,3) = -Sc(3)
		res2_c(2,3) = Sc(1)
		res2_c(1,4) = Tc(3)
		res2_c(2,4) = Tc(1)
		res2_c(3,4) = Tc(2)
	else
		write(6,*) 'Problem in hidden1: edges defined from',
     1		' a single point'
		stop
	endif
c
	r_11 = res(0,1)*res(0,1)
	r_22 = res(0,2)*res(0,2)
	r_33 = res(0,3)*res(0,3)
	diff = res(0,3)*res(1,2)-res(0,2)*res(1,3)
c
c	Check attachement with vertex C
c
	d0 = -2*res(0,1)*(r_11+r_22+r_33)
c
	d5 = res(0,1)*(res(0,1)*res2_c(1,4)+res(0,2)*res2_c(2,4)
     1		+res(0,3)*res2_c(3,4)-2*(res(1,3)*res2_c(1,3)
     2		+res(1,2)*res2_c(1,2))) + 2*res2_c(2,3)*diff
c
	dtest = d0*d5
c
c	if(ia.eq.92.and.ib.eq.98) then
c		write(6,*) 'ia,ib,ic,dtest:',ia,ib,ic,dtest
c	endif
c	if(ia.eq.98.and.ib.eq.92) then
c		write(6,*) 'ia,ib,ic,dtest:',ia,ib,ic,dtest
c	endif
	if(abs(dtest).lt.eps) then
		call edge_attach_gmp(ia,ib,ic,val,memory)
		memory = 1
		if(val.eq.1) testa = .true.
		return
	endif
c
c	If no problem, set testa to true if t < 0
c
	if(dtest.lt.0) testa = .true.
c
	return
	end
c
c	triangle_attach.f	Version 1 3/30/2007	Patrice Koehl
c
	subroutine triangle_attach(ia,ib,ic,id,S,T,Dabc,d,testa,eps,
     1				memory)
c
c	Input:
c
c	For the three points a,b,c that form the triangles, the program
c	needs as input the following determinants:
c
c	  S(i,j) = Minor(a,b,c,i,j,0)= det | a(i)  a(j)  1 |
c					   | b(i)  b(j)  1 |
c					   | c(i)  c(j)  1 |
c	for all i in [1,3], j in [i+1,4]
c
c	T(i,j) = M(a,b,c,i,j,4)    = Det | a(i) a(j) a(4) |
c		       		         | b(i) b(j) b(4) |
c					 | c(i) c(j) c(4) |
c
c	for all i in [1,2] and all j in [i+1,3]
c
c	Dabc = Det | a(1) a(2) a(3) |
c		   | b(1) b(2) b(3) |
c		   | c(1) c(2) c(3) |
c
c	and the coordinates of the fourth vertex d
c
c	Output:
c
c	testa	: flag set to 1 if the fourth point d 
c		  is inside the circumsphere of {a,b,c}
c
c	The program tests for problem with floating points (i.e. some
c	results smaller than EPS, in which case the sign of the
c	expression cannot be defined). If problems, IERR returns as 1,
c	and the program will then switch to GMP
c
	integer	iattach,ia,ib,ic,id,val
	integer	memory
c
	real*8	test,eps,sum2
	real*8	Dabc
	real*8	Det1, Det2, Det3, Deter
	real*8	S(3,4),T(2,3),d(4)
c
	logical testa
c
	save
c
	testa = .false.
c
c	We need to compute:
c
c	Det1 = Minor(a,b,c,d,2,3,4,0)
c	Det2 = Minor(a,b,c,d,1,3,4,0)
c	Det3 = Minor(a,b,c,d,1,2,4,0)
c	Deter= Minor(a,b,c,d,1,2,3,0)
c
	Det1  = - d(2)*S(3,4)+d(3)*S(2,4)-d(4)*S(2,3)+T(2,3)
	Det2  = - d(1)*S(3,4)+d(3)*S(1,4)-d(4)*S(1,3)+T(1,3)
	Det3  = - d(1)*S(2,4)+d(2)*S(1,4)-d(4)*S(1,2)+T(1,2)
	Deter = -d(1)*S(2,3) +d(2)*S(1,3)-d(3)*S(1,2)+Dabc
c
c	sums2 = S(1,2)*S(1,2)+S(1,3)*S(1,3)+S(2,3)*S(2,3)
c
c	check if the face is "attached" to the fourth vertex of the
c	parent tetrahedron
c	
c	test =  sums2 *(Det1*S(2,3)+Det2*S(1,3)+Det3*S(1,2)-2*Deter*Dabc) 
	test =  Det1*S(2,3)+Det2*S(1,3)+Det3*S(1,2)-2*Deter*Dabc 
c
c	Check for problems, in which case should be GMP
c
	if(abs(test).lt.eps) then
		call triangle_attach_gmp(ia,ib,ic,id,val,memory)
		memory = 1
		if(val.eq.1) testa = .true.
		return
	endif
c
c	If no problem, set testa to true if test > 0
c
	if(test.gt.0) then
		testa = .TRUE.
	endif
c
	return
	end
c
c	triangle_radius.f	Version 1 6/17/2005	Patrice Koehl
c
c	This subroutine computes the radius of the smallest circumsphere to
c	a triangle
c
	subroutine triangle_radius(ia,ib,ic,S,T,Dabc,testr,alpha,eps,
     1		scale,memory)
c
c	Input:
c
c	For the three points a,b,c that form the triangles, the program
c	needs as input the following determinants:
c
c	S(i,j)   = Minor(a,b,c,i,j,0)= det | a(i)  a(j)  1 |
c					   | b(i)  b(j)  1 |
c					   | c(i)  c(j)  1 |
c
c	for i in [1,3] and j in [i+1,4]
c
c	and:
c
c	T(i,j) = Minor(a,b,c,i,j,4)=det | a(i) a(j) a(4) |
c					| b(i) b(j) b(4) |
c					| c(i) c(j) c(4) |
c
c	and
c
c	Dabc  = Minor(a,b,c,1,2,3)
c
c	Output:
c
c	testr	: flag set to 1 if ALPHA is larger than rho, the radius
c		  of the circumsphere of the triangle
c
c	The program tests for problem with floating points (i.e. some
c	results smaller than EPS, in which case the sign of the
c	expression cannot be defined). If problems, IERR returns as 1,
c	and the program will then switch to GMP
c
	integer	ia,ib,ic,val
	integer	memory
c
	real*8	Dabc
	real*8	d0,d1,d2,d3,d4
	real*8	eps,alpha,scale
	real*8	sums2,num
	real*8	S(3,4),T(2,3)
c
	logical testr
c
	save
c
	testr = .false.
c
	sums2 = S(1,2)*S(1,2)+S(1,3)*S(1,3)+S(2,3)*S(2,3)
c
	d0 = sums2
c
	d1 = S(1,3)*S(3,4)+S(1,2)*S(2,4) -2*Dabc*S(2,3)
	d2 = S(1,2)*S(1,4)-S(2,3)*S(3,4) -2*Dabc*S(1,3)
	d3 = S(2,3)*S(2,4)+S(1,3)*S(1,4) +2*Dabc*S(1,2)
	d4 = S(1,2)*T(1,2)+S(1,3)*T(1,3) + S(2,3)*T(2,3) - 
     1		2*Dabc*Dabc
c
	num  = 4*(d1*d1+d2*d2+d3*d3) + 16*d0*d4
c
	if(abs(alpha-num).lt.eps) then
		call triangle_radius_gmp(ia,ib,ic,val,scale,alpha,memory)
		memory = 1
		if(val.eq.1) testr = .true.
		return
	endif
c
	if(alpha.gt.num) testr = .true.
c
	return
	end
c
c	vertex_attach.f		Version 1 6/17/2005	Patrice Koehl
c
c	This subroutine tests if a vertex is attached to another vertex.
c	The computation is done both way.
c
c       Let S be a simplex, and y_S the center of the ball orthogonal
c       to all balls in S. A point p is attached to S iff
c       pi(y_S, p) < 0, where pi is the power distance between the two
c       weighted points y_S and p.
c
c       Let S = {a}, with a of weight ra**2. Then y_S is the ball centered
c       at a, but with weight -ra**2.
c       The power distance between y_S and a point b is:
c
c       pi(y_S, b) = dist(a,b)**2 +ra**2 -rb**2

	subroutine vertex_attach(ia,ib,testa,testb,eps)
c
	integer	npointmax
c
	parameter (npointmax=MAX_POINT)
c
	integer	i,ia,ib,tst1,tst2
c
	real*8	a,b,ra,rb,ra2,rb2,dist2,eps,test1,test2
	real*8	Dab(3)
	real*8	coord(3*npointmax),radius(npointmax),weight(npointmax)
c
	logical testa,testb
c
	common /xyz_vertex/ coord,radius,weight
c
	save
c
	testa = .false.
	testb = .false.
c
	do 100 i = 1,3
		a = coord(3*(ia-1)+i)
		b = coord(3*(ib-1)+i)
		Dab(i) = a - b
100	continue
	ra = radius(ia)
	rb = radius(ib)
c
	ra2 = ra*ra
	rb2 = rb*rb
c
	dist2 = Dab(1)*Dab(1) + Dab(2)*Dab(2) + Dab(3)*Dab(3)
c
	test1 = dist2 + ra2 - rb2
	test2 = dist2 - ra2 + rb2
c
	if(abs(test1).lt.eps.or.abs(test2).lt.eps) then
		call vertex_attach_gmp(ia,ib,tst1,tst2)
		if(tst1.eq.1) testa = .true.
		if(tst2.eq.1) testb = .true.
		return
	endif
c
	if(test1.lt.0) testa = .true.
	if(test2.lt.0) testb = .true.
c
	return
	end
