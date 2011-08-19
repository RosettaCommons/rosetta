c	Write_simplices.f
c
c	This subroutine writes all simplices in the output file
c
c	Copyright (C) 2002 Patrice Koehl
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
	subroutine write_del(fname)
c
	integer	ntetra_max
	parameter (ntetra_max=MAX_TETRA)
c
	integer	i,j,idx,iorient
	integer	ntetra,nkeep
c
	integer*1 tetra_info(ntetra_max)
	integer*1 tetra_nindex(ntetra_max)
c
	integer tetra(4,ntetra_max)
	integer tetra_neighbour(4,ntetra_max)
c
	character	fname*100
c
	common	/tetra_zone/	ntetra,tetra,tetra_neighbour
	common	/tetra_stat/	tetra_info,tetra_nindex
c
1	format('Tetra  ',i6,'.',2x,4(i5,2x),i2,2x,5(i5,2x))
5	format('REMARK',/,
     1	       'REMARK  Tetrahedra in the Regular triangulation',/,
     2	       'REMARK',/,
     3         'REMARK  Number of tetrahedra :', i6,/,
     4	       'REMARK',/,
     5	       'REMARK  i,j,k,l             : the four vertices',
     s	       ' of the tetrahedron',/,
     s	       'REMARK  ccw                 : 1 if positive',
     s	       ' orientation, -1 otherwise',/,
     s	       'REMARK',/,
     s         'REMARK    #        i      j      k      l  ccw  ',/,
     s	       'REMARK')
c
	open(unit=1,file=fname,status='unknown')
c
c	First writes all tetrahedra
c
	nkeep = 0
	do 100 i = 1,ntetra
		if(btest(tetra_info(i),1)) then
			nkeep = nkeep + 1
		endif
100	continue
c
	write(1,5) nkeep
c
	nkeep = 0
c
	do 200 i = 1,ntetra
c
		if(.not.btest(tetra_info(i),1)) goto 200
c
		if(btest(tetra_info(i),0)) then
			iorient = 1
		else
			iorient = -1
		endif
		nkeep = nkeep + 1
		write(1,1) nkeep,(tetra(j,i)-4,j=1,4),
     1		iorient
c
200	continue
c
	close(unit=1)
c
	write(6,*) ' '
	write(6,*) 'In Delaunay:'
	write(6,*) 'Number of tetrahedron :',nkeep
	write(6,*) ' '
c
	return
	end
c
	subroutine write_alpha(fname)
c
	integer	ntetra_max,npointmax
	parameter (ntetra_max=MAX_TETRA)
	parameter (npointmax=MAX_POINT)
c
	integer	i,j,k,l,idx
	integer	ia,ib,ic
	integer	ntetra,nkeep
	integer	ntrig,nedge
	integer	iorient,itrig,iedge
	integer	npoints,nvertex,nvert,npass,ipair
	integer	itetra,jtetra,ktetra
	integer	trig1,trig2,trig_in,trig_out,triga,trigb
c
	integer	other3(3,4)
	integer	face_info(2,6),face_pos(2,6)
	integer	pair(2,6)
c
	integer*1 ival
c
	integer*1 tetra_info(ntetra_max)
	integer*1 tetra_nindex(ntetra_max)
	integer*1 tetra_edge(ntetra_max)
	integer*1 tetra_mask(ntetra_max)
c
	integer tetra(4,ntetra_max)
	integer tetra_neighbour(4,ntetra_max)
c
	integer*1	vertex_info(npointmax)
c
	character	fname*100
c
	common	/tetra_zone/	ntetra,tetra,tetra_neighbour
	common	/tetra_stat/	tetra_info,tetra_nindex
	common  /alp_zone/	tetra_edge
	common  /vertex_zone/	npoints,nvertex,vertex_info
c
	data other3 /2,3,4,1,3,4,1,2,4,1,2,3/
	data face_info/1,2,1,3,1,4,2,3,2,4,3,4/
	data face_pos/2,1,3,1,4,1,3,2,4,2,4,3/
	data pair/3,4,2,4,2,3,1,4,1,3,1,2/
c
	save
c
1	format('Tetra  ',i6,'.',2x,4(i5,2x),i2,2x,5(i5,2x))
2	format('Trig   ',i6,'.',2x,8(i5,2x))
3	format('Edge   ',i6,'.',2x,3(i5,2x))
4	format('Vertex ',i6,'.',2x,2(i5,2x))
5	format('REMARK',/,
     1	       'REMARK  Tetrahedra in the Alpha Complex',/,
     2	       'REMARK',/,
     3         'REMARK  Number of tetrahedra :', i6,/,
     4	       'REMARK',/,
     5	       'REMARK  i,j,k,l             : the four vertices',
     s	       ' of the tetrahedron',/,
     s	       'REMARK  ccw                 : 1 if positive',
     s	       ' orientation, -1 otherwise',/,
     s	       'REMARK',/,
     s         'REMARK    #        i      j      k      l  ccw  ',/,
     s	       'REMARK')
6	format('REMARK',/,
     1	       'REMARK  Triangles in the Alpha Complex',/,
     2	       'REMARK',/,
     3         'REMARK  Number of triangles :', i6,/,
     4	       'REMARK',/,
     5	       'REMARK  i,j,k             : the three vertices ',
     s	       'of the triangle',/,
     s	       'REMARK',/,
     s         'REMARK    #        i      j      k   ',/,
     s	       'REMARK')
7	format('REMARK',/,
     1	       'REMARK  Edges in the Alpha Complex',/,
     2	       'REMARK',/,
     3         'REMARK  Number of edges :', i6,/,
     4	       'REMARK',/,
     5	       'REMARK  i,j   : the two vertices ',
     s	       'of the edge',/,
     s	       'REMARK',/,
     s         'REMARK    #        i      j    ',/,
     s	       'REMARK')
8	format('REMARK',/,
     1	       'REMARK  Vertices in the Alpha Complex',/,
     2	       'REMARK',/,
     3         'REMARK  Number of vertices :', i6,/,
     4	       'REMARK',/,
     5	       'REMARK  i   : index of vertex ',
     s	       'in original list of points',/,
     s	       'REMARK',/,
     s         'REMARK    #        i      ',/,
     s	       'REMARK')
c
	open(unit=1,file=fname,status='unknown')
c
c	First writes all tetrahedra
c
	nkeep = 0
	do 100 i = 1,ntetra
		if(btest(tetra_info(i),7)) then
			nkeep = nkeep + 1
		endif
100	continue
c
	write(1,5) nkeep
c
	nkeep = 0
c
	do 200 i = 1,ntetra
c
		if(.not.btest(tetra_info(i),7)) goto 200
c
		if(btest(tetra_info(i),0)) then
			iorient = 1
		else
			iorient = -1
		endif
		nkeep = nkeep + 1
		write(1,1) nkeep,(tetra(j,i)-4,j=1,4),
     1		iorient
c
200	continue
c
c	Now writes all triangles
c
	ntrig = 0
	do 400 i = 1,ntetra
c
		if(.not.btest(tetra_info(i),1)) goto 400
		do 300 j = 1,4
			jtetra = tetra_neighbour(j,i) 
			if(jtetra.eq.0.or.jtetra.gt.i) then
				if(btest(tetra_info(i),2+j)) 
     1				ntrig = ntrig + 1
			endif
300		continue
c
400	continue
c
	write(1,6) ntrig
c
	ntrig = 0
	do 600 idx = 1,ntetra
c
		if(.not.btest(tetra_info(idx),1)) goto 600
c
		i = tetra(1,idx)
		j = tetra(2,idx)
		k = tetra(3,idx)
		l = tetra(4,idx)
c
		do 500 itrig = 1,4
c
			jtetra = tetra_neighbour(itrig,idx) 
			if(jtetra.eq.0.or.jtetra.gt.idx) then
c
				if(btest(tetra_info(idx),2+itrig)) then
c
					ntrig = ntrig + 1
					if(itrig.eq.1) then
						ia = j
						ib = k
						ic = l
					elseif(itrig.eq.2) then
						ia = i
						ib = k
						ic = l
					elseif(itrig.eq.3) then
						ia = i
						ib = j
						ic = l
					else
						ia = i
						ib = j
						ic = k
					endif
c
					write(1,2) ntrig,ia-4,ib-4,ic-4
				endif
c
			endif
c
500		continue
c
600	continue
c
c	Now write all edges
c
	do 900 i = 1,ntetra
		tetra_mask(i) = 0
900	continue
c
	nedge = 0
	do 1400 itetra = 1,ntetra
c
		if(.not.btest(tetra_info(itetra),1)) goto 1400
c
		do 1300 iedge = 1,6
c
			if(btest(tetra_mask(itetra),iedge-1)) goto 1300
c
			i = tetra(pair(1,iedge),itetra)
			j = tetra(pair(2,iedge),itetra)
c
			trig1 = face_info(1,iedge)
			trig2 = face_info(2,iedge)	
c
			ktetra = itetra
			npass = 1
			trig_out = trig1
			jtetra = tetra_neighbour(trig_out,ktetra)
c
1000			continue
c
			if(jtetra.eq.0) goto 1100
c
			if(jtetra.eq.itetra) goto 1200
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
			tetra_mask(jtetra) = ibset(tetra_mask(jtetra),
     1			ipair-1)
c
			ival = ibits(tetra_nindex(ktetra),
     1				2*(trig_out-1),2)
			trig_in = ival + 1
c
			triga = face_info(1,ipair)
			trigb = face_info(2,ipair)
c
			trig_out = triga
			if(trig_in.eq.triga) then
				trig_out = trigb
			endif
c
			ktetra = jtetra
			jtetra = tetra_neighbour(trig_out,ktetra)
c
			goto 1000
c
1100			continue
c
			if(npass.eq.2) goto 1200
			npass = npass + 1
			ktetra = itetra
			trig_out = trig2
			jtetra = tetra_neighbour(trig_out,ktetra)
			goto 1000
c
1200			continue
c
			if(btest(tetra_edge(itetra),iedge-1)) then
				nedge = nedge + 1
			endif
1300		continue
1400	continue
c
	write(1,7) nedge
c
	do 1500 i = 1,ntetra
		tetra_mask(i) = 0
1500	continue
c
	nedge = 0
	do 2000 itetra = 1,ntetra
c
		if(.not.btest(tetra_info(itetra),1)) goto 2000
c
		do 1900 iedge = 1,6
c
			if(btest(tetra_mask(itetra),iedge-1)) goto 1900
c
			i = tetra(pair(1,iedge),itetra)
			j = tetra(pair(2,iedge),itetra)
c
			trig1 = face_info(1,iedge)
			trig2 = face_info(2,iedge)	
c
			ktetra = itetra
			npass = 1
			trig_out = trig1
			jtetra = tetra_neighbour(trig_out,ktetra)
c
1600			continue
c
			if(jtetra.eq.0) goto 1700
c
			if(jtetra.eq.itetra) goto 1800
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
			tetra_mask(jtetra) = ibset(tetra_mask(jtetra),
     1			ipair-1)
c
			ival = ibits(tetra_nindex(ktetra),
     1				2*(trig_out-1),2)
			trig_in = ival + 1
c
			triga = face_info(1,ipair)
			trigb = face_info(2,ipair)
c
			trig_out = triga
			if(trig_in.eq.triga) then
				trig_out = trigb
			endif
c
			ktetra = jtetra
			jtetra = tetra_neighbour(trig_out,ktetra)
c
			goto 1600
c
1700			continue
c
			if(npass.eq.2) goto 1800
			npass = npass + 1
			ktetra = itetra
			trig_out = trig2
			jtetra = tetra_neighbour(trig_out,ktetra)
			goto 1600
c
1800			continue
c
			if(btest(tetra_edge(itetra),iedge-1)) then
				nedge = nedge + 1
				write(1,3) nedge,i-4,j-4
			endif
1900		continue
2000	continue
c
c	Now write all non redundant vertices
c
	nvert = 0
	do 2100 i = 1,nvertex
		if(btest(vertex_info(i),7)) nvert = nvert + 1
2100	continue
c
	write(1,8) nvert
c
	nvert = 0
	do 2200 i = 1,nvertex
		if(.not.btest(vertex_info(i),7)) goto 2200
		nvert = nvert + 1
		write(1,4) nvert,i-4
2200	continue
c
	close(unit=1)
c
	write(6,*) ' '
	write(6,*) 'In alpha complex:'
	write(6,*) 'Number of tetrahedron :',nkeep
	write(6,*) 'Number of triangles   :',ntrig
	write(6,*) 'Number of edges       :',nedge
	write(6,*) 'Number of vertices    :',nvert
	write(6,*) ' '
c
	return
	end
