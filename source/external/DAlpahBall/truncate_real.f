C	Truncate_real.f		Version 1 10/15/2002	Patrice Koehl
c
c	This small program truncates a real number to a given accuracy
c	level (i.e. number of digits after the decimal point)
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
	subroutine truncate_real(x_in,x_out,ndigit)
c
c	Input:
c		x_in	: real number before truncation
c		ndigit	: number of digits to be kept
c	Output:
c		x_out	: real number after truncation
c
	real*8	x_in,x_out,y
	real*8	fact
c
	integer	i,mantissa
	integer	ndigit
	integer	digit(16)
c
	mantissa = int(x_in)
c
	y = x_in - mantissa
c
	x_out = mantissa
c
	fact = 1
	do 100 i = 1,ndigit
		fact = fact*10.d0
		digit(i) = nint(y*10)
		y = 10*(y-digit(i)/10.d0)
		x_out = x_out + digit(i)/fact
100	continue
c
	return
	end
