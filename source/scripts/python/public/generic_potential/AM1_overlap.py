#!/usr/bin/env python3
#
# (c) Copyright Rosetta Commons Member Institutions.
# (c) This file is part of the Rosetta software suite and is made available under license.
# (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
# (c) For more information, see http://www.rosettacommons.org. Questions about this can be
# (c) addressed to University of Washington CoMotion, email: license@uw.edu.
'''
Helper function for AM1 charge calculation.

Author: Frank DiMaio 
'''
import os, sys
import numpy
import math
from numba import jit

FACT = numpy.array([
    1, 1, 2, 6, 24, 120, 720, 5040, 40320,
    362880, 3628800, 39916800, 479001600,
    6227020800, 87178291200, 1307674368000,
    20922789888000, 355687428096000, 6402373705728000,
    121645100408832000, 2432902008176640000], dtype='float64')

# overlap integrals for elements up to atno 18
#        orbitals: 1S, 2S, 3S, 2P, 3P
# (diat2 from MOPAC)
@jit(nopython=True)
def overlap_params(s1,s2,ni,nj,rab,ii):
    def aint(x,k):
        a = numpy.zeros(k)
        a[0] = numpy.exp(-x) / x
        for i in range(1,k):
            a[i] = (a[i-1]*i+numpy.exp(-x)) / x
        return a

    def bint(x,k):
        absx = abs(x)
        b = numpy.zeros(k)
        if (
            (absx > 3.) or
            ((absx > 2.) and (k<=10)) or
            ((absx > 1.) and (k<=7)) or
            ((absx > 0.5) and (k<=5))
        ):
            expx = numpy.exp(x);
            expmx = 1. / expx;
            b[0] = (expx - expmx) / x;
            for i in range(1,k):
                sgn = 1 if (i%2==0) else -1
                b[i] = (i*b[i-1] + sgn * expx - expmx) / x;
        elif (absx <= 1e-6):
            for i in range(k):
                b[i] = (2*((i+1)%2)) / (i+1.0)
        else:
            if (absx > 2):
                last = 16
            elif (absx > 1):
                last = 13
            elif (absx > 0.5):
                last = 8
            else:
                last = 7

            for i in range(k):
                y=0
                for j in range(last):
                   xf = FACT[j]
                   y += numpy.power(-x,j) * (2*((i+j+1)%2)) / (xf*(i+j+1.));
                b[i] = y;
        return b

    if (ni > nj):
        idx = 1
        sa = s2
        sb = s1
    else:
        idx = 2
        sa = s1
        sb = s2

    j = ii + 2
    if (ii > 3):
        j-= 1

    alpha = rab * .5 * (sa + sb)
    beta = rab * .5 * (sb - sa)
    a = aint(alpha, j)
    b = bint(beta, j)
    return (idx,sa,sb,alpha,beta,a,b)


# 'diat2' from MOPAC
@jit(nopython=True)
def inline_Sij_ref(ni,esa,epa,rab,nj,esb,epb):
    # S holds 5 nonzero elts:
    #  1) sb/sa
    #  2) pb/sa
    #  3) sb/pa
    #  4) pb/pa (sigma)
    #  5) pb/pa (pi)
    Sij = numpy.zeros((5))

    # figure out interaction type (ni,nj==valence shell)
    nmin = min(ni,nj)
    nmax = max(ni,nj)
    if (nmin<2 and nmax<2):
        ii=1 #1st:1st
    elif (nmin<2 and nmax<10):
        ii=2 #1st:2nd
    elif (nmin<2 and nmax<18):
        ii=3 #1st:3rd
    elif (nmin<10 and nmax<10):
        ii=4 #2nd:2nd
    elif (nmin<10 and nmax<18):
        ii=5 #2nd:3rd
    else:  #(nmin <18 and nmax <18):
        ii=6 #3rd:3rd

    (idx,sa,sb,alpha,beta,a,b) = overlap_params ( esa,esb,ni,nj,rab,ii )

    if (ii==1): #1st:1st
        Sij[0]=0.25*numpy.sqrt((sa*sb*rab*rab)**3) * (a[2]*b[0]-b[2]*a[0])
    elif (ii==2): #1st:2nd
        W=0.125*numpy.sqrt((sa**3)*(sb**5))*(rab**4)
        Sij[0] = W * numpy.sqrt(1.0/3.0) * (
            a[3]*b[0]-b[3]*a[0]+a[2]*b[1]-b[2]*a[1])
        if (ni > 1):
            (idx,sa,sb,alpha,beta,a,b) = overlap_params ( epa,esb,ni,nj,rab,ii )
        else:
            (idx,sa,sb,alpha,beta,a,b) = overlap_params ( esa,epb,ni,nj,rab,ii )
        W=0.125*numpy.sqrt((sa**3)*(sb**5))*(rab**4)
        Sij[idx]= W * (
            a[2]*b[0]-b[2]*a[0]+a[3]*b[1]-b[3]*a[1])
    elif (ii==3):
        W=0.0625*numpy.sqrt((sa**3)*(sb**7)/7.5)*(rab**5)
        Sij[0] = W * numpy.sqrt(1.0/3.0) * (
            a[4]*b[0]-b[4]*a[0]+12.0*(a[3]*b[1]-b[3]*a[1]))
        if (ni > 1):
            (idx,sa,sb,alpha,beta,a,b) = overlap_params ( epa,esb,ni,nj,rab,ii )
        else:
            (idx,sa,sb,alpha,beta,a,b) = overlap_params ( esa,epb,ni,nj,rab,ii )
        W=0.0625*numpy.sqrt((sa**3)*(sb**7)/7.5)*(rab**5)
        Sij[idx]= W * (
            a[3]*(b[0]+b[2])-b[3]*(a[0]+a[2])+b[1]*(a[2]+a[4])-a[1]*(b[2]+b[4]))
    elif (ii==4):
        W=0.0625*numpy.sqrt((sa*sb)**5)*(rab**5)
        Sij[0] = W / 3.0 * (
            a[4]*b[0]+b[4]*a[0]-2.0*a[2]*b[2])
        if (ni>nj):
            (idx,sa,sb,alpha,beta,a,b) = overlap_params ( epa,esb,ni,nj,rab,ii )
        else:
            (idx,sa,sb,alpha,beta,a,b) = overlap_params ( esa,epb,ni,nj,rab,ii )
        W=0.0625*numpy.sqrt((sa*sb)**5)*(rab**5)
        D=a[3]*(b[0]-b[2])-a[1]*(b[2]-b[4])
        E=b[3]*(a[0]-a[2])-b[1]*(a[2]-a[4])
        Sij[idx]=W * numpy.sqrt(1.0/3.0) * (D+E)
        if (ni>nj):
            (idx,sa,sb,alpha,beta,a,b) = overlap_params ( esa,epb,ni,nj,rab,ii )
        else:
            (idx,sa,sb,alpha,beta,a,b) = overlap_params ( epa,esb,ni,nj,rab,ii )
        W=0.0625*numpy.sqrt((sa*sb)**5)*(rab**5)
        D=a[3]*(b[0]-b[2])-a[1]*(b[2]-b[4])
        E=b[3]*(a[0]-a[2])-b[1]*(a[2]-a[4])
        Sij[3-idx]=-W * numpy.sqrt(1.0/3.0) * (E-D)
        (idx,sa,sb,alpha,beta,a,b) = overlap_params ( epa,epb,ni,nj,rab,ii )
        W=0.0625*numpy.sqrt((sa*sb)**5)*(rab**5)
        Sij[3] = -W*(b[2]*(a[4]+a[0])-a[2]*(b[4]+b[0]))
        Sij[4] = 0.5*W*(
            a[4]*(b[0]-b[2])-b[4]*(a[0]-a[2])-a[2]*b[0]+b[2]*a[0])

    elif (ii==5):
        W=0.03125*numpy.sqrt((sa**5)*(sb**7)/7.5)*(rab**6)
        Sij[0]=W * numpy.sqrt(1.0/27.0) * (
            a[5]*b[0]+a[4]*b[1]-2*(a[3]*b[2]+a[2]*b[3])+a[1]*b[4]+a[0]*b[5])
        if (ni>nj):
            (idx,sa,sb,alpha,beta,a,b) = overlap_params ( epa,esb,ni,nj,rab,ii )
        else:
            (idx,sa,sb,alpha,beta,a,b) = overlap_params ( esa,epb,ni,nj,rab,ii )
        W=0.03125*numpy.sqrt((sa**5)*(sb**7)/7.5)*(rab**6)
        Sij[idx]=W * numpy.sqrt(1.0/3.0) * (
            a[5]*b[1]+a[4]*b[0]-2.0*(a[3]*b[3]+a[2]*b[2])+a[1]*b[5]+a[0]*b[4])
        if (ni>nj):
            (idx,sa,sb,alpha,beta,a,b) = overlap_params ( esa,epb,ni,nj,rab,ii )
        else:
            (idx,sa,sb,alpha,beta,a,b) = overlap_params ( epa,esb,ni,nj,rab,ii )
        W=0.03125*numpy.sqrt((sa**5)*(sb**7)/7.5)*(rab**6)
        Sij[3-idx]=-W * numpy.sqrt(1.0/3.0) * (
            a[4]*(2.0*b[3]-b[0])-b[4]*(2.0*a[2]-a[0])-a[1]
            *(b[5]-2.0*b[3])+b[1]*(a[5]-2.0*a[3]))
        (idx,sa,sb,alpha,beta,a,b) = overlap_params ( epa,epb,ni,nj,rab,ii )
        W=0.03125*numpy.sqrt((sa**5)*(sb**7)/7.5)*(rab**6)
        Sij[3]=-W*(
            b[3]*(a[0]+a[4])-a[3]*(b[0]+b[4])+b[2]*(a[1]+a[5])-a[2]*(b[1]+b[5]))
        Sij[4]=0.5*W*(
            a[5]*(b[0]-b[2])-b[5]*(a[0]-a[2])+a[4]*(b[1]-b[3])-b[4]*(a[1]-a[3])
            -a[3]*b[0]+b[3]*a[0]-a[2]*b[1]+b[2]*a[1])
    else: #ii==6
        W=numpy.sqrt((sa*sb*rab*rab)**7)/480.0
        Sij[0] = W/3.0 * (
            a[6]*b[0]-3.0*(a[4]*b[2]-a[2]*b[4])-a[0]*b[6] )
        if (ni>nj):
            (idx,sa,sb,alpha,beta,a,b) = overlap_params ( epa,esb,ni,nj,rab,ii )
        else:
            (idx,sa,sb,alpha,beta,a,b) = overlap_params ( esa,epb,ni,nj,rab,ii )
        W=numpy.sqrt((sa*sb*rab*rab)**7)/480.0
        D=a[5]*(b[0]-b[2])-2.0*a[3]*(b[2]-b[4])+a[1]*(b[4]-b[6])
        E=b[5]*(a[0]-a[2])-2.0*b[3]*(a[2]-a[4])+b[1]*(a[4]-a[6])
        Sij[idx]=W * numpy.sqrt(1.0/3.0) * (D-E)
        if (ni>nj):
            (idx,sa,sb,alpha,beta,a,b) = overlap_params ( esa,epb,ni,nj,rab,ii )
        else:
            (idx,sa,sb,alpha,beta,a,b) = overlap_params ( epa,esb,ni,nj,rab,ii )
        W=numpy.sqrt((sa*sb*rab*rab)**7)/480.0
        D=a[5]*(b[0]-b[2])-2.0*a[3]*(b[2]-b[4])+a[1]*(b[4]-b[6])
        E=b[5]*(a[0]-a[2])-2.0*b[3]*(a[2]-a[4])+b[1]*(a[4]-a[6])
        Sij[3-idx]=-W * numpy.sqrt(1.0/3.0) * (-D-E)
        (idx,sa,sb,alpha,beta,a,b) = overlap_params ( epa,epb,ni,nj,rab,ii )
        W=numpy.sqrt((sa*sb*rab*rab)**7)/480.0
        Sij[3]=-W*(
            a[2]*(b[6]+2.0*b[2])-a[4]*(b[0]+2.0*b[4])-b[4]*a[0]+a[6]*b[2])
        Sij[4]=0.5*W*(
            a[6]*(b[0]-b[2])+b[6]*(a[0]-a[2])+a[4]*(b[4]-b[2]-b[0])
            +b[4]*(a[4]-a[2]-a[0])+2.0*a[2]*b[2])

    return Sij

@jit(nopython=True)
def inline_Sij(ni,esa,epa,nj,esb,epb,R,Rvec):
    def cross(a, b):
        result = numpy.zeros(3)
        result[0] = a[1] * b[2] - a[2] * b[1]
        result[1] = a[2] * b[0] - a[0] * b[2]
        result[2] = a[0] * b[1] - a[1] * b[0]
        return result

    norbsi = 1 if (ni==1) else 4
    norbsj = 1 if (nj==1) else 4

    # locally aligned
    Sij = inline_Sij_ref(ni,esa,epa,R,nj,esb,epb)

    # globally aligned
    E = numpy.eye(3)
    F = numpy.eye(3)
    F[2,0],F[2,1],F[2,2] = Rvec[0],Rvec[1],Rvec[2]
    if (numpy.abs(F[2,2]) < 0.999):
        F[0,:] = cross(F[2,:], E[2,:])
    else:
        F[0,:] = cross(F[2,:], E[1,:])
    F[0,:] = F[0,:] / numpy.sqrt( numpy.dot(F[0,:],F[0,:] ) )
    F[1,:] = cross (F[2,:], F[0,:])

    S = numpy.zeros((norbsi,norbsj))
    S[0,0] = Sij[0]
    for i in range(1,norbsi):
        S[i,0] = -Rvec[i-1]*Sij[1] #- Rvec[i-1]*Sij[1]
    for j in range(1,norbsj):
        S[0,j] = Rvec[j-1]*Sij[2] #- Rvec[j-1]*Sij[1]
    for i in range(1,norbsi):
        for j in range(1,norbsj):
            scaleXX = ( F[0,i-1]*F[0,j-1] + F[1,i-1]*F[1,j-1] )
            scaleZZ = ( F[2,i-1]*F[2,j-1] )
            S[i,j] = -scaleZZ*Sij[3] + scaleXX*Sij[4]
            S[j,i] = S[i,j]
    return S

@jit(nopython=True)
def inline_S(atomno, S, ns, zetaSs, zetaPs, Rvecs, Rs, norbs, idxorbs):
    natms = norbs.shape[0]
    for i in range( natms ):
        for j in range( i+1, natms ):
            if (atomno==-1 or i==atomno or j==atomno):
                Sij = inline_Sij(
                    ns[i],zetaSs[i],zetaPs[i],ns[j],zetaSs[j],zetaPs[j],Rs[i,j],Rvecs[i,j])

                for k in range(norbs[i]):
                    for l in range(norbs[j]):
                        S[idxorbs[i]+k,idxorbs[j]+l] = Sij[k,l]
                        S[idxorbs[j]+l,idxorbs[i]+k] = Sij[k,l]
