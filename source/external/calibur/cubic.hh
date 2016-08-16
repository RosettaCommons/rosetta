// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/kalngyk/cubic.hh
/// @author YK Ng & SC Li (kalngyk@gmail.com)

#ifndef apps_pilot_kalngyk_cubic_HH
#define apps_pilot_kalngyk_cubic_HH

#include <time.h>
#ifndef __WIN32__
#include <sys/resource.h>
#endif
#include <stdlib.h>
#include <iostream>
#include <cmath>


using namespace std;

#define cubic_roots(a,b,c,z) cubic_roots2(a,b,c,z)

#define PI 3.14159265358979323846

//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -

#ifdef __TEST_CUBIC__
#ifndef __WIN32__
static double
__timeval_difference(struct timeval * x, struct timeval * y)
{
    double elapsed;
    if (x->tv_usec < y->tv_usec)
    {
        int nsec = (y->tv_usec - x->tv_usec) / 1000000 + 1;
        y->tv_usec -= 1000000 * nsec;
        y->tv_sec += nsec;
    }
    if (x->tv_usec - y->tv_usec > 1000000)
    {
        int nsec = (x->tv_usec - y->tv_usec) / 1000000;
        y->tv_usec += 1000000 * nsec;
        y->tv_sec -= nsec;
    }

    elapsed = x->tv_sec - y->tv_sec;
    elapsed += (x->tv_usec - y->tv_usec)/(double)1000000;
    return elapsed;
}

static double
_get_elapsed(int set_start)
{
    static struct rusage last;
    struct rusage now;
    double elapsed = 0;
    if (set_start)
        getrusage(RUSAGE_SELF, &last);
    else
    {
        getrusage(RUSAGE_SELF, &now);
        elapsed += __timeval_difference(&(now.ru_utime), &(last.ru_utime));
        elapsed += __timeval_difference(&(now.ru_stime), &(last.ru_stime));
        last = now;
    }
    return elapsed;
}
#endif
#endif

//- = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = - = -


void cubic_roots1(double a2, double a1, double a0, double * z)
{
    double p_over_3 = (3*a1 - a2*a2)/9.;
    double q_over_2 = (9*a2*a1 - 27*a0 - 2*a2*a2*a2) /54.;
    if (p_over_3 >= 0)
    {
        // can proceed to obtain root using sinh(C)
        z[0] = z[1] = z[2] = 0;
        return;
    }
    else // p is negative
    {
        p_over_3 = -p_over_3;
        double C = q_over_2 / sqrt(p_over_3 * p_over_3 * p_over_3);
        if (C >= -1 && C <= 1)
        {
            // Then, cos(3*theta) = C. This gives 3 solutions, respectively
            // cos(acos(C)/3), cos((acos(C)+2pi)/3), cos((acos(C)+4pi)/3)
            double theta = acos(C);
            double a2_over_3 = a2/3.;
            z[0] = z[1] = z[2] = 2 * sqrt(p_over_3);
            z[0] = z[0] * cos( theta          / 3.) - a2_over_3;
            z[1] = z[1] * cos((theta + 2.*PI) / 3.) - a2_over_3;
            z[2] = z[2] * cos((theta + 4.*PI) / 3.) - a2_over_3;
        }
        else // only one root is real
        {
            // can proceed to obtain root using cosh(C)
            z[0] = z[1] = z[2] = 0;
        }
    }
}

/**
 * This is faster but leaves no option for obtaining imaginary roots.
 * (After "-O" in GCC, the speed gain is not much, only ~4%)
 */
void cubic_roots2(double a2, double a1, double a0, double* z)
{
    double Q = (3*a1 - a2*a2)/9.;
    double R = (9*a2*a1 - 27*a0 - 2*a2*a2*a2) /54.;
    double Q3 = Q*Q*Q;
    if (Q3 + R*R < 0) // condition for 3 real roots
    {
        // Q3+R2<0 implies Q3<0, and sqrt(-Q3) is real
        double negsqrtQ3 = sqrt(-Q3);
        double theta = acos(R/negsqrtQ3);
        double a2_over_3 = a2/3.;
        z[0] = z[1] = z[2] = 2 * sqrt(-Q);
        z[0] = z[0] * cos( theta          /3.) - a2_over_3;
        z[1] = z[1] * cos((theta + 2.*PI) /3.) - a2_over_3;
        z[2] = z[2] * cos((theta + 4.*PI) /3.) - a2_over_3;
    }
    else // only one root is real
    {
        z[0] = z[1] = z[2] = 0; // give up...
    }
}


#ifdef __TEST_CUBIC__
int main(void)
{
    double z[3];
    double x, time;
    _get_elapsed(1);
    for (int i=-300; i < 300; i++)
        for (int j=-300; j < 300; j++)
            for (int k=-300; k < 300; k++)
            {
                cubic_roots1((double)i, (double)j, (double)k, z);
                //cout << z[0] << "; " << z[1] << "; " << z[2] << endl;
                //get_roots_of_cubic2((double)i, (double)j, (double)k, z);
                //cout << z[0] << "; " << z[1] << "; " << z[2] << endl;
            }
    time = _get_elapsed(0);
    cout << "Used " << time << " seconds." << endl;

    _get_elapsed(1);
    for (int i=-300; i < 300; i++)
        for (int j=-300; j < 300; j++)
            for (int k=-300; k < 300; k++)
            {
                cubic_roots2((double)i, (double)j, (double)k, z);
            }
    time = _get_elapsed(0);
    cout << "Used " << time << " seconds." << endl;
}
#endif

#endif
