/* Author		Patrice Koehl
 * Revision #		1
 * Date			6/10/2005
 * 
 * This file contains all changeable constants for Alphavol. In particular,
 * Alphavol uses static memory allocation (for speed issue), and the size
 * of all arrays it needs are given in this file.
 */
#ifndef __DEFINES__
#define __DEFINES__

/**********************************************************************
 * general dimensions                                               
 *
 **********************************************************************/

#ifndef MAX_LENGTH
#define MAX_LENGTH 300
#endif

/***********************************************************************
 * Constants for Conversion to arbitrary precision arithmetics (APA)
 * Alphavol computes all geometric predicates first in floating points,
 * and if the results is imprecise, it switches to APA, using the GMP package.
 * Any value below a cutoff EPS triggers the switch to APA.
 * When a float is converted to APA, it is first multiplied by SCALE,
 * such that all its significant digits (NDIGIT) are taken care of.
 * 
 * EPS and SCALE are set up here.
 * The value of EPS should be quite safe. If speed is an issue for you, you 
 * you might want to decrease it (by 10, or even 100).
 * SCALE is set up to take into account 8 digits; this can also be changed
 *
 ***********************************************************************/

#ifndef EPS_APA
#define EPS_APA 0.0001
#endif

#ifndef SCALE_APA
#define SCALE_APA 100000000.
#endif

#ifndef NDIGIT_APA
#define NDIGIT_APA 8
#endif

/***********************************************************************
 * ARRAY sizes: these are constants you might want to change to suit 
 * your preferences.
 * The original AlphaBall is set for up to 1000000 atoms
 *
 ***********************************************************************/

#ifndef FSIZE
#define FSIZE 100
#endif
/*40290*/

#ifndef MAX_ATOM
#define MAX_ATOM 1000000
#endif

#ifndef MAX_COORD
#define MAX_COORD (3 * MAX_ATOM)
#endif

#ifndef MAX_POINT
#define MAX_POINT MAX_ATOM
#endif

/* #ifndef MAX_PAIR
 * #define MAX_PAIR (MAX_ATOM * (MAX_ATOM-1))/2
 * #endif
 */

/* 
 * The number of tetrahedra is arbitrarily set to 10 times the number of atoms.
 * This works well for molecules, but it may fail for other objects
 */

#ifndef MAX_TETRA
#define MAX_TETRA (10 * MAX_POINT)
#endif


/***********************************************************************
 * Size of auxiliary arrays: you should not have to change these values
 **********************************************************************/

#ifndef MAX_NEW
#define MAX_NEW 10000
#endif

#ifndef MAX_RED
#define MAX_RED 10000
#endif

#ifndef MAX_LINK
#define MAX_LINK 10000
#endif

#ifndef MAX_FACET
#define MAX_FACET 10000
#endif

#ifndef MAX_FREE
#define MAX_FREE 10000
#endif

/**********************************************************************
 * Define pi and multiples
 **********************************************************************/

#ifndef PI
#define PI 3.141592653589793238462643
#endif

#ifndef TWOPI
#define TWOPI (2 * PI)
#endif

#endif /* __DEFINES__ */

