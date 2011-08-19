// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
#ifndef typedefs_H
#define typedefs_H


typedef float real;
typedef real rvec[3];

#ifndef CPLUSPLUS 
   typedef int bool;
#endif 

#define TRUE 1
#define FALSE 0
#define STRLEN 2000
#define CONTINUE    '\\'
#define COMMENTSIGN ';'

#ifndef min
#define min(a,b) (((a) < (b)) ? (a) : (b) )
#endif
#ifndef max
#define max(a,b) (((a) > (b)) ? (a) : (b) )
#endif
#ifndef even
#define even(a) ( ( (a+1) / 2) == (a / 2) )
#endif

#endif
