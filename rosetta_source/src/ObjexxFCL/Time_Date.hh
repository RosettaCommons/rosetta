#ifndef INCLUDED_ObjexxFCL_Time_Date_hh
#define INCLUDED_ObjexxFCL_Time_Date_hh


// Time and Date Functions
//
// Project: Objexx Fortran Compatibility Library (ObjexxFCL)
//
// Version: 3.0.0
//
// Language: C++
//
// Copyright (c) 2000-2009 Objexx Engineering, Inc. All Rights Reserved.
// Use of this source code or any derivative of it is restricted by license.
// Licensing is available from Objexx Engineering, Inc.:  http://objexx.com  Objexx@objexx.com


// ObjexxFCL Headers
#include <ObjexxFCL/FArray1.fwd.hh>
#include <ObjexxFCL/Fstring.fwd.hh>


namespace ObjexxFCL {


/// @brief Current Time: HH, MM, SS
void
itime( FArray1_int & timearray );


/// @brief Current Date: DD, MM, YYYY
void
idate( FArray1_int & datearray );


/// @brief Current Date Numeric (Not Y2K Compliant): MM, DD, YY
void
idate( int & month, int & day, int & year );


/// @brief Current Date String (Not Y2K Compliant): DD-MMM-YY
void
date( Fstring & day );


} // namespace ObjexxFCL


#endif // INCLUDED_ObjexxFCL_Time_Date_HH
