// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
/*                                                                            */
/*                           ----  SPARTA   ----                              */
/*     Shifts Prediction from Analogue of Residue type and Torsion Angle      */
/*                           Yang Shen and Ad Bax                             */
/*                    J. Biomol. NMR, 38, 289-302 (2007)                      */
/*                 NIH, NIDDK, Laboratory of Chemical Physics                 */
/*                     version, 1.01 (build 2009.0928.17)                     */
/*                                                                            */
/*                      for any problem, please contact                       */
/*                          shenyang@niddk.nih.gov                            */
/*                                                                            */
/******************************************************************************/

#ifndef protocols_sparta_util_hh
#define protocols_sparta_util_hh

#include <utility/vector0.hh>

#include <vector>
#include <string>

namespace protocols {
namespace sparta {
typedef utility::vector0< std::string> StringList;
//strings functions
int contains( const std::string &str, const std::string &c );
int contains( const std::string &str, const char &c );


bool isDigit( const char &c );
bool isSpace( const char &c );

StringList split( const char sep, const std::string &str );
StringList split( const std::string &sep, const std::string &str );
//splits the std::string str into std::strings wherever a separator 'sep' occurs, and returns the list of those std::strings.
StringList split_WhiteSpace(const std::string &str);

char * section( const std::string &str, const char &seq, char *buff, int start, int end= 0xffffffff );
//returns a section of the std::string, each section is defined by char 'sep', numbers of start and end are the index number (begin with 0)

std::string simplifyWhiteSpace( const std::string &str );
//Returns a std::string that has whitespace removed from the start and the end,
//and which has each sequence of internal whitespace replaced with a single space.

}
}

#endif
