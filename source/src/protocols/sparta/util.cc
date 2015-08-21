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
#include <protocols/sparta/util.hh>

#include <string>
#include <vector>
#include <sstream>
#include <iterator>
#include <utility/exit.hh>

#include <stdio.h>

#include <utility/vector0.hh>


namespace protocols {
namespace sparta {

using namespace std;
std::string simplifyWhiteSpace(const std::string &str)
{
	if ( str.empty() ) {    // nothing to do
		return str;
	}

	char to[ 20000 ];
	runtime_assert( 20000 > str.size() );
	const char *from = str.c_str();
	const char *fromend = from+str.length();
	int outc=0;
	//  char *to = &(result[0]);

	while ( true ) {
		while ( from!=fromend && isSpace(*from) ) from++;
		while ( from!=fromend && !isSpace(*from) ) to[outc++] = *(from++);
		if ( from!=fromend ) {
			to[outc++] = ' ';
		} else {
			break;
		}
	}

	if ( outc > 0 && to[outc-1] == ' ' ) {
		outc--;
	}
	to[outc]='\0';
	std::string result(to);
	return result;
}


int contains( const std::string &str, const std::string &c )
{
	int count = 0;

	int n = str.length();
	int len = c.length();

	std::string temp = str.substr(0,len);
	while ( n-- ) {    // counts overlapping stringbs
		if ( temp == c ) {
			count++;
		}
		temp = str.substr(str.length()-n,len);
	}

	return count;
}


int contains( const string &str, const char &c )
{
	int count = 0;

	int n = str.length();
	while ( n-- ) {
		if ( str[n] == c ) {
			count++;
		}
	}
	return count;
}


bool isDigit( const char &c )
{
	if ( c >= '0' && c <= '9' ) return true;
	return false;
}


bool isSpace( const char &c )
{
	if ( (c >= 9 && c <= 13) || c == ' ' ) return true;
	return false;
}


StringList split(const char sep, const string &str)
{
	string sp = " "; sp[0] = sep;
	return split( sp, str );
}


StringList split(const string &sep, const string &str)
{
	StringList lst;

	int j = 0;
	int i = str.find( sep, j );

	while ( i != -1 ) {
		if ( str.substr(j, i - j ).length() > 0 ) {
			lst.push_back( str.substr( j, i - j ) );
		}
		j = i + sep.length();
		i = str.find( sep, j );
	}

	int l = str.length() - 1;
	if ( str.substr( j, l - j + 1 ).length() > 0 ) {
		lst.push_back( str.substr( j, l - j + 1 ) );
	}

	return lst;
}


StringList split_WhiteSpace(const string &str)
{
	StringList lst;
	std::stringstream stream( str );
	copy( istream_iterator< std::string >(stream), istream_iterator< std::string >(), back_inserter( lst ) );
	//   if ( str.empty() )    // nothing to do
	//     {
	//       lst.push_back( str );
	//       return lst;
	//     }

	//   string result;
	//   result.resize(str.length()); //.setLength( length() );
	//   const char *from = str.c_str();
	//   const char *fromend = from+str.length();
	//   int outc=0;
	//   char *to = &(result[0]);

	//   while ( true ) {
	//     //while ( from!=fromend && isSpace(*from) ) //(c >= 9 && c <= 13) || c == ' '
	//     while ( from!=fromend && ((*from >= 9 && *from <= 13) || *from == ' ') )
	//       from++;
	//     while ( from!=fromend && !isSpace(*from) )
	//       //while ( from!=fromend && !((*from >= 9 && *from <= 13) || *from == ' ') )
	//       to[outc++] = *from++;

	//     if ( from!=fromend )
	//       {
	//  to[outc++] = '\0';
	//  lst.push_back( to );
	//  to = &(result[0]);
	//  outc=0;
	//       }
	//     else
	//       break;
	//   }
	//   if ( outc > 0 )
	//     {
	//       to[outc++] = '\0';
	//       lst.push_back( to );
	//     }

	//   // if ( outc > 0 && to[outc-1] == ' ' )
	//   //  outc--;
	//   //    result.resize( outc );
	//   //    return result;
	return lst;
}


//returns a section of the string, each section is defined by char 'sep', numbers of start and end are the index number (begin with 0)
char *section( const string &str, const char &sep, char *buff, int start, int end ) {
	StringList fields = split(sep, str);

	string temp = "";
	if ( start < (int) fields.size() ) {
		for ( int i = start; i <= end; i++ ) {
			if ( i >= (int) fields.size() ) break;
			temp += fields[i];
			if ( i!=end ) temp+=sep;
		}
	}

	sprintf( buff, "%s", temp.c_str() );
	return buff;
}

}
}
