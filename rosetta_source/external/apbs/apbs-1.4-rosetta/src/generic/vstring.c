/**
 *  @file    vstring.c
 *  @ingroup Vstring
 *  @author  Nathan Baker
 *  @brief   Class Vstring methods
 *  @version
 *  @attention
 *  @verbatim
 *
 * APBS -- Adaptive Poisson-Boltzmann Solver
 *
 *  Nathan A. Baker (nathan.baker@pnnl.gov)
 *  Pacific Northwest National Laboratory
 *
 *  Additional contributing authors listed in the code documentation.
 *
 * Copyright (c) 2010-2012 Battelle Memorial Institute. Developed at the
 * Pacific Northwest National Laboratory, operated by Battelle Memorial
 * Institute, Pacific Northwest Division for the U.S. Department of Energy.
 *
 * Portions Copyright (c) 2002-2010, Washington University in St. Louis.
 * Portions Copyright (c) 2002-2010, Nathan A. Baker.
 * Portions Copyright (c) 1999-2002, The Regents of the University of
 * California.
 * Portions Copyright (c) 1995, Michael Holst.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * Neither the name of the developer nor the names of its contributors may be
 * used to endorse or promote products derived from this software without
 * specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
 * THE POSSIBILITY OF SUCH DAMAGE.
 *
 * @endverbatim
 */

#include "vstring.h"

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vstring_strcasecmp
//
//           Copyright (c) 1988-1993 The Regents of the University of
//                         California.
//           Copyright (c) 1995-1996 Sun Microsystems, Inc.
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vstring_strcasecmp(const char *s1, const char *s2) {

#if !defined(HAVE_STRCASECMP)
    unsigned char charmap[] = {
      0x00, 0x01, 0x02, 0x03, 0x04, 0x05, 0x06, 0x07,
      0x08, 0x09, 0x0a, 0x0b, 0x0c, 0x0d, 0x0e, 0x0f,
      0x10, 0x11, 0x12, 0x13, 0x14, 0x15, 0x16, 0x17,
      0x18, 0x19, 0x1a, 0x1b, 0x1c, 0x1d, 0x1e, 0x1f,
      0x20, 0x21, 0x22, 0x23, 0x24, 0x25, 0x26, 0x27,
      0x28, 0x29, 0x2a, 0x2b, 0x2c, 0x2d, 0x2e, 0x2f,
      0x30, 0x31, 0x32, 0x33, 0x34, 0x35, 0x36, 0x37,
      0x38, 0x39, 0x3a, 0x3b, 0x3c, 0x3d, 0x3e, 0x3f,
      0x40, 0x61, 0x62, 0x63, 0x64, 0x65, 0x66, 0x67,
      0x68, 0x69, 0x6a, 0x6b, 0x6c, 0x6d, 0x6e, 0x6f,
      0x70, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76, 0x77,
      0x78, 0x79, 0x7a, 0x5b, 0x5c, 0x5d, 0x5e, 0x5f,
      0x60, 0x61, 0x62, 0x63, 0x64, 0x65, 0x66, 0x67,
      0x68, 0x69, 0x6a, 0x6b, 0x6c, 0x6d, 0x6e, 0x6f,
      0x70, 0x71, 0x72, 0x73, 0x74, 0x75, 0x76, 0x77,
      0x78, 0x79, 0x7a, 0x7b, 0x7c, 0x7d, 0x7e, 0x7f,
      0x80, 0x81, 0x82, 0x83, 0x84, 0x85, 0x86, 0x87,
      0x88, 0x89, 0x8a, 0x8b, 0x8c, 0x8d, 0x8e, 0x8f,
      0x90, 0x91, 0x92, 0x93, 0x94, 0x95, 0x96, 0x97,
      0x98, 0x99, 0x9a, 0x9b, 0x9c, 0x9d, 0x9e, 0x9f,
      0xa0, 0xa1, 0xa2, 0xa3, 0xa4, 0xa5, 0xa6, 0xa7,
      0xa8, 0xa9, 0xaa, 0xab, 0xac, 0xad, 0xae, 0xaf,
      0xb0, 0xb1, 0xb2, 0xb3, 0xb4, 0xb5, 0xb6, 0xb7,
      0xb8, 0xb9, 0xba, 0xbb, 0xbc, 0xbd, 0xbe, 0xbf,
      0xc0, 0xe1, 0xe2, 0xe3, 0xe4, 0xc5, 0xe6, 0xe7,
      0xe8, 0xe9, 0xea, 0xeb, 0xec, 0xed, 0xee, 0xef,
      0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7,
      0xf8, 0xf9, 0xfa, 0xdb, 0xdc, 0xdd, 0xde, 0xdf,
      0xe0, 0xe1, 0xe2, 0xe3, 0xe4, 0xe5, 0xe6, 0xe7,
      0xe8, 0xe9, 0xea, 0xeb, 0xec, 0xed, 0xee, 0xef,
      0xf0, 0xf1, 0xf2, 0xf3, 0xf4, 0xf5, 0xf6, 0xf7,
      0xf8, 0xf9, 0xfa, 0xfb, 0xfc, 0xfd, 0xfe, 0xff,
    };

    unsigned char u1, u2;

    for ( ; ; s1++, s2++) {
    u1 = (unsigned char) *s1;
    u2 = (unsigned char) *s2;
    if ((u1 == '\0') || (charmap[u1] != charmap[u2])) {
        break;
    }
    }
    return charmap[u1] - charmap[u2];

#else

    return strcasecmp(s1, s2);

#endif

}

/* ///////////////////////////////////////////////////////////////////////////
// Routine:  Vstring_isdigit
//
//           Improves upon sscanf to see if a token is an int or not
//
//           Returns isdigit: 1 if a digit, 0 otherwise
/////////////////////////////////////////////////////////////////////////// */
VPUBLIC int Vstring_isdigit(const char *tok) {
    int i, isdigit, ti;
    char checkchar[1];
    char name[VMAX_BUFSIZE];
    strcpy(name,tok);
    isdigit = 1;
    for(i=0; ; i++){
      checkchar[0] =  name[i];
      if (name[i] == '\0'){
        break;
      }
      if (sscanf(checkchar, "%d", &ti) != 1){
        isdigit = 0;
        break;
      }
    }
    return isdigit;
}


/** Creates a wrapped and indented string from an input string
 *  @author Tucker Beck
 *  @ingroup Vstring
 *  @note    This funciton allocates a new string, so be sure to free it!
 */
char* Vstring_wrappedtext(const char* str, int right_margin, int left_padding)
{
    int span = right_margin - left_padding;
    int i = 0;
    int k = 0;
    int j = 0;
    int line_len = 0;
    int hyphenate = 0;
    char* wrap_str;
    int   wrap_len;
    int len = strlen( str );

    if( len == 0 )
        return VNULL;

    wrap_str = (char*)malloc( len * sizeof(char) );
    wrap_len = len;

    do
    {
        if( str[i] == ' ' )
        {
            i++;
        }
        else
        {
            /** @note:  The +2 is for the newline character
              *         and the null-terminating character;
              */
            if( k + right_margin + 2 > wrap_len )
            {
                wrap_len += right_margin + 2;
                wrap_str = (char*)realloc( wrap_str, wrap_len * sizeof( char ) );
            }


            if( i + span >= len )
            {
                hyphenate = 0;
                line_len = len - i;
            }
            else
            {
                j = span;
                do
                {
                    if( str[ i + j ] == ' ' )
                    {
                        hyphenate = 0;
                        line_len = j;
                        break;
                    }
                    else if( j == 0 )
                    {
                        hyphenate = 1;
                        line_len = span - 1;
                        break;
                    }
                    else
                    {
                        j--;
                    }
                } while( 1 );
            }

            memset( wrap_str + k, ' ', left_padding * sizeof( char ) );
            k += left_padding;

            memcpy( wrap_str + k, str + i, line_len * sizeof( char ) );
            k += line_len;
            i += line_len;

            if( hyphenate )
                wrap_str[k++] = '-';

            wrap_str[k++] = '\n';

            wrap_str[k] = '\0';
        }
    } while( i < len );

    return wrap_str;
}

