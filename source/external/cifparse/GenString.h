/*
FILE:     GenString.h
*/
/*
VERSION:  7.105
*/
/*
DATE:     10/21/2013
*/
/*
  Comments and Questions to: sw-help@rcsb.rutgers.edu
*/
/*
COPYRIGHT 1999-2013 Rutgers - The State University of New Jersey

This software is provided WITHOUT WARRANTY OF MERCHANTABILITY OR
FITNESS FOR A PARTICULAR PURPOSE OR ANY OTHER WARRANTY, EXPRESS OR
IMPLIED.  RUTGERS MAKE NO REPRESENTATION OR WARRANTY THAT THE
SOFTWARE WILL NOT INFRINGE ANY PATENT, COPYRIGHT OR OTHER
PROPRIETARY RIGHT.

The user of this software shall indemnify, hold harmless and defend
Rutgers, its governors, trustees, officers, employees, students,
agents and the authors against any and all claims, suits,
losses, liabilities, damages, costs, fees, and expenses including
reasonable attorneys' fees resulting from or arising out of the
use of this software.  This indemnification shall include, but is
not limited to, any and all claims alleging products liability.
*/
/*
               RCSB PDB SOFTWARE LICENSE AGREEMENT

BY CLICKING THE ACCEPTANCE BUTTON OR INSTALLING OR USING 
THIS "SOFTWARE, THE INDIVIDUAL OR ENTITY LICENSING THE  
SOFTWARE ("LICENSEE") IS CONSENTING TO BE BOUND BY AND IS 
BECOMING A PARTY TO THIS AGREEMENT.  IF LICENSEE DOES NOT 
AGREE TO ALL OF THE TERMS OF THIS AGREEMENT
THE LICENSEE MUST NOT INSTALL OR USE THE SOFTWARE.

1. LICENSE AGREEMENT

This is a license between you ("Licensee") and the Protein Data Bank (PDB) 
at Rutgers, The State University of New Jersey (hereafter referred to 
as "RUTGERS").   The software is owned by RUTGERS and protected by 
copyright laws, and some elements are protected by laws governing 
trademarks, trade dress and trade secrets, and may be protected by 
patent laws. 

2. LICENSE GRANT

RUTGERS grants you, and you hereby accept, non-exclusive, royalty-free 
perpetual license to install, use, modify, prepare derivative works, 
incorporate into other computer software, and distribute in binary 
and source code format, or any derivative work thereof, together with 
any associated media, printed materials, and on-line or electronic 
documentation (if any) provided by RUTGERS (collectively, the "SOFTWARE"), 
subject to the following terms and conditions: (i) any distribution 
of the SOFTWARE shall bind the receiver to the terms and conditions 
of this Agreement; (ii) any distribution of the SOFTWARE in modified 
form shall clearly state that the SOFTWARE has been modified from 
the version originally obtained from RUTGERS.  

2. COPYRIGHT; RETENTION OF RIGHTS.  

The above license grant is conditioned on the following: (i) you must 
reproduce all copyright notices and other proprietary notices on any 
copies of the SOFTWARE and you must not remove such notices; (ii) in 
the event you compile the SOFTWARE, you will include the copyright 
notice with the binary in such a manner as to allow it to be easily 
viewable; (iii) if you incorporate the SOFTWARE into other code, you 
must provide notice that the code contains the SOFTWARE and include 
a copy of the copyright notices and other proprietary notices.  All 
copies of the SOFTWARE shall be subject to the terms of this Agreement.  

3. NO MAINTENANCE OR SUPPORT; TREATMENT OF ENHANCEMENTS 

RUTGERS is under no obligation whatsoever to: (i) provide maintenance 
or support for the SOFTWARE; or (ii) to notify you of bug fixes, patches, 
or upgrades to the features, functionality or performance of the 
SOFTWARE ("Enhancements") (if any), whether developed by RUTGERS 
or third parties.  If, in its sole discretion, RUTGERS makes an 
Enhancement available to you and RUTGERS does not separately enter 
into a written license agreement with you relating to such bug fix, 
patch or upgrade, then it shall be deemed incorporated into the SOFTWARE 
and subject to this Agreement. You are under no obligation whatsoever 
to provide any Enhancements to RUTGERS or the public that you may 
develop over time; however, if you choose to provide your Enhancements 
to RUTGERS, or if you choose to otherwise publish or distribute your 
Enhancements, in source code form without contemporaneously requiring 
end users or RUTGERS to enter into a separate written license agreement 
for such Enhancements, then you hereby grant RUTGERS a non-exclusive,
royalty-free perpetual license to install, use, modify, prepare
derivative works, incorporate into the SOFTWARE or other computer
software, distribute, and sublicense your Enhancements or derivative
works thereof, in binary and source code form.

4. FEES.  There is no license fee for the SOFTWARE.  If Licensee
wishes to receive the SOFTWARE on media, there may be a small charge
for the media and for shipping and handling.  Licensee is
responsible for any and all taxes.

5. TERMINATION.  Without prejudice to any other rights, Licensor
may terminate this Agreement if Licensee breaches any of its terms
and conditions.  Upon termination, Licensee shall destroy all
copies of the SOFTWARE.

6. PROPRIETARY RIGHTS.  Title, ownership rights, and intellectual
property rights in the Product shall remain with RUTGERS.  Licensee 
acknowledges such ownership and intellectual property rights and will 
not take any action to jeopardize, limit or interfere in any manner 
with RUTGERS' ownership of or rights with respect to the SOFTWARE.  
The SOFTWARE is protected by copyright and other intellectual 
property laws and by international treaties.  Title and related 
rights in the content accessed through the SOFTWARE is the property 
of the applicable content owner and is protected by applicable law.  
The license granted under this Agreement gives Licensee no rights to such
content.

7. DISCLAIMER OF WARRANTY.  THE SOFTWARE IS PROVIDED FREE OF 
CHARGE, AND, THEREFORE, ON AN "AS IS" BASIS, WITHOUT WARRANTY OF 
ANY KIND, INCLUDING WITHOUT LIMITATION THE WARRANTIES THAT IT 
IS FREE OF DEFECTS, MERCHANTABLE, FIT FOR A PARTICULAR PURPOSE 
OR NON-INFRINGING.  THE ENTIRE RISK AS TO THE QUALITY AND 
PERFORMANCE OF THE SOFTWARE IS BORNE BY LICENSEE.  SHOULD THE 
SOFTWARE PROVE DEFECTIVE IN ANY RESPECT, THE LICENSEE AND NOT 
LICENSOR ASSUMES THE ENTIRE COST OF ANY SERVICE AND REPAIR.  
THIS DISCLAIMER OF WARRANTY CONSTITUTES AN ESSENTIAL PART OF 
THIS AGREEMENT.  NO USE OF THE PRODUCT IS AUTHORIZED HEREUNDER 
EXCEPT UNDER THIS DISCLAIMER.

8. LIMITATION OF LIABILITY.  TO THE MAXIMUM EXTENT PERMITTED BY
APPLICABLE LAW,  IN NO EVENT WILL LICENSOR BE LIABLE FOR ANY 
INDIRECT, SPECIAL, INCIDENTAL OR CONSEQUENTIAL DAMAGES ARISING 
OUT OF THE USE OF OR INABILITY TO USE THE SOFTWARE, INCLUDING, 
WITHOUT LIMITATION, DAMAGES FOR LOSS OF GOODWILL, WORK 
STOPPAGE, COMPUTER FAILURE OR MALFUNCTION, OR ANY AND ALL 
OTHER COMMERCIAL DAMAGES OR LOSSES, EVEN IF ADVISED OF THE
POSSIBILITY THEREOF. 
*/


#ifndef GENSTRING_H
#define GENSTRING_H


#include <string>
#include <functional>


/**
 ** \class Char
 **
 ** \brief Generic character class that contains character related methods.
 **
 ** This class is a static class that contains generic character related utility
 ** methods.
 */
class Char
{
  public:
    enum eCompareType
    {
        eCASE_SENSITIVE = 0,
        eCASE_INSENSITIVE,
        eWS_INSENSITIVE,  // But case-sensitive
        eAS_INTEGER
    };

    static char ToLower(const char c);
    static char ToUpper(const char c);

    static bool IsCiLess(const char c1, const char c2);

    static bool IsWhiteSpace(const char c);
    static bool IsDigit(const char c);
};


/**
 ** \class CharLess
 **
 ** \brief Public class that encapsulates character comparison.
 **
 ** This class encapsulates character comparison. It supports the following
 ** compare types: case-sensitive and case-insensitive.
 */
class CharLess
{
  public:
    CharLess(Char::eCompareType compareType = Char::eCASE_SENSITIVE);

    CharLess& operator=(const CharLess& in);

    bool operator()(const char c1, const char c2) const;

    inline Char::eCompareType GetCompareType();

  private:
    Char::eCompareType _compareType;
};


/**
 ** \class CharEqualTo
 **
 ** \brief Public class that encapsulates generic character equal_to functor.
 **
 ** This class is equal_to functor for generic character. It supports the
 ** following compare types: case-sensitive and case-insensitive.
 */
class CharEqualTo : public std::binary_function<char, char, bool>
{
  public:
    CharEqualTo(Char::eCompareType compareType = Char::eCASE_SENSITIVE);

    CharEqualTo& operator=(const CharEqualTo& in);

    bool operator()(const char c1, const char c2) const;

    inline Char::eCompareType GetCompareType();

  private:
    Char::eCompareType _compareType;
};


class WhiteSpace : public std::unary_function<char, bool>
{
  public:
    bool operator()(const char c) const;
    bool operator()(const char c1, const char c2) const;
};


/**
 ** \class StringLess
 **
 ** \brief Public class that encapsulates string comparison.
 **
 ** This class encapsulates string comparison. It supports the following
 ** compare types: case-sensitive, case-insensitive and as-integer.
 */
class StringLess
{
  public:
    StringLess(Char::eCompareType compareType = Char::eCASE_SENSITIVE);

    StringLess& operator=(const StringLess& in);

    bool operator()(const std::string& s1, const std::string& s2) const;

    inline Char::eCompareType GetCompareType();

  private:
    Char::eCompareType _compareType;
};


/**
 ** \class StringEqualTo
 **
 ** \brief Public class that encapsulates generic string equal_to functor.
 **
 ** This class is equal_to functor for generic strings. It supports the
 ** following compare types: case-sensitive, case-insensitive and as-integer.
 */
class StringEqualTo : public std::binary_function<std::string, std::string,
  bool>
{
  public:
    StringEqualTo(Char::eCompareType compareType = Char::eCASE_SENSITIVE);

    StringEqualTo& operator=(const StringEqualTo& in);

    bool operator()(const std::string& s1, const std::string& s2) const;

    inline Char::eCompareType GetCompareType();

  private:
    Char::eCompareType _compareType;
};


/**
 ** \class String
 **
 ** \brief Generic string class that contains string related utility methods.
 **
 ** This class is a static class that contains generic string related utility
 ** methods, such as: converting string to uppercase/lowercase, removing
 ** whitespaces, converting strings to/from integers/real numbers, determining
 ** if string a number, determining whether strings are equal, escaping and
 ** unescaping.
 */
class String
{
  public:
    static void LowerCase(const std::string& inString, std::string& outString);
    static void LowerCase(std::string& inOutString);
    static void UpperCase(const std::string& inString, std::string& outString);
    static void UpperCase(std::string& inOutString);

    static void RemoveWhiteSpace(const std::string& inString,
      std::string& outString);

    static std::string IntToString(int inInteger);
    static std::string DoubleToString(double inDouble);
    static int StringToInt(const std::string& inString);
    static double StringToDouble(const std::string& inString);
    static bool IsScientific(const std::string& number);
    static void ToFixedFormat(std::string& fixedFormat,
      const std::string& number);
    static bool StringToBoolean(const std::string& inString);

    static bool IsNumber(const std::string& inString);

    static bool IsCiEqual(const std::string& firstString,
      const std::string& secondString);
    static bool IsEqual(const std::string& firstString,
      const std::string& secondString,
      const Char::eCompareType compareType);

    static void StripLeadingWs(std::string& resString);
    static void StripTrailingWs(std::string& resString);
    static void StripAndCompressWs(std::string& resString);
    static void rcsb_clean_string(std::string& theString);

    static void UnEscape(std::string& outStr, const std::string& inStr);

    static void Replace(std::string& resString, const std::string& fromStr,
      const std::string& toStr);

  private:
    static std::string::const_iterator GetExpValue(int& expValue,
      const std::string::const_iterator& beg,
      const std::string::const_iterator& end);
    static void GetMantissa(std::string& mantissa, int& addExpValue,
      const std::string::const_iterator& beg,
      const std::string::const_iterator& end);
    static void ScientificNumberToFixed(std::string& fixed,
      const bool isPositive, const std::string& mantissa, const int exponent);
};


inline Char::eCompareType StringLess::GetCompareType()
{
    return (_compareType);
}

inline Char::eCompareType StringEqualTo::GetCompareType()
{
    return (_compareType);
}

#endif
