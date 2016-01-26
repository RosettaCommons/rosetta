/*
FILE:     GenString.C
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


#include <stdexcept>
#include <functional>
#include <algorithm>
#include <iterator>
#include <string>
#include <sstream>

#include "Exceptions.h"
#include "GenString.h"


using std::exception;
using std::out_of_range;
using std::runtime_error;
using std::not_equal_to;
using std::equal_to;
using std::not1;
using std::lexicographical_compare;
using std::find_if;
using std::remove_copy_if;
using std::replace_if;
using std::unique;
using std::transform;
using std::back_inserter;
using std::bind2nd;
using std::string;
using std::ios;
using std::istringstream;
using std::ostringstream;


char Char::ToLower(const char c)
{
    if (!((c >= 'A') && (c <= 'Z')))
    {
        // If not an uppercase letter, just return that character.
        return c;
    }

    return (c + ('a' - 'A'));
}


char Char::ToUpper(const char c)
{
    if (!((c >= 'a') && (c <= 'z')))
    {
        // If not a lowercase letter, just return that character.
        return c;
    }

    return (c - ('a' - 'A'));
}


bool Char::IsCiLess(const char c1, const char c2)
{
    return (Char::ToLower(c1) < Char::ToLower(c2));
}


bool Char::IsWhiteSpace(const char c)
{
    return ((c == ' ') || (c == '\n') || (c == '\t') || (c == '\f') ||
      (c == '\v') || (c == '\r'));
}


bool Char::IsDigit(const char c)
{
    return ((c >= '0') && (c <= '9'));
}


CharLess::CharLess(Char::eCompareType compareType) :
  _compareType(compareType)
{

}


CharLess& CharLess::operator=(const CharLess& in)
{
    _compareType = in._compareType;

    return (*this);
}


bool CharLess::operator()(const char c1, const char c2) const
{
    switch (_compareType)
    {
        case Char::eCASE_SENSITIVE:
        {
            return (c1 < c2);
            break;
        }
        case Char::eCASE_INSENSITIVE:
        {
            return (Char::IsCiLess(c1, c2));
            break;
        }
        case Char::eWS_INSENSITIVE:
        {
            if (Char::IsWhiteSpace(c1) && Char::IsWhiteSpace(c2))
            {
                return (false);
            }
            else
            {
                return (c1 < c2);
            }
            break;
        }
        default:
        {
            throw out_of_range("Invalid compare type in CharLess::operator()");
            break;
        }
    }

    return (true);
}


CharEqualTo::CharEqualTo(Char::eCompareType compareType) :
  _compareType(compareType)
{

}


CharEqualTo& CharEqualTo::operator=(const CharEqualTo& in)
{
    _compareType = in._compareType;

    return (*this);
}


bool CharEqualTo::operator()(const char c1, const char c2) const
{
    CharLess ciLess(_compareType);

    return (!ciLess(c1, c2) && !ciLess(c2, c1));
}


bool WhiteSpace::operator()(const char c) const
{
    return (Char::IsWhiteSpace(c));
}


bool WhiteSpace::operator()(const char c1, const char c2) const
{
    return ((Char::IsWhiteSpace(c1)) && (Char::IsWhiteSpace(c2)));
}


StringLess::StringLess(Char::eCompareType compareType) :
  _compareType(compareType)
{

}


StringLess& StringLess::operator=(const StringLess& in)
{
    _compareType = in._compareType;

    return (*this);
}


bool StringLess::operator()(const string& s1, const string& s2) const
{
    switch (_compareType)
    {
        case Char::eCASE_SENSITIVE:
        {
            return (s1 < s2);
            break;
        }
        case Char::eCASE_INSENSITIVE:
        {
            return (lexicographical_compare(s1.begin(), s1.end(),
              s2.begin(), s2.end(), Char::IsCiLess));
            break;
        }
        case Char::eAS_INTEGER:
        {
            int s1AsInt = String::StringToInt(s1);
            int s2AsInt = String::StringToInt(s2);
            return (s1AsInt < s2AsInt);
            break;
        }
        default:
        {
            throw out_of_range("Invalid compare type in "\
              "StringLess::operator()");
            break;
        }
    }

    return (true);
}


StringEqualTo::StringEqualTo(Char::eCompareType compareType) :
  _compareType(compareType)
{

}


StringEqualTo& StringEqualTo::operator=(const StringEqualTo& in)
{
    _compareType = in._compareType;

    return (*this);
}


bool StringEqualTo::operator()(const string& s1, const string& s2) const
{
    StringLess ciLess(_compareType);

    return (!ciLess(s1, s2) && !ciLess(s2, s1));
}


void String::LowerCase(const string& inString, string& outString)
{
    outString.clear();

    transform(inString.begin(), inString.end(), back_inserter(outString),
      Char::ToLower);
}


void String::LowerCase(string& resString)
{
    transform(resString.begin(), resString.end(), resString.begin(),
      Char::ToLower);
}


void String::UpperCase(const string& inString, string& outString)
{
    outString.clear();

    transform(inString.begin(), inString.end(), back_inserter(outString),
      Char::ToUpper);
}


void String::UpperCase(string& resString)
{
    transform(resString.begin(), resString.end(), resString.begin(),
      Char::ToUpper);
}


void String::RemoveWhiteSpace(const string& inString, string& outString)
{
    outString = inString;

    outString.erase(std::remove_if(outString.begin(), outString.end(),
      Char::IsWhiteSpace), outString.end());
}


string String::IntToString(int inInteger)
{
    ostringstream outStringStream;

    outStringStream << inInteger;

    return (outStringStream.str());
}


string String::DoubleToString(double inDouble)
{
    ostringstream outStringStream;

    outStringStream << inDouble;

    return (outStringStream.str());
}


int String::StringToInt(const string& inString)
{
    // No need to check if empty string. It will be picked up by the code
    // below.

    try
    {
        istringstream inStringStream(inString);

        inStringStream.exceptions(ios::badbit | ios::failbit);

        int ret;

        inStringStream >> ret;

        if (!inStringStream.eof())
        {
            throw runtime_error("Could not convert \"" + inString +
              "\" to a number in String::StringToInt");
        }

        return (ret);
    }
    catch (exception& e)
    {
        throw runtime_error("Could not convert \"" + inString +
          "\" to a number in String::StringToInt");
    }
}


double String::StringToDouble(const string& inString)
{
    // No need to check if empty string. It will be picked up by the code
    // below.

    try
    {
        istringstream inStringStream(inString);

        inStringStream.exceptions(ios::badbit | ios::failbit);

        double ret;

        inStringStream >> ret;

        if (!inStringStream.eof())
        {
            throw runtime_error("Could not convert \"" + inString +
              "\" to a number in String::StringToDouble");
        }

        return (ret);
    }
    catch (exception& e)
    {
        throw runtime_error("Could not convert \"" + inString +
          "\" to a number in String::StringToDouble");
    }
}


void String::ToFixedFormat(string& fixedFormat, const string& number)
{
    // Verify if the input is a valid number
    String::StringToDouble(number);

    bool isPositive = true;
    string::size_type numStrartInd = 0;

    if (number[0] == '-')
    {
        isPositive = false;
        numStrartInd = 1;
    }

    int expValue;

    string::const_iterator expIter = GetExpValue(expValue,
      number.begin() + numStrartInd, number.end());

    string mantissa;
    int addExpValue;
    GetMantissa(mantissa, addExpValue, number.begin() + numStrartInd,
      expIter);

    expValue += addExpValue;

    String::ScientificNumberToFixed(fixedFormat, isPositive, mantissa,
      expValue);
}


bool String::IsScientific(const string& number)
{
    // Find exponent letter
    string::const_iterator expIter = find_if(number.begin(), number.end(),
      bind2nd(CharEqualTo(Char::eCASE_INSENSITIVE), 'E'));

    if (expIter == number.end())
        return (false);
    else
        return (true);
}


string::const_iterator String::GetExpValue(int& expValue,
  const string::const_iterator& beg, const string::const_iterator& end)
{
    // Find exponent letter
    string::const_iterator expIter = find_if(beg, end,
      bind2nd(CharEqualTo(Char::eCASE_INSENSITIVE), 'E'));

    if (expIter == end)
    {
        expValue = 0;
    }
    else
    {
        expValue = StringToInt(string(expIter + 1, end));
    }

    return (expIter);
}


void String::GetMantissa(string& mantissa, int& addExpValue,
  const string::const_iterator& beg, const string::const_iterator& end)
{
    mantissa.clear();

    // Find the first non-zero digit.
    string::const_iterator firstNonZeroIter = find_if(beg, end,
      bind2nd(not_equal_to<char>(), '0'));

    if (*firstNonZeroIter == '.')
    {
        firstNonZeroIter = find_if(firstNonZeroIter + 1, end,
          bind2nd(not_equal_to<char>(), '0'));
    }

    if (firstNonZeroIter == end)
    {
        // Check if the mantissa is 0
        if (!(StringToDouble(string(beg, end)) == 0))
        {
            throw EmptyValueException("No mantissa", "String::GetMantissa");
        }
        else
        {
            mantissa.push_back('0');
            mantissa.push_back('.');
            addExpValue = 0;

            return;
        }
    }

    // Find the period
    string::const_iterator dotIter = find_if(beg, end,
      bind2nd(equal_to<char>(), '.'));

    mantissa.push_back(*firstNonZeroIter);
    mantissa.push_back('.');
    remove_copy_if(firstNonZeroIter + 1, end, back_inserter(mantissa),
      bind2nd(equal_to<char>(), '.'));

    // Strip trailing zeros from mantissa (from the end to the period)
    string::reverse_iterator lastNonZeroRevIter = find_if(mantissa.rbegin(),
      mantissa.rend(), bind2nd(not_equal_to<char>(), '0'));

    string::iterator lastNonZeroNormIter(lastNonZeroRevIter.base());

    mantissa.erase(lastNonZeroNormIter, mantissa.end());

    // IMPLEMENT - Simplify using ? operator 
    if (firstNonZeroIter < dotIter)
        addExpValue = dotIter - firstNonZeroIter - 1;
    else
        addExpValue = dotIter - firstNonZeroIter;
}
 

void String::ScientificNumberToFixed(string& fixed, const bool isPositive,
  const string& mantissa, const int exponent)
{
    // Check validity of numbers
    String::StringToDouble(mantissa);

    fixed.clear();

    if (!isPositive)
        fixed.push_back('-');

    // Build a fixed float format
    if (exponent < 0)
    {
        fixed.push_back('0');
        fixed.push_back('.');
        fixed.append(abs(exponent) - 1, '0');
        remove_copy_if(mantissa.begin(), mantissa.end(),
          back_inserter(fixed), bind2nd(equal_to<char>(), '.'));
    }
    else
    {
        fixed.push_back(mantissa[0]);
        int decPlaces = mantissa.size() - 2;
        if (decPlaces > exponent)
        {
            fixed.append(mantissa.begin() + 2, mantissa.begin() + 2 + exponent);
            fixed.push_back('.');
            fixed.append(mantissa.begin() + 2 + exponent, mantissa.end());
        }
        else
        {
            fixed.append(mantissa.begin() + 2, mantissa.begin() + 2 +
              decPlaces);
            fixed.append(abs(exponent) - decPlaces, '0');
            fixed.push_back('.');
        }
    }
}


bool String::StringToBoolean(const string& inString)
{
    return (StringToInt(inString) != 0);
}


bool String::IsNumber(const string& inString)
{
    try
    {
        String::StringToDouble(inString);

        return (true);
    }
    catch (exception)
    {
        return (false);
    }
}


bool String::IsCiEqual(const string& firstString, const string& secondString)
{
    return (IsEqual(firstString, secondString, Char::eCASE_INSENSITIVE));
}


bool String::IsEqual(const string& firstString, const string& secondString,
  const Char::eCompareType compareType)
{
    StringEqualTo stringEqualTo(compareType);

    return (stringEqualTo(firstString, secondString));
}


void String::StripLeadingWs(string& resString)
{
    string::iterator nonWhiteIter = find_if(resString.begin(),
      resString.end(), not1(WhiteSpace()));

    resString.erase(resString.begin(), nonWhiteIter);
}


void String::StripTrailingWs(string& resString)
{
    string::reverse_iterator nonWhiteRevIter = find_if(resString.rbegin(),
      resString.rend(), not1(WhiteSpace()));

    string::iterator nonWhiteIter = nonWhiteRevIter.base();

    resString.erase(nonWhiteIter, resString.end());
}


void String::StripAndCompressWs(string& resString)
{
    StripLeadingWs(resString);

    StripTrailingWs(resString);

    // unique() algorithm moves duplicate elements at the end and returns
    // logical end of the unique sequence.
    string::iterator endUnique = unique(resString.begin(), resString.end(),
      WhiteSpace());

    // Remove duplicate elements at the end
    resString.erase(endUnique, resString.end());

    replace_if(resString.begin(), resString.end(), WhiteSpace(), ' ');
}


void String::rcsb_clean_string(string& resString)
{
    StripLeadingWs(resString);

    StripTrailingWs(resString);

    replace_if(resString.begin(), resString.end(), WhiteSpace(), ' ');
}


void String::UnEscape(string& outStr, const string& inStr)
{
    outStr = inStr;

    Replace(outStr, "\\n", "\n");
    Replace(outStr, "\\t", "\t");
    Replace(outStr, "\\r", "\r");
}


void String::Replace(string& resString, const string& fromStr,
  const string& toStr)
{
    if (resString.empty() || fromStr.empty() || toStr.empty())
    {
        return;
    }

    string::size_type currInd = 0;
    string::size_type fromStrInd = 0;

    while ((fromStrInd = resString.find(fromStr, currInd)) != string::npos)
    {
        resString.replace(fromStrInd, fromStr.size(), toStr);
        currInd = fromStrInd + toStr.size();
    }
}

