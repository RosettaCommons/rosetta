/*
FILE:     mapped_vector.C
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


#ifndef MAPPED_VECTOR_C
#define MAPPED_VECTOR_C


#include <stdexcept>
#include <vector>

#include "mapped_vector.h"


using std::out_of_range;
using std::vector;


template <typename T, typename StringCompareT>
mapped_vector<T, StringCompareT>::mapped_vector()
{

    _current.first.clear();
    _current.second = 0;

}


template <typename T, typename StringCompareT>
mapped_vector<T, StringCompareT>::mapped_vector(const StringCompareT& cmp)
  : _index(cmp)
{

    _current.first.clear();
    _current.second = 0;

}


template <typename T, typename StringCompareT>
mapped_vector<T, StringCompareT>::mapped_vector(
  const mapped_vector& inMappedVector)
{

    _index = inMappedVector._index;
    _vector = inMappedVector._vector;

}


template <typename T, typename StringCompareT>
mapped_vector<T, StringCompareT>::~mapped_vector()
{

    clear();

}


template <typename T, typename StringCompareT>
void mapped_vector<T, StringCompareT>::push_back(const T& inT)
{

    _vector.push_back(inT);

    typename tIndex::value_type valuePair(inT, _vector.size() - 1);

    _index.insert(valuePair);

    _current.first = inT;
    _current.second = _vector.size() - 1;

}


template <typename T, typename StringCompareT>
unsigned int mapped_vector<T, StringCompareT>::size() const
{

    return(_vector.size());

}


template <typename T, typename StringCompareT>
bool mapped_vector<T, StringCompareT>::empty() const
{

    return(_vector.empty());

}


template <typename T, typename StringCompareT>
void mapped_vector<T, StringCompareT>::operator=(
  const mapped_vector& inMappedVector)
{

    _index = inMappedVector._index;
    _vector = inMappedVector._vector;

}


template <typename T, typename StringCompareT>
void mapped_vector<T, StringCompareT>::operator=(const vector<T>& inVector)
{

    clear();

    for (unsigned int index = 0; index < inVector.size(); ++index)
    {
        push_back(inVector[index]);
    }

}


template <typename T, typename StringCompareT>
bool mapped_vector<T, StringCompareT>::operator==(
  const mapped_vector& inMappedVector)
{

    return(_vector == inMappedVector._vector);

}


template <typename T, typename StringCompareT>
bool mapped_vector<T, StringCompareT>::operator!=(
  const mapped_vector& inMappedVector)
{

    return(!operator==(inMappedVector));

}


template <typename T, typename StringCompareT>
const T& mapped_vector<T, StringCompareT>::operator[](unsigned int index) const
{

    if (index >= size())
    {
        throw out_of_range("Invalid index in mapped_vector::operator[]");
    }

    return(_vector[index]);

}


template <typename T, typename StringCompareT>
const vector<T>& mapped_vector<T, StringCompareT>::get_vector() const
{

    return(_vector);

}


template <typename T, typename StringCompareT>
vector<T>& mapped_vector<T, StringCompareT>::get_vector()
{

    return(_vector);

}


template <typename T, typename StringCompareT>
void mapped_vector<T, StringCompareT>::erase(const T& inT)
{

    unsigned int index = get_index(inT);
    if (index >= size())
    {
        throw out_of_range("Element not found in mapped_vector::erase");
    }

    if (is_equal(_current.first, _vector[index]))
    {
        _current.first.clear();
        _current.second = 0;
    }

    _vector.erase(_vector.begin() + index);

    _index.erase(inT);

    for (typename tIndex::iterator pos = _index.begin();
      pos != _index.end(); ++pos)
    {
        if (pos->second >= index)
        {
            --(pos->second);
        }
    }

}


template <typename T, typename StringCompareT>
void mapped_vector<T, StringCompareT>::insert(const unsigned int index,
  const T& inT)
{

    unsigned int existingIndex = get_index(inT);
    if (existingIndex != size())
    {
        throw out_of_range("Element exists in mapped_vector::insert");
    }

    _current.first = inT;
    _current.second = index;

    _vector.insert(_vector.begin() + index, inT);

    for (typename tIndex::iterator pos = _index.begin(); pos != _index.end();
      ++pos)
    {
        if (pos->second >= index)
        {
            ++(pos->second);
        }
    }

    typename tIndex::value_type valuePair(inT, index);

    _index.insert(valuePair);

}


template <typename T, typename StringCompareT>
void mapped_vector<T, StringCompareT>::index_it()
{

    for (unsigned int index = 0; index < _vector.size(); ++index)
    {
        typename tIndex::value_type valuePair(_vector[index], index);

        _index.insert(valuePair);
    }

}


template <typename T, typename StringCompareT>
void mapped_vector<T, StringCompareT>::clear()
{

    _index.clear();
    _vector.clear();

    _current.first.clear();
    _current.second = 0;

}


template <typename T, typename StringCompareT>
unsigned int mapped_vector<T, StringCompareT>::find(const T& inT) const
{

    return(get_index(inT));

}


template <typename T, typename StringCompareT>
unsigned int mapped_vector<T, StringCompareT>::get_index(const T& inT) const
{

    if (is_equal(_current.first, inT))
    {
        return(_current.second);
    }

    // Return index of found value or invalid index
    typename tIndex::const_iterator pos = _index.find(inT);
    if (pos != _index.end())
    {
        // Found
        _current.first = inT;
        _current.second = pos->second;

        return(pos->second);
    }
    else
    {
        // Not found. Return invalid index.
        return(_vector.size());
    }

}


template <typename T, typename StringCompareT>
bool mapped_vector<T, StringCompareT>::is_equal(const T& firstT,
  const T& secondT) const
{

    typename tIndex::key_compare keyComp = _index.key_comp();

    return(!(keyComp(firstT, secondT) || keyComp(secondT, firstT)));

}


#endif

