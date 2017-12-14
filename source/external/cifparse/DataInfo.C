/*
FILE:     DataInfo.C
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


#include <iostream>

#include "Exceptions.h"
#include "rcsb_types.h"
#include "GenCont.h"
#include "RcsbFile.h"
#include "CifString.h"
#include "DataInfo.h"


using std::string;
using std::vector;
using std::cerr;
using std::endl;


#ifndef VLAD_ATOM_SITES_ALT_ID_IGNORE
string CIF_ITEM;
#endif


DataInfo::DataInfo()
= default;


DataInfo::~DataInfo()
= default;


bool DataInfo::IsUnknownValueAllowed(const string& catName, 
  const string& attribName)
{
    return (true);
}


bool DataInfo::AreItemsValuesValid(const string& catName,
  const vector<string>& columnNames, const vector<unsigned int>& columnIndices,
  const vector<bool>& allowedNullAttribs, const vector<string>& row,
  const Char::eCompareType compareType)
{
    if (row.empty())
    {
        cerr << "  Empty row detected." << endl;
        return (false);
    }

    for (unsigned int i = 0; i < columnNames.size(); ++i)
    {
        if (IsKeyItem(catName, columnNames[i], compareType))
        {
            if (CifString::IsEmptyValue(row[columnIndices[i]]))
            {
                string cifItem;
                CifString::MakeCifItem(cifItem, catName, columnNames[i]);

#ifndef VLAD_ATOM_SITES_ALT_ID_IGNORE
                if ((cifItem == "_atom_sites_alt.id") &&
                  (row[columnIndices[i]] == CifString::InapplicableValue))
                {
                    CIF_ITEM = cifItem;
                    return (false);
                }
#endif
                cerr << "  Key item \"" << cifItem << "\" has invalid value \""
                  << row[columnIndices[i]] << "\"" << endl;

                return (false);
            }
        }
        else
        {
            if (CifString::IsUnknownValue(row[columnIndices[i]]))
            {
                if (!allowedNullAttribs[i])
                {
                    string cifItem;
                    CifString::MakeCifItem(cifItem, catName, columnNames[i]);

                    cerr << "  Non-key item \"" << cifItem << "\", that must "\
                      "not have unknown value, has invalid value \"" <<
                      row[columnIndices[i]] << "\"" << endl;

                    return (false);
                }
            }
        }
    }

    return (true);
}


bool DataInfo::IsKeyItem(const string& catName,
  const string& attribName, const Char::eCompareType compareType)
{
    if (attribName.empty())
        return(false);

    const vector<string>& keys = GetCatKeys(catName);

    string itemName;
    CifString::MakeCifItem(itemName, catName, attribName);

    return (GenCont::IsInVector(itemName, keys));
}


bool DataInfo::AreAllKeyItems(const string& catName,
  const vector<string>& attributes)
{
    // Return true for cases where all attributes are key.

    if (attributes.empty())
        return(false);

    const vector<string>& keys = GetCatKeys(catName);
 
    unsigned int keyCount = 0;

    for (const auto & attribute : attributes)
    {
        if (IsKeyItem(catName, attribute))
        {
            keyCount++;
            continue;
        }
        else
        {
            keyCount = 0;
            break;
        }
    }

    if (!keys.empty() && (keyCount == keys.size()))
    {
        return (true);
    }
    else
    {
        if (keyCount != 0)
        {
            cerr << "Missing keys in category: " << catName << endl;
        }
        return (false);
    }
}


bool DataInfo::MustConvertItem(const string& catName,
  const string& itemName)
{
    string cifItem;
    CifString::MakeCifItem(cifItem, catName, itemName);

    return (IsItemDefined(cifItem));
}


void DataInfo::GetItemsTypes(vector<eTypeCode>& itemsTypes,
  const string& catName, const vector<string>& attribsNames)
{
    itemsTypes.clear();

    for (const auto & attribsName : attribsNames)
    {
        string cifItem;
        CifString::MakeCifItem(cifItem, catName, attribsName);

        if (!MustConvertItem(catName, attribsName))
        {
            itemsTypes.push_back(eTYPE_CODE_NONE);

            continue;
        }

        itemsTypes.push_back(_GetDataType(cifItem));
    }
}


eTypeCode DataInfo::_GetDataType(const string& itemName) 
{
    const vector<string>& primitiveCode = GetItemAttribute(itemName,
      CifString::CIF_DDL_CATEGORY_ITEM_TYPE_LIST,
      CifString::CIF_DDL_ITEM_PRIMITIVE_CODE);

    const vector<string>& dataType = GetItemAttribute(itemName,
      CifString::CIF_DDL_CATEGORY_ITEM_TYPE,
      CifString::CIF_DDL_ITEM_CODE);

    if (primitiveCode.empty() || dataType.empty())
    {
        cerr << "Warning - ITEM_TYPE_CODE: Item \"" + itemName +
          "\" has no type code defined. Using default string type." << endl;
#ifdef VLAD_TEST
        throw EmptyValueException("Item \"" + itemName + "\" has no type "\
          "code defined.", "DataInfo::_GetDataType");
#endif

        // Missing item primitive or type code, assume string type
        return (eTYPE_CODE_STRING);
    }

    if (primitiveCode[0].empty() || dataType[0].empty())
    {
        cerr << "Warning - ITEM_TYPE_CODE: Item \"" + itemName +
          "\" has no type code defined. Using default string type." << endl;
#ifdef VLAD_TEST
        throw EmptyValueException("Item \"" + itemName + "\" has no type "\
          "code defined.", "DataInfo::_GetDataType");
#endif
        // Missing item primitive or type code, assume string type
        return (eTYPE_CODE_STRING);
    }

    if (String::IsCiEqual(primitiveCode[0], "char") ||
      String::IsCiEqual(primitiveCode[0], "uchar"))
    {
        if (String::IsCiEqual(dataType[0], "yyyy-mm-dd"))
        {
            return (eTYPE_CODE_DATETIME);
        }
        else if (String::IsCiEqual(dataType[0], "text"))
        {
            return (eTYPE_CODE_TEXT);
        }
        else
        {
            // Assume string type
            return (eTYPE_CODE_STRING);
        }
    }
    else if (String::IsCiEqual(primitiveCode[0], "numb"))
    {
        if (String::IsCiEqual(dataType[0], "float"))
        {
            return (eTYPE_CODE_FLOAT);
        }
        else if (String::IsCiEqual(dataType[0], "double"))
        {
            return (eTYPE_CODE_FLOAT);
        }
        else if (String::IsCiEqual(dataType[0], "int"))
        {
            return (eTYPE_CODE_INT);
        }
        else if (String::IsCiEqual(dataType[0], "numb"))
        {
            return (eTYPE_CODE_INT);
        }
        else
        {
            // Assume string type
            return (eTYPE_CODE_STRING);
        }
    }
    else
    {
        throw EmptyValueException("Item \"" + itemName + "\" has invalid "\
          "primitive code \"" + primitiveCode[0] + "\"",
          "DataInfo::_GetDataType");
    }
}


bool DataInfo::IsItemMandatory(const string& catName,
  const string& attribName)
{
    string itemName;
    CifString::MakeCifItem(itemName, catName, attribName);

    return (IsItemMandatory(itemName));
}


void DataInfo::GetMandatoryItems(vector<string>& mandItemsNames,
  const string& catName)
{
    mandItemsNames.clear();

    const vector<string>& itemsNames = GetItemsNames();
    for (const auto & itemName : itemsNames)
    {
        string currCatName;
        CifString::GetCategoryFromCifItem(currCatName, itemName);
        if (currCatName != catName)
        {
            continue;
        }

        if (IsItemMandatory(itemName))
        {
            mandItemsNames.push_back(itemName);
        }
    }
}


bool DataInfo::IsItemMandatory(const string& itemName)
{
    const vector<string>& mCode = GetItemAttribute(itemName,
      CifString::CIF_DDL_CATEGORY_ITEM,
      CifString::CIF_DDL_ITEM_MANDATORY_CODE);

    if (String::IsCiEqual(mCode[0], "Y") || String::IsCiEqual(mCode[0], "YES"))
    {
        return (true);
    }
    else
    {
        return (false);
    }
}


bool DataInfo::IsSimpleDataType(const string& itemName)
{
    bool simpleTyping = false;

    bool iRange = false;
    bool iUnits = false;

    const vector<string>& enums = GetItemAttribute(itemName,
      CifString::CIF_DDL_CATEGORY_ITEM_ENUMERATION,
      CifString::CIF_DDL_ITEM_VALUE);

    const vector<string>& units = GetItemAttribute(itemName,
      CifString::CIF_DDL_CATEGORY_ITEM_UNITS,
      CifString::CIF_DDL_ITEM_CODE);

    if (!units.empty())
    {
      iUnits = true;
    }

    const vector<string>& rangeMin = GetItemAttribute(itemName,
      CifString::CIF_DDL_CATEGORY_ITEM_RANGE,
      CifString::CIF_DDL_ITEM_MINIMUM);

    const vector<string>& rangeMax = GetItemAttribute(itemName,
      CifString::CIF_DDL_CATEGORY_ITEM_RANGE,
      CifString::CIF_DDL_ITEM_MAXIMUM);

    if ((!rangeMin.empty()) && (!rangeMax.empty()) &&
      (rangeMin.size() == rangeMax.size()))
    {
        iRange = true;
    }

    if (!iRange && enums.empty() && !iUnits)
    {
        simpleTyping = true;
    }

    return (simpleTyping);
}


void DataInfo::StandardizeEnumItem(string& value,
  const string& catName, const string& attribName)
{
    // If value is empty, no need to go and search. Return.
    if (value.empty())
        return;

    string itemName;
    CifString::MakeCifItem(itemName, catName, attribName);

    const vector<string>& enums = GetItemAttribute(itemName,
      CifString::CIF_DDL_CATEGORY_ITEM_ENUMERATION,
      CifString::CIF_DDL_ITEM_VALUE);

    for (const auto & i : enums)
    {
        if (String::IsCiEqual(value, i))
        {
            value = i;
            break;
        }
    }
}

