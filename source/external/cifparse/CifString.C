/*
FILE:     CifString.C
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


#include <string>
#include <vector>

#include "Exceptions.h"
#include "CifString.h"


using std::string;
using std::vector;


const string CifString::CIF_DDL_CATEGORY_BLOCK("datablock");
const string CifString::CIF_DDL_CATEGORY_DATABLOCK("datablock");
const string CifString::CIF_DDL_CATEGORY_DATABLOCK_METHODS("datablock_methods");
const string CifString::CIF_DDL_CATEGORY_ITEM("item");
const string CifString::CIF_DDL_CATEGORY_ITEM_LINKED("item_linked");
const string CifString::CIF_DDL_CATEGORY_PDBX_ITEM_LINKED_GROUP("pdbx_item_linked_group");
const string CifString::CIF_DDL_CATEGORY_PDBX_ITEM_LINKED_GROUP_LIST("pdbx_item_linked_group_list");
const string CifString::CIF_DDL_CATEGORY_CATEGORY("category");
const string CifString::CIF_DDL_CATEGORY_CATEGORY_EXAMPLES("category_examples");
const string CifString::CIF_DDL_CATEGORY_NDB_CATEGORY_EXAMPLES("ndb_category_examples");
const string CifString::CIF_DDL_CATEGORY_CATEGORY_KEY("category_key");
const string CifString::CIF_DDL_CATEGORY_CATEGORY_GROUP("category_group");
const string CifString::CIF_DDL_CATEGORY_CATEGORY_GROUP_LIST("category_group_list");
const string CifString::CIF_DDL_CATEGORY_CATEGORY_METHODS("category_methods");
const string CifString::CIF_DDL_CATEGORY_SUB_CATEGORY("sub_category");
const string CifString::CIF_DDL_CATEGORY_SUB_CATEGORY_EXAMPLES("sub_category_examples");
const string CifString::CIF_DDL_CATEGORY_SUB_CATEGORY_METHODS("sub_category_methods");
const string CifString::CIF_DDL_CATEGORY_ITEM_SUB_CATEGORY("item_sub_category");
const string CifString::CIF_DDL_CATEGORY_ITEM_TYPE("item_type");
const string CifString::CIF_DDL_CATEGORY_ITEM_TYPE_CONDITIONS("item_type_conditions");
const string CifString::CIF_DDL_CATEGORY_ITEM_METHODS("item_methods");
const string CifString::CIF_DDL_CATEGORY_ITEM_TYPE_LIST("item_type_list");
const string CifString::CIF_DDL_CATEGORY_ITEM_STRUCTURE("item_structure");
const string CifString::CIF_DDL_CATEGORY_ITEM_STRUCTURE_LIST("item_structure_list");
const string CifString::CIF_DDL_CATEGORY_ITEM_DESCRIPTION("item_description");
const string CifString::CIF_DDL_CATEGORY_NDB_ITEM_DESCRIPTION("ndb_item_description");
const string CifString::CIF_DDL_CATEGORY_NDB_CATEGORY_DESCRIPTION("ndb_category_description");
const string CifString::CIF_DDL_CATEGORY_ITEM_EXAMPLES("item_examples");
const string CifString::CIF_DDL_CATEGORY_NDB_ITEM_EXAMPLES("ndb_item_examples");
const string CifString::CIF_DDL_CATEGORY_ITEM_DEPENDENT("item_dependent");
const string CifString::CIF_DDL_CATEGORY_ITEM_RELATED("item_related");
const string CifString::CIF_DDL_CATEGORY_ITEM_RANGE("item_range");
const string CifString::CIF_DDL_CATEGORY_ITEM_ENUMERATION("item_enumeration");
const string CifString::CIF_DDL_CATEGORY_NDB_ITEM_ENUMERATION("ndb_item_enumeration");
const string CifString::CIF_DDL_CATEGORY_ITEM_DEFAULT("item_default");
const string CifString::CIF_DDL_CATEGORY_ITEM_ALIASES("item_aliases");
const string CifString::CIF_DDL_CATEGORY_DICTIONARY("dictionary");
const string CifString::CIF_DDL_CATEGORY_DICTIONARY_HISTORY("dictionary_history");
const string CifString::CIF_DDL_CATEGORY_ITEM_UNITS("item_units");
const string CifString::CIF_DDL_CATEGORY_ITEM_UNITS_LIST("item_units_list");
const string CifString::CIF_DDL_CATEGORY_ITEM_UNITS_CONVERSION("item_units_conversion");
const string CifString::CIF_DDL_CATEGORY_METHOD_LIST("method_list");

const string CifString::CIF_DDL_ITEM_ID("id");
const string CifString::CIF_DDL_ITEM_CATEGORY_ID("category_id");
const string CifString::CIF_DDL_ITEM_SUB_CATEGORY_ID("sub_category_id");
const string CifString::CIF_DDL_ITEM_METHOD_ID("method_id");
const string CifString::CIF_DDL_ITEM_PARENT_NAME("parent_name");
const string CifString::CIF_DDL_ITEM_CHILD_NAME("child_name");
const string CifString::CIF_DDL_ITEM_CHILD_CATEGORY_ID("child_category_id");
const string CifString::CIF_DDL_ITEM_PARENT_CATEGORY_ID("parent_category_id");
const string CifString::CIF_DDL_ITEM_LINK_GROUP_ID("link_group_id");
const string CifString::CIF_DDL_ITEM_LABEL("label");
const string CifString::CIF_DDL_ITEM_CONTEXT("context");
const string CifString::CIF_DDL_ITEM_CONDITION_ID("condition_id");
const string CifString::CIF_DDL_ITEM_ALIAS_NAME("alias_name");
const string CifString::CIF_DDL_ITEM_DICTIONARY("dictionary");
const string CifString::CIF_DDL_ITEM_TITLE("title");
const string CifString::CIF_DDL_ITEM_VERSION("version");
const string CifString::CIF_DDL_ITEM_NAME("name");
const string CifString::CIF_DDL_ITEM_CODE("code");
const string CifString::CIF_DDL_ITEM_PRIMITIVE_CODE("primitive_code");
const string CifString::CIF_DDL_ITEM_CONSTRUCT("construct");
const string CifString::CIF_DDL_ITEM_ORGANIZATION("organization");
const string CifString::CIF_DDL_ITEM_INDEX("index");
const string CifString::CIF_DDL_ITEM_DIMENSION("dimension");
const string CifString::CIF_DDL_ITEM_DATABLOCK_ID("datablock_id");
const string CifString::CIF_DDL_ITEM_DESCRIPTION("description");
const string CifString::CIF_DDL_ITEM_NDB_DESCRIPTION("ndb_description");
const string CifString::CIF_DDL_ITEM_CASE("case");
const string CifString::CIF_DDL_ITEM_MANDATORY_CODE("mandatory_code");
const string CifString::CIF_DDL_ITEM_DETAIL("detail");
const string CifString::CIF_DDL_ITEM_MAXIMUM("maximum");
const string CifString::CIF_DDL_ITEM_MINIMUM("minimum");
const string CifString::CIF_DDL_ITEM_VALUE("value");
const string CifString::CIF_DDL_ITEM_DEPENDENT_NAME("dependent_name");
const string CifString::CIF_DDL_ITEM_RELATED_NAME("related_name");
const string CifString::CIF_DDL_ITEM_FUNCTION_CODE("function_code");
const string CifString::CIF_DDL_ITEM_OFFSET("file_offset");
const string CifString::CIF_DDL_ITEM_OPERATOR("operator");
const string CifString::CIF_DDL_ITEM_FACTOR("factor");
const string CifString::CIF_DDL_ITEM_FROM_CODE("from_code");
const string CifString::CIF_DDL_ITEM_TO_CODE("to_code");
const string CifString::CIF_DDL_ITEM_UPDATE("update");
const string CifString::CIF_DDL_ITEM_REVISION("revision");
const string CifString::CIF_DDL_ITEM_INLINE("item_inline");
const string CifString::CIF_DDL_ITEM_LANGUAGE("item_language");
const string CifString::CIF_DDL_ITEM_PARENT_ID("parent_id");


const string CifString::UnknownValue("?");
const string CifString::InapplicableValue(".");


void CifString::MakeCifItem(string& cifItem, const string& categoryName,
  const string& attribName)
{
    if (categoryName.empty())
    {
        throw EmptyValueException("Empty category name",
          "CifString::MakeCifItem");
    }

    if (attribName.empty())
    {
        throw EmptyValueException("Empty attribute name",
          "CifString::MakeCifItem");
    }

    cifItem = PREFIX_CHAR + categoryName + JOIN_CHAR + attribName;
}


void CifString::MakeCifItems(vector<string>& cifItems,
  const string& categoryName, const vector<string>& attribsNames)
{
    cifItems.assign(attribsNames.size(), string());

    for (unsigned int attribI = 0; attribI < attribsNames.size(); ++attribI)
    {
        MakeCifItem(cifItems[attribI], categoryName, attribsNames[attribI]);
    }
}


void CifString::GetItemFromCifItem(string& itemName, const string& cifItem)
{
    /* ----------------------------------------------------------------------
     Purpose: CifString::GetItemFromCifItem(itemName, cifItem)

     Get the itemName part of an item name (ie. _<category>.<itemName>)
     Return 1 for success or 0 otherwise.

     *--------------------------------------------------------------------- */

    unsigned int k, ilen;

    if (cifItem.empty() || (cifItem[0] != '_') || (ilen = cifItem.size()) < 4)
    {
        throw EmptyValueException("Invalid CIF item \"" + cifItem + "\"",
          "CifString::GetItemFromCifItem");
    }

    for (k = 1; k < ilen; ++k)
    {
        if (cifItem[k] == JOIN_CHAR)
            break;
    }

    if ((k == 1) || (k == (ilen - 1)))
    {
        throw EmptyValueException("Invalid CIF item \"" + cifItem + "\"",
          "CifString::GetItemFromCifItem");
    }

    itemName.clear();

    for (unsigned int i = k + 1; i < ilen; ++i)
    {
        itemName.push_back(cifItem[i]);
    }
}


void CifString::GetCategoryFromCifItem(string& categoryName,
    const string& cifItem)
/* ----------------------------------------------------------------------
 Purpose: CifString::GetCategoryFromCifItem(categoryName, cifItem)

 Get the category part of an item name (ie. _<category>.<keyword>)
 Return 1 for success or 0 otherwise.

 * ---------------------------------------------------------------------- */
{
    if (cifItem.empty() || (cifItem[0] != '_'))
    {
        throw EmptyValueException("Invalid CIF item \"" + cifItem + "\"",
          "CifString::GetCategoryFromCifItem");
    }

    // Skip the first char and search for a join character
    string::size_type dotIndex = cifItem.find(JOIN_CHAR, 1);

    if (dotIndex == string::npos)
    {
        throw EmptyValueException("Invalid CIF item \"" + cifItem + "\"",
          "CifString::GetCategoryFromCifItem");
#ifdef VLAD_DEL
        // If join char not found, set its index to be the end of the string.
        dotIndex = cifItem.size();
#endif
    }

    categoryName.clear();

    // Copy from the first char to the character prior to the dot char.
    categoryName.assign(cifItem, 1, dotIndex - 1);
}


bool CifString::IsEmptyValue(const string& value)
{
    if (value.empty() || (value == InapplicableValue) ||
      (value == UnknownValue))
    {
        return (true);
    }
    else
    {
        return (false);
    }
}


bool CifString::IsUnknownValue(const string& value)
{
    if (value.empty() || (value == UnknownValue))
    {
        return (true);
    }
    else
    {
        return (false);
    }
}


bool CifString::IsSpecialChar(const char charValue)
{
    switch (charValue)
    {
        case '(':
        case ')':
        case '[':
        case ']':
        case '{':
        case '}':
        {
            return (true);
            break;
        }
        default:
        {
            return (false);
            break;
        }
    }
}


bool CifString::IsSpecialFirstChar(const char charValue)
{
    switch (charValue)
    {
        case '$':
        case '#':
        case '_':
        case ';':
        case '(':
        case ')':
        case '[':
        case ']':
        case '{':
        case '}':
        {
            return (true);
            break;
        }
        default:
        {
            return (false);
            break;
        }
    }
}

