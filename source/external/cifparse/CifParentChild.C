/*
FILE:     CifParentChild.C
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


/*!
** \file CifParentChild.C
**
** \brief Implementation file for CifFile class.
*/


#include <string>
#include <vector>
#include <sstream>
#include <iostream>
#include <algorithm>

#include "GenString.h"
#include "GenCont.h"
#include "CifString.h"
#include "ISTable.h"
#include "TableFile.h"
#include "CifParentChild.h"


using std::string;
using std::vector;
using std::map;
using std::ostringstream;
using std::cout;
using std::endl;
using std::sort;


CifParentChild::CifParentChild(Block& block) : ParentChild(),
  _parChildTableP(nullptr), _inParChildGroupP(nullptr), _inParChildGroupListP(nullptr)
{
    _parChildTableP = block.GetTablePtr("item_linked");

    _inParChildGroupP = block.GetTablePtr("pdbx_item_linked_group");
    _inParChildGroupListP = block.GetTablePtr("pdbx_item_linked_group_list");

    if ((_inParChildGroupP == nullptr) && (_inParChildGroupListP != nullptr))
    {
        throw EmptyValueException("Empty item linked group table",
          "CifParentChild::CifParentChild");
    }

    if ((_inParChildGroupP != nullptr) && (_inParChildGroupListP == nullptr))
    {
        throw EmptyValueException("Empty item linked group list table",
          "CifParentChild::CifParentChild");
    }

    if (_inParChildGroupP != nullptr)
        *_groupTableP = *_inParChildGroupP;
    if (_inParChildGroupListP != nullptr)
        *_groupListTableP = *_inParChildGroupListP;

    Init(block);
}


CifParentChild::CifParentChild(Block& block, ISTable* parChildTableP) :
  ParentChild(), _parChildTableP(parChildTableP)
{
    // This method is used when it is needed to build group and grup list
    // tables out of item linked tables.

    ISTable* categoryP = block.GetTablePtr("category");
    ISTable* itemP = block.GetTablePtr("item");

    vector<string> cats;
    for (unsigned int rowI = 0; rowI < categoryP->GetNumRows(); ++rowI)
    {
        cats.push_back((*categoryP)(rowI, "id"));
    }

    sort(cats.begin(), cats.end());

    vector<vector<string> > items;

    for (const auto & cat : cats)
    {
        vector<string> searchValue;
        searchValue.push_back(cat);
   
        vector<string> searchCols;
        searchCols.emplace_back("category_id");
  
        vector<unsigned int> found;
        itemP->Search(found, searchValue, searchCols);

        vector<string> cifItemNames;
        for (unsigned int foundI : found)
        {
           cifItemNames.push_back((*itemP)(foundI, "name"));
        } 

        items.push_back(cifItemNames);
    }

    BuildOldTables(cats, items);

    AddParentCategoryToItemLinkedGroup(*_groupTableP, *_groupListTableP);

    CreateAllRelations(*_groupTableP, *_groupListTableP);
}


void CifParentChild::Init(Block& block)
{
    ISTable* categoryP = block.GetTablePtr("category");
    ISTable* itemP = block.GetTablePtr("item");

    vector<string> cats;
    for (unsigned int rowI = 0; rowI < categoryP->GetNumRows(); ++rowI)
    {
        cats.push_back((*categoryP)(rowI, "id"));
    }

    sort(cats.begin(), cats.end());

    vector<vector<string> > items;

    for (const auto & cat : cats)
    {
        vector<string> searchValue;
        searchValue.push_back(cat);
   
        vector<string> searchCols;
        searchCols.emplace_back("category_id");
  
        vector<unsigned int> found;
        itemP->Search(found, searchValue, searchCols);

        vector<string> cifItemNames;
        for (unsigned int foundI : found)
        {
           cifItemNames.push_back((*itemP)(foundI, "name"));
        } 

        items.push_back(cifItemNames);
    }

    AddParentCategoryToItemLinkedGroup(*_groupTableP, *_groupListTableP);

    BuildNewTables(cats, items);

    CreateAllRelations(*_groupTableP, *_groupListTableP);
}


CifParentChild::~CifParentChild()
= default;


void CifParentChild::BuildOldTables(const vector<string>& cats, 
  const vector<vector<string> >& items)
{
    for (unsigned int catI = 0; catI < cats.size(); ++catI)
    {
        unsigned int groupNum = 1;

        const vector<string>& cifItemNames = items[catI];

        map<string, unsigned int> maxKeyGroups;
        ISTable* keysTableP = CreateKeysTableOld(cifItemNames, maxKeyGroups);

        if (keysTableP->GetNumRows() == 0)
        {
            // No parents for this child
            delete (keysTableP);
            continue;
        }

        for (auto & maxKeyGroup : maxKeyGroups)
        {
            // pos->first is parent name
            // pos->second is number of parent's groups
            for (unsigned int groupI = 0; groupI < maxKeyGroup.second; ++groupI,
              ++groupNum)
            {
                vector<string> newGroupRow;
                newGroupRow.push_back(cats[catI]);
                newGroupRow.push_back(String::IntToString(groupNum));
                newGroupRow.push_back(cats[catI] + ":" + maxKeyGroup.first + ":" +
                  String::IntToString(groupNum));
                newGroupRow.push_back(CifString::InapplicableValue);
                newGroupRow.push_back(CifString::InapplicableValue);
                _groupTableP->AddRow(newGroupRow);

                vector<string> searchCol;
                searchCol.emplace_back("keyGroup");
                searchCol.emplace_back("parCategory");

                vector<string> searchVal;
                searchVal.push_back(String::IntToString(groupI + 1));
                searchVal.push_back(maxKeyGroup.first);

                vector<unsigned int> found;
                keysTableP->Search(found, searchVal, searchCol);

                for (unsigned int foundI : found)
                {
                    vector<string> newGroupListRow;
                    newGroupListRow.push_back(cats[catI]);
                    newGroupListRow.\
                      push_back(String::IntToString(groupNum));
                    newGroupListRow.push_back\
                      ((*keysTableP)(foundI, "childKeyCifItem"));
                    newGroupListRow.push_back\
                      ((*keysTableP)(foundI, "parKeyCifItem"));
                    newGroupListRow.push_back(maxKeyGroup.first);
                    _groupListTableP->AddRow(newGroupListRow);
                } // (for every parent/child pair in a group)
            } // for (all group numbers in a parent)
        } // for (all parents)

        delete (keysTableP); 
    } // for (all categories in the dictionary, acting as child categories)
}


void CifParentChild::BuildNewTables(const vector<string>& cats, 
  const vector<vector<string> >& items)
{
    if (_parChildTableP == nullptr)
    {
        return;
    }

    // Maps a pair of child category/parent category to the last group
    // number used.

    map<string, unsigned int> maxGroupNum;

    for (unsigned int catI = 0; catI < cats.size(); ++catI)
    {
        const vector<string>& cifItemNames = items[catI];

        vector<string> childSearchCol;
        childSearchCol.emplace_back("child_name");

        for (const auto & cifItemName : cifItemNames)
        {
            // Check if item has already been processed
            vector<string> childTarget;
            childTarget.push_back(cifItemName);

            if (_groupListTableP != nullptr)
            {
                unsigned int row =
                  _groupListTableP->FindFirst(childTarget,
                  childSearchCol);

                if (row != _groupListTableP->GetNumRows())
                    continue;
            }

            vector<unsigned int> parLoc;
            _parChildTableP->Search(parLoc, childTarget, childSearchCol);

            for (unsigned int parI : parLoc)
            {
                const string& parentItem = (*_parChildTableP)(parI,
                  "parent_name");

                string parCatName;
                CifString::GetCategoryFromCifItem(parCatName, parentItem);

                if (maxGroupNum[cats[catI]] == 0)
                {
                    unsigned int lastGroupNum = LastGroupNum(cats[catI]);
                    maxGroupNum[cats[catI]] = lastGroupNum + 1;
                }
                else
                    maxGroupNum[cats[catI]] += 1;

                cout << "Info: Creating a new group \"" <<
                  String::IntToString(maxGroupNum[cats[catI]]) <<
                  "\" for child category \"" << cats[catI] <<
                  "\" for child item \"" << cifItemName <<
                  "\" and parent item \"" << parentItem <<
                  "\", from \"item_linked\" table, as these "\
                  "are not defined in group tables." << endl;

                vector<string> newGroupRow;
                newGroupRow.push_back(cats[catI]);
                newGroupRow.push_back(String::IntToString(
                  maxGroupNum[cats[catI]]));
                newGroupRow.push_back(cats[catI] + ":" + parCatName + ":" +
                  String::IntToString(maxGroupNum[cats[catI]]));
                newGroupRow.push_back(CifString::InapplicableValue);
                newGroupRow.push_back(CifString::InapplicableValue);
                newGroupRow.push_back(parCatName);
                _groupTableP->AddRow(newGroupRow);

                vector<string> row;
                row.push_back(cats[catI]);
                row.push_back(String::IntToString(
                  maxGroupNum[cats[catI]]));
                row.push_back(cifItemName);
                row.push_back(parentItem);
                row.push_back(parCatName);
                _groupListTableP->AddRow(row);
            } // for (all child's parent items)
        } // for (all child items)
    }
}


int CifParentChild::CheckParentChild(Block& block, ISTable& catTable,
  ostringstream& log)
{
    int ret = 1;

    const string& childCatName = catTable.GetName();

    vector<vector<string> > parParKeys; 
    vector<vector<string> > comboComboKeys;

    GetParents(parParKeys, comboComboKeys, childCatName);


    if (parParKeys.empty())
    {
#ifdef JW_DEBUG
        cout << "No parents to check for category " <<  catTable.GetName() << endl;
#endif
        // No parents. Just return.
        return(ret);
    }

    vector<string> cifItemNames;
    CifString::MakeCifItems(cifItemNames, catTable.GetName(),
      catTable.GetColumnNames());

    FilterMissingItems(parParKeys, comboComboKeys, cifItemNames);

    for (unsigned int allParI = 0; allParI < parParKeys.size(); ++allParI)
    {
        vector<string>& parKeys = parParKeys[allParI];
        vector<string>& comboKeys = comboComboKeys[allParI];

        string parCatName;
        CifString::GetCategoryFromCifItem(parCatName, parKeys[0]);

#if JW_DEBUG
        cout << "Current category "    <<  catTable.GetName() << endl;
        cout << "Parent category  "    << parCatName << endl;
        cout << "Max group number "    << pos->second << endl;
        //
        for ( unsigned int ii =0; ii < parKeys.size(); ii++)
        {
            cout << "parent  # ii " << ii << " " << parKeys[ii] << endl;
        }

        for ( unsigned int ii =0; ii < comboKeys.size(); ii++)
        {
            for (unsigned int jj =0; jj < comboKeys[ii].size(); jj++)
            {
                cout << " comboKeys[][] " << ii << " " << jj <<  " child = "
                 << comboKeys[ii][jj] << endl;
            }
        }
#endif

        vector<string> parKeyItems;
        for (const auto & parKey : parKeys)
        {
            string tmpItem;
            CifString::GetItemFromCifItem(tmpItem, parKey);
            parKeyItems.push_back(tmpItem);
        }

#ifdef VLAD_PERF
        cout << "Processing parent \"" << parCatName << "\" combo key: " <<
          comboKeyI << " of " << comboKeys.size() << endl;
#endif
#ifdef JW_DEBUG
        cout << "Processing parent \"" << parCatName << "\" combo key: "
             << comboKeyI + 1 << " of " << comboKeys.size() <<  endl;

        cout << "In category " << catTable.GetName() << " row count  "
             <<  catTable.GetNumRows() << endl;
#endif
#ifdef JW_HACK
        int iKeyDif;
        bool allowMissing = false;
#endif
        for (unsigned int l = 0; l < catTable.GetNumRows(); ++l)
        {
#ifdef JW_DEBUG
            cout << "++ In category " << catTable.GetName() << " row: " << l
              << " of " << catTable.GetNumRows() << endl;
#endif
#ifdef VLAD_PERF
            cout << "  Processing combo key row: " << l << " of " <<
              catTable.GetNumRows() << endl;
#endif
            vector<string> parKeyNonEmptyValues;
            vector<string> parKeyAllValues;
            vector<string> parKeyNonEmptyItems;
            vector<string> childKeyNonEmptyItems;
#ifdef JW_HACK
            vector<string> skippedParKeyItems;
            iKeyDif=0;
            allowMissing=false;
#endif
            for (unsigned int keyI = 0; keyI < comboKeys.size(); ++keyI)
            {
                string itemName;
                CifString::GetItemFromCifItem(itemName, comboKeys[keyI]);

                // loop for all rows
                const string& parKeyValue = catTable(l, itemName);
                parKeyAllValues.push_back(parKeyValue);
                if (CifString::IsEmptyValue(parKeyValue))
                {
#ifdef JW_HACK
                    // Skip search for null/na/missing values
                    skippedParKeyItems.push_back(parKeyItems[keyI]);
#endif
                }
                else
                {
                    parKeyNonEmptyValues.push_back(parKeyValue);
                    parKeyNonEmptyItems.push_back(parKeyItems[keyI]);
                    childKeyNonEmptyItems.push_back(comboKeys[keyI]);
                }
            } // For every key item in a combo key

            if (parKeyNonEmptyValues.empty())
                continue;

#ifndef VLAD_FIX_SINGLE_KEY_UNKNOWN_VALUE_NOT_REPORT_PARENT
            if (!block.IsTablePresent(parCatName))
            {
                // Parent table does not exist. Ignore this fact.
#ifdef JW_HACK
                if (parCatName != "chem_comp_atom" &&
                 parCatName != "atom_sites_alt")
                {
                   log << "++ERROR - In block \"" << block.GetName() <<
                     "\", parent category \"" << parCatName <<
                     "\", of category \"" << catTable.GetName() <<
                     "\", is missing." << endl;
                }
#endif
                break;
            }

            ISTable* parentCatTableP = block.GetTablePtr(parCatName);
            if (parentCatTableP->GetNumRows() == 0)
            {
#ifdef JW_HACK
                log << "WARN - In block " << block.GetName() << " category " <<
                  catTable.GetName() << " parent table empty " << parCatName <<
                 endl;
#endif
                 // Parent table has no rows. Ignore this fact.
                 break;
            }

            // Check existence in the parent
            for (unsigned int parKeyI = 0; parKeyI < parKeyItems.size();
              ++parKeyI)
            {
               if (!CifString::IsEmptyValue(parKeyAllValues[parKeyI]) &&
                 !parentCatTableP->IsColumnPresent(parKeyItems[parKeyI]))
               {
                   log << "ERROR - In block \"" << block.GetName() <<
                     "\", for child item \"" <<
                     comboKeys[parKeyI] <<
                     "\", parent key attribute \"" <<
                     parKeyItems[parKeyI] <<
                     "\" is not defined in parent category \"" <<
                     parCatName << "\"" << endl;
                   ret = 0;

                   return(ret);
               }
           }
#endif

#ifdef JW_HACK
            iKeyDif = comboKeys.size() - parKeyNonEmptyValues.size();
            //
            // HACK - to deal with real missing values - 
            //
            if (iKeyDif == 1)
            {
                if (
		      ((catTable.GetName() == "atom_site") && (parCatName == "entity_poly_seq") && (skippedParKeyItems[0] == "num")) ||
		      ((catTable.GetName() == "atom_site") && (parCatName == "pdbx_poly_seq_scheme") && (skippedParKeyItems[0] == "seq_id")) 

		      )
                {
                    allowMissing = true;
		}

#ifdef JW_DEBUG
		if (allowMissing)
                {
                    cout << "++INFO - Allowed exception (keydif=" << iKeyDif << ") " 
			 << " parent " << parCatName 
			 << " skipped parent item(s) " << skippedParKeyItems[0]
			 << " in child " 
			 << catTable.GetName() << " row: " << l+1 << " of " 
			 << catTable.GetNumRows() << endl;
		  }
#endif

            }

            if (iKeyDif == 2)
            {
                  if ((catTable.GetName() == "pdbx_poly_seq_scheme") && (parCatName == "atom_site") &&
		      ((skippedParKeyItems[0] == "auth_comp_id" &&  skippedParKeyItems[1] == "pdbx_PDB_ins_code") ||
		       (skippedParKeyItems[1] == "auth_comp_id" &&   skippedParKeyItems[0] == "pdbx_PDB_ins_code")))
		    {
		      allowMissing = true;
		    } else if ((catTable.GetName() == "atom_site") && (parCatName == "pdbx_poly_seq_scheme") &&
		      ((skippedParKeyItems[0] == "seq_id" &&   skippedParKeyItems[1] == "pdb_ins_code") ||
		       (skippedParKeyItems[1] == "seq_id" &&   skippedParKeyItems[0] == "pdb_ins_code"))) 
		    {
		      allowMissing = true;
		    }
		  
#ifdef JW_DEBUG
		    if (allowMissing) {
		      cout << "++INFO - Allowed exception (keydif=" << iKeyDif << ") " 
			   << " parent " << parCatName
			   << " skipped parent item(s) " << skippedParKeyItems[0]
			   << " and " << skippedParKeyItems[1]
			   << " in child " 
			   << catTable.GetName() << " row: " << l+1 << " of " 
			   << catTable.GetNumRows() << endl;
		    }
#endif
            }

            if (iKeyDif == 3)
            { 
                  if ( (catTable.GetName() == "pdbx_poly_seq_scheme") && (parCatName == "atom_site") &&
		       ((skippedParKeyItems[0] == "auth_seq_id" &&
			 skippedParKeyItems[1] == "auth_comp_id" &&
			 skippedParKeyItems[2] == "pdbx_PDB_ins_code") ||
			(skippedParKeyItems[0] == "auth_seq_id"  &&
			 skippedParKeyItems[1] == "auth_comp_id" &&
			 skippedParKeyItems[2] == "pdbx_PDB_ins_code")))  {
		    allowMissing = true;
		  } else if ((catTable.GetName() == "struct_ref_seq_dif") && (parCatName == "pdbx_poly_seq_scheme") &&
		       ((skippedParKeyItems[0] == "mon_id" &&
			 skippedParKeyItems[1] == "seq_id" &&
			 skippedParKeyItems[2] == "pdb_ins_code") ||
			(skippedParKeyItems[1] == "mon_id"  &&
			 skippedParKeyItems[0] == "seq_id" &&
			 skippedParKeyItems[2] == "pdb_ins_code")))  {
		    allowMissing = true;
		  }
#ifdef JW_DEBUG
		  if (allowMissing) {		   
                    cout << "++INFO - Allowed exception (keydif=" << iKeyDif << ") " 
			 << " parent " << parCatName 
			 << " skipped parent item(s) " << skippedParKeyItems[0]
			 << " and " << skippedParKeyItems[1]
			 << " and " << skippedParKeyItems[2]
			 << " in child " 
			 << catTable.GetName() << " row: " << l+1 << " of " 
			 << catTable.GetNumRows() << endl;
		  }
#endif
            }

            // END HACK -


#ifdef JW_DEBUG
            if (iKeyDif != 0 && !allowMissing)
            {
                cout << "++INFO - child category (keydif=" << iKeyDif << ") " << catTable.GetName() <<
                      " row: " << l << " of " << catTable.GetNumRows() << endl;
            }

            if (iKeyDif > 0)
            {
                if (allowMissing)
                {
                    for ( unsigned int ii =0; ii < skippedParKeyItems.size(); ii++)
                    {
                        cout <<
                          "++INFO - Allowed incomplete parent key in " << parCatName << " "
                          << ii+1 << " of " << skippedParKeyItems.size() <<
                          " -> " << skippedParKeyItems[ii] << endl;
                    }
                }
                else
                {
                    for (unsigned int ii =0;
                      ii < skippedParKeyItems.size(); ii++)
                    {
                        cout << "++INFO - Incomplete parent key in " 
                          << parCatName << " " 
                          << ii+1 << " of " << skippedParKeyItems.size() <<
                          " -> " << skippedParKeyItems[ii] << endl;
                    }
                }
            }
#endif
#endif

#ifdef JW_HACK
            if ((iKeyDif == 0) || ((iKeyDif > 0) && !allowMissing))
#endif

            {
                    unsigned int searchIndex =
                      parentCatTableP->FindFirst(parKeyNonEmptyValues,
                        parKeyNonEmptyItems);
                    if (searchIndex == parentCatTableP->GetNumRows())
                    {
#ifdef OLD_IMPL
                        log << "ERROR - In block \"" << block.GetName()
                          << " unmatched value in " << catTable.GetName()
                          << " row " << l+1 <<  " in parent category " <<
                          parCatName  << endl;

                        for (unsigned int keyValI = 0;
                          keyValI < parKeyNonEmptyValues.size(); ++keyValI)
                        {
                            log << "  \"" << childKeyNonEmptyItems[keyValI]
                              << "\" -> \"_" << parCatName << "." <<
                              parKeyNonEmptyItems[keyValI] << "\" value =\"" <<
                              parKeyNonEmptyValues[keyValI] << "\"" << endl;
                        }
#else
                        string linkGroupIdLabel;
                        GetLinkGroupIdLabel(linkGroupIdLabel, parKeys,
                          comboKeys);
                        if (linkGroupIdLabel.empty())
                        {
                            log << "BIG TROUBLE; linkGroupId not found!" <<
                              endl;
                        }

                        log << "ERROR PARCHILD \"" << linkGroupIdLabel <<
                          "\" - In block \"" << block.GetName() <<
                          "\", in category \"" << childCatName <<
                          "\", in row# " << l + 1 <<
                          ", unmatched value in the parent \"" << parCatName <<
                          "\"";

                        log << endl;

#ifdef VLAD_MORE_CORRECTIONS
                        log <<  "(";
    
                        for (unsigned int keyValI = 0;
                          keyValI < parKeyNonEmptyValues.size(); ++keyValI)
                        {
                            log << "\"" << childKeyNonEmptyItems[keyValI] 
                              << "\" -> \"_" << parCatName << "." <<
                              parKeyNonEmptyItems[keyValI] << "\" value = \"" <<
                              parKeyNonEmptyValues[keyValI] << "\"";
                            if (keyValI != parKeyNonEmptyValues.size() - 1)
                                log << ", ";
                        }

                        log << ")" << endl;
#endif

                        for (unsigned int keyValI = 0;
                          keyValI < parKeyNonEmptyValues.size(); ++keyValI)
                        {
                            log << "  \"" << childKeyNonEmptyItems[keyValI]
                              << "\" -> \"_" << parCatName << "." <<
                              parKeyNonEmptyItems[keyValI] << "\" value =\"" <<
                              parKeyNonEmptyValues[keyValI] << "\"" << endl;
                        }
#endif


#ifdef JW_DEBUG
			cout << "ERROR - In block \"" << block.GetName() 
			     << " unsatisfied child value in " << catTable.GetName() 
			     << " row " << l+1 <<  " parent category " << parCatName 
			     << endl;
			
			cout << "ERROR - full key \"" << comboKeys.size()
			     << " realized key size " <<
			  parKeyNonEmptyValues.size() << endl;
			
			for (unsigned int keyValI = 0; keyValI < parKeyNonEmptyValues.size(); ++keyValI)
			  {
			    cout << "  \"" << childKeyNonEmptyItems[keyValI] 
				 << "\"->\"" << parCatName << "." <<
			      parKeyNonEmptyItems[keyValI] << "\" value =\"" <<
			      parKeyNonEmptyValues[keyValI] << "\"" << endl;
			  }

#endif // JW_DEBUG
                    }
                    ret = 0;
            }
        } // For every row in child
    } // For every parent

    return(ret);
}


void CifParentChild::WriteGroupTables(Block& block)
{
    auto* tmpGroupTableP = new ISTable();
    *tmpGroupTableP = *_groupTableP;

    auto* tmpGroupListTableP = new ISTable();
    *tmpGroupListTableP = *_groupListTableP;

    tmpGroupTableP->DeleteColumn("parent_category_id");

    if (tmpGroupTableP->GetNumRows() != 0)
        block.WriteTable(tmpGroupTableP);
    if (tmpGroupListTableP->GetNumRows() != 0)
        block.WriteTable(tmpGroupListTableP);
}


ISTable* CifParentChild::CreateKeysTableOld(const vector<string>& cifItemNames,
  map<string, unsigned int>& maxKeyGroups)
{
    // Prepare the table that will contain all the needed keys information
    // This table will have four columns: child key CIF item, group number
    // parent key CIF item and parent category name

    auto* keysTableP = new ISTable();

    keysTableP->AddColumn("childKeyCifItem");
    keysTableP->AddColumn("keyGroup");
    keysTableP->AddColumn("parKeyCifItem");
    keysTableP->AddColumn("parCategory");

    FillKeysTableOld(*keysTableP, cifItemNames, maxKeyGroups);

    return(keysTableP);
}


void CifParentChild::FillKeysTableOld(ISTable& keysTable,
  const vector<string>& cifItemNames, map<string, unsigned int>& maxKeyGroups)
{
    vector<string> keyList;
    keyList.emplace_back("parKeyCifItem");

    for (const auto & cifItemName : cifItemNames)
    {
        // cifItemNames are child's items
        vector<string> parCifItems;
        GetParentCifItems(parCifItems, cifItemName);

        sort(parCifItems.begin(), parCifItems.end());

        for (const auto & parCifItem : parCifItems)
        {
            string parCatName;
            CifString::GetCategoryFromCifItem(parCatName, parCifItem);

            vector<string> keyTarget;
            keyTarget.push_back(parCifItem);

            vector<unsigned int> parents;
            keysTable.Search(parents, keyTarget, keyList);

            unsigned int newKeyGroupInt = 1;
            if (!parents.empty())
            {
                newKeyGroupInt = String::StringToInt(keysTable(
                  parents[parents.size() - 1], "keyGroup"));
                ++newKeyGroupInt;
            }

            if (newKeyGroupInt > maxKeyGroups[parCatName])
                maxKeyGroups[parCatName] = newKeyGroupInt;

            string newKeyGroup = String::IntToString(newKeyGroupInt);

            vector<string> row;

            row.push_back(cifItemName);
            row.push_back(newKeyGroup);
            row.push_back(parCifItem);
            row.push_back(parCatName);

            keysTable.AddRow(row);
        } // For all item's parent items (that can be in different parents)
    } // For all items in the category
}


void CifParentChild::GetParentCifItems(vector<string>& parCifItems,
  const string& cifItemName)
{
    parCifItems.clear();

    vector<string> childCifItem;
    childCifItem.push_back(cifItemName);

    vector<string> childNameCol;
    childNameCol.emplace_back("child_name");

    vector<unsigned int> parLoc;
    _parChildTableP->Search(parLoc, childCifItem, childNameCol);

    for (unsigned int parLocI : parLoc)
    {
        parCifItems.push_back((*_parChildTableP)(parLocI,
          "parent_name"));
    }
}


unsigned int CifParentChild::LastGroupNum(const string& childCat)
{
    vector<string> searchCol;
    searchCol.emplace_back("category_id");

    vector<string> searchVal;
    searchVal.push_back(childCat);

    vector<unsigned int> found;

    unsigned int lastGroupNum = 0;

    _groupTableP->Search(found, searchVal, searchCol);

    for (unsigned int foundI : found)
    {
        unsigned int currGroupNum =
          String::StringToInt((*_groupTableP)(foundI,
          "link_group_id"));

        if (currGroupNum > lastGroupNum)
            lastGroupNum = currGroupNum;
    }

    return (lastGroupNum);
}


void CifParentChild::FilterMissingItems(vector<vector<string> >& parParKeys,
  vector<vector<string> >& comboComboKeys, const vector<string>& cifItemNames)
{
    // Identify missing items indices

    for (unsigned int allParI = 0; allParI < parParKeys.size(); ++allParI)
    {
        vector<string>& parKeys = parParKeys[allParI];
        vector<string>& comboKeys = comboComboKeys[allParI];

        vector<unsigned int> remIndices;
        for (unsigned int keysI = 0; keysI < comboKeys.size(); ++keysI)
        {
            if (!GenCont::IsInVector(comboKeys[keysI], cifItemNames))
                remIndices.push_back(keysI);
        }
 
        sort (remIndices.begin(), remIndices.end());

        for (unsigned int indI = 0; indI < remIndices.size(); ++indI)
        {
            // Delete elements from end to begin
            parKeys.erase(parKeys.begin() + remIndices[remIndices.size() -
              indI - 1]);
            comboKeys.erase(comboKeys.begin() + remIndices[remIndices.size() -
              indI - 1]);
        }
    }
}

