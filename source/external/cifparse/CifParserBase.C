/*
FILE:     CifParserBase.C
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
** \file CifParserBase.C
**
** \brief Implementation file for CifParser class.
*/


/* 
  PURPOSE:    A DDL 2.1 compliant CIF file parser.
*/

#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "Exceptions.h"
#include "GenString.h"
#include "RcsbFile.h"
#include "CifString.h"
#include "CifScannerBase.h"
#include "CifParserInt.h"
#include "CifParserBase.h"

using CIFPARSER__BUFFER_STATE = struct cifparser__buffer_state *;

// WIN32 - move inside extern "C"
extern "C" FILE* cifparser_in;

extern "C"
{
    int cifparser_parse();
    void cifparser_restart(FILE*);
    CIFPARSER__BUFFER_STATE cifparser__scan_string(const char*);
    void cifparser__delete_buffer(CIFPARSER__BUFFER_STATE);

// WIN32 - move inside extern "C"

char* Glob_tBufKeyword;
char* Glob_pBufValue;
char* Glob_dataBlockName = nullptr;

}

CifParser* CifParserP = nullptr;

#ifdef PDB_TRACE
using std::cout;
#endif
using std::endl;

CifParser::CifParser(CifFile *fo, bool verbose)
{

    if (CifParserP != nullptr)
    {
        // Attempting to create a new parser, during the lifetime of
        // an existing parser.
        throw AlreadyExistsException("Cannot create a new parser, since "\
          "one already exists.", "CifParser::CifParser");
    }

    if (fo != nullptr)
        _fobj = fo;
    else
        throw EmptyValueException("fo is a NULL pointer",
          "CifParser::CifParser");

    _verbose=verbose;

    Clear();

    CifParserP = this;
}


CifParser::CifParser(CifFile *fo, CifFileReadDef readDef, bool verbose)
{

    if (CifParserP != nullptr)
    {
        // Attempting to create a new parser, during the lifetime of
        // an existing parser.
        throw AlreadyExistsException("Cannot create a new parser, since "\
          "one already exists.", "CifParser::CifParser");
    }

    if (fo != nullptr)
        _fobj = fo;
    else
        throw EmptyValueException("fo is a NULL pointer",
          "CifParser::CifParser");

    _verbose=verbose;

    Clear();
    _readDef = readDef;

    CifParserP = this;
}

void CifParser::Parse(const string& fileName, string& diagnostics,
  const string& parseLogFileName)
{
    diagnostics.clear();

    FILE* cifIn;

    if ((cifIn = fopen(fileName.c_str(), "r")) == nullptr)
    {
      diagnostics = "Unable to open file.";

      throw NotFoundException("File \"" + fileName + "\" cannot be opened",
        "CifParser::Parse");
    }

    string logFileName;
    if (!parseLogFileName.empty())
    {
        logFileName = parseLogFileName;
    }
    else
    {
        RcsbFile::RelativeFileName(logFileName, fileName);

        logFileName += "-parser.log";
    }

    OpenLog(logFileName, _verbose);

    cifparser_in = cifIn;

    int ret;

    cifparser_restart(cifparser_in);

    ret = cifparser_parse();
    if (ret != 0)
    {
        int b = 0;
        b++;
    }

    fclose(cifIn);

    if (RcsbFile::IsEmpty(log))
    {
        log.close();
        RcsbFile::Delete(logFileName);
    }
    else
    {
        log.close();
    }

    if (this->errorLog.size() > 0)
    {
        diagnostics = this->errorLog;
    }
}

void CifParser::ParseString(const string& cifString, string& diagnostics)
{

    diagnostics.clear();

    CIFPARSER__BUFFER_STATE bufferState =
      cifparser__scan_string(cifString.c_str());

    int ret;

    ret = cifparser_parse();

    cifparser__delete_buffer(bufferState);

    if (ret != 0)
    {
        int b = 0;
        b++;
    }

    if (this->errorLog.size() > 0)
    {
        diagnostics = this->errorLog;
    }

}

void cifparser_error(const char *s)
{
    CifParserP->Error(s);
}

void CifParser::Error(const char* s)
/* 
 * Purpose:  yyerror() Print errors in CifParsr.log.
 */
{
  errorLog += "ERROR - ";
  errorLog += s;
  errorLog += " at line ";
  errorLog += String::IntToString(NDBlineNo);
#if DEBUG
  errorLog += "   ";
  errorLog += yytext;
#endif
  errorLog += '\n';
  log <<"ERROR - "<< s << " at line " <<  NDBlineNo << endl;
  _err++;
}

void CifParser::Reset() 
{
  if (_tbl && (_curItemNo > 0)) { // write the current table / management of _tbl by _fobj
      _ComplexWriteTable();
  }
  _fieldList.clear();
}

void CifParser::Clear()
{
  _nTablesInBlock=0;
  _numDataBlocks=0;
  _afterLoop = false;
  _curItemNo=0;
  _curRow=0;
  _curValueNo=0;
  _fieldListAlloc=100;
  _fieldList.reserve(_fieldListAlloc);
  _pBufValue.clear();
  _tbl=nullptr;
  _err=0;
  _warn=0;
  _tBufKeyword.clear();
  _curCategoryName.clear();
  _curDataBlockName.clear();
  _prevDataBlockName = _curDataBlockName;
}


int CifParser::ProcessLoopDeclaration(void)
/* ----------------------------------------------------------------------
     Purpose: CifParser::ProcessLoopDeclaration(void)

              Handles initialization for a new loop, by creating a new 
              category and adding the current item name to this category.
   ---------------------------------------------------------------------- */
{
  string categoryName;

#if DEBUG
  if (_verbose != 0)
  {
    log << "Processing loop declaration at line " << NDBlineNo <<
      " value " << _tBufKeyword << endl;
  }
#endif

  _afterLoop = true;

  // Write the previously processed table, if it exists
  if (_tbl != nullptr)
  {
    _ComplexWriteTable();
    // VLAD: Should _tbl be destroyed here first, prior to assigning NULL?
    _tbl = nullptr;
    _curCategoryName.clear();
  }

  CifString::GetCategoryFromCifItem(categoryName, _tBufKeyword);
  if (categoryName.empty())
  {
    log << "Warning -  Error in category name at line " << NDBlineNo <<
      " value " << _tBufKeyword << endl;
    errorLog += "Warning - Error in category name at line ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += " value ";
    errorLog += _tBufKeyword;
    errorLog += '\n';
    _warn++;
    return(0);
  }

  if (_readDef.AreAllCatsRead() == true)
  {
    return(STOP_PARSING);
  }

  if (!(_readDef.Category_OK(categoryName) &&
    _readDef.Datablock_OK(_curDataBlockName)))
  {
    return(0);
  }
  else
  {
    _readDef.IncreaseNumReadCats();
  }

  // If category aready exists, log the warning
  bool tablePresent = false;

  if (_fobj->IsBlockPresent(_curDataBlockName))
  {
      Block& block = _fobj->GetBlock(_curDataBlockName);
      tablePresent = block.IsTablePresent(categoryName);
  }

  if (tablePresent)
  {
    log << "Warning - Duplicate category name " << categoryName << 
      " at line " << NDBlineNo << endl;
    errorLog += "Warning - Duplicate category name ";
    errorLog += categoryName;
    errorLog += " at line ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += '\n';
    _warn++;
    return(0);
  }

  _curCategoryName = categoryName;

  _tbl = new ISTable(categoryName, _fobj->GetCaseSensitivity());

  _curRow = 0;

  ProcessItemNameList();

  return(0);

} // End of  CifParser::ProcessLoopDeclaration()


int CifParser::ProcessItemNameList(void)
/* ----------------------------------------------------------------------
   Purpose: CifParser::ProcessItemNameList(void)

            Registers the item keyword for the the current item in the 
	    current category.  Maintains an index array of "valid" keyword 
	    names in fieldList[].  This array is used to indirectly 
	    reference between keywords and values ...  
 * ----------------------------------------------------------------------*/
{

  string keywordName;
  string categoryName;

#if DEBUG
  if (_verbose) log << "Processing item name list at line " <<  NDBlineNo << " keyword " <<  _tBufKeyword << endl;
#endif

  // If no table exists, it means that item name list of a duplicate
  // category are being processed, and they should be ignored. Just return.
  if (_tbl == nullptr)
  {
    return(0);
  }

  CifString::GetCategoryFromCifItem(categoryName, _tBufKeyword);
  if (categoryName.empty())
  {
    log << "Warning - Error in category name at line " << NDBlineNo <<
      " value " << _tBufKeyword << endl;
    errorLog += "Warning - Error in category name at line ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += " value ";
    errorLog += _tBufKeyword;
    errorLog += '\n';
    _warn++;
#ifdef VLAD_FIX
    return(0);
#endif
  }

  if (_afterLoop == false)
  {
    if (_readDef.AreAllCatsRead() == true)
    {
      return(STOP_PARSING);
    }

    if (!(_readDef.Category_OK(categoryName) &&
      _readDef.Datablock_OK(_curDataBlockName)))
    {
      return(0);
    }
    else
    {
      _readDef.IncreaseNumReadCats();
    }
  }

  if (_curItemNo > _fieldListAlloc - 1)
  {
    _fieldListAlloc = _curItemNo + _fieldListAlloc; 
    _fieldList.reserve(_fieldListAlloc);
  }

  if (!String::IsEqual(categoryName, _curCategoryName,
    _fobj->GetCaseSensitivity()))
  {
    // If parsed category is different than current category in the loop,
    // just log and return.
    log << "Warning - Wrong category at line " << NDBlineNo <<
      " at item " << _tBufKeyword << endl;
    errorLog += "Warning - Wrong category at line ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += " item ";
    errorLog += _tBufKeyword;
    errorLog += '\n';
    _fieldList.clear();
    _warn++;
    _curItemNo = 0;
    return(0);
  }

  CifString::GetItemFromCifItem(keywordName, _tBufKeyword);
  if (keywordName.empty())
  {
    log << "Warning - Syntax error at line " << NDBlineNo <<
      " at item " << _tBufKeyword << endl;
    errorLog += "Warning - Syntax error at line ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += " item ";
    errorLog += _tBufKeyword;
    errorLog += '\n';
    _fieldList[_curItemNo].clear();
    _warn++;
    _curItemNo++;
    return(0);
  }

  if (!_tbl->IsColumnPresent(keywordName))
  {
#if DEBUG
    if (_verbose) log << "Line " << NDBlineNo << " keyword is " << keywordName << endl;
#endif
    // Add the new item
    _tbl->AddColumn(keywordName);
    if (_curItemNo >= (int)_fieldList.size())
      _fieldList.push_back(keywordName);
    else
      _fieldList[_curItemNo] = keywordName;
  } // Item not in the category
  else
  {
    log << "Warning - Duplicate item name " << _tBufKeyword  <<
      " at line " << NDBlineNo << endl;
    errorLog += "Warning - Duplicate item name ";
    errorLog += _tBufKeyword;
    errorLog += " at line ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += '\n';
    _warn++;
    return(0);
  } // Item found in the category

  _curItemNo++;

  return(0);

} // End of CifParser::ProcessItemNameList()


int CifParser::ProcessValueList(void)
/* ----------------------------------------------------------------------
     Purpose:  CifParser::ProcessValueList(void)

               Add the current value to the appropriate column in the 
               the current row.  Start a new row if necessary.
 * ----------------------------------------------------------------------*/
{

  // If no table exists, it means that value list of a duplicate
  // category is being processed, and it should be ignored. Just return.
  if (_tbl == nullptr)
  {
    return(0);
  }

  if (_afterLoop == false)
  {
    if (_readDef.AreAllCatsRead() == true)
    {
      return(STOP_PARSING);
    }

    if (!(_readDef.Category_OK(_curCategoryName) &&
      _readDef.Datablock_OK(_curDataBlockName)))
    {
      return(0);
    }
    else
    {
      _readDef.IncreaseNumReadCats();
    }
  }

#if DEBUG
  if (_verbose)
  {
    if (!_pBufValue.empty())
      log << "Processing value at line " << NDBlineNo << " value [" << _pBufValue << "]" << endl;
  }
#endif

  if (_fieldList.empty())
  {
      log << "ERROR: Parsing error at line " << NDBlineNo << " value [" << _pBufValue << "]" << endl;

      errorLog += "ERROR: Parsing error at line ";
      errorLog += String::IntToString(NDBlineNo);
      errorLog += " value [";
      errorLog += _pBufValue;
      errorLog += "]";
      errorLog += '\n';

      return(0);
  }

  if (!_fieldList[_curValueNo].empty())
  {
    // If it is a value for a valid item, process it
    //if (_tbl->GetNumColumns() == 1)
    if (_curValueNo == 0)
    {
      // If this is the very first value for item list, create a new row
      vector<string> rowBuf(_tbl->GetNumColumns(), CifString::UnknownValue);
      _curRow++;
      _tbl->AddRow(rowBuf);
    }

    if (!_pBufValue.empty())
    {
#ifdef PDB_TRACE
        cout << "TRACE: value in list \"" << _pBufValue << "\"" << endl;
#endif
	_tbl->UpdateCell(_curRow - 1, _fieldList[_curValueNo], _pBufValue);
    }
    else
    {
#ifdef PDB_TRACE
        cout << "DEBUG: value in list \"" << CifString::UnknownValue << "\"" << endl;
#endif
        _tbl->UpdateCell(_curRow - 1, _fieldList[_curValueNo], CifString::UnknownValue);
    }
  }

  _curValueNo++;

  if (_curValueNo == _curItemNo)
  {
    // If the number of values reached the number of items
#if DEBUG
    if (_verbose)
    {
      log << "Loading row " << _curRow -1 << " with " <<  rowBuf.size() << " elements" << endl;
      for (int i=0; i < rowBuf.size(); i++)
      {
        log << "Column [" << i << "] value "<< rowBuf[i] << endl;
      }
    }
#endif

    _curValueNo = 0;
  }

  return(0);

} // End of CifParser::ProcessValueList()


int CifParser::ProcessItemValuePair(void)
/* ----------------------------------------------------------------------
      Purpose: CifParser::ProcessItemValuePair()

               Assign the current value to its associated item name.
 * ----------------------------------------------------------------------*/
{

  string categoryName;
  string keywordName;

  _curItemNo  = 1;
  _curValueNo = 0;

#if DEBUG
  if (_verbose)
  {
    if (!_pBufValue.empty()) 
      log << "Processing " << _tBufKeyword << " at " <<  NDBlineNo << " value " << _pBufValue  << endl;

  }
#endif

  CifString::GetCategoryFromCifItem(categoryName, _tBufKeyword);
  if (categoryName.empty())
  {
    log << "Warning - Error in category name at line " << NDBlineNo <<
      " value " << _tBufKeyword << endl;
    errorLog += "Warning - Error in category name at line ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += " value ";
    errorLog += _tBufKeyword;
    errorLog += '\n';
    _warn++;
    return(0);
  }

  CifString::GetItemFromCifItem(keywordName, _tBufKeyword);
  if (keywordName.empty())
  {
    log << "Warning - Syntax error line at " << NDBlineNo << " item " << _tBufKeyword << endl;
    errorLog += "Warning - Syntax error at line ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += " item ";
    errorLog += _tBufKeyword;
    errorLog += '\n';
    _warn++;
    return(0);
  }

  if (!String::IsEqual(categoryName, _curCategoryName,
    _fobj->GetCaseSensitivity()))
  {
    /* If the new category is different than the current category. */
    if (_readDef.AreAllCatsRead() == true)
    {
      return(STOP_PARSING);
    }

    if (!(_readDef.Category_OK(categoryName) &&
      _readDef.Datablock_OK(_curDataBlockName)))
    {
      return(0);
    }
    else
    {
      _readDef.IncreaseNumReadCats();
    }
  }

  int tableWritten = 0;

  bool tablePresent = false;

  if ( _fobj->IsBlockPresent(_curDataBlockName))
  {
      Block& block = _fobj->GetBlock(_curDataBlockName);
      tablePresent = block.IsTablePresent(categoryName);
  }

  if (tablePresent)
  { //  duplicates a persistent table?
    log << "Warning - Duplicate category name " << categoryName << " at line " << NDBlineNo << endl;
    errorLog += "Warning - Duplicate category name ";
    errorLog += categoryName;
    errorLog += " at line ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += '\n';
    _warn++;

    tableWritten = 1;
  }

  if (String::IsEqual(categoryName, _curCategoryName,
    _fobj->GetCaseSensitivity()) && _afterLoop)
  { //  duplicates a persistent table-immediately after 
    // loop goes pair item-value for same category
    log << "Warning - Duplicate category name " << categoryName << " at line " << NDBlineNo << endl;
    errorLog += "Warning - Duplicate category name ";
    errorLog += categoryName;
    errorLog += " at line ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += '\n';
    _warn++;
  }

  if (!String::IsEqual(categoryName,_curCategoryName,
    _fobj->GetCaseSensitivity()))
  {
    if (_tbl != nullptr)
    {
      _ComplexWriteTable();

      // VLAD: Should _tbl be destroyed here first, prior to assigning NULL?
      _tbl = nullptr;
      _curCategoryName.clear();
    }

    if (tableWritten == 1)
    {
        // If table that is just read was already written, get it.

        Block& block = _fobj->GetBlock(_curDataBlockName);

        _tbl = block.GetTablePtr(categoryName);
    }
    else
    {
        // If table that is just read was not already written. Create a new
        // one.
        _tbl = new ISTable(categoryName, _fobj->GetCaseSensitivity());
    }

    _curCategoryName = categoryName;
  }

  _afterLoop=false;

  if (_tbl->IsColumnPresent(keywordName))
  {
    log << "Warning - Duplicate item name " <<_tBufKeyword  << " at line " << NDBlineNo << endl;
    errorLog += "Warning - Duplicate item name ";
    errorLog += _tBufKeyword;
    errorLog += " at line ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += '\n';
    _warn++;
  
    unsigned int numRows = _tbl->GetNumRows();
    if (!_pBufValue.empty())
    {
        for (unsigned int i = 0; i < numRows; i++)
        {
#ifdef PDB_TRACE
            cout << "DEBUG: value in pair \"" << _pBufValue << "\"" << endl;
#endif
            _tbl->UpdateCell(i, keywordName, _pBufValue);
        }
    }
    else
    {
        for (unsigned int i = 0; i < numRows; i++)
        {
#ifdef PDB_TRACE
            cout << "DEBUG: value in pair \"" << CifString::UnknownValue << "\"" << endl;
#endif
            _tbl->UpdateCell(i, keywordName, CifString::UnknownValue);
        }
    }
  }
  else
  {

#if DEBUG
  if (_verbose) log << "Line " << NDBlineNo << " keyword is " << _tBufKeyword << endl;
#endif

    _tbl->AddColumn(keywordName);

    if (_tbl->GetNumRows() == 0)
        _tbl->AddRow();

    if (!_pBufValue.empty())
    {
        for (unsigned int i = 0; i < _tbl->GetNumRows(); i++)
        {
#ifdef PDB_TRACE
          cout << "DEBUG: value in pair \"" << _pBufValue << "\"" << endl;
#endif
          _tbl->UpdateCell(i, keywordName, _pBufValue);
        }
    }
    else
    {
        for (unsigned int i = 0; i < _tbl->GetNumRows(); i++)
        {
#ifdef PDB_TRACE
            cout << "DEBUG: value in pair \"" << CifString::UnknownValue << "\"" << endl;
#endif
            _tbl->UpdateCell(i, keywordName, CifString::UnknownValue);
        }
    }
  }

  return(0);

} // End of CifParser::ProcessItemValuePair()


void CifParser::ProcessAssignments()
{
  if (_curValueNo!=0) {
  log <<"ERROR - Number of data values is not exact multiples of the number of data names   (look above line "<<NDBlineNo<<")"<<endl;
  errorLog += "ERROR - Number of data values is not exact multiples of the number of data names   (look above line ";
  errorLog += String::IntToString(NDBlineNo);
  errorLog += ")";
  errorLog += '\n';
  _warn++;
  }
}

void CifParser::ProcessLoop()
{
    _curItemNo = 0;
    _curValueNo = 0;
}

void CifParser::ProcessItemName()
{
    _tBufKeyword = Glob_tBufKeyword;
}

void CifParser::ProcessItemValue()
{
    _pBufValue = Glob_pBufValue;
}

void CifParser::ProcessLsItemValue()
{
    _pBufValue = Glob_pBufValue;
}

void CifParser::ProcessUnknownValue()
{
   _pBufValue = CifString::InapplicableValue;
}

void CifParser::ProcessMissingValue()
{
   _pBufValue = CifString::UnknownValue;
}

void CifParser::ProcessDataBlockName()
{

  string newDataBlock;


  // Write category from previous datablock
  if ((_tbl != nullptr) && (_curItemNo > 0) && (!_curDataBlockName.empty()) &&
    (!_curCategoryName.empty()))
  {
#if DEBUG    
    if (_verbose) log << " Save category " << _curCategoryName << " in " << _curDataBlockName << endl;
#endif
    _ComplexWriteTable();

    // VLAD: Should _tbl be destroyed here first, prior to assigning NULL?
    _tbl = nullptr;

    _curCategoryName.clear();
  }

  _prevDataBlockName = _curDataBlockName;

  _numDataBlocks++;
  _nTablesInBlock = 0;

  // Set current datablock name
  if ((Glob_dataBlockName != nullptr) &&
    (strlen(Glob_dataBlockName) > strlen(DATA_TAG)))
  {
    _curDataBlockName = &(Glob_dataBlockName)[strlen(DATA_TAG)];

    newDataBlock = _fobj->AddBlock(_curDataBlockName);

    if (newDataBlock != _curDataBlockName)
    {
      // Issue a warning on duplicate data block
      log << "Warning - Duplicate datablock name " << _curDataBlockName << 
        " at line " << NDBlineNo << " replaced with " << newDataBlock <<
        endl;
      errorLog += "Warning - Duplicate datablock name ";
      errorLog += _curDataBlockName;
      errorLog += " at line ";
      errorLog += String::IntToString(NDBlineNo);
      errorLog += " replaced with ";
      errorLog += newDataBlock;
      errorLog += '\n';
      _warn++;
    }
  }
  else
  {
    newDataBlock = _fobj->AddBlock("");

    // Issue a warning on empty data block
    log << "Warning - Empty datablock name " <<  "at line " <<
      NDBlineNo << " replaced with " << newDataBlock << endl;
    errorLog += "Warning - Empty datablock name ";
    errorLog += " at line ";
    errorLog += String::IntToString(NDBlineNo);
    errorLog += " replaced with ";
    errorLog += newDataBlock;
    errorLog += '\n';
    _warn++;
  }

  _curDataBlockName = newDataBlock;

  _curCategoryName.clear();

#if DEBUG
  if (_verbose) log << " Previous data block is now  " <<  _prevDataBlockName << endl;
  if (_verbose) log << " Current  data block is now  " <<  _curDataBlockName << endl;
#endif

} // End of CifParser::ProcessDataBlockName()

CifParser::~CifParser()
{
    Reset();
    CifParserP = nullptr;
}


void CifParser::_ComplexWriteTable()
{

    if (_curDataBlockName.empty())
    {
        string newDataBlock = _fobj->AddBlock(_curDataBlockName);

        _curDataBlockName = newDataBlock;
    }

    // If current block is empty add it.
    Block& block = _fobj->GetBlock(_curDataBlockName);

    block.WriteTable(_tbl);

}


void ProcessAssignmentsFromParser()
{
    CifParserP->ProcessAssignments();
}

void ProcessLoopFromParser()
{
    CifParserP->ProcessLoop();
}

int ProcessItemValuePairFromParser()
{
    return(CifParserP->ProcessItemValuePair());
}

int ProcessLoopDeclarationFromParser()
{
    return(CifParserP->ProcessLoopDeclaration());
}

int ProcessItemNameListFromParser()
{
    return(CifParserP->ProcessItemNameList());
}

int ProcessValueListFromParser()
{
    return(CifParserP->ProcessValueList());
}

void ProcessItemNameFromParser()
{
    CifParserP->ProcessItemName();
}

void ProcessItemValueFromParser()
{
    CifParserP->ProcessItemValue();
}

void ProcessLsItemValueFromParser()
{
    CifParserP->ProcessLsItemValue();
}

void ProcessUnknownValueFromParser()
{
    CifParserP->ProcessUnknownValue();
}

void ProcessMissingValueFromParser()
{
    CifParserP->ProcessMissingValue();
}

void ProcessDataBlockNameFromParser()
{
    CifParserP->ProcessDataBlockName();
}

