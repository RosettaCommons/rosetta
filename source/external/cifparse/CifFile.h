/*
FILE:     CifFile.h
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
** \file CifFile.h
**
** \brief Header file for CifFile class.
*/


/* 
  PURPOSE:    Base class for read/write cif files
*/


#ifndef CIFFILE_H
#define CIFFILE_H


#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <sstream>

#include "GenString.h"
#include "CifString.h"
#include "TableFile.h"
#include "CifParentChild.h"



/**
**  \class CifFile
**
**  \brief Public class that represents a CIF file, composed of blocks with
**    tables.
**
**  This class represents a CIF file. In addition to inherited methods from
**  \e TableFile class, this class provides methods for writing the data to
**  a text file, along with methods for controlling how the data is written,
**  and a method for verifying the CIF file against dictionary.
*/
class CifFile : public TableFile
{
  public:
    std::string _parsingDiags;
    std::string _checkingDiags;

    static const unsigned int STD_CIF_LINE_LENGTH = 80;

    enum eQuoting
    {
        eSINGLE = 0,
        eDOUBLE
    };

    /**
    **  Constructs a CIF file.
    **
    **  \param[in] fileMode - CIF file mode. Possible values are
    **    read-only, create, update and virtual. Detailed description of
    **    file mode is given in \e TableFile documentation.
    **  \param[in] fileName - relative or absolute name of the file
    **    where object persistency is maintained. If \e fileMode specifies
    **    virtual mode, this parameter is ignored.
    **  \param[in] verbose - optional parameter that indicates whether
    **    logging should be turned on (if true) or off (if false).
    **    If \e verbose is not specified, logging is turned off.
    **  \param[in] caseSense - optional parameter that indicates case
    **    sensitivity of table names. Possible values are case sensitive
    **    and case in-sensitive. If not specified, case sensitive table
    **    names are assumed.
    **  \param[in] maxLineLength - optional parameter that indicates the
    **    maximum number of written characters in one line in the written
    **    text file. If not specified, \e STD_CIF_LINE_LENGTH is used.
    **  \param[in] nullValue - optional parameter that indicates the
    **    character that is to be used to denote unknown value in the written
    **    CIF file. If not specified, \e CifString::UnknownValue is used.
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    CifFile(const eFileMode fileMode, const std::string& fileName,
      const bool verbose = false, const Char::eCompareType
      caseSense = Char::eCASE_SENSITIVE,
      const unsigned int maxLineLength = STD_CIF_LINE_LENGTH,
      const std::string& nullValue = CifString::UnknownValue);

    /**
    **  Constructs a CIF file in virtual mode.
    **
    **  \param[in] verbose - optional parameter that indicates whether
    **    logging should be turned on (if true) or off (if false).
    **    If \e verbose is not specified, logging is turned off.
    **  \param[in] caseSense - optional parameter that indicates case
    **    sensitivity of table names. Possible values are case sensitive
    **    and case in-sensitive. If not specified, case sensitive table
    **    names are assumed.
    **  \param[in] maxLineLength - optional parameter that indicates the
    **    maximum number of written characters in one line in the written
    **    text file. If not specified, \e STD_CIF_LINE_LENGTH is used.
    **  \param[in] nullValue - optional parameter that indicates the
    **    character that is to be used to denote unknown value in the written
    **    CIF file. If not specified, \e CifString::UnknownValue is used.
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    CifFile(const bool verbose = false, const Char::eCompareType
      caseSense = Char::eCASE_SENSITIVE,
      const unsigned int maxLineLength = STD_CIF_LINE_LENGTH,
      const std::string& nullValue = CifString::UnknownValue);

    /**
    **  Destructs a CIF file, by releasing all consumed resources.
    **
    **  \param: Not applicable
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    ~CifFile();

    /**
    **  Sets file name of a file that was the source of the object data.
    **
    **  \param srcFileName - The name of the source data file.
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void SetSrcFileName(const std::string& srcFileName);


    /**
    **  Retrieves source file name.
    **
    **  \param: None
    **
    **  \return - source file name
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    const std::string& GetSrcFileName();


    /**
    **  Retrieves logging option.
    **
    **  \param: None
    **
    **  \return true - if logging is turned on
    **  \return false - if logging is turned off
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline bool GetVerbose();

    /**
    **  Sets smart printing option. Smart printing is used to beautify the
    **  output of a written text file.
    **
    **  \param smartPrint - smart printing. If false, smart printing is
    **    disabled. If true, smart printing is enabled. If not specified,
    **    smart printing is enabled.
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline void SetSmartPrint(bool smartPrint = true);


    /**
    **  Retrieves smart printing option.
    **
    **  \param: None
    **
    **  \return true - if smart printing is enabled
    **  \return false - if smart printing is disabled
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline bool IsSmartPrint();

    /**
    **  Sets quoting option. This option is used in order to
    **  select the type of quoting to be used in the written text file.
    **
    **  \param quoting - type of quoting. If \e eSINGLE, single quotes are
    **    used. If \e eDOUBLE, double quotes are used.
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void SetQuoting(eQuoting quoting);

    /**
    **  Retrieves quoting option.
    **
    **  \param: None
    **
    **  \return eSINGLE - if single quotes are used
    **  \return eDOUBLE - if double quotes are used
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    unsigned int GetQuoting();
 
    /**
    **  This method is used in order to control how single row categories are
    **  written: in form of a "loop_" construct or as an item-value pair.
    **
    **  \param catName - category name
    **  \param looping - category looping option. If false and the
    **    category is a single row category, that category will not be
    **    written with "loop_" construct. Otherwise, if true, single row
    **    category will be written with "loop_" construct.
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void SetLooping(const std::string& catName, bool looping = false);

    /**
    **  Retrieves looping option of a category.
    **
    **  \param catName - category name
    **
    **  \return - category looping option, as described in SetLooping() method.
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    bool GetLooping(const std::string& catName);
 
    /**
    **  Writes the data out to a text file.
    **
    **  \param[in] cifFileName - relative or absolute name of the text file
    **    to which the data from \e CifFile object is to be written to.
    **  \param[in] sortTables - optional parameter that indicates whether
    **    written tables should be sorted (if true) or not sorted (if false).
    **    If \e sortTables is not specified, tables are not sorted prior to
    **    writing them.
    **  \param[in] writeEmptyTables - optional parameter that indicates
    **    whether empty tables (0 rows) are to be written (if true) or not
    **    written (if false). If \e writeEmptyTables is not specified, empty
    **    tables are not written.
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void Write(const std::string& cifFileName, const bool sortTables = false,
      const bool writeEmptyTables = false);

    /**
    **  Writes the data out to a text file.
    **
    **  \param[in] cifFileName - relative or absolute name of the text file
    **    to which the data from \e CifFile object is to be written to.
    **  \param[in] tableOrder - vector of table names that indicates the
    **    order of written tables.
    **  \param[in] writeEmptyTables - optional parameter that indicates
    **    whether empty tables (0 rows) are to be written (if true) or not
    **    written (if false). If \e writeEmptyTables is not specified, empty
    **    tables are not written.
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void Write(const std::string& cifFileName,
      const std::vector<std::string>& tableOrder,
      const bool writeEmptyTables = false);

    /**
    **  Writes the data out to an output stream.
    **
    **  \param[in] outStream - a reference to the output stream
    **  \param[in] sortTables - optional parameter that indicates whether
    **    written tables should be sorted (if true) or not sorted (if false).
    **    If \e sortTables is not specified, tables are not sorted prior to
    **    writing them.
    **  \param[in] writeEmptyTables - optional parameter that indicates
    **    whether empty tables (0 rows) are to be written (if true) or not
    **    written (if false). If \e writeEmptyTables is not specified, empty
    **    tables are not written.
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void Write(std::ostream& outStream, const bool sortTables = false,
      const bool writeEmptyTables = false);

    /**
    **  Writes the data out to a text file in NMR-STAR format.
    **
    **  \param[in] nmrStarFileName - relative or absolute name of the text file
    **    to which the data from \e CifFile object is to be written to.
    **  \param[in] globalBlockName - the name of the global NMR-STAR block.
    **  \param[in] sortTables - optional parameter that indicates whether
    **    written tables should be sorted (if true) or not sorted (if false).
    **    If \e sortTables is not specified, tables are not sorted prior to
    **    writing them.
    **  \param[in] writeEmptyTables - optional parameter that indicates
    **    whether empty tables (0 rows) are to be written (if true) or not
    **    written (if false). If \e writeEmptyTables is not specified, empty
    **    tables are not written.
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void WriteNmrStar(const std::string& nmrStarFileName,
      const std::string& globalBlockName,  const bool sortTables = false,
      const bool writeEmptyTables = false);

    /**
    **  Checks a CIF file (all blocks in it) against the dictionary.
    **
    **  \param[in] dicRef - reference to a dictionary file. The check is
    **    done against the first block in the dictionary file.
    **  \param[in] diagFileName - relative or absolute name of the file,
    **    where diagnostic information is stored.
    **  \param[in] extraDictChecks - optional parameter that indicates whether
    **    to perform additional, non-standard, dictionary checks. If not
    **    specified, those checks are not performed.
    **  \param[in] extraCifChecks - optional parameter that indicates whether
    **    to perform additional, non-standard, CIF checks. If not specified,
    **    those checks are not performed.
    **
    **  \return 0 - if all checks passed
    **  \return different than 0 - if checks failed
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    int DataChecking(CifFile& dicRef, const std::string& diagFileName,
      const bool extraDictChecks = false, const bool extraCifChecks = false);

    /**
    **  Checks a block of CIF file against the specified reference block.
    **
    **  \param[in] block - reference to a block that is to be checked
    **  \param[in] refBlock - reference to a reference block against which
    **    \e block is to be checked
    **  \param[out] diagBuf - diagnostics buffer that holds checking results
    **  \param[in] extraDictChecks - optional parameter that indicates whether
    **    to perform additional, non-standard, checks. If not specified, those
    **    checks are not performed.
    **  \param[in] extraCifChecks - optional parameter that indicates whether
    **    to perform additional, non-standard, CIF checks. If not specified,
    **    those checks are not performed.
    **
    **  \return 0 - if all checks passed
    **  \return different than 0 - if checks failed
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    int DataChecking(Block& block, Block& refBlock, std::ostringstream& buf,
      const bool extraDictChecks = false, const bool extraCifChecks = false);

    /**
    **  Sets enumerations checking option for case-insensitive types.
    **
    **  \param caseSense - case sensitivity of enumeration values check. If
    **    false, enumeration values of case-insensitive types will be checked
    **    as case-insensitive. If true, enumeration values of case-insensitive
    **    types will be checked as case-sensitive.
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void SetEnumCheck(bool caseSense = false);

    /**
    **  Retrieves enumerations checking option for case-insensitive types.
    **
    **  \param: None
    **
    **  \return true - if case-sensitive enumeration check is enabled
    **  \return false - if case-insensitive enumeration check is enabled
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    bool GetEnumCheck();

    /**
    **  Gets parsing diagnostics.
    **
    **  \param: None
    **
    **  \return - reference to parsing diagnostics
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    const std::string& GetParsingDiags();


    /**
    **  Finds indices of rows that contain all CIF null values. A CIF null
    **  value is defined as a "?" or "".
    **
    **  \param[out] nullRowsIndices - vector of null rows indices.
    **  \param[in] isTable - table reference
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void FindCifNullRows(std::vector<unsigned int>& nullRowsIndices,
      const ISTable& isTable);

    void GetAttributeValue(std::string& attribVal, const std::string& blockId,
      const std::string& category, const std::string& attribute);
    void GetAttributeValueIf(std::string& attribVal, const std::string& blockId,
      const std::string& category, const std::string& attributeA,
      const std::string& attributeB, const std::string& valB);
    bool IsAttributeValueDefined(const std::string& blockId,
      const std::string& category, const std::string& attribute);

    void SetAttributeValue(const std::string& blockId,
      const std::string& category,
      const std::string& attribute, const std::string& value,
      const bool create = false);
    void SetAttributeValueIf(const std::string& blockId,
      const std::string& category, const std::string& attributeA,
      const std::string& valA,
      const std::string& attributeB, const std::string& valB,
       const bool create = false);
    void SetAttributeValueIfNull(const std::string& blockId,
      const std::string& category, const std::string& attribute,
      const std::string& value);

    void GetAttributeValues(std::vector<std::string>& strings,
      const std::string& blockId,
      const std::string& category, const std::string& attribute);
    void GetAttributeValuesIf(std::vector<std::string>& strings, 
      const std::string& blockId, const std::string& category,
      const std::string& attributeA, 
      const std::string& attributeB, const std::string& valB);

    void SetAttributeValues(const std::string& blockId,
      const std::string& category, const std::string& attribute,
      const std::vector<std::string>& values);

#ifdef VLAD_TO_CIF_FILE_NOT_USED
    void del_attribute_value_where(CifFile *fobj, const char *blockId,
      const char *category, const char *attributeB, const char *valB);
#endif // VLAD_TO_CIF_FILE_NOT_USED not defined

    int CheckCategories(Block& block, Block& refBlock, std::ostringstream& log);
    void CheckCategoryKey(Block& block, std::ostringstream& log);
    void CheckItemsTable(Block& block, std::ostringstream& log);
    int CheckItems(Block& block, Block& refBlock, std::ostringstream& log);


  protected:
    static const unsigned int STD_PRINT_SPACING = 3;
    static const unsigned int SMART_PRINT_SPACING = 1;
    static const unsigned int HEADER_SPACING = 40;

    enum eIdentType
    {
        eNONE = 0,
        eLEFT,
        eRIGHT
    };

    std::string _beginDataKeyword;
    std::string _endDataKeyword;

    std::string _beginLoopKeyword;
    std::string _endLoopKeyword;

    unsigned int _maxCifLineLength;
    std::string _nullValue;
    bool _verbose;
    bool _smartPrint;
    std::string _quotes;
    std::map<std::string, bool> _looping;
    bool _enumCaseSense;

    int _IsQuotableText(const std::string& itemValue);
    eIdentType _FindPrintType(const std::vector<std::string>& values);

    void _PrintItemIdent(std::ostream& cifo, unsigned int& linePos);
    void _PrintItemName(std::ostream& cifo, const std::string& category,
      const std::string& itemName, unsigned int& linePos);
    void _PrintPostItemSeparator(std::ostream& cifo, unsigned int& linePos,
      const bool ident = false, const unsigned int numSpaces = 1);

    int _PrintItemValue(std::ostream& cifo, const std::string& itemValue,
      unsigned int& linePos, const eIdentType identType = eNONE,
      const unsigned int width = 0);

    int _PrintItemNameInHeader(std::ostream& cifo, const std::string& itemValue,
      unsigned int& linePos, const eIdentType identType = eNONE,
      const unsigned int width = 0);

    void _PrintHeaderedItems(std::ostream& cifo,
      const std::vector<std::string>& colNames,
      const std::vector<unsigned int>& colWidths,
      const std::vector<eIdentType> colPrintType);

    void Write(std::ostream& cifo, const std::vector<std::string>& catOrder,
      const bool writeEmptyTables = false);

    void Write(std::ostream& cifo, std::vector<unsigned int>& tables,
      const bool writeEmptyTables = false);


  private:
    std::string _srcFileName;

    bool _extraDictChecks;
    bool _extraCifChecks;

    void Init();

    bool IsCatDefinedInRef(const std::string& catName, ISTable& catTable);
    bool IsItemDefinedInRef(const std::string& catName,
      const std::string& itemName, ISTable& refItemTable);
    void CheckKeyItems(const std::string& blockName, ISTable& catTable,
      ISTable& keyTable, std::ostringstream& log);
    void CheckKeyValues(const std::vector<std::string>& keyItems,
      ISTable& catTable, std::ostringstream& log);

    void GetKeyAttributes(std::vector<std::string>& keyAttributes,
      const std::string& catTableName, ISTable& catKeyTable);
    void CheckKeyItems(const std::string& blockName, ISTable& catTable,
      const std::vector<std::string>& keyAttributes, std::ostringstream& log);

    void CheckMandatoryItems(const std::string& blockName, ISTable& catTable,
      ISTable& refItemTable, const std::vector<std::string>& keyItems,
      std::ostringstream& log);

    void CheckAndRectifyItemTypeCode(Block& block, std::ostringstream& log);
    void RectifyItemTypeCode(std::string& retItemTypeCode,
      std::ostringstream& log, Block& block, CifParentChild& cifParentChild,
      const std::string& cifItemName);

    int CheckRegExpRangeEnum(Block& block, ISTable& catTable,
      const std::string& attribName, ISTable& itemTypeTable,
      ISTable& itemTypeListTable, ISTable* itemRangeTableP,
      ISTable* itemEnumTableP, ISTable& parChildTable, ISTable* itemAliasesP,
      std::ostringstream& log);

    int CheckCellRange(const std::string& cell, const std::string& typeCode,
      const std::vector<std::string>& minlist,
      const std::vector<std::string>& maxlist);

    int CheckCellEnum(const std::string& cell, const std::string& typeCode,
      const std::string& primCode, const std::vector<std::string>& enumlist);

    int CheckCellFloatRange(const std::string& cell,
      const std::vector<std::string>& minlist,
      const std::vector<std::string>& maxlist);

    int CheckCellIntRange(const std::string& cell,
      const std::vector<std::string>& minlist,
      const std::vector<std::string>& maxlist);

    int CheckCellFloatEnum(const std::string& cell,
      const std::vector<std::string>& enumlist);

    int CheckCellIntEnum(const std::string& cell,
      const std::vector<std::string>& enumlist);

    int CheckCellOtherEnum(const std::string& cell, const std::string& primCode,
      const std::vector<std::string>& enumlist);

    void GetItemTypeCode(std::string& typeCode, const std::string& cifItemName,
      ISTable& itemTypeTable);

    void ConvertEscapedString(const std::string& inString,
      std::string& outString);
};


inline bool CifFile::GetVerbose()
{
    return(_verbose);
}


inline void CifFile::SetSmartPrint(bool smartPrint)
{
    _smartPrint = smartPrint;
}


inline bool CifFile::IsSmartPrint()
{
    return(_smartPrint);
}


#endif
