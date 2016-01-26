/*
FILE:     TableFile.h
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
** \file TableFile.h
**
** \brief Header file for Block and TableFile classes.
*/


#ifndef TABLEFILE_H
#define TABLEFILE_H


#include <vector>
#include <set>

#include "rcsb_types.h"
#include "mapped_ptr_vector.h"
#include "mapped_ptr_vector.C"
#include "GenString.h"
#include "ISTable.h"
#include "Serializer.h"


/**
**  \class Block
**
**  \brief Public class that represents a data block, that contains tables.
**
**  This class represents a data block, that can come from DDL,
**  dictionary or CIF files. Data block is a container of tables.
**  This class provides methods for construction and destruction, tables
**  manipulation (addition, retrieval, deleting, writing), data blocks
**  comparison.
*/
class Block
{
  public:
    mapped_ptr_vector<ISTable, StringLess> _tables;

    /**
    **  Utility method, not part of users public API, and will soon be removed.
    **
    **  Constructs a data block.
    **
    **  \param[in] name - the name of the data block
    **  \param[in] serP - pointer to the File Navigator object
    **  \param[in] fileMode - optional parameter that indicates data block
    **    mode. Possible values are read-only, create, update and virtual.
    **  \param[in] caseSense - optional parameter that indicates case
    **    sensitivity of table names. Possible values are case sensitive
    **    and case in-sensitive. If not specified, case sensitive table
    **    names are assumed.
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    Block(const string& name, Serializer* serP,
      const eFileMode fileMode = READ_MODE, const Char::eCompareType
      caseSense = Char::eCASE_SENSITIVE);

    /**
    **  Utility method, not part of users public API, and will soon be removed.
    **
    **  Destructs a data block.
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
    ~Block();

    /**
    **  Compares a data block to another data block.
    **
    **  \param[in] inBlock - reference to input data block
    **
    **  \return vector of pairs, where the first value in a pair is a
    **    table name and the second value in a pair is one of the following
    **    indicators of table differences: \n \n
    **  eMISSING - table exists only in the input block and not in this
    **    block \n
    **  eEXTRA - table exists only in this block and not in the input block \n
    **  eCASE_SENSE - table exists in both blocks, but with different column
    **    name case sensitivity \n
    **  eMORE_COLS - table exists in both blocks, but the table in this block
    **    has more columns than the table in the input block \n
    **  eLESS_COLS - table exists in both blocks, but the table in this block
    **    has less columns than the table in the input block \n
    **  eCOL_NAMES - table exists in both blocks, but tables have different
    **    column names \n
    **  eMORE_ROWS - table exists in both blocks, but the table in this block
    **    has more rows than the table in the input block \n
    **  eLESS_ROWS - table exists in both blocks, but the table in this block
    **    has less rows than the table in the input block \n
    **  eCELLS - table exists in both blocks, but tables have different
    **    content \n
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    vector<pair<string, ISTable::eTableDiff> > operator==(Block& inBlock);

    /**
    **  Utility method, not part of users public API, and will soon be removed.
    **
    **  Sets the name of a data block.
    **
    **  \param[in] name - the name of the data block
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline void SetName(const string& name);

    /**
    **  Retrieves data block name.
    **
    **  \param: None
    **
    **  \return Constant reference to a string that contains data block name.
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline const string& GetName() const;

    /**
    **  Adds a table to the block. If a table with the specified name
    **  already exists, it will be overwritten.
    **
    **  \param[in] name - optional parameter that indicates the name of the
    **    table to be added
    **  \param[in] colCaseSense - optional parameter that indicates case
    **    sensitivity of column names. Possible values are case sensitive and
    **    case in-sensitive. If not specified, a table with case sensitive
    **    column names is constructed.
    **
    **  \return Reference to the table
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    ISTable& AddTable(const std::string& name = string(),
      const Char::eCompareType colCaseSense = Char::eCASE_SENSITIVE);

    /**
    **  Utility method, not part of users public API, and will soon be removed.
    */
    void AddTable(const string& name, const int indexInFile = 0,
      ISTable* isTableP = NULL);

    /**
    **  Changes the name of a table in the data block.
    **
    **  \param[in] oldName - the name of the table which is to be renamed
    **  \param[in] newName - the new table name
    **
    **  \return None
    **
    **  \pre \e oldName must be non-empty
    **  \pre Table with name \e oldName must be present
    **  \pre \e newName must be non-empty
    **  \pre Table with name \e newName must not be present
    **  \pre Block must be in create or update mode
    **
    **
    **  \post None
    **
    **  \exception EmptyValueException - if \e oldName is empty
    **  \exception NotFoundException - if table with name \e oldName does
    **    not exist
    **  \exception EmptyValueException - if \e newName is empty
    **  \exception AlreadyExistsException - if table with name \e newName
    **    already exists
    **  \exception FileModeException - if block is not in create or
    **    update mode
    */
    void RenameTable(const string& oldName, const string& newName);

    /**
    **  Retrieves names of all tables in a data block.
    **
    **  \param[out] tableNames - retrieved table names
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void GetTableNames(vector<string>& tableNames);

    /**
    **  Checks for table presence in the data block.
    **
    **  \param[in] tableName - table name
    **
    **  \return true - if table exists
    **  \return false - if table does not exist
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    bool IsTablePresent(const string& tableName);

    /**
    **  Retrieves a table reference.
    **
    **  \param[in] tableName - table name
    **
    **  \return Reference to the table, if table was found
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception NotFoundException - if table with name \e tableName
    **    does not exist
    */
    ISTable& GetTable(const string& tableName);

    /**
    **  Retrieves a pointer to the table.
    **
    **  \param[in] tableName - table name
    **
    **  \return Pointer to the table, if table was found
    **  \return NULL, if table was not found
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    ISTable* GetTablePtr(const string& tableName);

    /**
    **  Deletes a table from a data block.
    **
    **  \param[in] tableName - table name
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void DeleteTable(const string& tableName);

    /**
    **  Writes a table to the data block. In this context, writing means
    **  adding it (if it does not already exist) or updating it (if it
    **  already exists).
    **
    **  \param[in] isTable - reference to the table
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void WriteTable(ISTable& isTable);

    /**
    **  Writes a table to the data block. In this context, writing means
    **  adding it (if it does not already exist) or updating it (if it
    **  already exists).
    **
    **  \param[in] isTableP - pointer to the table
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void WriteTable(ISTable* isTableP);

    /**
    **  Utility method, not part of users public API, and will soon be removed.
    */
    void Print();

  private:
    string _name;
    eFileMode _fileMode;
    Serializer* _ser;

    Block(const Block& t);
    Block& operator=(const Block& inBlock);

    ISTable* _GetTablePtr(const unsigned int tableIndex);
};


/**
**  \class TableFile
**
**  \brief Public class that represents a file composed of blocks with tables.
**
**  This class represents an ordered container of data blocks. Data blocks can
**  come from DDL, dictionary or CIF files, where each data block is a
**  container of tables. This class provides methods for construction and
**  destruction, data blocks manipulation (addition, retrieval, renaming.).
**  The class does in-memory management of data blocks, as well as
**  serialization and de-serialization to and from a file. The class supports
**  the following file modes: read-only, create, update and virtual. In
**  read-only mode, blocks and tables can only be read (from an existing table
**  file that has been previously serialized to a file) and cannot be
**  modified. Create mode is used to create a table file from scratch and
**  add/update blocks and tables in it and serialize it to a file. Update mode
**  is used to update an existing table file (that has been previously
**  serialized to a file). Virtual mode only provides in-memory management of
**  data blocks and is used when object persistency is not needed. Hence, all
**  modes except virtual mode provide association between in-memory data
**  blocks and persistent data blocks stored in a file.
*/
class TableFile
{
  public:
    enum eStatusInd
    {
        eCLEAR_STATUS = 0x0000,
        eDUPLICATE_BLOCKS = 0x0001,
        eUNNAMED_BLOCKS = 0x0002
    };

    /**
    **  Constructs a table file.
    **
    **  \param[in] caseSense - optional parameter that indicates case
    **    sensitivity of table names in blocks. Possible values are case
    **    sensitive and case in-sensitive. If not specified, case sensitive
    **    table names are assumed.
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post Constructed table file has virtual file mode.
    **
    **  \exception: None
    */
    TableFile(const Char::eCompareType caseSense = Char::eCASE_SENSITIVE);

    /**
    **  Constructs a table file.
    **
    **  \param[in] fileMode - table file mode. Possible values are
    **    read-only, create, update and virtual.
    **  \param[in] fileName - relative or absolute name of the file
    **    where object persistency is maintained. If \e fileMode specifies
    **    virtual mode, this parameter is ignored.
    **  \param[in] caseSense - optional parameter that indicates case
    **    sensitivity of table names in blocks. Possible values are case
    **    sensitive and case in-sensitive. If not specified, case sensitive
    **    table names are assumed.
    **
    **  \return Not applicable
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    TableFile(const eFileMode fileMode, const string& fileName,
      const Char::eCompareType caseSense = Char::eCASE_SENSITIVE);

    /**
    **  Destructs a table file, by first flushing all the modified tables in
    **  data blocks (for create mode or update mode) and then releasing all
    **  in-memory objects.
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
    virtual ~TableFile();

    /**
    **  Retrieves the name of the file that persistently holds data blocks
    **  and their tables.
    **
    **  \param: None
    **
    **  \return String that contains the file name.
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline string GetFileName(void);

    /**
    **  Retrieves table file mode.
    **
    **  \param: None
    **
    **  \return READ_MODE - if read-only mode
    **  \return CREATE_MODE - if create mode
    **  \return UPDATE_MODE - if update mode
    **  \return VIRTUAL_MODE - if virtual mode
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline eFileMode GetFileMode(void);

    /**
    **  Retrieves case sensitivity of table names in blocks.
    **
    **  \param: None
    **
    **  \return eCASE_SENSITIVE - if case sensitive
    **  \return eCASE_INSENSITIVE - if case in-sensitive
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline Char::eCompareType GetCaseSensitivity(void);

    /**
    **  Retrieves table file status in form of one or more flags.
    **
    **  \param: None
    **
    **  \return One or more of these flags: \n
    **    eCLEAR_STATUS - no flag value indicates that there are no flags set \n
    **    eDUPLICATE_BLOCKS - flag that indicates existence of blocks with
    **    the same name, which are internally stored with different names \n
    **    eUNNAMED_BLOCKS - flag that indicates existence of blocks with
    **    empty names
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline unsigned int GetStatusInd(void);

    /**
    **  Retrieves the number of data blocks in the table file.
    **
    **  \param: None
    **
    **  \return The number of data blocks in the table file.
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    inline unsigned int GetNumBlocks();

    /**
    **  Retrieves data block names.
    **
    **  \param[out] blockNames - retrieved data block names
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void GetBlockNames(vector<string>& blockNames);

    /**
    **  Retrieves the name of the first data block.
    **
    **  \param: None
    **
    **  \return String that contains the name of the first data block.
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    string GetFirstBlockName();

    /**
    **  Checks for data block existence.
    **
    **  \param[in] blockName - the name of the data block
    **
    **  \return true - if data block exists
    **  \return false - if data block does not exist
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    bool IsBlockPresent(const string& blockName);

    /**
    **  Adds a block to the table file. If a block with the specified name
    **  already exists, table file stores it under different internal name,
    **  that is obtained by appending a "#" symbol and the current block
    **  count. After writing blocks, with these kinds of block names,
    **  to an ASCII file, "#" symbol becomes a comment and the text after
    **  it is ignored. This enables the preservation of all duplicate blocks
    **  in the written file.
    **
    **  \param[in] blockName - the name of the data block
    **
    **  \return String that contains the internally assigned data block name.
    **    This value is different from \e blockName, if data block with
    **    the name \e blockName, already exists when this method is invoked. 
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    string AddBlock(const string& blockName);

    /**
    **  Retrieves a reference to the data block in the table file.
    **
    **  \param[in] blockName - the name of the data block
    **
    **  \return Reference to the data block in the table file.
    **
    **  \pre Data block with name \e blockName must be present
    **
    **  \post None
    **
    **  \exception NotFoundException - if data block with name \e blockName
    **    does not exist
    */
    Block& GetBlock(const string& blockName);

    /**
    **  Changes the data block name.
    **
    **  \param[in] oldBlockName - the name of the data block which is to
    **    be renamed
    **  \param[in] newBlockName - the new data block name
    **
    **  \return String that contains the internally assigned data block name.
    **    This value is different from \e newBlockName, if data block with
    **    the name \e newBlockName, already exists when this method is invoked. 
    **
    **  \pre Table file must have at least one data block.
    **  \pre Data block with name \e oldBlockName must be present
    **
    **  \post None
    **
    **  \exception EmptyContainerException - if table file has no data blocks
    **  \exception NotFoundException - if data block with name \e oldBlockName
    **    does not exist
    */
    string RenameBlock(const string& oldBlockName, const string& newBlockName);

    /**
    **  Changes the name of the first data block in table file.
    **
    **  \param[in] newBlockName - the new data block name
    **
    **  \return String that contains the internally assigned data block name.
    **    This value is different from \e newBlockName, if data block with
    **    the name \e newBlockName, already exists when this method is invoked. 
    **
    **  \pre Table file must have at least one data block.
    **
    **  \post None
    **
    **  \exception EmptyContainerException - if table file has no data blocks
    */
    inline string RenameFirstBlock(const string& newBlockName);

    /**
    **  Writes only the new or modified tables in data blocks to the
    **  associated persistent storage file (specified at table file
    **  construction time).
    **
    **  \param: None
    **
    **  \return None
    **
    **  \pre Table file must be in create or update mode
    **
    **  \post None
    **
    **  \exception FileModeException - if table file is not in create or
    **    update mode
    */
    void Flush();

    /**
    **  Writes all the data blocks and their tables in the specified file.
    **  The inteded purpose of this method is to write to a file different
    **  than the one specified at construction time.
    **
    **  \param[in] fileName - relative or absolute name of the file
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void Serialize(const string& fileName);

    /**
    **  Flushes the table file (if in create or update mode) and closes
    **  the associated persistent storage file.
    **
    **  \param: None
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void Close();

  protected:
    string _fileName;

    eFileMode _fileMode;

    // Indicates case sensitivity of identifiers
    Char::eCompareType _caseSense;

    // Indicates the current status of the object
    unsigned int _statusInd;  

    mapped_ptr_vector<Block, StringLess> _blocks;

    Serializer* _f;

    void _SetStatusInd(const string& blockName);

    void _AddBlock(const string& blockName, Serializer* serP);

    void _GetNumTablesInBlocks(vector<UInt32>& numTablesInBlocks);

    ISTable* _GetTablePtr(const unsigned int blockIndex,
      const unsigned int tableIndex);
    void _GetAllTables();

    unsigned int GetTotalNumTables();
    void GetTableNames(vector<string>& tableNames);

    void GetTablesIndices(vector<unsigned int>& tablesIndices);
    void GetSortedTablesIndices(vector<unsigned int>& tablesIndices);

    void _ReadFileIndex();
    void _ReadFileIndexVersion0();
    void _ReadFileIndexVersion1();
    void _WriteFileIndex(Serializer* serP,
      const vector<unsigned int>& tableLocs);

  private:
    static const string _version;
    void Init();
    void Open(const string& fileName, const eFileMode fileMode);
    unsigned int GetBlockIndexFromTableId(const string& tableId);
    string GetTableNameFromTableId(const string& tableId);
    string MakeInternalBlockName(const string& blockName,
      const unsigned int blockIndex);
    void PrintHeaderInfo();
};


inline void Block::SetName(const string& name)
{
    _name = name;
}


inline const string& Block::GetName() const
{
    return _name;
} 


inline string TableFile::GetFileName(void)
{
    return _fileName;
}


inline eFileMode TableFile::GetFileMode(void)
{
    return _fileMode;
}


inline Char::eCompareType TableFile::GetCaseSensitivity(void)
{
    return(_caseSense);
}


inline unsigned int TableFile::GetStatusInd(void)
{
    return _statusInd;
}

   
inline unsigned int TableFile::GetNumBlocks()
{
    return _blocks.size();
}


inline string TableFile::RenameFirstBlock(const string& newBlockName)
{
    return(RenameBlock(GetFirstBlockName(), newBlockName));
}


#endif // TABLEFILE_H
