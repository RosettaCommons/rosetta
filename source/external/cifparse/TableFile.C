/*
FILE:     TableFile.C
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
** \file TableFile.C
**
** \brief Implementation file for Block and TableFile classes. 
*/


#include <string>
#include <iostream>

#include "Exceptions.h"
#include "GenString.h"
#include "TableFile.h"


using std::string;
using std::vector;
using std::make_pair;
using std::cout;
using std::endl;


const string TableFile::_version("V1");


Block::Block(const string& name, Serializer* ser,
  const eFileMode fileMode, const Char::eCompareType caseSense) :
  _tables(StringLess(caseSense)), _fileMode(fileMode), _ser(ser)
{

    _name = name;

}


Block::Block(const Block& block)
{
 
    _name = block._name;
    _fileMode = block._fileMode;
    _ser = block._ser;
    _tables = block._tables;

}


Block::~Block()
{

    _name.clear();
    _ser = nullptr;
    _tables.clear();

}


vector<pair<string, ISTable::eTableDiff> > Block::operator==(Block& inBlock)
{

    vector<pair<string, ISTable::eTableDiff> > diff;

    vector<string> firstBlockTableNames, secondBlockTableNames;

    GetTableNames(firstBlockTableNames);
    inBlock.GetTableNames(secondBlockTableNames);

    for (auto & firstBlockTableName : firstBlockTableNames)
    {
        if (!inBlock.IsTablePresent(firstBlockTableName))
        {
            diff.emplace_back(firstBlockTableName,
              ISTable::eEXTRA);
        }
        else
        {
            ISTable* firstBlockTableP =
              GetTablePtr(firstBlockTableName);
            ISTable* secondBlockTableP =
              inBlock.GetTablePtr(firstBlockTableName);

            ISTable::eTableDiff tableDiff = ((*firstBlockTableP) ==
              (*secondBlockTableP));
            if (tableDiff != ISTable::eNONE)
            {
                diff.emplace_back(firstBlockTableName,
                  tableDiff);
            }
        }
    }

    for (auto & secondBlockTableName : secondBlockTableNames)
    {
        if (!IsTablePresent(secondBlockTableName))
        {
            diff.emplace_back(secondBlockTableName,
              ISTable::eMISSING);
        }
    }

    return(diff);

}


ISTable& Block::AddTable(const std::string& name,
  const Char::eCompareType colCaseSense)
{
    ISTable* isTableP = new ISTable(name, colCaseSense);

    WriteTable(isTableP);

    return (*isTableP);
}


void Block::AddTable(const string& name, const int indexInFile,
  ISTable* isTableP)
{

    if (name.empty())
    {
        throw EmptyValueException("Empty table name",
          "Block::AddTable");
    }

    _tables.push_back(name, indexInFile);

    if (isTableP != nullptr)
        _tables.set(isTableP);

}


void Block::RenameTable(const string& oldName, const string& newName)
{

    if (_fileMode == READ_MODE)
    {
        throw FileModeException("Rename table in read-only block",
          "Block::RenameTable");
    }

    if (oldName.empty())
    {
        // Empty old table name.
        throw EmptyValueException("Empty old table name",
          "Block::RenameTable");
    }

    if (newName.empty())
    {
        // Empty new table name.
        throw EmptyValueException("Empty new table name",
          "Block::RenameTable");
    }

    if (IsTablePresent(newName))
    {
        // Table with new name already exists
        throw AlreadyExistsException("Found table with new name",
          "Block::RenameTable");
    }

    if (!IsTablePresent(oldName))
    {
        // Table with old name does not exist
        throw NotFoundException("Table with old name not found",
          "Block::RenameTable");
    }

    // Indicate that the table has been modified
    ISTable* isTableP = GetTablePtr(oldName);
    isTableP->SetModified(true);

    // Change the name in the container
    _tables.rename(oldName, newName);

}


void Block::GetTableNames(vector<string>& names)
{

    names.clear();

    for (unsigned int tableI = 0; tableI < _tables.size(); ++tableI)
    {
        names.push_back(_tables.get_name(tableI));
    }

}


bool Block::IsTablePresent(const string& name)
{

    if (name.empty())
    {
        return(false);
    }

    unsigned int tableIndex = _tables.find(name);
    if (tableIndex != _tables.size())
    {
        // Found it.
        return(true);
    }
    else
    {
        // Not found.
        return(false);
    }

}


ISTable& Block::GetTable(const string& name)
{
    ISTable* isTableP = GetTablePtr(name);

    if (isTableP == nullptr)
    {
        throw NotFoundException("Table \"" + name + "\" not found.",
          "Block::GetTable");
    }

    return (*isTableP);
}


ISTable* Block::GetTablePtr(const string& name)
{

    // VLAD - IMPROVEMENT
    // SHOULD WE THROW EXCEPTION NotFoundException rather than returning NULL

    if (name.empty())
    {
        return(nullptr);
    }

    unsigned int tableIndex = _tables.find(name);
    if (tableIndex == _tables.size())
    {
        return(nullptr);
    }
 
    return(_GetTablePtr(tableIndex));

}


void Block::DeleteTable(const string& name)
{

    if (_fileMode == READ_MODE)
    {
        throw FileModeException("Delete table in read-only block",
          "Block::DeleteTable");
    }

    if (name.empty())
    {
        throw EmptyValueException("Empty table name",
          "Block::DeleteTable");
    }

    unsigned int tableIndex = _tables.find(name);
    if (tableIndex == _tables.size())
    {
        throw NotFoundException("Table not found",
          "Block::DeleteTable");
    }

    ISTable* isTableP = &(_tables[tableIndex]);

    if (isTableP != nullptr)
    {
        delete isTableP;
    }

    _tables.erase(name);

}


void Block::WriteTable(ISTable& isTable)
{
    WriteTable(&isTable);
}


void Block::WriteTable(ISTable* isTableP)
{

    if (_fileMode == READ_MODE)
    {
        throw FileModeException("Write table in read-only block",
          "Block::WriteTable");
    }

    if (isTableP == nullptr)
    {
        throw EmptyValueException("NULL ISTable pointer",
          "Block::WriteTable");
    }

    if ((isTableP->GetName()).empty())
    {
        throw EmptyValueException("Empty table name",
          "Block::WriteTable");
    }

    isTableP->SetModified(true);

    unsigned int tableIndex = _tables.find(isTableP->GetName());
    if (tableIndex != _tables.size())
    {
        // Found it
        ISTable* currIsTableP = &(_tables[tableIndex]);

        if (currIsTableP != isTableP)
        {
            delete currIsTableP;
            _tables.set(isTableP);  
        }
    }
    else 
    {
        // Not found. Add new table to the block
        AddTable(isTableP->GetName(), 0, isTableP);
    }

}


void Block::Print()
{

    cout << _name << endl;

}


ISTable* Block::_GetTablePtr(const unsigned int tableIndex)
{

    if (&(_tables[tableIndex]) != nullptr)
    {
        // Table is already in memory. Just return the pointer to it.
        return &(_tables[tableIndex]);
    }
    else
    {
        // For write or virtual mode, there are no tables on the disk.
        if ((_fileMode == CREATE_MODE) || (_fileMode == VIRTUAL_MODE))
            return(nullptr);

        string name = _tables.get_name(tableIndex);

        ISTable* isTableP = new ISTable(name);

        isTableP->SetSerializer(_ser);

        _tables.set(isTableP);

        _tables.read(name);

        return(isTableP);
    }

}


TableFile::TableFile(const Char::eCompareType caseSense) :
_caseSense(caseSense), _blocks(StringLess(caseSense))
{

    Init();

}


TableFile::TableFile(const eFileMode fileMode, const string& fileName,
  const Char::eCompareType caseSense) : _caseSense(caseSense),
  _blocks(StringLess(caseSense))
{

    Init();

    _fileName = fileName;
    _fileMode = fileMode;

    Open(fileName, fileMode);

}


TableFile::~TableFile()
{

    Close();

    _fileName.clear();

    for (unsigned int blockI = 0; blockI < _blocks.size(); ++blockI)
    {
        for (unsigned int tableI = 0; tableI < _blocks[blockI]._tables.size();
          ++tableI)
        {
            if (&(_blocks[blockI]._tables[tableI]) != nullptr)
            {
                delete &(_blocks[blockI]._tables[tableI]);
            }
        }
        delete &(_blocks[blockI]);
    }

    _blocks.clear();

}


string TableFile::AddBlock(const string& blockName)
{

    string internalBlockName = MakeInternalBlockName(blockName,
      _blocks.size() + 1);

    _AddBlock(internalBlockName, _f);

    return(internalBlockName);

} // End of TableFile::AddBlock()


bool TableFile::IsBlockPresent(const string& blockName)
{

    if (blockName.empty())
        return(false);

    unsigned int index = _blocks.find(blockName);
    if (index != _blocks.size())
    {
        // Found it.
        return(true);
    }
    else
    {
        // Not found.
        return(false);
    }

} // End of TableFile::IsBlockPresent()


string TableFile::GetFirstBlockName()
{

    try
    {
        return(_blocks[0].GetName());
    }
    catch (out_of_range)
    {
        throw EmptyContainerException("TableFile empty",
          "TableFile::GetFirstBlockName");
    }

#ifdef VLAD_DELETED
    // VLAD - EXCEPTIONS HERE
    if (_blocks.empty())
        return(string());

    return(_blocks[0].GetName());
#endif

} // End of TableFile::GetFirstBlockName()


void TableFile::GetBlockNames(vector<string>& blockNames)
{

    blockNames.clear();

    for (unsigned int blockI = 0; blockI < _blocks.size(); ++blockI)
    {
        blockNames.push_back(_blocks[blockI].GetName());
    }

}


string TableFile::RenameBlock(const string& oldBlockName,
  const string& newBlockName)
{

    // VLAD - Carefully examine all scenarios here. The returned value
    // may not be correct.

    if (oldBlockName == newBlockName)
        return(oldBlockName);

    if (_blocks.empty())
    {
        // No data blocks.
        throw EmptyContainerException("No data blocks in table file",
          "TableFile::RenameBlock");
    }

    unsigned int blockIndex = _blocks.find(oldBlockName);
    if (blockIndex == _blocks.size())
    {
        // Block not found.
        throw NotFoundException("Block not found",
          "TableFile::RenameBlock");
    }

    string internalBlockName = MakeInternalBlockName(newBlockName,
      blockIndex);

    _blocks.rename(oldBlockName, internalBlockName);

    return(internalBlockName);

} // End of TableFile::RenameBlock()


Block& TableFile::GetBlock(const string& blockName)
{

    return(_blocks[blockName]);

} // End of TableFile::GetBlock()


void TableFile::Flush()
{

    if ((_fileMode == READ_MODE) || (_fileMode == VIRTUAL_MODE))
    {
        throw FileModeException("Flush table file in inappropriate block mode",
          "TableFile::Flush");
    }

    // Write all modified tables 
    for (unsigned int blockI = 0; blockI < _blocks.size(); ++blockI)
    {
        for (unsigned int tableI = 0; tableI < _blocks[blockI]._tables.size();
          ++tableI)
        {
            if (&(_blocks[blockI]._tables[tableI]) == nullptr)
            {
                continue;
            }

            if (_blocks[blockI]._tables[tableI].GetModified())
            {
                _blocks[blockI]._tables[tableI].SetSerializer(_f);

                _blocks[blockI]._tables.write(_blocks[blockI].\
                  _tables[tableI].GetName());

                _blocks[blockI]._tables[tableI].SetModified(false);
            }
        }
    }

    vector<unsigned int> tableLocs;

    pair<unsigned int, unsigned int> indices;
    string name;

    for (unsigned int blockI = 0; blockI < _blocks.size(); ++blockI)
    {
        for (unsigned int tableI = 0; tableI < _blocks[blockI]._tables.size();
          ++tableI)
        {
            name = _blocks[blockI]._tables.get_name(tableI);
            indices = _blocks[blockI]._tables.get_indices(name);
 
            tableLocs.push_back(indices.second);
        }
    }

    _WriteFileIndex(_f, tableLocs);

    tableLocs.clear();

} // End of TableFile::Flush()


void TableFile::Serialize(const string& fileName)
{

    // Do we need to check if this is the same as attached file
    if (fileName == _fileName)
    {
        Flush();
        return;
    }

    if (fileName.empty())
    {
        throw EmptyValueException("Empty file name",
          "TableFile::Serialize");
    }

    auto* ser = new Serializer(fileName, CREATE_MODE);

    // Read in all the missing tables in memory
    for (unsigned int blockI = 0; blockI < _blocks.size(); ++blockI)
    {
        for (unsigned int tableI = 0; tableI < _blocks[blockI]._tables.size();
          ++tableI)
        {
            if (&(_blocks[blockI]._tables[tableI]) == nullptr)
            {
                _GetTablePtr(blockI, tableI);
            }
        }
    }

    // Write all the tables. Do not touch the flag, since this
    // flag is only related to the attached binary file and not
    // the new one.

    vector<unsigned int> tableLocs;
    int loc;

    for (unsigned int blockI = 0; blockI < _blocks.size(); ++blockI)
    {
        for (unsigned int tableI = 0; tableI < _blocks[blockI]._tables.size();
          ++tableI)
        {
            _blocks[blockI]._tables[tableI].SetSerializer(ser);
            loc = _blocks[blockI]._tables.write(_blocks[blockI].\
              _tables[tableI].GetName());
            tableLocs.push_back(loc);
        }
    }

    _WriteFileIndex(ser, tableLocs);

    tableLocs.clear();

    delete ser;

} // End of TableFile::Serialize()


void TableFile::Close()
{

    if (_fileMode == NO_MODE)
    {
        // Already closed.
        return;
    }

    if (_fileMode == VIRTUAL_MODE)
    {
        _fileMode = NO_MODE;
        return;
    }

    if (_fileMode != READ_MODE)
    {
        Flush();
    }

    if (_f != nullptr)
    {
        delete _f;
        _f = nullptr;
    }

    _fileMode = NO_MODE;

}


void TableFile::PrintHeaderInfo()
{

    cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"\
      "++++++++++++++++++++++++++" << endl;
    cout << "TableFile::PrintHeaderInfo*TableFile::PrintHeaderInfo*"\
      "TableFile::PrintHeaderInfo" << endl;
    cout << "File name: " << _fileName << endl;
    cout << "File mode: " << _fileMode << endl;
    cout << "numBlocks = " << _blocks.size() << endl;
    cout << "numTables = " << GetTotalNumTables() << endl;

    for (unsigned int i = 0; i < _blocks.size(); i++)
    {
        cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++"\
          "+++++++++++++++++" << endl;
        cout << "Block [" << i << "] " << _blocks[i].GetName() << endl; 
    }

    for (unsigned int blockI = 0; blockI < _blocks.size(); ++blockI)
    {
        for (unsigned int i = 0; i < _blocks[blockI]._tables.size(); i++)
        {
            cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++"\
              "+++++++++++++++++" << endl;
            cout << "Table [" << i << "] " << endl; 
            cout << "block = " << _blocks[blockI].GetName() << endl; 
            cout << "name = " << _blocks[blockI]._tables[i].GetName() << endl; 
            cout << "indexInMemory = " << i << endl; 
        }
    }
    cout << "TableFile::PrintHeaderInfo*TableFile::PrintHeaderInfo*"\
      "TableFile::PrintHeaderInfo" << endl;
    cout << "eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+"\
      "eod+eod+eod+eod+eod+eod" << endl;

}




void TableFile::_SetStatusInd(const string& blockName)
{

    string::size_type hashSymbolIndex = blockName.find('#');
    if (hashSymbolIndex != string::npos)
    {
        // Found a block with a hash symbol in its name
        if (hashSymbolIndex > 1)
        {
            // Set the flag that duplicate blocks are found
            _statusInd |= eDUPLICATE_BLOCKS;
        } // More than a space symbol prior to a hash 
        else
        {
            // Set the flag that empty block name is found
            _statusInd |= eUNNAMED_BLOCKS;
        } // Only a space symbol prior to a hash
    }

}


void TableFile::_AddBlock(const string& blockName, Serializer* ser)
{

    auto* blockP = new Block(blockName, ser, _fileMode, _caseSense);

    _blocks.push_back(blockP);

} // End of TableFile::_AddBlock()


// VLAD - FIX THIS TO USE UInt32 and fix Serializer API correspondingly
void TableFile::_GetNumTablesInBlocks(vector<UInt32>& numTablesInBlocks)
{

    numTablesInBlocks.clear();

    for (unsigned int blockI = 0; blockI < _blocks.size(); ++blockI)
    {
        numTablesInBlocks.push_back(_blocks[blockI]._tables.size());
    }

}


ISTable* TableFile::_GetTablePtr(const unsigned int blockIndex,
  const unsigned int tableIndex)
{

    if (&(_blocks[blockIndex]._tables[tableIndex]) != nullptr)
    {
        // Table is already in memory. Just return the pointer to it.
        return &(_blocks[blockIndex]._tables[tableIndex]);
    }
    else
    {
        // For write or virtual mode, there are no tables on the disk.
        if ((_fileMode == CREATE_MODE) || (_fileMode == VIRTUAL_MODE))
            return(nullptr);

        string name = _blocks[blockIndex]._tables.get_name(tableIndex);

        ISTable* tableP = new ISTable(name);

        tableP->SetSerializer(_f);

        _blocks[blockIndex]._tables.set(tableP);

        _blocks[blockIndex]._tables.read(name);

        return(tableP);
    }

}


void TableFile::_GetAllTables()
{

    for (unsigned int blockI = 0; blockI < _blocks.size(); ++blockI)
    {
        for (unsigned int tableI = 0; tableI < _blocks[blockI]._tables.size();
          ++tableI)
        {
            _GetTablePtr(blockI, tableI);
        }
    }

}


unsigned int TableFile::GetTotalNumTables()
{

    unsigned int totalNumTables = 0;

    for (unsigned int blockI = 0; blockI < _blocks.size(); ++blockI)
    {
        totalNumTables += _blocks[blockI]._tables.size();
    }

    return(totalNumTables);

}


void TableFile::GetTableNames(vector<string>& tableNames)
{
    // This method gets all the table names in all blocks. It may have
    // duplicate table names if tables with the same name are present in
    // multiple blocks.
    tableNames.clear();

    for (unsigned int blockI = 0; blockI < _blocks.size(); ++blockI)
    {
        for (unsigned int tableI = 0; tableI < _blocks[blockI]._tables.size();
          ++tableI)
        {
            tableNames.push_back(_blocks[blockI]._tables.get_name(tableI));
        }
    }

}


void TableFile::GetTablesIndices(vector<unsigned int>& tablesIndices)
{

    tablesIndices.clear();

    for (unsigned int blockI = 0; blockI < _blocks.size(); ++blockI)
    {
        for (unsigned int tableI = 0; tableI < _blocks[blockI]._tables.size();
          ++tableI)
        {
            tablesIndices.push_back(tableI);
        }
    }

}


void TableFile::GetSortedTablesIndices(vector<unsigned int>& tablesIndices)
{

    tablesIndices.clear();

    vector<unsigned int> partIndices;

    for (unsigned int blockI = 0; blockI < _blocks.size(); ++blockI)
    {
        _blocks[blockI]._tables.get_sorted_indices(partIndices);

        tablesIndices.insert(tablesIndices.end(), partIndices.begin(),
          partIndices.end());
    }

}


void TableFile::_ReadFileIndex()
{

    unsigned int where = _f->GetNumDataIndices() - 1;

    string version;

    _f->ReadString(version, where);

    if (version == _version)
    {
        _ReadFileIndexVersion1();
    }
    else
    {
        _ReadFileIndexVersion0();
    }

}


void TableFile::_ReadFileIndexVersion1()
{

    unsigned int where = _f->GetNumDataIndices();
    if (where <  5)
    {
        throw FileException("Read file header size is inconsistent",
          "TableFile::_ReadFileIndexVersion1");
    }

    where -= 5; 

    vector<string> blockNames;
    _f->ReadStrings(blockNames, where);
    ++where;

    vector<UInt32> numTablesInBlocks;
    _f->ReadUInt32s(numTablesInBlocks, where);
    ++where;

    vector<string> tableNames;
    _f->ReadStrings(tableNames, where);
    ++where;

    vector<UInt32> locs;
    _f->ReadUInt32s(locs, where);
    ++where;

    for (unsigned int blockI = 0, tableI = 0; blockI < blockNames.size();
      ++blockI)
    {
        _AddBlock(blockNames[blockI], _f);

        _SetStatusInd(blockNames[blockI]);

         for (unsigned int tableNumI = 0; tableNumI <
           (unsigned int)numTablesInBlocks[blockI]; ++tableNumI, ++tableI)
         {
             _blocks[blockI].AddTable(tableNames[tableI], locs[tableI]);
         }
    } // For all blocks

}


void TableFile::_ReadFileIndexVersion0()
{
    unsigned int i, numTables = 0, numBlocks = 0;

    unsigned int where = _f->GetNumDataIndices();
    if (where < 7)
    {
        throw FileException("Read file header size is inconsistent",
          "TableFile::_ReadFileIndexVersion0");
    }

    where -= 7;

    vector<string> tableIds;
    _f->ReadStrings(tableIds, where);
    numTables = tableIds.size();

    vector<string> tableNames;
    _f->ReadStrings(tableNames, where + 1);
    numTables = tableNames.size();

    vector<string> tableBlockNames;
    _f->ReadStrings(tableBlockNames, where + 2);
    numTables = tableBlockNames.size();

    vector<string> blockNames;
    _f->ReadStrings(blockNames, where + 3);
    numBlocks = blockNames.size();

    vector<UInt32> tOrder2;
    _f->ReadUInt32s(tOrder2, where + 4);

    vector<UInt32> bOrder2;
    _f->ReadUInt32s(bOrder2, where + 5);

    vector<UInt32> locs;
    _f->ReadUInt32s(locs, where + 6);

    if (bOrder2.empty())
    {
        for (i = 0; i < numBlocks; i++)
        {
            _AddBlock(blockNames[i], _f);

            _SetStatusInd(blockNames[i]);
        }
    }
    else
    {
        auto* bOrder = new int[numBlocks];
        for (i = 0; i < numBlocks; i++)
        {
            bOrder[bOrder2[i]] = i;
        }
        for (i = 0; i < numBlocks; i++)
        {
            _AddBlock(blockNames[bOrder[i]], _f);

            _SetStatusInd(blockNames[bOrder[i]]);
        }

        if (bOrder)
            delete[] bOrder; 
        bOrder = nullptr;
    }

    if (tOrder2.empty())
    {
        for (i = 0; i < numTables; i++)
        {
            /* Find block index from table Id. */ 
            unsigned int blockIndex = GetBlockIndexFromTableId(tableIds[i]);

            /* Find table name from table Id. */ 
            string tableName = GetTableNameFromTableId(tableIds[i]);

            _blocks[blockIndex].AddTable(tableName, locs[i]);
        }
    }
    else
    {
        auto* tOrder = new int[numTables];
	for (i = 0; i < numTables; i++)
        {
            tOrder[tOrder2[i]]=i;
	}
	for (i = 0; i < numTables; i++)
        {
            /* Find block index from table Id. */ 
            unsigned int blockIndex = GetBlockIndexFromTableId(tableIds[i]);

            /* Find table name from table Id. */ 
            string tableName = GetTableNameFromTableId(tableIds[tOrder[i]]);

            _blocks[blockIndex].AddTable(tableName, locs[tOrder[i]]);
        }
	if (tOrder) delete[] tOrder; 
	tOrder = nullptr;
    }

}


void TableFile::_WriteFileIndex(Serializer* ser,
  const vector<unsigned int>& tableLocs)
{

    if (ser == nullptr)
    {
        throw EmptyValueException("NULL ser pointer",
          "TableFile::_WriteFileIndex");
    }

    vector<string> blockNames;
    GetBlockNames(blockNames);
 
    vector<UInt32> numTablesInBlocks;
    _GetNumTablesInBlocks(numTablesInBlocks);

    vector<string> tableNames;
    GetTableNames(tableNames);

    ser->WriteStrings(blockNames);
    ser->WriteUInt32s(numTablesInBlocks);
    ser->WriteStrings(tableNames);
    ser->WriteUInt32s(tableLocs);

    ser->WriteString(_version);

}




void TableFile::Init()
{

    _fileMode = VIRTUAL_MODE;
    _statusInd = eCLEAR_STATUS;
    _f = nullptr;

}


void TableFile::Open(const string& fileName, const eFileMode fileMode)
{

    if (_fileMode == VIRTUAL_MODE)
    {
        // No file is to be open in virtual mode. In virtual mode
        // file will never be opened.
        return;
    }

    _f = new Serializer(fileName, fileMode);

    if (_fileMode != CREATE_MODE)
        _ReadFileIndex();

}


unsigned int TableFile::GetBlockIndexFromTableId(const string& tableId)
{

    string::size_type firstUnderscoreIndex = tableId.find('_');
    string::size_type secondUnderscoreIndex = tableId.find('_',
      firstUnderscoreIndex + 1);

    unsigned int ret = String::StringToInt(tableId.substr(
      firstUnderscoreIndex + 1, secondUnderscoreIndex - 1));

    return(ret);

}


string TableFile::GetTableNameFromTableId(const string& tableId)
{

    string::size_type firstUnderscoreIndex = tableId.find('_');
    string::size_type secondUnderscoreIndex = tableId.find('_',
      firstUnderscoreIndex + 1);

    return tableId.substr(secondUnderscoreIndex + 1, tableId.size() - 1 -
      secondUnderscoreIndex);

}

string TableFile::MakeInternalBlockName(const string& blockName,
  const unsigned int blockIndex)
{

    if (blockName.empty())
    {
        // Set the flag that empty block name is found
        _statusInd |= eUNNAMED_BLOCKS;
    }

    string internalBlockName = blockName;

    if (blockName.empty() || IsBlockPresent(blockName))
    {
        // This is either a nameless block or a block that is already present.

        // Set the flag that duplicate blocks are found
        if ((blockName.empty() && (_blocks.size() > 1)) ||
          IsBlockPresent(blockName))
        {
            _statusInd |= eDUPLICATE_BLOCKS;
        }

        // Create a new block name.
        internalBlockName += " # ";
        internalBlockName += String::IntToString(blockIndex);
    }

    return(internalBlockName);

}
