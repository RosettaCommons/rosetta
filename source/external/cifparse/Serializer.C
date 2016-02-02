/*
FILE:     Serializer.C
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


#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>

// Windows work around for Posix commands
#ifndef WIN32
#include <unistd.h>
#else // On Win32
#include <io.h>
#define S_IRUSR _S_IREAD
#define S_IWUSR _S_IWRITE
#endif

#include <iostream>
#include <iomanip>
#include <stdexcept>
#include <string>

#include "rcsb_types.h"
#include "Exceptions.h"
#include "GenString.h"
#include "RcsbPlatform.h"
#include "Serializer.h"


using std::string;
using std::vector;
using std::ios;
using std::endl;
using std::cerr;
using std::setw;
using std::out_of_range;


bool Serializer::_littleEndian = RcsbPlatform::IsLittleEndian();


Serializer::Serializer(const string& fileName,
  const eFileMode fileMode, const bool verbose)
{
    if (fileName.empty())
    {
        // Empty file name.
        throw EmptyValueException("Empty file name",
          "BlockIO::BlockIO");
    }

    if ((fileMode != READ_MODE) && (fileMode != CREATE_MODE) &&
      (fileMode != UPDATE_MODE))
    {
        throw FileModeException(string("Invalid file mode: ") + \
          String::IntToString(fileMode), "BlockIO::BlockIO");
    }

    _numBlocksIO = 0;
    _currentBlockIO = 0;

    OpenFileIO(fileName, fileMode);

    Init();

    _fileName = fileName;

    _mode = fileMode;

    _verbose = verbose;
    if (_verbose)
    {
        _log.open("Serializer.log", ios::out | ios::app);
    }

    _theBlock.AssociateBuffer(&_buffer);

    if (GetNumBlocksIO() >= 3)
    {
        // Read file header if it exists. It will not exist for files
        // opened in create mode, where there are fewer than 3 blocks.
        _ReadFileHeader();
    }

}


Serializer::~Serializer()
{
    if (_mode != READ_MODE)
    {
        // Finish writing data in the current block
        WriteBlock(_currentBlock);

        _WriteFileHeader();
    }

    if (_verbose)
        _log.close();

    CloseFileIO();
}


UInt32 Serializer::ReadUInt32(const UInt32 index)
{

    if (_verbose)
        _log << "ReadUInt32() index = " << index << endl;

    if ((_mode != READ_MODE) && (_mode != UPDATE_MODE))
    {
        throw FileModeException("Read attempt in write-only file",
          "Serializer::ReadUInt32");
    }

    if (index >= _indices.size())
    {
        throw out_of_range("Invalid index in Serializer::ReadUInt32");
    }

    if (_indices[index].blockNumber == 0)
    {
        throw InvalidStateException("Attempt to read deleted index",
          "Serializer::ReadUInt32");
    }

    if ((_indices[index].dataType != UWORD_TYPE) && 
      (_indices[index].dataType != WORD_TYPE))
    {
        throw InvalidStateException("Attempt to read non-UInt32",
          "Serializer::ReadUInt32");
    }

    _currentBlock = _indices[index].blockNumber;
    UInt32 bytesRead = ReadBlock(_indices[index].blockNumber);
    if (bytesRead < _indices[index].offset + UINT32_SIZE)
    {
        throw FileException("Reading less that than needed for UInt32",
          "Serializer::ReadUInt32");
    }

    char* temp = _buffer + _indices[index].offset;
    if ((temp - _buffer) > (int)(BLKSIZE - UINT32_SIZE))
    {
        throw FileException("Reading less that than needed for UInt32",
          "Serializer::ReadUInt32");
    }

    return (_GetUInt32(temp));

}

void Serializer::ReadUInt32s(vector<UInt32>& UInt32s, const UInt32 index)
{
    UInt32s.clear();

    if (_verbose)
        _log << "ReadUInt32s() index = " << index << endl;

    if ((_mode != READ_MODE) && (_mode != UPDATE_MODE))
    {
        throw FileModeException("Read attempt in write-only file",
          "Serializer::ReadUInt32s");
    }

    if (index >= _indices.size())
    {
        throw out_of_range("Invalid index in Serializer::ReadUInt32s");
    }

    if (_indices[index].blockNumber == 0)
    {
        throw InvalidStateException("Attempt to read deleted index",
          "Serializer::ReadUInt32s");
    }

    if ((_indices[index].dataType != UWORDS_TYPE) && 
      (_indices[index].dataType != WORDS_TYPE))

    {
        throw InvalidStateException("Attempt to read non-UInt32s",
          "Serializer::ReadUInt32s");
    }

    _currentBlock = _indices[index].blockNumber;
    UInt32 bytesRead = ReadBlock(_indices[index].blockNumber);
    if (bytesRead < _indices[index].offset + UINT32_SIZE)
    {
        throw FileException("Reading less that than needed for UInt32s",
          "Serializer::ReadUInt32s");
    }

    char* temp = _buffer + _indices[index].offset;

    UInt32 numWords = _GetUInt32(temp);

    temp += UINT32_SIZE;
    if (numWords == 0)
    {
        return;
    }

    UInt32 wordsToRead = numWords;
    if (((numWords + 1)*UINT32_SIZE + _indices[index].offset) > bytesRead
        && bytesRead != BLKSIZE)
    {
        throw FileException("Reading less that than needed for UInt32s",
          "Serializer::ReadUInt32s");
    }

    if (((numWords + 1)*UINT32_SIZE) != _indices[index].length)
    {
        throw InvalidStateException("Invalid index length",
          "Serializer::ReadUInt32s");
    }

    UInt32 blockSpan = (_indices[index].length + _indices[index].offset - 1)
        / BLKSIZE + 1;
    UInt32 boundary = BLKSIZE / UINT32_SIZE;
    UInt32 endBlockNumber = _indices[index].blockNumber + blockSpan - 1;

    while (blockSpan--)
    {
        if (_currentBlock == endBlockNumber)
        {
            if (bytesRead < wordsToRead*UINT32_SIZE)
            {
                numWords = 0;
                throw FileException("Reading less that than needed for UInt32s",
                  "Serializer::ReadUInt32s");
            }
            else
                boundary = wordsToRead;
        }
        else
        {
            boundary = (_buffer + BLKSIZE - temp) / UINT32_SIZE;
        }
        for (UInt32 i = 0; i < boundary; i++)
        {
            UInt32s.push_back(_GetUInt32(temp));

            temp += UINT32_SIZE;
        }
        wordsToRead -= boundary;
        if (blockSpan)
        {
            bytesRead = ReadBlock(++_currentBlock);
            temp = _buffer;
        }
    }

    if (wordsToRead)
    {
        throw InvalidStateException("Invalid number of UInt32s",
          "Serializer::ReadUInt32s");
    }

}


void Serializer::ReadString(string& retString, const UInt32 index)
{

    retString.clear();

    if (_verbose)
        _log << "ReadString() index = " << index << endl;

    if ((_mode != READ_MODE) && (_mode != UPDATE_MODE))
    {
        throw FileModeException("Read attempt in write-only file",
          "Serializer::ReadString");
    }

    if (index >= _indices.size())
    {
        throw out_of_range("Invalid index in Serializer::ReadString");
    }

    if (_indices[index].blockNumber == 0)
    {
        throw InvalidStateException("Attempt to read deleted index",
          "Serializer::ReadString");
    }

    if (_indices[index].dataType != STRING_TYPE)
    {
        throw InvalidStateException("Attempt to read non-string",
          "Serializer::ReadString");
    }

    _currentBlock = _indices[index].blockNumber;
    UInt32 bytesRead = ReadBlock(_indices[index].blockNumber);

    if (bytesRead < _indices[index].offset + UINT32_SIZE)
    {
        if (_verbose)
        {
            _log << " -X0---------------------------------------------------"
                << endl;
            _log << "UINT32_SIZE   =              " << UINT32_SIZE << endl;
            _log << "BLKSIZE     =             " << BLKSIZE << endl;
            _log << "_indices[index].offset  = " << _indices[index].offset
                << endl;
            _log << "bytesRead   =             " << bytesRead << endl;
        }

        throw FileException("Reading less that than needed for String",
          "Serializer::ReadString");
    }

    char* temp = _buffer + _indices[index].offset;

    UInt32 stringSize = _GetUInt32(temp);

    temp += UINT32_SIZE;
    UInt32 bytesToRead = stringSize;
    if ((stringSize + UINT32_SIZE + _indices[index].offset) > bytesRead
        && bytesRead != BLKSIZE)
    {
        throw FileException("Reading less that than needed for String",
          "Serializer::ReadString");
    }
    if ((stringSize + UINT32_SIZE) != _indices[index].length)
    {
        throw InvalidStateException("Invalid index length",
          "Serializer::ReadString");
    }

    UInt32 blockSpan = (_indices[index].length + _indices[index].offset - 1)
        / BLKSIZE + 1;
    UInt32 boundary = BLKSIZE;
    UInt32 endBlockNumber = _indices[index].blockNumber + blockSpan - 1;

    while (blockSpan--)
    {
        if (_currentBlock == endBlockNumber)
        {
            if (bytesRead < bytesToRead)
            {
                throw FileException("Reading less that than needed for String",
                  "Serializer::ReadString");
            }
            else
                boundary = bytesToRead;
        }
        else
        {
            boundary = _buffer + BLKSIZE - temp;
        }
        for (UInt32 i = 0; i < boundary; i++)
            retString.push_back(*temp++);

        bytesToRead -= boundary;
        if (blockSpan)
        {
            bytesRead = ReadBlock(++_currentBlock);
            temp = _buffer;
        }
    }

    if (bytesToRead)
    {
        retString.clear();

        throw InvalidStateException("Invalid number of String octets",
          "Serializer::ReadString");
    }

}

void Serializer::ReadStrings(vector<string>& theStrings, const UInt32 index)
{

    theStrings.clear();

    if (_verbose)
        _log << "ReadStrings() index = " << index << endl;

    if ((_mode != READ_MODE) && (_mode != UPDATE_MODE))
    {
        throw FileModeException("Read attempt in write-only file",
          "Serializer::ReadStrings");
    }

    if (index >= _indices.size())
    {
        throw out_of_range("Invalid index in Serializer::ReadStrings");
    }

    if (_indices[index].blockNumber == 0)
    {
        throw InvalidStateException("Attempt to read deleted index",
          "Serializer::ReadStrings");
    }
    if (_indices[index].dataType != STRINGS_TYPE)
    {
        throw InvalidStateException("Attempt to read non-Strings",
          "Serializer::ReadStrings");
    }

    vector<UInt32> stringSizes;

    _currentBlock = _indices[index].blockNumber;
    UInt32 bytesRead = ReadBlock(_indices[index].blockNumber);

    char* temp = _buffer + _indices[index].offset;

    UInt32 numStrings = _GetUInt32(temp);

    temp += UINT32_SIZE;

    if (numStrings == 0)
    {
        return;
    }

    UInt32 wordsLeftInBlock = (_buffer + BLKSIZE - temp) / UINT32_SIZE;
    for (UInt32 i = 0; i < numStrings; ++i)
    {
        if (wordsLeftInBlock == 0)
        {
            bytesRead = ReadBlock(++_currentBlock);
            temp = _buffer;
            wordsLeftInBlock = (_buffer + BLKSIZE - temp) / UINT32_SIZE;
            if ((bytesRead < (_indices[index].length - (i+2)*UINT32_SIZE))
                && (bytesRead != BLKSIZE))
            {
                throw FileException("Reading less that than needed for Strings",
                  "Serializer::ReadStrings");
            }
        }

        stringSizes.push_back(_GetUInt32(temp));

        temp += UINT32_SIZE;
        wordsLeftInBlock--;
    }

    for (UInt32 i = 0; i < numStrings; ++i)
    {
        theStrings.push_back("");

#ifndef VLAD_REPORTED_MEMORY_ACCESS_IN_PURIFY_IF_THIS_CODE_IS_ENABLED
        theStrings[i].reserve(stringSizes[i]);
#endif
        UInt32 j = stringSizes[i];
        if ((bytesRead < (j + temp - _buffer + BLKSIZE)) && (bytesRead
            != BLKSIZE))
        {
            theStrings.clear();

            throw FileException("Reading less that than needed for Strings",
              "Serializer::ReadStrings");
        }

        while (j--)
        {
            if (temp >= _buffer + BLKSIZE)
            {
                bytesRead = ReadBlock(++_currentBlock);
                temp = _buffer;
                if ((bytesRead < j) && (bytesRead != BLKSIZE))
                {
                    theStrings.clear();
                    throw FileException("Reading less that than needed for Strings",
                      "Serializer::ReadStrings");
                }
            }
            theStrings[i] += (*temp++);
        }
    }

}


UInt32 Serializer::WriteUInt32(const UInt32 theWord)
{
    UInt32 temp = _indices.size();

    WriteUInt32AtIndex(theWord, _indices.size());

    SetVirtualLength(temp);

    return(temp);
}


UInt32 Serializer::WriteUInt32s(const vector<UInt32>& theWords)
{
    UInt32 temp = _indices.size();

    WriteUInt32sAtIndex(theWords, _indices.size());

    SetVirtualLength(temp);
    
    return(temp);
}

UInt32 Serializer::WriteString(const string& theString)
{
    UInt32 temp = _indices.size();

    WriteStringAtIndex(theString, _indices.size());

    SetVirtualLength(temp);

    return (temp);
}

UInt32 Serializer::WriteStrings(const vector<string>& theStrings)
{
    UInt32 temp = _indices.size();

    WriteStringsAtIndex(theStrings, _indices.size());

    SetVirtualLength(temp);

    return (temp);
}

UInt32 Serializer::UpdateUInt32(const UInt32 theWord, const UInt32 oldIndex)
{
    if (_mode == READ_MODE)
    {
        throw FileModeException("Update attempt in read-only file",
          "Serializer::UpdateUInt32");
    }

    if (oldIndex >= _indices.size())
    {
        throw out_of_range("Invalid index in Serializer::UpdateUInt32");
    }

    if ((_indices[oldIndex].dataType != UWORD_TYPE) &&
      (_indices[oldIndex].dataType != WORD_TYPE))
    {
        throw InvalidStateException("Attempt to update non-UInt32",
          "Serializer::UpdateUInt32");
    }

    if (_indices[oldIndex].blockNumber == 0)
    {
        throw InvalidStateException("Attempt to update deleted index",
          "Serializer::UpdateUInt32");
    }

    UInt32 newIndex = oldIndex;  // This used to be -1

    if (_indices[oldIndex].vLength < UINT32_SIZE)
    {
        Delete(oldIndex);
        newIndex = WriteUInt32(theWord);
    }
    else
    {
        WriteUInt32AtIndex(theWord, oldIndex);
    }

    return (newIndex);
}

UInt32 Serializer::UpdateUInt32s(const vector<UInt32>& theWords,
  const UInt32 oldIndex)
{
    if (_mode == READ_MODE)
    {
        throw FileModeException("Update attempt in read-only file",
          "Serializer::UpdateUInt32s");
    }

    if (oldIndex >= _indices.size())
    {
        throw out_of_range("Invalid index in Serializer::UpdateUInt32s");
    }

    if ((_indices[oldIndex].dataType != UWORDS_TYPE) &&
      (_indices[oldIndex].dataType != WORDS_TYPE))
    {
        throw InvalidStateException("Attempt to update non-UInt32s",
          "Serializer::UpdateUInt32s");
    }

    if (_indices[oldIndex].blockNumber == 0)
    {
        throw InvalidStateException("Attempt to update deleted index",
          "Serializer::UpdateUInt32s");
    }

    if (theWords.empty())
    {
        throw EmptyValueException("Empty words vector",
          "Serializer::UpdateUInt32s");
        // return UWORDS_NULL_BEFORE_WRITE;
    }

    UInt32 newIndex = oldIndex; // VLAD - This used to be -1

    UInt32 totalLength = (theWords.size() + 1) * UINT32_SIZE;

    if (_indices[oldIndex].vLength < totalLength)
    {
        Delete(oldIndex);
        newIndex = WriteUInt32s(theWords);
    }
    else
    {
        WriteUInt32sAtIndex(theWords, oldIndex);
    }

    return (newIndex);
}

UInt32 Serializer::UpdateString(const string& theString,
  const UInt32 oldIndex)
{
    if (_mode == READ_MODE)
    {
        throw FileModeException("Update attempt in read-only file",
          "Serializer::UpdateString");
    }

    if (oldIndex >= _indices.size())
    {
        throw out_of_range("Invalid index in Serializer::UpdateString");
    }

    if (_indices[oldIndex].dataType != STRING_TYPE)
    {
        throw InvalidStateException("Attempt to update non-String",
          "Serializer::UpdateString");
    }

    if (_indices[oldIndex].blockNumber == 0)
    {
        throw InvalidStateException("Attempt to update deleted index",
          "Serializer::UpdateString");
    }

    //  if (!theString) return STRING_NULL_BEFORE_WRITE;

    UInt32 newIndex = oldIndex;  // VLAD - This used to be -1

    if (_indices[oldIndex].vLength < (theString.size() + UINT32_SIZE))
    {
        Delete(oldIndex);
        newIndex = WriteString(theString);
    }
    else
    {
        WriteStringAtIndex(theString, oldIndex);
    }

    return (newIndex);
}

UInt32 Serializer::UpdateStrings(const vector<string>& theStrings,
  const UInt32 oldIndex)
{
    if (_mode == READ_MODE)
    {
        throw FileModeException("Update attempt in read-only file",
          "Serializer::UpdateStrings");
    }

    if (oldIndex >= _indices.size())
    {
        throw out_of_range("Invalid index in Serializer::UpdateStrings");
    }

    if (_indices[oldIndex].dataType != STRINGS_TYPE)
    {
        throw InvalidStateException("Attempt to update non-String",
          "Serializer::UpdateStrings");
    }

    if (_indices[oldIndex].blockNumber == 0)
    {
        throw InvalidStateException("Attempt to update deleted index",
          "Serializer::UpdateStrings");
    }

    if (theStrings.empty())
    {
        throw EmptyValueException("Empty strings vector",
          "Serializer::UpdateStrings");
        // return STRING_NULL_BEFORE_WRITE;
    }

    UInt32 newIndex = oldIndex; // VLAD - This used to be -1

    UInt32 totalLength = 0;

    for (UInt32 i = 0; i < theStrings.size(); i++)
    {
        totalLength += theStrings[i].size();
    }

    totalLength += (theStrings.size() + 1) * UINT32_SIZE;

    if (_indices[oldIndex].vLength < totalLength)
    {
        Delete(oldIndex);
        newIndex = WriteStrings(theStrings);
    }
    else
    {
        WriteStringsAtIndex(theStrings, oldIndex);
    }

    return (newIndex);
}

void Serializer::Delete(const UInt32 index)
{
    if (_mode == READ_MODE)
    {
        throw FileModeException("Delete index attempt in read-only file",
          "Serializer::Delete");
    }

    if (index >= _indices.size())
    {
        throw out_of_range("Invalid index in Serializer::Delete");
    }

    _indices[index].blockNumber = 0;

    // VLAD - WE BETTER ADD CODE TO MOVE BACK THE POINTER AND REUSE THE
    // SPACE AT THE END OF THE FILE. IN OTHER WORDS IF THE DELETED INDEX
    // EQUALS _indices.size() - 1, do special processing
}

void Serializer::WriteUInt32AtIndex(const UInt32 theWord, const UInt32 index)
{
    if (_mode == READ_MODE)
    {
        throw FileModeException("Write attempt in read-only file",
          "Serializer::WriteUInt32AtIndex");
    }

    if (index > _indices.size())
    {
        throw out_of_range("Invalid index in Serializer::WriteUInt32AtIndex");
    }

    char* temp = GetWritingPoint(index);

    _indices[index].blockNumber = _currentBlock;
    _indices[index].offset = _currentOffset;
    _indices[index].length = UINT32_SIZE;
    _indices[index].dataType = UWORD_TYPE;

    _PutUInt32(theWord, temp);

    temp += UINT32_SIZE;

    WriteLast(temp);
}

void Serializer::WriteUInt32sAtIndex(const vector<UInt32>& theWords,
  const UInt32 index)
{
    if (_mode == READ_MODE)
    {
        throw FileModeException("Write attempt in read-only file",
          "Serializer::WriteUInt32sAtIndex");
    }

    if (index > _indices.size())
    {
        throw out_of_range("Invalid index in "\
          "Serializer::WriteUInt32sAtIndex");
    }

    char* temp = GetWritingPoint(index);

    // Number of words + the size of each word
    UInt32 totalLength = UINT32_SIZE * (theWords.size() + 1);

    _indices[index].blockNumber = _currentBlock;
    _indices[index].offset = _currentOffset;
    _indices[index].length = totalLength;
    _indices[index].dataType = UWORDS_TYPE;

    // First write the number of words
    _PutUInt32(theWords.size(), temp);

    temp += UINT32_SIZE;

    UInt32 wordsLeft = (_buffer + BLKSIZE - temp) / UINT32_SIZE;

    for (UInt32 i = 0; i < theWords.size(); i++)
    {
        if (wordsLeft == 0)
        {
            WriteBlock(_currentBlock++);

            temp = _buffer;
            wordsLeft = BLKSIZE / UINT32_SIZE;
        }

        // Write the current word
        _PutUInt32(theWords[i], temp);

        temp += UINT32_SIZE;
        wordsLeft--;
    }

    WriteLast(temp);
}

void Serializer::WriteStringAtIndex(const string& theString,
  const UInt32 index)
{
    if (_mode == READ_MODE)
    {
        throw FileModeException("Write attempt in read-only file",
          "Serializer::WriteStringAtIndex");
    }

    if (index > _indices.size())
    {
        throw out_of_range("Invalid index in "\
          "Serializer::WriteStringAtIndex");
    }

    char* temp = GetWritingPoint(index);

    _indices[index].blockNumber = _currentBlock;
    _indices[index].offset = _currentOffset;
    _indices[index].length = UINT32_SIZE + theString.size();
    _indices[index].dataType = STRING_TYPE;

    // First write the string size
    _PutUInt32(theString.size(), temp);

    temp += UINT32_SIZE;

    UInt32 octetsLeft = _buffer + BLKSIZE - temp;

    for (UInt32 i = 0; i < theString.size(); i++)
    {
        if (octetsLeft == 0)
        {
            WriteBlock(_currentBlock++);

            temp = _buffer;
            octetsLeft = BLKSIZE;
        }

        // Write the current character
        *temp++ = theString[i];
        octetsLeft--;
    }

    WriteLast(temp);
}

void Serializer::WriteStringsAtIndex(const vector<string>& theStrings,
  const UInt32 index)
{
    if (_mode == READ_MODE)
    {
        throw FileModeException("Write attempt in read-only file",
          "Serializer::WriteStringsAtIndex");
    }

    if (index > _indices.size())
    {
        throw out_of_range("Invalid index in "\
          "Serializer::WriteStringsAtIndex");
    }

    char* temp = GetWritingPoint(index);

    UInt32 totalLength = 0;
    for (UInt32 i = 0; i < theStrings.size(); i++)
    {
        totalLength += theStrings[i].size();
    }

    // number of strings + size of each string
    totalLength += UINT32_SIZE * (theStrings.size() + 1);

    _indices[index].blockNumber = _currentBlock;
    _indices[index].offset = _currentOffset;
    _indices[index].length = totalLength;
    _indices[index].dataType = STRINGS_TYPE;

    // First write the number of strings
    _PutUInt32(theStrings.size(), temp);

    temp += UINT32_SIZE;

    // Write the size of each string
    UInt32 wordsLeft = (_buffer + BLKSIZE - temp) / UINT32_SIZE;
    for (UInt32 i = 0; i < theStrings.size(); i++)
    {
        if (wordsLeft == 0)
        {
            WriteBlock(_currentBlock++);

            temp = _buffer;
            wordsLeft = BLKSIZE / UINT32_SIZE;
        }

        // Write the current string size
        _PutUInt32(theStrings[i].size(), temp);

        temp += UINT32_SIZE;
        wordsLeft--;
    }

    UInt32 octetsLeft = _buffer + BLKSIZE - temp;
    for (UInt32 i = 0; i < theStrings.size(); i++)
    {
        for (UInt32 j = 0; j < theStrings[i].size(); j++)
        {
            if (octetsLeft == 0)
            {
                WriteBlock(_currentBlock++);

                temp = _buffer;
                octetsLeft = BLKSIZE;
            }
            *temp++ = theStrings[i][j]; //see if incrementing is faster
            octetsLeft--;
        }
    }

    WriteLast(temp);
}

void Serializer::SwapHeader(tFileHeader& out, const tFileHeader& in)
{
    out.fileIndexBlock = SwapUInt32(in.fileIndexBlock);
    out.fileIndexNumBlocks = SwapUInt32(in.fileIndexNumBlocks);
    out.fileIndexLength = SwapUInt32(in.fileIndexLength);
    out.numIndices = SwapUInt32(in.numIndices);
    out.reserved[0] = SwapUInt32(in.reserved[0]);
    out.reserved[1] = SwapUInt32(in.reserved[1]);
    out.reserved[2] = SwapUInt32(in.reserved[2]);
    out.version = SwapUInt32(in.version);
}

void Serializer::SwapIndex(EntryIndex& out, const EntryIndex& in)
{
    out.blockNumber = SwapUInt32(in.blockNumber);
    out.offset = SwapUInt32(in.offset);
    out.length = SwapUInt32(in.length);
    out.dataType = SwapUInt32(in.dataType);
    out.vLength = SwapUInt32(in.vLength);
    out.reserved[0] = SwapUInt32(in.reserved[0]);
    out.reserved[1] = SwapUInt32(in.reserved[1]);
    out.reserved[2] = SwapUInt32(in.reserved[2]);
}

UInt32 Serializer::SwapUInt32(const UInt32 theWord)
{
    UInt32 r;
    char * sp, * dp;

    sp = (char*)&theWord;
    dp = (char*)&r;

    for (unsigned int i=1; i<sizeof(UInt32); i++)
        *dp++=sp[sizeof(UInt32)-i];
    *dp=sp[0];
    return r;
}

void Serializer::_GetHeader(const char* where)
{
    tFileHeader& fileHeader = *((tFileHeader*)where);

    if (_littleEndian)
        _fileHeader = fileHeader;
    else
        SwapHeader(_fileHeader, fileHeader);
}

void Serializer::_PutHeader(char* where)
{
    tFileHeader& fileHeader = *((tFileHeader*)where);

    if (_littleEndian)
        fileHeader = _fileHeader;
    else
        SwapHeader(fileHeader, _fileHeader);
}

void Serializer::_GetIndex(EntryIndex& outIndex, const char* where)
{
    EntryIndex& index = *((EntryIndex*)where);

    if (_littleEndian)
        outIndex = index;
    else
        SwapIndex(outIndex, index);
}

void Serializer::_PutIndex(const EntryIndex& inIndex, char* where)
{
    EntryIndex& index = *((EntryIndex*)where);

    if (_littleEndian)
        index = inIndex;
    else
        SwapIndex(index, inIndex);
}

UInt32 Serializer::_GetUInt32(const char* where)
{
    UInt32& word = *((UInt32*)where);

    if (_littleEndian)
        return (word);
    else
        return (SwapUInt32(word));
}

void Serializer::_PutUInt32(const UInt32 inWord, char* where)
{
    UInt32& word = *((UInt32*)where);

    if (_littleEndian)
        word = inWord;
    else
        word = SwapUInt32(inWord);
}

void Serializer::Init()
{
    _fileHeader.fileIndexBlock = 2;
    _fileHeader.fileIndexNumBlocks = 1;
    _fileHeader.fileIndexLength = 0;
    _fileHeader.numIndices = 0;
    _fileHeader.reserved[0] = 0;
    _fileHeader.reserved[1] = 0;
    _fileHeader.reserved[2] = 0;
    _fileHeader.version = _version;

    _indices.reserve(INDEX_INCREMENT);

    _currentBlock = 1;
    _currentOffset = 0;
}

void Serializer::PrintIndexPosition(const UInt32 position)
{
    if (_verbose)
        _log << "Serializer Position " << position << " blockNumber "
            << _indices[position].blockNumber << " offset      "
            << _indices[position].offset << " length      "
            << _indices[position].length << " vLength     "
            << _indices[position].vLength << " dataType    "
            << _indices[position].dataType << endl;
}

void Serializer::PrintIndex()
{
    if (_verbose)
    {
        _log
            << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
            << endl;
        _log
            << "Serializer::PrintIndex*Serializer::PrintIndex*Serializer::PrintIndex"
            << endl;
        _log << "File Navigator Index" << endl;
        _log << "_fileName:         " << _fileName << endl;
        _log << "_mode:             " << _mode << endl;
        _log << "_currentBlock:     " << _currentBlock << endl;
        _log << "_currentOffset:    " << _currentOffset << endl;
        _log << "_fileHeader.fileIndexBlock:     "
            << _fileHeader.fileIndexBlock << endl;
        _log << "_fileHeader.fileIndexNumBlocks: "
            << _fileHeader.fileIndexNumBlocks << endl;
        _log << "_fileHeader.fileIndexLength:    "
            << _fileHeader.fileIndexLength << endl;
        _log << "------------------------------------------------------"
            << endl;
    }

    for (unsigned int i = 0; i < _indices.size(); ++i)
    {
        PrintIndexPosition(i);
    }

    if (_verbose)
    {
        _log
            << "Serializer::PrintIndex*Serializer::PrintIndex*Serializer::PrintIndex"
            << endl;
        _log
            << "eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod"
            << endl << endl;
    }
}

void Serializer::DumpFile()
{
    for (UInt32 position = 0; position < _indices.size(); ++position)
    {
        if (_verbose)
        {
            _log
                << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++"
                << endl;
            _log << "File Position [" << position << "] blockNumber "
                << _indices[position].blockNumber << " dataType    "
                << _indices[position].dataType << " offset      "
                << _indices[position].offset << " length      "
                << _indices[position].length << " vLength     "
                << _indices[position].vLength << endl;
        }

        if (_indices[position].blockNumber == 0)
        {
            if (_verbose)
                _log << "Data deleted" << endl;
            continue;
        }

        switch (_indices[position].dataType)
        {
            case STRINGS_TYPE:
            {
                vector<string> ss;
                ReadStrings(ss, position);
                UInt32 num = ss.size();
                if (_verbose)
                    _log << "   Strings = " << num << endl;
                for (UInt32 j = 0; j < num; ++j)
                {
                    if (_verbose)
                        _log << "   [" << j << "] " << ss[j] << endl;
                }
                ss.clear();
                break;
            }
            case STRING_TYPE:
            {
                string s;
                ReadString(s, position);
                if (_verbose)
                {
                    if (s.empty())
                        _log << "   String is NULL" << endl;
                    else
                        _log << "   String[" << s.size() << "] = " << s
                            << endl;
                }
                break;
            }
#ifdef VLAD_DECIDE_LATER
            case WORD_TYPE:
                Word w = 0;
                w = GetWord(position, err);
                if (_verbose)
                    _log << "   Word = " << w << endl;
                break;

            case WORDS_TYPE:
                vector<Word> ws;
                GetWords(position, ws, err);
                if (_verbose)
                    _log << "   Words = " << ws.size() << endl;
                for (UInt32 j = 0; j < ws.size(); ++j)
                {
                    if (_verbose)
                        _log << "   [" << j << "] " << ws[j] << endl;
                }

                ws.clear();
                break;
#endif

            case UWORD_TYPE:
            {
                UInt32 uw = ReadUInt32(position);
                if (_verbose)
                    _log << "   UInt32 = " << uw << endl;
                break;
            }
            case UWORDS_TYPE:
            {
                vector<UInt32> uws;
                ReadUInt32s(uws, position);
                if (_verbose)
                    _log << "   UInt32s = " << uws.size() << endl;
                for (UInt32 j = 0; j < uws.size(); ++j)
                {
                    if (_verbose)
                        _log << "   [" << j << "] " << uws[j] << endl;
                }
                if (!uws.empty())
                    uws.clear();
                break;
            }
            default:
            {
                if (_verbose)
                    _log << "No method for data type "
                        << _indices[position].dataType << endl;
                break;
            }
        }
    }

    if (_verbose)
        _log
            << "+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+eod+"
            << endl;

}

void Serializer::_WriteFileHeader()
{
    _fileHeader.numIndices = _indices.size();

    _fileHeader.fileIndexNumBlocks = (_fileHeader.numIndices
        / _indicesPerBlock) + 1;
    if (_fileHeader.numIndices % _indicesPerBlock == 0)
    {
        _fileHeader.fileIndexNumBlocks--;
    }

    _fileHeader.fileIndexLength = _fileHeader.numIndices * sizeof(EntryIndex);

    GetLastDataBuffer();

    if (_currentOffset != 0)
    {
        _currentOffset = 0;
        _currentBlock++;
    }

    _fileHeader.fileIndexBlock = _currentBlock;

    char* temp = _buffer;

    UInt32 indicesInBlock = _indicesPerBlock;

    for (UInt32 i = 0; i < _fileHeader.fileIndexNumBlocks; ++i)
    {
        if (i == (_fileHeader.fileIndexNumBlocks - 1))
        {
            indicesInBlock = _fileHeader.numIndices - (_indicesPerBlock * i);
        }

        for (UInt32 j = 0; j < indicesInBlock; ++j)
        {
            _PutIndex(_indices[i * _indicesPerBlock + j], temp);

            temp += sizeof(EntryIndex);
        }

        WriteBlock(_currentBlock++);
        temp = _buffer;
    }

    // Set the buffer to all zeroes
    memset(_buffer, 0, BLKSIZE);

    _PutHeader(_buffer);

    WriteBlock(0);
}

void Serializer::_ReadFileHeader()
{
    UInt32 bytesRead = ReadBlock(0);
    if (bytesRead < sizeof(tFileHeader))
    {
        throw FileException("Read file header is too short",
          "Serializer::_ReadFileHeader");
    }

    _GetHeader(_buffer);

    if (_fileHeader.numIndices == 0)
    {
        return;
    }
    
    if (((_fileHeader.fileIndexLength/sizeof(EntryIndex)) !=
      _fileHeader.numIndices) || (_fileHeader.fileIndexBlock < 2)
        || (_fileHeader.fileIndexNumBlocks < 1))
    {
        // Not enough bytes read
        throw FileException("Read file header size is inconsistent",
          "Serializer::_ReadFileHeader");
    }

    UInt32 indicesInBlock = _indicesPerBlock;

    for (UInt32 i = 0; i < _fileHeader.fileIndexNumBlocks; ++i)
    {
        if (i == (_fileHeader.fileIndexNumBlocks - 1))
        {
            indicesInBlock = _fileHeader.numIndices - (_indicesPerBlock * i);
        }

        bytesRead = ReadBlock(_fileHeader.fileIndexBlock + i);

        if (bytesRead < (indicesInBlock * sizeof(EntryIndex)))
        {
            throw FileException("File header content is inconsistent",
              "Serializer::_ReadFileHeader");
        }

        char* temp = _buffer;

        for (UInt32 j = 0; j < indicesInBlock; ++j, temp += sizeof(EntryIndex))
        {
            EntryIndex tmpindex;
            _GetIndex(tmpindex, temp);

            if (tmpindex.blockNumber != 0)
            {
                // Only store non-deleted indices.
                _indices.push_back(tmpindex);
            }
        }
    }

    _currentBlock = _indices[_indices.size() - 1].blockNumber;

    _currentOffset = _indices[_indices.size() - 1].offset +
      _indices[_indices.size() - 1].vLength;
}

void Serializer::GetLastDataBuffer(void)
{
    if (_mode == READ_MODE)
    {
        throw FileModeException("Write attempt in read-only file",
          "Serializer::GetLastDataBuffer");
    }

    if (_indices.empty())
    {
        // The first item is always written in block number 1 at offset 0
        _currentBlock = 1;
        _currentOffset = 0;

        return;
    } 

    UInt32 prevBlock = _currentBlock;

    _currentBlock = _indices[_indices.size() - 1].blockNumber;

    _currentOffset = _indices[_indices.size() - 1].offset
        + _indices[_indices.size() - 1].vLength;

    int a1 = _buffer + _currentOffset - (char*)0; // word align the offset

    UInt32 n = a1 % UINT32_SIZE;

    if (n != 0)
        _currentOffset += (UINT32_SIZE - n);

    while (_currentOffset >= BLKSIZE)
    {
        _currentOffset -= BLKSIZE;
        _currentBlock++;

        a1 = _buffer + _currentOffset - (char*)0; // word align the offset

        n = a1 % UINT32_SIZE;

        if (n != 0)
            _currentOffset += (UINT32_SIZE - n);
    }

    if (prevBlock != _currentBlock)
    {
        WriteBlock(prevBlock);

        UInt32 bytesRead = ReadBlock(_currentBlock);
        if (bytesRead != BLKSIZE)
        {
            if (_verbose)
                _log << "**Quitting with bytesRead = " << bytesRead << endl;

            throw FileException("Reading less than block size",
              "Serializer::GetLastDataBuffer");
            //return NOT_ENOUGH_BYTES_READ;
        }
    }
}

void Serializer::GetDataBufferAtIndex(const UInt32 index)
{
    if (_mode == READ_MODE)
    {
        throw FileModeException("Write attempt in read-only file",
          "Serializer::GetDataBufferAtIndex");
    }

    if (index >= _indices.size())
    {
        throw out_of_range("Invalid index in Serializer::"\
          "GetDataBufferAtIndex");
    }

    UInt32 currBlock = 0;
    if (_indices[index].blockNumber == 0)
        currBlock = _indices[_indices.size() - 1].blockNumber;
    else
        currBlock = _indices[index].blockNumber;

    UInt32 currOffset = _indices[index].offset;

    if (currBlock != _currentBlock)
    {
        WriteBlock(_currentBlock);

        _currentBlock = currBlock;
        _currentOffset = currOffset;

        UInt32 bytesRead = ReadBlock(_currentBlock);
        if (bytesRead != BLKSIZE)
        {
            throw FileException("Reading less than block size",
              "Serializer::GetDataBufferAtIndex");
        }
    }
}

void Serializer::OpenFileIO(const string& fileName, const eFileMode fileMode)
{
    int openMode = O_RDONLY;

    switch (fileMode)
    {
        case READ_MODE:
            openMode = O_RDONLY;
            break;
        case CREATE_MODE:
            openMode = O_RDWR | O_CREAT | O_TRUNC;
            break;
        case UPDATE_MODE:
            openMode = O_RDWR | O_CREAT;
            break;
        default:
            openMode = O_RDONLY;
            break;
    }

    if ((_fd = open(fileName.c_str(), openMode, S_IRUSR|S_IWUSR)) < 0)
    {
        throw FileException("Could not open file: " + fileName,
          "BlockIO::OpenFile");
    }

    SInt32 lseekoff = lseek(_fd, 0L, SEEK_END);

    if (lseekoff == -1)
    {
        throw FileException("Could not seek to the end of the file: " +
          fileName, "BlockIO::OpenFile");
    }

    _numBlocksIO = lseekoff / BLKSIZE + 1;

    if ((lseekoff % BLKSIZE) == 0)
        _numBlocksIO--;

    if (lseek(_fd, 0L, SEEK_SET) == -1)
    {
        throw FileException("Could not seek in the file: " + fileName,
          "BlockIO::OpenFile");
    }

    _currentBlockIO = 0;
}

void Serializer::CloseFileIO()
{
    close(_fd);
}

unsigned int Serializer::ReadBlock(const UInt32 blockNum)
{
    if ((blockNum != 0) && (_currentBlockIO == blockNum))
        return BLKSIZE;

    if (blockNum >= _numBlocksIO)
    {
        // Wrong block number.
        throw out_of_range("Invalid block number in BlockIO::ReadBlock");
    }

    _currentBlockIO = blockNum;

    return(_theBlock.ReadBlock(_fd, blockNum));
}

unsigned int Serializer::WriteBlock(const UInt32 blockNum)
{
    if (blockNum >= _numBlocksIO)
        _numBlocksIO = blockNum + 1;

    _currentBlockIO = blockNum;

    return(_theBlock.WriteBlock(_fd, blockNum));
}

void Serializer::AllocateIndices(const UInt32 index)
{
    if (index == _indices.size())
    {
        EntryIndex tmpEntryIndex;

        tmpEntryIndex.blockNumber = 0;
        tmpEntryIndex.offset = 0;
        tmpEntryIndex.length = 0;
        tmpEntryIndex.dataType = 0;
        tmpEntryIndex.vLength = 0; 
        tmpEntryIndex.reserved[0] = 0; 
        tmpEntryIndex.reserved[1] = 0; 
        tmpEntryIndex.reserved[2] = 0; 

        _indices.push_back(tmpEntryIndex);
    }
}

void Serializer::WriteLast(const char* const where)
{
    _currentOffset = where - _buffer;
    UInt32 n = _currentOffset % UINT32_SIZE;
    if (n != 0)
        _currentOffset += (UINT32_SIZE - n);
    if (_currentOffset + UINT32_SIZE > BLKSIZE)
    {
        WriteBlock(_currentBlock++);
        _currentOffset = 0;
    }
}

char* Serializer::GetWritingPoint(const UInt32 index)
{
    if (index != _indices.size())
    {
        GetDataBufferAtIndex(index);
    }
    else
    {
        GetLastDataBuffer();
    }

    AllocateIndices(index);

    return ((char*)(_buffer + _currentOffset));
}

void Serializer::SetVirtualLength(const UInt32 index)
{
    UInt32 n = _indices[index].length % UINT32_SIZE;

    if (n == 0)
        n = UINT32_SIZE;

    _indices[index].vLength = _indices[index].length + (UINT32_SIZE - n);
}
