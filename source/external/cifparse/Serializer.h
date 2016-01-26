/*
FILE:     Serializer.h
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


#ifndef SERIALIZER_H
#define SERIALIZER_H


#include <string>
#include <vector>
#include <fstream>

#include "rcsb_types.h"
#include "BlockIO.h"


const int NO_TYPE = 0; // This is reserved
const unsigned int STRINGS_TYPE = 1;
const unsigned int STRING_TYPE = 2;
const int INT_TYPE = 3;
const int LONG_TYPE = 4;
const int FLOAT_TYPE = 5;
const int DOUBLE_TYPE = 6;
const unsigned int WORD_TYPE = 7;
const unsigned int WORDS_TYPE = 8;
const unsigned int UWORD_TYPE = 9;
const unsigned int UWORDS_TYPE = 10;

const int INDEX_INCREMENT = 1024;

enum eFileMode
{
    NO_MODE = 0,
    READ_MODE,
    CREATE_MODE,
    UPDATE_MODE,
    VIRTUAL_MODE
};


class Serializer
{
  public:
    // Constructors and destructor
    Serializer(const std::string& fileName, const eFileMode fileMode,
      const bool verbose = false);
    ~Serializer();

    inline unsigned int GetNumDataIndices();

    // Read methods
    UInt32 ReadUInt32(const UInt32 index);
    void ReadUInt32s(std::vector<UInt32>& UInt32s, const UInt32 index);
    void ReadString(std::string& retString, const UInt32 index);
    void ReadStrings(std::vector<std::string>& theStrings, const UInt32 index);

    // Write methods 
    UInt32 WriteUInt32(const UInt32 theWord);
    UInt32 WriteUInt32s(const std::vector<UInt32>& theWords);
    UInt32 WriteString(const std::string& theString);
    UInt32 WriteStrings(const std::vector<std::string>& theStrings);

    // Update methods
    UInt32 UpdateUInt32(const UInt32 theWord, const UInt32 oldIndex);
    UInt32 UpdateUInt32s(const std::vector<UInt32>& theWords,
      const UInt32 oldIndex);
    UInt32 UpdateString(const std::string& theString, const UInt32 oldIndex);
    UInt32 UpdateStrings(const std::vector<std::string>& theStrings,
      const UInt32 oldIndex);

  private:
    typedef struct
    {
        // Block number of the start of the indices info
        UInt32 fileIndexBlock;

        // Number of blocks that hold the indices
        UInt32 fileIndexNumBlocks;

        // Total size in bytes of all indices
        UInt32 fileIndexLength;

        // Number of indices
        UInt32 numIndices;

        // Reserved information
        UInt32 reserved[3];

        // File version
        UInt32 version;
    } tFileHeader;

    // Represents an entry index. Entry is a value of some type that has
    // been stored in to the file. Index shows entry location (in which
    // block and offset from the start of the block), its length, dataType.
    typedef struct
    {
        UInt32 blockNumber; // Block in which data is located
        UInt32 offset;      // Data offset in the block
        UInt32 length;      // The length of the data 
        UInt32 dataType;    // Type of data
        UInt32 vLength;     // Virtual length (length adjusted for word size)
        UInt32 reserved[3];
    } EntryIndex;

    static bool _littleEndian;

    static const UInt32 _version = 1;
    static const UInt32 _indicesPerBlock = BLKSIZE / sizeof(EntryIndex);

    // An array of index entries (i.e., these are indices)
    std::string _fileName;

    // Stored in block 0, this holds the info about the index
    tFileHeader _fileHeader;

    std::vector<EntryIndex> _indices;

    std::ofstream _log;

    bool _verbose;

    UInt32 _currentBlock;  // The current block number of the current buffer
    UInt32 _currentOffset; // The offset into the current buffer

    char* _buffer;

    eFileMode _mode;

    BlockIO _theBlock; // A block for doing read/write a block at a time

    void Init();

    void WriteUInt32AtIndex(const UInt32 theWord, const UInt32 index);
    void WriteUInt32sAtIndex(const std::vector<UInt32>& Words,
      const UInt32 index);
    void WriteStringAtIndex(const std::string& theString, const UInt32 index);
    void WriteStringsAtIndex(const std::vector<std::string>& theStrings,
      const UInt32 index);

    void Delete(const UInt32 index);

    void GetLastDataBuffer(void);
    void GetDataBufferAtIndex(const UInt32 index);

    void _GetHeader(const char* where);
    void _PutHeader(char* where);

    void _GetIndex(EntryIndex& outIndex, const char* where);
    void _PutIndex(const EntryIndex& inIndex, char* where);

    UInt32 _GetUInt32(const char* where);
    void _PutUInt32(const UInt32 inWord, char* where);

    void SwapHeader(tFileHeader& out, const tFileHeader& in);
    void SwapIndex(EntryIndex& out, const EntryIndex& in);
    UInt32 SwapUInt32(const UInt32 theWord);

    void _ReadFileHeader();
    void _WriteFileHeader();

    void AllocateIndices(const UInt32 index);

    UInt32 ReadBlock(const UInt32 blockNum);
    UInt32 WriteBlock(const UInt32 blockNum);

    char* GetWritingPoint(const UInt32 index);
    void WriteLast(const char* const where);

    void SetVirtualLength(const UInt32 index);

    int _fd; // The file descriptor of the file that is opened, -1 if unopened

    UInt32 _numBlocksIO; // The number of blocks in the currently opened file
    UInt32 _currentBlockIO; // The block that is currently read into buffer

    void OpenFileIO(const std::string& filename, const eFileMode fileMode);
    void CloseFileIO();

    inline UInt32 GetCurrentBlockNumberIO() const;
    inline UInt32 GetNumBlocksIO() const;

    void PrintIndex();
    void PrintIndexPosition(const UInt32 position);
    void DumpFile();
    void PrintBuffer();
};


inline UInt32 Serializer::GetNumDataIndices()
{
    return (_indices.size());
}


inline UInt32 Serializer::GetCurrentBlockNumberIO() const
{
    return (_currentBlockIO);
}


inline UInt32 Serializer::GetNumBlocksIO() const
{
    return (_numBlocksIO);
}

#endif

