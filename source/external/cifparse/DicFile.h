/*
FILE:     DicFile.h
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
** \file DicFile.h
**
** \brief Header file for DicFile class.
*/


#ifndef DICFILE_H
#define DICFILE_H


#include <string>
#include <iostream>

#include "GenString.h"
#include "ISTable.h"
#include "CifFile.h"


/**
**  \class DicFile
**
**  \brief Public class that represents a dictionary file, composed of
**    blocks with tables.
**
**  This class represents a dictionary file. In addition to inherited methods
**  from \e CifFile class, this class provides a method for writing the
**  content of "item_aliases" table to a text file.
*/
class DicFile : public CifFile
{
  public:
    using CifFile::Write;

    /**
    **  Constructs a dictionary file.
    **
    **  \param[in] fileMode - dictionary file mode. Possible values are
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
    DicFile(const eFileMode fileMode, const std::string& objFileName,
      const bool verbose = false, const Char::eCompareType
      caseSense = Char::eCASE_SENSITIVE,
      const unsigned int maxLineLength = STD_CIF_LINE_LENGTH,
      const std::string& nullValue = CifString::UnknownValue);

    /**
    **  Constructs a dictionary file in virtual mode.
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
    DicFile(const bool verbose = false, const Char::eCompareType
      caseSense = Char::eCASE_SENSITIVE,
      const unsigned int maxLineLength = STD_CIF_LINE_LENGTH,
      const std::string& nullValue = CifString::UnknownValue);

    /**
    **  Destructs a dictionary file, by releasing all consumed resources.
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
    ~DicFile();

    /**
    **  Writes the content of "item_aliases" table to a text file.
    **
    **  \param[in] fileName - relative or absolute name of the text file
    **    to which the content of "item_aliases" table is to be written to.
    **
    **  \return None
    **
    **  \pre None
    **
    **  \post None
    **
    **  \exception: None
    */
    void WriteItemAliases(const std::string& fileName);


    /**
    **  Method, not currently part of users public API, and will soon be
    **  re-examined.
    */
    ISTable* GetFormatTable();
 
    /**
    **  Method, not currently part of users public API, and will soon be
    **  re-examined.
    */
    int WriteFormatted(const std::string& cifFileName, ISTable* formatP = NULL);

    /**
    **  Method, not currently part of users public API, and will soon be
    **  re-examined.
    */
    int WriteFormatted(const std::string& cifFileName, TableFile* ddl,
      ISTable* formatP = NULL);

    /**
    **  Method, not currently part of users public API, and will soon be
    **  re-examined.
    */
    void Compress(CifFile* ddl);

    CifFile* GetRefFile();

  protected:
    ISTable* _formatP;

    int WriteFormatted(std::ostream& cifo, ISTable* formatP);
    int WriteFormatted(std::ostream& cifo, TableFile* ddl, ISTable* formatP);

    void WriteItemAliases(std::ostream& cifo);

  private:
    void AddRefRow(ISTable& table, const char* first, const char* second, 
      const char* third);
};

#endif
