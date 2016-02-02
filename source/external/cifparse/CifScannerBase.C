/*
FILE:     CifScannerBase.C
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
** \file CifScannerBase.C
**
** \brief Implementation file for CifScanner class.
*/


/* 
  PURPOSE:    A DDL 2.1 compliant CIF file parser.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CifScannerBase.h"
#include "CifScannerInt.h"
#include "CifParserBase.h"
#include "CifParser.h"

// WIN32 - all extern as "C"

extern "C" int cifparser_leng;
extern "C" char* cifparser_text;
extern "C" YYSTYPE cifparser_lval;
extern "C" FILE* cifparser_in;

#define yyleng cifparser_leng
#define yytext cifparser_text
#define yylval cifparser_lval
#define yyin cifparser_in

// WIN32 - extern C
extern "C" void cif_yy_less(int i);

extern CifParser* CifParserP;

using std::ios;
using std::endl;

CifScanner::CifScanner() {
  Clear();
  _verbose=false;
  // VLAD: WATCH HERE
  _tBuf = new string;
  _tBuf->clear();
 if (_verbose) log << "Default constructor called" << endl;
}

void CifScanner::OpenLog(const string& logName, bool verboseLevel)
{
  _verbose = verboseLevel;
  if (!logName.empty()) log.open(logName.c_str(),ios::out|ios::trunc);
}

/*
CifScanner::CifScanner(istream *in) {
  Clear();
  _verbose=false;
  _tBuf = new CifString(1025,512);
  _tBuf->Clear();
  yyin=in;
 if (_verbose) log << "CifScanner::CifScanner(istream *in, int verbose) constructor called" << endl;
}
*/

void CifScanner::Clear(void) {
  _tBuf=NULL;
  _isText = false;
  errorLog.clear();
  NDBlineNo=1;
}


void CifScanner::Reset(void) {
if (_tBuf) delete _tBuf;
  _tBuf=NULL;
  Clear();
}

int CifScanner::yylex()
{
   return(0);
}


int CifScanner::ProcessNone()
{

#if DEBUG
          log << "LS0: line "<<  NDBlineNo <<  " length " << yyleng << " yytext=" << yytext << endl;
#endif  
          NDBlineNo++;
          if (_isText == true) {          /* end of text value */
#ifdef VLAD_DEBUG
             cout << "DEBUG - In CifScanner::ProcessNone() - end of multi-line text \"" << yytext << "\"" << endl;
#endif

             // Check if the first character is semicolon followed by a
             // non-newline. This is considered invalid.
             if ((yyleng > 1) && ((yytext[0] == ';') && (yytext[1] != '\n')))
             {
                 log << "ERROR - Invalid syntax. Improperly placed semicolon in line " << NDBlineNo - 1 <<  endl;

                 errorLog += "ERROR - Invalid syntax. Improperly placed semicolon in line ";
                 errorLog += String::IntToString(NDBlineNo - 1);
                 errorLog += '\n';
             }

             for (_i=yyleng-1; _i >= 0; _i--) {
               if ( yytext[_i] == ' ' || yytext[_i] == '\t' ||  yytext[_i] == '\n' || yytext[_i] == '\r') {
                  yytext[_i]='\0';
               } else if ( yytext[_i] == ';') {
                    yytext[_i]='\0';
                    break;
               } else
                  break;
             }
             (*_tBuf)+=yytext;
          _tBuf->erase(strlen(_tBuf->c_str())-1,1);
             yylval.cBuf=(char*)_tBuf->c_str();
             _isText = false;
#if DEBUG
          log << "LS1: String[" <<  strlen(yylval.cBuf) << "] " << yylval.cBuf << endl;
#endif
             return(LSITEMVALUE_CIF);
          } else {  /* text value begins */
#ifdef VLAD_DEBUG
             cout << "DEBUG - In CifScanner::ProcessNone() - begin of multi-line text \"" << yytext << "\"" << endl;
#endif
             _isText = true;

             for (_i=0; _i < yyleng; _i++) {
                 if (yytext[_i] == ';') {  break; }
             }
             _tBuf->clear();
             string tmpP;
             for (unsigned int tmpI = 0; tmpI < strlen(&yytext[_i+1]); tmpI++)
             {
                 if (yytext[_i+1+tmpI] != '\r')
                 {
                     tmpP.push_back(yytext[_i+1+tmpI]);
                 }
             }
       
             (*_tBuf) += tmpP;
          }
    return (0);
}


void CifScanner::ProcessWhiteSpace()
{
         for (_i=0; _i < yyleng; _i++)
           if (yytext[_i] == '\n') NDBlineNo++;
    if (_isText)
        (*_tBuf) += yytext;

}

int CifScanner::ProcessData()
{
      if (_isText)
        (*_tBuf) += yytext;
      else {
        yylval.cBuf=yytext;
        return (DATABLOCK_CIF);
      }

      return (0);
}

int CifScanner::ProcessLoopScanner()
{
      if (_isText)
        (*_tBuf) += yytext;
      else
        return (LOOP_CIF);

      return (0);
}

void CifScanner::ProcessStop()
{
      if (_isText)
        (*_tBuf) += yytext;
}

int CifScanner::ProcessDot()
{
      if (_isText)
        (*_tBuf) += yytext;
      else
        return (UNKNOWN_CIF);

      return (0);
}

int CifScanner::ProcessQuestion()
{
      if (_isText)
        (*_tBuf) += yytext;
      else
        return (MISSING_CIF);

      return (0);
}

void CifScanner::ProcessComment()
{
      if (_isText)
        (*_tBuf) += yytext;
}

int CifScanner::ProcessUnderscore()
{
      if (_isText) {
        (*_tBuf) += yytext;
      } else {
        /* If the beginning of text is in buffer yytext */
         for (_i=0; _i<yyleng; _i++)
           if (yytext[_i] == '_') break;
         yylval.cBuf=yytext;
         return(ITEMNAME_CIF);
      }

      return (0);
}

int CifScanner::ProcessBadStrings()
{
     if (!_isText) {
        /*
        ** The string is not part of a multiline CIF value, but a CIF value
        ** anywhere else, in a loop or out of the loop.
        */
        _j=0;
        yylval.cBuf=&yytext[_j];
#if DEBUG
        log << "UQ: String " << yylval.cBuf << endl;
#endif

#ifdef REPORT_EMBEDDED_QUOTES
        unsigned int cBufLen = strlen(yylval.cBuf);
        for (unsigned int i = 0; i < cBufLen; ++i)
        {
            if ((yylval.cBuf[i] == '\'') || (yylval.cBuf[i] == '\"'))
            {
                log << "ERROR - Invalid character at line " <<
                  String::IntToString(NDBlineNo) << " in CIF value " <<
                  yylval.cBuf << endl;

                errorLog += "ERROR - Invalid character at line ";
                errorLog += String::IntToString(NDBlineNo);
                errorLog += " in CIF value ";
                errorLog += yylval.cBuf;
                errorLog += '\n';
            }
        }
#endif // REPORT_EMBEDDED_QUOTES

        return(ITEMVALUE_CIF);
     }
     else {
        /*
        ** The string is part of a multiline CIF value. It is processed as is.
        */
#if DEBUG
          log << "UQx: String " << yytext<< endl;
#endif
        (*_tBuf) += yytext;
     }

     return (0);
}

int CifScanner::ProcessSQuotedStrings()
{
char * p;
     if (!_isText) {
        p=yytext;
                  p++;
        while ((p=strchr(p,'\''))) {
          p++;
          if ( p[0] == ' ' || p[0] == '\t' || p[0] == '\n' || p[0] == '\r') {
             _i=yyleng-strlen(p);
             cif_yy_less(_i);
             p=&yytext[yyleng];
          }
        }
        yylval.cBuf=&yytext[1];
        yylval.cBuf[_i-2]='\0';
        return(ITEMVALUE_CIF);
     }
     else {
        if (yytext[yyleng-1] == '\n')
        {
            if (yytext[yyleng-2] == '\r')
            {
               cif_yy_less(yyleng-2);
            }
            else
            {
               cif_yy_less(yyleng-1);
            }
        }

        (*_tBuf) += yytext;
     }

     return (0);
}

int CifScanner::ProcessDQuotedStrings()
{
char * p;
     if (!_isText) {
        p=yytext;
                  p++;
        while ((p=strchr(p,'\"'))) {
          p++;
          if ( p[0] == ' ' || p[0] == '\t' || p[0] == '\n' || p[0] == '\r') {
             _i=yyleng-strlen(p);
             cif_yy_less(_i);
             p=&yytext[yyleng];
          }
        }
        yylval.cBuf=&yytext[1];
        yylval.cBuf[_i-2]='\0';
        return(ITEMVALUE_CIF);
     }
     else {
        if (yytext[yyleng-1] == '\n')
        {
            if (yytext[yyleng-2] == '\r')
            {
               cif_yy_less(yyleng-2);
            }
            else
            {
               cif_yy_less(yyleng-1);
            }
        }

        (*_tBuf) += yytext;
     }

     return (0);
}

int CifScanner::ProcessEof()
{
   if (_isText == true) {
      _isText=false;

          errorLog += "ERROR - String is not finished above line ";
          errorLog += String::IntToString(NDBlineNo);
          errorLog += '\n';

           log<<"ERROR - String is not finished above line "<<NDBlineNo<< endl;
           return(1);
        }
        else
           return(0);
}

int ProcessNoneFromScanner()
{
    return (CifParserP->ProcessNone());
}

void ProcessWhiteSpaceFromScanner()
{
    CifParserP->ProcessWhiteSpace();
}

int ProcessDataFromScanner()
{
    return (CifParserP->ProcessData());
}

int ProcessLoopFromScanner()
{
    return (CifParserP->ProcessLoopScanner());
}

void ProcessStopFromScanner()
{
    CifParserP->ProcessStop();
}

int ProcessDotFromScanner()
{
    return (CifParserP->ProcessDot());
}


int ProcessQuestionFromScanner()
{
    return (CifParserP->ProcessQuestion());
}


void ProcessCommentFromScanner()
{
    CifParserP->ProcessComment();
}


int ProcessUnderscoreFromScanner()
{
    return (CifParserP->ProcessUnderscore());
}


int ProcessBadStringsFromScanner()
{
    return (CifParserP->ProcessBadStrings());
}


int ProcessSQuotedStringsFromScanner()
{
    return (CifParserP->ProcessSQuotedStrings());
}


int ProcessDQuotedStringsFromScanner()
{
    return (CifParserP->ProcessDQuotedStrings());
}


int ProcessEofFromScanner()
{
    return (CifParserP->ProcessEof());
}

