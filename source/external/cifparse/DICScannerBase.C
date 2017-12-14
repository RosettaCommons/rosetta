/*
FILE:     DICScannerBase.C
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
** \file DICScannerBase.C
**
** \brief Implementation file for DICScanner class.
*/


/* 
  PURPOSE:    DDL 2.1 compliant CIF file lexer ...
*/


#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "DICScannerBase.h"
#include "DICScannerInt.h"
#include "DICParserBase.h"
#include "DICParser.h"

extern "C" int dicparser_leng;
extern "C" char* dicparser_text;
extern "C" YYSTYPE dicparser_lval;
extern "C" FILE* dicparser_in;

#define yyleng dicparser_leng
#define yytext dicparser_text
#define yylval dicparser_lval
#define yyin dicparser_in

extern "C" void dic_yy_less(int i);

extern DICParser* DICParserP;


#if 0
int yywrap(void)
{
   return(1);
}
#endif

using std::ios;
using std::endl;

DICScanner::DICScanner() {
  Clear();
  _verbose=false;
   // VLAD: WATCH HERE
  _tBuf = new string;
  _tBuf->clear();
 if (_verbose) log << "Default constructor called" << endl;
}

void DICScanner::OpenLog(const string& logName, bool verboseLevel)
{
  _verbose = verboseLevel;
//  if (_verbose && logName) log.open(logName,ios::out|ios::trunc);
  if (!logName.empty()) log.open(logName.c_str(),ios::out|ios::trunc);
}

/*
DICScanner::DICScanner(istream *in) {
  Clear();
  _verbose=false;
  // VLAD: WATCH HERE
  _tBuf = new string(1025,512);
  _tBuf->Clear();
  yyin=in;
 if (_verbose) log << "DICScanner::DICScanner(istream *in, int verbose) constructor called" << endl;
}
*/

void DICScanner::Clear(void) {
//if (_tBuf) delete _tBuf;
  _tBuf=nullptr;
  _isText = false;
  NDBlineNo=1;
}


void DICScanner::Reset(void) {
if (_tBuf) delete _tBuf;
  _tBuf=nullptr;
  Clear();
}

int DICScanner::yylex()
{
   return(0);
}

int DICScanner::ProcessNone()
{
#if DEBUG
          log << "LS0: line "<<  NDBlineNo <<  " length " << yyleng << " yytext=" << yytext << endl;
#endif
          NDBlineNo++;
          if (_isText == true) {          /* end of text value */
             for (_i=yyleng-1; _i >= 0; _i--) {
               if ( yytext[_i] == ' ' || yytext[_i] == '\t' ||  yytext[_i] == '\n') {
                  yytext[_i]='\0';
               } else if ( yytext[_i] == ';') {
                    yytext[_i]='\0';
                    break;
               } else
                  break;
             }
             (*_tBuf)+=yytext;
//          _tBuf->InsertAt(strlen(_tBuf->Text())-1,"\0");
          _tBuf->erase(strlen(_tBuf->c_str())-1,1);
//printf("%d   \n",strlen(_tBuf->Text()),_tBuf->Text());
//cout<<strlen(_tBuf->Text());
//_tBuf->Print();
//cout<<endl;
             yylval.cBuf=(char*)_tBuf->c_str();
             _isText = false;
#if DEBUG
          log << "LS1: String[" <<  strlen(yylval.cBuf) << "] " << yylval.cBuf << endl;
#endif
             return(LSITEMVALUE_DIC);
          } else {  /* text value begins */
             _isText = true;
             for (_i=0; _i < yyleng; _i++) {
                 if (yytext[_i] == ';') {  break; }
             }
             _tBuf->clear();
             (*_tBuf) += &yytext[_i+1];

          }

          return(0);
}

void DICScanner::ProcessWhiteSpace()
{
         for (_i=0; _i < yyleng; _i++)
           if (yytext[_i] == '\n') NDBlineNo++;
    if (_isText)
        (*_tBuf) += yytext;

/*      if (_isText)
        (*_tBuf) += yytext;
       else {
         for (_i=0; _i < yyleng; _i++)
           if (yytext[_i] == '\n') NDBlineNo++;
       }*/
}

int DICScanner::ProcessData()
{
      if (_isText)
        (*_tBuf) += yytext;
      else {
        yylval.cBuf=yytext;
        return (DATABLOCK_DIC);
      }

      return(0);
}

int DICScanner::ProcessItemSaveBegin()
{

    if (_isText)
        (*_tBuf) += yytext;
    else {
#if DEBUG
   log<< "SF: line "<<NDBlineNo<<", Starting save frame "<<&yytext[5]<<endl;
#endif
       if (isSave)
          log<< "Syntax error line "<< NDBlineNo<<" with "<< &yytext[5]<<", end of save expected"<<endl;
       else {
          yylval.cBuf=yytext;
          isSave=2;
          return(SAVE_BEGIN_DIC);
       }
    }

    return(0);
}

int DICScanner::ProcessCategorySaveBegin()
{

   if (_isText)
      (*_tBuf) += yytext;
   else {
#if DEBUG
   log<< "SF: line "<<NDBlineNo<<" Starting save frame "<<&yytext[5]<<endl;
#endif
       if (isSave)
          log<< "Syntax error line "<< NDBlineNo<<" with "<< &yytext[5]<<", end of save expected"<<endl;
       else {
          yylval.cBuf=yytext;
          isSave=1;
          return(SAVE_BEGIN_DIC);
          }
       }

    return(0);
}

int DICScanner::ProcessSaveEndScanner()
{
    if (_isText)
        (*_tBuf) += yytext;
    else {
    #if DEBUG
    log<<"SF: line "<<NDBlineNo<<" Ending save frame"<<endl;
    #endif
       if (!isSave)
          log<< "Syntax error line "<< NDBlineNo<<" no open save frame "<<endl;
       else {
                    yylval.cBuf=yytext;
                    isSave=0;
          return(SAVE_END_DIC);
       }
    }

    return(0);
}

int DICScanner::ProcessLoopScanner()
{
      if (_isText)
        (*_tBuf) += yytext;
      else
        return (LOOP_DIC);

      return(0);
}

void DICScanner::ProcessStop()
{
      if (_isText)
        (*_tBuf) += yytext;
}

int DICScanner::ProcessDot()
{
      if (_isText)
        (*_tBuf) += yytext;
      else
        return (UNKNOWN_DIC);

      return(0);
}

int DICScanner::ProcessQuestion()
{
      if (_isText)
        (*_tBuf) += yytext;
      else
        return (MISSING_DIC);

      return(0);
}

void DICScanner::ProcessComment()
{
      if (_isText)
        (*_tBuf) += yytext;
/*      else
        NDBlineNo++;*/
}

int DICScanner::ProcessItemNameScanner()
{
      if (_isText) {
        (*_tBuf) += yytext;
      } else {
        /* If the beginning of text is in buffer yytext */
         for (_i=0; _i<yyleng; _i++)
           if (yytext[_i] == '_') break;
         yylval.cBuf=yytext;
         return(ITEMNAME_DIC);
      }

      return(0);
}

int DICScanner::ProcessUnquotedString()
{
     if (!_isText) {
        _j=0;
        /*
        ** The string is not part of a multiline CIF value, but a CIF value
        ** anywhere else, in a loop or out of the loop.
        */
/*        for (_i=yyleng-1; _i >= 0; _i--) {
            if ( yytext[_i] == '\'' || yytext[_i] == '\"') {
               yytext[_i]='\0';
               break;
            } else
             break;
        }
        for (_i=0; _i < yyleng; _i++) {
             if (yytext[_i] == '\'' || yytext[_i] == '\"') {
                _j++;
                break;
            } else
               break;
        }
*/
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
            }
        }
#endif // REPORT_EMBEDDED_QUOTES

        return(ITEMVALUE_DIC);
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

     return(0);
}

int DICScanner::ProcessSQuotedString()
{
char * p;
     if (!_isText) {
        p=yytext;
                  p++;
        while ((p=strchr(p,'\''))) {
          p++;
          if ( p[0] == ' ' || p[0] == '\t' || p[0] == '\n') {
             _i=yyleng-strlen(p);
             dic_yy_less(_i);
             p=&yytext[yyleng];
          }
        }
        yylval.cBuf=&yytext[1];
        yylval.cBuf[_i-2]='\0';
#if DEBUG
//        ndb_log_message_text(NDB_MSG_DEBUG,"SQ String:%s",cifp_lval.TempBuffer);
#endif
        return(ITEMVALUE_DIC);
     }
     else {
        if (yytext[yyleng-1] == '\n') dic_yy_less(yyleng-1);
        (*_tBuf) += yytext;

     }

     return(0);
}

int DICScanner::ProcessDQuotedString()
{
char * p;
     if (!_isText) {
        p=yytext;
                  p++;
        while ((p=strchr(p,'\"'))) {
          p++;
          if ( p[0] == ' ' || p[0] == '\t' || p[0] == '\n') {
             _i=yyleng-strlen(p);
             dic_yy_less(_i);
             p=&yytext[yyleng];
          }
        }
        yylval.cBuf=&yytext[1];
        yylval.cBuf[_i-2]='\0';
#if DEBUG
//        ndb_log_message_text(NDB_MSG_DEBUG,"SQ String:%s",cifp_lval.TempBuffer);
#endif
        return(ITEMVALUE_DIC);
     }
     else {
        if (yytext[yyleng-1] == '\n') dic_yy_less(yyleng-1);
        (*_tBuf) += yytext;

     }

     return(0);
}

int DICScanner::ProcessEof()
{
   if (_isText == true) {
      _isText=false;
           log<<"String is not not finish at line "<<NDBlineNo<< endl;
           return(1);
        }
        else
           return(0);
}


int ProcessNoneFromDICScanner()
{
    return (DICParserP->ProcessNone());
}

void ProcessWhiteSpaceFromDICScanner()
{
    DICParserP->ProcessWhiteSpace();
}

int ProcessDataFromDICScanner()
{
    return (DICParserP->ProcessData());
}

int ProcessItemSaveBeginFromDICScanner()
{
    return (DICParserP->ProcessItemSaveBegin());
}

int ProcessCategorySaveBeginFromDICScanner()
{
    return (DICParserP->ProcessCategorySaveBegin());
}

int ProcessSaveEndFromDICScanner()
{
    return (DICParserP->ProcessSaveEndScanner());
}

int ProcessLoopFromDICScanner()
{
    return (DICParserP->ProcessLoopScanner());
}

void ProcessStopFromDICScanner()
{
    DICParserP->ProcessStop();
}

int ProcessDotFromDICScanner()
{
    return (DICParserP->ProcessDot());
}

int ProcessQuestionFromDICScanner()
{
    return (DICParserP->ProcessQuestion());
}

void ProcessCommentFromDICScanner()
{
    DICParserP->ProcessComment();
}

int ProcessItemNameFromDICScanner()
{
    return (DICParserP->ProcessItemNameScanner());
}

int ProcessUnquotedStringFromDICScanner()
{
    return (DICParserP->ProcessUnquotedString());
}

int ProcessSQuotedStringFromDICScanner()
{
    return (DICParserP->ProcessSQuotedString());
}

int ProcessDQuotedStringFromDICScanner()
{
    return (DICParserP->ProcessDQuotedString());
}

int ProcessEofFromDICScanner()
{
    return (DICParserP->ProcessEof());
}

