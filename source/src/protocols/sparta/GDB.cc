// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
/*                                                                            */
/*                           ----  SPARTA   ----                              */
/*     Shifts Prediction from Analogue of Residue type and Torsion Angle      */
/*                           Yang Shen and Ad Bax                             */
/*                    J. Biomol. NMR, 38, 289-302 (2007)                      */
/*                 NIH, NIDDK, Laboratory of Chemical Physics                 */
/*                     version, 1.00 (build 2010.0607.00)                     */
/*                                                                            */
/*                      for any problem, please contact                       */
/*                          shenyang@niddk.nih.gov                            */
/*                                                                            */
/******************************************************************************/


/* GDB.cpp: class for a simple generic database */


#include <fstream>
#include <protocols/sparta/GDB.hh>
#include <protocols/sparta/util.hh>
#include <utility/exit.hh>
// Utility headers
#include <basic/Tracer.hh>
//#include <boost/unordered_map.hpp>

#include <utility/vector0.hh>
#include <utility/vector1.hh>



static basic::Tracer tr("protocols.sparta");

namespace protocols {
namespace sparta {

using namespace std;

GDB::GDB()
{
  VarsNumber = 0;
  plain_text = false;
}


GDB::GDB(const string &fileName)
{
  GDBfileName = fileName;

  VarsNumber = 0;

  loadGDB(fileName);
  plain_text = false;
}


void GDB::loadGDB(const string &fileName)
{
  GDBfileName = fileName;

  ifstream file(fileName.c_str());
  if (! file.is_open() ){
    tr.Error << "\tCan't open file " << fileName << " for reading" << endl;
    exit(0);
  }

  int entry_count=1;
  string str;
  bool DATA_START = false;

  while (! file.eof() ) {
    getline (file, str);
    if (str.empty() || str.size() == 0) continue;

    //StringList fields = split(" ", simplifyWhiteSpace(str)); split_WhiteSpace
    StringList fields = split_WhiteSpace(str);
    if ( fields.size() == 0) continue;

    //read entries
    if ( DATA_START && VarsNumber > 0 && (int) fields.size() == VarsNumber )
      {
	for(int i = 0; i < VarsNumber; i++)
	  {
	    Entries[ entry_count ][ VARS[i] ] = fields[i];
	  }
	entry_count++;
      }

    //skip the "REMARK"
    if ( fields[0] == "REMARK" )
      {
	REMARKS.push_back( str.substr( fields[0].length()+1, str.length() ) );
      }
    else if ( fields[0] == "DATA" ) {
			string text;
			text = str.substr( fields[0].length()+1, str.length() );
			DATA.push_back( text );
		} else if ( fields[0] == "VARS" ) { //read the VARS
			string text;
			text = str.substr(fields[0].length(), str.length());
			VARS_str = simplifyWhiteSpace(text);
			VARS_str_parser(VARS_str);
		} else if ( fields[0].compare("FORMAT") == 0 ) {//read the FORMAT
			FORMAT_str = str.substr(fields[0].length(), str.length());
			FORMAT_str_parser(FORMAT_str);
			DATA_START = true;
		}
  }

  file.close();

  //re-format the SEQUENCE if there exist SEQUENCE DATA
  string seq = getData("SEQUENCE");
  string firstResS = getData("FIRST_RESID");

  firstResID = (firstResS.length() > 0)? atoi( firstResS.c_str() ):1;
  if ( seq.length() > 0 )
    {
      for(int i = 0; i < (int) seq.length(); i++)
	{
	  if ( !(seq[i]>='A' && seq[i]<'Z') && seq[i]!='c' && seq[i]!='p' && seq[i]!='?') continue;
	  residList[residList.size()+firstResID] = seq.substr(i,1); // insert one-letter amino acid code
	}
    }

}



void GDB::saveGDB(const string &fileName)
{
  ofstream out(fileName.c_str(), ios::trunc);

  if ( !out ){
    tr.Error << "\tCan't save file " << fileName.c_str() << endl;
  }
	showGDB( out );
}

void GDB::showGDB( std::ostream & out ) {

  if (!plain_text) {
		for(int i = 0; i < (int)REMARKS.size(); i++)	{
			out << "REMARK " << REMARKS[i] << endl;
		}
		out << endl;

		for(int i = 0; i < (int)DATA.size(); i++)	{
			int pos = DATA[i].find_first_of(' ');
			string name = DATA[i].substr(0,pos);

			if ( name == "SEQUENCE" ) {
	      string seq = simplifyWhiteSpace( DATA[i].substr( pos+1,DATA[i].length()-pos-1 ) );
	      int len = seq.length();
	      if (len <= 0) continue;

	      for(int i = 0; i<= len/55; i++)	{
					string temp = seq.substr(i*55, 55);
					out << "DATA SEQUENCE " << temp.c_str() << endl;
				}
	    } else out << "DATA " << DATA[i] << endl;
		}
		out << endl;

		VarList::iterator iterR;
		out << "VARS   ";
		for ( iterR = VARS.begin(); iterR != VARS.end(); iterR++ )
			out << iterR->second.c_str() << " ";
		out << endl;

		out << "FORMAT ";
		for ( iterR = FORMAT.begin(); iterR != FORMAT.end(); iterR++ )
			out << iterR->second.c_str();
		out << endl << endl;
	}

  // write the entries
  for ( EntryList::iterator it = Entries.begin(); it != Entries.end(); it++ ) {
		GDB_Entry ent = it->second;

		for(int i = 0; i < (int) VARS.size(); i++) {
			if ( contains(FORMAT[i],'s') == 1)	{
				sprintf(buf, FORMAT[i].c_str(), ent[ VARS[i] ].c_str() ) ;
				out << buf;
			}	else if ( contains(FORMAT[i],'d') == 1) {
				sprintf(buf, FORMAT[i].c_str(), atoi(ent[ VARS[i] ].c_str() )) ;
				out << buf;
			}	else if ( contains(FORMAT[i],'f') == 1) {
				sprintf(buf, FORMAT[i].c_str(), atof(ent[ VARS[i] ].c_str() )) ;
				out << buf ;
			}
		}
		out << endl;
	}

}



GDB::GDB_Entry GDB::getEntry(int number)
{
  return Entries[number];
}


// get the index-th entry with VName=VVal, default return the first satisfied entry
GDB::GDB_Entry GDB::getEntry(const string &VName, const string &VVal, int index)
{
  GDB_Entry ent;

  int count = 0;
	EntryList::iterator it;
  for ( it = Entries.begin(); it != Entries.end(); it++ )
    {
      ent = it->second;
      if ( ent[VName] == VVal ) count++;

      if ( count == index && index > 0 ) return ent;
    }

  ent.clear();

  return ent;
}


// get the index-th entry with VName1=VVal1 and VName2=VVal2, default return the first satisfied entry
GDB::GDB_Entry GDB::getEntry(const string &VName1, const string &VVal1, const string &VName2, const string &VVal2, int index)
{
  GDB_Entry ent;
  int count = 0;
  EntryList::iterator it;
  for ( it = Entries.begin(); it != Entries.end(); it++ ) {
		ent = it->second;
		if ( ent[VName1] == VVal1 && ent[VName2] == VVal2 ) count++;
		if ( count == index && index > 0 ) return ent;
	}
  ent.clear();
  return ent;
}



string GDB::getResidName(int rNum)
{
  return residList[rNum];
}



int GDB::getEntryCount() // return size of current entries
{
  return Entries.size();
}



//set value to a variable of a given entry
void GDB::setEntry(
	int index,
	const string & VarName,
	const string & VarVal
) {
	VarList::const_iterator it_V;
	for( it_V = VARS.begin(); it_V != VARS.end(); it_V++ ) {
		if ( it_V->second == VarName ) break;
	}
	if ( it_V != VARS.end() ) {
		Entries[ index ][VarName.c_str()] = VarVal.c_str();
	} else {
		string const msg( "\tInvalid variable name '" + VarName + "'\n" );
		//tr.Error << "\tInvalid varible name '" << VarName << "'" << endl;
		tr.Error << msg;
		utility_exit_with_message(msg);
	}
}



//add a new entry to the end of entries list
void GDB::addEntry(const string &VarName, const string &VarVal) {
  setEntry( Entries.size()+1, VarName, VarVal);
}



//add one VAR with given FORMAT to the end of VARS list
void GDB::addVAR(const string &VAR_Name, const string &FORMAT_Name) {
  int size = VARS.size();

  if ( !checkFormat(FORMAT_Name) )
    tr.Error << "\tBad format syntax '" << FORMAT_Name <<  "'" << endl;
  else if (contains(VAR_Name, ' ') > 0)
    tr.Error << "\tInvalid varible name '" << VAR_Name <<  "' (with space)" << endl;
  else {
    if (size == 0)
      {
	VARS[0] = VAR_Name;
	FORMAT[0] = FORMAT_Name;
	VarsNumber++;
      }
    else {
      int i;
      for(i=0; i< size; i++)
	if (VARS[i] == VAR_Name) break;

      VARS[i] = VAR_Name;
      FORMAT[i] = FORMAT_Name;

      if (i >= size)  // VAR is not exist, add the VAR and its FORMAT
	VarsNumber++;
    }
  }
}


//re-set one VAR with given FORMAT, 'index' number starts from 1 and can't larger than current size + 1
//if 'index' equals to current VARS size + 1, add new VAR to the end of VARS list
void GDB::setVAR(int index, const string &VAR_Name, const string &FORMAT_Name) {
  int size = VARS.size();
  if ( !checkFormat(FORMAT_Name) )
    tr.Error << "\tBad format syntax '" << FORMAT_Name <<  "'" << endl;
  else if (contains(VAR_Name, ' ') > 0)
    tr.Error << "\tInvalid varible name '" << VAR_Name <<  "' (with space)" << endl;
  else {
    if (index > size+1 || index < 1)
      tr.Error << "\tPlease use number 1-" << size+1 << " as index for varible '" << VAR_Name << "' with format '" << FORMAT_Name << "'" << endl ;
    else {
      VARS[index-1] = VAR_Name;
      FORMAT[index-1] = FORMAT_Name;
      if (index == size+1) VarsNumber++;
    }
  }
}



//add a new REMARK entry
void GDB::addRemark(const string &str) {
  REMARKS.push_back(str);
}



void GDB::setData(const string &DataName, const string &DataVal) {
  DATA.push_back( DataName+" "+DataVal );
}



string GDB::getData(const string &DataName) {
  string data="";
  for(int i = 0; i < (int) DATA.size(); i++) {
		int pos = DATA[i].find_first_of(' ');
		if (DATA[i].substr(0,pos) == DataName )
			data = data+ DATA[i].substr( pos+1,DATA[i].length()-pos-1 );
	}
	return data;
}



bool GDB::isVarFloat(int index) {
  if (contains(FORMAT[index],'f') == 1) return true;
  return false;
}



bool GDB::isVarFloat(const string &VarName) {
  VarList::const_iterator it_V;
  for( it_V = VARS.begin(); it_V != VARS.end(); it_V++ ) {
		if ( it_V->second == VarName ) return isVarFloat(it_V->first);
	}
  return false;
}



bool GDB::isVarInt(int index) {
  if (contains(FORMAT[index],'d') == 1) return true;
  return false;
}



bool GDB::isVarInt(const string &VarName) {
  VarList::const_iterator it_V;
  for( it_V = VARS.begin(); it_V != VARS.end(); it_V++ )  {
		if ( it_V->second == VarName ) return isVarInt(it_V->first);
	}
  return false;
}



bool GDB::isVarString(int index) {
  if (contains(FORMAT[index],'s') == 1) return true;
  return false;
}



bool GDB::isVarString(const string &VarName) {
  VarList::const_iterator it_V;

  for( it_V = VARS.begin(); it_V != VARS.end(); it_V++ )
    {
      if ( it_V->second == VarName ) return isVarString(it_V->first);
    }

  return false;
}



// check if f is a valid FORMAT
bool GDB::checkFormat(const string& f)
{
  string str = simplifyWhiteSpace(f);
  int last = str.length()-1;

  if ( contains(str, '%') != 1 || str[0] != '%') return false;
  if ( contains(str, 's') != 1 && contains(str, 'd') != 1 && contains(str, 'f') != 1 ) return false;
  if ( str[last] != 's' && str[last] != 'd' && str[last] != 'f' ) return false;

  for(int i=0; i< (int)str.length();i++)
    {
      if ( str[i] != '%' && str[i] != 's' && str[i] != 'd' && str[i] != 'f' && str[i] != '-' && str[i] != '.' && !isDigit(str[i]) )
	return false;
      else if ( (str[i] == 's' || str[i] == 'd' || str[i] == 'f') && i != last)
	return false;
    }

  return true;
}


//parse the 'VARS' string and store to VARS list
void GDB::VARS_str_parser(const string &str)
{
  VARS.clear();
  Entries.clear();

  StringList V_Fields = split(" ", str);
  VarsNumber = V_Fields.size();

  for(int i = 0; i < VarsNumber; i++)
    VARS[i] = V_Fields[i];
}


//parse the 'FORMAT' string and store to FORMAT list
void GDB::FORMAT_str_parser(const string &str)
{
  FORMAT.clear();
  Entries.clear();

  string temp = str;
  temp = " " + temp;

  if (contains(temp,'%') == VarsNumber) {
      for(int i = 0; i < VarsNumber; i++) {
				string f_str = "%" + (string) section(temp,'%',buf,i+1, i+1);

				if ( checkFormat(f_str) )
					FORMAT[i] = f_str;
				else {
					tr.Error << "\tBad format syntax '" << f_str <<  "'" << endl;
					exit(0);
				}
			}
	}
}


// pre-set the VARS and FORMAT for a specified Class
void GDB::presetClass(const string &ClassName)
{
  if (ClassName == "TALOS_SHIFT" ){
    ClassType = "TALOS_SHIFT";

    VARS_str_parser("RESID RESNAME ATOMNAME SHIFT");
    FORMAT_str_parser("%4d %1s %4s  %8.3f");
  }
  else if (ClassName == "TALOS_PRED" ){
    ClassType = "TALOS_PRED";

    VARS_str_parser("INDEX PHI PSI DIST W R1 R2 R3 SOURCE");
    FORMAT_str_parser("%2d %9.3f %9.3f %8.3f %5.3f %-4s %-4s %-4s %s");
  }
}


void GDB::set_plaintext()
{
  plain_text = true;
}



}
}
