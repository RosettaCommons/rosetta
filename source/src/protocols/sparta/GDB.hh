// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
/*                                                                            */
/*                     ---- TALOS+ C++ verison ----                           */
/*           TALOS: Torsion Angle Likeness Optimized By Shifts.               */
/*        Yang Shen, Gabriel Cornilescu, Frank Delaglio, and Ad Bax           */
/*                   NIH Laboratory of Chemical Physics                       */
/*                         version, x.2010.0607.00                            */
/*                                                                            */
/*                        for any problem, please contact                     */
/*                           shenyang@niddk.nih.gov                           */
/*                                                                            */
/******************************************************************************/

#ifndef GDB_H
#define GDB_H

// AUTO-REMOVED #include <string>
#include <vector>
// AUTO-REMOVED #include <map>
#include <boost/unordered_map.hpp>
#include <utility/vector0.hh>
#include <utility/vector1.hh>
#include <map>

namespace protocols {
namespace sparta {




class GDB
{
  std::string ClassType;
  bool plain_text; //only print the data matrix
  char buf[5000];

public:

  std::string GDBfileName;

  std::string FORMAT_str, VARS_str;
  int VarsNumber;
  int firstResID;
	typedef boost::unordered_map< std::string, std::string > GDB_Entry;
  GDB_Entry EMPTY;

	typedef std::map<int, std::string> ResidueList;
	ResidueList residList;

  typedef utility::vector0< std::string > StringList;
	StringList REMARKS;
  StringList DATA;
	typedef std::map<int, std::string> VarList;
	typedef std::map<int, std::string> FormatList;
	VarList VARS;
	FormatList FORMAT;

	typedef std::map< int, GDB_Entry > EntryList;
	EntryList Entries;

  GDB();
  GDB(const std::string& fileName);

  void loadGDB(const std::string &fileName);
  void saveGDB(const std::string &fileName);

  void showGDB(std::ostream& os);

  //add one VAR with given FORMAT
  void addVAR(const std::string &VAR_Name, const std::string &FORMAT_Name);
  //re-set one VAR with given FORMAT
  void setVAR(int index, const std::string &VAR_Name, const std::string &FORMAT_Name);

  //set the DATA with 'DataVal'
  void setData(const std::string &DataName, const std::string &DataVal);
  //get the DATA value with 'DataName'
  std::string getData(const std::string &DataName);


  GDB_Entry getEntry(int number);
  GDB_Entry getEntry(const std::string &VName, const std::string &VVal, int index);
  GDB_Entry getEntry(const std::string &VName1, const std::string &VVal1, const std::string &VName2, const std::string &VVal2, int index);
  //re-set the 'index'-th data entry with new 'VarVal' for given 'VarName'
  void setEntry(int index, const std::string &VarName, const std::string &VarVal);
  //add one data entry to the end of current Entries
  void addEntry(const std::string &VarName, const std::string &VarVal);


  void addRemark(const std::string &str); // add one REMARK

  std::string getResidName(int rNum);
  int getEntryCount(); // return size of current entries

  void presetClass(const std::string &ClassName); // pre-set the VARS and FORMAT

  bool checkFormat(const std::string& f); // check if f is a valid FORMAT

  bool isVarFloat(int index);
  bool isVarInt(int index);
  bool isVarString(int index);

  bool isVarFloat(const std::string &VarName);
  bool isVarInt(const std::string &VarName);
  bool isVarString(const std::string &VarName);

  void VARS_str_parser(const std::string &str);
  void FORMAT_str_parser(const std::string &str);

  void set_plaintext();


};

}
}
#endif
