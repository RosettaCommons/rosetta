// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/Field.cc
///
/// @brief Each line of a PDB file is a Record which is divided into Fields
/// @author Matthew O'Meara (mattjomeara@gmail.com)

#include <core/io/pdb/Field.hh>

// Utitility headers
#include <utility/tools/make_map.hh>


namespace core {
namespace io {
namespace pdb {


using std::string;
using std::ostream;

/// @brief various constructors - only for convinience.
Field::Field() :
  type(""),
  value(""),
  start(0),
  end(0)
{}

Field::Field(Size s, Size e) :
  type(""),
  value(""),
  start(s),
  end(e)
{}

Field::Field(string t, Size s, Size e) :
  type(t),
  value(""),
  start(s),
  end(e)
{}

/// @brief read field value from given string.
void
Field::getValueFrom(string source) {
  value = string(source.begin()+start-1, source.begin()+end);
}


/// @brief Debug output.
std::ostream&
operator <<(
  std::ostream &os,
  Field const & F
) {
  os << "[" << F.start << ", " <<  F.end << "]=" << F.value << "";
  return os;
}

/// @details static holder for collection of Fields.
RecordRef Field::PDB_Records_;

/// @details check if records table was init, init table otherwise
/// return reference to private static records collection
RecordRef & Field::getRecordCollection()
{
	if( PDB_Records_.size() == 0 ) {
		PDB_Records_ = utility::tools::make_map<string, Record>(
			"MODEL ", utility::tools::make_map<string, Field>(
				"type",			Field( 1,    6),
				"serial",		Field( 7,   80) ),

			"HEADER", utility::tools::make_map<string, Field>(
				"type",           Field( 1,  6),
				"classification", Field(11, 50),
				"depDate",        Field(51, 59),
				"idCode",         Field(63, 66)	),

			"TITLE ", utility::tools::make_map<string, Field>(
				"type",         Field( 1,  6),
				"continuation", Field( 9, 10),
				"title",        Field(11, 70) ),

			"COMPND", utility::tools::make_map<string, Field>(
				"type",         Field( 1,  6),
				"continuation", Field( 9, 10),
				"compound",     Field(11, 70) ),

			"KEYWDS", utility::tools::make_map<string, Field>(
				"type",         Field( 1,  6),
				"continuation", Field( 9, 10),
				"keywords",     Field(11, 70) ),

			"EXPDTA", utility::tools::make_map<string, Field>(
				"type",         Field( 1,  6),
				"continuation", Field( 9, 10),
				"technique",    Field(11, 70) ),

			"REMARK", utility::tools::make_map<string, Field>(
				"type",      Field( 1,  6),
				"remarkNum", Field( 8, 10),
				"value",     Field(12, 70) ), // non standard name.

			"ATOM  ", utility::tools::make_map<string, Field>(
				"type",       Field( 1,  6),
				"serial",     Field( 7, 11), /// Integer
				"name",       Field(13, 16), /// Atom
				"altLoc",     Field(17, 17), /// Character
				"resName",    Field(18, 20), /// Residue name
				"chainID",    Field(22, 22), /// Character
				"resSeq",     Field(23, 26), /// Integer
				"iCode",      Field(27, 27), /// AChar
				"x",          Field(31, 38), /// Real(8.3)
				"y",          Field(39, 46), /// Real(8.3)
				"z",          Field(47, 54), /// Real(8.3)
				"occupancy",  Field(55, 60), /// Real(6.2)
				"tempFactor", Field(61, 66), /// Real(6.2)
				///"segID",     Field(73, 76),
				"element",    Field(77, 78), /// LString(2)
				"charge",     Field(79, 80)  /// LString(2)
				),
			"HETATM", utility::tools::make_map<string, Field>(
				"type",       Field( 1,  6),
				"serial",     Field( 7, 11), /// Integer
				"name",       Field(13, 16), /// Atom
				"altLoc",     Field(17, 17), /// Character
				"resName",    Field(18, 20), /// Residue name
				"chainID",    Field(22, 22), /// Character
				"resSeq",     Field(23, 26), /// Integer
				"iCode",      Field(27, 27), /// AChar
				"x",          Field(31, 38), /// Real(8.3)
				"y",          Field(39, 46), /// Real(8.3)
				"z",          Field(47, 54), /// Real(8.3)
				"occupancy",  Field(55, 60), /// Real(6.2)
				"tempFactor", Field(61, 66), /// Real(6.2)
				///"segID",     Field(73, 76),
				"element",    Field(77, 78), /// LString(2)
				"charge",     Field(79, 80)  /// LString(2)
				),
			"SSBOND", utility::tools::make_map<string, Field>(
				"type",     Field( 1,  6),
				"serNum",   Field( 8, 10), /// Integer
				"CYS",      Field(12, 14),
				"chainID1", Field(16, 16), /// Character
				"seqNum1",  Field(18, 21), /// Integer
				"icode1",   Field(22, 22), /// AChar
				"CYS",      Field(26, 28),
				"chainID2", Field(30, 30), /// Character
				"seqNum2",  Field(32, 35), /// Integer
				"icode2",   Field(36, 36), /// AChar
				"sym1",     Field(60, 65), /// SymOP
				"sym2",     Field(67, 72)  /// SymOP
				),
			"TER   ", utility::tools::make_map<string, Field>(
				"type",       Field( 1,  6),
				"serial",     Field( 7, 11),
				"resName",    Field(18, 20),
				"chainID",    Field(22, 22),
				"resSeq",     Field(23, 26),
				"iCode",      Field(27, 27) ),

			"UNKNOW", utility::tools::make_map<string, Field>(
				"type", Field( 1,  6),
				"info", Field( 7, 80) )
			);
	}
	return PDB_Records_;
}


/// @details Debug printing, serialazing to Tracer like object.
std::ostream&
operator <<(std::ostream &os,Record const & R) {
  for(Record::const_iterator p=R.begin(); p!=R.end(); p++ ) {
    os << "<Record>{" << p->first << ":" << p->second
    << "}\n";
  }

  return os;
}


} // namespace
} // namespace
} // namespace
