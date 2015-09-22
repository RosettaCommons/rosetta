// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/Field.cc
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

/// @brief various constructors - only for convenience.
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

/// @details check if records table was init, init table otherwise
/// return reference to private static records collection
/// @remarks See http://www.wwpdb.org/docs.html#format for details.
RecordRef & Field::getRecordCollection()
{
	static RecordRef records;
	if ( records.size() == 0 ) {
		records = utility::tools::make_map<string, Record>(
			// Title Section
			"HEADER", utility::tools::make_map<string, Field>(
			"type",           Field( 1,  6),
			"classification", Field(11, 50),
			"depDate",        Field(51, 59),
			"idCode",         Field(63, 66) ),

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
			"value",     Field(12, 70) ), // non-standard name

			// Primary Structure Section

			// Heterogen Section
			"HETNAM", utility::tools::make_map<string, Field>(
			"type",         Field( 1,  6),
			"continuation", Field( 9, 10),
			"hetID",        Field(12, 14), // LString(3): HET identifier, right-justified
			"text",         Field(16, 70)  // String: chemical name
			),

			// Secondary Structure Section

			// Connectivity Annotation Section
			"SSBOND", utility::tools::make_map<string, Field>(
			"type",     Field( 1,  6),
			"serNum",   Field( 8, 10), // Integer
			"resName1", Field(12, 14),
			"chainID1", Field(16, 16), // Character
			"resSeq1",  Field(18, 21), // Integer
			"iCode1",   Field(22, 22), // AChar
			"resName2", Field(26, 28),
			"chainID2", Field(30, 30), // Character
			"resSeq2",  Field(32, 35), // Integer
			"iCode2",   Field(36, 36), // AChar
			"sym1",     Field(60, 65), // SymOP
			"sym2",     Field(67, 72), // SymOP
			"length",   Field(74, 78)  // Real(5.2) 3.3 standard
			),
			"LINK  ", utility::tools::make_map<string, Field>(
			"type",     Field( 1,  6),
			"name1",    Field(13, 16), // Atom
			"altLoc1",  Field(17, 17), // Character
			"resName1", Field(18, 20), // Residue name
			"chainID1", Field(22, 22), // Character
			"resSeq1",  Field(23, 26), // Integer
			"iCode1",   Field(27, 27), // AChar
			"name2",    Field(43, 46), // Atom
			"altLoc2",  Field(47, 47), // Character
			"resName2", Field(48, 50), // Residue name
			"chainID2", Field(52, 52), // Character
			"resSeq2",  Field(53, 56), // Integer
			"iCode2",   Field(57, 57), // AChar
			"sym1",     Field(60, 65), // SymOP
			"sym2",     Field(67, 72), // SymOP
			"length",   Field(74, 78)  // Real(5.2)
			),

			// Miscellaneous Features Section

			// Crystallographic & Coordinate Transformation Section

			// Coordinate Section
			"MODEL ", utility::tools::make_map<string, Field>(
			"type",   Field( 1,  6),
			"serial", Field( 7, 80) ),

			"CRYST1", utility::tools::make_map<string, Field>(
			"type",       Field( 1,  6),
			"a",          Field( 7, 15), // Real(9.3)
			"b",          Field(16, 24), // Real(9.3)
			"c",          Field(25, 33), // Real(9.3)
			"alpha",      Field(34, 40), // Real(7.2)
			"beta",       Field(41, 47), // Real(7.2)
			"gamma",      Field(48, 54), // Real(7.2)
			"spacegroup", Field(56, 66), // LString
			"z",          Field(67, 70)  // Integer
			),
			"ATOM  ", utility::tools::make_map<string, Field>(
			"type",       Field( 1,  6),
			"serial",     Field( 7, 11), // Integer
			"name",       Field(13, 16), // Atom
			"altLoc",     Field(17, 17), // Character
			"resName",    Field(18, 20), // Residue name
			"chainID",    Field(22, 22), // Character
			"resSeq",     Field(23, 26), // Integer
			"iCode",      Field(27, 27), // AChar
			"x",          Field(31, 38), // Real(8.3)
			"y",          Field(39, 46), // Real(8.3)
			"z",          Field(47, 54), // Real(8.3)
			"occupancy",  Field(55, 60), // Real(6.2)
			"tempFactor", Field(61, 66), // Real(6.2)
			//"segID",     Field(73, 76),
			"element",    Field(77, 78), // LString(2)
			"charge",     Field(79, 80)  // LString(2)
			),
			"HETATM", utility::tools::make_map<string, Field>(
			"type",       Field( 1,  6),
			"serial",     Field( 7, 11), // Integer
			"name",       Field(13, 16), // Atom
			"altLoc",     Field(17, 17), // Character
			"resName",    Field(18, 20), // Residue name
			"chainID",    Field(22, 22), // Character
			"resSeq",     Field(23, 26), // Integer
			"iCode",      Field(27, 27), // AChar
			"x",          Field(31, 38), // Real(8.3)
			"y",          Field(39, 46), // Real(8.3)
			"z",          Field(47, 54), // Real(8.3)
			"occupancy",  Field(55, 60), // Real(6.2)
			"tempFactor", Field(61, 66), // Real(6.2)
			//"segID",     Field(73, 76),
			"element",    Field(77, 78), // LString(2)
			"charge",     Field(79, 80)  // LString(2)
			),
			"TER   ", utility::tools::make_map<string, Field>(
			"type",    Field( 1,  6),
			"serial",  Field( 7, 11),
			"resName", Field(18, 20),
			"chainID", Field(22, 22),
			"resSeq",  Field(23, 26),
			"iCode",   Field(27, 27) ),

			// Connectivity Section

			// Bookkeeping Section

			// Unknown Type
			"UNKNOW", utility::tools::make_map<string, Field>(
			"type", Field( 1,  6),
			"info", Field( 7, 80) )
		);
	}
	return records;
}


/// @details Debug printing, serializing to Tracer like object.
std::ostream&
operator <<(std::ostream &os,Record const & R) {
	for ( Record::const_iterator p=R.begin(), end=R.end(); p!=end; ++p ) {
		os << "<Record>{" << p->first << ":" << p->second
			<< "}\n";
	}

	return os;
}


} // pdb
} // io
} // core
