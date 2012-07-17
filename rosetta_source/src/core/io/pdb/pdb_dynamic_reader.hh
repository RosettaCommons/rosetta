// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/pdb_dynamic_reader.hh
///
/// @brief
/// @author Sergey Lyskov (Sergey.Lyskov@jhu.edu)

#ifndef INCLUDED_core_io_pdb_pdb_dynamic_reader_hh
#define INCLUDED_core_io_pdb_pdb_dynamic_reader_hh


// Unit headers
#include <core/io/pdb/file_data.fwd.hh>
#include <core/io/pdb/pdb_dynamic_reader_options.fwd.hh>

// Utility headers
#include <utility/Show.hh>
// AUTO-REMOVED #include <utility/vector1.hh>

// C++ headers
//#include <cstdlib>
#ifdef WIN32
#include <string>
#endif

#include <map>
#include <vector>
#include <iostream>

#include <utility/vector1.fwd.hh>


namespace core {
namespace io {
namespace pdb {


// Forward
class Field;

/// Record type, represent one line in pdb file
typedef std::string String;
typedef std::map<String, Field> Record;

/// @brief Debug printing, serialazing to Tracer like object.
std::ostream& operator <<(std::ostream &os,Record const & R);

/// ------------------------------------------------------------------------------------------------
/// @brief Data type Class to represent one field in PDB file.
class Field : public utility::Show
{
public:
	/// @brief various constructors - only for convinience.
	Field() : type(""), value(""), start(-1), end(-1) {};
	Field(int s, int e) { start=s; end=e; };
	Field(String type_, int s, int e) { type=type_; start=s; end=e; };

	/// @brief read field value from given string.
	void getValueFrom(String source) { value = String(source.begin()+start-1, source.begin()+end); }


public:  /// This class is intended to be just 'data' type class - no need to make it private.

	/// @brief string value of field, type of the field.
	String type, value;

	/// @brief begining position in line, ending postion in line
	int start, end;


	/// @brief Debug output.
	friend std::ostream& operator <<(std::ostream &os, Field const & F) {
		//os << "type[" << F.start << ", " <<  F.end << "]="<<  F.type << "  value=" <<  F.value << "\n";
		os << "[" << F.start << ", " <<  F.end << "]="<< F.value << "";
		return os;
	};
};



/// @brief PDB Reader it self, D - for dynamic approch of type handling
class PDB_DReader
{
public:
	/// @brief creating record from given string. Also, read Field values from string.
	static Record mapStringToRecord(const String & s);


	/// @brief Reverse opearation - create PDB string from given Record
	static String createPDBString(const Record & R);


	/// @brief Parse whole PDB string and return vector of records in order they was in PDB.
	static std::vector<Record> parse(const String &);


	/// @brief create File data sturcture from array of Records. 
	static FileData createFileData(std::vector<Record> &);

	/// @brief create File data sturcture from array of Records and a set of options. 
	static FileData createFileData(std::vector<Record> &, PDB_DReaderOptions const & options);

	/// @brief create File data sturcture from string containing PDB information.
	static FileData createFileData(const String & data);

	/// @brief create File data sturcture from string containing PDB information and a set of options.
	static FileData createFileData(const String & data, PDB_DReaderOptions const & options);
	
	/// @brief create PDB-like string to represent given FileData object
	static String createPDBData(FileData const &fd);

	/// @brief create PDB-like vector of string to represent given FileData object.
	static utility::vector1< std::string > createPDBData_vector(FileData const & fd );


	/// @brief create vector of records for given FileData object.
	static std::vector<Record> createRecords(FileData const & fd);


private:
	/// @brief collection of all possible records (line types), that can exist in PDB file.
	typedef std::map<String, Record> RecordRef;

	/// @brief static holder for collection.
	static RecordRef PDB_Records_;


	/// @brief collection builder
	static RecordRef & getRecordCollection();
};

} // namespace pdb
} // namespace io
} // namespace core


#endif // INCLUDED_core_io_pdb_pdb_dynamic_reader_HH
