// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    core/io/pdb/record_def_io.hh
/// @brief   Database input/output function definitions for PDB record definition data.
/// @author  Labonte <JWLabonte@jhu.edu>


// Unit header
#include <core/io/pdb/record_def_io.hh>
#include <core/io/pdb/Field.hh>

// Project header
#include <core/types.hh>

// Utility header
#include <utility/io/util.hh>

// Basic header
#include <basic/Tracer.hh>

// C++ header
#include <sstream>


// Construct tracer.
static THREAD_LOCAL basic::Tracer TR( "core.io.pdb.record_def_io" );


namespace core {
namespace io {
namespace pdb {

// Return a mapping of PDB record types to record definitions.
RecordRef
read_record_definitions_from_file( std::string const & filename )
{
	using namespace std;
	using namespace utility;

	vector1< string > const lines( utility::io::get_lines_from_file_data( filename ) );
	RecordRef record_defs;

	Size const n_lines( lines.size() );
	for ( uint i( 1 ); i <= n_lines; ++i ) {
		istringstream line_word_by_word( lines[ i ] );
		string record_type;  // The outer map key is the 6-letter record type name.
		string field_name;  // The name of the record field.
		core::uint start_col, end_col;  // Field column definitions.
		string data_type_string;  // A data type as defined by the PDB.
		PDBDataType data_type;

		// The first five columns must be present.
		line_word_by_word >> record_type >> field_name >> start_col >> end_col >> data_type_string;
		data_type = get_pdb_data_type_from_string( data_type_string );
		Record record;
		while ( ! line_word_by_word.fail() ) {
			record[ field_name ] = Field( start_col, end_col, data_type );
			line_word_by_word >> field_name >> start_col >> end_col >> data_type_string;
			data_type = get_pdb_data_type_from_string( data_type_string );
		}

		record_defs[ record_type ] = record;
	}

	if ( TR.Debug.visible() ) {
		TR.Debug << "Read " << record_defs.size() << " PDB record definitions from " << filename << '.' << endl;
	}

	return record_defs;
}

}  // namespace pdb
}  // namespace io
}  // namespace core
