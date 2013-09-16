// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/pdb_dynamic_reader_options.hh
/// @brief  Declarations for PDB_DReaderOptions.
/// @author Brian D. Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_core_io_pdb_pdb_dynamic_reader_options_HH
#define INCLUDED_core_io_pdb_pdb_dynamic_reader_options_HH


// Unit headers
#include <core/io/pdb/pdb_dynamic_reader_options.fwd.hh>
#include <core/io/pdb/file_data_options.hh>

// C++ headers
#include <string>


namespace core {
namespace io {
namespace pdb {

class PDB_DReaderOptions : public FileDataOptions
{
public:
	PDB_DReaderOptions();
	
	virtual ~PDB_DReaderOptions();
	
	virtual
	void parse_my_tag( utility::tag::TagPtr tag );
	
	virtual
	std::string type() const;
	
	// accessors
	bool new_chain_order() const;
	bool obey_ENDMDL() const;
	bool read_pdb_header() const;
	bool read_link_records() const;
	
	// mutators
	void set_new_chain_order( bool setting );
	void set_obey_ENDMDL( bool setting );
	void set_read_pdb_header( bool setting );
	void set_read_link_records( bool setting );

private:
	/// @brief Assigns user specified values to primitive members using command line options
	void init_from_options();

private:
	bool new_chain_order_;
	bool obey_ENDMDL_;
	bool read_pdb_header_;
	bool read_link_records_;
};

} // namespace pdb
} // namespace io
} // namespace core

#endif // INCLUDED_core_io_pdb_pdb_dynamic_reader_options_HH
