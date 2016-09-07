// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/io/pdb/StructFileRepOptions.hh
/// @brief  Declarations for StructFileReaderOptions.
/// @author Brian D. Weitzner brian.weitzner@gmail.com

#ifndef INCLUDED_core_io_StructFileReaderOptions_HH
#define INCLUDED_core_io_StructFileReaderOptions_HH


// Unit headers
#include <core/io/StructFileReaderOptions.fwd.hh>
#include <core/io/StructFileRepOptions.hh>

// Basic headers
#include <utility/vector1.hh>
#include <basic/resource_manager/ResourceOptions.hh>

// C++ headers
#include <string>


namespace core {
namespace io {

class StructFileReaderOptions : public StructFileRepOptions
{
public:
	/// @brief Constructor that takes default values from the global OptionCollection object, basic::options::option.
	StructFileReaderOptions();
	/// @brief Constructor that takes default values from a provided OptionCollection object
	StructFileReaderOptions( utility::options::OptionCollection const & options );

	~StructFileReaderOptions() override;


	void parse_my_tag( utility::tag::TagCOP tag ) override;


	std::string type() const override;

	// accessors
	bool new_chain_order() const;
	bool obey_ENDMDL() const;
	bool read_pdb_header() const;
	bool glycam_pdb_format() const;

	// mutators
	void set_new_chain_order( bool setting );
	void set_obey_ENDMDL( bool setting );
	void set_read_pdb_header( bool setting );
	void set_glycam_pdb_format( bool setting );

	/// @brief Declare the list of options that are read in the process of reading a PDB (or SDF) and converting
	/// it into a Pose.
	static
	void
	list_options_read( utility::options::OptionKeyList & read_options );

private:
	/// @brief Assigns user specified values to primitive members using command line options
	void init_from_options( utility::options::OptionCollection const & options );

private:
	bool new_chain_order_;
	bool obey_ENDMDL_;
	bool read_pdb_header_;
	bool glycam_pdb_format_;
};

} // namespace io
} // namespace core

#endif // INCLUDED_core_io_pdb_StructFileRepOptions_HH
