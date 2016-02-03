// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/StructFileReaderOptions.cc
/// @brief  Definitions for StructFileReaderOptions.
/// @author Brian D. Weitzner brian.weitzner@gmail.com

// Unit headers
#include <core/io/StructFileReaderOptions.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/run.OptionKeys.gen.hh>
#include <basic/options/keys/carbohydrates.OptionKeys.gen.hh>

// Utility headers
#include <utility/tag/Tag.hh>


namespace core {
namespace io {

StructFileReaderOptions::StructFileReaderOptions() { init_from_options(); }

StructFileReaderOptions::~StructFileReaderOptions() {}

void StructFileReaderOptions::parse_my_tag( utility::tag::TagCOP tag )
{
	StructFileRepOptions::parse_my_tag( tag );

	set_new_chain_order(  tag->getOption< bool >( "new_chain_order", 0 ) );
	set_obey_ENDMDL(  tag->getOption< bool >( "obey_ENDMDL", 0 ) );
	set_read_pdb_header( tag->getOption< bool >( "preserve_header", 0 ) );
	set_read_link_records( tag->getOption< bool >( "read_pdb_link_records", 0 ) );
	set_glycam_pdb_format( tag->getOption< bool >( "glycam_pdb_format", 0 ) );
}

std::string StructFileReaderOptions::type() const { return "StructFileReaderOptions"; }

// accessors
bool StructFileReaderOptions::new_chain_order() const { return new_chain_order_; }
bool StructFileReaderOptions::obey_ENDMDL() const { return obey_ENDMDL_; }
bool StructFileReaderOptions::read_pdb_header() const { return read_pdb_header_; }
bool StructFileReaderOptions::read_link_records() const { return read_link_records_; }
bool StructFileReaderOptions::glycam_pdb_format() const { return glycam_pdb_format_; }

// mutators
void StructFileReaderOptions::set_new_chain_order( bool setting ) { new_chain_order_ = setting; }
void StructFileReaderOptions::set_obey_ENDMDL( bool setting ) { obey_ENDMDL_ = setting; }
void StructFileReaderOptions::set_read_pdb_header( bool setting ) { read_pdb_header_ = setting; }
void StructFileReaderOptions::set_read_link_records( bool setting ) { read_link_records_ = setting; }
void StructFileReaderOptions::set_glycam_pdb_format( bool setting ) { glycam_pdb_format_ = setting; }

void StructFileReaderOptions::init_from_options()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;

	set_new_chain_order( option[ in::file::new_chain_order ]() );
	set_obey_ENDMDL( option[ in::file::obey_ENDMDL ].value() );
	set_read_pdb_header( option[ run::preserve_header ]() );
	set_read_link_records( option[ in::file::read_pdb_link_records ]() );
	set_glycam_pdb_format( option[ carbohydrates::glycam_pdb_format ]() );
}

} // namespace io
} // namespace core
