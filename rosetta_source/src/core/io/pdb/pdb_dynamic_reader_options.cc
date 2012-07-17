// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/io/pdb/pdb_dynamic_reader_options.cc
///
/// @brief
/// @author Brian D. Weitzner brian.weitzner@gmail.com

// Unit headers
#include <core/io/pdb/pdb_dynamic_reader_options.hh>

// Basic headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

// Utility headers
#include <utility/tag/Tag.hh>

// C++ headers

namespace core {
namespace io {
namespace pdb {

PDB_DReaderOptions::PDB_DReaderOptions() { init_from_options(); }

PDB_DReaderOptions::~PDB_DReaderOptions() {}

void PDB_DReaderOptions::parse_my_tag( utility::tag::TagPtr tag )
{
	FileDataOptions::parse_my_tag( tag );
	
	set_new_chain_order(  tag->getOption< bool >( "new_chain_order", 0 ));
	set_obey_ENDMDL(  tag->getOption< bool >( "obey_ENDMDL", 0 ));
}

void PDB_DReaderOptions::init_from_options()
{
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	
	set_new_chain_order( option[ in::file::new_chain_order ]() );
	set_obey_ENDMDL( option[ in::file::obey_ENDMDL ].value() );	
}

} // namespace pdb
} // namespace io
} // namespace core
