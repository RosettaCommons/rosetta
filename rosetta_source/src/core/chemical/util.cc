// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file core/chemical/util.cc
/// @brief Utilities for modifying and utilizing Residues and other core::chemical classes.


// Unit headers
#include <core/chemical/util.hh>

// Package Headers
// AUTO-REMOVED #include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/ChemicalManager.hh>

// Project Headers
#include <core/types.hh>
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/vector1.hh>

namespace core {
namespace chemical {

core::chemical::ResidueTypeSetCAP
rsd_set_from_cmd_line() {
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;

	std::string const type_set_name( option[ in::file::residue_type_set ]() );
	ResidueTypeSetCAP set = ChemicalManager::get_instance()->residue_type_set(
		type_set_name
	);

	return set;
}

} // namespace chemical
} // namespace core
