// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   protocols/multistate_design/util.cc
/// @brief  collection of useful utilities for setting up and using a GeneticAlgorithm with multistate design
/// @author Justin Ashworth wrote this, but Andrew Leaver-Fay copy-and-pasted it into this file.

#ifndef INCLUDED_protocols_multistate_design_util_hh
#define INCLUDED_protocols_multistate_design_util_hh

// AUTO-REMOVED #include <protocols/genetic_algorithm/Entity.hh>
#include <core/pack/task/PackerTask.fwd.hh>

#include <core/types.hh>
#include <protocols/genetic_algorithm/Entity.fwd.hh>
#include <utility/vector1.hh>


namespace protocols {
namespace multistate_design {

protocols::genetic_algorithm::EntityElements
list_amino_acid_options(
	core::Size residue_index,
	core::pack::task::ResidueLevelTask const & rtask
);

/// @breif Creates a set of PosType entity-elements from a string of 1-letter AA codes
protocols::genetic_algorithm::EntityElements
entity_elements_from_1letterstring(
	std::string const & input
);


}
}

#endif

