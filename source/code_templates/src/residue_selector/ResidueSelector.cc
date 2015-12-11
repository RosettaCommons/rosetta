// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pack/task/residue_selector/--class--.hh
/// @brief  The --class-- selects residues in a given proximity of set focus residues
/// @author Robert Lindner (rlindner@mpimf-heidelberg.mpg.de)

// Unit headers
#include <--path--/--class--.hh>
#include --res_sel_creator--

// Basic Headers
#include <basic/datacache/DataMap.hh>

// Package headers
#include <core/pose/selection.hh>
#include <core/conformation/Residue.hh>
#include <core/select/residue_selector/ResidueSelectorFactory.hh>

// Utility Headers
#include <utility/tag/Tag.hh>
#include <basic/Tracer.hh>

// C++ headers
#include <utility/assert.hh>

static THREAD_LOCAL basic::Tracer TR( "--namespace_dot--.--class--" );


--namespace--

using namespace core::select::residue_selector;


--class--::--class--(){

}


--class--::~--class--() {}

ResidueSubset
--class--::apply( core::pose::Pose const & pose ) const
{


}

void
--class--::parse_my_tag(
	utility::tag::TagCOP tag,
	basic::datacache::DataMap & datamap)


}

std::string --class--::get_name() const {
	return --class--::class_name();
}

std::string --class--::class_name() {
	return "--class--";
}

ResidueSelectorOP
--class--Creator::create_residue_selector() const {
	return ResidueSelectorOP( new --class-- );
}

std::string
--class--Creator::keyname() const {
	return --class--::class_name();
}


--end_namespace--



