// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @author Alex Ford (fordas@uw.edu)
//
#pragma once

#include <boost/range.hpp>
#include <ndarray.h>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <core/conformation/Conformation.fwd.hh>
#include <core/conformation/Residue.fwd.hh>
#include <protocols/indexed_structure_store/Datatypes.hh>

namespace protocols { namespace indexed_structure_store {

ResidueEntry
extract_residue_entry(
	core::conformation::Residue const & res);

ndarray::Array<ResidueEntry, 1>
extract_residue_entries(
	core::pose::Pose const & pose,
	bool ignore_non_protein = false
);

void
apply_residue_entries_to_pose(
	ndarray::Array<ResidueEntry, 1> residue_entries,
	core::pose::Pose & pose,
	core::Size start_residue = 1,
	bool apply_bb=true,
	bool apply_sidechain=true,
	bool apply_orient=true
);

core::pose::PoseOP
residue_entries_to_pose(
	ndarray::Array<ResidueEntry, 1> residue_entries,
	std::string residue_type = "fa_standard", bool auto_termini=true
);

template<typename ResidueRange>
core::pose::PoseOP
residue_entries_to_pose(ResidueRange residue_entries, std::string residue_type = "fa_standard", bool auto_termini=true){
	ndarray::Array<ResidueEntry, 1> residues((boost::size(residue_entries)));
	std::copy( boost::begin(residue_entries), boost::end(residue_entries), residues.begin());
	return residue_entries_to_pose(residues, residue_type, auto_termini);
}


core::pose::PoseOP
initial_pose_for_residues(
	ndarray::Array<ResidueEntry, 1> residue_entries,
	std::string residue_type = "fa_standard", bool auto_termini=true
);

class PoseUtilityPlaceholder {};

} }
