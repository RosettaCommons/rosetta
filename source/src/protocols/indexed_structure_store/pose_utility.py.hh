// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @author Alex Ford (fordas@uw.edu)
//

#pragma once
#include "stdint.h"
#include <iostream>

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"

#include "ndarray.h"
#include "ndarray/pybind11.h"

#include <core/conformation/Residue.hh>
#include <core/pose/Pose.hh>

#include <protocols/indexed_structure_store/pose_utility.hh>
#include <protocols/indexed_structure_store/Datatypes.py.hh>
#include <protocols/indexed_structure_store/StructureStore.hh>

namespace protocols { namespace indexed_structure_store {

core::pose::PoseOP
pose_from_store(StructureStore & store, uint32_t structure_id) {
	StructureEntry s = store.get_structure(structure_id);
	if ( s.id != structure_id ) {
		throw std::invalid_argument("structure_id not found in store.");
	}

	return residue_entries_to_pose(store.get_residues(structure_id));
}


template<typename Module> void bind_pose_utility(Module & m){

	namespace py = pybind11;

	m.def("extract_residue_entry", extract_residue_entry,
		"Extract residue entry from residue.",
		py::arg("res")
	);

	m.def("extract_residue_entries", &extract_residue_entries,
		"Extract residue entry array from pose.",
		py::arg("pose"),
		py::arg("ignore_non_protein") = false
	);

	m.def("apply_residue_entries_to_pose", &apply_residue_entries_to_pose,
		"Apply components residue entry array to pose conformation.",
		py::arg("residue_entries"),
		py::arg("pose"),
		py::arg("start_residue") = 1,
		py::arg("apply_bb") = true,
		py::arg("apply_sidechain") = true,
		py::arg("apply_orient") = true
	);

	m.def("residue_entries_to_pose", &residue_entries_to_pose<ndarray::Array<ResidueEntry, 1>>,
		"Create pose from residue entry array.",
		py::arg("residue_entries"),
		py::arg("residue_type") = "fa_standard",
		py::arg("auto_termini") = true
	);

	m.def("initial_pose_for_residues", &initial_pose_for_residues,
		"Create initial pose from residue array sequence.",
		py::arg("residue_entries"),
		py::arg("residue_type") = "fa_standard",
		py::arg("auto_termini") = true
	);

	m.def("pose_from_store", &pose_from_store,
		"Get pose for given structure in store.",
		py::arg("store"),
		py::arg("structure_id")
	);
}

} }
