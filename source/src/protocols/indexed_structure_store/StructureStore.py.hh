// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/indexed_structure_store/StructureStore.py.hh
/// @brief Python bindings for a POD-based database of residue-level protein structures.
/// @author Alex Ford (fordas@uw.edu)
//
#pragma once

#include "stdint.h"
#include <iostream>

#include "pybind11/pybind11.h"
#include "pybind11/numpy.h"
#include "pybind11/eval.h"

#include "ndarray.h"
#include "ndarray/pybind11.h"

#include <protocols/indexed_structure_store/StructureStore.hh>
#include <protocols/indexed_structure_store/Datatypes.hh>
#include <protocols/indexed_structure_store/Datatypes.py.hh>

namespace protocols { namespace indexed_structure_store {

template<typename Module> void bind_StructureStore(Module & m)
{
	namespace py = pybind11;

	// Guard by numpy import, aborting binding definition if numpy interface
	// isn't available or numpy version < 1.9.
	try {
		auto numpy = py::module::import("numpy");
		std::string version_check = "[int(x) for x in __version__.split('.')[:2]] < [int(x) for x in '1.9'.split('.')]";
		if ( py::eval(version_check, numpy.attr("__dict__")).cast<bool>() ) {
			return;
		}
	} catch (...) {
		return;
	}



	PYBIND11_NUMPY_DTYPE(StructureEntry, id, name);
	PYBIND11_NUMPY_DTYPE(ResidueBackboneEntry, phi, psi, omega);

	::pybind11::detail::npy_format_descriptor<ResidueSidechainEntry>::register_dtype({
		PYBIND11_FIELD_DESCRIPTOR(ResidueSidechainEntry, chi1),
		PYBIND11_FIELD_DESCRIPTOR(ResidueSidechainEntry, chi2),
		PYBIND11_FIELD_DESCRIPTOR(ResidueSidechainEntry, chi3),
		PYBIND11_FIELD_DESCRIPTOR(ResidueSidechainEntry, chi4),
		::pybind11::detail::field_descriptor{
		"aa",
		offsetof(ResidueSidechainEntry, aa),
		sizeof(decltype(std::declval<ResidueSidechainEntry>().aa)),
		::pybind11::format_descriptor<char[1]>::format(),
		::pybind11::detail::npy_format_descriptor<char[1]>::dtype()
		}
		});

	PYBIND11_NUMPY_DTYPE(ResidueOrientEntry, N, C, CA, O);
	PYBIND11_NUMPY_DTYPE(ResidueEntry, structure_id, residue_id, bb, sc, orient, chain_ending);

	m.attr("structure_entry") = py::dtype::of<StructureEntry>();
	m.attr("residue_entry_dtype") = py::dtype::of<ResidueEntry>();

	py::class_< StructureStore, std::shared_ptr<StructureStore> >(m, "StructureStore")
		.def(py::init<>())
		.def("get_structure", &StructureStore::get_structure)
		.def("get_residues", &StructureStore::get_residues)
		.def("check_entries", &StructureStore::check_entries)
		.def("resize", &StructureStore::resize)
		.def_readwrite("structure_entries", &StructureStore::structure_entries)
		.def_readwrite("residue_entries", &StructureStore::residue_entries);
}

}}
