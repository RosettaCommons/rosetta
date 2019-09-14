// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/indexed_structure_store/FragmentStore.py.hh
/// @brief Python bindings for a database of residue-level fragments.
/// @author Bobby Langan (robert.langan@gmail.com)
//
#pragma once

#include "stdint.h"
#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <numeric/types.hh>

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/stl_bind.h"

#include <protocols/indexed_structure_store/FragmentStore.hh>

namespace protocols { namespace indexed_structure_store {

template<typename Module> void bind_FragmentStore(Module & m)
{

	namespace py = pybind11;

	// These lines correctly bind the member variables of FragmentStore to
	// python bindings so that they can be accessed/changed in the
	// H5PyFragmentStoreProvider
	py::class_< FragmentStore, std::shared_ptr<FragmentStore> >(m, "FragmentStore")
		.def(py::init<>())
		.def_readwrite("int64_groups", &FragmentStore::int64_groups)
		.def_readwrite("real_groups", &FragmentStore::real_groups)
		.def_readwrite("realVector_groups", &FragmentStore::realVector_groups)
		.def_readwrite("string_groups", &FragmentStore::string_groups);
}

}}
