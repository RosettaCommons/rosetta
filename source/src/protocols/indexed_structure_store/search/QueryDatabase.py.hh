// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/indexed_structure_store/search/QueryDatabase.py.hh
/// @brief Database for generic alignment-based queries over structure stores.
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

#include <protocols/indexed_structure_store/search/QueryDatabase.hh>

#include <protocols/indexed_structure_store/Datatypes.hh>
#include <protocols/indexed_structure_store/Datatypes.py.hh>

namespace protocols { namespace indexed_structure_store { namespace search {

template<typename Module> void bind_QueryDatabase(Module & m)
{
	namespace py = pybind11;

	// Guard by numpy import, aborting binding definition if numpy interface
	// isn't available or numpy version < 1.9.
	try {
		auto numpy = py::module::import("numpy");
		std::string version_check = "__version__.split('.') < '1.9'.split('.')";
		if ( py::eval(version_check, numpy.attr("__dict__")).cast<bool>() ) {
			return;
		}
	} catch (...) {
		return;
	}

	PYBIND11_NUMPY_DTYPE(StructurePairQueryResult, fragment_a_start, fragment_b_start, result_rmsd);
	PYBIND11_NUMPY_DTYPE(StructureSingleQueryResult, fragment_start, result_rmsd);

	py::class_< StructureDatabase >(m, "StructureDatabase")
		.def(py::init<>())
		.def("initialize",
		(void (StructureDatabase::*)(
		ndarray::Array<SearchReal, 3, 3> coordinate_buffer,
		ndarray::Array<Index, 1> structure_endpoints,
		ndarray::Array<Index, 1> chain_endpoints))
		&StructureDatabase::initialize )
		.def("initialize",
		(void (StructureDatabase::*)(
		ndarray::Array<ResidueEntry, 1>))
		&StructureDatabase::initialize )
		.def("get_structure_data", &StructureDatabase::get_structure_data, py::return_value_policy::reference_internal);

	py::class_<StructureData >(m, "StructureData")
		.def_readonly("structure_offset",  &StructureData::structure_offset)
		.def_readonly("coordinate_buffer",  &StructureData::coordinate_buffer)
		.def_readonly("chain_endpoints",  &StructureData::chain_endpoints)
		.def("get_fragment_indicies",  &StructureData::get_fragment_indicies_copy)
		.def("get_fragment_centers_of_mass",  &StructureData::get_fragment_centers_of_mass_copy);

	py::class_< StructurePairQuery >(m, "StructurePairQuery")
		.def(
		py::init< ndarray::Array<SearchReal, 3, 3>, ndarray::Array<SearchReal, 3, 3>, SearchReal>())
		.def(
		py::init< ndarray::Array<SearchReal, 3, 3>, ndarray::Array<SearchReal, 3, 3>, SearchReal, int, int >())
		.def_readonly("n_entry_a",  &StructurePairQuery::n_entry_a)
		.def_readonly("n_entry_b",  &StructurePairQuery::n_entry_b)
		.def_readonly("q_buffer_a",  &StructurePairQuery::q_buffer_a )
		.def_readonly("q_buffer_b",  &StructurePairQuery::q_buffer_b )
		.def_readonly("c_per_entry",  &StructurePairQuery::c_per_entry)
		.def_readonly("rmsd_tolerance",  &StructurePairQuery::rmsd_tolerance);

	py::class_< PairQueryExecutor >(m, "PairQueryExecutor")
		.def(py::init<StructurePairQuery const &>())
		.def("execute", &PairQueryExecutor::execute )
		.def("execute_structure", &PairQueryExecutor::execute_structure )
		.def_property_readonly("query_results",
		[](PairQueryExecutor & src) {
			return py::array_t<StructurePairQueryResult>(
			src.query_results.size(), &src.query_results[0]);
		});

	py::class_<PairQuerySummaryStatistics >(m, "PairQuerySummaryStatistics")
		.def_readonly("structures_considered", &PairQuerySummaryStatistics::structures_considered)
		.def_readonly("fragments_considered", &PairQuerySummaryStatistics::fragments_considered)
		.def_readonly("fragments_expanded", &PairQuerySummaryStatistics::fragments_expanded)
		.def_readonly("pairs_considered", &PairQuerySummaryStatistics::pairs_considered)
		.def_readonly("pairs_aligned", &PairQuerySummaryStatistics::pairs_aligned)
		.def_readonly("result_count", &PairQuerySummaryStatistics::result_count);

	py::class_< StructureSingleQuery >(m, "StructureSingleQuery")
		.def(py::init< ndarray::Array<SearchReal, 3, 3>, SearchReal >())
		.def_readonly("n_entry",  &StructureSingleQuery::n_entry)
		.def_readonly("q_buffer",  &StructureSingleQuery::q_buffer)
		.def_readonly("c_per_entry",  &StructureSingleQuery::c_per_entry)
		.def_readonly("rmsd_tolerance",  &StructureSingleQuery::rmsd_tolerance);

	py::class_< SingleQueryExecutor >(m, "SingleQueryExecutor")
		.def(py::init<StructureSingleQuery const &>())
		.def("execute", &SingleQueryExecutor::execute )
		.def("execute_structure", &SingleQueryExecutor::execute_structure )
		.def_property_readonly("query_results",
		[](SingleQueryExecutor & src) {
			return py::array_t<StructureSingleQueryResult>(
			src.query_results.size(), &src.query_results[0]);
		});

	py::class_<SingleQuerySummaryStatistics >(m, "SingleQuerySummaryStatistics")
		.def_readonly("structures_considered", &SingleQuerySummaryStatistics::structures_considered)
		.def_readonly("fragments_considered", &SingleQuerySummaryStatistics::fragments_considered)
		.def_readonly("result_count", &SingleQuerySummaryStatistics::result_count);
}

}}}
