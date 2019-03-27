// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/indexed_structure_store/JSONStructureStoreBackend
/// @brief JSON-backed protein structure store.
/// @author Alex Ford <fordas@uw.edu>
//
#include <basic/Tracer.hh>
#include <boost/format.hpp>

#include <utility/file/file_sys_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.hh>
#include <utility/io/ozstream.hh>
#include <utility/exit.hh>

#include <iostream>
#include <fstream>

#include <protocols/indexed_structure_store/JSONStructureStoreBackend.hh>
#include <protocols/indexed_structure_store/Datatypes.hh>
#include <protocols/indexed_structure_store/Datatypes.json.hh>
#include <protocols/indexed_structure_store/StructureStore.hh>
#include <json.hpp>

#include <string>
#include <stdexcept>

namespace protocols
{
namespace indexed_structure_store
{
using nlohmann::json;

static basic::Tracer TR("core.indexed_structure_store.JSONStructureStoreBackend");

struct StructureRecord {
	StructureEntry structure;
	std::vector<ResidueEntry> residues;
};

inline void to_json(json& j, const StructureRecord& i) {
	j = json{
		{"structure" , i.structure},
		{"residues" , i.residues},
		};
}

inline void from_json(const json& j, StructureRecord& i) {
	i.structure = j["structure"].get<StructureEntry>();
	i.residues = j["residues"].get<std::vector<ResidueEntry>>();
}

bool JSONStructureStoreBackend::can_load(std::string store_path) {
	if ( utility::file::file_extension(store_path) == "json" ) {
		return true;
		//} else if ( utility::file::file_extension(store_path) == "msgpack" ) {
		//return true;
	} else {
		return false;
	}
}

StructureStoreOP JSONStructureStoreBackend::load_store(std::string store_path)
{
	TR.Debug << "load_store: " << store_path<< std::endl;

	std::ifstream data;
	data.open(store_path);

	std::vector<StructureEntry> structure_entries;
	std::vector<ResidueEntry> residue_entries;

	// json parser does not properly handle streaming reads on larger files
	// Update to 3.* release series should resolve this issue, but until
	// then we need to read-by-line
	// https://github.com/nlohmann/json/issues/367
	int line = -1;
	for ( std::string record_line; std::getline(data, record_line); ++line ) {
		TR.Debug << boost::format("load_store loading record line: %i pos: %i") % line % data.tellg() << std::endl;

		StructureRecord record = json::parse(record_line);
		TR.Debug << "load_store loaded record: " << json(record.structure) << std::endl;

		structure_entries.push_back(record.structure);
		std::copy(
			record.residues.begin(), record.residues.end(),
			std::back_inserter(residue_entries)
		);
	}

	StructureStoreOP structure_store(new StructureStore());

	TR.Debug << "load_store structures size: " << structure_entries.size() << std::endl;
	structure_store->structure_entries = ndarray::allocate(structure_entries.size());
	std::copy(
		structure_entries.begin(), structure_entries.end(),
		structure_store->structure_entries.begin());

	TR.Debug << "get_residues_store residues size: " << residue_entries.size() << std::endl;
	structure_store->residue_entries = ndarray::allocate(residue_entries.size());
	std::copy(
		residue_entries.begin(), residue_entries.end(),
		structure_store->residue_entries.begin());

	TR.Debug << "Validating entries." << std::endl;
	structure_store->check_entries();
	TR.Debug << "Entries valid." << std::endl;

	return structure_store;
}

bool JSONStructureStoreBackend::can_write(std::string store_path) {
	if ( utility::file::file_extension(store_path) == "json" ) {
		return true;
		//} else if ( utility::file::file_extension(store_path) == "msgpack" ) {
		//return true;
	} else {
		return false;
	}
}

void JSONStructureStoreBackend::write_store(std::string store_path, StructureStore & store)
{
	TR << "write_store: " << store_path<< std::endl;

	utility::io::ozstream outstream( store_path.c_str() );

	auto res_end = store.residue_entries.begin();
	for ( auto & structure : store.structure_entries ) {
		TR.Debug << "writing structure: " << json(structure) << std::endl;
		auto res_start = res_end;
		if ( res_start->structure_id != structure.id ) {
			TR.Debug << "structure: " << json(structure) << " residue: " << json(*res_start) << std::endl;
			utility_exit_with_message("Invalid store residue ordering.");
		}
		while ( res_end != store.residue_entries.end() && res_end->structure_id == structure.id ) {
			++res_end;
		}

		StructureRecord outrecord = {
			structure, std::vector<ResidueEntry>(res_start, res_end)};

		outstream << json(outrecord) << std::endl;
	}
}

}
}
