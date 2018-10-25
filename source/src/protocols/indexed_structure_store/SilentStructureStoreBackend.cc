// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file ./src/protocols/indexed_structure_store/SilentStructureStoreBackend.cc
/// @brief SilentFile-backed protein structure store.
/// @author Alex Ford (fordas@uw.edu)
//
#include <basic/Tracer.hh>

#include <utility/file/file_sys_util.hh>
#include <utility/file/FileName.hh>
#include <utility/file/PathName.hh>
#include <utility/io/izstream.hh>
#include <utility/exit.hh>

#include <protocols/indexed_structure_store/SilentStructureStoreBackend.hh>
#include <protocols/indexed_structure_store/Datatypes.hh>
#include <protocols/indexed_structure_store/StructureStore.hh>
#include <protocols/indexed_structure_store/pose_utility.hh>

#include <core/pose/Pose.hh>

#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentStruct.hh>

#include <string>
#include <stdexcept>

namespace protocols
{
namespace indexed_structure_store
{

static basic::Tracer TR("core.indexed_structure_store.SilentStructureStoreBackend");

bool SilentStructureStoreBackend::can_load(std::string store_path) {
	if ( utility::file::file_extension(store_path) == "silent" ) {
		return true;
	} else if ( utility::file::file_extension(store_path) == "out" ) {
		return true;
	} else {
		return false;
	}
}

StructureStoreOP SilentStructureStoreBackend::load_store(std::string store_path)
{
	TR.Debug << "load_store: " << store_path << std::endl;

	// SilentFileData::read_file throws exceptions on errors while loading silent file.
	core::io::silent::SilentFileOptions file_options;
	file_options.read_through_errors(false);
	core::io::silent::SilentFileData file_data(file_options);
	file_data.read_file(store_path);

	TR.Debug << "load_store loading silent structures: " << file_data.size() << std::endl;

	std::vector<StructureEntry> structure_entries;
	std::vector<ResidueEntry> residue_entries;

	uint32_t structure_id = 0;

	core::pose::PoseOP work_pose(new core::pose::Pose());

	for ( std::string const & tag : file_data.tags() ) {
		TR.Debug << "load_store loading silent tag: " << tag << std::endl;

		StructureEntry structure_entry;
		structure_entry.id = structure_id;

		std::string structure_name(tag);
		structure_name.resize(32);
		TR.Debug << "load_store loading name: " << structure_name << std::endl;
		structure_name.copy(structure_entry.name, 32, 0);
		TR.Debug << "load_store chararray name: " << std::string(structure_entry.name) << std::endl;

		structure_entries.push_back(structure_entry);

		file_data[tag]->fill_pose(*work_pose);
		ndarray::Array<ResidueEntry, 1> pose_residue_entries = extract_residue_entries(*work_pose, true);
		for ( auto & r : pose_residue_entries ) {
			r.structure_id = structure_id;
		}
		TR.Debug << "load_store loaded residues: " << pose_residue_entries.getSize<0>() << std::endl;

		std::copy(
			pose_residue_entries.begin(), pose_residue_entries.end(),
			std::back_inserter(residue_entries)
		);

		structure_id += 1;
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

bool SilentStructureStoreBackend::can_write(std::string) {
	return false;
}

void SilentStructureStoreBackend::write_store(std::string, StructureStore &)
{
	utility_exit_with_message("SilentStructureStoreBackend does not support writes.");
}

}
}
