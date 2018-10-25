// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.
//
/// @file protocols/indexed_structure_store/H5StructureStoreBackend
/// @brief HDF5-backed protein structure store.
/// @author Alex Ford <fordas@uw.edu>
//

#ifdef USEHDF5

#include <basic/Tracer.hh>

#include <utility/file/PathName.hh>
#include <utility/file/file_sys_util.hh>

#include <protocols/indexed_structure_store/H5StructureStoreBackend.hh>
#include <protocols/indexed_structure_store/Datatypes.hh>
#include <protocols/indexed_structure_store/StructureStore.hh>

#include <string>
#include <stdexcept>

#include "H5Cpp.h"

namespace protocols
{
namespace indexed_structure_store
{

static basic::Tracer TR("core.indexed_structure_store.H5StructureStoreBackend");

H5::DataType H5StructureStoreBackend::StructureEntry_datatype()
{
	using namespace H5;
	CompType type(sizeof(StructureEntry));
	type.insertMember("id"   , HOFFSET(StructureEntry , id)   , PredType::NATIVE_UINT32);
	type.insertMember("name" , HOFFSET(StructureEntry , name) , StrType(H5::PredType::C_S1 , 32));

	return type;
}

H5::DataType H5StructureStoreBackend::ResidueEntry_datatype()
{
	using namespace H5;
	CompType type(sizeof(ResidueEntry));

	type.insertMember("id",           HOFFSET(ResidueEntry , structure_id) , PredType::NATIVE_UINT32);
	type.insertMember("resn",         HOFFSET(ResidueEntry , residue_id)   , PredType::NATIVE_UINT32);
	type.insertMember("bb",           HOFFSET(ResidueEntry , bb)           , ResidueBackboneEntry_datatype());
	type.insertMember("sc",           HOFFSET(ResidueEntry , sc)           , ResidueSidechainEntry_datatype());
	type.insertMember("orient",       HOFFSET(ResidueEntry , orient)       , ResidueOrientEntry_datatype());
	type.insertMember("chain_ending", HOFFSET(ResidueEntry , chain_ending) , PredType::NATIVE_HBOOL);

	return type;
}

H5::DataType H5StructureStoreBackend::ResidueBackboneEntry_datatype()
{
	using namespace H5;
	CompType type(sizeof(ResidueBackboneEntry));
	type.insertMember("phi"   , HOFFSET(ResidueBackboneEntry , phi)   , PredType::NATIVE_FLOAT);
	type.insertMember("psi"   , HOFFSET(ResidueBackboneEntry , psi)   , PredType::NATIVE_FLOAT);
	type.insertMember("omega" , HOFFSET(ResidueBackboneEntry , omega) , PredType::NATIVE_FLOAT);

	return type;
}

H5::DataType H5StructureStoreBackend::ResidueSidechainEntry_datatype()
{
	using namespace H5;
	CompType type(sizeof(ResidueSidechainEntry));
	type.insertMember("chi1" , HOFFSET(ResidueSidechainEntry , chi1) , PredType::NATIVE_FLOAT);
	type.insertMember("chi2" , HOFFSET(ResidueSidechainEntry , chi2) , PredType::NATIVE_FLOAT);
	type.insertMember("chi3" , HOFFSET(ResidueSidechainEntry , chi3) , PredType::NATIVE_FLOAT);
	type.insertMember("chi4" , HOFFSET(ResidueSidechainEntry , chi4) , PredType::NATIVE_FLOAT);
	type.insertMember("aa"   , HOFFSET(ResidueSidechainEntry , aa)   , PredType::FORTRAN_S1);

	return type;
}

H5::DataType H5StructureStoreBackend::ResidueOrientEntry_datatype()
{
	using namespace H5;
	CompType type(sizeof(ResidueOrientEntry));

	hsize_t dimension = 3;
	type.insertMember("N"  , HOFFSET(ResidueOrientEntry , N)  , ArrayType(PredType::NATIVE_FLOAT , 1 , &dimension));
	type.insertMember("C"  , HOFFSET(ResidueOrientEntry , C)  , ArrayType(PredType::NATIVE_FLOAT , 1 , &dimension));
	type.insertMember("CA" , HOFFSET(ResidueOrientEntry , CA) , ArrayType(PredType::NATIVE_FLOAT , 1 , &dimension));
	type.insertMember("O" , HOFFSET(ResidueOrientEntry  , O)  , ArrayType(PredType::NATIVE_FLOAT , 1 , &dimension));

	return type;
}

bool H5StructureStoreBackend::can_load(std::string store_path) {
	return utility::file::file_extension(store_path) == "h5";
}

StructureStoreOP H5StructureStoreBackend::load_store(std::string store_path)
{
	using namespace H5;
	TR.Info << "Loading backend: " << store_path << std::endl;

	H5::H5File target_file;
	target_file.openFile(store_path, H5F_ACC_RDONLY);

	StructureStoreOP structure_store(new StructureStore());

	TR.Debug << "Opening: /structures" << std::endl;
	DataSet structure_dataset(target_file.openDataSet("structures"));
	DataSpace structure_dataspace(structure_dataset.getSpace());

	TR.Debug << "Loading: /structures size:" << structure_dataspace.getSimpleExtentNpoints() << std::endl;
	structure_store->structure_entries = ndarray::allocate(structure_dataspace.getSimpleExtentNpoints());
	structure_dataset.read(structure_store->structure_entries.begin(), StructureEntry_datatype());

	TR.Debug << "Opening: /residues" << std::endl;
	DataSet residue_dataset(target_file.openDataSet("residues"));
	DataSpace residue_dataspace(residue_dataset.getSpace());

	TR.Debug << "Loading: /residues size:" << residue_dataspace.getSimpleExtentNpoints() << std::endl;
	structure_store->residue_entries = ndarray::allocate(residue_dataspace.getSimpleExtentNpoints());
	residue_dataset.read(structure_store->residue_entries.begin(), ResidueEntry_datatype());

	TR.Debug << "Validating entries." << std::endl;
	structure_store->check_entries();
	TR.Debug << "Entries valid." << std::endl;

	return structure_store;
}

bool H5StructureStoreBackend::can_write(std::string) {
	return false;
}

void H5StructureStoreBackend::write_store(std::string, StructureStore &)
{
	utility_exit_with_message("H5StructureStoreBackend does not support writes.");
}

}
}

#endif
