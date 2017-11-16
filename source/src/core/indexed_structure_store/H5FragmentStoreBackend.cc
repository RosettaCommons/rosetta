// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file
/// @brief
/// @author Alex Ford <fordas@uw.edu>

#ifdef USEHDF5

#include <string>
#include <map>
#include <stdexcept>

#include <boost/algorithm/string.hpp>

#include "H5Cpp.h"

#include <basic/Tracer.hh>

#include <numeric/xyzVector.hh>

#include <core/indexed_structure_store/FragmentStore.hh>
#include <core/indexed_structure_store/H5FragmentStoreBackend.hh>

namespace core
{
namespace indexed_structure_store
{

static basic::Tracer TR( "core.indexed_structure_store.H5FragmentStoreBackend" );

#ifdef ROSETTA_FLOAT
#define RealPredType PredType::NATIVE_FLOAT
#else
#define RealPredType PredType::NATIVE_DOUBLE
#endif
#define Int64PredType PredType::NATIVE_UINT64

H5::DataType H5FragmentStoreBackend::FragmentThresholdDistanceEntryDatatype()
{
	using namespace H5;
	CompType entry_type(sizeof(numeric::Real));
	entry_type.insertMember("threshold_distance", 0, RealPredType);
	return entry_type;
}

H5::DataType H5FragmentStoreBackend::FragmentCoordinateEntryDatatype(FragmentSpecification fragment_spec)
{
  using namespace H5;

	CompType entry_type(sizeof(numeric::xyzVector<numeric::Real>) * fragment_spec.coordinates_per_fragment());

  hsize_t array_length[2];
	array_length[0] = fragment_spec.coordinates_per_fragment();
	array_length[1] = 3;

	entry_type.insertMember("coordinates", 0, ArrayType(RealPredType, 2, array_length));

	return entry_type;
}

H5::DataType H5FragmentStoreBackend::FragmentInt64GroupEntryDatatype(std::string group_field)
{
  using namespace H5;

	CompType entry_type(sizeof(numeric::Size));
	entry_type.insertMember(group_field, 0, Int64PredType);

	return entry_type;
}

H5::DataType H5FragmentStoreBackend::FragmentRealGroupEntryDatatype(std::string group_field)
{
  using namespace H5;

	CompType entry_type(sizeof(numeric::Size));
	entry_type.insertMember(group_field, 0, RealPredType);

	return entry_type;
}
H5::DataType H5FragmentStoreBackend::FragmentReal1PerResGroupEntryDatatype(std::string group_field,FragmentSpecification fragment_spec)
{
	using namespace H5;
	CompType entry_type(sizeof(numeric::Size)*fragment_spec.fragment_length);
	hsize_t array_length[1];
	array_length[0] = fragment_spec.fragment_length ;
	entry_type.insertMember(group_field, 0, ArrayType(RealPredType, 1, array_length));
	return entry_type;
}


H5::DataType H5FragmentStoreBackend::FragmentString1PerResGroupEntryDatatype(std::string group_field,FragmentSpecification fragment_spec)
{
	using namespace H5;
	CompType entry_type(sizeof(char)*fragment_spec.fragment_length);
	hsize_t array_length[1];
	array_length[0] = fragment_spec.fragment_length ;
	H5::StrType string_type(H5::PredType::C_S1, 1);
	entry_type.insertMember(group_field, 0, ArrayType(string_type, 1, array_length));
	return entry_type;
}

H5::DataType H5FragmentStoreBackend::FragmentString5PerResGroupEntryDatatype(std::string group_field,FragmentSpecification fragment_spec)
{
	using namespace H5;
	CompType entry_type(sizeof(char)*fragment_spec.fragment_length*5);
	hsize_t array_length[1];
	array_length[0] = fragment_spec.fragment_length ;
	H5::StrType string_type(H5::PredType::C_S1, 5);
	entry_type.insertMember(group_field, 0, ArrayType(string_type, 1, array_length));
	return entry_type;
}


H5FragmentStoreBackend::H5FragmentStoreBackend(std::string target_filename)
{
	using namespace H5;
	TR.Info << "Loading backend: " << target_filename << std::endl;

	target_filename_ = target_filename;
	target_file_.openFile(target_filename, H5F_ACC_RDONLY);
}

FragmentStoreOP H5FragmentStoreBackend::get_fragment_store(std::string store_name)
{
	using namespace H5;

	TR.Debug << "get_fragment_store: " << target_filename_ << " store: " << store_name << std::endl;

	std::string store_path = "/fragments/" + store_name;

	TR.Debug << "Opening: " << store_path << std::endl;
	DataSet store_dataset(target_file_.openDataSet(store_path));

	FragmentSpecification fragment_spec;

  std::string raw_atoms ("");
	TR.Debug << "Reading attribute: fragment_atoms" << std::endl;
	Attribute atoms_attribute = store_dataset.openAttribute("fragment_atoms");
	atoms_attribute.read(atoms_attribute.getStrType(), raw_atoms);
	TR.Debug << "Read attribute: fragment_atoms value:" << raw_atoms << std::endl;
	fragment_spec.fragment_atoms.resize(0);
  boost::split(fragment_spec.fragment_atoms, raw_atoms, boost::is_any_of(","));
	for (numeric::Size i = 0; i < fragment_spec.fragment_atoms.size(); i++)
	{
		TR.Debug << "Parsed fragment atom:" << fragment_spec.fragment_atoms[i] << std::endl;
	}

	TR.Debug << "Reading attribute: fragment_length" << std::endl;
	Attribute length_attribute = store_dataset.openAttribute("fragment_length");
	unsigned int fragment_length;
	length_attribute.read(PredType::NATIVE_UINT, &fragment_length);
	TR.Debug << "Read attribute: fragment_length value:" << fragment_length << std::endl;

	fragment_spec.fragment_length = fragment_length;
	TR.Debug << "Loaded attribute: fragment_length value:" << fragment_spec.fragment_length << std::endl;

	DataSpace store_dataspace(store_dataset.getSpace());
	TR.Debug << "Loading: " << store_path << " size:" << store_dataspace.getSimpleExtentNpoints() << std::endl;

	FragmentStoreOP fragment_store = FragmentStoreOP(new FragmentStore(fragment_spec, store_dataspace.getSimpleExtentNpoints()));

	store_dataset.read(&fragment_store->fragment_threshold_distances[0],FragmentThresholdDistanceEntryDatatype());

	store_dataset.read(
      &fragment_store->fragment_coordinates[0],
      FragmentCoordinateEntryDatatype(fragment_spec));
	return fragment_store;
}

void H5FragmentStoreBackend::append_to_fragment_store(FragmentStoreOP fragment_store, std::string store_name, std::string group_field, std::string group_type){
	using namespace H5;
	std::string store_path = "/fragments/" + store_name;
	TR.Debug << "Re-Opening: " << store_path << std::endl;
	DataSet store_dataset(target_file_.openDataSet(store_path));
	DataSpace store_dataspace(store_dataset.getSpace());
	TR.Debug <<"Loading group type: " << group_field << " of type " << group_type << std::endl;
	if(group_type != "int64" && group_type != "real" && group_type != "char_per_residue" && group_type != "real_per_residue" && group_type != "five_char_per_residue")
		utility_exit_with_message(group_type + " is not a valid entry in the fragment store. Currently only int64,real,char and 5char are implemented");
	if(group_type =="int64"){
		std::vector<numeric::Size> int64_group;
		int64_group.resize(store_dataspace.getSimpleExtentNpoints());
		store_dataset.read(
		&int64_group[0],
		FragmentInt64GroupEntryDatatype(group_field));
		fragment_store->int64_groups.insert(std::pair<std::string,std::vector <numeric::Size> > (group_field,int64_group));

		}
	if(group_type == "real"){
		std::vector<numeric::Real> real_group;
		real_group.resize(store_dataspace.getSimpleExtentNpoints());
		store_dataset.read(
		&real_group[0],
		FragmentRealGroupEntryDatatype(group_field));
		fragment_store->real_groups.insert(std::pair<std::string,std::vector <numeric::Real> > (group_field,real_group));
		}
	if(group_type == "real_per_residue"){
		std::vector<numeric::Real> real_group;
		std::vector<std::vector <numeric::Real> > real_group_processed;
		real_group.resize(store_dataspace.getSimpleExtentNpoints()*fragment_store->fragment_specification.coordinates_per_fragment());
		store_dataset.read(
		&real_group[0],
		FragmentReal1PerResGroupEntryDatatype(group_field,fragment_store->fragment_specification));
		std::vector<numeric::Real>::iterator begin_itr,end_itr;
		begin_itr = real_group.begin();
		for(numeric::Size ii=0; ii<(numeric::Size)store_dataspace.getSimpleExtentNpoints(); ++ii){
			end_itr=begin_itr+fragment_store->fragment_specification.coordinates_per_fragment();
			std::vector<numeric::Real> numericSplit(begin_itr,end_itr);
			real_group_processed.push_back(numericSplit);
			begin_itr=begin_itr+fragment_store->fragment_specification.coordinates_per_fragment();
		}
		fragment_store->realVector_groups.insert(std::pair<std::string,std::vector<std::vector<numeric::Real> > > (group_field,real_group_processed));
	}
	if(group_type == "char_per_residue"){
		//1 residue string with same length as the number of frags.
		std::vector<char> char_per_residue_group;
		std::vector<std::string> char_per_residue_group_processed;
		char_per_residue_group.resize(store_dataspace.getSimpleExtentNpoints()*fragment_store->fragment_specification.coordinates_per_fragment());
		store_dataset.read(&char_per_residue_group[0],FragmentString1PerResGroupEntryDatatype(group_field,fragment_store->fragment_specification));
		std::vector<char>::iterator begin_itr,end_itr;
		begin_itr = char_per_residue_group.begin();
		for(numeric::Size ii=0; ii<(numeric::Size)store_dataspace.getSimpleExtentNpoints(); ++ii){
			end_itr=begin_itr+fragment_store->fragment_specification.coordinates_per_fragment();
			std::string aa(begin_itr,end_itr);
			char_per_residue_group_processed.push_back(aa);
			begin_itr=begin_itr+fragment_store->fragment_specification.coordinates_per_fragment();
		}
		fragment_store->string_groups.insert(std::pair<std::string,std::vector<std::string> > (group_field,char_per_residue_group_processed));
	}
	if(group_type == "five_char_per_residue"){
		//1 residue string with same length as the number of frags.
		std::vector<char> char_per_residue_group;
		std::vector<std::string> char_per_residue_group_processed;
		char_per_residue_group.resize(store_dataspace.getSimpleExtentNpoints()*fragment_store->fragment_specification.coordinates_per_fragment()*5);
		store_dataset.read(&char_per_residue_group[0],FragmentString5PerResGroupEntryDatatype(group_field,fragment_store->fragment_specification));
		std::vector<char>::iterator begin_itr,end_itr;
		begin_itr = char_per_residue_group.begin();
		for(numeric::Size ii=0; ii<(numeric::Size)store_dataspace.getSimpleExtentNpoints(); ++ii){
			end_itr=begin_itr+5;
			std::string name(begin_itr,end_itr);
			char_per_residue_group_processed.push_back(name);
			begin_itr=begin_itr+fragment_store->fragment_specification.coordinates_per_fragment()*5;//skip the other 8 residues
		}
		fragment_store->string_groups.insert(std::pair<std::string,std::vector<std::string> > (group_field,char_per_residue_group_processed));
		}
	}

}
}

#endif
