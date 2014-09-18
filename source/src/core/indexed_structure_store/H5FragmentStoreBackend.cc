// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Alex Ford <fordas@uw.edu>

#ifdef USEHDF5

#include <string>
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

static thread_local basic::Tracer TR( "core.indexed_structure_store.H5FragmentStoreBackend" );

#ifdef ROSETTA_FLOAT
#define RealPredType PredType::NATIVE_FLOAT
#else
#define RealPredType PredType::NATIVE_DOUBLE
#endif

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

  FragmentStoreOP fragment_store = new FragmentStore(fragment_spec, store_dataspace.getSimpleExtentNpoints());

	store_dataset.read(
      &fragment_store->fragment_threshold_distances[0],
      FragmentThresholdDistanceEntryDatatype());
	store_dataset.read(
      &fragment_store->fragment_coordinates[0],
      FragmentCoordinateEntryDatatype(fragment_spec));

	return fragment_store;
}

}
}

#endif
