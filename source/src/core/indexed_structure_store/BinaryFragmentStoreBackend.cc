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

#include <fstream>
#include <string>
#include <stdexcept>

#include <boost/algorithm/string.hpp>

#include <utility/exit.hh>
#include <utility/json_spirit/json_spirit.h>
#include <utility/json_spirit/json_spirit_tools.hh>

#include <basic/Tracer.hh>

#include <numeric/xyzVector.hh>

#include <core/indexed_structure_store/FragmentStore.hh>
#include <core/indexed_structure_store/BinaryFragmentStoreBackend.hh>

namespace core
{
namespace indexed_structure_store
{


static THREAD_LOCAL basic::Tracer TR( "core.indexed_structure_store.BinaryFragmentStoreBackend" );

BinaryFragmentStoreBackend::BinaryFragmentStoreBackend(std::string target_path) :
	target_path_(target_path)
{
	TR.Info << "Loading backend: " << target_path << std::endl;
}

FragmentStoreOP BinaryFragmentStoreBackend::get_fragment_store(std::string store_name)
{
	TR.Debug << "get_fragment_store: " << target_path_ << " store: " << store_name << std::endl;

	std::string store_path = target_path_ + "/fragments/" + store_name;
	TR.Debug << "Opening: " << store_path << std::endl;

	std::string metadata_path = store_path + "/metadata.json";
	utility::json_spirit::mValue md_doc;
	std::ifstream md_stream(metadata_path.c_str(), std::ios::in | std::ios::binary);
	std::string md_string( (std::istreambuf_iterator<char>(md_stream)), std::istreambuf_iterator<char>() );
	TR.Debug << "Read store metadata: " << metadata_path << "\n" << md_string << std::endl;

	utility::json_spirit::mObject md_object = utility::json_spirit::read_mObject( md_string );

	FragmentSpecification fragment_spec;

	utility::json_spirit::mArray atoms_src = utility::json_spirit::get_mArray(md_object, "fragment_atoms");
	fragment_spec.fragment_atoms.resize(atoms_src.size());
	for ( numeric::Size i = 0; i < atoms_src.size(); i++ ) {
		fragment_spec.fragment_atoms[i] = atoms_src[i].get_str();
	}

	fragment_spec.fragment_length = utility::json_spirit::get_int(md_object, "fragment_length");

	TR.Debug << "Loaded fragment specification: " << fragment_spec << std::endl;

	numeric::Size num_entries = utility::json_spirit::get_int(md_object, "num_entries");
	TR.Debug << "Loading: " << store_path << " size:" << num_entries << std::endl;

	FragmentStoreOP fragment_store( new FragmentStore(fragment_spec, num_entries) );

	std::string coordinates_path = store_path + "/" + utility::json_spirit::get_string(md_object, "coordinates_file");
	std::fstream coordinates_file;
	coordinates_file.exceptions(std::ifstream::failbit | std::ifstream::badbit | std::ifstream::eofbit);
	coordinates_file.open(coordinates_path.c_str(),std::ios::in|std::ios::binary);

	for ( auto & fragment_coordinate : fragment_store->fragment_coordinates ) {
		// Stored precision may not be double
		double read_coordinate[3];
		coordinates_file.read((char*)&read_coordinate[0], sizeof(double) * 3);

		fragment_coordinate.x() = read_coordinate[0];
		fragment_coordinate.y() = read_coordinate[1];
		fragment_coordinate.z() = read_coordinate[2];
	}
	coordinates_file.close();

	std::string threshold_distance_path = store_path + "/" + utility::json_spirit::get_string(md_object, "threshold_distance_file");
	std::fstream threshold_distance_file;
	threshold_distance_file.exceptions(std::ifstream::failbit | std::ifstream::badbit | std::ifstream::eofbit);
	threshold_distance_file.open(threshold_distance_path.c_str(),std::ios::in|std::ios::binary);

	for ( double & fragment_threshold_distance : fragment_store->fragment_threshold_distances ) {
		// Stored precision may not be double
		double read_coordinate;
		threshold_distance_file.read((char*)&read_coordinate, sizeof(double));

		fragment_threshold_distance = read_coordinate;
	}
	threshold_distance_file.close();

	return fragment_store;
}

}
}
