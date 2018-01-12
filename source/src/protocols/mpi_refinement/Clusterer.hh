// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/protocols/mpi_refinement/Clusterer.cc
/// @brief
/// @author Hahnbeom Park

#ifndef INCLUDED_protocols_mpi_refinement_Clusterer_hh
#define INCLUDED_protocols_mpi_refinement_Clusterer_hh

#include <protocols/mpi_refinement/util.hh>
#include <protocols/wum/SilentStructStore.hh>

#include <core/io/silent/SilentStruct.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/pose/Pose.hh>

namespace protocols {
namespace mpi_refinement {

class Clusterer
{
public:

	Clusterer();
	~Clusterer();

	void set_defaults();

	protocols::wum::SilentStructStore
	apply( protocols::wum::SilentStructStore structs,
		core::Size const ncluster,
		core::Real const dist_cut ) const;

private:

	bool
	get_distance( core::io::silent::SilentStructOP ss1,
		core::io::silent::SilentStructOP ss2,
		core::Real &distance,
		core::Real const dist_cut ) const;

	protocols::wum::SilentStructStore
	energy_sort_cluster( protocols::wum::SilentStructStore structs,
		core::Size const ncluster,
		core::Real const dist_cut = 2.0 ) const;

private:

	//core::Real simlimit_;
	std::string similarity_method_;
	std::string similarity_measure_;
	core::Real simtol_;
	std::string sim_replace_obj_;
	std::string method_;

}; // class Clusterer
}
}

#endif
