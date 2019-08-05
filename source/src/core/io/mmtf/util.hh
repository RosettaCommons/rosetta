// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file core/io/mmtf/util.hh
/// @brief Functions for MMTF writing/reading
/// @author Daniel Farrell (danpf@uw.edu)


#ifndef INCLUDED_core_io_mmtf_util_HH
#define INCLUDED_core_io_mmtf_util_HH

// Package headers
#include <core/io/AtomInformation.hh>
#include <utility/vector0.hh>
#include <core/types.hh>

#include <ostream>
#include <sstream>
#include <tuple>
#include <mmtf.hpp>

namespace core {
namespace io {
namespace mmtf {



//  AtomInformation Group/Chain/Pose
typedef utility::vector0<core::io::AtomInformation> aiGroup;
typedef utility::vector0<aiGroup> aiChain;
typedef utility::vector0<aiChain> aiPose;
typedef utility::vector0<aiPose> aiModels;

/// @brief storage for indicies to something in a mmtf struct
///        this makes it easier to iterate the struct in relation to
///        AtomInformation vectors
struct sd_index {
	sd_index() {
		model_index = -1;
		chain_index = -1;
		group_index = -1;
		group_atom_index = -1;
	}

	sd_index(int32_t mi, int32_t ci, int32_t gi, int32_t gai) {
		model_index = mi;
		chain_index = ci;
		group_index = gi;
		group_atom_index = gai;
	}

	bool operator==(sd_index const & rhs) const;

	std::string
	to_string() const;

	int32_t model_index, chain_index, group_index, group_atom_index;
};

std::ostream & operator<<( std::ostream & os, core::io::mmtf::sd_index const & sd_i );


/// @brief use the sd_index struct to make iterating StructureData easier
std::map<core::Size, sd_index>
make_atom_num_to_sd_map(::mmtf::StructureData const & sd);


/// @brief simple comparator for AtomInformation when grouping by residue/grp
struct ai_cmp {
	bool operator()(AtomInformation const & A, AtomInformation const & B) const {
		return std::tie(A.chainID, A.resSeq) < std::tie(B.chainID, B.resSeq);
	}
};


}  // namespace mmtf
}  // namespace io
}  // namespace core

#endif  // INCLUDED_core_io_mmtf_util_HH
