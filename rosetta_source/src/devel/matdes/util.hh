// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.
//
/// @file 
/// @brief 
/// @author Neil King ( neilking@uw.edu )
/// @author Javier Castellanos ( javiercv@uw.edu )
/// @author Jacob Bale (balej@uw.edu)

#ifndef INCLUDED_devel_matdes_util_hh
#define INCLUDED_devel_matdes_util_hh

// Package headers

// project headers
#include <core/pack/task/PackerTask.fwd.hh>
#include <core/pose/Pose.fwd.hh>
#include <core/types.hh>
#include <utility/vector1.hh>

// STL headers
#include <set>
#include <sstream>
#include <algorithm>
#include <iterator>

namespace devel {
namespace matdes {
using core::pose::Pose;
using core::Real;
using core::Size;

template<typename T>
std::string stringify_iterable(const T& iterable, std::string sep=" "){
		std::ostringstream ss;
		const char* c_sep = sep.c_str();
		std::copy(iterable.begin(), iterable.end(), std::ostream_iterator<typename T::value_type>(ss, c_sep));
		return ss.str();
}

void add_native_bias_constraints(Pose & pose, core::Real cst_weight, const std::set<Size>& design_pos);


utility::vector1<Real>
sidechain_sasa(Pose const & pose, Real probe_radius = 2.2);

std::set<Size>
pick_design_position(core::pose::Pose const & pose, Size nsub_bblock = 1, Real contact_dist=10.0, Real bblock_dist = 5.0, Real probe_radius = 2.2);

core::pose::Pose
get_neighbor_subs(core::pose::Pose const &pose, utility::vector1<Size> intra_subs);

utility::vector1<Size>
get_neighbor_sub_resis(core::pose::Pose const &pose, Real contact_dist=10.0, std::string sym_dof_name = "" );

} // devel
} // matdes
#endif
