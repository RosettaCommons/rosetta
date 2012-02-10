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

#ifndef INCLUDED_devel_matdes_util_hh
#define INCLUDED_devel_matdes_util_hh

// Package headers

// project headers

namespace devel {
namespace matdes {

PackerTaskOP
make_interface_design_packertask(core:pose:Pose & pose);

void add_native_bias_constraints(Pose & pose, Real cst_weight, const std::set<Size>& design_pos);

utility::vector1<Real>
sidechain_sasa(Pose const & pose, Real probe_radius = 2.2);

std::set<Size>
pick_design_position(core::pose::Pose & pose, Real contact_dist=10.0, Real bblock_dist = 5.0, Real probe_radius = 2.2);

} // devel
} // matdes
#endif
