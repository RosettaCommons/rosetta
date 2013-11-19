// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet;
//
// This file is made available under the Rosetta Commons license.
// See http://www.rosettacommons.org/license
// (C) 199x-2007 University of Washington
// (C) 199x-2007 University of California Santa Cruz
// (C) 199x-2007 University of California San Francisco
// (C) 199x-2007 Johns Hopkins University
// (C) 199x-2007 University of North Carolina, Chapel Hill
// (C) 199x-2007 Vanderbilt University

/// @brief Declaration of a filter to extract constellations shared by multiple chains
/// @author Andrea Bazzoli

#ifndef INCLUDED_InterfaceFilter_hh
#define INCLUDED_InterfaceFilter_hh

#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <core/types.hh>

using core::pose::Pose;
using core::Size;

namespace devel {
namespace constel {

bool at_interface(Pose const& ps, utility::vector1<Size> const& cnl);

} // constel
} // devel

#endif
