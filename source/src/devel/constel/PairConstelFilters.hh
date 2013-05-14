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

/// @brief Filters for constellations formed by a pair of residues.
/// @author Andrea Bazzoli

#ifndef INCLUDED_PairConstelFilters_hh
#define INCLUDED_PairConstelFilters_hh

#include <core/pose/Pose.fwd.hh>
#include <utility/vector1.fwd.hh>
#include <core/types.hh>

using core::pose::Pose;
using core::Size;

namespace devel {
namespace constel {

///
/// @brief Class to filter out constellations that cannot (putatively) be rescued
///		by compounds containing indole and a carboxylic group.
///
class FilterByIndoleCOO {

  public:
  static bool sat(Pose const& ps, utility::vector1<Size> const& cnl);
};


///
/// @brief Class to filter out constellations that cannot (putatively) be rescued
/// 	by tryptamine.
///
class FilterByTryptamine {

  public:
  static bool sat(Pose const& ps, utility::vector1<Size> const& cnl);
};


///
/// @brief Class to filter out constellations that cannot (putatively) be rescued
/// 	by amphetamine.
///
class FilterByAmphetamine {

  public:
  static bool sat(Pose const& ps, utility::vector1<Size> const& cnl);
};


///
/// @brief Class to filter out constellations that cannot (putatively) be rescued
/// 	by histamine.
///
class FilterByHistamine {

  public:
  static bool sat(Pose const& ps, utility::vector1<Size> const& cnl);
};


} // constel
} // devel 

#endif
