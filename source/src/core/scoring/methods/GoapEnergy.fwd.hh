// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/methods/CartesianBondedEnergy.cc
/// @brief  C++ implementaion of GOAP(Generalized Orientation-dependent, All-atom statistical Potential)
///         by Zhou H & Skolnick J, Biophys J 2011, 101(8):2043-52.
/// @author Hahnbeom Park

#ifndef INCLUDED_core_scoring_methods_GoapEnergy_FWD_HH
#define INCLUDED_core_scoring_methods_GoapEnergy_FWD_HH

#include <utility/pointer/owning_ptr.hh>
#include <core/chemical/ResidueType.hh>

namespace core {
namespace scoring {
namespace methods {

  class GoapRsdType;
  typedef utility::pointer::shared_ptr< GoapRsdType > GoapRsdTypeOP;
  typedef utility::pointer::shared_ptr< GoapRsdType const > GoapRsdTypeCOP;

  class GoapEnergy;
  typedef utility::pointer::shared_ptr< GoapEnergy > GoapEnergyOP;

  //typedef std::map< chemical::AA const *, GoapRsdTypeOP > GoapRsdTypeMap;
  typedef std::map< std::string, GoapRsdTypeOP > GoapRsdTypeMap;

} // methods
} // scoring
} // core

#endif
