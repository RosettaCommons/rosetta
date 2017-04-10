// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    ResHashMap.fwd.hh

/// @brief   Forward declarations for ResHashMap - a data structure that stores Residue hash keys and corresponding rot_ind values
/// @author  Aliza Rubenstein (aliza.rubenstein@gmail.com)

#ifndef INCLUDED_protocols_mean_field_ResHashMap_FWD_HH
#define INCLUDED_protocols_mean_field_ResHashMap_FWD_HH

// Utility header
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace mean_field {

/// @brief hashmap used to uniquely identify rotamer across different backbones (i.e. backbone-independent rotamers)
class ResHashMap;

typedef utility::pointer::shared_ptr<ResHashMap> ResHashMapOP;
typedef utility::pointer::shared_ptr<ResHashMap const> ResHashMapCOP;

}  // namespace mean_field
}  // namespace protocols

#endif  // INCLUDED_protocols_mean_field_ResHashMap_FWD_HH
