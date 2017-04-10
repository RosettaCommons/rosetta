// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is protocolsoped by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file    FlexBBDesignMeanField.hh

/// @brief   Forward declarations for FlexBBDesignMeanField.
/// @author  arubenstein

#ifndef INCLUDED_protocols_mean_field_FlexBBDesignMeanField_FWD_HH
#define INCLUDED_protocols_mean_field_FlexBBDesignMeanField_FWD_HH

// Utility header
#include <utility/pointer/owning_ptr.hh>

namespace protocols {
namespace mean_field {

/// @brief
class FlexBBDesignMeanField;

typedef utility::pointer::shared_ptr<FlexBBDesignMeanField> FlexBBDesignMeanFieldOP;
typedef utility::pointer::shared_ptr<FlexBBDesignMeanField const> FlexBBDesignMeanFieldCOP;

}  // namespace mean_field
}  // namespace protocols

#endif  // INCLUDED_protocols_mean_field_FlexBBDesignMeanField_FWD_HH
