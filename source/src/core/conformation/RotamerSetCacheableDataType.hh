// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/conformation/RotamerSetCacheableDataType.hh
/// @brief  enum for the DataCache within the RotamerSet class;
/// EnergyMethods in core/scoring forward-declare themselves here if they intend to cache
/// data within the the RotamerSet.
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_conformation_RotamerSetCacheableDataType_hh
#define INCLUDED_core_conformation_RotamerSetCacheableDataType_hh


namespace core {
namespace conformation {

// hold the enum within a descriptive namespace to avoid name collisions
namespace RotamerSetCacheableDataType {

/// @brief enum indexing the types for the DataCache within the RotamerSet class
enum Enum {
	LK_BALL_ROTAMER_SET_INFO = 1,
	FACTS_ROTAMER_SET_INFO,
	GEN_BORN_ROTAMER_SET_INFO,
	MULTIPOLE_ELEC_ROTAMER_SET_INFO,
	// *** IMPORTANT ***
	// The 'num_cacheable_data_types' below must be the last enum, and must
	// always be set equal to the (last-1) enum.  If you append a new enum
	// to the list, remember to change the value below!
	num_cacheable_data_types = MULTIPOLE_ELEC_ROTAMER_SET_INFO
};

} // namespace RotamerSetCacheableDataType
} // namespace conformation
} // namespace core


#endif
