// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/scoring/EnergiesCacheableDataType.hh
/// @brief  enum for the DataCache within the Energies class
/// @author

#ifndef INCLUDED_core_scoring_EnergiesCacheableDataType_hh
#define INCLUDED_core_scoring_EnergiesCacheableDataType_hh


namespace core {
namespace scoring {


// hold the enum within a descriptive namespace to avoid name collisions
namespace EnergiesCacheableDataType {

/// @brief enum indexing the types for the DataCache within the Energies class
enum Enum {
	ETABLE_NBLIST = 1, // for indexing into vector1
	ELEC_NBLIST,
    GEOM_SOLV_NBLIST,
    LK_POLARNONPOLAR_NBLIST,
	MM_LJ_INTER_NBLIST,
	MM_LJ_INTRA_NBLIST,
	HBOND_SET,
	H2O_HBOND_SET,
	ETABLE_TRIE_COLLECTION,
	HBOND_TRIE_COLLECTION,
	ELEC_TRIE_COLLECTION,
	MM_LJ_TRIE_COLLECTION,
	// *** IMPORTANT ***
	// The 'num_cacheable_data_types' below must be the last enum, and must
	// always be set equal to the (last-1) enum.  If you append a new enum
	// to the list, remember to change the value below!
	num_cacheable_data_types = MM_LJ_TRIE_COLLECTION
};

} // namespace EnergiesCacheableDataType


} // namespace scoring
} // namespace core


#endif /* INCLUDED_core_scoring_EnergiesCacheableDataType_HH */
